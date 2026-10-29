"""CSD entry adapter for the release-scale CoordRep audit.

The adapter deliberately never invents oxidation, charge, spin, disorder
variants, or lattice translations.  It emits the structural information that
is present in the CSD molecular graph and assigns an explicit partial or
ambiguous status when the source object does not support a unique state.
"""

from __future__ import annotations

import hashlib
import math
from collections import Counter, defaultdict, deque
from dataclasses import dataclass, replace
from itertools import combinations
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple

import numpy as np

from coordrep.geometry.shape import ShapeCalculator

from .metal_policy import is_target_metal_symbol, metal_blocks
from .records import (
    CoordinationIncidence,
    DonorGroupRelation,
    DonorPairRelation,
    DonorGroup,
    EntryRecord,
    ExternalCenterRelation,
    ExternalDonorRelation,
    MetalRelation,
    MetalSiteRecord,
    validate_entry_record,
)


HAPTIC_MARKERS = ("pi", "deloc")
MAX_CONFIRMED_COLLECTIVE_SITE_MEMBERS = 12
PROTOCOL_ID = "csd-all-metal-structure-audit-20260825-v1"


class CSDRecordError(RuntimeError):
    """A controlled failure that must be written to the failure ledger."""


@dataclass(frozen=True)
class _Bond:
    a: int
    b: int
    bond_type: str
    is_haptic: bool


@dataclass
class _GroupDraft:
    kind: str
    atoms: Tuple[int, ...]
    targets: Set[int]
    bond_types: Set[str]
    incidences: Dict[int, Dict[int, Set[str]]] = None
    incidence_roles: Dict[int, Dict[int, Set[str]]] = None
    evidence: str = "CSD_native_bond"
    coordination_role: str = "sigma"
    issues: Set[str] = None

    def __post_init__(self) -> None:
        if self.issues is None:
            self.issues = set()
        if self.incidences is None:
            self.incidences = defaultdict(lambda: defaultdict(set))
        if self.incidence_roles is None:
            self.incidence_roles = defaultdict(lambda: defaultdict(set))


def _bool_attr(obj: Any, name: str) -> bool:
    try:
        return bool(getattr(obj, name))
    except Exception:
        return False


def _is_metal_atom(atom: Any) -> bool:
    """Return membership in the frozen CCDC-derived all-metal domain."""

    return is_target_metal_symbol(getattr(atom, "atomic_symbol", "?"))


def _bond_type(bond: Any) -> str:
    try:
        value = getattr(bond, "bond_type")
    except Exception:
        value = None
    text = "unknown" if value is None else str(value)
    return " ".join(text.split()) or "unknown"


def _is_haptic_type(value: str) -> bool:
    lowered = value.lower()
    return any(marker in lowered for marker in HAPTIC_MARKERS)


def _formal_charge(atom: Any) -> str:
    for name in ("formal_charge", "charge"):
        try:
            value = getattr(atom, name)
        except Exception:
            continue
        if value is not None:
            return str(value)
    return "unk"


def _coordinates(atom: Any) -> Optional[np.ndarray]:
    # Fractional coordinates are never a Cartesian fallback: using them as
    # angstrom coordinates corrupts distances, angles, and CShM in non-orthogonal
    # cells.  Fractional values are handled only by the PBC helpers below.
    for name in ("coordinates",):
        try:
            value = getattr(atom, name)
        except Exception:
            continue
        if value is None:
            continue
        try:
            if all(hasattr(value, axis) for axis in ("x", "y", "z")):
                arr = np.asarray([value.x, value.y, value.z], dtype=float)
            else:
                arr = np.asarray(list(value), dtype=float)
            if arr.shape == (3,) and np.all(np.isfinite(arr)):
                return arr
        except Exception:
            continue
    return None


def _fractional_coordinates(atom: Any) -> Optional[np.ndarray]:
    try:
        value = atom.fractional_coordinates
    except Exception:
        return None
    if value is None:
        return None
    try:
        if all(hasattr(value, axis) for axis in ("x", "y", "z")):
            arr = np.asarray([value.x, value.y, value.z], dtype=float)
        else:
            arr = np.asarray(list(value), dtype=float)
    except Exception:
        return None
    if arr.shape != (3,) or not np.all(np.isfinite(arr)):
        return None
    return arr


def _relative_image_delta(metal: Any, donor: Any) -> Optional[Tuple[int, int, int]]:
    metal_fractional = _fractional_coordinates(metal)
    donor_fractional = _fractional_coordinates(donor)
    if metal_fractional is None or donor_fractional is None:
        return None
    epsilon = 1.0e-8
    metal_lift = np.floor(metal_fractional + epsilon).astype(int)
    donor_lift = np.floor(donor_fractional + epsilon).astype(int)
    return tuple(int(value) for value in donor_lift - metal_lift)


def _periodic_orbit_key(indices: Sequence[int], atoms: Sequence[Any]) -> str:
    payload = []
    for index in indices:
        fractional = _fractional_coordinates(atoms[index])
        if fractional is None:
            return ""
        wrapped = np.mod(fractional, 1.0)
        payload.append((
            str(getattr(atoms[index], "atomic_symbol", "?")),
            tuple(round(float(value), 6) for value in wrapped),
        ))
    return _short_hash(tuple(sorted(payload)))


def _short_hash(value: Any, length: int = 20) -> str:
    return hashlib.sha256(repr(value).encode("utf-8")).hexdigest()[:length]


def _collect_graph(molecule: Any) -> Tuple[List[Any], List[_Bond], List[List[Tuple[int, str]]]]:
    atoms = list(molecule.atoms)
    if not atoms:
        raise CSDRecordError("EMPTY_MOLECULE")

    by_object = {id(atom): index for index, atom in enumerate(atoms)}
    by_source_index: Dict[int, List[int]] = defaultdict(list)
    by_label: Dict[str, List[int]] = defaultdict(list)
    for index, atom in enumerate(atoms):
        try:
            by_source_index[int(atom.index)].append(index)
        except Exception:
            pass
        by_label[str(getattr(atom, "label", ""))].append(index)

    def atom_index(atom: Any) -> Optional[int]:
        found = by_object.get(id(atom))
        if found is not None:
            return found
        try:
            source_candidates = by_source_index.get(int(atom.index), [])
        except Exception:
            source_candidates = []
        if len(source_candidates) == 1:
            return source_candidates[0]
        label = str(getattr(atom, "label", ""))
        candidates = by_label.get(label, [])
        if len(candidates) == 1:
            return candidates[0]
        return None

    raw_bonds: Iterable[Any]
    try:
        raw_bonds = list(molecule.bonds)
    except Exception:
        seen_objects: Dict[int, Any] = {}
        for atom in atoms:
            try:
                for bond in atom.bonds:
                    seen_objects[id(bond)] = bond
            except Exception:
                pass
        raw_bonds = list(seen_objects.values())

    bonds: List[_Bond] = []
    seen: Set[Tuple[int, int, str]] = set()
    for bond in raw_bonds:
        try:
            endpoints = list(bond.atoms)
        except Exception:
            continue
        if len(endpoints) != 2:
            continue
        a = atom_index(endpoints[0])
        b = atom_index(endpoints[1])
        if a is None or b is None or a == b:
            continue
        bt = _bond_type(bond)
        key = (min(a, b), max(a, b), bt)
        if key in seen:
            continue
        seen.add(key)
        bonds.append(_Bond(key[0], key[1], bt, _is_haptic_type(bt)))

    adjacency: List[List[Tuple[int, str]]] = [[] for _ in atoms]
    for bond in bonds:
        adjacency[bond.a].append((bond.b, bond.bond_type))
        adjacency[bond.b].append((bond.a, bond.bond_type))
    return atoms, bonds, adjacency


def _wl_colours(
    atoms: Sequence[Any], adjacency: Sequence[Sequence[Tuple[int, str]]], rounds: int = 4
) -> List[str]:
    colours = [
        _short_hash((
            str(getattr(atom, "atomic_symbol", "?")),
            _formal_charge(atom),
            len(adjacency[index]),
        ))
        for index, atom in enumerate(atoms)
    ]
    for _ in range(rounds):
        colours = [
            _short_hash((
                colours[index],
                tuple(sorted((bond_type, colours[other]) for other, bond_type in neighbours)),
            ))
            for index, neighbours in enumerate(adjacency)
        ]
    return colours


def _nonmetal_components(
    atoms: Sequence[Any],
    adjacency: Sequence[Sequence[Tuple[int, str]]],
    colours: Sequence[str],
) -> Tuple[Dict[int, int], Dict[int, str]]:
    nonmetals = {
        index
        for index, atom in enumerate(atoms)
        if not _is_metal_atom(atom)
    }
    component_by_atom: Dict[int, int] = {}
    component_keys: Dict[int, str] = {}
    component_id = 0
    for start in sorted(nonmetals):
        if start in component_by_atom:
            continue
        queue = deque([start])
        members: List[int] = []
        component_by_atom[start] = component_id
        while queue:
            current = queue.popleft()
            members.append(current)
            for other, _ in adjacency[current]:
                if other not in nonmetals or other in component_by_atom:
                    continue
                component_by_atom[other] = component_id
                queue.append(other)
        elements = sorted(str(getattr(atoms[i], "atomic_symbol", "?")) for i in members)
        internal_edges = []
        member_set = set(members)
        for i in members:
            for j, bt in adjacency[i]:
                if j in member_set and i < j:
                    endpoint_a = (
                        str(getattr(atoms[i], "atomic_symbol", "?")),
                        colours[i],
                    )
                    endpoint_b = (
                        str(getattr(atoms[j], "atomic_symbol", "?")),
                        colours[j],
                    )
                    left, right = sorted((endpoint_a, endpoint_b))
                    internal_edges.append((left, bt, right))
        component_keys[component_id] = _short_hash((tuple(elements), tuple(sorted(internal_edges))))
        component_id += 1
    return component_by_atom, component_keys


def _connected_components(nodes: Set[int], adjacency: Sequence[Sequence[Tuple[int, str]]]) -> List[Tuple[int, ...]]:
    components: List[Tuple[int, ...]] = []
    unseen = set(nodes)
    while unseen:
        start = min(unseen)
        unseen.remove(start)
        queue = deque([start])
        members = [start]
        while queue:
            current = queue.popleft()
            for other, _ in adjacency[current]:
                if other in unseen:
                    unseen.remove(other)
                    queue.append(other)
                    members.append(other)
        components.append(tuple(sorted(members)))
    return components


def _build_group_drafts(
    atoms: Sequence[Any],
    bonds: Sequence[_Bond],
    adjacency: Sequence[Sequence[Tuple[int, str]]],
    metal_indices: Sequence[int],
) -> Tuple[
    List[_GroupDraft],
    Set[Tuple[int, int, str]],
    Set[Tuple[int, int, str]],
]:
    metal_set = set(metal_indices)
    sigma: Dict[int, _GroupDraft] = {}
    pi_by_metal: Dict[int, Set[int]] = defaultdict(set)
    pi_bond_type: Dict[Tuple[int, int], Set[str]] = defaultdict(set)
    direct_mm: Set[Tuple[int, int, str]] = set()
    external_mm: Set[Tuple[int, int, str]] = set()

    for bond in bonds:
        a_metal = bond.a in metal_set
        b_metal = bond.b in metal_set
        if a_metal and b_metal:
            direct_mm.add((min(bond.a, bond.b), max(bond.a, bond.b), bond.bond_type))
            continue
        if a_metal == b_metal:
            continue
        metal_index = bond.a if a_metal else bond.b
        donor_index = bond.b if a_metal else bond.a
        if _is_metal_atom(atoms[donor_index]):
            external_mm.add((metal_index, donor_index, bond.bond_type))
            continue
        if bond.is_haptic:
            pi_by_metal[metal_index].add(donor_index)
            pi_bond_type[(metal_index, donor_index)].add(bond.bond_type)
        else:
            draft = sigma.setdefault(
                donor_index,
                _GroupDraft(
                    kind="atom",
                    atoms=(donor_index,),
                    targets=set(),
                    bond_types=set(),
                    coordination_role="sigma",
                ),
            )
            draft.targets.add(metal_index)
            draft.bond_types.add(bond.bond_type)
            draft.incidences[metal_index][donor_index].add(bond.bond_type)
            draft.incidence_roles[metal_index][donor_index].add("sigma")

    pi: Dict[Tuple[int, ...], _GroupDraft] = {}
    for metal_index, donor_indices in pi_by_metal.items():
        for component in _connected_components(set(donor_indices), adjacency):
            if len(component) == 1:
                donor_index = component[0]
                draft = sigma.setdefault(
                    donor_index,
                    _GroupDraft(
                        kind="atom",
                        atoms=(donor_index,),
                        targets=set(),
                        bond_types=set(),
                        coordination_role="pi",
                    ),
                )
                draft.targets.add(metal_index)
                draft.bond_types.update(pi_bond_type[(metal_index, donor_index)])
                draft.incidences[metal_index][donor_index].update(
                    pi_bond_type[(metal_index, donor_index)]
                )
                draft.incidence_roles[metal_index][donor_index].add("pi")
                continue
            is_confirmed_collective = (
                len(component) <= MAX_CONFIRMED_COLLECTIVE_SITE_MEMBERS
            )
            site_kind = (
                "pi_fragment" if is_confirmed_collective else "candidate_atom_set"
            )
            draft = pi.setdefault(
                component,
                _GroupDraft(
                    kind=site_kind,
                    atoms=component,
                    targets=set(),
                    bond_types=set(),
                    coordination_role=("pi" if is_confirmed_collective else "pi_candidate"),
                ),
            )
            if not is_confirmed_collective:
                draft.issues.add("COLLECTIVE_SITE_MEMBERSHIP_OUT_OF_DOMAIN")
            draft.targets.add(metal_index)
            for donor_index in component:
                draft.bond_types.update(pi_bond_type[(metal_index, donor_index)])
                if pi_bond_type[(metal_index, donor_index)]:
                    draft.incidences[metal_index][donor_index].update(
                        pi_bond_type[(metal_index, donor_index)]
                    )
                    draft.incidence_roles[metal_index][donor_index].add("pi")

    return list(sigma.values()) + list(pi.values()), direct_mm, external_mm


def _metal_sort_key(
    metal_index: int,
    atoms: Sequence[Any],
    colours: Sequence[str],
    drafts: Sequence[_GroupDraft],
    external_mm: Sequence[Tuple[int, int, str]],
) -> Tuple[Any, ...]:
    local = [draft for draft in drafts if metal_index in draft.targets]
    donor_signature = []
    for draft in local:
        donor_signature.append((
            draft.kind,
            tuple(sorted(str(getattr(atoms[i], "atomic_symbol", "?")) for i in draft.atoms)),
            len(draft.atoms),
            len(draft.targets),
            tuple(sorted(colours[i] for i in draft.atoms)),
        ))
    external_signature = tuple(sorted(
        (
            str(getattr(atoms[partner], "atomic_symbol", "?")),
            colours[partner],
            bond_type,
        )
        for target, partner, bond_type in external_mm
        if target == metal_index
    ))
    return (
        str(getattr(atoms[metal_index], "atomic_symbol", "?")),
        colours[metal_index],
        len(local),
        tuple(sorted(donor_signature)),
        external_signature,
    )


def _group_coordinate(
    draft: _GroupDraft, atoms: Sequence[Any], target: Optional[int] = None
) -> Optional[np.ndarray]:
    atom_indices = draft.atoms
    if target is not None and target in draft.incidences:
        atom_indices = tuple(sorted(draft.incidences[target]))
    points = [_coordinates(atoms[index]) for index in atom_indices]
    if not points or any(point is None for point in points):
        return None
    return np.mean(np.asarray(points, dtype=float), axis=0)


def _local_geometry(
    metal_index: int,
    local_drafts: Sequence[Tuple[str, _GroupDraft]],
    atoms: Sequence[Any],
) -> Dict[str, Any]:
    metal_coord = _coordinates(atoms[metal_index])
    points = [_group_coordinate(draft, atoms, metal_index) for _, draft in local_drafts]
    if metal_coord is None or any(point is None for point in points):
        return {
            "geometry_kind": "topology-only",
            "shape_best": "",
            "shape_values": (),
            "radial_profile": (),
            "angular_profile": (),
            "issues": ("COORDINATES_UNAVAILABLE",),
        }
    if not points:
        return {
            "geometry_kind": "topology-only",
            "shape_best": "",
            "shape_values": (),
            "radial_profile": (),
            "angular_profile": (),
            "issues": ("NO_DONOR_GROUPS",),
        }

    vectors = np.asarray(points, dtype=float) - metal_coord
    radii = np.linalg.norm(vectors, axis=1)
    radial = tuple(sorted(round(float(value), 4) for value in radii))
    angular: List[float] = []
    for i, j in combinations(range(len(vectors)), 2):
        denominator = float(radii[i] * radii[j])
        if denominator <= 1e-12:
            continue
        cosine = float(np.dot(vectors[i], vectors[j]) / denominator)
        angular.append(round(max(-1.0, min(1.0, cosine)), 6))
    angular_profile = tuple(sorted(angular))

    if 2 <= len(points) <= 6:
        try:
            result = ShapeCalculator(rounding_decimals=4).compute(
                np.asarray(points, dtype=float), metal_coord
            )
            values = tuple(float(value) for value in result.values_rounded)
            if values:
                best_index = min(range(len(values)), key=values.__getitem__)
                return {
                    "geometry_kind": "CShM",
                    "shape_best": str(result.ref_shapes[best_index]),
                    "shape_values": values,
                    "radial_profile": radial,
                    "angular_profile": angular_profile,
                    "issues": (),
                }
        except Exception:
            pass

    issue = "REFERENCE_SET_UNAVAILABLE" if len(points) > 6 else "CSHM_UNAVAILABLE"
    return {
        "geometry_kind": "distance-angle-profile",
        "shape_best": "",
        "shape_values": (),
        "radial_profile": radial,
        "angular_profile": angular_profile,
        "issues": (issue,),
    }


def adapt_csd_entry(
    entry: Any,
    source_index: int,
    *,
    source_release: str,
    protocol_id: str = PROTOCOL_ID,
) -> EntryRecord:
    """Convert one in-domain CSD entry into a typed audit record.

    Raises
    ------
    CSDRecordError
        If the entry is not a three-dimensional metal-containing record or a
        structural record cannot be constructed.  Callers must write such
        failures separately and must not count them as emitted records.
    """
    if not _bool_attr(entry, "has_3d_structure"):
        raise CSDRecordError("NO_3D_STRUCTURE")
    try:
        molecule = entry.molecule
    except Exception as exc:
        raise CSDRecordError("MOLECULE_ACCESS_FAILED") from exc
    if molecule is None:
        raise CSDRecordError("NO_MOLECULE")

    atoms, bonds, adjacency = _collect_graph(molecule)
    metal_indices = [
        index
        for index, atom in enumerate(atoms)
        if _is_metal_atom(atom)
    ]
    if not metal_indices:
        raise CSDRecordError("NO_IN_DOMAIN_METAL")

    colours = _wl_colours(atoms, adjacency)
    component_by_atom, component_keys = _nonmetal_components(
        atoms, adjacency, colours
    )
    drafts, direct_mm, external_mm = _build_group_drafts(
        atoms, bonds, adjacency, metal_indices
    )
    is_polymeric = bool(
        _bool_attr(entry, "is_polymeric") or _bool_attr(molecule, "is_polymeric")
    )

    metal_set = set(metal_indices)
    external_donor_edge_types: Dict[Tuple[int, int], Set[str]] = defaultdict(set)
    for bond in bonds:
        endpoint_pairs = ((bond.a, bond.b), (bond.b, bond.a))
        for donor_index, partner_index in endpoint_pairs:
            if donor_index in metal_set or partner_index in metal_set:
                continue
            if _is_metal_atom(atoms[donor_index]):
                continue
            if not _is_metal_atom(atoms[partner_index]):
                continue
            external_donor_edge_types[(donor_index, partner_index)].add(
                bond.bond_type
            )
            break

    def external_donor_context(donor_index: int) -> Tuple[Any, ...]:
        donor_coord = _coordinates(atoms[donor_index])
        contexts = []
        for (candidate_donor, partner_index), bond_types in external_donor_edge_types.items():
            if candidate_donor != donor_index:
                continue
            partner_coord = _coordinates(atoms[partner_index])
            distance = None
            if donor_coord is not None and partner_coord is not None:
                distance = round(float(np.linalg.norm(partner_coord - donor_coord)), 4)
            image_delta = (
                _relative_image_delta(atoms[donor_index], atoms[partner_index])
                if is_polymeric
                else (0, 0, 0)
            )
            contexts.append((
                str(getattr(atoms[partner_index], "atomic_symbol", "?")),
                colours[partner_index],
                tuple(sorted(bond_types)),
                image_delta is not None,
                image_delta if image_delta is not None else (0, 0, 0),
                distance,
            ))
        return tuple(sorted(contexts, key=repr))

    keyed_metals = [
        (_metal_sort_key(index, atoms, colours, drafts, tuple(external_mm)), index)
        for index in metal_indices
    ]
    metal_key_counts = Counter(key for key, _ in keyed_metals)
    ordered_metals = [index for _, index in sorted(keyed_metals, key=lambda item: (item[0], item[1]))]
    metal_labels = {index: f"M{rank + 1}" for rank, index in enumerate(ordered_metals)}
    exact_metal_order = len(ordered_metals) == 1 or all(count == 1 for count in metal_key_counts.values())

    def draft_intrinsic_key(draft: _GroupDraft) -> Tuple[Any, ...]:
        donors = tuple(sorted(
            (
                str(getattr(atoms[i], "atomic_symbol", "?")),
                colours[i],
                external_donor_context(i),
            )
            for i in draft.atoms
        ))
        targets = tuple(sorted(metal_labels[i] for i in draft.targets))
        target_geometry = []
        for target in sorted(draft.targets, key=lambda index: metal_labels[index]):
            metal_coord = _coordinates(atoms[target])
            group_coord = _group_coordinate(draft, atoms, target)
            group_distance = None
            if metal_coord is not None and group_coord is not None:
                group_distance = round(float(np.linalg.norm(group_coord - metal_coord)), 4)
            donor_incidence = []
            for donor, bond_types in draft.incidences.get(target, {}).items():
                donor_coord = _coordinates(atoms[donor])
                distance = None
                if metal_coord is not None and donor_coord is not None:
                    distance = round(float(np.linalg.norm(donor_coord - metal_coord)), 4)
                image_delta = (
                    _relative_image_delta(atoms[target], atoms[donor])
                    if is_polymeric
                    else (0, 0, 0)
                )
                donor_incidence.append((
                    str(getattr(atoms[donor], "atomic_symbol", "?")),
                    colours[donor],
                    tuple(sorted(bond_types)),
                    tuple(sorted(draft.incidence_roles[target][donor])),
                    image_delta,
                    distance,
                ))
            target_geometry.append((
                metal_labels[target],
                group_distance,
                tuple(sorted(donor_incidence, key=repr)),
            ))
        component_ids = {component_by_atom.get(index) for index in draft.atoms}
        component_ids.discard(None)
        component_signature = tuple(sorted(component_keys[cid] for cid in component_ids))
        return (
            draft.kind,
            donors,
            targets,
            tuple(sorted(draft.bond_types)),
            component_signature,
            tuple(target_geometry),
        )

    draft_keys = [draft_intrinsic_key(draft) for draft in drafts]
    draft_key_counts = Counter(draft_keys)
    exact_group_order = all(count == 1 for count in draft_key_counts.values())
    ordered_drafts = sorted(
        drafts,
        key=lambda draft: (repr(draft_intrinsic_key(draft)), tuple(draft.atoms)),
    )
    draft_labels = {id(draft): f"G{rank + 1}" for rank, draft in enumerate(ordered_drafts)}
    donor_groups: List[DonorGroup] = []
    incidences: List[CoordinationIncidence] = []
    donor_ordinals: Dict[int, Dict[int, int]] = {}
    exact_donor_order = True
    for draft in ordered_drafts:
        def donor_intrinsic_key(index: int) -> Tuple[Any, ...]:
            donor_coord = _coordinates(atoms[index])
            incidence_context = []
            for target in sorted(draft.targets, key=lambda value: metal_labels[value]):
                if index not in draft.incidences.get(target, {}):
                    continue
                metal_coord = _coordinates(atoms[target])
                distance = None
                if metal_coord is not None and donor_coord is not None:
                    distance = round(float(np.linalg.norm(donor_coord - metal_coord)), 4)
                image_delta = (
                    _relative_image_delta(atoms[target], atoms[index])
                    if is_polymeric
                    else (0, 0, 0)
                )
                incidence_context.append((
                    metal_labels[target],
                    tuple(sorted(draft.incidences[target][index])),
                    tuple(sorted(draft.incidence_roles[target][index])),
                    image_delta,
                    distance,
                ))
            internal_geometry = []
            for other in draft.atoms:
                if other == index:
                    continue
                other_coord = _coordinates(atoms[other])
                separation = None
                if donor_coord is not None and other_coord is not None:
                    separation = round(float(np.linalg.norm(donor_coord - other_coord)), 4)
                internal_geometry.append((
                    str(getattr(atoms[other], "atomic_symbol", "?")),
                    colours[other],
                    separation,
                ))
            return (
                str(getattr(atoms[index], "atomic_symbol", "?")),
                colours[index],
                tuple(incidence_context),
                external_donor_context(index),
                tuple(sorted(internal_geometry, key=repr)),
            )

        donor_keys = [donor_intrinsic_key(index) for index in draft.atoms]
        if any(count > 1 for count in Counter(donor_keys).values()):
            exact_donor_order = False
        ordered_donor_atoms = tuple(sorted(
            draft.atoms,
            key=lambda index: (repr(donor_intrinsic_key(index)), index),
        ))
        donor_elements = tuple(
            str(getattr(atoms[index], "atomic_symbol", "?")) for index in ordered_donor_atoms
        )
        topology_keys = tuple(colours[index] for index in ordered_donor_atoms)
        target_labels = tuple(sorted(metal_labels[index] for index in draft.targets))
        component_ids = {component_by_atom.get(index) for index in draft.atoms}
        component_ids.discard(None)
        ligand_key = ""
        if len(component_ids) == 1:
            ligand_key = component_keys[next(iter(component_ids))]
        elif len(component_ids) > 1:
            ligand_key = _short_hash(tuple(sorted(component_keys[cid] for cid in component_ids)))
            draft.issues.add("DONOR_GROUP_SPANS_COMPONENTS")

        distances = []
        group_coord = _group_coordinate(draft, atoms)
        if group_coord is not None:
            for target in draft.targets:
                metal_coord = _coordinates(atoms[target])
                if metal_coord is not None:
                    distances.append(round(float(np.linalg.norm(group_coord - metal_coord)), 4))

        donor_groups.append(DonorGroup(
            label=draft_labels[id(draft)],
            kind=draft.kind,
            donor_elements=donor_elements,
            topology_keys=topology_keys,
            target_metals=target_labels,
            hapticity=len(draft.atoms) if draft.coordination_role == "pi" else 1,
            bridge_degree=len(target_labels),
            external_donor_contact_count=sum(
                len(external_donor_context(index)) for index in draft.atoms
            ),
            ligand_component_key=ligand_key,
            periodic_orbit_key=(
                _periodic_orbit_key(ordered_donor_atoms, atoms)
                if is_polymeric
                else ""
            ),
            bond_types=tuple(sorted(draft.bond_types)),
            distance_profile=tuple(sorted(distances)),
            evidence=draft.evidence,
            issues=tuple(sorted(draft.issues)),
        ))
        ordinal_by_atom = {
            atom_index: ordinal for ordinal, atom_index in enumerate(ordered_donor_atoms)
        }
        donor_ordinals[id(draft)] = ordinal_by_atom
        for target_index, donor_map in draft.incidences.items():
            target_label = metal_labels[target_index]
            target_coord = _coordinates(atoms[target_index])
            for donor_index, bond_types in donor_map.items():
                donor_coord = _coordinates(atoms[donor_index])
                distance = None
                if target_coord is not None and donor_coord is not None:
                    distance = round(float(np.linalg.norm(donor_coord - target_coord)), 4)
                role = "+".join(sorted(draft.incidence_roles[target_index][donor_index]))
                image_delta = (0, 0, 0)
                image_delta_resolved = True
                if is_polymeric:
                    resolved_delta = _relative_image_delta(
                        atoms[target_index], atoms[donor_index]
                    )
                    if resolved_delta is None:
                        draft.issues.add("PBC_IMAGE_DELTA_UNAVAILABLE")
                        image_delta_resolved = False
                    else:
                        image_delta = resolved_delta
                incidences.append(CoordinationIncidence(
                    metal=target_label,
                    donor_group=draft_labels[id(draft)],
                    donor_ordinal=ordinal_by_atom[donor_index],
                    donor_element=str(getattr(atoms[donor_index], "atomic_symbol", "?")),
                    bond_role=role,
                    bond_type="+".join(sorted(bond_types)),
                    image_delta=image_delta,
                    image_delta_resolved=image_delta_resolved,
                    distance=distance,
                ))
        if donor_groups[-1].issues != tuple(sorted(draft.issues)):
            donor_groups[-1] = replace(
                donor_groups[-1], issues=tuple(sorted(draft.issues))
            )

    external_donor_relations: List[ExternalDonorRelation] = []
    for draft in ordered_drafts:
        group_label = draft_labels[id(draft)]
        for donor_index in draft.atoms:
            for (candidate_donor, partner_index), bond_types in external_donor_edge_types.items():
                if candidate_donor != donor_index:
                    continue
                donor_coord = _coordinates(atoms[donor_index])
                partner_coord = _coordinates(atoms[partner_index])
                distance = None
                if donor_coord is not None and partner_coord is not None:
                    distance = round(float(np.linalg.norm(partner_coord - donor_coord)), 4)
                image_delta = (0, 0, 0)
                image_delta_resolved = True
                if is_polymeric:
                    resolved_delta = _relative_image_delta(
                        atoms[donor_index], atoms[partner_index]
                    )
                    if resolved_delta is None:
                        image_delta_resolved = False
                        draft.issues.add(
                            "PBC_EXTERNAL_DONOR_IMAGE_DELTA_UNAVAILABLE"
                        )
                    else:
                        image_delta = resolved_delta
                external_donor_relations.append(ExternalDonorRelation(
                    donor_group=group_label,
                    donor_ordinal=donor_ordinals[id(draft)][donor_index],
                    partner_element=str(
                        getattr(atoms[partner_index], "atomic_symbol", "?")
                    ),
                    partner_topology_key=colours[partner_index],
                    bond_type="+".join(sorted(bond_types)),
                    distance=distance,
                    image_delta=image_delta,
                    image_delta_resolved=image_delta_resolved,
                ))
        group_position = next(
            index for index, group in enumerate(donor_groups)
            if group.label == group_label
        )
        if donor_groups[group_position].issues != tuple(sorted(draft.issues)):
            donor_groups[group_position] = replace(
                donor_groups[group_position], issues=tuple(sorted(draft.issues))
            )

    external_donor_relations.sort(key=lambda relation: (
        relation.donor_group,
        relation.donor_ordinal,
        relation.partner_element,
        relation.partner_topology_key,
        relation.bond_type,
        not relation.image_delta_resolved,
        relation.image_delta,
        relation.distance is None,
        relation.distance if relation.distance is not None else 0.0,
    ))

    incidences.sort(key=lambda item: (
        item.metal,
        item.donor_group,
        item.donor_ordinal,
        item.image_delta,
        item.bond_role,
        item.bond_type,
    ))

    external_center_relations: List[ExternalCenterRelation] = []
    external_metal_neighbours: Dict[int, Set[int]] = defaultdict(set)

    def external_relation_sort_key(
        item: Tuple[int, int, str]
    ) -> Tuple[Any, ...]:
        target_index, partner_index, bond_type = item
        target_coord = _coordinates(atoms[target_index])
        partner_coord = _coordinates(atoms[partner_index])
        distance = None
        if target_coord is not None and partner_coord is not None:
            distance = round(float(np.linalg.norm(partner_coord - target_coord)), 4)
        image_delta = (
            _relative_image_delta(atoms[target_index], atoms[partner_index])
            if is_polymeric
            else (0, 0, 0)
        )
        return (
            metal_labels[target_index],
            str(getattr(atoms[partner_index], "atomic_symbol", "?")),
            colours[partner_index],
            bond_type,
            image_delta is None,
            image_delta if image_delta is not None else (0, 0, 0),
            distance is None,
            distance if distance is not None else 0.0,
            partner_index,
        )

    for target_index, partner_index, bond_type in sorted(
        external_mm,
        key=external_relation_sort_key,
    ):
        target_coord = _coordinates(atoms[target_index])
        partner_coord = _coordinates(atoms[partner_index])
        distance = None
        if target_coord is not None and partner_coord is not None:
            distance = round(float(np.linalg.norm(partner_coord - target_coord)), 4)
        image_delta = (0, 0, 0)
        image_delta_resolved = True
        if is_polymeric:
            resolved_delta = _relative_image_delta(
                atoms[target_index], atoms[partner_index]
            )
            if resolved_delta is None:
                image_delta_resolved = False
            else:
                image_delta = resolved_delta
        external_center_relations.append(ExternalCenterRelation(
            metal=metal_labels[target_index],
            partner_element=str(getattr(atoms[partner_index], "atomic_symbol", "?")),
            partner_topology_key=colours[partner_index],
            bond_type=bond_type,
            distance=distance,
            image_delta=image_delta,
            image_delta_resolved=image_delta_resolved,
        ))
        external_metal_neighbours[target_index].add(partner_index)

    exact_label_free = exact_metal_order and exact_group_order and exact_donor_order
    canonical_status = "exact_label_free" if exact_label_free else "source_order_noncanonical"

    local_sites: List[MetalSiteRecord] = []
    direct_metal_neighbours: Dict[int, Set[int]] = defaultdict(set)
    for a, b, _ in direct_mm:
        direct_metal_neighbours[a].add(b)
        direct_metal_neighbours[b].add(a)
    high_cn = False
    any_haptic = any(group.kind == "pi_fragment" and group.hapticity > 1 for group in donor_groups)
    any_pi_candidate = any(group.kind == "candidate_atom_set" for group in donor_groups)
    try:
        crystal = entry.crystal
    except Exception:
        crystal = None
    crystal_disorder_object = None
    if crystal is not None:
        try:
            crystal_disorder_object = crystal.disorder
        except Exception:
            crystal_disorder_object = None
    has_disorder = bool(
        _bool_attr(entry, "has_disorder")
        or (crystal is not None and _bool_attr(crystal, "has_disorder"))
        or crystal_disorder_object is not None
    )
    periodic_translation_complete = bool(
        not is_polymeric
        or (
            bool(incidences or direct_mm or external_mm)
            and all(
                group.periodic_orbit_key
                and "PBC_IMAGE_DELTA_UNAVAILABLE" not in group.issues
                for group in donor_groups
            )
            and all(incidence.image_delta_resolved for incidence in incidences)
            and all(
                _relative_image_delta(atoms[a], atoms[b]) is not None
                for a, b, _ in direct_mm
            )
            and all(
                relation.image_delta_resolved
                for relation in external_center_relations
            )
            and all(
                relation.image_delta_resolved
                for relation in external_donor_relations
            )
        )
    )

    for metal_index in ordered_metals:
        label = metal_labels[metal_index]
        local_pairs = [
            (draft_labels[id(draft)], draft)
            for draft in ordered_drafts
            if metal_index in draft.targets
        ]
        local_incidences = [incidence for incidence in incidences if incidence.metal == label]
        local_group_labels = {incidence.donor_group for incidence in local_incidences}
        local_group_records = [group for group in donor_groups if group.label in local_group_labels]
        cn_site = len(local_group_records)
        cn_atom = len(local_incidences)
        eta_sum = (
            None
            if any(group.kind == "candidate_atom_set" for group in local_group_records)
            else sum(group.hapticity for group in local_group_records)
        )
        composition = Counter(
            incidence.donor_element for incidence in local_incidences
        )
        geometry = _local_geometry(metal_index, local_pairs, atoms)
        geometry_kind = geometry["geometry_kind"]
        site_issues = set(geometry["issues"])
        if cn_site > 6:
            high_cn = True
        if is_polymeric and not periodic_translation_complete:
            site_issues.add("PBC_TRANSLATION_UNRESOLVED")
        if has_disorder:
            site_issues.add("DISORDER_VARIANTS_UNRESOLVED")
        has_external_center = bool(external_metal_neighbours.get(metal_index))
        local_external_bridge_count = sum(
            relation.donor_group in local_group_labels
            for relation in external_donor_relations
        )
        has_external_bridge = local_external_bridge_count > 0
        if has_external_center:
            site_issues.add("OUT_OF_DOMAIN_METAL_COCENTER")
        if has_external_bridge:
            site_issues.add("DONOR_SHARED_WITH_OUT_OF_DOMAIN_METAL")
        if (
            not local_group_records
            and (
                direct_metal_neighbours.get(metal_index)
                or external_metal_neighbours.get(metal_index)
            )
        ):
            record_level = "topology"
            resolution = "ambiguous" if has_disorder else "partial"
            if direct_metal_neighbours.get(metal_index):
                site_issues.add("DIRECT_METAL_RELATION_ONLY")
            if external_metal_neighbours.get(metal_index):
                site_issues.add("EXTERNAL_METAL_RELATION_ONLY")
            geometry_kind = "metal-graph-topology"
        elif not local_group_records:
            record_level = "audit_only"
            resolution = "unresolved"
            site_issues.add("NO_NATIVE_DONOR_EDGE")
        elif any(group.kind == "candidate_atom_set" for group in local_group_records):
            record_level = "topology"
            resolution = "ambiguous" if has_disorder else "partial"
            site_issues.add("AMBIGUOUS_COLLECTIVE_PI_SITE")
        elif any(group.kind == "pi_fragment" and group.hapticity > 1 for group in local_group_records):
            record_level = "site_object"
            resolution = "ambiguous" if has_disorder else (
                "partial"
                if (is_polymeric and not periodic_translation_complete)
                or has_external_center
                or has_external_bridge
                else "resolved"
            )
        elif geometry["geometry_kind"] == "CShM":
            record_level = "local_state"
            resolution = "ambiguous" if has_disorder else (
                "partial"
                if (is_polymeric and not periodic_translation_complete)
                or has_external_center
                or has_external_bridge
                else "resolved"
            )
        else:
            record_level = "topology"
            resolution = "ambiguous" if has_disorder else "partial"

        local_sites.append(MetalSiteRecord(
            label=label,
            element=str(getattr(atoms[metal_index], "atomic_symbol", "?")),
            donor_groups=tuple(sorted(group.label for group in local_group_records)),
            cn_site=cn_site,
            cn_atom=cn_atom,
            eta_sum=eta_sum,
            metal_degree=len(direct_metal_neighbours.get(metal_index, set())),
            external_metal_degree=len(external_metal_neighbours.get(metal_index, set())),
            external_bridge_count=local_external_bridge_count,
            donor_composition=tuple(sorted(composition.items())),
            geometry_kind=geometry_kind,
            shape_best=geometry["shape_best"],
            shape_values=geometry["shape_values"],
            radial_profile=geometry["radial_profile"],
            angular_profile=geometry["angular_profile"],
            record_level=record_level,
            resolution=resolution,
            canonical_status=canonical_status,
            issues=tuple(sorted(site_issues)),
        ))

    donor_pair_relations: List[DonorPairRelation] = []
    for metal_index in ordered_metals:
        metal_label = metal_labels[metal_index]
        metal_coord = _coordinates(atoms[metal_index])
        local_pairs = sorted(
            (
                draft_labels[id(draft)],
                draft,
                _group_coordinate(draft, atoms, metal_index),
            )
            for draft in ordered_drafts
            if metal_index in draft.targets
        )
        for (label_a, _, point_a), (label_b, _, point_b) in combinations(local_pairs, 2):
            angle = None
            centroid_distance = None
            angle_class = "unavailable"
            if point_a is not None and point_b is not None:
                centroid_distance = round(float(np.linalg.norm(point_a - point_b)), 4)
            if metal_coord is not None and point_a is not None and point_b is not None:
                vector_a = point_a - metal_coord
                vector_b = point_b - metal_coord
                denominator = float(np.linalg.norm(vector_a) * np.linalg.norm(vector_b))
                if denominator > 1.0e-12:
                    cosine = max(-1.0, min(1.0, float(np.dot(vector_a, vector_b) / denominator)))
                    angle = round(float(np.degrees(np.arccos(cosine))), 3)
                    if angle >= 165.0:
                        angle_class = "trans"
                    elif 75.0 <= angle <= 105.0:
                        angle_class = "cis"
                    else:
                        angle_class = "other"
            donor_pair_relations.append(DonorPairRelation(
                metal=metal_label,
                donor_group_a=label_a,
                donor_group_b=label_b,
                angle_deg=angle,
                angle_class=angle_class,
                centroid_distance=centroid_distance,
            ))

    donor_group_relations: List[DonorGroupRelation] = []
    overlap_draft_pairs: List[Tuple[_GroupDraft, _GroupDraft]] = []
    for draft_a, draft_b in combinations(ordered_drafts, 2):
        shared_atoms = sorted(set(draft_a.atoms).intersection(draft_b.atoms))
        if not shared_atoms:
            continue
        label_a = draft_labels[id(draft_a)]
        label_b = draft_labels[id(draft_b)]
        if label_a > label_b:
            draft_a, draft_b = draft_b, draft_a
            label_a, label_b = label_b, label_a
        member_pairs = tuple(sorted(
            (
                donor_ordinals[id(draft_a)][atom_index],
                donor_ordinals[id(draft_b)][atom_index],
            )
            for atom_index in shared_atoms
        ))
        donor_group_relations.append(DonorGroupRelation(
            donor_group_a=label_a,
            donor_group_b=label_b,
            relation="overlaps_on_member",
            shared_member_ordinals=member_pairs,
        ))
        overlap_draft_pairs.append((draft_a, draft_b))

    relation_map: Dict[Tuple[str, str, str, str], Set[str]] = defaultdict(set)
    relation_geometry: Dict[
        Tuple[str, str, str, str],
        Tuple[Optional[float], Tuple[int, int, int], bool],
    ] = {}
    for a, b, bond_type in direct_mm:
        la, lb = sorted((metal_labels[a], metal_labels[b]))
        key = (la, lb, "direct_metal_bond", bond_type)
        relation_map[key]
        coord_a = _coordinates(atoms[a])
        coord_b = _coordinates(atoms[b])
        distance = None
        if coord_a is not None and coord_b is not None:
            distance = round(float(np.linalg.norm(coord_a - coord_b)), 4)
        image_delta = (0, 0, 0)
        image_delta_resolved = True
        if is_polymeric:
            forward_a, forward_b = (a, b) if metal_labels[a] == la else (b, a)
            resolved_delta = _relative_image_delta(atoms[forward_a], atoms[forward_b])
            if resolved_delta is not None:
                image_delta = resolved_delta
            else:
                image_delta_resolved = False
        relation_geometry[key] = (distance, image_delta, image_delta_resolved)
    for group in donor_groups:
        if group.bridge_degree < 2:
            continue
        for a, b in combinations(group.target_metals, 2):
            la, lb = sorted((a, b))
            relation_map[(la, lb, "shared_donor", "")].add(group.label)
    for draft_a, draft_b in overlap_draft_pairs:
        through = {
            draft_labels[id(draft_a)],
            draft_labels[id(draft_b)],
        }
        for target_a in draft_a.targets:
            for target_b in draft_b.targets:
                if target_a == target_b:
                    continue
                la, lb = sorted((metal_labels[target_a], metal_labels[target_b]))
                relation_map[(la, lb, "shared_donor_member", "")].update(through)
    relations = tuple(
        MetalRelation(
            metal_a=key[0],
            metal_b=key[1],
            relation=key[2],
            through_groups=tuple(sorted(groups)),
            bond_type=key[3],
            distance=relation_geometry.get(key, (None, (0, 0, 0), True))[0],
            image_delta=relation_geometry.get(key, (None, (0, 0, 0), True))[1],
            image_delta_resolved=relation_geometry.get(
                key, (None, (0, 0, 0), True)
            )[2],
        )
        for key, groups in sorted(relation_map.items())
    )

    scope = "mono_eta1"
    if len(ordered_metals) > 1 and (any_haptic or any_pi_candidate):
        scope = "multi_haptic"
    elif len(ordered_metals) > 1:
        scope = "multi_eta1"
    elif any_haptic or any_pi_candidate:
        scope = "mono_haptic"

    flags = []
    issues = []
    if len(ordered_metals) > 1:
        flags.append("multimetal")
    entry_metal_blocks = metal_blocks(
        str(getattr(atoms[index], "atomic_symbol", "?"))
        for index in ordered_metals
    )
    flags.extend(f"metal_block_{block}" for block in entry_metal_blocks)
    if len(entry_metal_blocks) > 1:
        flags.append("mixed_metal_blocks")
    if any_haptic:
        flags.append("haptic")
    if any_pi_candidate:
        flags.append("ambiguous_pi_candidate")
        issues.append("AMBIGUOUS_COLLECTIVE_PI_SITE")
    if any(
        incidence.donor_element == "H" for incidence in incidences
    ):
        flags.append("metal_hydrogen")
    if external_center_relations:
        flags.append("heterometal_co_center")
        issues.append("OUT_OF_DOMAIN_METAL_COCENTER")
    if external_donor_relations:
        flags.append("external_metal_shared_donor")
        issues.append("DONOR_SHARED_WITH_OUT_OF_DOMAIN_METAL")
    if has_disorder:
        flags.append("disorder")
        issues.append("DISORDER_VARIANTS_UNRESOLVED")
    if is_polymeric:
        flags.append("polymeric")
        if periodic_translation_complete:
            flags.append("translation_resolved_source_graph")
        else:
            issues.append("PBC_TRANSLATION_UNRESOLVED")
    if len(ordered_metals) > 12:
        flags.append("high_nuclearity")
    if high_cn:
        flags.append("high_cn")
    if not exact_label_free:
        flags.append("source_order_noncanonical")
        issues.append("NONCANONICAL_SOURCE_TIEBREAK")
    if not exact_metal_order:
        flags.append("metal_order_tie")
    if not exact_group_order:
        flags.append("donor_group_order_tie")
    if not exact_donor_order:
        flags.append("within_group_donor_order_tie")

    if all(site.record_level == "audit_only" for site in local_sites):
        entry_resolution = "EMITTED_AUDIT_ONLY"
    elif has_disorder:
        entry_resolution = "EMITTED_AMBIGUOUS"
    elif any(site.resolution != "resolved" for site in local_sites):
        entry_resolution = "EMITTED_PARTIAL"
    else:
        entry_resolution = "EMITTED_RESOLVED"

    record = EntryRecord(
        source_index=int(source_index),
        refcode=str(entry.identifier),
        scope=scope,
        scope_flags=tuple(sorted(flags)),
        expected_metal_sites=len(ordered_metals),
        metal_sites=tuple(local_sites),
        donor_groups=tuple(donor_groups),
        incidences=tuple(incidences),
        metal_relations=relations,
        external_center_relations=tuple(external_center_relations),
        external_donor_relations=tuple(external_donor_relations),
        donor_group_relations=tuple(sorted(
            donor_group_relations,
            key=lambda relation: (
                relation.donor_group_a,
                relation.donor_group_b,
                relation.relation,
            ),
        )),
        donor_pair_relations=tuple(donor_pair_relations),
        entry_resolution=entry_resolution,
        issue_codes=tuple(sorted(set(issues))),
        source_release=source_release,
        protocol_id=protocol_id,
    ).with_digest()
    validation_issues = validate_entry_record(record)
    if validation_issues:
        raise CSDRecordError("SCHEMA_VALIDATION:" + ",".join(validation_issues))
    return record
