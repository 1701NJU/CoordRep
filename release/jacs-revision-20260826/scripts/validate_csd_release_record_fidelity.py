#!/usr/bin/env python3
"""Independent source-topology validation for emitted CSD audit records.

This validator intentionally does not import the candidate CSD adapter.  It
reconstructs metal--donor incidences and site grouping directly from CCDC
objects and compares chemically meaningful multisets with the emitted JSON.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import sys
from collections import Counter, defaultdict, deque
from itertools import combinations
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Sequence, Set, Tuple

import numpy as np
import ccdc
import scipy
from ccdc.io import EntryReader


BASE = Path(__file__).resolve().parent.parent
if str(BASE) not in sys.path:
    sys.path.insert(0, str(BASE))

from coordrep.audit.metal_policy import (
    METAL_POLICY_ID,
    METAL_POLICY_SHA256,
    is_target_metal_symbol,
)


MAX_CONFIRMED_COLLECTIVE_SITE_MEMBERS = 12


def _is_pi(value: Any) -> bool:
    text = str(value).lower() if value is not None else ""
    return "pi" in text or "deloc" in text


def _is_source_metal(atom: Any) -> bool:
    return is_target_metal_symbol(getattr(atom, "atomic_symbol", "?"))


def _formal_charge(atom: Any) -> str:
    for name in ("formal_charge", "charge"):
        try:
            value = getattr(atom, name)
        except Exception:
            continue
        if value is not None:
            return str(value)
    return "unk"


def _short_hash(value: Any, length: int = 20) -> str:
    return hashlib.sha256(repr(value).encode("utf-8")).hexdigest()[:length]


def _wl_colours(
    atoms_by_index: Mapping[int, Any],
    adjacency: Mapping[int, Sequence[Tuple[int, str]]],
    rounds: int = 4,
) -> Dict[int, str]:
    colours = {
        index: _short_hash((
            str(getattr(atom, "atomic_symbol", "?")),
            _formal_charge(atom),
            len(adjacency.get(index, ())),
        ))
        for index, atom in atoms_by_index.items()
    }
    for _ in range(rounds):
        colours = {
            index: _short_hash((
                colours[index],
                tuple(sorted(
                    (bond_type, colours[other])
                    for other, bond_type in adjacency.get(index, ())
                )),
            ))
            for index in atoms_by_index
        }
    return colours


def _fractional(atom: Any):
    try:
        value = atom.fractional_coordinates
        array = np.asarray([value.x, value.y, value.z], dtype=float)
    except Exception:
        return None
    return array if array.shape == (3,) and np.all(np.isfinite(array)) else None


def _cartesian(atom: Any):
    try:
        value = atom.coordinates
        array = np.asarray([value.x, value.y, value.z], dtype=float)
    except Exception:
        return None
    return array if array.shape == (3,) and np.all(np.isfinite(array)) else None


def _angle_relation(
    metal: Any, members_a: Sequence[Any], members_b: Sequence[Any]
) -> Tuple[Any, str, Any]:
    metal_coord = _cartesian(metal)
    points_a = [_cartesian(atom) for atom in members_a]
    points_b = [_cartesian(atom) for atom in members_b]
    if (
        metal_coord is None
        or not points_a
        or not points_b
        or any(point is None for point in points_a + points_b)
    ):
        return None, "unavailable", None
    point_a = np.mean(np.asarray(points_a), axis=0)
    point_b = np.mean(np.asarray(points_b), axis=0)
    centroid_distance = round(float(np.linalg.norm(point_a - point_b)), 4)
    vector_a = point_a - metal_coord
    vector_b = point_b - metal_coord
    denominator = float(np.linalg.norm(vector_a) * np.linalg.norm(vector_b))
    if denominator <= 1.0e-12:
        return None, "unavailable", centroid_distance
    cosine = max(-1.0, min(1.0, float(np.dot(vector_a, vector_b) / denominator)))
    angle = round(float(np.degrees(np.arccos(cosine))), 3)
    angle_class = "trans" if angle >= 165.0 else (
        "cis" if 75.0 <= angle <= 105.0 else "other"
    )
    return angle, angle_class, centroid_distance


def _image_delta_fact(
    metal: Any, donor: Any, polymeric: bool
) -> Tuple[Tuple[int, int, int], bool]:
    if not polymeric:
        return (0, 0, 0), True
    fm = _fractional(metal)
    fd = _fractional(donor)
    if fm is None or fd is None:
        return (0, 0, 0), False
    return tuple(
        int(value)
        for value in np.floor(fd + 1.0e-8).astype(int)
        - np.floor(fm + 1.0e-8).astype(int)
    ), True


def _canonical_direct_signature(
    element_a: str,
    element_b: str,
    bond_type: str,
    distance: Any,
    image_delta: Sequence[int],
    image_delta_resolved: bool,
) -> Tuple[Any, ...]:
    delta = tuple(int(value) for value in image_delta)
    reversed_delta = tuple(-value for value in delta)
    oriented = min(
        (element_a, element_b, delta),
        (element_b, element_a, reversed_delta),
    )
    return oriented + (bool(image_delta_resolved), bond_type, distance)


def _geometry_profiles(
    metal: Any,
    groups: Sequence[Tuple[str, Tuple[int, ...], str]],
    by_index: Mapping[int, Any],
) -> Tuple[Tuple[float, ...], Tuple[float, ...]]:
    metal_coord = _cartesian(metal)
    points = []
    for _, members, _ in groups:
        member_points = [_cartesian(by_index[index]) for index in members]
        if not member_points or any(point is None for point in member_points):
            return (), ()
        points.append(np.mean(np.asarray(member_points), axis=0))
    if metal_coord is None or not points:
        return (), ()
    vectors = np.asarray(points) - metal_coord
    radii = np.linalg.norm(vectors, axis=1)
    radial = tuple(sorted(round(float(value), 4) for value in radii))
    angular = []
    for left in range(len(vectors)):
        for right in range(left + 1, len(vectors)):
            denominator = float(radii[left] * radii[right])
            if denominator <= 1.0e-12:
                continue
            cosine = float(np.dot(vectors[left], vectors[right]) / denominator)
            angular.append(round(max(-1.0, min(1.0, cosine)), 6))
    return radial, tuple(sorted(angular))


def _source_group_descriptor(
    group: Mapping[str, Any], by_index: Mapping[int, Any]
) -> Tuple[Any, ...]:
    return (
        group["kind"],
        tuple(sorted(by_index[index].atomic_symbol for index in group["members"])),
        len(group["members"]),
        tuple(sorted(by_index[index].atomic_symbol for index in group["targets"])),
    )


def _record_group_descriptor(
    group: Mapping[str, Any], metal_by_label: Mapping[str, str]
) -> Tuple[Any, ...]:
    return (
        group["kind"],
        tuple(sorted(group["donor_elements"])),
        len(group["donor_elements"]),
        tuple(sorted(metal_by_label[label] for label in group["target_metals"])),
    )


def _source_facts(entry: Any) -> Dict[str, Any]:
    molecule = entry.molecule
    atoms = list(molecule.atoms)
    by_index = {int(atom.index): atom for atom in atoms}
    if len(by_index) != len(atoms):
        raise ValueError("source atom.index is not unique")
    metals = [atom for atom in atoms if is_target_metal_symbol(atom.atomic_symbol)]
    metal_indices = {int(atom.index) for atom in metals}
    polymeric = bool(getattr(entry, "is_polymeric", False) or getattr(molecule, "is_polymeric", False))

    adjacency: Dict[int, Set[int]] = defaultdict(set)
    topology_adjacency: Dict[int, List[Tuple[int, str]]] = defaultdict(list)
    donor_edges: List[Tuple[int, int, str, str]] = []
    direct_mm: List[Tuple[int, int, str]] = []
    external_edges: List[Tuple[int, int, str]] = []
    external_donor_edge_types: Dict[Tuple[int, int], Set[str]] = defaultdict(set)
    donor_targets: Dict[int, Set[int]] = defaultdict(set)
    # Atomic coordination sites are keyed by source atom identity, irrespective
    # of whether CSD labels a particular metal--atom edge Single or
    # Pi/Delocalised.  A Pi/Delocalised donor is promoted to a collective
    # pi_fragment only when it belongs to a connected eta>1 set for that metal.
    atomic_targets: Dict[int, Set[int]] = defaultdict(set)
    pi_donors_by_metal: Dict[int, Set[int]] = defaultdict(set)

    seen_bonds: Set[Tuple[int, int, str]] = set()
    for bond in molecule.bonds:
        a, b = list(bond.atoms)
        ai, bi = int(a.index), int(b.index)
        raw_bond_text = str(bond.bond_type) if bond.bond_type is not None else "unknown"
        bond_text = " ".join(raw_bond_text.split()) or "unknown"
        bond_key = (min(ai, bi), max(ai, bi), bond_text)
        if bond_key in seen_bonds:
            continue
        seen_bonds.add(bond_key)
        adjacency[ai].add(bi)
        adjacency[bi].add(ai)
        topology_adjacency[ai].append((bi, bond_text))
        topology_adjacency[bi].append((ai, bond_text))
        a_metal, b_metal = ai in metal_indices, bi in metal_indices
        if a_metal and b_metal:
            direct_mm.append((min(ai, bi), max(ai, bi), bond_text))
            continue
        if not a_metal and not b_metal:
            a_source_metal = _is_source_metal(a)
            b_source_metal = _is_source_metal(b)
            if a_source_metal != b_source_metal:
                donor_index = bi if a_source_metal else ai
                partner_index = ai if a_source_metal else bi
                external_donor_edge_types[(donor_index, partner_index)].add(
                    bond_text
                )
            continue
        if a_metal == b_metal:
            continue
        metal_index = ai if a_metal else bi
        donor_index = bi if a_metal else ai
        donor = by_index[donor_index]
        if _is_source_metal(donor):
            external_edges.append((metal_index, donor_index, bond_text))
            continue
        role = "pi" if _is_pi(bond.bond_type) else "sigma"
        donor_edges.append((metal_index, donor_index, role, bond_text))
        donor_targets[donor_index].add(metal_index)
        if role == "pi":
            pi_donors_by_metal[metal_index].add(donor_index)

    topology_colours = _wl_colours(by_index, topology_adjacency)

    local_site_counts: Dict[int, int] = {}
    haptic_signatures: Counter = Counter()
    candidate_signatures: Counter = Counter()
    pi_component_targets: Dict[Tuple[int, ...], Set[int]] = defaultdict(set)
    local_group_members: Dict[int, List[Tuple[str, Tuple[int, ...], str]]] = {}
    for metal in metals:
        metal_index = int(metal.index)
        local_edges = [edge for edge in donor_edges if edge[0] == metal_index]
        sigma_atoms = {edge[1] for edge in local_edges if edge[2] == "sigma"}
        unseen = set(pi_donors_by_metal[metal_index])
        pi_components = []
        while unseen:
            start = min(unseen)
            unseen.remove(start)
            queue = deque([start])
            component = [start]
            while queue:
                current = queue.popleft()
                for other in adjacency[current]:
                    if other in unseen:
                        unseen.remove(other)
                        queue.append(other)
                        component.append(other)
            pi_components.append(component)
            if len(component) > 1:
                pi_component_targets[tuple(sorted(component))].add(metal_index)
                if len(component) <= MAX_CONFIRMED_COLLECTIVE_SITE_MEMBERS:
                    haptic_signatures[(metal.atomic_symbol, len(component))] += 1
                else:
                    candidate_signatures[(metal.atomic_symbol, len(component))] += 1
        singleton_pi_atoms = {
            component[0] for component in pi_components if len(component) == 1
        }
        atomic_atoms = sigma_atoms | singleton_pi_atoms
        for donor_index in atomic_atoms:
            atomic_targets[donor_index].add(metal_index)

        local_atomic_groups = []
        for donor_index in sorted(atomic_atoms):
            target_roles = {
                edge[2]
                for edge in local_edges
                if edge[1] == donor_index
            }
            local_atomic_groups.append((
                "atom",
                (donor_index,),
                "+".join(sorted(target_roles)),
            ))
        local_pi_groups = [
            (
                (
                    "pi_fragment"
                    if len(component) <= MAX_CONFIRMED_COLLECTIVE_SITE_MEMBERS
                    else "candidate_atom_set"
                ),
                tuple(component),
                "pi",
            )
            for component in pi_components
            if len(component) > 1
        ]
        local_group_members[metal_index] = local_atomic_groups + local_pi_groups
        local_site_counts[metal_index] = len(local_atomic_groups) + len(local_pi_groups)

    incidence_signature = Counter()
    incidence_edges: Dict[Tuple[int, int], Dict[str, Set[str]]] = defaultdict(
        lambda: {"roles": set(), "bond_types": set()}
    )
    for metal_index, donor_index, role, bond_text in donor_edges:
        incidence_edges[(metal_index, donor_index)]["roles"].add(role)
        incidence_edges[(metal_index, donor_index)]["bond_types"].add(bond_text)
    for (metal_index, donor_index), edge_data in incidence_edges.items():
        metal = by_index[metal_index]
        donor = by_index[donor_index]
        image_delta, image_delta_resolved = _image_delta_fact(
            metal, donor, polymeric
        )
        incidence_signature[(
            metal.atomic_symbol,
            donor.atomic_symbol,
            "+".join(sorted(edge_data["roles"])),
            "+".join(sorted(edge_data["bond_types"])),
            image_delta,
            image_delta_resolved,
            (
                round(float(np.linalg.norm(_cartesian(metal) - _cartesian(donor))), 4)
                if _cartesian(metal) is not None and _cartesian(donor) is not None
                else None
            ),
        )] += 1

    donor_pair_signature = Counter()
    geometry_profile_signature = Counter()
    for metal in metals:
        metal_index = int(metal.index)
        groups = local_group_members[metal_index]
        radial_profile, angular_profile = _geometry_profiles(
            metal, groups, by_index
        )
        geometry_profile_signature[(
            metal.atomic_symbol,
            radial_profile,
            angular_profile,
        )] += 1
        for left in range(len(groups)):
            for right in range(left + 1, len(groups)):
                kind_a, members_a, role_a = groups[left]
                kind_b, members_b, role_b = groups[right]
                descriptor_a = (
                    kind_a,
                    tuple(sorted(by_index[index].atomic_symbol for index in members_a)),
                    role_a,
                )
                descriptor_b = (
                    kind_b,
                    tuple(sorted(by_index[index].atomic_symbol for index in members_b)),
                    role_b,
                )
                descriptors = tuple(sorted((descriptor_a, descriptor_b), key=repr))
                angle, angle_class, separation = _angle_relation(
                    metal,
                    [by_index[index] for index in members_a],
                    [by_index[index] for index in members_b],
                )
                donor_pair_signature[(
                    metal.atomic_symbol,
                    descriptors,
                    angle,
                    angle_class,
                    separation,
                )] += 1

    bridge_signature = Counter()
    for donor_index, targets in atomic_targets.items():
        if len(targets) > 1:
            bridge_signature[(
                "atom",
                (by_index[donor_index].atomic_symbol,),
                1,
                len(targets),
            )] += 1
    for component, targets in pi_component_targets.items():
        if len(targets) > 1:
            elements = tuple(sorted(by_index[index].atomic_symbol for index in component))
            is_confirmed_collective = (
                len(component) <= MAX_CONFIRMED_COLLECTIVE_SITE_MEMBERS
            )
            bridge_signature[(
                "pi_fragment" if is_confirmed_collective else "candidate_atom_set",
                elements,
                len(component) if is_confirmed_collective else 1,
                len(targets),
            )] += 1

    source_groups: Dict[Tuple[str, Tuple[int, ...]], Dict[str, Any]] = {}
    for donor_index, targets in atomic_targets.items():
        key = ("atom", (donor_index,))
        source_groups[key] = {
            "kind": "atom",
            "members": (donor_index,),
            "targets": set(targets),
        }
    for component, targets in pi_component_targets.items():
        kind = (
            "pi_fragment"
            if len(component) <= MAX_CONFIRMED_COLLECTIVE_SITE_MEMBERS
            else "candidate_atom_set"
        )
        key = (kind, tuple(component))
        source_groups[key] = {
            "kind": kind,
            "members": tuple(component),
            "targets": set(targets),
        }

    external_donor_signature = Counter()
    external_group_contact_signature = Counter()
    external_contacts_by_group: Dict[
        Tuple[str, Tuple[int, ...]], int
    ] = defaultdict(int)
    used_external_donor_edges: Set[Tuple[int, int]] = set()
    for group_key, group in source_groups.items():
        group_descriptor = _source_group_descriptor(group, by_index)
        for donor_index in group["members"]:
            for (candidate_donor, partner_index), bond_types in external_donor_edge_types.items():
                if candidate_donor != donor_index:
                    continue
                donor = by_index[donor_index]
                partner = by_index[partner_index]
                donor_coord, partner_coord = _cartesian(donor), _cartesian(partner)
                distance = (
                    round(float(np.linalg.norm(partner_coord - donor_coord)), 4)
                    if donor_coord is not None and partner_coord is not None
                    else None
                )
                image_delta, image_delta_resolved = _image_delta_fact(
                    donor, partner, polymeric
                )
                external_donor_signature[(
                    group_descriptor,
                    donor.atomic_symbol,
                    topology_colours[donor_index],
                    partner.atomic_symbol,
                    topology_colours[partner_index],
                    "+".join(sorted(bond_types)),
                    distance,
                    image_delta,
                    image_delta_resolved,
                )] += 1
                external_contacts_by_group[group_key] += 1
                used_external_donor_edges.add((donor_index, partner_index))
        external_group_contact_signature[(
            group_descriptor,
            external_contacts_by_group[group_key],
        )] += 1

    external_bridge_site_signature = Counter()
    for metal in metals:
        metal_index = int(metal.index)
        relation_count = sum(
            count
            for group_key, count in external_contacts_by_group.items()
            if metal_index in source_groups[group_key]["targets"]
        )
        external_bridge_site_signature[(metal.atomic_symbol, relation_count)] += 1

    shared_donor_metal_map: Dict[
        Tuple[int, int], Set[Tuple[str, Tuple[int, ...]]]
    ] = defaultdict(set)
    for group_key, group in source_groups.items():
        for metal_a, metal_b in combinations(sorted(group["targets"]), 2):
            shared_donor_metal_map[(metal_a, metal_b)].add(group_key)
    shared_donor_metal_signature = Counter()
    for (metal_a, metal_b), through_keys in shared_donor_metal_map.items():
        metal_elements = tuple(sorted((
            by_index[metal_a].atomic_symbol,
            by_index[metal_b].atomic_symbol,
        )))
        through_descriptors = tuple(sorted(
            (
                _source_group_descriptor(source_groups[key], by_index)
                for key in through_keys
            ),
            key=repr,
        ))
        shared_donor_metal_signature[(metal_elements, through_descriptors)] += 1

    overlap_signature = Counter()
    shared_member_metal_map: Dict[
        Tuple[int, int], Set[Tuple[str, Tuple[int, ...]]]
    ] = defaultdict(set)
    source_group_items = list(source_groups.items())
    for left in range(len(source_group_items)):
        key_a, group_a = source_group_items[left]
        for right in range(left + 1, len(source_group_items)):
            key_b, group_b = source_group_items[right]
            shared_members = set(group_a["members"]).intersection(group_b["members"])
            if not shared_members:
                continue
            descriptor_a = _source_group_descriptor(group_a, by_index)
            descriptor_b = _source_group_descriptor(group_b, by_index)
            descriptors = tuple(sorted((descriptor_a, descriptor_b), key=repr))
            overlap_signature[(
                descriptors,
                tuple(sorted(by_index[index].atomic_symbol for index in shared_members)),
            )] += 1
            for target_a in group_a["targets"]:
                for target_b in group_b["targets"]:
                    if target_a == target_b:
                        continue
                    metal_pair = tuple(sorted((target_a, target_b)))
                    shared_member_metal_map[metal_pair].update((key_a, key_b))

    shared_member_metal_signature = Counter()
    for (metal_a, metal_b), through_keys in shared_member_metal_map.items():
        metal_elements = tuple(sorted((
            by_index[metal_a].atomic_symbol,
            by_index[metal_b].atomic_symbol,
        )))
        through_descriptors = tuple(sorted(
            (
                _source_group_descriptor(source_groups[key], by_index)
                for key in through_keys
            ),
            key=repr,
        ))
        shared_member_metal_signature[(metal_elements, through_descriptors)] += 1

    direct_signature = Counter()
    for a, b, bond_type in direct_mm:
        atom_a, atom_b = by_index[a], by_index[b]
        coord_a, coord_b = _cartesian(atom_a), _cartesian(atom_b)
        distance = (
            round(float(np.linalg.norm(coord_a - coord_b)), 4)
            if coord_a is not None and coord_b is not None
            else None
        )
        image_delta, image_delta_resolved = _image_delta_fact(
            atom_a, atom_b, polymeric
        )
        direct_signature[_canonical_direct_signature(
            atom_a.atomic_symbol,
            atom_b.atomic_symbol,
            bond_type,
            distance,
            image_delta,
            image_delta_resolved,
        )] += 1
    direct_degree: Counter = Counter()
    for a, b, _ in direct_mm:
        direct_degree[a] += 1
        direct_degree[b] += 1

    external_signature = Counter()
    external_neighbours: Dict[int, Set[int]] = defaultdict(set)
    for metal_index, partner_index, bond_type in external_edges:
        metal = by_index[metal_index]
        partner = by_index[partner_index]
        metal_coord, partner_coord = _cartesian(metal), _cartesian(partner)
        distance = (
            round(float(np.linalg.norm(partner_coord - metal_coord)), 4)
            if metal_coord is not None and partner_coord is not None
            else None
        )
        image_delta, image_delta_resolved = _image_delta_fact(
            metal, partner, polymeric
        )
        external_signature[(
            metal.atomic_symbol,
            partner.atomic_symbol,
            topology_colours[partner_index],
            bond_type,
            distance,
            image_delta,
            image_delta_resolved,
        )] += 1
        external_neighbours[metal_index].add(partner_index)
    translation_resolved_flag = bool(
        polymeric
        and (incidence_edges or direct_mm or external_edges)
        and all(
            _image_delta_fact(by_index[metal_index], by_index[donor_index], True)[1]
            for metal_index, donor_index in incidence_edges
        )
        and all(
            _image_delta_fact(by_index[metal_a], by_index[metal_b], True)[1]
            for metal_a, metal_b, _ in direct_mm
        )
        and all(
            _image_delta_fact(by_index[metal_index], by_index[partner_index], True)[1]
            for metal_index, partner_index, _ in external_edges
        )
        and all(
            _image_delta_fact(by_index[donor_index], by_index[partner_index], True)[1]
            for donor_index, partner_index in used_external_donor_edges
        )
    )
    return {
        "metal_elements": Counter(atom.atomic_symbol for atom in metals),
        "metal_count": len(metals),
        "incidence_signature": incidence_signature,
        "incidence_count": len(incidence_edges),
        "local_cn_atom": Counter(
            (
                by_index[index].atomic_symbol,
                sum(metal_index == index for metal_index, _ in incidence_edges),
            )
            for index in metal_indices
        ),
        "local_cn_site": Counter(
            (by_index[index].atomic_symbol, local_site_counts[index]) for index in metal_indices
        ),
        "haptic_signature": haptic_signatures,
        "candidate_signature": candidate_signatures,
        "bridge_signature": bridge_signature,
        "shared_donor_metal_signature": shared_donor_metal_signature,
        "overlap_signature": overlap_signature,
        "shared_member_metal_signature": shared_member_metal_signature,
        "direct_signature": direct_signature,
        "metal_degree_signature": Counter(
            (by_index[index].atomic_symbol, direct_degree[index]) for index in metal_indices
        ),
        "external_signature": external_signature,
        "external_degree_signature": Counter(
            (by_index[index].atomic_symbol, len(external_neighbours[index]))
            for index in metal_indices
        ),
        "external_donor_signature": external_donor_signature,
        "external_group_contact_signature": external_group_contact_signature,
        "external_bridge_site_signature": external_bridge_site_signature,
        "local_eta_sum": Counter(
            (
                by_index[index].atomic_symbol,
                (
                    None
                    if any(
                        kind == "candidate_atom_set"
                        for kind, _, _ in local_group_members[index]
                    )
                    else sum(
                        len(members) if kind == "pi_fragment" else 1
                        for kind, members, _ in local_group_members[index]
                    )
                ),
            )
            for index in metal_indices
        ),
        "donor_pair_signature": donor_pair_signature,
        "geometry_profile_signature": geometry_profile_signature,
        "polymeric": polymeric,
        "translation_resolved_flag": translation_resolved_flag,
    }


def _record_facts(record: Mapping[str, Any]) -> Dict[str, Any]:
    metal_by_label = {site["label"]: site["element"] for site in record["metal_sites"]}
    group_by_label = {group["label"]: group for group in record["donor_groups"]}
    incidence_signature = Counter()
    for incidence in record["incidences"]:
        incidence_signature[(
            metal_by_label[incidence["metal"]],
            incidence["donor_element"],
            incidence["bond_role"],
            incidence["bond_type"],
            tuple(incidence["image_delta"]),
            bool(incidence.get("image_delta_resolved", True)),
            incidence.get("distance"),
        )] += 1
    haptic_signature = Counter()
    candidate_signature = Counter()
    bridge_signature = Counter()
    for group in record["donor_groups"]:
        if group["kind"] == "pi_fragment" and int(group["hapticity"]) > 1:
            for target in group["target_metals"]:
                haptic_signature[(metal_by_label[target], int(group["hapticity"]))] += 1
        if group["kind"] == "candidate_atom_set":
            for target in group["target_metals"]:
                candidate_signature[(
                    metal_by_label[target],
                    len(group["donor_elements"]),
                )] += 1
        if int(group["bridge_degree"]) > 1:
            bridge_signature[(
                group["kind"],
                tuple(sorted(group["donor_elements"])),
                int(group["hapticity"]),
                int(group["bridge_degree"]),
            )] += 1

    overlap_signature = Counter()
    for relation in record.get("donor_group_relations", []):
        group_a = group_by_label[relation["donor_group_a"]]
        group_b = group_by_label[relation["donor_group_b"]]
        descriptor_a = _record_group_descriptor(group_a, metal_by_label)
        descriptor_b = _record_group_descriptor(group_b, metal_by_label)
        descriptors = tuple(sorted((descriptor_a, descriptor_b), key=repr))
        shared_elements = []
        for ordinal_a, ordinal_b in relation.get("shared_member_ordinals", []):
            element_a = group_a["donor_elements"][int(ordinal_a)]
            element_b = group_b["donor_elements"][int(ordinal_b)]
            shared_elements.append(
                element_a if element_a == element_b else f"{element_a}!={element_b}"
            )
        overlap_signature[(descriptors, tuple(sorted(shared_elements)))] += 1

    shared_member_metal_signature = Counter()
    for relation in record["metal_relations"]:
        if relation["relation"] != "shared_donor_member":
            continue
        metal_elements = tuple(sorted((
            metal_by_label[relation["metal_a"]],
            metal_by_label[relation["metal_b"]],
        )))
        through_descriptors = tuple(sorted(
            (
                _record_group_descriptor(group_by_label[label], metal_by_label)
                for label in relation.get("through_groups", [])
            ),
            key=repr,
        ))
        shared_member_metal_signature[(metal_elements, through_descriptors)] += 1

    shared_donor_metal_signature = Counter()
    for relation in record["metal_relations"]:
        if relation["relation"] != "shared_donor":
            continue
        metal_elements = tuple(sorted((
            metal_by_label[relation["metal_a"]],
            metal_by_label[relation["metal_b"]],
        )))
        through_descriptors = tuple(sorted(
            (
                _record_group_descriptor(group_by_label[label], metal_by_label)
                for label in relation.get("through_groups", [])
            ),
            key=repr,
        ))
        shared_donor_metal_signature[(metal_elements, through_descriptors)] += 1

    direct_signature = Counter()
    for relation in record["metal_relations"]:
        if relation["relation"] != "direct_metal_bond":
            continue
        direct_signature[_canonical_direct_signature(
            metal_by_label[relation["metal_a"]],
            metal_by_label[relation["metal_b"]],
            relation["bond_type"],
            relation.get("distance"),
            relation.get("image_delta", (0, 0, 0)),
            bool(relation.get("image_delta_resolved", True)),
        )] += 1

    external_signature = Counter()
    for relation in record.get("external_center_relations", []):
        external_signature[(
            metal_by_label[relation["metal"]],
            relation["partner_element"],
            relation["partner_topology_key"],
            relation["bond_type"],
            relation.get("distance"),
            tuple(relation.get("image_delta", (0, 0, 0))),
            bool(relation.get("image_delta_resolved", True)),
        )] += 1
    external_donor_signature = Counter()
    for relation in record.get("external_donor_relations", []):
        group = group_by_label[relation["donor_group"]]
        donor_ordinal = int(relation["donor_ordinal"])
        external_donor_signature[(
            _record_group_descriptor(group, metal_by_label),
            group["donor_elements"][donor_ordinal],
            group["topology_keys"][donor_ordinal],
            relation["partner_element"],
            relation["partner_topology_key"],
            relation["bond_type"],
            relation.get("distance"),
            tuple(relation.get("image_delta", (0, 0, 0))),
            bool(relation.get("image_delta_resolved", True)),
        )] += 1
    external_group_contact_signature = Counter(
        (
            _record_group_descriptor(group, metal_by_label),
            int(group.get("external_donor_contact_count", 0)),
        )
        for group in record["donor_groups"]
    )
    group_roles = defaultdict(set)
    for incidence in record["incidences"]:
        group_roles[(
            incidence["donor_group"],
            incidence["metal"],
        )].add(incidence["bond_role"])
    donor_pair_signature = Counter()
    for relation in record.get("donor_pair_relations", []):
        group_a = group_by_label[relation["donor_group_a"]]
        group_b = group_by_label[relation["donor_group_b"]]
        descriptor_a = (
            group_a["kind"],
            tuple(sorted(group_a["donor_elements"])),
            "+".join(sorted(group_roles[(group_a["label"], relation["metal"])])),
        )
        descriptor_b = (
            group_b["kind"],
            tuple(sorted(group_b["donor_elements"])),
            "+".join(sorted(group_roles[(group_b["label"], relation["metal"])])),
        )
        descriptors = tuple(sorted((descriptor_a, descriptor_b), key=repr))
        donor_pair_signature[(
            metal_by_label[relation["metal"]],
            descriptors,
            relation.get("angle_deg"),
            relation["angle_class"],
            relation.get("centroid_distance"),
        )] += 1
    return {
        "metal_elements": Counter(metal_by_label.values()),
        "metal_count": len(metal_by_label),
        "incidence_signature": incidence_signature,
        "incidence_count": len(record["incidences"]),
        "local_cn_atom": Counter((site["element"], int(site["cn_atom"])) for site in record["metal_sites"]),
        "local_cn_site": Counter((site["element"], int(site["cn_site"])) for site in record["metal_sites"]),
        "haptic_signature": haptic_signature,
        "candidate_signature": candidate_signature,
        "bridge_signature": bridge_signature,
        "shared_donor_metal_signature": shared_donor_metal_signature,
        "overlap_signature": overlap_signature,
        "shared_member_metal_signature": shared_member_metal_signature,
        "direct_signature": direct_signature,
        "metal_degree_signature": Counter(
            (site["element"], int(site.get("metal_degree", 0)))
            for site in record["metal_sites"]
        ),
        "external_signature": external_signature,
        "external_degree_signature": Counter(
            (site["element"], int(site.get("external_metal_degree", 0)))
            for site in record["metal_sites"]
        ),
        "external_donor_signature": external_donor_signature,
        "external_group_contact_signature": external_group_contact_signature,
        "external_bridge_site_signature": Counter(
            (site["element"], int(site.get("external_bridge_count", 0)))
            for site in record["metal_sites"]
        ),
        "local_eta_sum": Counter(
            (site["element"], site.get("eta_sum"))
            for site in record["metal_sites"]
        ),
        "donor_pair_signature": donor_pair_signature,
        "geometry_profile_signature": Counter(
            (
                site["element"],
                tuple(site.get("radial_profile", ())),
                tuple(site.get("angular_profile", ())),
            )
            for site in record["metal_sites"]
        ),
        "polymeric": "polymeric" in record["scope_flags"],
        "translation_resolved_flag": (
            "translation_resolved_source_graph" in record["scope_flags"]
        ),
    }


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--database", default="CSD")
    parser.add_argument("--records", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    reader = EntryReader(args.database)
    records = [json.loads(line) for line in args.records.open(encoding="utf-8") if line.strip()]
    comparison_fields = (
        "metal_elements",
        "metal_count",
        "incidence_signature",
        "incidence_count",
        "local_cn_atom",
        "local_cn_site",
        "haptic_signature",
        "candidate_signature",
        "bridge_signature",
        "shared_donor_metal_signature",
        "overlap_signature",
        "shared_member_metal_signature",
        "direct_signature",
        "metal_degree_signature",
        "external_signature",
        "external_degree_signature",
        "external_donor_signature",
        "external_group_contact_signature",
        "external_bridge_site_signature",
        "local_eta_sum",
        "donor_pair_signature",
        "geometry_profile_signature",
        "polymeric",
        "translation_resolved_flag",
    )
    rows = []
    failures = Counter()
    for position, record in enumerate(records, 1):
        refcode = record["refcode"]
        issues: List[str] = []
        try:
            source = _source_facts(reader.entry(refcode))
            emitted = _record_facts(record)
            for field in comparison_fields:
                if source[field] != emitted[field]:
                    issues.append(field)
                    failures[field] += 1
        except Exception as exc:
            issues.append("validator_exception:" + type(exc).__name__)
            failures["validator_exception"] += 1
        rows.append({
            "refcode": refcode,
            "pass": not issues,
            "mismatched_fields": issues,
            "record_checksum_sha256": record.get("record_checksum_sha256", ""),
        })
        if position % 200 == 0:
            print(json.dumps({"validated": position, "passed": sum(row["pass"] for row in rows)}), flush=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="\n") as handle:
        for row in rows:
            handle.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")
    passed = sum(row["pass"] for row in rows)
    summary = {
        "records": len(rows),
        "passed": passed,
        "failed": len(rows) - passed,
        "pass_rate": passed / len(rows) if rows else None,
        "failure_fields": dict(failures),
        "records_sha256": hashlib.sha256(args.records.read_bytes()).hexdigest(),
        "validation_sha256": hashlib.sha256(args.output.read_bytes()).hexdigest(),
        "validator_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "metal_policy_id": METAL_POLICY_ID,
        "metal_policy_sha256": METAL_POLICY_SHA256,
        "database": {
            "argument": args.database,
            "entry_count": len(reader),
        },
        "runtime": {
            "python": sys.version.replace("\n", " "),
            "implementation": platform.python_implementation(),
            "platform": platform.platform(),
            "ccdc": str(getattr(ccdc, "__version__", "unknown")),
            "numpy": str(np.__version__),
            "scipy": str(scipy.__version__),
        },
    }
    summary_path = args.output.with_suffix(args.output.suffix + ".summary.json")
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
