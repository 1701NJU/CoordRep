"""Conservative distance-completion layer for CoordRep CSD audits.

This experimental v2 module is deliberately separate from the native-bond
serializer.  It never edits, deletes, or retypes a CSD-native edge.  Instead,
it records every non-native metal--neighbour observation as an accepted,
candidate, or rejected *distance decision* with an explicit reason.  Only the
accepted decisions may be used to construct an auxiliary, distance-augmented
first sphere; they do not change the v1 structural-record coverage statistic.

The default radii are the Cordero et al. single-bond covalent radii (angstrom):
B. Cordero et al., Dalton Trans. 2008, 2832--2838,
DOI 10.1039/B801115J.  Values are frozen here so that a release does not depend
on a mutable chemistry-toolkit table.
"""

from __future__ import annotations

import hashlib
import json
import math
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple


PROTOCOL_ID = "CoordRep-DistanceFirstSphere/2.0.0-20260825"

CONTACT_STATUSES = {"native", "accepted", "candidate", "rejected"}
METAL_BLOCKS = {"s", "p", "d", "f"}
DISORDER_STATES = {"clear", "compatible", "ambiguous", "incompatible"}

EVIDENCE_NATIVE = "N0_CSD_NATIVE_EDGE"
EVIDENCE_DISTANCE_STRONG = "D3_DISTANCE_STRONG"
EVIDENCE_DISTANCE_SHELL = "D2_DISTANCE_AND_SHELL"
EVIDENCE_DISTANCE_CANDIDATE = "D1_DISTANCE_CANDIDATE"
EVIDENCE_PROTOCOL_REJECTED = "D0_PROTOCOL_REJECTED"

# Atomic-number order.  Cordero values are available through Cm (Z = 96).
# Later values are intentionally None rather than silently substituted.
_ELEMENT_SYMBOLS: Tuple[str, ...] = (
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
    "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
    "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
)

_CORDER0_RADII: Tuple[Optional[float], ...] = (
    0.31, 0.28, 1.28, 0.96, 0.84, 0.76, 0.71, 0.66, 0.57, 0.58,
    1.66, 1.41, 1.21, 1.11, 1.07, 1.05, 1.02, 1.06, 2.03, 1.76,
    1.70, 1.60, 1.53, 1.39, 1.39, 1.32, 1.26, 1.24, 1.32, 1.22,
    1.22, 1.20, 1.19, 1.20, 1.20, 1.16, 2.20, 1.95, 1.90, 1.75,
    1.64, 1.54, 1.47, 1.46, 1.42, 1.39, 1.45, 1.44, 1.42, 1.39,
    1.39, 1.38, 1.39, 1.40, 2.44, 2.15, 2.07, 2.04, 2.03, 2.01,
    1.99, 1.98, 1.98, 1.96, 1.94, 1.92, 1.92, 1.89, 1.90, 1.87,
    1.87, 1.75, 1.70, 1.62, 1.51, 1.44, 1.41, 1.36, 1.36, 1.32,
    1.45, 1.46, 1.48, 1.40, 1.50, 1.50, 2.60, 2.21, 2.15, 2.06,
    2.00, 1.96, 1.90, 1.87, 1.80, 1.69,
    None, None, None, None, None, None, None, None, None, None, None,
    None, None, None, None, None, None, None, None, None, None, None,
)

if len(_ELEMENT_SYMBOLS) != len(_CORDER0_RADII):  # pragma: no cover
    raise RuntimeError("Frozen Cordero radius table is misaligned")

COVALENT_RADII_ANGSTROM: Mapping[str, Optional[float]] = dict(
    zip(_ELEMENT_SYMBOLS, _CORDER0_RADII)
)


@dataclass(frozen=True)
class DistancePolicy:
    """Frozen decision policy for non-native first-sphere observations."""

    policy_id: str = PROTOCOL_ID
    radii_source: str = "Cordero2008_DOI:10.1039/B801115J"
    strict_scale_s: float = 1.18
    strict_scale_p: float = 1.15
    strict_scale_d: float = 1.15
    strict_scale_f: float = 1.20
    candidate_scale_s: float = 1.45
    candidate_scale_p: float = 1.30
    candidate_scale_d: float = 1.30
    candidate_scale_f: float = 1.45
    ionic_shell_scale_s: float = 1.35
    ionic_shell_scale_p: float = 1.15
    ionic_shell_scale_d: float = 1.15
    ionic_shell_scale_f: float = 1.35
    shell_gap_absolute_A: float = 0.35
    shell_gap_relative: float = 0.15
    minimum_distance_A: float = 0.50
    search_cap_A: float = 4.50
    full_occupancy_floor: float = 0.98
    distance_round_digits: int = 4
    auto_accept_donor_elements: Tuple[str, ...] = (
        "N", "O", "F", "P", "S", "Cl", "Se", "Br", "Te", "I"
    )
    ionic_extension_donor_elements: Tuple[str, ...] = (
        "N", "O", "F", "S", "Cl", "Se", "Br", "Te", "I"
    )
    ionic_extension_blocks: Tuple[str, ...] = ("s", "f")

    def strict_scale(self, block: str) -> float:
        return float(getattr(self, f"strict_scale_{block}"))

    def candidate_scale(self, block: str) -> float:
        return float(getattr(self, f"candidate_scale_{block}"))

    def ionic_shell_scale(self, block: str) -> float:
        return float(getattr(self, f"ionic_shell_scale_{block}"))

    def canonical_payload(self) -> Dict[str, object]:
        payload = asdict(self)
        payload["covalent_radii_A"] = {
            key: COVALENT_RADII_ANGSTROM[key] for key in sorted(COVALENT_RADII_ANGSTROM)
        }
        return payload

    def digest(self) -> str:
        encoded = json.dumps(
            self.canonical_payload(), sort_keys=True, separators=(",", ":"), allow_nan=False
        )
        return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


DEFAULT_POLICY = DistancePolicy()


@dataclass(frozen=True)
class GeometryAtom:
    """Minimal atom view used by the adapter-independent neighbour enumerator.

    Coordinates are Cartesian for molecular enumeration and fractional for
    periodic enumeration.  Metal membership and block must be supplied by the
    caller's frozen metal policy rather than re-inferred from the element name.
    """

    key: str
    element: str
    coordinates: Tuple[float, float, float]
    is_metal: bool = False
    metal_block: str = ""
    orbit_key: str = ""
    occupancy: Optional[float] = 1.0
    disorder_assembly: str = ""
    disorder_group: str = ""
    ligand_component_key: str = ""
    pi_component_key: str = ""

    def __post_init__(self) -> None:
        if len(self.coordinates) != 3:
            raise ValueError("coordinates must contain three numbers")
        if self.is_metal and self.metal_block not in METAL_BLOCKS:
            raise ValueError("a target metal atom requires an s/p/d/f metal_block")


@dataclass(frozen=True)
class NativeEdgeKey:
    """Oriented source edge from a target metal to one donor physical image."""

    metal_key: str
    donor_key: str
    image_delta: Tuple[int, int, int] = (0, 0, 0)


@dataclass(frozen=True)
class ContactObservation:
    """One metal--neighbour observation in a specific physical image.

    ``donor_orbit_key`` identifies the wrapped crystallographic atom orbit.
    ``image_delta`` identifies the physical periodic image.  The same orbit at
    two translations is therefore represented by two observations.
    """

    metal_key: str
    metal_element: str
    metal_block: str
    donor_key: str
    donor_element: str
    distance_A: Optional[float]
    native_edge: bool = False
    donor_is_metal: bool = False
    periodic: bool = False
    donor_orbit_key: str = ""
    image_delta: Optional[Tuple[int, int, int]] = (0, 0, 0)
    disorder_state: str = "clear"
    metal_occupancy: Optional[float] = 1.0
    donor_occupancy: Optional[float] = 1.0
    ligand_component_key: str = ""
    pi_component_key: str = ""

    def __post_init__(self) -> None:
        if self.metal_block not in METAL_BLOCKS:
            raise ValueError(f"metal_block must be one of {sorted(METAL_BLOCKS)}")
        if self.disorder_state not in DISORDER_STATES:
            raise ValueError(f"unsupported disorder_state: {self.disorder_state}")
        if self.image_delta is not None and len(self.image_delta) != 3:
            raise ValueError("image_delta must contain three integers or be None")

    @property
    def physical_key(self) -> Tuple[object, ...]:
        orbit = self.donor_orbit_key or self.donor_key
        image = tuple(self.image_delta) if self.image_delta is not None else ("?", "?", "?")
        return (self.metal_key, orbit, image)


def _euclidean(vector: Sequence[float]) -> float:
    return math.sqrt(sum(float(value) ** 2 for value in vector))


def _pair_disorder_state(metal: GeometryAtom, donor: GeometryAtom) -> str:
    metal_tag = (metal.disorder_assembly, metal.disorder_group)
    donor_tag = (donor.disorder_assembly, donor.disorder_group)
    if metal_tag == ("", "") and donor_tag == ("", ""):
        return "clear"
    if metal_tag == donor_tag and metal_tag != ("", ""):
        return "compatible"
    if (
        metal.disorder_assembly
        and metal.disorder_assembly == donor.disorder_assembly
        and metal.disorder_group
        and donor.disorder_group
        and metal.disorder_group != donor.disorder_group
    ):
        return "incompatible"
    return "ambiguous"


def _make_observation(
    metal: GeometryAtom,
    donor: GeometryAtom,
    distance_A: Optional[float],
    *,
    native_edge: bool,
    periodic: bool,
    image_delta: Optional[Tuple[int, int, int]],
) -> ContactObservation:
    return ContactObservation(
        metal_key=metal.key,
        metal_element=metal.element,
        metal_block=metal.metal_block,
        donor_key=donor.key,
        donor_element=donor.element,
        distance_A=distance_A,
        native_edge=native_edge,
        donor_is_metal=donor.is_metal,
        periodic=periodic,
        donor_orbit_key=donor.orbit_key or donor.key,
        image_delta=image_delta,
        disorder_state=_pair_disorder_state(metal, donor),
        metal_occupancy=metal.occupancy,
        donor_occupancy=donor.occupancy,
        ligand_component_key=donor.ligand_component_key,
        pi_component_key=donor.pi_component_key,
    )


def enumerate_molecular_observations(
    atoms: Iterable[GeometryAtom],
    native_edges: Iterable[NativeEdgeKey] = (),
    policy: DistancePolicy = DEFAULT_POLICY,
) -> Tuple[ContactObservation, ...]:
    """Enumerate molecular metal--atom pairs without altering native edges."""

    atom_list = tuple(atoms)
    atom_by_key = {atom.key: atom for atom in atom_list}
    if len(atom_by_key) != len(atom_list):
        raise ValueError("GeometryAtom keys must be unique")
    native_set = set(native_edges)
    invalid = [edge for edge in native_set if edge.image_delta != (0, 0, 0)]
    if invalid:
        raise ValueError("molecular native edges must use image_delta=(0,0,0)")

    output: List[ContactObservation] = []
    seen_native = set()
    for metal in atom_list:
        if not metal.is_metal:
            continue
        for donor in atom_list:
            if donor.key == metal.key:
                continue
            edge = NativeEdgeKey(metal.key, donor.key, (0, 0, 0))
            vector = tuple(
                float(donor.coordinates[index]) - float(metal.coordinates[index])
                for index in range(3)
            )
            distance = _euclidean(vector)
            if edge in native_set:
                seen_native.add(edge)
            if distance <= policy.search_cap_A or edge in native_set:
                output.append(
                    _make_observation(
                        metal,
                        donor,
                        distance,
                        native_edge=edge in native_set,
                        periodic=False,
                        image_delta=(0, 0, 0),
                    )
                )

    missing = native_set - seen_native
    if missing:
        unresolved = sorted((edge.metal_key, edge.donor_key) for edge in missing)
        raise ValueError(f"native molecular edges could not be enumerated: {unresolved}")
    return tuple(output)


def _determinant3(matrix: Sequence[Sequence[float]]) -> float:
    a, b, c = matrix
    return (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    )


def _reciprocal_rows(cell: Sequence[Sequence[float]]) -> Tuple[Tuple[float, ...], ...]:
    """Return rows of A^-1 for row-vector cell convention."""

    a, b, c = tuple(tuple(float(value) for value in row) for row in cell)
    determinant = _determinant3((a, b, c))
    if abs(determinant) < 1e-12:
        raise ValueError("periodic cell is singular")
    return (
        (
            (b[1] * c[2] - b[2] * c[1]) / determinant,
            (a[2] * c[1] - a[1] * c[2]) / determinant,
            (a[1] * b[2] - a[2] * b[1]) / determinant,
        ),
        (
            (b[2] * c[0] - b[0] * c[2]) / determinant,
            (a[0] * c[2] - a[2] * c[0]) / determinant,
            (a[2] * b[0] - a[0] * b[2]) / determinant,
        ),
        (
            (b[0] * c[1] - b[1] * c[0]) / determinant,
            (a[1] * c[0] - a[0] * c[1]) / determinant,
            (a[0] * b[1] - a[1] * b[0]) / determinant,
        ),
    )


def _fractional_to_cartesian(
    fractional: Sequence[float], cell: Sequence[Sequence[float]]
) -> Tuple[float, float, float]:
    return tuple(
        sum(float(fractional[row]) * float(cell[row][column]) for row in range(3))
        for column in range(3)
    )  # type: ignore[return-value]


def enumerate_periodic_observations(
    atoms: Iterable[GeometryAtom],
    cell_vectors_A: Sequence[Sequence[float]],
    native_edges: Iterable[NativeEdgeKey] = (),
    policy: DistancePolicy = DEFAULT_POLICY,
) -> Tuple[ContactObservation, ...]:
    """Enumerate physical periodic neighbours as orbit + integer translation.

    ``GeometryAtom.coordinates`` are fractional. ``cell_vectors_A`` contains
    the three Cartesian lattice vectors as rows.  Search bounds are derived
    from reciprocal-vector norms, so non-orthogonal and short cells are not
    truncated to a fixed [-1, 1] image cube.
    """

    atom_list = tuple(atoms)
    atom_by_key = {atom.key: atom for atom in atom_list}
    if len(atom_by_key) != len(atom_list):
        raise ValueError("GeometryAtom keys must be unique")
    cell = tuple(tuple(float(value) for value in row) for row in cell_vectors_A)
    if len(cell) != 3 or any(len(row) != 3 for row in cell):
        raise ValueError("cell_vectors_A must be a 3x3 matrix")
    inverse = _reciprocal_rows(cell)
    reciprocal_columns = tuple(
        tuple(inverse[row][column] for row in range(3)) for column in range(3)
    )
    bounds = tuple(
        int(math.ceil(policy.search_cap_A * _euclidean(column))) + 1
        for column in reciprocal_columns
    )
    native_set = set(native_edges)
    seen_native = set()
    output: List[ContactObservation] = []

    for metal in atom_list:
        if not metal.is_metal:
            continue
        for donor in atom_list:
            for i in range(-bounds[0], bounds[0] + 1):
                for j in range(-bounds[1], bounds[1] + 1):
                    for k in range(-bounds[2], bounds[2] + 1):
                        image = (i, j, k)
                        if donor.key == metal.key and image == (0, 0, 0):
                            continue
                        delta_fractional = tuple(
                            float(donor.coordinates[index])
                            + image[index]
                            - float(metal.coordinates[index])
                            for index in range(3)
                        )
                        distance = _euclidean(_fractional_to_cartesian(delta_fractional, cell))
                        edge = NativeEdgeKey(metal.key, donor.key, image)
                        if edge in native_set:
                            seen_native.add(edge)
                        if distance <= policy.search_cap_A or edge in native_set:
                            output.append(
                                _make_observation(
                                    metal,
                                    donor,
                                    distance,
                                    native_edge=edge in native_set,
                                    periodic=True,
                                    image_delta=image,
                                )
                            )

    missing = native_set - seen_native
    if missing:
        unknown = sorted(
            (edge.metal_key, edge.donor_key, edge.image_delta) for edge in missing
            if edge.metal_key not in atom_by_key or edge.donor_key not in atom_by_key
        )
        if unknown:
            raise ValueError(f"native edges reference unknown atoms: {unknown}")
        invalid_centres = sorted(
            edge.metal_key for edge in missing if not atom_by_key[edge.metal_key].is_metal
        )
        if invalid_centres:
            raise ValueError(f"native edges use non-target metal centres: {invalid_centres}")
        # A native edge can lie beyond the image bounds derived from the search
        # cap.  It is still retained, because native provenance outranks v2.
        for edge in sorted(missing, key=lambda item: (item.metal_key, item.donor_key, item.image_delta)):
            metal = atom_by_key[edge.metal_key]
            donor = atom_by_key[edge.donor_key]
            delta_fractional = tuple(
                float(donor.coordinates[index])
                + edge.image_delta[index]
                - float(metal.coordinates[index])
                for index in range(3)
            )
            distance = _euclidean(_fractional_to_cartesian(delta_fractional, cell))
            output.append(
                _make_observation(
                    metal,
                    donor,
                    distance,
                    native_edge=True,
                    periodic=True,
                    image_delta=edge.image_delta,
                )
            )
    output.sort(key=lambda item: (item.metal_key, item.donor_key, item.image_delta or (999, 999, 999)))
    return tuple(output)


@dataclass(frozen=True)
class ContactDecision:
    metal_key: str
    metal_element: str
    metal_block: str
    donor_key: str
    donor_element: str
    donor_orbit_key: str
    image_delta: Optional[Tuple[int, int, int]]
    distance_A: Optional[float]
    normalized_distance: Optional[float]
    status: str
    evidence_level: str
    reason: str
    periodic: bool
    disorder_state: str
    pi_component_key: str
    ligand_component_key: str
    policy_id: str

    def __post_init__(self) -> None:
        if self.status not in CONTACT_STATUSES:
            raise ValueError(f"unsupported status: {self.status}")

    def to_dict(self) -> Dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class CollectiveSiteCandidate:
    """Distance-supported multiatom site without an inferred eta value."""

    metal_key: str
    pi_component_key: str
    member_keys: Tuple[str, ...]
    member_elements: Tuple[str, ...]
    member_distances_A: Tuple[float, ...]
    status: str = "candidate"
    hapticity: Optional[int] = None
    reason: str = "HAPTIC_ASSIGNMENT_WITHHELD"
    evidence_level: str = EVIDENCE_DISTANCE_CANDIDATE

    def to_dict(self) -> Dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class DistanceFirstSphereResult:
    policy_id: str
    policy_digest: str
    native_edges: Tuple[ContactDecision, ...]
    inferred_contacts: Tuple[ContactDecision, ...]
    collective_site_candidates: Tuple[CollectiveSiteCandidate, ...]

    @property
    def accepted(self) -> Tuple[ContactDecision, ...]:
        return tuple(item for item in self.inferred_contacts if item.status == "accepted")

    @property
    def candidates(self) -> Tuple[ContactDecision, ...]:
        return tuple(item for item in self.inferred_contacts if item.status == "candidate")

    @property
    def rejected(self) -> Tuple[ContactDecision, ...]:
        return tuple(item for item in self.inferred_contacts if item.status == "rejected")

    def to_dict(self) -> Dict[str, object]:
        return {
            "policy_id": self.policy_id,
            "policy_digest": self.policy_digest,
            "native_edges": [item.to_dict() for item in self.native_edges],
            "inferred_contacts": [item.to_dict() for item in self.inferred_contacts],
            "collective_site_candidates": [
                item.to_dict() for item in self.collective_site_candidates
            ],
        }


def _finite_distance(value: Optional[float]) -> bool:
    return value is not None and math.isfinite(float(value))


def _occupancy_requires_adjudication(observation: ContactObservation, floor: float) -> bool:
    values = (observation.metal_occupancy, observation.donor_occupancy)
    return any(value is None or float(value) < floor for value in values)


def _radius_sum(observation: ContactObservation) -> Optional[float]:
    metal_radius = COVALENT_RADII_ANGSTROM.get(observation.metal_element)
    donor_radius = COVALENT_RADII_ANGSTROM.get(observation.donor_element)
    if metal_radius is None or donor_radius is None:
        return None
    return float(metal_radius) + float(donor_radius)


def _shell_cutoffs(
    observations: Sequence[ContactObservation], policy: DistancePolicy
) -> Dict[str, Optional[float]]:
    """Return a conservative first-gap cutoff for each metal site.

    A native edge anchors the shell: no inferred gap is allowed inside the
    outermost native contact.  Distances beyond the block-specific candidate
    gate still participate as possible second-shell witnesses up to the frozen
    global search cap.
    """

    grouped: Dict[str, List[ContactObservation]] = defaultdict(list)
    for observation in observations:
        if _finite_distance(observation.distance_A):
            distance = float(observation.distance_A)
            if policy.minimum_distance_A <= distance <= policy.search_cap_A:
                grouped[observation.metal_key].append(observation)

    output: Dict[str, Optional[float]] = {}
    for metal_key, group in grouped.items():
        distances = sorted({round(float(item.distance_A), 6) for item in group})
        native_distances = [
            float(item.distance_A) for item in group if item.native_edge and item.distance_A is not None
        ]
        anchor = max(native_distances) if native_distances else distances[0]
        cutoff: Optional[float] = None
        for left, right in zip(distances, distances[1:]):
            if left + 1e-9 < anchor:
                continue
            required = max(policy.shell_gap_absolute_A, policy.shell_gap_relative * left)
            if right - left >= required:
                cutoff = (left + right) / 2.0
                break
        output[metal_key] = cutoff
    return output


def _decision(
    observation: ContactObservation,
    policy: DistancePolicy,
    *,
    status: str,
    evidence: str,
    reason: str,
    normalized_distance: Optional[float],
) -> ContactDecision:
    distance = observation.distance_A
    if distance is not None and math.isfinite(float(distance)):
        distance = round(float(distance), policy.distance_round_digits)
    else:
        distance = None
    if normalized_distance is not None:
        normalized_distance = round(float(normalized_distance), 6)
    return ContactDecision(
        metal_key=observation.metal_key,
        metal_element=observation.metal_element,
        metal_block=observation.metal_block,
        donor_key=observation.donor_key,
        donor_element=observation.donor_element,
        donor_orbit_key=observation.donor_orbit_key or observation.donor_key,
        image_delta=observation.image_delta,
        distance_A=distance,
        normalized_distance=normalized_distance,
        status=status,
        evidence_level=evidence,
        reason=reason,
        periodic=observation.periodic,
        disorder_state=observation.disorder_state,
        pi_component_key=observation.pi_component_key,
        ligand_component_key=observation.ligand_component_key,
        policy_id=policy.policy_id,
    )


def classify_distance_first_sphere(
    observations: Iterable[ContactObservation],
    policy: DistancePolicy = DEFAULT_POLICY,
) -> DistanceFirstSphereResult:
    """Classify native and non-native observations without changing v1 data.

    Rejected means "not accepted by this frozen protocol", not proof that the
    atom pair is chemically non-interacting.  Candidate decisions likewise do
    not increment the auxiliary accepted coordination number.
    """

    items = tuple(observations)
    shell_cutoffs = _shell_cutoffs(items, policy)

    pi_members: Dict[Tuple[str, str], set] = defaultdict(set)
    for item in items:
        if not item.native_edge and item.pi_component_key:
            pi_members[(item.metal_key, item.pi_component_key)].add(item.donor_key)
    collective_keys = {key for key, members in pi_members.items() if len(members) >= 2}

    native: List[ContactDecision] = []
    inferred: List[ContactDecision] = []

    for item in items:
        radius_sum = _radius_sum(item)
        normalized = None
        if radius_sum is not None and _finite_distance(item.distance_A):
            normalized = float(item.distance_A) / radius_sum

        if item.native_edge:
            native.append(
                _decision(
                    item,
                    policy,
                    status="native",
                    evidence=EVIDENCE_NATIVE,
                    reason="NATIVE_EDGE_IMMUTABLE",
                    normalized_distance=normalized,
                )
            )
            continue

        if not _finite_distance(item.distance_A):
            inferred.append(
                _decision(
                    item, policy, status="rejected", evidence=EVIDENCE_PROTOCOL_REJECTED,
                    reason="COORDINATE_DISTANCE_UNAVAILABLE", normalized_distance=None,
                )
            )
            continue

        distance = float(item.distance_A)
        if distance < policy.minimum_distance_A:
            inferred.append(
                _decision(
                    item, policy, status="rejected", evidence=EVIDENCE_PROTOCOL_REJECTED,
                    reason="GEOMETRIC_OVERLAP", normalized_distance=normalized,
                )
            )
            continue
        if distance > policy.search_cap_A:
            inferred.append(
                _decision(
                    item, policy, status="rejected", evidence=EVIDENCE_PROTOCOL_REJECTED,
                    reason="BEYOND_GLOBAL_SEARCH_CAP", normalized_distance=normalized,
                )
            )
            continue
        if item.donor_is_metal:
            inferred.append(
                _decision(
                    item, policy, status="candidate", evidence=EVIDENCE_DISTANCE_CANDIDATE,
                    reason="METAL_METAL_INFERENCE_WITHHELD", normalized_distance=normalized,
                )
            )
            continue
        if item.periodic and item.image_delta is None:
            inferred.append(
                _decision(
                    item, policy, status="candidate", evidence=EVIDENCE_DISTANCE_CANDIDATE,
                    reason="PBC_IMAGE_UNRESOLVED", normalized_distance=normalized,
                )
            )
            continue
        if item.disorder_state == "incompatible":
            inferred.append(
                _decision(
                    item, policy, status="rejected", evidence=EVIDENCE_PROTOCOL_REJECTED,
                    reason="DISORDER_ALTERNATIVES_INCOMPATIBLE", normalized_distance=normalized,
                )
            )
            continue
        if item.disorder_state == "ambiguous" or _occupancy_requires_adjudication(
            item, policy.full_occupancy_floor
        ):
            inferred.append(
                _decision(
                    item, policy, status="candidate", evidence=EVIDENCE_DISTANCE_CANDIDATE,
                    reason="DISORDER_ADJUDICATION_REQUIRED", normalized_distance=normalized,
                )
            )
            continue
        if (item.metal_key, item.pi_component_key) in collective_keys:
            inferred.append(
                _decision(
                    item, policy, status="candidate", evidence=EVIDENCE_DISTANCE_CANDIDATE,
                    reason="HAPTIC_ASSIGNMENT_WITHHELD", normalized_distance=normalized,
                )
            )
            continue
        if radius_sum is None:
            inferred.append(
                _decision(
                    item, policy, status="rejected", evidence=EVIDENCE_PROTOCOL_REJECTED,
                    reason="FROZEN_RADIUS_UNAVAILABLE", normalized_distance=None,
                )
            )
            continue
        if normalized is None or normalized > policy.candidate_scale(item.metal_block):
            inferred.append(
                _decision(
                    item, policy, status="rejected", evidence=EVIDENCE_PROTOCOL_REJECTED,
                    reason="BEYOND_BLOCK_CANDIDATE_CUTOFF", normalized_distance=normalized,
                )
            )
            continue

        shell_cutoff = shell_cutoffs.get(item.metal_key)
        if shell_cutoff is not None and distance > shell_cutoff:
            inferred.append(
                _decision(
                    item, policy, status="rejected", evidence=EVIDENCE_PROTOCOL_REJECTED,
                    reason="OUTSIDE_FIRST_SHELL_GAP", normalized_distance=normalized,
                )
            )
            continue

        if item.donor_element not in policy.auto_accept_donor_elements:
            inferred.append(
                _decision(
                    item, policy, status="candidate", evidence=EVIDENCE_DISTANCE_CANDIDATE,
                    reason="DONOR_CHEMISTRY_ADJUDICATION_REQUIRED", normalized_distance=normalized,
                )
            )
            continue

        if normalized <= policy.strict_scale(item.metal_block):
            inferred.append(
                _decision(
                    item, policy, status="accepted", evidence=EVIDENCE_DISTANCE_STRONG,
                    reason="WITHIN_BLOCK_STRICT_CUTOFF", normalized_distance=normalized,
                )
            )
            continue

        if (
            item.metal_block in policy.ionic_extension_blocks
            and item.donor_element in policy.ionic_extension_donor_elements
            and normalized <= policy.ionic_shell_scale(item.metal_block)
            and shell_cutoff is not None
        ):
            inferred.append(
                _decision(
                    item, policy, status="accepted", evidence=EVIDENCE_DISTANCE_SHELL,
                    reason="IONIC_EXTENSION_WITH_FIRST_SHELL_GAP", normalized_distance=normalized,
                )
            )
            continue

        inferred.append(
            _decision(
                item, policy, status="candidate", evidence=EVIDENCE_DISTANCE_CANDIDATE,
                reason="WITHIN_CANDIDATE_BAND_ONLY", normalized_distance=normalized,
            )
        )

    collective: List[CollectiveSiteCandidate] = []
    for metal_key, pi_key in sorted(collective_keys):
        members = [
            item
            for item in items
            if not item.native_edge
            and item.metal_key == metal_key
            and item.pi_component_key == pi_key
            and _finite_distance(item.distance_A)
        ]
        members.sort(key=lambda item: (item.donor_key, float(item.distance_A)))
        collective.append(
            CollectiveSiteCandidate(
                metal_key=metal_key,
                pi_component_key=pi_key,
                member_keys=tuple(item.donor_key for item in members),
                member_elements=tuple(item.donor_element for item in members),
                member_distances_A=tuple(
                    round(float(item.distance_A), policy.distance_round_digits) for item in members
                ),
            )
        )

    sort_key = lambda decision: (
        decision.metal_key,
        decision.donor_orbit_key,
        decision.image_delta if decision.image_delta is not None else (999, 999, 999),
        decision.donor_key,
    )
    return DistanceFirstSphereResult(
        policy_id=policy.policy_id,
        policy_digest=policy.digest(),
        native_edges=tuple(sorted(native, key=sort_key)),
        inferred_contacts=tuple(sorted(inferred, key=sort_key)),
        collective_site_candidates=tuple(collective),
    )


def summarize_distance_results(
    results: Iterable[DistanceFirstSphereResult],
) -> Dict[str, object]:
    """Return the minimum release statistics required for this v2 layer."""

    native_count = 0
    decisions: List[ContactDecision] = []
    collective_count = 0
    for result in results:
        native_count += len(result.native_edges)
        decisions.extend(result.inferred_contacts)
        collective_count += len(result.collective_site_candidates)

    by_status = Counter(item.status for item in decisions)
    by_block_status = Counter((item.metal_block, item.status) for item in decisions)
    by_reason = Counter(item.reason for item in decisions)
    by_evidence = Counter(item.evidence_level for item in decisions)
    flag_counts = {
        "periodic_decisions": sum(bool(item.periodic) for item in decisions),
        "disorder_affected_decisions": sum(item.disorder_state != "clear" for item in decisions),
        "haptic_member_decisions": sum(bool(item.pi_component_key) for item in decisions),
        "collective_site_candidates": collective_count,
    }
    return {
        "native_edges_unchanged": native_count,
        "distance_observations": len(decisions),
        "by_status": dict(sorted(by_status.items())),
        "by_block_and_status": {
            f"{block}:{status}": count
            for (block, status), count in sorted(by_block_status.items())
        },
        "by_evidence": dict(sorted(by_evidence.items())),
        "by_reason": dict(sorted(by_reason.items())),
        "flags": flag_counts,
    }


__all__ = [
    "COVALENT_RADII_ANGSTROM",
    "ContactDecision",
    "ContactObservation",
    "CollectiveSiteCandidate",
    "DEFAULT_POLICY",
    "DistanceFirstSphereResult",
    "DistancePolicy",
    "GeometryAtom",
    "NativeEdgeKey",
    "PROTOCOL_ID",
    "classify_distance_first_sphere",
    "enumerate_molecular_observations",
    "enumerate_periodic_observations",
    "summarize_distance_results",
]
