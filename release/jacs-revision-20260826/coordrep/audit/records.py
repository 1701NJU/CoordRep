"""Typed records for the release-scale CSD audit.

This module defines a conservative envelope around the scope-native CoordRep
serializers.  A record may be resolved, partial, or ambiguous; those statuses
are chemical outcomes, not exceptions.  Failed input and failed pipeline
states are kept outside :class:`EntryRecord` and therefore cannot inflate the
structural-record emission rate.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass, field, replace
from typing import Any, Dict, List, Optional, Sequence, Tuple


SCHEMA_VERSION = "CoordRep-Record/3"

ENTRY_RESOLUTIONS = {
    "EMITTED_RESOLVED",
    "EMITTED_PARTIAL",
    "EMITTED_AMBIGUOUS",
    "EMITTED_AUDIT_ONLY",
}
SITE_RESOLUTIONS = {"resolved", "partial", "ambiguous", "unresolved"}
RECORD_LEVELS = {"local_state", "topology", "site_object", "audit_only"}
DONOR_KINDS = {"atom", "pi_fragment", "candidate_atom_set"}
CANONICAL_STATUSES = {"exact_label_free", "source_order_noncanonical"}
ANGLE_CLASSES = {"cis", "trans", "other", "unavailable"}


def _normalise(value: Any) -> Any:
    """Return JSON-ready data with deterministic mapping and tuple handling."""
    if hasattr(value, "to_dict"):
        value = value.to_dict()
    if isinstance(value, dict):
        return {str(key): _normalise(value[key]) for key in sorted(value)}
    if isinstance(value, (list, tuple)):
        return [_normalise(item) for item in value]
    return value


def canonical_json(value: Any) -> str:
    """Serialize *value* to the canonical JSON used for hashes and reruns."""
    return json.dumps(
        _normalise(value),
        ensure_ascii=True,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    )


def record_digest(value: Any) -> str:
    """Return the complete SHA-256 digest of a structural record body."""
    return hashlib.sha256(canonical_json(value).encode("utf-8")).hexdigest()


@dataclass(frozen=True)
class DonorGroup:
    """One atom-resolved or collective coordination site.

    ``topology_keys`` are label-free graph-refinement witnesses.  They retain
    source connectivity without exporting CSD atom labels or coordinates.
    ``target_metals`` refers to serializer-local metal labels assigned by the
    adapter; exact canonical status is reported separately.
    """

    label: str
    kind: str
    donor_elements: Tuple[str, ...]
    topology_keys: Tuple[str, ...]
    target_metals: Tuple[str, ...]
    hapticity: int = 1
    bridge_degree: int = 1
    external_donor_contact_count: int = 0
    ligand_component_key: str = ""
    periodic_orbit_key: str = ""
    bond_types: Tuple[str, ...] = ()
    distance_profile: Tuple[float, ...] = ()
    evidence: str = "CSD_native_bond"
    issues: Tuple[str, ...] = ()

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class CoordinationIncidence:
    """One physical metal--donor contact in the source record.

    Molecular records use ``image_delta=(0, 0, 0)``.  Periodic adapters must
    retain the crystallographic image delta rather than merging contacts to
    the same wrapped donor orbit.
    """

    metal: str
    donor_group: str
    donor_ordinal: int
    donor_element: str
    bond_role: str
    bond_type: str
    image_delta: Tuple[int, int, int] = (0, 0, 0)
    image_delta_resolved: bool = True
    distance: Optional[float] = None

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class MetalRelation:
    metal_a: str
    metal_b: str
    relation: str
    through_groups: Tuple[str, ...] = ()
    bond_type: str = ""
    distance: Optional[float] = None
    image_delta: Tuple[int, int, int] = (0, 0, 0)
    image_delta_resolved: bool = True
    evidence: str = "CSD_native_bond"

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class DonorPairRelation:
    """Geometry-resolved relation between two coordination-site objects."""

    metal: str
    donor_group_a: str
    donor_group_b: str
    angle_deg: Optional[float]
    angle_class: str
    centroid_distance: Optional[float] = None
    evidence: str = "coordinate_geometry"

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class DonorGroupRelation:
    """Explicit membership relation between coordination-site objects.

    Collective haptic sites may share one or more member atoms with a second
    atom-resolved or collective site directed to another metal.  The ordinal
    pairs retain that eta/sigma or eta/eta overlap without merging chemically
    distinct, target-specific coordination-site objects.
    """

    donor_group_a: str
    donor_group_b: str
    relation: str
    shared_member_ordinals: Tuple[Tuple[int, int], ...]
    evidence: str = "CSD_native_bond_membership"

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class ExternalCenterRelation:
    """A source bond from an in-domain metal site to a metal outside policy.

    The partner is retained as a typed co-center rather than being treated as
    a ligand donor.  Under the frozen all-metal policy this class is reserved
    for forward compatibility with symbols outside the frozen API definition.
    """

    metal: str
    partner_element: str
    partner_topology_key: str
    bond_type: str
    distance: Optional[float] = None
    image_delta: Tuple[int, int, int] = (0, 0, 0)
    image_delta_resolved: bool = True
    evidence: str = "CSD_native_bond"

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class ExternalDonorRelation:
    """A donor-group member shared with a metal outside the audit domain."""

    donor_group: str
    donor_ordinal: int
    partner_element: str
    partner_topology_key: str
    bond_type: str
    distance: Optional[float] = None
    image_delta: Tuple[int, int, int] = (0, 0, 0)
    image_delta_resolved: bool = True
    evidence: str = "CSD_native_bond"

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class MetalSiteRecord:
    label: str
    element: str
    donor_groups: Tuple[str, ...]
    cn_site: int
    cn_atom: int
    eta_sum: Optional[int]
    metal_degree: int
    external_metal_degree: int
    external_bridge_count: int
    donor_composition: Tuple[Tuple[str, int], ...]
    geometry_kind: str
    shape_best: str = ""
    shape_values: Tuple[float, ...] = ()
    radial_profile: Tuple[float, ...] = ()
    angular_profile: Tuple[float, ...] = ()
    record_level: str = "topology"
    resolution: str = "partial"
    canonical_status: str = "source_order_noncanonical"
    issues: Tuple[str, ...] = ()
    coordrep_state: str = ""

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class EntryRecord:
    """One source entry with complete in-domain metal-site accounting."""

    source_index: int
    refcode: str
    scope: str
    scope_flags: Tuple[str, ...]
    expected_metal_sites: int
    metal_sites: Tuple[MetalSiteRecord, ...]
    donor_groups: Tuple[DonorGroup, ...]
    incidences: Tuple[CoordinationIncidence, ...]
    metal_relations: Tuple[MetalRelation, ...]
    external_center_relations: Tuple[ExternalCenterRelation, ...]
    external_donor_relations: Tuple[ExternalDonorRelation, ...]
    donor_group_relations: Tuple[DonorGroupRelation, ...]
    donor_pair_relations: Tuple[DonorPairRelation, ...]
    entry_resolution: str
    issue_codes: Tuple[str, ...] = ()
    schema_version: str = SCHEMA_VERSION
    source_database: str = "CSD"
    source_release: str = ""
    protocol_id: str = ""
    record_checksum_sha256: str = ""

    def structural_body(self) -> Dict[str, Any]:
        """Return fields that define the structural record digest.

        Source index, refcode, release, protocol, and issue bookkeeping are
        provenance rather than chemical identity and are excluded.
        """
        return {
            "schema_version": self.schema_version,
            "scope": self.scope,
            "scope_flags": self.scope_flags,
            "metal_sites": self.metal_sites,
            "donor_groups": self.donor_groups,
            "incidences": self.incidences,
            "metal_relations": self.metal_relations,
            "external_center_relations": self.external_center_relations,
            "external_donor_relations": self.external_donor_relations,
            "donor_group_relations": self.donor_group_relations,
            "donor_pair_relations": self.donor_pair_relations,
        }

    def with_digest(self) -> "EntryRecord":
        return replace(
            self,
            record_checksum_sha256=record_digest(self.structural_body()),
        )

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


def validate_entry_record(record: EntryRecord) -> List[str]:
    """Return validation issue codes; an empty list means schema-valid."""
    issues: List[str] = []
    if record.schema_version != SCHEMA_VERSION:
        issues.append("SCHEMA_VERSION_MISMATCH")
    if record.entry_resolution not in ENTRY_RESOLUTIONS:
        issues.append("INVALID_ENTRY_RESOLUTION")
    if record.expected_metal_sites != len(record.metal_sites):
        issues.append("METAL_SITE_ACCOUNTING_MISMATCH")
    if record.expected_metal_sites < 1:
        issues.append("NO_IN_DOMAIN_METAL_SITE")

    metal_labels = [site.label for site in record.metal_sites]
    if len(metal_labels) != len(set(metal_labels)):
        issues.append("DUPLICATE_METAL_LABEL")
    metal_label_set = set(metal_labels)

    group_labels = [group.label for group in record.donor_groups]
    if len(group_labels) != len(set(group_labels)):
        issues.append("DUPLICATE_DONOR_GROUP_LABEL")
    group_label_set = set(group_labels)
    group_by_label = {group.label: group for group in record.donor_groups}

    incidence_keys = [
        (
            incidence.metal,
            incidence.donor_group,
            incidence.donor_ordinal,
            incidence.image_delta,
            incidence.bond_role,
        )
        for incidence in record.incidences
    ]
    if len(incidence_keys) != len(set(incidence_keys)):
        issues.append("DUPLICATE_COORDINATION_INCIDENCE")
    for incidence in record.incidences:
        if incidence.metal not in metal_label_set:
            issues.append("UNKNOWN_INCIDENCE_METAL")
        if incidence.donor_group not in group_label_set:
            issues.append("UNKNOWN_INCIDENCE_GROUP")
        group = next((g for g in record.donor_groups if g.label == incidence.donor_group), None)
        if group is not None:
            if not 0 <= incidence.donor_ordinal < len(group.donor_elements):
                issues.append("INVALID_DONOR_ORDINAL")
            elif incidence.donor_element != group.donor_elements[incidence.donor_ordinal]:
                issues.append("INCIDENCE_DONOR_ELEMENT_MISMATCH")
        if len(incidence.image_delta) != 3 or not all(
            isinstance(value, int) for value in incidence.image_delta
        ):
            issues.append("INVALID_IMAGE_DELTA")
        if not isinstance(incidence.image_delta_resolved, bool):
            issues.append("INVALID_IMAGE_DELTA_RESOLUTION")

    for group in record.donor_groups:
        if group.kind not in DONOR_KINDS:
            issues.append("INVALID_DONOR_KIND")
        if group.hapticity < 1 or group.bridge_degree < 1:
            issues.append("INVALID_COORDINATION_MULTIPLICITY")
        if group.external_donor_contact_count < 0:
            issues.append("INVALID_EXTERNAL_DONOR_CONTACT_COUNT")
        if group.bridge_degree != len(set(group.target_metals)):
            issues.append("BRIDGE_DEGREE_MISMATCH")
        if not set(group.target_metals).issubset(metal_label_set):
            issues.append("UNKNOWN_DONOR_TARGET")
        if len(group.donor_elements) != len(group.topology_keys):
            issues.append("DONOR_TOPOLOGY_LENGTH_MISMATCH")
        incidence_targets = {
            incidence.metal
            for incidence in record.incidences
            if incidence.donor_group == group.label
        }
        if incidence_targets != set(group.target_metals):
            issues.append("DONOR_TARGET_INCIDENCE_MISMATCH")
        external_relation_count = sum(
            relation.donor_group == group.label
            for relation in record.external_donor_relations
        )
        if group.external_donor_contact_count != external_relation_count:
            issues.append("EXTERNAL_DONOR_CONTACT_COUNT_MISMATCH")

    for site in record.metal_sites:
        if site.record_level not in RECORD_LEVELS:
            issues.append("INVALID_RECORD_LEVEL")
        if site.resolution not in SITE_RESOLUTIONS:
            issues.append("INVALID_SITE_RESOLUTION")
        if site.canonical_status not in CANONICAL_STATUSES:
            issues.append("INVALID_CANONICAL_STATUS")
        if not set(site.donor_groups).issubset(group_label_set):
            issues.append("UNKNOWN_LOCAL_DONOR_GROUP")
        local_incidences = [i for i in record.incidences if i.metal == site.label]
        local_group_labels = {i.donor_group for i in local_incidences}
        local_groups = [g for g in record.donor_groups if g.label in local_group_labels]
        if site.cn_site != len(local_groups):
            issues.append("CN_SITE_MISMATCH")
        if set(site.donor_groups) != local_group_labels:
            issues.append("LOCAL_DONOR_GROUP_LIST_MISMATCH")
        if site.cn_atom != len(local_incidences):
            issues.append("CN_ATOM_MISMATCH")
        has_candidate_group = any(
            group.kind == "candidate_atom_set" for group in local_groups
        )
        expected_eta_sum = sum(group.hapticity for group in local_groups)
        if has_candidate_group:
            if site.eta_sum is not None:
                issues.append("AMBIGUOUS_ETA_SUM_MUST_BE_NULL")
        elif site.eta_sum != expected_eta_sum:
            issues.append("ETA_SUM_MISMATCH")
        expected_composition: Dict[str, int] = {}
        for incidence in local_incidences:
            expected_composition[incidence.donor_element] = (
                expected_composition.get(incidence.donor_element, 0) + 1
            )
        if tuple(sorted(expected_composition.items())) != site.donor_composition:
            issues.append("DONOR_COMPOSITION_MISMATCH")
        direct_degree = sum(
            relation.relation == "direct_metal_bond"
            and site.label in (relation.metal_a, relation.metal_b)
            for relation in record.metal_relations
        )
        if site.metal_degree != direct_degree:
            issues.append("METAL_DEGREE_MISMATCH")
        external_degree = sum(
            relation.metal == site.label
            for relation in record.external_center_relations
        )
        if site.external_metal_degree != external_degree:
            issues.append("EXTERNAL_METAL_DEGREE_MISMATCH")
        external_bridge_count = sum(
            relation.donor_group in local_group_labels
            for relation in record.external_donor_relations
        )
        if site.external_bridge_count != external_bridge_count:
            issues.append("EXTERNAL_BRIDGE_COUNT_MISMATCH")

        local_pair_relations = [
            relation
            for relation in record.donor_pair_relations
            if relation.metal == site.label
        ]
        expected_pairs = site.cn_site * (site.cn_site - 1) // 2
        if len(local_pair_relations) != expected_pairs:
            issues.append("DONOR_PAIR_RELATION_COUNT_MISMATCH")

    pair_keys = []
    local_group_sets = {
        site.label: set(site.donor_groups) for site in record.metal_sites
    }
    for relation in record.donor_pair_relations:
        if relation.metal not in metal_label_set:
            issues.append("UNKNOWN_DONOR_PAIR_METAL")
            continue
        if relation.donor_group_a >= relation.donor_group_b:
            issues.append("NONCANONICAL_DONOR_PAIR_ORDER")
        if relation.angle_class not in ANGLE_CLASSES:
            issues.append("INVALID_DONOR_PAIR_ANGLE_CLASS")
        if relation.angle_class == "unavailable" and relation.angle_deg is not None:
            issues.append("UNAVAILABLE_DONOR_PAIR_HAS_ANGLE")
        if relation.angle_class != "unavailable" and relation.angle_deg is None:
            issues.append("DONOR_PAIR_ANGLE_MISSING")
        local_groups = local_group_sets[relation.metal]
        if not {relation.donor_group_a, relation.donor_group_b}.issubset(local_groups):
            issues.append("DONOR_PAIR_NOT_LOCAL_TO_METAL")
        pair_keys.append((
            relation.metal,
            relation.donor_group_a,
            relation.donor_group_b,
        ))
    if len(pair_keys) != len(set(pair_keys)):
        issues.append("DUPLICATE_DONOR_PAIR_RELATION")

    for relation in record.metal_relations:
        if relation.metal_a not in metal_label_set or relation.metal_b not in metal_label_set:
            issues.append("UNKNOWN_METAL_RELATION_ENDPOINT")
        if relation.metal_a >= relation.metal_b:
            issues.append("NONCANONICAL_METAL_RELATION_ORDER")
        if len(relation.image_delta) != 3 or not all(
            isinstance(value, int) for value in relation.image_delta
        ):
            issues.append("INVALID_METAL_RELATION_IMAGE_DELTA")
        if not isinstance(relation.image_delta_resolved, bool):
            issues.append("INVALID_METAL_RELATION_IMAGE_DELTA_RESOLUTION")

    for relation in record.external_center_relations:
        if relation.metal not in metal_label_set:
            issues.append("UNKNOWN_EXTERNAL_CENTER_RELATION_METAL")
        if not relation.partner_element or not relation.partner_topology_key:
            issues.append("INVALID_EXTERNAL_CENTER_PARTNER")
        if len(relation.image_delta) != 3 or not all(
            isinstance(value, int) for value in relation.image_delta
        ):
            issues.append("INVALID_EXTERNAL_CENTER_IMAGE_DELTA")
        if not isinstance(relation.image_delta_resolved, bool):
            issues.append("INVALID_EXTERNAL_CENTER_IMAGE_DELTA_RESOLUTION")

    for relation in record.external_donor_relations:
        if relation.donor_group not in group_label_set:
            issues.append("UNKNOWN_EXTERNAL_DONOR_GROUP")
            continue
        group = group_by_label[relation.donor_group]
        if not 0 <= relation.donor_ordinal < len(group.donor_elements):
            issues.append("INVALID_EXTERNAL_DONOR_ORDINAL")
        if not relation.partner_element or not relation.partner_topology_key:
            issues.append("INVALID_EXTERNAL_DONOR_PARTNER")
        if len(relation.image_delta) != 3 or not all(
            isinstance(value, int) for value in relation.image_delta
        ):
            issues.append("INVALID_EXTERNAL_DONOR_IMAGE_DELTA")
        if not isinstance(relation.image_delta_resolved, bool):
            issues.append("INVALID_EXTERNAL_DONOR_IMAGE_DELTA_RESOLUTION")

    group_relation_keys = []
    for relation in record.donor_group_relations:
        if relation.donor_group_a not in group_label_set or relation.donor_group_b not in group_label_set:
            issues.append("UNKNOWN_DONOR_GROUP_RELATION_ENDPOINT")
            continue
        if relation.donor_group_a >= relation.donor_group_b:
            issues.append("NONCANONICAL_DONOR_GROUP_RELATION_ORDER")
        if relation.relation != "overlaps_on_member":
            issues.append("INVALID_DONOR_GROUP_RELATION")
        group_a = group_by_label[relation.donor_group_a]
        group_b = group_by_label[relation.donor_group_b]
        if not relation.shared_member_ordinals:
            issues.append("EMPTY_DONOR_GROUP_OVERLAP")
        for ordinal_a, ordinal_b in relation.shared_member_ordinals:
            if not 0 <= ordinal_a < len(group_a.donor_elements):
                issues.append("INVALID_DONOR_GROUP_OVERLAP_ORDINAL")
            if not 0 <= ordinal_b < len(group_b.donor_elements):
                issues.append("INVALID_DONOR_GROUP_OVERLAP_ORDINAL")
        group_relation_keys.append((
            relation.donor_group_a,
            relation.donor_group_b,
            relation.relation,
        ))
    if len(group_relation_keys) != len(set(group_relation_keys)):
        issues.append("DUPLICATE_DONOR_GROUP_RELATION")

    if record.entry_resolution != "EMITTED_AUDIT_ONLY":
        has_direct_metal_relation = any(
            relation.relation == "direct_metal_bond" for relation in record.metal_relations
        )
        if (
            not record.donor_groups
            and not has_direct_metal_relation
            and not record.external_center_relations
        ):
            issues.append("STRUCTURAL_RECORD_WITHOUT_DONORS")
        if all(site.record_level == "audit_only" for site in record.metal_sites):
            issues.append("STRUCTURAL_RECORD_HAS_ONLY_AUDIT_SITES")

    expected_digest = record_digest(record.structural_body())
    if (
        record.record_checksum_sha256
        and record.record_checksum_sha256 != expected_digest
    ):
        issues.append("RECORD_DIGEST_MISMATCH")
    return sorted(set(issues))
