from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from coordrep.audit.csd_adapter import CSDRecordError, adapt_csd_entry
from coordrep.audit.records import canonical_json, validate_entry_record


@dataclass(frozen=True)
class FakeCoordinates:
    x: float
    y: float
    z: float


class FakeAtom:
    def __init__(self, label, element, xyz):
        self.label = label
        self.atomic_symbol = element
        self.coordinates = FakeCoordinates(*xyz)
        self.fractional_coordinates = FakeCoordinates(*xyz)
        self.formal_charge = None
        self.bonds = []


class FakeBond:
    def __init__(self, a, b, bond_type="Single"):
        self.atoms = (a, b)
        self.bond_type = bond_type
        a.bonds.append(self)
        b.bonds.append(self)


class FakeMolecule:
    def __init__(self, atoms, bonds):
        self.atoms = atoms
        self.bonds = bonds


class FakeEntry:
    def __init__(self, refcode, atoms, bonds, *, disorder=False, polymer=False):
        self.identifier = refcode
        self.has_3d_structure = True
        self.has_disorder = disorder
        self.is_polymeric = polymer
        self.molecule = FakeMolecule(atoms, bonds)


def _entry(refcode, atom_specs, bond_specs, **flags):
    atoms = [FakeAtom(*spec) for spec in atom_specs]
    for index, atom in enumerate(atoms):
        atom.index = index
    bonds = [FakeBond(atoms[a], atoms[b], kind) for a, b, kind in bond_specs]
    return FakeEntry(refcode, atoms, bonds, **flags)


def test_mononuclear_eta1_emits_resolved_structural_record():
    entry = _entry(
        "FAKEPT",
        [
            ("Pt1", "Pt", (0, 0, 0)),
            ("N1", "N", (1, 0, 0)),
            ("N2", "N", (-1, 0, 0)),
            ("N3", "N", (0, 1, 0)),
            ("N4", "N", (0, -1, 0)),
        ],
        [(0, i, "Single") for i in range(1, 5)],
    )
    record = adapt_csd_entry(entry, 7, source_release="test")
    assert record.scope == "mono_eta1"
    assert record.entry_resolution == "EMITTED_RESOLVED"
    assert len(record.metal_sites) == 1
    assert len(record.donor_groups) == 4
    assert len(record.incidences) == 4
    assert record.metal_sites[0].cn_site == 4
    assert record.metal_sites[0].cn_atom == 4
    assert record.metal_sites[0].geometry_kind == "CShM"
    assert not validate_entry_record(record)


def test_bridge_is_one_group_with_two_metal_incidences():
    entry = _entry(
        "FAKECU2",
        [
            ("Cu1", "Cu", (-1, 0, 0)),
            ("Cu2", "Cu", (1, 0, 0)),
            ("Cl1", "Cl", (0, 1, 0)),
            ("Cl2", "Cl", (0, -1, 0)),
        ],
        [
            (0, 2, "Single"), (1, 2, "Single"),
            (0, 3, "Single"), (1, 3, "Single"),
        ],
    )
    record = adapt_csd_entry(entry, 8, source_release="test")
    assert record.scope == "multi_eta1"
    assert len(record.donor_groups) == 2
    assert len(record.incidences) == 4
    assert all(group.bridge_degree == 2 for group in record.donor_groups)
    assert all(site.cn_site == 2 for site in record.metal_sites)
    assert sum(rel.relation == "shared_donor" for rel in record.metal_relations) == 1
    assert not validate_entry_record(record)


def test_haptic_fragments_preserve_site_and_atom_counts():
    atoms = [("Fe1", "Fe", (0, 0, 0))]
    bonds = []
    for ring in range(2):
        z = -1.5 if ring == 0 else 1.5
        offset = len(atoms)
        for i in range(5):
            angle = 2 * np.pi * i / 5
            atoms.append((f"C{ring}_{i}", "C", (np.cos(angle), np.sin(angle), z)))
            bonds.append((0, offset + i, "Pi"))
            bonds.append((offset + i, offset + (i + 1) % 5, "Aromatic"))
    entry = _entry("FAKECP2", atoms, bonds)
    record = adapt_csd_entry(entry, 9, source_release="test")
    assert record.scope == "mono_haptic"
    assert [group.hapticity for group in record.donor_groups] == [5, 5]
    assert len(record.incidences) == 10
    assert record.metal_sites[0].cn_site == 2
    assert record.metal_sites[0].cn_atom == 10
    assert record.metal_sites[0].record_level == "site_object"
    assert not validate_entry_record(record)


def test_roqmer_like_eta13_is_an_ambiguous_candidate_not_confirmed_haptic():
    atom_specs = [("Cr1", "Cr", (0, 0, 0))]
    bond_specs = []
    for index in range(13):
        angle = 2 * np.pi * index / 13
        atom_specs.append((
            f"N{index + 1}",
            "N",
            (2.0 * np.cos(angle), 2.0 * np.sin(angle), 0),
        ))
        bond_specs.append((0, index + 1, "Pi"))
        bond_specs.append((index + 1, (index + 1) % 13 + 1, "Single"))

    record = adapt_csd_entry(
        _entry(
            "FAKEROQMER04",
            atom_specs,
            bond_specs,
            disorder=True,
            polymer=True,
        ),
        29,
        source_release="test",
    )

    assert len(record.donor_groups) == 1
    group = record.donor_groups[0]
    assert group.kind == "candidate_atom_set"
    assert group.hapticity == 1
    assert "COLLECTIVE_SITE_MEMBERSHIP_OUT_OF_DOMAIN" in group.issues
    assert "ambiguous_pi_candidate" in record.scope_flags
    assert "haptic" not in record.scope_flags
    assert record.metal_sites[0].record_level == "topology"
    assert record.metal_sites[0].resolution == "ambiguous"
    assert "AMBIGUOUS_COLLECTIVE_PI_SITE" in record.metal_sites[0].issues
    assert not validate_entry_record(record)


def test_multimetal_haptic_disorder_polymer_is_not_forced_resolved():
    entry = _entry(
        "FAKEMIX",
        [
            ("Fe1", "Fe", (-1, 0, 0)),
            ("Fe2", "Fe", (1, 0, 0)),
            ("C1", "C", (0, 1, 0)),
            ("C2", "C", (0, 2, 0)),
        ],
        [
            (0, 2, "Pi"), (0, 3, "Pi"),
            (1, 2, "Pi"), (1, 3, "Pi"),
            (2, 3, "Single"),
        ],
        disorder=True,
        polymer=True,
    )
    record = adapt_csd_entry(entry, 10, source_release="test")
    assert record.scope == "multi_haptic"
    assert record.entry_resolution == "EMITTED_AMBIGUOUS"
    assert {"disorder", "polymeric", "haptic", "multimetal"}.issubset(record.scope_flags)
    assert all(site.resolution == "ambiguous" for site in record.metal_sites)
    assert not validate_entry_record(record)


def test_cn_above_reference_library_emits_partial_topology_record():
    atom_specs = [("La1", "La", (0, 0, 0))]
    for i in range(7):
        angle = 2 * np.pi * i / 7
        atom_specs.append((f"O{i}", "O", (np.cos(angle), np.sin(angle), 0.2 * (-1) ** i)))
    entry = _entry(
        "FAKECN7",
        atom_specs,
        [(0, i, "Single") for i in range(1, 8)],
    )
    record = adapt_csd_entry(entry, 11, source_release="test")
    assert record.entry_resolution == "EMITTED_PARTIAL"
    assert "high_cn" in record.scope_flags
    assert record.metal_sites[0].cn_atom == 7
    assert record.metal_sites[0].record_level == "topology"
    assert "REFERENCE_SET_UNAVAILABLE" in record.metal_sites[0].issues
    assert not validate_entry_record(record)


def test_zero_native_donors_is_audit_only_not_structural_coverage():
    entry = _entry("FAKEEMPTY", [("Zn1", "Zn", (0, 0, 0))], [])
    record = adapt_csd_entry(entry, 12, source_release="test")
    assert record.entry_resolution == "EMITTED_AUDIT_ONLY"
    assert record.metal_sites[0].record_level == "audit_only"
    assert not validate_entry_record(record)


def test_native_metal_hydride_is_retained_as_a_coordination_incidence():
    entry = _entry(
        "FAKEHYDRIDE",
        [("Fe1", "Fe", (0, 0, 0)), ("H1", "H", (1.5, 0, 0))],
        [(0, 1, "Single")],
    )
    record = adapt_csd_entry(entry, 15, source_release="test")
    assert record.entry_resolution == "EMITTED_PARTIAL"
    assert record.metal_sites[0].cn_atom == 1
    assert record.donor_groups[0].donor_elements == ("H",)
    assert record.incidences[0].donor_element == "H"
    assert not validate_entry_record(record)


def test_bridging_metal_bound_hydrogen_is_one_shared_site():
    entry = _entry(
        "FAKEBRIDGEH",
        [
            ("Ru1", "Ru", (-1, 0, 0)),
            ("Ru2", "Ru", (1, 0, 0)),
            ("H1", "H", (0, 0, 0)),
        ],
        [(0, 2, "Single"), (1, 2, "Single")],
    )
    record = adapt_csd_entry(entry, 17, source_release="test")
    assert len(record.donor_groups) == 1
    assert record.donor_groups[0].donor_elements == ("H",)
    assert record.donor_groups[0].bridge_degree == 2
    assert len(record.incidences) == 2
    assert "metal_hydrogen" in record.scope_flags
    assert not validate_entry_record(record)


def test_mixed_single_delocalised_singleton_is_one_shared_atom_site():
    entry = _entry(
        "FAKEMIXEDBRIDGE",
        [
            ("Fe1", "Fe", (-1, 0, 0)),
            ("Fe2", "Fe", (1, 0, 0)),
            ("O1", "O", (0, 0, 0)),
        ],
        [(0, 2, "Single"), (1, 2, "Delocalised")],
    )
    record = adapt_csd_entry(entry, 21, source_release="test")
    assert len(record.donor_groups) == 1
    group = record.donor_groups[0]
    assert group.kind == "atom"
    assert group.donor_elements == ("O",)
    assert group.bridge_degree == 2
    assert set(group.target_metals) == {site.label for site in record.metal_sites}
    assert {incidence.bond_role for incidence in record.incidences} == {"pi", "sigma"}
    assert all(site.cn_site == 1 for site in record.metal_sites)
    shared = [
        relation
        for relation in record.metal_relations
        if relation.relation == "shared_donor"
    ]
    assert len(shared) == 1
    assert shared[0].through_groups == (group.label,)
    assert not validate_entry_record(record)


def test_duplicate_endpoint_bond_labels_are_one_union_incidence():
    entry = _entry(
        "FAKEMIXEDENDPOINT",
        [("Fe1", "Fe", (0, 0, 0)), ("O1", "O", (1.8, 0, 0))],
        [(0, 1, "Single"), (0, 1, "Delocalised")],
    )
    record = adapt_csd_entry(entry, 23, source_release="test")
    assert len(record.donor_groups) == 1
    assert len(record.incidences) == 1
    incidence = record.incidences[0]
    assert incidence.bond_role == "pi+sigma"
    assert incidence.bond_type == "Delocalised+Single"
    assert record.metal_sites[0].cn_atom == 1
    assert record.metal_sites[0].cn_site == 1
    assert not validate_entry_record(record)


def test_eta5_site_and_sigma_atom_site_retain_shared_member_relation():
    atom_specs = [
        ("Fe1", "Fe", (0, 0, 0)),
        ("Cu1", "Cu", (2.5, 0, 1.5)),
    ]
    bond_specs = []
    ring_start = len(atom_specs)
    for index in range(5):
        angle = 2 * np.pi * index / 5
        atom_specs.append((
            f"C{index + 1}",
            "C",
            (np.cos(angle), np.sin(angle), 1.5),
        ))
        bond_specs.append((0, ring_start + index, "Pi"))
        bond_specs.append((
            ring_start + index,
            ring_start + (index + 1) % 5,
            "Aromatic",
        ))
    bond_specs.append((1, ring_start, "Single"))

    record = adapt_csd_entry(
        _entry("FAKEETA5SIGMA", atom_specs, bond_specs),
        22,
        source_release="test",
    )
    assert record.scope == "multi_haptic"
    assert len(record.donor_groups) == 2
    pi_group = next(group for group in record.donor_groups if group.kind == "pi_fragment")
    atom_group = next(group for group in record.donor_groups if group.kind == "atom")
    assert pi_group.hapticity == 5
    assert atom_group.hapticity == 1
    assert len(record.donor_group_relations) == 1
    overlap = record.donor_group_relations[0]
    assert overlap.relation == "overlaps_on_member"
    assert len(overlap.shared_member_ordinals) == 1
    ordinal_pi, ordinal_atom = (
        overlap.shared_member_ordinals[0]
        if overlap.donor_group_a == pi_group.label
        else tuple(reversed(overlap.shared_member_ordinals[0]))
    )
    assert pi_group.donor_elements[ordinal_pi] == "C"
    assert atom_group.donor_elements[ordinal_atom] == "C"
    shared_member_relations = [
        relation
        for relation in record.metal_relations
        if relation.relation == "shared_donor_member"
    ]
    assert len(shared_member_relations) == 1
    assert set(shared_member_relations[0].through_groups) == {
        pi_group.label,
        atom_group.label,
    }
    assert not validate_entry_record(record)


def test_direct_metal_bond_only_is_a_structural_topology_record():
    entry = _entry(
        "FAKEMM",
        [("Pt1", "Pt", (-1, 0, 0)), ("Pt2", "Pt", (1, 0, 0))],
        [(0, 1, "Single")],
    )
    record = adapt_csd_entry(entry, 16, source_release="test")
    assert record.entry_resolution == "EMITTED_PARTIAL"
    assert all(site.record_level == "topology" for site in record.metal_sites)
    assert all("DIRECT_METAL_RELATION_ONLY" in site.issues for site in record.metal_sites)
    assert len(record.metal_relations) == 1
    assert record.metal_relations[0].relation == "direct_metal_bond"
    assert not validate_entry_record(record)


def test_aluminium_is_promoted_to_an_internal_metal_site():
    entry = _entry(
        "FAKEFEAL",
        [("Fe1", "Fe", (0, 0, 0)), ("Al1", "Al", (2.3, 0, 0))],
        [(0, 1, "Single")],
    )
    record = adapt_csd_entry(entry, 24, source_release="test")
    assert record.entry_resolution == "EMITTED_PARTIAL"
    assert not record.donor_groups
    assert not record.incidences
    assert {site.element for site in record.metal_sites} == {"Al", "Fe"}
    assert len(record.metal_relations) == 1
    relation = record.metal_relations[0]
    assert relation.relation == "direct_metal_bond"
    assert relation.distance == 2.3
    assert not record.external_center_relations
    assert all(site.metal_degree == 1 for site in record.metal_sites)
    assert all(site.record_level == "topology" for site in record.metal_sites)
    assert {"metal_block_p", "metal_block_d", "mixed_metal_blocks"}.issubset(
        record.scope_flags
    )
    assert not validate_entry_record(record)


def test_donor_shared_with_aluminium_is_an_internal_bridge():
    record = adapt_csd_entry(
        _entry(
            "FAKETINAL",
            [
                ("Ti1", "Ti", (0, 0, 0)),
                ("N1", "N", (1.9, 0, 0)),
                ("Al1", "Al", (3.8, 0, 0)),
            ],
            [(0, 1, "Single"), (1, 2, "Single")],
        ),
        30,
        source_release="test",
    )

    assert record.entry_resolution == "EMITTED_PARTIAL"
    assert len(record.donor_groups) == 1
    assert len(record.incidences) == 2
    assert not record.external_center_relations
    assert not record.external_donor_relations
    assert record.donor_groups[0].bridge_degree == 2
    assert record.donor_groups[0].external_donor_contact_count == 0
    assert set(record.donor_groups[0].target_metals) == {
        site.label for site in record.metal_sites
    }
    assert any(
        relation.relation == "shared_donor" for relation in record.metal_relations
    )
    assert all(site.external_bridge_count == 0 for site in record.metal_sites)
    assert not validate_entry_record(record)


def test_s_and_f_block_metals_emit_structure_records():
    for refcode, metal in (("FAKEMG", "Mg"), ("FAKECE", "Ce")):
        record = adapt_csd_entry(
            _entry(
                refcode,
                [(f"{metal}1", metal, (0, 0, 0)), ("O1", "O", (2.2, 0, 0))],
                [(0, 1, "Single")],
            ),
            27,
            source_release="test",
        )
        assert record.expected_metal_sites == 1
        assert record.metal_sites[0].element == metal
        assert len(record.incidences) == 1
        expected_flag = "metal_block_s" if metal == "Mg" else "metal_block_f"
        assert expected_flag in record.scope_flags
        assert not validate_entry_record(record)


def test_polymeric_no_edge_is_not_translation_resolved():
    record = adapt_csd_entry(
        _entry("FAKEPBCEMPTY", [("Zn1", "Zn", (0.2, 0.3, 0.4))], [], polymer=True),
        25,
        source_release="test",
    )
    assert record.entry_resolution == "EMITTED_AUDIT_ONLY"
    assert "polymeric" in record.scope_flags
    assert "translation_resolved_source_graph" not in record.scope_flags
    assert "PBC_TRANSLATION_UNRESOLVED" in record.issue_codes
    assert not validate_entry_record(record)


def test_fractional_coordinates_are_not_used_as_cartesian_geometry():
    entry = _entry(
        "FAKEFRACTIONALONLY",
        [
            ("Pt1", "Pt", (0.1, 0.2, 0.3)),
            ("N1", "N", (0.2, 0.2, 0.3)),
            ("N2", "N", (0.0, 0.2, 0.3)),
        ],
        [(0, 1, "Single"), (0, 2, "Single")],
    )
    for atom in entry.molecule.atoms:
        atom.coordinates = None
    record = adapt_csd_entry(entry, 26, source_release="test")
    site = record.metal_sites[0]
    assert record.entry_resolution == "EMITTED_PARTIAL"
    assert site.geometry_kind == "topology-only"
    assert not site.radial_profile
    assert not site.angular_profile
    assert not site.shape_values
    assert "COORDINATES_UNAVAILABLE" in site.issues
    assert not validate_entry_record(record)


def test_polymeric_source_graph_retains_translation_label():
    entry = _entry(
        "FAKEPBC",
        [
            ("Cr1", "Cr", (0.8, 0.2, 0.2)),
            ("O1", "O", (0.9, -0.1, 0.2)),
            ("O2", "O", (0.7, 0.3, 0.2)),
        ],
        [(0, 1, "Delocalised"), (0, 2, "Single")],
        polymer=True,
    )
    record = adapt_csd_entry(entry, 14, source_release="test")
    assert record.entry_resolution == "EMITTED_RESOLVED"
    assert "translation_resolved_source_graph" in record.scope_flags
    assert {incidence.image_delta for incidence in record.incidences} == {
        (0, -1, 0),
        (0, 0, 0),
    }
    assert all(group.periodic_orbit_key for group in record.donor_groups)
    assert not validate_entry_record(record)


def test_record_digest_and_canonical_json_are_repeatable():
    entry = _entry(
        "FAKEREPEAT",
        [("Ni1", "Ni", (0, 0, 0)), ("N1", "N", (1, 0, 0)), ("N2", "N", (-1, 0, 0))],
        [(0, 1, "Single"), (0, 2, "Single")],
    )
    first = adapt_csd_entry(entry, 13, source_release="test")
    second = adapt_csd_entry(entry, 13, source_release="test")
    assert canonical_json(first) == canonical_json(second)
    assert first.record_checksum_sha256 == second.record_checksum_sha256


def test_label_free_group_order_is_invariant_to_atom_list_permutation():
    specs = [
        ("Pt1", "Pt", (0, 0, 0)),
        ("N1", "N", (1.0, 0, 0)),
        ("N2", "N", (-1.2, 0, 0)),
        ("Cl1", "Cl", (0, 1.4, 0)),
        ("Cl2", "Cl", (0, -1.6, 0)),
    ]
    bonds = [(0, 1, "Single"), (0, 2, "Single"), (0, 3, "Single"), (0, 4, "Single")]
    permuted = [specs[0], specs[2], specs[1], specs[4], specs[3]]
    first = adapt_csd_entry(_entry("FAKEPERM", specs, bonds), 18, source_release="test")
    second = adapt_csd_entry(_entry("FAKEPERM", permuted, bonds), 18, source_release="test")
    assert first.record_checksum_sha256 == second.record_checksum_sha256
    assert all(site.canonical_status == "exact_label_free" for site in first.metal_sites)
    assert all(site.canonical_status == "exact_label_free" for site in second.metal_sites)


def test_heteroatom_ligand_component_key_is_atom_permutation_invariant():
    specs = [
        ("Pt1", "Pt", (0, 0, 0)),
        ("N1", "N", (1, 0, 0)),
        ("C1", "C", (2, 0, 0)),
        ("O1", "O", (3, 0, 0)),
    ]
    bonds = [
        (0, 1, "Single"),
        (1, 2, "Single"),
        (2, 3, "Double"),
    ]
    permutation = [0, 2, 1, 3]
    old_to_new = {
        old_index: new_index
        for new_index, old_index in enumerate(permutation)
    }
    permuted_specs = [specs[index] for index in permutation]
    permuted_bonds = [
        (old_to_new[a], old_to_new[b], bond_type)
        for a, b, bond_type in bonds
    ]

    first = adapt_csd_entry(
        _entry("FAKEHETEROLIG", specs, bonds),
        28,
        source_release="test",
    )
    second = adapt_csd_entry(
        _entry("FAKEHETEROLIG", permuted_specs, permuted_bonds),
        28,
        source_release="test",
    )

    assert first.donor_groups[0].ligand_component_key == (
        second.donor_groups[0].ligand_component_key
    )
    assert first.record_checksum_sha256 == second.record_checksum_sha256
    assert all(site.canonical_status == "exact_label_free" for site in first.metal_sites)
    assert all(site.canonical_status == "exact_label_free" for site in second.metal_sites)


def test_donor_pair_map_distinguishes_cis_and_trans_isomers():
    cis_entry = _entry(
        "FAKECIS",
        [
            ("Pt1", "Pt", (0, 0, 0)),
            ("N1", "N", (1, 0, 0)),
            ("N2", "N", (0, 1, 0)),
            ("Cl1", "Cl", (-1, 0, 0)),
            ("Cl2", "Cl", (0, -1, 0)),
        ],
        [(0, i, "Single") for i in range(1, 5)],
    )
    trans_entry = _entry(
        "FAKETRANS",
        [
            ("Pt1", "Pt", (0, 0, 0)),
            ("N1", "N", (1, 0, 0)),
            ("N2", "N", (-1, 0, 0)),
            ("Cl1", "Cl", (0, 1, 0)),
            ("Cl2", "Cl", (0, -1, 0)),
        ],
        [(0, i, "Single") for i in range(1, 5)],
    )
    cis = adapt_csd_entry(cis_entry, 19, source_release="test")
    trans = adapt_csd_entry(trans_entry, 20, source_release="test")
    assert len(cis.donor_pair_relations) == 6
    assert len(trans.donor_pair_relations) == 6
    groups_cis = {group.label: group.donor_elements for group in cis.donor_groups}
    groups_trans = {group.label: group.donor_elements for group in trans.donor_groups}
    cis_same_element_trans = [
        relation
        for relation in cis.donor_pair_relations
        if relation.angle_class == "trans"
        and groups_cis[relation.donor_group_a] == groups_cis[relation.donor_group_b]
    ]
    trans_same_element_trans = [
        relation
        for relation in trans.donor_pair_relations
        if relation.angle_class == "trans"
        and groups_trans[relation.donor_group_a] == groups_trans[relation.donor_group_b]
    ]
    assert not cis_same_element_trans
    assert len(trans_same_element_trans) == 2
    assert cis.record_checksum_sha256 != trans.record_checksum_sha256
    assert not validate_entry_record(cis)
    assert not validate_entry_record(trans)
