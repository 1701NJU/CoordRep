"""Red-team tests for simultaneous duplicate-ligand groups in the core."""

from __future__ import annotations

import itertools
import random

import pytest

from coordrep.core import (
    AssemblyGraph,
    ConstraintSet,
    CoordComplex,
    DonorSite,
    LigandModule,
    MetalState,
)


LIGAND_TYPES = {
    "A": ("N", "N", 0),
    "B": ("[Cl-]", "Cl", -1),
    "C": ("P", "P", 0),
}


def make_complex(instance_order, trans_edges, element="Co") -> CoordComplex:
    ligands = []
    donor_element = {}
    for source_id in instance_order:
        ligand_type = source_id[0]
        smiles, donor, charge = LIGAND_TYPES[ligand_type]
        donor_element[source_id] = donor
        ligands.append(LigandModule(
            lig_id=source_id,
            smiles=smiles,
            attach_atoms=[0],
            donor_elements=[donor],
            dent=1,
            attachment_set_key=f"[{donor}:1]",
            donor_attachment_keys=[f"[{donor}:2]"],
            charge=charge,
        ))

    trans_pairs = [
        (
            DonorSite(a, donor_element[a], 1),
            DonorSite(b, donor_element[b], 1),
        )
        for a, b in trans_edges
    ]
    return CoordComplex(
        metal=MetalState(element=element, oxidation=3),
        ligands=ligands,
        graph=AssemblyGraph(
            metal_idx=0,
            donor_indices=list(range(len(ligands))),
            lig_assignments={},
            bond_orders={},
        ),
        shape=None,
        constraints=ConstraintSet(trans_pairs=trans_pairs),
    )


def test_cn4_two_duplicate_groups_require_global_cartesian_minimum():
    # This is a direct counterexample to the former group-at-a-time search:
    # reversing both A and B simultaneously used to retain L4 in one variant.
    edges = [("A0", "B0")]
    outputs = set()
    for a_order in itertools.permutations(("A0", "A1")):
        for b_order in itertools.permutations(("B0", "B1")):
            cc = make_complex(a_order + b_order, edges, element="Pt")
            outputs.add(cc.canonicalize().to_string())
    assert len(outputs) == 1


def test_cn6_three_duplicate_groups_all_traversals_collapse():
    edges = [("A0", "B0"), ("A1", "C0"), ("B1", "C1")]
    outputs = set()
    for a_order in itertools.permutations(("A0", "A1")):
        for b_order in itertools.permutations(("B0", "B1")):
            for c_order in itertools.permutations(("C0", "C1")):
                cc = make_complex(a_order + b_order + c_order, edges)
                outputs.add(cc.canonicalize().to_string())
    assert len(outputs) == 1


def test_cn6_random_global_ligand_list_orders_and_source_labels_collapse():
    base_instances = ["A0", "A1", "B0", "B1", "C0", "C1"]
    base_edges = [("A0", "B0"), ("A1", "C0"), ("B1", "C1")]
    expected = make_complex(base_instances, base_edges).canonicalize().to_string()

    for seed in range(100):
        rng = random.Random(seed)
        order = list(base_instances)
        rng.shuffle(order)
        source_map = {
            instance: f"{instance[0]}_arbitrary_source_label_{seed}_{idx}"
            for idx, instance in enumerate(base_instances)
        }
        renamed_order = [source_map[instance] for instance in order]
        renamed_edges = [
            (source_map[a], source_map[b]) for a, b in base_edges
        ]
        # make_complex reads the first character as the ligand type, so the
        # arbitrary labels deliberately retain only that chemical type prefix.
        candidate = make_complex(renamed_order, renamed_edges).canonicalize().to_string()
        assert candidate == expected


def test_out_of_scope_factorial_search_fails_explicitly_not_approximately():
    instances = [f"A{idx}" for idx in range(9)]
    cc = make_complex(instances, [])
    with pytest.raises(ValueError, match="outside the validated CN<=6 core"):
        cc.canonicalize()


def test_internal_symmetric_donor_rank_flip_uses_declared_orbits():
    """Whole-record minimization resolves legal internal donor ties."""
    ligand_ids = ["E0", "E1", "E2"]
    base_relations = [
        ("E0", 1, "E1", 2),
        ("E0", 2, "E2", 1),
        ("E1", 1, "E2", 2),
    ]
    outputs = set()
    for flips in itertools.product((False, True), repeat=3):
        ligands = [
            LigandModule(
                lig_id=ligand_id,
                smiles="NCCN",
                attach_atoms=[0, 1],
                donor_elements=["N", "N"],
                dent=2,
                attachment_set_key="[NH2:1]CC[NH2:1]",
                donor_attachment_keys=["en-donor-orbit", "en-donor-orbit"],
                donor_rank_orbits=[[('N', 1), ('N', 2)]],
            )
            for ligand_id in ligand_ids
        ]
        relations = []
        for a, rank_a, b, rank_b in base_relations:
            if flips[ligand_ids.index(a)]:
                rank_a = 3 - rank_a
            if flips[ligand_ids.index(b)]:
                rank_b = 3 - rank_b
            relations.append((
                DonorSite(a, "N", rank_a),
                DonorSite(b, "N", rank_b),
            ))
        cc = CoordComplex(
            metal=MetalState(element="Co", oxidation=3),
            ligands=ligands,
            graph=AssemblyGraph(0, list(range(6)), {}, {}),
            shape=None,
            constraints=ConstraintSet(trans_pairs=relations),
        )
        outputs.add(cc.canonicalize().to_string())

    assert len(outputs) == 1
