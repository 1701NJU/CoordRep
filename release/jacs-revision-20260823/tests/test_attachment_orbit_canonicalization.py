"""Regression tests for ligand-local attachment and donor-orbit identity."""

from __future__ import annotations

import random
from pathlib import Path

import networkx as nx
import numpy as np
import pytest

from coordrep.canonical.canonicalize import _ligand_equivalence_key
from coordrep.core import (
    AssemblyGraph,
    ConstraintSet,
    CoordComplex,
    DonorSite,
    LigandModule,
    MetalState,
)
from coordrep.encode import _load_molecule, encode_molecule
from coordrep.graph.ligand_module import LigandExtractor
from coordrep.graph.donor_sites import assign_donor_ranks
from coordrep.io.tmqm_reader import Atom, RawMolecule


FIXTURE_DIR = Path(__file__).parent / "fixtures" / "xyz"
ABAZAH_FIXTURE = FIXTURE_DIR / "ABAZAH.xyz"
TOHTIV_FIXTURE = FIXTURE_DIR / "TOHTIV.xyz"


def _complex(ligands, trans_pairs=(), cn=None):
    if cn is None:
        cn = sum(lig.dent for lig in ligands)
    return CoordComplex(
        metal=MetalState("Ir", oxidation=3),
        ligands=list(ligands),
        graph=AssemblyGraph(
            metal_idx=0,
            donor_indices=list(range(cn)),
            lig_assignments={},
            bond_orders={},
        ),
        shape=None,
        constraints=ConstraintSet(trans_pairs=list(trans_pairs)),
    )


def _monodentate(lig_id, attach_atom, payload, element, attachment_key):
    return LigandModule(
        lig_id=lig_id,
        smiles=payload,
        attach_atoms=[attach_atom],
        donor_elements=[element],
        dent=1,
        attachment_set_key=attachment_key,
        donor_attachment_keys=[f"{attachment_key}|target"],
    )


def test_equivalent_ligands_ignore_distinct_global_attach_indices():
    left = _monodentate("source_left", 7, "N", "N", "[NH3:1]")
    right = _monodentate("source_right", 913, "N", "N", "[NH3:1]")
    chloride = _monodentate("source_cl", 44, "[Cl-]", "Cl", "[Cl-:1]")

    assert _ligand_equivalence_key(left) == _ligand_equivalence_key(right)

    first = _complex(
        [left, right, chloride],
        [(DonorSite("source_left", "N", 1), DonorSite("source_cl", "Cl", 1))],
    )
    second = _complex(
        [
            _monodentate("arbitrary_b", 3, "N", "N", "[NH3:1]"),
            _monodentate("arbitrary_cl", 1001, "[Cl-]", "Cl", "[Cl-:1]"),
            _monodentate("arbitrary_a", 81, "N", "N", "[NH3:1]"),
        ],
        [(DonorSite("arbitrary_a", "N", 1), DonorSite("arbitrary_cl", "Cl", 1))],
    )

    assert first.canonicalize().to_string() == second.canonicalize().to_string()


def _cyclohexane_attachment_descriptors(donors):
    atoms = [
        Atom(index=i, element="C", x=float(i), y=0.0, z=0.0)
        for i in range(6)
    ]
    graph = nx.Graph()
    graph.add_nodes_from((i, {"element": "C"}) for i in range(6))
    for i in range(6):
        graph.add_edge(i, (i + 1) % 6, bond_order=1.0)
    return LigandExtractor()._generate_attachment_descriptors(
        atoms, graph, set(range(6)), donors
    )


def test_same_payload_nonorbit_attachment_sets_do_not_merge_or_serialize_equal():
    adjacent_key, adjacent_donors = _cyclohexane_attachment_descriptors([0, 1])
    opposite_key, opposite_donors = _cyclohexane_attachment_descriptors([0, 3])
    assert adjacent_key != opposite_key

    adjacent = LigandModule(
        lig_id="Lsrc",
        smiles="C1CCCCC1",
        attach_atoms=[100, 101],
        donor_elements=["C", "C"],
        dent=2,
        attachment_set_key=adjacent_key,
        donor_attachment_keys=[adjacent_donors[0], adjacent_donors[1]],
    )
    opposite = LigandModule(
        lig_id="Lsrc",
        smiles="C1CCCCC1",
        attach_atoms=[400, 900],
        donor_elements=["C", "C"],
        dent=2,
        attachment_set_key=opposite_key,
        donor_attachment_keys=[opposite_donors[0], opposite_donors[3]],
    )

    assert _ligand_equivalence_key(adjacent) != _ligand_equivalence_key(opposite)
    assert (
        _complex([adjacent]).canonicalize().to_string()
        != _complex([opposite]).canonicalize().to_string()
    )


def test_repeated_payload_without_attachment_identity_fails_closed():
    unresolved = [
        LigandModule(
            lig_id=f"U{idx}",
            smiles="c1ccncc1",
            attach_atoms=[10 + idx],
            donor_elements=["N"],
            dent=1,
        )
        for idx in range(2)
    ]
    with pytest.raises(ValueError, match="attachment_set_key"):
        _complex(unresolved).canonicalize()


def test_internal_symmetric_donor_tie_collapses_by_whole_record_minimum():
    ligand_ids = ["E0", "E1", "E2"]
    base_relations = [
        ("E0", 1, "E1", 2),
        ("E0", 2, "E2", 1),
        ("E1", 1, "E2", 2),
    ]
    outputs = set()

    for mask in range(8):
        ligands = [
            LigandModule(
                lig_id=lig_id,
                smiles="NCCN",
                attach_atoms=[50 + 2 * idx, 51 + 2 * idx],
                donor_elements=["N", "N"],
                dent=2,
                attachment_set_key="[NH2:1]CC[NH2:1]",
                donor_attachment_keys=["en-orbit", "en-orbit"],
                donor_rank_orbits=[[('N', 1), ('N', 2)]],
            )
            for idx, lig_id in enumerate(ligand_ids)
        ]
        relations = []
        for a, rank_a, b, rank_b in base_relations:
            if mask & (1 << ligand_ids.index(a)):
                rank_a = 3 - rank_a
            if mask & (1 << ligand_ids.index(b)):
                rank_b = 3 - rank_b
            relations.append((
                DonorSite(a, "N", rank_a),
                DonorSite(b, "N", rank_b),
            ))
        outputs.add(_complex(ligands, relations).canonicalize().to_string())

    assert len(outputs) == 1


def test_adapter_declares_only_same_orbit_same_geometry_donor_ties():
    atoms = [
        Atom(0, "Co", 0.0, 0.0, 0.0),
        Atom(1, "N", -2.0, 0.0, 0.0),
        Atom(2, "N", 2.0, 0.0, 0.0),
        Atom(3, "C", -3.2, 0.0, 0.0),
        Atom(4, "C", 3.2, 0.0, 0.0),
    ]
    bond_orders = np.zeros((5, 5))
    for a, b in ((0, 1), (0, 2), (1, 3), (2, 4)):
        bond_orders[a, b] = bond_orders[b, a] = 1.0
    molecule = RawMolecule("SYMMETRIC_TIE", atoms, bond_orders=bond_orders)
    ligand = LigandModule(
        lig_id="E0",
        smiles="NCCN",
        attach_atoms=[1, 2],
        donor_elements=["N", "N"],
        dent=2,
        attachment_set_key="[NH2:1]CC[NH2:1]",
        donor_attachment_keys=["en-orbit", "en-orbit"],
    )

    sites = assign_donor_ranks(molecule, ligand, metal_idx=0)

    assert {site.donor_rank for site in sites.values()} == {1, 2}
    assert ligand.donor_rank_orbits == [[("N", 1), ("N", 2)]]


def _renumber_raw_molecule(mol, permutation):
    atoms = [
        Atom(
            index=new_index,
            element=mol.atoms[old_index].element,
            x=mol.atoms[old_index].x,
            y=mol.atoms[old_index].y,
            z=mol.atoms[old_index].z,
        )
        for new_index, old_index in enumerate(permutation)
    ]
    bond_orders = None
    if mol.bond_orders is not None:
        bond_orders = mol.bond_orders[permutation][:, permutation]
    return RawMolecule(
        mol_id=mol.mol_id,
        atoms=atoms,
        bond_orders=bond_orders,
        properties=mol.properties,
    )


@pytest.mark.skipif(
    not ABAZAH_FIXTURE.is_file(),
    reason="public release omits coordinate fixtures; hash and retrieval provenance are documented",
)
def test_abazah_raw_molecule_atom_permutations_are_bit_identical():
    source = _load_molecule(str(ABAZAH_FIXTURE))
    outputs = set()
    variant_count = 16

    for seed in range(variant_count):
        permutation = list(range(source.n_atoms))
        random.Random(seed).shuffle(permutation)
        variant = _renumber_raw_molecule(source, permutation)
        encoded = encode_molecule(variant)
        assert len(encoded.graph.donor_indices) == 6
        assert len(encoded.ligands) == 3
        outputs.add(encoded.canonicalize().to_string())

    assert len(outputs) == 1


@pytest.mark.skipif(
    not TOHTIV_FIXTURE.is_file(),
    reason="public release omits coordinate fixtures; hash and retrieval provenance are documented",
)
def test_tohtiv_nonidentical_equal_mass_ligands_have_stable_order():
    """Equal-mass nonidentical ligands must not be ordered by float MW."""
    source = _load_molecule(str(TOHTIV_FIXTURE))
    outputs = set()

    for seed in range(20):
        permutation = list(range(source.n_atoms))
        random.Random(seed).shuffle(permutation)
        variant = _renumber_raw_molecule(source, permutation)
        encoded = encode_molecule(variant)
        assert len(encoded.graph.donor_indices) == 5
        assert len(encoded.ligands) == 3
        outputs.add(encoded.canonicalize().to_string())

    assert len(outputs) == 1
