"""Regression fixtures for the narrowly scoped 1.1.1-candidate repairs."""

from __future__ import annotations

import numpy as np
import networkx as nx

from coordrep.core import (
    AssemblyGraph,
    ConstraintSet,
    CoordComplex,
    DonorSite,
    LigandModule,
    MetalState,
    ShapeVector,
)
from coordrep.identity.identity_keys import extract_identity_keys
from coordrep.io.tmqm_reader import Atom
from coordrep.serialize.parse_string import parse_constraint_block


def _site(ligand: str, element: str, rank: int) -> DonorSite:
    return DonorSite(ligand, element, rank)


def _ligand(
    lig_id: str,
    payload: str,
    provenance: str = "SMILES",
    connectivity_status: str = "supported_bond_order_graph",
) -> LigandModule:
    return LigandModule(
        lig_id=lig_id,
        smiles=payload,
        attach_atoms=[int(lig_id[1:]) - 1],
        donor_elements=["N"],
        dent=1,
        payload_provenance=provenance,
        connectivity_status=connectivity_status,
    )


def _complex(
    *,
    cn: int = 4,
    ligands: list[LigandModule] | None = None,
    constraints: ConstraintSet | None = None,
) -> CoordComplex:
    ligands = ligands or [_ligand("L1", "N"), _ligand("L2", "Cl")]
    donors = list(range(cn))
    refs = ["Td", "SP"] if cn == 4 else []
    values = np.array([0.4, 8.0]) if cn == 4 else np.array([])
    return CoordComplex(
        metal=MetalState("Pt", oxidation=2, dcount=8),
        ligands=ligands,
        graph=AssemblyGraph(
            metal_idx=cn,
            donor_indices=donors,
            lig_assignments={idx: ligands[idx % len(ligands)].lig_id for idx in donors},
            bond_orders={idx: 1.0 for idx in donors},
        ),
        shape=ShapeVector(
            cn=cn,
            ref_shapes=refs,
            values=values,
            values_rounded=values.copy(),
        ),
        constraints=constraints or ConstraintSet(),
    )


def test_cis_relations_are_serialized_and_losslessly_parsed():
    constraints = ConstraintSet(
        trans_pairs=[(_site("L1", "N", 1), _site("L2", "Cl", 1))],
        cis_pairs=[
            (_site("L1", "N", 1), _site("L3", "O", 1)),
            (_site("L2", "Cl", 1), _site("L3", "O", 1)),
        ],
        notes={"fac_mer": "fac"},
    ).canonicalize()
    cc = _complex(constraints=constraints)
    rendered = cc.to_string()
    parsed = parse_constraint_block(rendered)

    assert len(parsed.trans) == 1
    assert len(parsed.cis) == 2
    assert parsed.fac_mer == ("fac",)
    assert parsed.to_string() in rendered
    keys = extract_identity_keys(rendered)
    assert keys.constraint_signature == "T1_C2_fac"


def test_cn7_shape_is_explicitly_unsupported_and_validator_does_not_crash():
    cc = _complex(cn=7)
    rendered = cc.to_string()
    issues = cc.validate()
    keys = extract_identity_keys(rendered)

    assert (
        "<ShapeStatus:unsupported|Reason:no-reference-implementation|CN:7>"
        in rendered
    )
    assert any(issue.code == "SHAPE_REFERENCE_UNSUPPORTED" for issue in issues)
    assert keys.best_shape == "?"
    assert not keys.shape_supported
    assert keys.shape_status == "unsupported:no-reference-implementation"


def test_formula_payload_is_not_claimed_as_connectivity():
    ligands = [
        _ligand(
            "L1",
            "C2H3N",
            provenance="FORMULA",
            connectivity_status="unsupported_formula_fallback",
        )
    ]
    rendered = _complex(ligands=ligands).to_string()
    keys = extract_identity_keys(rendered)

    assert "|L1=FORMULA:C2H3N|" in rendered
    assert not keys.L3_connectivity_supported
    assert keys.L3_connectivity_status == "unsupported_formula_fallback"
    assert "CONNECTIVITY:UNSUPPORTED_FORMULA_FALLBACK" in keys.L3_ConnID


def test_formula_fallback_provenance_is_set_by_extractor(monkeypatch):
    import coordrep.graph.ligand_module as ligand_module

    monkeypatch.setattr(ligand_module, "RDKIT_AVAILABLE", False)
    atoms = [
        Atom(index=0, element="C", x=0.0, y=0.0, z=0.0),
        Atom(index=1, element="N", x=1.2, y=0.0, z=0.0),
    ]
    graph = nx.Graph()
    graph.add_nodes_from([0, 1])
    graph.add_edge(0, 1, bond_order=1.0)
    result = ligand_module.LigandExtractor()._generate_smiles(
        atoms,
        graph,
        {0, 1},
        has_bond_orders=False,
    )

    payload, _, provenance, status, issue = result
    assert payload == "CN"
    assert provenance == "FORMULA"
    assert status == "unsupported_formula_fallback"
    assert issue.code == "FORMULA_FALLBACK"


def test_smiles_payload_has_explicit_provenance():
    rendered = _complex(ligands=[_ligand("L1", "NCCN")]).to_string()
    keys = extract_identity_keys(rendered)

    assert "|L1=SMILES:NCCN|" in rendered
    assert keys.L3_connectivity_supported
    assert keys.ligand_provenance == ("SMILES",)


def test_legacy_untagged_string_remains_parseable_but_not_provenance_supported():
    legacy = (
        "[Metal:Pt|ox:+2|d:d8|CN:4]"
        "<ShapeBest:SP|Class:ideal|Delta:2|V:8.00,0.40>"
        "{cis:L1:N:1--L2:Cl:1}|L1=N|L2=Cl|"
    )
    keys = extract_identity_keys(legacy)

    assert keys.constraint_signature == "C1"
    assert keys.ligand_provenance == (
        "LEGACY_UNSPECIFIED",
        "LEGACY_UNSPECIFIED",
    )
    assert not keys.L3_connectivity_supported
    assert keys.L3_connectivity_status == "legacy_provenance_unspecified"
