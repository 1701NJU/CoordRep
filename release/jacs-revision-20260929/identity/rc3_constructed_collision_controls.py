"""Run constructed stereoisomer and chelate-pointer controls through public rc3.

    python release/jacs-revision-20260929/identity/rc3_constructed_collision_controls.py

This exercises the published 1.1.2rc3 shape, relation, exact canonicalizer,
serializer, and L0--L3 identity-key code. The constructed shells are not
corpus records and make no claim about collision prevalence.
"""

from __future__ import annotations

import csv
import hashlib
import json
import re
import sys
from collections import Counter
from copy import deepcopy
from pathlib import Path

import numpy as np
from rdkit import Chem, rdBase


ROOT = Path(__file__).resolve().parents[3]
RELEASE = ROOT / "release/jacs-revision-20260826"
assert RELEASE.is_dir(), f"Published rc3 copy not found: {RELEASE}"
assert 'version = "1.1.2rc3"' in (RELEASE / "pyproject.toml").read_text()
sys.path.insert(0, str(RELEASE))

from coordrep.core import (  # noqa: E402
    AssemblyGraph, CoordComplex, DonorSite, LigandModule, MetalState, ShapeVector,
)
from coordrep.geometry.rel_config import ConstraintDetector  # noqa: E402
from coordrep.geometry.shape import ShapeCalculator  # noqa: E402
from coordrep.identity.identity_keys import extract_identity_keys_from_cc  # noqa: E402
from coordrep.io.tmqm_reader import Atom, RawMolecule  # noqa: E402
from coordrep.serialize.parse_string import parse_constraint_block  # noqa: E402


SP = np.asarray([[1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]], float)
OH = np.asarray(
    [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]],
    float,
)
CHELATE_OH = np.asarray(
    [[1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]],
    float,
)
SPECS = {
    "Pt_cis": ("Pt", 2, 8, SP, ("Cl", "Cl", "N", "N")),
    "Pt_trans": ("Pt", 2, 8, SP, ("Cl", "N", "Cl", "N")),
    "Co_fac": ("Co", 3, 6, OH, ("Cl", "N", "Cl", "N", "Cl", "N")),
    "Co_mer": ("Co", 3, 6, OH, ("Cl", "Cl", "Cl", "N", "N", "N")),
}
PAIR_RE = re.compile(r"L\d+:([A-Za-z]+):\d+--L\d+:([A-Za-z]+):\d+")


def sha(data: str | bytes) -> str:
    if isinstance(data, str):
        data = data.encode("utf-8")
    return hashlib.sha256(data).hexdigest()


def attachment_descriptors(smiles: str, donor_atom_indices: list[int]):
    """Mirror rc3's atom-map-1 set and atom-map-2 singled-donor semantics."""
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None and Chem.MolToSmiles(mol) == smiles

    def coloured(target: int | None) -> str:
        copied = Chem.Mol(mol)
        for idx in donor_atom_indices:
            copied.GetAtomWithIdx(idx).SetAtomMapNum(1)
        if target is not None:
            copied.GetAtomWithIdx(target).SetAtomMapNum(2)
        return Chem.MolToSmiles(copied, canonical=True, isomericSmiles=True)

    return coloured(None), [coloured(idx) for idx in donor_atom_indices]


def ligand(lig_id: str, smiles: str, indices: list[int], elements: list[str], charge: int, donor_atoms: list[int]):
    attachment_set_key, donor_attachment_keys = attachment_descriptors(smiles, donor_atoms)
    return LigandModule(
        lig_id=lig_id,
        smiles=smiles,
        attach_atoms=indices,
        donor_elements=elements,
        dent=len(indices),
        charge=charge,
        payload_provenance="SMILES",
        connectivity_status="nominal_constructed_smiles",
        attachment_set_key=attachment_set_key,
        donor_attachment_keys=donor_attachment_keys,
    )


def build(name: str):
    if name == "Co_glycinato_tetraammine":
        metal, oxidation, dcount, geometry = "Co", 3, 6, CHELATE_OH
        elements = ("N", "O", "N", "N", "N", "N")
        gly = Chem.MolFromSmiles("NCC(=O)[O-]")
        assert gly is not None
        donor_atoms = [
            next(a.GetIdx() for a in gly.GetAtoms() if a.GetSymbol() == "N"),
            next(a.GetIdx() for a in gly.GetAtoms() if a.GetSymbol() == "O" and a.GetFormalCharge() == -1),
        ]
        ligands = [ligand("L1", "NCC(=O)[O-]", [1, 2], ["N", "O"], -1, donor_atoms)]
        ligands.extend(
            ligand(f"L{i - 1}", "N", [i], ["N"], 0, [0])
            for i in range(3, 7)
        )
        assignment = {1: "L1", 2: "L1", **{i: f"L{i - 1}" for i in range(3, 7)}}
        sites = {1: DonorSite("L1", "N", 1), 2: DonorSite("L1", "O", 1)}
        sites.update({i: DonorSite(f"L{i - 1}", "N", 1) for i in range(3, 7)})
    else:
        metal, oxidation, dcount, geometry, elements = SPECS[name]
        ligands = [
            ligand(
                f"L{i}", "[Cl-]" if label == "Cl" else "N", [i],
                [label], -1 if label == "Cl" else 0, [0],
            )
            for i, label in enumerate(elements, start=1)
        ]
        assignment = {i: f"L{i}" for i in range(1, len(elements) + 1)}
        sites = {i: DonorSite(f"L{i}", label, 1) for i, label in enumerate(elements, start=1)}
    cn = len(elements)
    atoms = [Atom(0, metal, 0, 0, 0)] + [
        Atom(i, el, *map(float, coord))
        for i, (el, coord) in enumerate(zip(elements, geometry), start=1)
    ]
    donor_indices = list(range(1, cn + 1))
    constraints = ConstraintDetector().detect_constraints_with_sites(
        RawMolecule(name, atoms), 0, donor_indices, ligands, sites
    )
    shape = ShapeCalculator(rounding_decimals=2).compute(geometry, np.zeros(3))
    assert shape.best_shape == ("SP" if cn == 4 else "Oh")
    cc = CoordComplex(
        metal=MetalState(metal, oxidation=oxidation, dcount=dcount),
        ligands=ligands,
        graph=AssemblyGraph(
            metal_idx=0, donor_indices=donor_indices,
            lig_assignments=assignment,
            bond_orders={i: 1.0 for i in donor_indices},
        ),
        shape=ShapeVector(cn=cn, ref_shapes=shape.ref_shapes, values=shape.values.copy()),
        constraints=constraints,
        source_id=f"constructed_{name}",
    )
    return cc, geometry, elements


def run_case(name: str) -> dict:
    cc, geometry, elements = build(name)
    keys = extract_identity_keys_from_cc(cc)
    assert keys.L0_StateKey == cc.canonicalize().to_string()
    reversed_cc = deepcopy(cc)
    reversed_cc.ligands.reverse()
    reversed_keys = extract_identity_keys_from_cc(reversed_cc)
    assert keys.L0_StateKey == reversed_keys.L0_StateKey, name
    parsed = parse_constraint_block(keys.L0_StateKey)
    assert len(parsed.trans) == len(elements) // 2
    assert len(parsed.trans) + len(parsed.cis) == len(elements) * (len(elements) - 1) // 2
    trans_element_counts = Counter()
    for relation in parsed.trans:
        m = PAIR_RE.fullmatch(relation)
        assert m, relation
        trans_element_counts["-".join(sorted(m.groups()))] += 1
    opposite_cl_pairs = sum(
        np.isclose(np.dot(geometry[i], geometry[j]), -1)
        for i in range(len(elements)) for j in range(i + 1, len(elements))
        if elements[i] == elements[j] == "Cl"
    )
    if name in SPECS:
        assert int(opposite_cl_pairs) == (0 if name.endswith(("cis", "fac")) else 1)
        assert int(trans_element_counts.get("Cl-Cl", 0)) == int(opposite_cl_pairs)
    result = {
        "name": name,
        "metal": cc.metal.element,
        "oxidation": cc.metal.oxidation,
        "cn": len(elements),
        "donor_elements_by_position": list(elements),
        "donor_vectors_angstrom": geometry.tolist(),
        "nominal_smiles_multiset": sorted(l.smiles for l in cc.ligands),
        "ligand_attachment_set_keys": {l.lig_id: l.attachment_set_key for l in cc.ligands},
        "ligand_donor_attachment_keys": {l.lig_id: l.donor_attachment_keys for l in cc.ligands},
        "opposite_Cl_Cl_pairs": int(opposite_cl_pairs),
        "trans_element_composition": dict(sorted(trans_element_counts.items())),
        "trans_relations": list(parsed.trans),
        "cis_relations": list(parsed.cis),
        "shape_reference_order": list(cc.shape.ref_shapes),
        "shape_values_unrounded": [float(v) for v in cc.shape.values],
        "shape_values_rounded": [round(float(v), 2) for v in cc.shape.values],
        "shape_bin": keys.binned_shape.token,
        "constraint_signature": keys.constraint_signature,
        "L0_StateKey": keys.L0_StateKey,
        "L0_sha256": sha(keys.L0_StateKey),
        "L1_ShapeID": keys.L1_ShapeID,
        "L1_sha256": sha(keys.L1_ShapeID),
        "L2_TopoID": keys.L2_TopoID,
        "L3_ConnID": keys.L3_ConnID,
        "reversed_ligand_traversal_L0_equal": keys.L0_StateKey == reversed_keys.L0_StateKey,
    }
    if name == "Co_glycinato_tetraammine":
        # L0 retains both pointers for the two donor atoms of one ligand.
        relation = "{cis:L1:N:1--L1:O:1}"
        assert relation in keys.L0_StateKey
        assert "L1:N:1" not in keys.L1_ShapeID and "L1:O:1" not in keys.L1_ShapeID
        result["same_ligand_cis_relation_in_L0"] = relation
        result["dropped_in_L1"] = ["L1:N:1", "L1:O:1", relation]
    return result


def main() -> None:
    names = ["Pt_cis", "Pt_trans", "Co_fac", "Co_mer", "Co_glycinato_tetraammine"]
    cases = {name: run_case(name) for name in names}
    pairs = {}
    for label, a, b in (("cis_trans", "Pt_cis", "Pt_trans"), ("fac_mer", "Co_fac", "Co_mer")):
        left, right = cases[a], cases[b]
        fixed = {
            field: left[field] == right[field]
            for field in (
                "metal", "oxidation", "cn", "donor_vectors_angstrom",
                "nominal_smiles_multiset", "shape_reference_order",
                "shape_values_rounded", "shape_bin",
            )
        }
        assert all(fixed.values()), (label, fixed)
        equality = {
            layer: left[layer] == right[layer]
            for layer in ("L0_StateKey", "L1_ShapeID", "L2_TopoID", "L3_ConnID")
        }
        assert equality == {"L0_StateKey": False, "L1_ShapeID": True, "L2_TopoID": True, "L3_ConnID": True}
        assert left["trans_element_composition"] != right["trans_element_composition"]
        pairs[label] = {
            "cases": [a, b], "fixed_fields": fixed,
            "layer_key_equal": equality,
            "L1_collision": True,
            "trans_element_composition_differs": True,
        }
    source_files = (
        "coordrep/core.py", "coordrep/geometry/shape.py", "coordrep/geometry/rel_config.py",
        "coordrep/canonical/canonicalize.py", "coordrep/serialize/to_string.py",
        "coordrep/identity/identity_keys.py",
    )
    report = {
        "package_version": "1.1.2rc3",
        "package_path": str(RELEASE.relative_to(ROOT)).replace("\\", "/"),
        "python_version": sys.version.split()[0],
        "rdkit_version": rdBase.rdkitVersion,
        "numpy_version": np.__version__,
        "source_sha256": {p: sha((RELEASE / p).read_bytes()) for p in source_files},
        "cases": cases,
        "pairs": pairs,
        "limitations": [
            "Constructed ideal donor shells with manually declared canonical ligand SMILES and ligand-local atom-map descriptors; not full structural ingestion.",
            "Runs the published rc3 production shape, relation, exact canonicalization, serialization, and L0-L3 key path.",
            "These pairwise controls establish the stated L1 collisions only, not a corpus collision frequency.",
            "For monodentate M(A)3(B)3, fac/mer differs in pairwise relations; the detector's explicit fm marker is limited to a denticity-three ligand.",
        ],
    }
    out = Path(__file__).with_name("rc3_constructed_collision_controls.json")
    out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    csv_out = Path(__file__).with_name("rc3_constructed_collision_controls.csv")
    with csv_out.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=(
            "case", "metal", "oxidation", "cn", "shape_bin", "constraint_signature",
            "trans_element_composition", "L0_StateKey", "L0_sha256", "L1_ShapeID", "L1_sha256",
            "reversed_ligand_traversal_L0_equal",
        ))
        writer.writeheader()
        for case in cases.values():
            writer.writerow({
                "case": case["name"],
                **{key: case[key] for key in writer.fieldnames if key in case},
                "trans_element_composition": json.dumps(case["trans_element_composition"], sort_keys=True),
            })
    print(f"Wrote {out} and {csv_out}")
    for label, pair in pairs.items():
        print(label, pair["layer_key_equal"])
    print("Chelate same-ligand L0 relation:", cases["Co_glycinato_tetraammine"]["same_ligand_cis_relation_in_L0"])
    print("Chelate L1:", cases["Co_glycinato_tetraammine"]["L1_ShapeID"])


if __name__ == "__main__":
    main()
