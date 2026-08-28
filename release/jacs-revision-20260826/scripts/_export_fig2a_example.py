#!/usr/bin/env python3
"""Export selected Fig 2A canonicalization example: GEFCON."""

import sys, os, json, copy
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from scipy.spatial.transform import Rotation
from coordrep.io.tmqm_reader import TMQMReader, Atom
from coordrep.encode import encode_molecule
from coordrep.identity.identity_keys import extract_identity_keys

TMQM_DIR = 'inputs/tmQM'
OUT_DIR = os.path.join(os.path.dirname(__file__), '..', 'revision_results', 'fig2a_example')
MOL_ID = 'ZENBID'


def rotate_mol(mol, seed=42):
    rng = np.random.RandomState(seed)
    R = Rotation.from_rotvec(rng.randn(3)).as_matrix()
    new_atoms = [Atom(a.index, a.element, *(R @ a.coords)) for a in mol.atoms]
    m = copy.deepcopy(mol)
    m.atoms = new_atoms
    m.mol_id = f"{mol.mol_id}_rot{seed}"
    return m


def permute_mol(mol, seed=42):
    rng = np.random.RandomState(seed)
    perm = rng.permutation(len(mol.atoms))
    new_atoms = [Atom(i, mol.atoms[perm[i]].element, mol.atoms[perm[i]].x,
                      mol.atoms[perm[i]].y, mol.atoms[perm[i]].z)
                 for i in range(len(mol.atoms))]
    m = copy.deepcopy(mol)
    m.atoms = new_atoms
    if mol.bond_orders is not None:
        n = len(mol.atoms)
        new_bo = np.zeros_like(mol.bond_orders)
        for i in range(n):
            for j in range(n):
                new_bo[np.where(perm == i)[0][0], np.where(perm == j)[0][0]] = mol.bond_orders[i, j]
        m.bond_orders = new_bo
    m.mol_id = f"{mol.mol_id}_perm{seed}"
    return m


def reverse_nonmetal(mol):
    from coordrep.io.tmqm_reader import TRANSITION_METALS
    metal_idx = next(i for i, a in enumerate(mol.atoms) if a.element in TRANSITION_METALS)
    others = [i for i in range(len(mol.atoms)) if i != metal_idx]
    others_rev = list(reversed(others))
    perm = list(range(len(mol.atoms)))
    for src, dst in zip(others, others_rev):
        perm[dst] = src
    new_atoms = [Atom(i, mol.atoms[perm[i]].element, mol.atoms[perm[i]].x,
                      mol.atoms[perm[i]].y, mol.atoms[perm[i]].z)
                 for i in range(len(mol.atoms))]
    m = copy.deepcopy(mol)
    m.atoms = new_atoms
    if mol.bond_orders is not None:
        n = len(mol.atoms)
        new_bo = np.zeros_like(mol.bond_orders)
        for i in range(n):
            for j in range(n):
                new_bo[i, j] = mol.bond_orders[perm[i], perm[j]]
        m.bond_orders = new_bo
    m.mol_id = f"{mol.mol_id}_ligtrav"
    return m


def write_xyz(atoms, mol_id, filepath):
    with open(filepath, 'w') as f:
        f.write(f"{len(atoms)}\n")
        f.write(f"CSD_code = {mol_id}\n")
        for a in atoms:
            f.write(f"{a.element:3s} {a.x:12.6f} {a.y:12.6f} {a.z:12.6f}\n")


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    reader = TMQMReader(TMQM_DIR)
    reader.load()
    mol = reader.get_molecule(MOL_ID)
    assert mol is not None, f"{MOL_ID} not found"

    # 1. Original
    cc = encode_molecule(mol)
    cc_can = cc.canonicalize()
    L0_orig = cc_can.to_string()
    keys = extract_identity_keys(L0_orig)

    # 2. Transformed variants
    mol_rot = rotate_mol(mol, seed=42)
    mol_perm = permute_mol(mol, seed=42)
    mol_lig = reverse_nonmetal(mol)

    cc_rot = encode_molecule(mol_rot).canonicalize()
    cc_perm = encode_molecule(mol_perm).canonicalize()
    cc_lig = encode_molecule(mol_lig).canonicalize()

    L0_rot = cc_rot.to_string()
    L0_perm = cc_perm.to_string()
    L0_lig = cc_lig.to_string()

    assert L0_rot == L0_orig, "Rotation invariance FAILED"
    assert L0_perm == L0_orig, "Permutation invariance FAILED"
    assert L0_lig == L0_orig, "Ligand traversal invariance FAILED"

    # 3. Write XYZ files for rendering
    write_xyz(mol.atoms, f"{MOL_ID}_original", os.path.join(OUT_DIR, f"{MOL_ID}_original.xyz"))
    write_xyz(mol_rot.atoms, f"{MOL_ID}_rotated", os.path.join(OUT_DIR, f"{MOL_ID}_rotated.xyz"))
    write_xyz(mol_perm.atoms, f"{MOL_ID}_permuted", os.path.join(OUT_DIR, f"{MOL_ID}_permuted.xyz"))

    # 4. Build simplified CoordRep block diagram
    block = []
    block.append("=== CoordRep Block Diagram ===")
    block.append("")
    block.append(f"Metal block:    {cc_can.metal.element}, CN={len(cc_can.graph.donor_indices)}")
    block.append(f"                [Metal:{cc_can.metal.element}|CN:{len(cc_can.graph.donor_indices)}]")
    block.append("")
    block.append(f"Shape block:    Best={keys.best_shape}, CShM=[{', '.join(f'{s}={v:.2f}' for s, v in zip(cc_can.shape.ref_shapes, cc_can.shape.values_rounded))}]")
    shape_str = L0_orig.split('<')[1].split('>')[0] if '<' in L0_orig else ''
    block.append(f"                <{shape_str}>")
    block.append("")
    block.append("Stereo block:")
    for tp in cc_can.constraints.trans_pairs:
        block.append(f"                {{trans:{tp[0]}--{tp[1]}}}")
    block.append("")
    block.append("Ligand dict:")
    for lig in cc_can.ligands:
        block.append(f"  {lig.lig_id}: dent={lig.dent}, donors={lig.donor_elements}, SMILES={lig.smiles}")
    block.append(f"                |{'|'.join(f'{lig.lig_id}={lig.smiles}' for lig in cc_can.ligands)}|")
    block_text = "\n".join(block)

    with open(os.path.join(OUT_DIR, "coordrep_block_diagram.txt"), 'w') as f:
        f.write(block_text)

    # 5. Write summary JSON
    summary = {
        "selected_example_id": MOL_ID,
        "source": "tmQM (public)",
        "metal": cc_can.metal.element,
        "cn": len(cc_can.graph.donor_indices),
        "n_atoms": mol.n_atoms,
        "n_ligands": len(cc_can.ligands),
        "denticity_list": [lig.dent for lig in cc_can.ligands],
        "donor_elements": [lig.donor_elements for lig in cc_can.ligands],
        "smiles": [lig.smiles for lig in cc_can.ligands],
        "best_shape": keys.best_shape,
        "shape_class": L0_orig.split('Class:')[1].split('|')[0] if 'Class:' in L0_orig else '?',
        "L0_StateKey": L0_orig,
        "L1_ShapeID": keys.L1_ShapeID,
        "L2_TopoID": keys.L2_TopoID,
        "L3_ConnID": keys.L3_ConnID,
        "L0_length": len(L0_orig),
        "invariance_verification": {
            "rotated_L0_match": L0_rot == L0_orig,
            "permuted_L0_match": L0_perm == L0_orig,
            "lig_traversal_L0_match": L0_lig == L0_orig,
        },
        "transformed_inputs": {
            "rotated": {"mol_id": mol_rot.mol_id, "L0": L0_rot},
            "permuted": {"mol_id": mol_perm.mol_id, "L0": L0_perm},
            "lig_traversal": {"mol_id": mol_lig.mol_id, "L0": L0_lig},
        },
        "files": [
            f"{MOL_ID}_original.xyz",
            f"{MOL_ID}_rotated.xyz",
            f"{MOL_ID}_permuted.xyz",
            "coordrep_block_diagram.txt",
            "fig2a_summary.json",
        ],
    }

    with open(os.path.join(OUT_DIR, "fig2a_summary.json"), 'w') as f:
        json.dump(summary, f, indent=2)

    # Print everything
    print("=" * 70)
    print("Fig 2A Selected Example")
    print("=" * 70)
    print(f"1. selected_example_id: {MOL_ID}")
    print(f"2. metal={cc_can.metal.element}, CN={len(cc_can.graph.donor_indices)}, denticity_list={[lig.dent for lig in cc_can.ligands]}")
    print(f"3. Full L0 CoordRep ({len(L0_orig)} chars):")
    print(f"   {L0_orig}")
    print()
    print("4. Transformed inputs:")
    print(f"   a) Rotated:      L0 match = {L0_rot == L0_orig}")
    print(f"   b) Permuted:     L0 match = {L0_perm == L0_orig}")
    print(f"   c) Lig traversal: L0 match = {L0_lig == L0_orig}")
    print()
    print("5. Invariance verification: ALL PASSED")
    print()
    print("6. Exported files:")
    for f_name in summary['files']:
        fpath = os.path.join(OUT_DIR, f_name)
        size = os.path.getsize(fpath)
        print(f"   {f_name:40s} ({size} bytes)")
    print()
    print("Block diagram:")
    print(block_text)
    print()
    print(f"Identity keys:")
    print(f"  L0 (StateKey): {keys.L0_StateKey[:80]}...")
    print(f"  L1 (ShapeID):  {keys.L1_ShapeID}")
    print(f"  L2 (TopoID):   {keys.L2_TopoID}")
    print(f"  L3 (ConnID):   {keys.L3_ConnID}")


if __name__ == '__main__':
    main()
