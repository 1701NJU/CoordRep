#!/usr/bin/env python3
"""
coordrep_examples_gen.py
========================
Task 3: Generate 3 fully annotated CoordRep examples from real pipeline.

A. Monodentate-only square-planar (cisplatin-like)
B. Bidentate chelating (en-like or bpy-like)
C. Octahedral fac/mer with multiple stereo constraints

All strings come from the real encode → canonicalize → serialize pipeline.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
import textwrap
from pathlib import Path

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np
from coordrep.io.tmqm_reader import Atom, RawMolecule
from coordrep.encode import encode_molecule
from coordrep.core import CoordRepConfig
from coordrep.canonical.canonicalize import canonicalize_complex
from coordrep.serialize.to_string import serialize_complex
from coordrep.identity import extract_identity_keys
from coordrep_tools.validate import is_valid_coordrep


def build_molecule(atoms_data, bonds=None):
    """Build a RawMolecule from atoms list and optional bonds."""
    atoms = []
    for i, (elem, x, y, z) in enumerate(atoms_data):
        atoms.append(Atom(index=i, element=elem, x=x, y=y, z=z))
    n = len(atoms)
    bo = np.zeros((n, n), dtype=np.float32)
    if bonds:
        for i, j, order in bonds:
            bo[i, j] = order
            bo[j, i] = order
    return RawMolecule(mol_id="example", atoms=atoms, bond_orders=bo)


def annotate_coordrep(coordrep_str, cc):
    """Create annotation dict for a CoordRep string."""
    # Parse identity keys
    try:
        keys = extract_identity_keys(coordrep_str)
    except Exception:
        keys = None

    # Parse metal block
    metal_m = re.search(r'\[([^\]]+)\]', coordrep_str)
    metal_block = metal_m.group(0) if metal_m else ""

    # Parse shape block
    shape_m = re.search(r'<([^>]+)>', coordrep_str)
    shape_block = shape_m.group(0) if shape_m else ""

    # Parse constraint blocks
    constraint_blocks = re.findall(r'\{[^}]+\}', coordrep_str)

    # Parse ligand blocks
    lig_blocks = re.findall(r'\|(L\d+=(?:[^|]*))', coordrep_str)

    # Build ligand annotations
    ligand_annots = []
    for lig in cc.ligands:
        ligand_annots.append({
            'lig_id': lig.lig_id,
            'smiles': lig.smiles,
            'denticity': lig.dent,
            'donor_elements': lig.donor_elements,
            'eta': lig.eta,
        })

    return {
        'full_string': coordrep_str,
        'metal_block': metal_block,
        'shape_block': shape_block,
        'constraint_blocks': constraint_blocks,
        'ligand_blocks': lig_blocks,
        'ligand_annotations': ligand_annots,
        'L0_StateKey': keys.L0_StateKey if keys else "",
        'L1_ShapeID': keys.L1_ShapeID if keys else "",
        'L2_TopoID': keys.L2_TopoID if keys else "",
        'L3_ConnID': keys.L3_ConnID if keys else "",
        'best_shape': keys.best_shape if keys else "",
        'is_boundary': keys.binned_shape.is_boundary if keys and keys.binned_shape else False,
    }


def pretty_print(coordrep_str):
    """Split CoordRep string into readable multiline format."""
    lines = []
    # Metal block
    m = re.match(r'(\[[^\]]+\])', coordrep_str)
    if m:
        lines.append(m.group(1))
        rest = coordrep_str[m.end():]
    else:
        rest = coordrep_str

    # Shape block
    m = re.match(r'(<[^>]+>)', rest)
    if m:
        lines.append(m.group(1))
        rest = rest[m.end():]

    # Constraint blocks
    while True:
        m = re.match(r'(\{[^}]+\})', rest)
        if not m:
            break
        lines.append(m.group(1))
        rest = rest[m.end():]

    # Ligand blocks
    lig_parts = re.findall(r'\|([^|]+)', rest)
    for lp in lig_parts:
        if lp.strip():
            lines.append(f"|{lp}|")

    return '\n'.join(lines)


def find_in_pipeline(jsonl_path, criteria_fn, max_scan=60000):
    """Find a complex in the pipeline output matching criteria."""
    with open(jsonl_path) as f:
        for i, line in enumerate(f):
            if i >= max_scan:
                break
            d = json.loads(line)
            if criteria_fn(d):
                return d
    return None


def generate_example_A(config, pipeline_path):
    """
    Example A: Monodentate-only, square-planar CN=4.
    Look for a Pt or Pd complex with CN=4 and all dent=1.
    """
    print("  Searching for monodentate square-planar example …")

    def criteria(d):
        return (d.get('metal') in ('Pt', 'Pd') and
                d.get('cn') == 4 and
                d.get('ligand_dents') and
                all(dd == 1 for dd in d['ligand_dents']) and
                d.get('coordrep') and
                '{trans:' in d.get('coordrep', ''))

    rec = find_in_pipeline(pipeline_path, criteria)
    if rec:
        return {
            'label': 'A_monodentate_square_planar',
            'description': f"Monodentate-only square-planar {rec['metal']}(II) CN=4, source: tmQM {rec['mol_id']}",
            'mol_id': rec['mol_id'],
            'coordrep': rec['coordrep'],
            'source': 'tmQM_pipeline',
        }

    # Fallback: build cisplatin manually
    print("    Not found in pipeline, building cis-[PtCl2(NH3)2] from coordinates …")
    # Approximate square-planar cisplatin geometry
    atoms_data = [
        ('Pt', 0.0, 0.0, 0.0),
        ('N',  2.05, 0.0, 0.0),     # NH3 trans to Cl
        ('N',  0.0, 2.05, 0.0),     # NH3 trans to Cl
        ('Cl', -2.30, 0.0, 0.0),    # Cl trans to N
        ('Cl', 0.0, -2.30, 0.0),    # Cl trans to N
        # NH3 hydrogens
        ('H', 2.45, 0.95, 0.0),
        ('H', 2.45, -0.48, 0.82),
        ('H', 2.45, -0.48, -0.82),
        ('H', 0.0, 2.45, 0.95),
        ('H', -0.48, 2.45, -0.48),
        ('H', 0.48, 2.45, -0.48),
    ]
    bonds = [
        (0, 1, 1.0), (0, 2, 1.0), (0, 3, 1.0), (0, 4, 1.0),
        (1, 5, 1.0), (1, 6, 1.0), (1, 7, 1.0),
        (2, 8, 1.0), (2, 9, 1.0), (2, 10, 1.0),
    ]
    mol = build_molecule(atoms_data, bonds)
    mol.mol_id = "cisplatin_manual"
    return {
        'label': 'A_monodentate_square_planar',
        'description': 'cis-[PtCl2(NH3)2] (cisplatin), manually built coordinates',
        'mol_id': 'cisplatin_manual',
        'raw_molecule': mol,
        'source': 'manual_construction',
    }


def generate_example_B(config, pipeline_path):
    """
    Example B: Bidentate chelating ligand (en, bpy, etc).
    Look for a complex with at least one dent=2 and stereo constraints.
    """
    print("  Searching for bidentate chelating example …")

    def criteria(d):
        dents = d.get('ligand_dents', [])
        return (d.get('cn') in (4, 5, 6) and
                any(dd == 2 for dd in dents) and
                d.get('coordrep') and
                len(d.get('coordrep', '')) < 300 and
                '{trans:' in d.get('coordrep', ''))

    rec = find_in_pipeline(pipeline_path, criteria)
    if rec:
        return {
            'label': 'B_bidentate_chelating',
            'description': f"Bidentate chelating {rec['metal']} CN={rec['cn']}, "
                          f"dents={rec['ligand_dents']}, source: tmQM {rec['mol_id']}",
            'mol_id': rec['mol_id'],
            'coordrep': rec['coordrep'],
            'source': 'tmQM_pipeline',
        }
    return None


def generate_example_C(config, pipeline_path):
    """
    Example C: Octahedral CN=6 with multiple stereo constraints.
    """
    print("  Searching for octahedral fac/mer example …")

    def criteria(d):
        cr = d.get('coordrep', '')
        return (d.get('cn') == 6 and
                cr.count('{trans:') >= 2 and
                len(cr) < 400 and
                d.get('coordrep'))

    rec = find_in_pipeline(pipeline_path, criteria)
    if rec:
        return {
            'label': 'C_octahedral_multi_stereo',
            'description': f"Octahedral {rec['metal']} CN=6, "
                          f"dents={rec.get('ligand_dents', [])}, "
                          f"multiple trans constraints, source: tmQM {rec['mol_id']}",
            'mol_id': rec['mol_id'],
            'coordrep': rec['coordrep'],
            'source': 'tmQM_pipeline',
        }
    return None


def process_example(example, config):
    """
    Process one example: if it has a raw_molecule, run pipeline;
    otherwise parse the existing coordrep string.
    """
    if 'raw_molecule' in example:
        mol = example['raw_molecule']
        cc = encode_molecule(mol, config)
        cc = canonicalize_complex(cc)
        coordrep_str = cc.to_string()
        example['coordrep'] = coordrep_str
    else:
        coordrep_str = example['coordrep']
        # We need cc for annotations — create a minimal one by parsing
        cc = None

    annotation = {}
    annotation['label'] = example['label']
    annotation['description'] = example['description']
    annotation['mol_id'] = example['mol_id']
    annotation['source'] = example.get('source', '')
    annotation['full_coordrep'] = coordrep_str
    annotation['pretty_coordrep'] = pretty_print(coordrep_str)
    annotation['is_valid'] = is_valid_coordrep(coordrep_str, strict=True)

    # Parse identity keys
    try:
        keys = extract_identity_keys(coordrep_str)
        annotation['L0_StateKey'] = keys.L0_StateKey
        annotation['L1_ShapeID'] = keys.L1_ShapeID
        annotation['L2_TopoID'] = keys.L2_TopoID
        annotation['L3_ConnID'] = keys.L3_ConnID
        annotation['best_shape'] = keys.best_shape
    except Exception:
        annotation['L0_StateKey'] = ""
        annotation['L1_ShapeID'] = ""
        annotation['L2_TopoID'] = ""
        annotation['L3_ConnID'] = ""
        annotation['best_shape'] = ""

    # Parse blocks
    metal_m = re.search(r'\[([^\]]+)\]', coordrep_str)
    annotation['metal_block'] = metal_m.group(0) if metal_m else ""

    shape_m = re.search(r'<([^>]+)>', coordrep_str)
    annotation['shape_block'] = shape_m.group(0) if shape_m else ""

    annotation['constraint_blocks'] = re.findall(r'\{[^}]+\}', coordrep_str)
    annotation['ligand_blocks'] = re.findall(r'\|(L\d+=[^|]*)', coordrep_str)

    # Parse donor sites
    donor_refs = re.findall(r'(L\d+):([A-Za-z]+):(\d+)', coordrep_str)
    annotation['donor_sites'] = [f"{lid}:{e}:{r}" for lid, e, r in donor_refs]

    # Parse CN
    cn_m = re.search(r'CN:(\d+)', coordrep_str)
    annotation['cn'] = int(cn_m.group(1)) if cn_m else 0

    # Parse metal
    metal_m2 = re.search(r'Metal:([A-Z][a-z]?)', coordrep_str)
    annotation['metal'] = metal_m2.group(1) if metal_m2 else "?"

    return annotation


def write_maintext_md(examples, path):
    """Write main-text ready markdown with annotated examples."""
    with open(path, "w") as f:
        f.write("# Full Annotated CoordRep Examples\n\n")
        f.write("These examples demonstrate the complete CoordRep representation\n")
        f.write("for real coordination complexes at varying complexity levels.\n\n")

        for ex in examples:
            f.write(f"## {ex['label']}\n\n")
            f.write(f"**Description:** {ex['description']}\n\n")
            f.write(f"**Metal:** {ex['metal']}  |  **CN:** {ex['cn']}  |  "
                    f"**Best shape:** {ex.get('best_shape', '?')}\n\n")
            f.write("**Full CoordRep string:**\n\n")
            f.write(f"```\n{ex['full_coordrep']}\n```\n\n")
            f.write("**Pretty-printed (multiline):**\n\n")
            f.write(f"```\n{ex['pretty_coordrep']}\n```\n\n")
            f.write("**Annotation table:**\n\n")
            f.write("| Field | Value |\n|-------|-------|\n")
            f.write(f"| Metal block | `{ex['metal_block']}` |\n")
            f.write(f"| Shape block | `{ex['shape_block']}` |\n")
            for i, cb in enumerate(ex.get('constraint_blocks', [])):
                f.write(f"| Stereo constraint {i+1} | `{cb}` |\n")
            for lb in ex.get('ligand_blocks', []):
                lid = lb.split('=')[0]
                smi = lb.split('=', 1)[1] if '=' in lb else ''
                f.write(f"| Ligand {lid} | `{smi}` |\n")
            f.write(f"| Donor sites | {', '.join(ex.get('donor_sites', []))} |\n")
            f.write(f"| Valid | {ex['is_valid']} |\n")
            f.write("\n**Identity keys:**\n\n")
            f.write(f"- **L0 (StateKey):** `{ex['L0_StateKey'][:80]}…`\n")
            f.write(f"- **L1 (ShapeID):** `{ex['L1_ShapeID'][:80]}…`\n")
            f.write(f"- **L2 (TopoID):** `{ex['L2_TopoID'][:80]}…`\n")
            f.write(f"- **L3 (ConnID):** `{ex['L3_ConnID'][:80]}…`\n")
            f.write("\n---\n\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pipeline",
                        default="inputs/tmqm/results.jsonl")
    parser.add_argument("--out", default="revision_results/examples")
    args = parser.parse_args()

    os.makedirs(args.out, exist_ok=True)
    config = CoordRepConfig.default()

    # Generate 3 examples
    examples_raw = []

    ex_a = generate_example_A(config, args.pipeline)
    if ex_a:
        examples_raw.append(ex_a)

    ex_b = generate_example_B(config, args.pipeline)
    if ex_b:
        examples_raw.append(ex_b)

    ex_c = generate_example_C(config, args.pipeline)
    if ex_c:
        examples_raw.append(ex_c)

    print(f"  Generated {len(examples_raw)} raw examples")

    # Process through pipeline
    examples = []
    for ex in examples_raw:
        ann = process_example(ex, config)
        examples.append(ann)
        print(f"  {ann['label']}: valid={ann['is_valid']}, cn={ann['cn']}, "
              f"metal={ann['metal']}, shape={ann.get('best_shape', '?')}")

    # ── Write outputs ─────────────────────────────────────
    # JSON (full data)
    p = os.path.join(args.out, "coordrep_examples_full.json")
    with open(p, "w") as f:
        json.dump(examples, f, indent=2)
    print(f"  {p}")

    # Main-text markdown
    p = os.path.join(args.out, "coordrep_examples_maintext.md")
    write_maintext_md(examples, p)
    print(f"  {p}")

    # SI markdown (same content for now)
    p = os.path.join(args.out, "coordrep_examples_si.md")
    write_maintext_md(examples, p)
    print(f"  {p}")

    # Plain strings
    p = os.path.join(args.out, "coordrep_example_strings.txt")
    with open(p, "w") as f:
        for ex in examples:
            f.write(f"# {ex['label']} ({ex['description']})\n")
            f.write(f"{ex['full_coordrep']}\n\n")
    print(f"  {p}")

    print(f"\n{'='*60}")
    print("COORDREP EXAMPLES GENERATED")
    print(f"{'='*60}")
    for ex in examples:
        print(f"\n  [{ex['label']}]")
        print(f"  {ex['description']}")
        print(f"  CoordRep: {ex['full_coordrep'][:100]}…")
    print()


if __name__ == "__main__":
    main()
