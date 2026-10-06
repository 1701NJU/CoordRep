#!/usr/bin/env python3
"""Phase 1: extract donor-bearing non-metal CSD components.

Run with the CSD Python API interpreter.  The output deliberately retains
source atom labels only as a temporary mapping layer for the independent RDKit
phase; canonical payloads and donor keys are assessed in phase 2.
"""

from __future__ import annotations

import ast
import json
import sys
from collections import defaultdict
from pathlib import Path

from ccdc.io import EntryReader


ROOT = Path(__file__).resolve().parents[1]
OUTDIR = ROOT / "revision_results" / "csd_ligand_component_feasibility"

TRANSITION_METALS = {
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
}


def literal_assignment(path: Path, name: str):
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    for node in tree.body:
        if isinstance(node, ast.AnnAssign) and getattr(node.target, "id", None) == name:
            return ast.literal_eval(node.value)
        if isinstance(node, ast.Assign):
            if any(getattr(target, "id", None) == name for target in node.targets):
                return ast.literal_eval(node.value)
    raise RuntimeError(f"assignment {name!r} not found in {path}")


def is_pi_bond(bond) -> bool:
    value = str(bond.bond_type) if bond.bond_type is not None else ""
    return "pi" in value.lower() or "deloc" in value.lower()


def safe_component_smiles(component):
    try:
        value = component.smiles
    except Exception as exc:
        return "", f"{type(exc).__name__}:{str(exc)[:160]}"
    if not value:
        return "", "empty_smiles"
    return str(value), ""


def safe_mol2(component):
    try:
        return component.to_string("mol2"), ""
    except Exception as exc:
        return "", f"{type(exc).__name__}:{str(exc)[:160]}"


def audit_entry(reader, scope: str, refcode: str, expected_metals):
    try:
        entry = reader.entry(refcode)
    except Exception as exc:
        return [], {
            "scope": scope,
            "refcode": refcode,
            "status": "entry_lookup_failed",
            "reason": f"{type(exc).__name__}:{str(exc)[:160]}",
        }
    if entry is None:
        return [], {
            "scope": scope,
            "refcode": refcode,
            "status": "entry_not_found",
            "reason": "",
        }

    molecule = entry.molecule
    if molecule is None:
        return [], {
            "scope": scope,
            "refcode": refcode,
            "status": "no_molecule",
            "reason": "",
        }

    metals = [atom for atom in molecule.atoms if atom.atomic_symbol in TRANSITION_METALS]
    observed_metals = sorted(atom.atomic_symbol for atom in metals)
    expected_metals = sorted(expected_metals)
    metal_scope_match = observed_metals == expected_metals
    donor_records = {}
    for metal in metals:
        for bond in metal.bonds:
            other = bond.atoms[0] if bond.atoms[1] == metal else bond.atoms[1]
            if other.atomic_symbol in TRANSITION_METALS:
                continue
            label = str(other.label)
            record = donor_records.setdefault(label, {
                "label": label,
                "element": other.atomic_symbol,
                "target_metals": [],
                "pi": False,
            })
            if str(metal.label) not in record["target_metals"]:
                record["target_metals"].append(str(metal.label))
            record["pi"] = bool(record["pi"] or is_pi_bond(bond))

    stripped = molecule.copy()
    stripped_metals = [
        atom for atom in stripped.atoms if atom.atomic_symbol in TRANSITION_METALS
    ]
    stripped.remove_atoms(stripped_metals)
    components = list(stripped.components)

    label_to_components = defaultdict(list)
    for component_idx, component in enumerate(components):
        for atom in component.atoms:
            label_to_components[str(atom.label)].append(component_idx)

    missing_donors = sorted(
        label for label in donor_records if len(label_to_components.get(label, [])) != 1
    )

    component_rows = []
    donor_component_ids = set()
    for component_idx, component in enumerate(components):
        donors = [
            donor_records[label]
            for label in donor_records
            if len(label_to_components.get(label, [])) == 1
            and component_idx == label_to_components[label][0]
        ]
        if not donors:
            continue
        donor_component_ids.add(component_idx)
        smiles, smiles_error = safe_component_smiles(component)
        mol2, mol2_error = safe_mol2(component)
        component_rows.append({
            "scope": scope,
            "refcode": refcode,
            "component_index": component_idx,
            "metal_count": len(metals),
            "metals": observed_metals,
            "expected_metals": expected_metals,
            "metal_scope_match": metal_scope_match,
            "component_atom_labels": [str(atom.label) for atom in component.atoms],
            "component_heavy_atom_labels": [
                str(atom.label) for atom in component.atoms if atom.atomic_symbol != "H"
            ],
            "component_formula_elements": sorted(
                atom.atomic_symbol for atom in component.atoms
            ),
            "donors": donors,
            "n_donors": len(donors),
            "n_pi_donors": sum(bool(donor["pi"]) for donor in donors),
            "csd_smiles": smiles,
            "csd_smiles_error": smiles_error,
            "mol2": mol2,
            "mol2_error": mol2_error,
        })

    mapped_donors = sum(
        1 for label in donor_records if len(label_to_components.get(label, [])) == 1
    )
    entry_row = {
        "scope": scope,
        "refcode": refcode,
        "status": "processed" if not missing_donors else "donor_mapping_incomplete",
        "reason": "",
        "metal_count": len(metals),
        "observed_metals": observed_metals,
        "expected_metals": expected_metals,
        "metal_scope_match": metal_scope_match,
        "donor_atoms": len(donor_records),
        "mapped_donor_atoms": mapped_donors,
        "missing_or_ambiguous_donor_labels": missing_donors,
        "nonmetal_components": len(components),
        "donor_bearing_components": len(donor_component_ids),
        "component_ligand_merge_gain": max(0, len(donor_records) - len(donor_component_ids)),
    }
    return component_rows, entry_row


def main() -> int:
    multi_cases = literal_assignment(
        ROOT / "coordrep" / "v2beta" / "cases_multi.py", "MULTI_CASES_RAW"
    )
    haptic_cases = literal_assignment(
        ROOT / "coordrep" / "v2beta" / "cases_haptic.py", "HAPTIC_CASES_RAW"
    )
    requested = [
        ("multi", case[1], [spec[0] for spec in case[3]]) for case in multi_cases
    ] + [
        ("haptic", case[1], [case[3][0]]) for case in haptic_cases
    ]

    reader = EntryReader("CSD")
    component_rows = []
    entry_rows = []
    for idx, (scope, refcode, expected_metals) in enumerate(requested, 1):
        rows, entry_row = audit_entry(reader, scope, refcode, expected_metals)
        component_rows.extend(rows)
        entry_rows.append(entry_row)
        print(f"[{idx:02d}/{len(requested)}] {scope} {refcode}: {entry_row['status']}")

    OUTDIR.mkdir(parents=True, exist_ok=True)
    components_path = OUTDIR / "phase1_donor_components.jsonl"
    with components_path.open("w", encoding="utf-8") as handle:
        for row in component_rows:
            handle.write(json.dumps(row, ensure_ascii=False) + "\n")
    entries_path = OUTDIR / "phase1_entries.json"
    entries_path.write_text(
        json.dumps(entry_rows, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    by_scope = {}
    for scope in ("multi", "haptic"):
        subset = [row for row in entry_rows if row["scope"] == scope]
        processed = [row for row in subset if row["status"] in {"processed", "donor_mapping_incomplete"}]
        components_subset = [row for row in component_rows if row["scope"] == scope]
        by_scope[scope] = {
            "requested_entries": len(subset),
            "found_entries": len(processed),
            "missing_or_lookup_failed": len(subset) - len(processed),
            "donor_mapping_complete_entries": sum(
                row["status"] == "processed" for row in processed
            ),
            "metal_scope_matched_entries": sum(
                bool(row.get("metal_scope_match")) for row in processed
            ),
            "donor_atoms": sum(row.get("donor_atoms", 0) for row in processed),
            "mapped_donor_atoms": sum(row.get("mapped_donor_atoms", 0) for row in processed),
            "donor_bearing_components": len(components_subset),
            "csd_smiles_nonempty_components": sum(bool(row["csd_smiles"]) for row in components_subset),
            "mol2_nonempty_components": sum(bool(row["mol2"]) for row in components_subset),
            "component_ligand_merge_gain": sum(
                row.get("component_ligand_merge_gain", 0) for row in processed
            ),
        }
    summary = {
        "ccdc_api_version": "3.4.1",
        "csd_entry_count": len(reader),
        "definition": "remove all transition-metal atoms and treat each resulting non-metal connected component containing a donor as one ligand object",
        "by_scope": by_scope,
    }
    (OUTDIR / "phase1_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
