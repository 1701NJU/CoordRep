#!/usr/bin/env python3
"""Phase 2: RDKit payload/atom-map and permutation feasibility audit."""

from __future__ import annotations

import json
import random
from collections import Counter, defaultdict
from pathlib import Path

from rdkit import Chem
from rdkit import RDLogger


ROOT = Path(__file__).resolve().parents[1]
OUTDIR = ROOT / "revision_results" / "csd_ligand_component_feasibility"
N_PERMUTATIONS = 20


def atom_name(atom) -> str:
    for prop in ("_TriposAtomName", "molFileAlias"):
        if atom.HasProp(prop):
            return atom.GetProp(prop)
    return ""


def donor_key_multiset(molecule, donor_labels):
    ranks = list(Chem.CanonicalRankAtoms(
        molecule,
        breakTies=False,
        includeChirality=True,
        includeIsotopes=True,
    ))
    names = {atom_name(atom): atom.GetIdx() for atom in molecule.GetAtoms()}
    if any(label not in names for label in donor_labels):
        return None
    return tuple(sorted(
        (
            molecule.GetAtomWithIdx(names[label]).GetSymbol(),
            int(ranks[names[label]]),
        )
        for label in donor_labels
    ))


def audit_component(row):
    result = {
        "scope": row["scope"],
        "refcode": row["refcode"],
        "component_index": row["component_index"],
        "metal_scope_match": bool(row.get("metal_scope_match")),
        "csd_smiles": row["csd_smiles"],
        "n_donors": row["n_donors"],
        "n_pi_donors": row["n_pi_donors"],
        "status": "",
        "canonical_smiles": "",
        "donor_key_multiset": None,
        "permutation_trials": 0,
        "permutation_passed": 0,
    }
    if not row["mol2"]:
        result["status"] = "mol2_missing"
        return result

    molecule = Chem.MolFromMol2Block(
        row["mol2"], sanitize=True, removeHs=False, cleanupSubstructures=True
    )
    if molecule is None:
        result["status"] = "rdkit_mol2_parse_or_sanitize_failed"
        return result

    donor_labels = [str(donor["label"]) for donor in row["donors"]]
    names = [atom_name(atom) for atom in molecule.GetAtoms()]
    if any(label not in names for label in donor_labels):
        result["status"] = "donor_atom_name_map_failed"
        return result

    try:
        canonical_smiles = Chem.MolToSmiles(
            molecule, canonical=True, isomericSmiles=True
        )
        donor_keys = donor_key_multiset(molecule, donor_labels)
    except Exception:
        result["status"] = "canonical_payload_or_rank_failed"
        return result
    if donor_keys is None:
        result["status"] = "donor_atom_name_map_failed"
        return result

    result["canonical_smiles"] = canonical_smiles
    result["donor_key_multiset"] = donor_keys
    rng = random.Random(f"{row['scope']}:{row['refcode']}:{row['component_index']}")
    for _ in range(N_PERMUTATIONS):
        order = list(range(molecule.GetNumAtoms()))
        rng.shuffle(order)
        permuted = Chem.RenumberAtoms(molecule, order)
        result["permutation_trials"] += 1
        try:
            permuted_smiles = Chem.MolToSmiles(
                permuted, canonical=True, isomericSmiles=True
            )
            permuted_keys = donor_key_multiset(permuted, donor_labels)
        except Exception:
            continue
        if permuted_smiles == canonical_smiles and permuted_keys == donor_keys:
            result["permutation_passed"] += 1

    result["status"] = (
        "full_payload_atom_map_invariant"
        if result["permutation_passed"] == result["permutation_trials"]
        else "permutation_invariance_failed"
    )
    return result


def main() -> int:
    RDLogger.DisableLog("rdApp.warning")
    RDLogger.DisableLog("rdApp.error")
    rows = [
        json.loads(line)
        for line in (OUTDIR / "phase1_donor_components.jsonl").read_text(
            encoding="utf-8"
        ).splitlines()
        if line.strip()
    ]
    audited = [audit_component(row) for row in rows]

    output_path = OUTDIR / "phase2_component_results.jsonl"
    with output_path.open("w", encoding="utf-8") as handle:
        for row in audited:
            handle.write(json.dumps(row, ensure_ascii=False) + "\n")

    status_counts = Counter(row["status"] for row in audited)
    matched = [row for row in audited if row["metal_scope_match"]]
    matched_status_counts = Counter(row["status"] for row in matched)
    matched_by_entry = defaultdict(list)
    for row in matched:
        matched_by_entry[(row["scope"], row["refcode"])].append(row)
    entry_success = {
        f"{scope}:{refcode}": bool(component_rows)
        and all(
            row["status"] == "full_payload_atom_map_invariant"
            for row in component_rows
        )
        for (scope, refcode), component_rows in matched_by_entry.items()
    }

    summary = {
        "definition": "CSD metal-stripped donor component -> Mol2 -> RDKit sanitized molecule -> canonical isomeric SMILES and donor symmetry-rank multiset",
        "permutations_per_parseable_component": N_PERMUTATIONS,
        "all_extracted_components": len(audited),
        "all_status_counts": dict(status_counts),
        "metal_scope_matched_components": len(matched),
        "metal_scope_matched_status_counts": dict(matched_status_counts),
        "metal_scope_matched_entries_with_components": len(matched_by_entry),
        "full_payload_success_entries": sum(entry_success.values()),
        "entry_success": entry_success,
        "permutation_trials_on_matched_components": sum(
            row["permutation_trials"] for row in matched
        ),
        "permutation_passed_on_matched_components": sum(
            row["permutation_passed"] for row in matched
        ),
        "interpretation": (
            "NO-GO for a full ligand-chemical-identity claim from the current 55/33 template set: "
            "most refcodes are absent or chemically mismatched, and eta5 Cp components fail neutral aromatic Mol2/SMILES sanitization."
        ),
    }
    (OUTDIR / "phase2_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
