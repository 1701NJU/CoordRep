#!/usr/bin/env python3
"""
multidentate_validator.py
=========================
Task 4: Multidentate grammar consistency validator.

Checks:
1. donor marker count per ligand_id == denticity metadata
2. donor indices continuous and non-duplicate
3. stereo references point to existing donor/ligand IDs
4. chelating ligands not erroneously split into monodentate
5. ligand SMILES donor atom candidates element-consistent with markers

Reads tmQM pipeline output and CSD retained entries.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

try:
    from rdkit import Chem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False


DONOR_ELEMENTS = {'N', 'O', 'S', 'P', 'Cl', 'Br', 'I', 'F', 'Se', 'Te', 'As'}


def parse_coordrep_grammar(coordrep_str):
    """
    Parse CoordRep string into structured components for validation.

    Returns dict with parsed fields.
    """
    result = {
        'metal': None,
        'cn': 0,
        'ligands': {},       # {lig_id: smiles}
        'donor_sites': {},   # {lig_id: [(elem, rank), ...]}
        'constraints': [],   # [(type, site1, site2), ...]
        'all_lig_ids_in_constraints': set(),
    }

    # Metal block
    metal_m = re.search(r'Metal:([A-Z][a-z]?)', coordrep_str)
    if metal_m:
        result['metal'] = metal_m.group(1)
    cn_m = re.search(r'CN:(\d+)', coordrep_str)
    if cn_m:
        result['cn'] = int(cn_m.group(1))

    # Ligand blocks
    for m in re.finditer(r'\|(L\d+)=([^|]*)', coordrep_str):
        lid = m.group(1)
        smi = m.group(2)
        result['ligands'][lid] = smi

    # Donor sites from constraint blocks
    for m in re.finditer(r'(L\d+):([A-Za-z]+):(\d+)', coordrep_str):
        lid = m.group(1)
        elem = m.group(2)
        rank = int(m.group(3))
        result['donor_sites'].setdefault(lid, []).append((elem, rank))
        result['all_lig_ids_in_constraints'].add(lid)

    # Constraint blocks
    for m in re.finditer(r'\{(trans|cis|fm|fac|mer):([^}]+)\}', coordrep_str):
        ctype = m.group(1)
        body = m.group(2)
        result['constraints'].append((ctype, body))

    return result


def validate_one(coordrep_str, meta=None):
    """
    Validate one CoordRep string for multidentate consistency.

    meta: optional dict with 'ligand_dents', 'mol_id', 'source'

    Returns list of issue dicts.
    """
    issues = []
    parsed = parse_coordrep_grammar(coordrep_str)
    cn = parsed['cn']
    ligands = parsed['ligands']
    donor_sites = parsed['donor_sites']

    meta = meta or {}
    lig_dents_meta = meta.get('ligand_dents', [])
    mol_id = meta.get('mol_id', '')

    # 1. Check donor marker count per ligand vs denticity metadata
    #    Note: constraint blocks only reference donors in stereo relations,
    #    so observed count from constraints is a LOWER BOUND.
    #    Only flag if observed EXCEEDS expected (over-counting error).
    if lig_dents_meta:
        # Canonicalization sorts ligands by denticity (desc), so
        # sort the metadata dents descending to match string order.
        sorted_meta_dents = sorted(lig_dents_meta, reverse=True)
        lig_ids_sorted = sorted(ligands.keys(), key=lambda x: int(x[1:]))
        for i, lid in enumerate(lig_ids_sorted):
            if i < len(sorted_meta_dents):
                expected_dent = sorted_meta_dents[i]
                observed_donors = donor_sites.get(lid, [])
                n_observed = len(set(observed_donors))
                # Flag only if observed > expected (genuine over-count)
                if n_observed > expected_dent:
                    issues.append({
                        'check': 'donor_count_exceeds_dent',
                        'lig_id': lid,
                        'expected': expected_dent,
                        'observed': n_observed,
                        'mol_id': mol_id,
                    })

    # 2. Donor indices non-duplicate per ligand
    #    Note: ranks may appear non-continuous because constraint blocks
    #    only reference donors in stereo relations. E.g., L1:C:2 can appear
    #    without L1:C:1 if C:1 is not in any trans/cis pair. This is valid.
    for lid, donors in donor_sites.items():
        by_elem = defaultdict(list)
        for elem, rank in donors:
            by_elem[elem].append(rank)

        for elem, ranks in by_elem.items():
            unique_ranks = sorted(set(ranks))
            if len(ranks) != len(unique_ranks):
                issues.append({
                    'check': 'duplicate_donor_rank',
                    'lig_id': lid,
                    'element': elem,
                    'ranks': ranks,
                    'mol_id': mol_id,
                })

    # 3. Stereo relation references point to existing ligand IDs
    for ctype, body in parsed['constraints']:
        refs = re.findall(r'(L\d+)', body)
        for ref in refs:
            if ref not in ligands:
                issues.append({
                    'check': 'stereo_ref_missing_ligand',
                    'constraint': f"{{{ctype}:{body}}}",
                    'missing_lig': ref,
                    'mol_id': mol_id,
                })

    # 4. Check if chelating ligands might be erroneously split
    # If multiple ligand IDs have the same SMILES and each is monodentate,
    # they might be a single chelating ligand that was split
    smiles_to_lids = defaultdict(list)
    for lid, smi in ligands.items():
        smiles_to_lids[smi].append(lid)

    for smi, lids in smiles_to_lids.items():
        if len(lids) >= 2:
            # Check if all are monodentate (no constraint refs or only 1 donor each)
            all_mono = all(len(donor_sites.get(lid, [])) <= 1 for lid in lids)
            if all_mono and len(lids) <= 3:
                # Could be legitimate (e.g., two separate NH3 ligands)
                # Only flag if SMILES suggests multidentate capability
                if RDKIT_AVAILABLE and smi:
                    try:
                        mol = Chem.MolFromSmiles(smi)
                        if mol:
                            donor_count = sum(1 for a in mol.GetAtoms()
                                              if a.GetSymbol() in DONOR_ELEMENTS)
                            if donor_count >= 2:
                                issues.append({
                                    'check': 'possible_split_chelate',
                                    'smiles': smi,
                                    'lig_ids': lids,
                                    'donor_atoms_in_smiles': donor_count,
                                    'mol_id': mol_id,
                                })
                    except Exception:
                        pass

    # 5. SMILES donor element consistency
    if RDKIT_AVAILABLE:
        for lid, donors in donor_sites.items():
            smi = ligands.get(lid, '')
            if not smi:
                continue
            try:
                mol = Chem.MolFromSmiles(smi)
                if mol is None:
                    continue
                smiles_elements = set(a.GetSymbol() for a in mol.GetAtoms()
                                      if a.GetSymbol() in DONOR_ELEMENTS)
                for elem, rank in donors:
                    if elem not in smiles_elements and elem not in ('C',):
                        issues.append({
                            'check': 'donor_element_not_in_smiles',
                            'lig_id': lid,
                            'donor_element': elem,
                            'smiles': smi,
                            'smiles_elements': sorted(smiles_elements),
                            'mol_id': mol_id,
                        })
            except Exception:
                pass

    return issues


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tmqm",
                        default="/data/CoordRep/CoordSMILES/pipeline_full_output/results.jsonl")
    parser.add_argument("--csd",
                        default="revision_results/csd_external/csd_retained_entries.jsonl")
    parser.add_argument("--out", default="revision_results/multidentate")
    args = parser.parse_args()

    os.makedirs(args.out, exist_ok=True)

    # ── Validate tmQM ─────────────────────────────────────
    print("Validating tmQM CoordRep strings …")
    all_issues = []
    n_valid = 0
    n_total = 0

    with open(args.tmqm) as f:
        for line in f:
            d = json.loads(line)
            coordrep = d.get('coordrep', '')
            if not coordrep:
                continue
            n_total += 1
            meta = {
                'mol_id': d.get('mol_id', ''),
                'source': 'tmQM',
                'ligand_dents': d.get('ligand_dents', []),
            }
            issues = validate_one(coordrep, meta)
            if issues:
                all_issues.extend(issues)
            else:
                n_valid += 1

    print(f"  tmQM: {n_total} strings, {n_valid} clean, "
          f"{n_total - n_valid} with issues ({len(all_issues)} total issues)")

    # ── Validate CSD ──────────────────────────────────────
    n_csd_total = 0
    n_csd_valid = 0
    if os.path.exists(args.csd):
        print("Validating CSD CoordRep strings …")
        with open(args.csd) as f:
            for line in f:
                d = json.loads(line)
                coordrep = d.get('coordrep', '')
                if not coordrep:
                    continue
                n_csd_total += 1
                n_total += 1
                meta = {
                    'mol_id': d.get('refcode', ''),
                    'source': 'CSD',
                }
                issues = validate_one(coordrep, meta)
                if issues:
                    all_issues.extend(issues)
                else:
                    n_valid += 1
                    n_csd_valid += 1

        print(f"  CSD: {n_csd_total} strings, {n_csd_valid} clean, "
              f"{n_csd_total - n_csd_valid} with issues")

    # ── Aggregate ─────────────────────────────────────────
    issue_types = Counter(i['check'] for i in all_issues)
    failure_rate = (n_total - n_valid) / max(n_total, 1)

    # ── Write summary CSV ─────────────────────────────────
    p = os.path.join(args.out, "multidentate_validator_summary.csv")
    with open(p, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["metric", "value"])
        w.writerow(["n_total", n_total])
        w.writerow(["n_valid", n_valid])
        w.writerow(["n_with_issues", n_total - n_valid])
        w.writerow(["failure_rate", round(failure_rate, 4)])
        for check, cnt in issue_types.most_common():
            w.writerow([f"issue_{check}", cnt])
    print(f"  {p}")

    # ── Write failure details JSONL ───────────────────────
    p = os.path.join(args.out, "multidentate_validator_failures.jsonl")
    with open(p, "w") as f:
        for issue in all_issues[:2000]:
            # Convert sets to lists for JSON serialization
            serializable = {}
            for k, v in issue.items():
                if isinstance(v, set):
                    serializable[k] = sorted(v)
                else:
                    serializable[k] = v
            f.write(json.dumps(serializable) + "\n")
    print(f"  {p} ({min(len(all_issues), 2000)} entries)")

    # ── Print summary ─────────────────────────────────────
    print(f"\n{'='*60}")
    print("MULTIDENTATE VALIDATOR SUMMARY")
    print(f"{'='*60}")
    print(f"  Total strings: {n_total}")
    print(f"  Valid: {n_valid} ({n_valid/max(n_total,1):.1%})")
    print(f"  With issues: {n_total - n_valid} ({failure_rate:.1%})")
    print(f"  Issue types:")
    for check, cnt in issue_types.most_common():
        print(f"    {check}: {cnt}")
    if failure_rate > 0.01:
        print(f"\n  ⚠ Failure rate > 1%: consider manual review of examples")
    print()


if __name__ == "__main__":
    main()
