#!/usr/bin/env python3
"""
run_csd_identity_baselines.py
=============================
Part 3: CSD family identity / duplicate linking baselines.

Compares:
- Canonical SMILES multiset key
- WL graph hash key
- ECFP Tanimoto ligand-set similarity
- CoordRep L0/L1/L2/L3 identity keys (from existing results)

Outputs: revision_results/gnn_baselines/csd_identity_*.csv
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Set, Tuple

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from coordrep_tools.gnn_baselines.non_neural import (
    canonical_smiles_multiset_key,
    ligand_set_wl_key,
    ligand_set_tanimoto,
    wl_hash,
    ecfp_fingerprint,
    tanimoto_similarity,
)


# ── Load CSD data ────────────────────────────────────────

CSD_PATH = "inputs/csd/csd_retained_entries.jsonl"
CSD_FAMILY_PATH = "inputs/csd/csd_family_identity_summary.csv"


def _extract_ligand_smiles_from_coordrep(coordrep: str) -> List[str]:
    """Extract ligand SMILES from |L1=...|L2=...| blocks in CoordRep string."""
    import re
    return re.findall(r'\|L\d+=([^|]+)\|', coordrep)


def load_csd_entries(path: str) -> List[dict]:
    entries = []
    with open(path) as f:
        for line in f:
            d = json.loads(line)
            # Fill in ligand_smiles from coordrep if missing
            if 'ligand_smiles' not in d and 'coordrep' in d:
                d['ligand_smiles'] = _extract_ligand_smiles_from_coordrep(d['coordrep'])
            # Fill in donor_elements from coordrep if missing
            if 'donor_elements' not in d and 'coordrep' in d:
                import re
                donors = re.findall(r':([A-Z][a-z]?):\d+', d['coordrep'])
                d['donor_elements'] = donors
            entries.append(d)
    return entries


def extract_families(entries: List[dict]) -> Dict[str, List[dict]]:
    """Group entries by refcode family (first 6 chars of refcode)."""
    families = defaultdict(list)
    for e in entries:
        refcode = e.get('refcode', e.get('mol_id', ''))
        if len(refcode) >= 6:
            family = refcode[:6]
        else:
            family = refcode
        families[family].append(e)
    return families


# ── Identity key generators ───────────────────────────────

def make_smiles_key(entry: dict) -> str:
    """Canonical SMILES multiset key."""
    smiles = entry.get('ligand_smiles', [])
    if not smiles:
        return ''
    return canonical_smiles_multiset_key(smiles)


def make_metal_smiles_key(entry: dict) -> str:
    """Metal + SMILES multiset key."""
    metal = entry.get('metal', '?')
    cn = entry.get('cn', 0)
    smiles = entry.get('ligand_smiles', [])
    smi_key = canonical_smiles_multiset_key(smiles) if smiles else ''
    return f"{metal}|{cn}|{smi_key}"


def make_wl_key(entry: dict) -> str:
    """WL graph hash key over ligand set."""
    smiles = entry.get('ligand_smiles', [])
    metal = entry.get('metal', '?')
    cn = entry.get('cn', 0)
    if not smiles:
        return ''
    return ligand_set_wl_key(smiles, metal, cn)


def make_coordrep_l0(entry: dict) -> str:
    """CoordRep L0: use pre-computed if available, else metal+CN."""
    if 'L0' in entry and entry['L0']:
        return entry['L0']
    return f"{entry.get('metal','?')}|{entry.get('cn',0)}"


def make_coordrep_l1(entry: dict) -> str:
    """CoordRep L1: use pre-computed if available."""
    if 'L1' in entry and entry['L1']:
        return entry['L1']
    donors = entry.get('donor_elements', [])
    donor_str = ','.join(sorted(donors))
    return f"{make_coordrep_l0(entry)}|{donor_str}"


def make_coordrep_l2(entry: dict) -> str:
    """CoordRep L2: use pre-computed if available."""
    if 'L2' in entry and entry['L2']:
        return entry['L2']
    # Fallback: L1 + sorted ligand SMILES
    smiles = entry.get('ligand_smiles', [])
    smi_key = canonical_smiles_multiset_key(smiles) if smiles else ''
    return f"{make_coordrep_l1(entry)}|{smi_key}"


def make_coordrep_l3(entry: dict) -> str:
    """CoordRep L3: use pre-computed if available."""
    if 'L3' in entry and entry['L3']:
        return entry['L3']
    return make_coordrep_l2(entry)


# ── Evaluation metrics ────────────────────────────────────

def evaluate_identity_method(families: Dict[str, List[dict]],
                              key_fn, method_name: str) -> dict:
    """
    Evaluate an identity key method:
    - within_family_recall: fraction of pairs in same family sharing same key
    - cross_family_precision: fraction of same-key pairs that are in same family
    - false_non_duplicate: same family, different key
    - false_merge: different family, same key
    """
    # Build key → entries mapping
    key_to_entries = defaultdict(list)
    entry_to_family = {}

    for family_id, entries in families.items():
        if len(entries) < 2:
            continue
        for e in entries:
            refcode = e.get('refcode', e.get('mol_id', ''))
            key = key_fn(e)
            key_to_entries[key].append(refcode)
            entry_to_family[refcode] = family_id

    # Within-family recall: for each family with >=2 members,
    # what fraction of pairs share the same key?
    within_pairs = 0
    within_matches = 0
    for family_id, entries in families.items():
        if len(entries) < 2:
            continue
        keys = [key_fn(e) for e in entries]
        for i in range(len(keys)):
            for j in range(i + 1, len(keys)):
                within_pairs += 1
                if keys[i] == keys[j]:
                    within_matches += 1

    within_recall = within_matches / max(within_pairs, 1)

    # Cross-family precision: among all entries sharing the same key,
    # what fraction are in the same family?
    cross_pairs = 0
    cross_correct = 0
    for key, refcodes in key_to_entries.items():
        if len(refcodes) < 2:
            continue
        for i in range(len(refcodes)):
            for j in range(i + 1, len(refcodes)):
                cross_pairs += 1
                if entry_to_family.get(refcodes[i]) == entry_to_family.get(refcodes[j]):
                    cross_correct += 1

    cross_precision = cross_correct / max(cross_pairs, 1)

    # Count unique keys
    all_keys = set()
    n_entries = 0
    for entries in families.values():
        for e in entries:
            all_keys.add(key_fn(e))
            n_entries += 1

    return {
        'method': method_name,
        'n_families': len([f for f in families.values() if len(f) >= 2]),
        'n_entries': n_entries,
        'n_unique_keys': len(all_keys),
        'within_family_recall': round(within_recall, 4),
        'cross_family_precision': round(cross_precision, 4),
        'within_pairs': within_pairs,
        'within_matches': within_matches,
    }


def evaluate_similarity_retrieval(families: Dict[str, List[dict]],
                                   max_families: int = 200) -> dict:
    """
    Evaluate ECFP Tanimoto similarity for top-k retrieval within families.
    """
    # Build a fingerprint database
    all_entries = []
    family_labels = []
    for fam_id, entries in families.items():
        if len(entries) < 2:
            continue
        for e in entries:
            all_entries.append(e)
            family_labels.append(fam_id)
        if len(set(family_labels)) >= max_families:
            break

    if len(all_entries) < 10:
        return {'method': 'ECFP_Tanimoto', 'recall_at_5': 0, 'recall_at_10': 0}

    # Compute pairwise Tanimoto on ligand-set level
    fps = []
    for e in all_entries:
        smiles = e.get('ligand_smiles', [])
        if smiles:
            combined = '.'.join(smiles)  # crude but fast
            fps.append(ecfp_fingerprint(combined))
        else:
            fps.append(None)

    # For each entry, retrieve top-k most similar and check family match
    n = len(all_entries)
    recall_at_5 = []
    recall_at_10 = []

    for i in range(n):
        if fps[i] is None:
            continue
        target_fam = family_labels[i]
        sims = []
        for j in range(n):
            if i == j:
                continue
            sim = tanimoto_similarity(fps[i], fps[j])
            sims.append((sim, j))
        sims.sort(key=lambda x: -x[0])

        top5_fams = [family_labels[j] for _, j in sims[:5]]
        top10_fams = [family_labels[j] for _, j in sims[:10]]
        recall_at_5.append(int(target_fam in top5_fams))
        recall_at_10.append(int(target_fam in top10_fams))

    return {
        'method': 'ECFP_Tanimoto_retrieval',
        'n_entries': n,
        'recall_at_5': round(float(np.mean(recall_at_5)), 4) if recall_at_5 else 0,
        'recall_at_10': round(float(np.mean(recall_at_10)), 4) if recall_at_10 else 0,
    }


# ── Main ──────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--csd", default=CSD_PATH)
    parser.add_argument("--out", default="revision_results/gnn_baselines")
    parser.add_argument("--max_families", type=int, default=300)
    args = parser.parse_args()

    os.makedirs(args.out, exist_ok=True)

    # Load CSD data
    print("Loading CSD entries …")
    if not os.path.exists(args.csd):
        print(f"  ERROR: {args.csd} not found. Skipping CSD identity.")
        return

    entries = load_csd_entries(args.csd)
    families = extract_families(entries)
    multi_fam = {k: v for k, v in families.items() if len(v) >= 2}
    print(f"  {len(entries)} entries, {len(families)} families, "
          f"{len(multi_fam)} with >=2 members")

    # Evaluate identity methods
    print(f"\n{'='*60}")
    print("Evaluating identity methods")
    print(f"{'='*60}")

    methods = [
        ('SMILES_multiset', make_smiles_key),
        ('Metal+SMILES', make_metal_smiles_key),
        ('WL_hash', make_wl_key),
        ('CoordRep_L0', make_coordrep_l0),
        ('CoordRep_L1', make_coordrep_l1),
        ('CoordRep_L2', make_coordrep_l2),
        ('CoordRep_L3', make_coordrep_l3),
    ]

    results = []
    for name, fn in methods:
        r = evaluate_identity_method(multi_fam, fn, name)
        results.append(r)
        print(f"  {name}: recall={r['within_family_recall']:.3f} "
              f"precision={r['cross_family_precision']:.3f} "
              f"unique_keys={r['n_unique_keys']}")

    # Similarity retrieval
    print("\nECFP Tanimoto retrieval …")
    sim_result = evaluate_similarity_retrieval(multi_fam, args.max_families)
    print(f"  recall@5={sim_result.get('recall_at_5', 0):.3f} "
          f"recall@10={sim_result.get('recall_at_10', 0):.3f}")

    # Write outputs
    print(f"\n{'='*60}")
    print("Writing outputs")
    print(f"{'='*60}")

    # Summary CSV
    p = os.path.join(args.out, "csd_identity_baselines.csv")
    if results:
        keys = list(results[0].keys())
        with open(p, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            for r in results:
                w.writerow(r)
    print(f"  {p}")

    # Precision-recall CSV
    p = os.path.join(args.out, "csd_identity_precision_recall.csv")
    with open(p, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["method", "within_family_recall", "cross_family_precision",
                     "n_unique_keys"])
        for r in results:
            w.writerow([r['method'], r['within_family_recall'],
                        r['cross_family_precision'], r['n_unique_keys']])
    print(f"  {p}")

    # Summary JSON
    summary = {
        'identity_methods': {r['method']: r for r in results},
        'similarity_retrieval': sim_result,
        'n_families': len(multi_fam),
        'n_entries': len(entries),
    }
    p = os.path.join(args.out, "csd_identity_summary.json")
    with open(p, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"  {p}")

    # Final table
    print(f"\n{'='*60}")
    print("CSD IDENTITY BASELINES")
    print(f"{'='*60}")
    print(f"  {'Method':<20s} {'Recall':>8s} {'Precision':>10s} {'Keys':>8s}")
    print(f"  {'-'*48}")
    for r in results:
        print(f"  {r['method']:<20s} {r['within_family_recall']:>8.3f} "
              f"{r['cross_family_precision']:>10.3f} {r['n_unique_keys']:>8d}")


if __name__ == "__main__":
    main()
