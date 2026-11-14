#!/usr/bin/env python3
"""
multidentate_coverage.py
========================
Task 1: Multidentate coverage statistics across tmQM, COD, and CSD.

Computes complex-level and ligand-level denticity distributions,
stratified by source, metal row, and CN.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

METAL_ROW = {}
for m in ['Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn']:
    METAL_ROW[m] = '3d'
for m in ['Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd']:
    METAL_ROW[m] = '4d'
for m in ['La','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg']:
    METAL_ROW[m] = '5d'


def parse_dent_from_coordrep(coordrep: str):
    """
    Parse denticity info from a CoordRep string.

    Returns dict: {lig_id: dent, ...}, cn, metal
    """
    # Parse CN
    cn_m = re.search(r'CN:(\d+)', coordrep)
    cn = int(cn_m.group(1)) if cn_m else 0

    # Parse metal
    metal_m = re.search(r'Metal:([A-Z][a-z]?)', coordrep)
    metal = metal_m.group(1) if metal_m else "?"

    # Parse ligand IDs from ligand blocks
    lig_ids = re.findall(r'\|(L\d+)=', coordrep)

    # Parse donor sites from constraint blocks
    donor_refs = re.findall(r'(L\d+):([A-Za-z]+):(\d+)', coordrep)
    lig_donor_sites = defaultdict(set)
    for lid, elem, rank in donor_refs:
        lig_donor_sites[lid].add((elem, rank))

    # Assign denticity
    dents = {}
    accounted = 0
    for lid in lig_ids:
        if lid in lig_donor_sites:
            d = len(lig_donor_sites[lid])
            dents[lid] = d
            accounted += d
        else:
            dents[lid] = 1  # at least monodentate
            accounted += 1

    # Distribute remaining CN to unconstrained ligands
    remaining = cn - accounted
    if remaining > 0:
        unconstrained = [lid for lid in lig_ids if lid not in lig_donor_sites]
        if unconstrained:
            extra_each = remaining // len(unconstrained)
            leftover = remaining % len(unconstrained)
            for i, lid in enumerate(unconstrained):
                dents[lid] += extra_each + (1 if i < leftover else 0)

    return dents, cn, metal


def load_tmqm_data(path: str):
    """Load tmQM pipeline output, return list of complex records."""
    records = []
    with open(path) as f:
        for line in f:
            d = json.loads(line)
            if not d.get('coordrep'):
                continue
            lig_dents = d.get('ligand_dents', [])
            if not lig_dents:
                continue
            records.append({
                'source': 'tmQM',
                'id': d['mol_id'],
                'metal': d.get('metal', '?'),
                'cn': d.get('cn', 0),
                'ligand_dents': lig_dents,
                'coordrep': d['coordrep'],
            })
    return records


def load_csd_data(path: str):
    """Load CSD retained entries, parse denticity from CoordRep strings."""
    records = []
    with open(path) as f:
        for line in f:
            d = json.loads(line)
            coordrep = d.get('coordrep', '')
            if not coordrep:
                continue
            dents, cn, metal = parse_dent_from_coordrep(coordrep)
            lig_dents = list(dents.values())
            records.append({
                'source': 'CSD',
                'id': d['refcode'],
                'metal': d.get('metal', metal),
                'cn': d.get('cn', cn),
                'ligand_dents': lig_dents,
                'coordrep': coordrep,
            })
    return records


def compute_stats(records):
    """Compute all multidentate coverage statistics."""
    stats = {
        'n_complexes': 0,
        'mono_only': 0,
        'has_bi': 0,
        'has_ge3': 0,
        'max_dent_dist': Counter(),
        'lig_dent_dist': Counter(),
        'n_ligands': 0,
        'by_source': defaultdict(lambda: {
            'n_complexes': 0, 'mono_only': 0, 'has_bi': 0, 'has_ge3': 0,
            'n_ligands': 0, 'lig_dent_dist': Counter(),
        }),
        'by_cn': defaultdict(lambda: {
            'n_complexes': 0, 'has_multidentate': 0,
        }),
        'by_metal_row': defaultdict(lambda: {
            'n_complexes': 0, 'has_multidentate': 0,
        }),
        'failures': [],
    }

    for rec in records:
        dents = rec['ligand_dents']
        source = rec['source']
        metal = rec['metal']
        cn = rec['cn']

        if not dents:
            continue

        stats['n_complexes'] += 1
        max_d = max(dents)
        is_mono_only = all(d == 1 for d in dents)
        has_bi = any(d == 2 for d in dents)
        has_ge3 = any(d >= 3 for d in dents)

        if is_mono_only:
            stats['mono_only'] += 1
        if has_bi:
            stats['has_bi'] += 1
        if has_ge3:
            stats['has_ge3'] += 1

        stats['max_dent_dist'][max_d] += 1

        for d in dents:
            stats['n_ligands'] += 1
            if d >= 4:
                stats['lig_dent_dist']['dent>=4'] += 1
            else:
                stats['lig_dent_dist'][f'dent={d}'] += 1

        # By source
        src = stats['by_source'][source]
        src['n_complexes'] += 1
        if is_mono_only:
            src['mono_only'] += 1
        if has_bi:
            src['has_bi'] += 1
        if has_ge3:
            src['has_ge3'] += 1
        src['n_ligands'] += len(dents)
        for d in dents:
            k = f'dent={d}' if d < 4 else 'dent>=4'
            src['lig_dent_dist'][k] += 1

        # By CN
        cn_key = cn if cn in (4, 5, 6) else 'other'
        stats['by_cn'][cn_key]['n_complexes'] += 1
        if not is_mono_only:
            stats['by_cn'][cn_key]['has_multidentate'] += 1

        # By metal row
        row = METAL_ROW.get(metal, 'other')
        stats['by_metal_row'][row]['n_complexes'] += 1
        if not is_mono_only:
            stats['by_metal_row'][row]['has_multidentate'] += 1

        # Check for segmentation issues
        sum_d = sum(dents)
        if cn > 0 and sum_d != cn:
            stats['failures'].append({
                'id': rec['id'],
                'source': source,
                'cn': cn,
                'ligand_dents': dents,
                'sum_dents': sum_d,
                'reason': f'dent_sum={sum_d} != cn={cn}',
            })

    return stats


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tmqm", default="inputs/tmqm/results.jsonl")
    parser.add_argument("--csd", default="revision_results/csd_external/csd_retained_entries.jsonl")
    parser.add_argument("--out", default="revision_results/multidentate")
    args = parser.parse_args()

    os.makedirs(args.out, exist_ok=True)

    # Load data
    print("Loading tmQM data …")
    tmqm = load_tmqm_data(args.tmqm)
    print(f"  {len(tmqm)} tmQM complexes")

    print("Loading CSD data …")
    csd = []
    if os.path.exists(args.csd):
        csd = load_csd_data(args.csd)
        print(f"  {len(csd)} CSD complexes")
    else:
        print("  CSD retained entries not found, skipping")

    all_records = tmqm + csd

    # Compute stats
    print("Computing statistics …")
    stats = compute_stats(all_records)

    n = stats['n_complexes']
    nl = stats['n_ligands']

    # ── summary.json ──────────────────────────────────────
    csd_src = stats['by_source'].get('CSD', {})
    csd_n = csd_src.get('n_complexes', 0)
    csd_multi = csd_n - csd_src.get('mono_only', 0)

    summary = {
        "n_complexes_total": n,
        "fraction_monodentate_only": round(stats['mono_only'] / max(n, 1), 4),
        "fraction_with_bidentate": round(stats['has_bi'] / max(n, 1), 4),
        "fraction_with_dent_ge3": round(stats['has_ge3'] / max(n, 1), 4),
        "n_ligands_total": nl,
        "ligand_fraction_dent1": round(stats['lig_dent_dist'].get('dent=1', 0) / max(nl, 1), 4),
        "ligand_fraction_dent2": round(stats['lig_dent_dist'].get('dent=2', 0) / max(nl, 1), 4),
        "ligand_fraction_dent_ge3": round(
            (stats['lig_dent_dist'].get('dent=3', 0) +
             stats['lig_dent_dist'].get('dent>=4', 0)) / max(nl, 1), 4),
        "csd_fraction_with_multidentate": round(csd_multi / max(csd_n, 1), 4),
        "n_segmentation_failures": len(stats['failures']),
    }

    with open(os.path.join(args.out, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2)

    # ── multidentate_coverage_summary.csv ─────────────────
    rows = [
        ["metric", "count", "fraction"],
        ["total_complexes", n, 1.0],
        ["monodentate_only", stats['mono_only'], round(stats['mono_only'] / max(n, 1), 4)],
        ["contains_bidentate", stats['has_bi'], round(stats['has_bi'] / max(n, 1), 4)],
        ["contains_dent_ge3", stats['has_ge3'], round(stats['has_ge3'] / max(n, 1), 4)],
        ["total_ligands", nl, ""],
    ]
    for k in sorted(stats['lig_dent_dist'].keys()):
        c = stats['lig_dent_dist'][k]
        rows.append([k, c, round(c / max(nl, 1), 4)])
    p = os.path.join(args.out, "multidentate_coverage_summary.csv")
    with open(p, "w", newline="") as f:
        csv.writer(f).writerows(rows)
    print(f"  {p}")

    # ── multidentate_by_source.csv ────────────────────────
    p = os.path.join(args.out, "multidentate_by_source.csv")
    with open(p, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["source", "n_complexes", "mono_only", "has_bidentate",
                     "has_dent_ge3", "frac_multidentate"])
        for src in ['tmQM', 'CSD']:
            d = stats['by_source'].get(src, {})
            nc = d.get('n_complexes', 0)
            mo = d.get('mono_only', 0)
            w.writerow([src, nc, mo, d.get('has_bi', 0), d.get('has_ge3', 0),
                         round((nc - mo) / max(nc, 1), 4)])
    print(f"  {p}")

    # ── multidentate_by_cn.csv ────────────────────────────
    p = os.path.join(args.out, "multidentate_by_cn.csv")
    with open(p, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["cn", "n_complexes", "has_multidentate", "frac_multidentate"])
        for cn in [4, 5, 6, 'other']:
            d = stats['by_cn'].get(cn, {})
            nc = d.get('n_complexes', 0)
            hm = d.get('has_multidentate', 0)
            w.writerow([cn, nc, hm, round(hm / max(nc, 1), 4)])
    print(f"  {p}")

    # ── multidentate_by_metal_row.csv ─────────────────────
    p = os.path.join(args.out, "multidentate_by_metal_row.csv")
    with open(p, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["metal_row", "n_complexes", "has_multidentate", "frac_multidentate"])
        for row in ['3d', '4d', '5d', 'other']:
            d = stats['by_metal_row'].get(row, {})
            nc = d.get('n_complexes', 0)
            hm = d.get('has_multidentate', 0)
            w.writerow([row, nc, hm, round(hm / max(nc, 1), 4)])
    print(f"  {p}")

    # ── Print summary ─────────────────────────────────────
    print(f"\n{'='*60}")
    print("MULTIDENTATE COVERAGE SUMMARY")
    print(f"{'='*60}")
    for k, v in summary.items():
        print(f"  {k}: {v}")
    print(f"  segmentation failures: {len(stats['failures'])}")
    if stats['failures']:
        print(f"  (first 5 failure examples:)")
        for f in stats['failures'][:5]:
            print(f"    {f['id']}: {f['reason']}")
    print()


if __name__ == "__main__":
    main()
