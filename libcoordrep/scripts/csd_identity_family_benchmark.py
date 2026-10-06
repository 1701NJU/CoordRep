#!/usr/bin/env python3
"""
csd_identity_family_benchmark.py
================================
Task 2: CSD refcode-family identity benchmark.

Reads retained entries from Task 1, groups by refcode family,
and tests within-family agreement at each identity layer (L0-L3).

Usage:
    python scripts/csd_identity_family_benchmark.py [--retained PATH]
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from collections import Counter, defaultdict
from typing import Dict, List

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))


def _pairwise_match_rate(values: List[str]) -> float:
    """Fraction of all ordered pairs that are identical."""
    n = len(values)
    if n < 2:
        return 1.0
    matches = sum(1 for i in range(n) for j in range(i + 1, n)
                  if values[i] == values[j])
    total = n * (n - 1) // 2
    return matches / total if total > 0 else 1.0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--retained",
                        default="revision_results/csd_external/csd_retained_entries.jsonl")
    parser.add_argument("--out", default="revision_results/csd_external")
    parser.add_argument("--min_family_size", type=int, default=2)
    parser.add_argument("--max_examples", type=int, default=50,
                        help="Max example families for figure JSONL")
    args = parser.parse_args()

    os.makedirs(args.out, exist_ok=True)

    # ── Load retained entries ─────────────────────────────────
    print("Loading retained entries …")
    entries = []
    with open(args.retained) as f:
        for line in f:
            entries.append(json.loads(line))
    print(f"  {len(entries)} entries")

    # ── Group by family ───────────────────────────────────────
    families: Dict[str, List[dict]] = defaultdict(list)
    for e in entries:
        families[e["family"]].append(e)

    # Filter to families with >= min_family_size
    fam_list = [(fam, members) for fam, members in families.items()
                if len(members) >= args.min_family_size]
    fam_list.sort(key=lambda x: -len(x[1]))
    print(f"  {len(fam_list)} families with size >= {args.min_family_size}")

    # ── Compute per-family identity metrics ───────────────────
    detail_rows = []
    layer_match_sums = {"L0": 0.0, "L1": 0.0, "L2": 0.0, "L3": 0.0}
    boundary_in_families = 0
    total_members_in_families = 0
    cn_mismatch_families = 0
    ligand_mismatch_families = 0
    stereo_mismatch_families = 0
    false_nondup_L0 = 0  # families where L0 disagrees but L3 agrees

    for fam, members in fam_list:
        n = len(members)
        total_members_in_families += n

        l0s = [m["L0"] for m in members]
        l1s = [m["L1"] for m in members]
        l2s = [m["L2"] for m in members]
        l3s = [m["L3"] for m in members]

        l0_match = _pairwise_match_rate(l0s)
        l1_match = _pairwise_match_rate(l1s)
        l2_match = _pairwise_match_rate(l2s)
        l3_match = _pairwise_match_rate(l3s)

        layer_match_sums["L0"] += l0_match
        layer_match_sums["L1"] += l1_match
        layer_match_sums["L2"] += l2_match
        layer_match_sums["L3"] += l3_match

        if l0_match < 1.0 and l3_match == 1.0:
            false_nondup_L0 += 1

        # boundary count
        n_boundary = sum(1 for m in members if m.get("is_boundary", False))
        boundary_in_families += n_boundary

        # CN variation
        cns = set(m["cn"] for m in members)
        if len(cns) > 1:
            cn_mismatch_families += 1

        # Metal variation (should all be same in a family)
        metals = set(m["metal"] for m in members)

        # Shape variation
        shapes = set(m.get("best_shape", "") for m in members)
        if len(shapes) > 1:
            stereo_mismatch_families += 1

        detail_rows.append({
            "family": fam,
            "size": n,
            "metals": ";".join(sorted(metals)),
            "cns": ";".join(str(c) for c in sorted(cns)),
            "L0_match": round(l0_match, 4),
            "L1_match": round(l1_match, 4),
            "L2_match": round(l2_match, 4),
            "L3_match": round(l3_match, 4),
            "n_boundary": n_boundary,
            "n_unique_L0": len(set(l0s)),
            "n_unique_L1": len(set(l1s)),
            "n_unique_L2": len(set(l2s)),
            "n_unique_L3": len(set(l3s)),
        })

    n_fam = len(fam_list)

    # ── Summary metrics ───────────────────────────────────────
    summary = {
        "n_families": n_fam,
        "mean_family_size": round(total_members_in_families / max(n_fam, 1), 2),
        "L0_within_family_match": round(layer_match_sums["L0"] / max(n_fam, 1), 4),
        "L1_within_family_match": round(layer_match_sums["L1"] / max(n_fam, 1), 4),
        "L2_within_family_match": round(layer_match_sums["L2"] / max(n_fam, 1), 4),
        "L3_within_family_match": round(layer_match_sums["L3"] / max(n_fam, 1), 4),
        "false_nonduplicate_L0": false_nondup_L0,
        "boundary_fraction_in_families": round(
            boundary_in_families / max(total_members_in_families, 1), 4),
        "cn_mismatch_families": cn_mismatch_families,
        "stereo_mismatch_families": stereo_mismatch_families,
    }

    # ── Layer retention ───────────────────────────────────────
    # How many families have perfect agreement at each layer
    layer_retention = []
    for layer in ["L0", "L1", "L2", "L3"]:
        n_perfect = sum(1 for r in detail_rows if r[f"{layer}_match"] == 1.0)
        layer_retention.append({
            "layer": layer,
            "n_perfect_families": n_perfect,
            "pct_perfect": round(100 * n_perfect / max(n_fam, 1), 2),
            "mean_match": round(layer_match_sums[layer] / max(n_fam, 1), 4),
        })

    # ── Write outputs ─────────────────────────────────────────
    print("\nWriting outputs …")

    # summary
    p = os.path.join(args.out, "csd_family_identity_summary.csv")
    with open(p, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["metric", "value"])
        for k, v in summary.items():
            w.writerow([k, v])
    print(f"  {p}")

    # detail
    p = os.path.join(args.out, "csd_family_identity_detail.csv")
    if detail_rows:
        cols = list(detail_rows[0].keys())
        with open(p, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=cols)
            w.writeheader()
            for r in detail_rows:
                w.writerow(r)
    print(f"  {p}")

    # layer retention
    p = os.path.join(args.out, "csd_identity_layer_retention.csv")
    with open(p, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["layer", "n_perfect_families",
                                          "pct_perfect", "mean_match"])
        w.writeheader()
        for r in layer_retention:
            w.writerow(r)
    print(f"  {p}")

    # examples for figure (top families)
    p = os.path.join(args.out, "csd_family_examples_for_figure.jsonl")
    examples = []
    # Pick families with interesting L0 variation
    varied = [r for r in detail_rows if r["L0_match"] < 1.0 and r["L3_match"] == 1.0]
    varied.sort(key=lambda r: r["size"], reverse=True)
    for r in varied[:args.max_examples]:
        fam = r["family"]
        members = families[fam]
        examples.append({
            "family": fam,
            "size": r["size"],
            "L0_match": r["L0_match"],
            "L3_match": r["L3_match"],
            "refcodes": [m["refcode"] for m in members],
            "metals": list(set(m["metal"] for m in members)),
            "cns": list(set(m["cn"] for m in members)),
            "shapes": list(set(m.get("best_shape", "") for m in members)),
        })
    with open(p, "w") as f:
        for ex in examples:
            f.write(json.dumps(ex) + "\n")
    print(f"  {p} ({len(examples)} examples)")

    # update summary.json
    sum_path = os.path.join(args.out, "summary.json")
    if os.path.exists(sum_path):
        with open(sum_path) as f:
            existing = json.load(f)
        existing.update(summary)
        summary = existing
    with open(sum_path, "w") as f:
        json.dump(summary, f, indent=2)

    # ── Print key numbers ─────────────────────────────────────
    print(f"\n{'='*60}")
    print("CSD IDENTITY FAMILY BENCHMARK")
    print(f"{'='*60}")
    for k, v in summary.items():
        print(f"  {k}: {v}")
    print("\n  Layer retention:")
    for r in layer_retention:
        print(f"    {r['layer']}: {r['pct_perfect']:.1f}% perfect, "
              f"mean_match={r['mean_match']:.3f}")
    print()


if __name__ == "__main__":
    main()
