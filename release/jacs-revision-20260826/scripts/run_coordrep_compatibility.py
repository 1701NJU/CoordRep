#!/usr/bin/env python3
"""
run_coordrep_compatibility.py
=============================
CoordRep-Score / Counterfactual Coordination Compatibility Benchmark.

Given real CoordRep complexes and a trained MLM checkpoint, generate
counterfactual decoys, score everything with pseudo-log-likelihood,
and output stratified results proving the model learns whole-complex
coordination compatibility — not just ligand→donor lookup.

Usage (1701 env):
  python scripts/run_coordrep_compatibility.py \\
      --checkpoint .../best_model.pt \\
      --test .../fig5_tasks.jsonl \\
      --decoys-per-type 10 \\
      --out revision_results/coordrep_score
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import sys
import time
from collections import Counter, defaultdict
from typing import Dict, List

import numpy as np

# ── project imports ──────────────────────────────────────────────
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from coordrep_tools.infer import load_model_and_tokenizer
from coordrep_tools.ablation_masking import _detect_blocks
from coordrep_tools.counterfactual_decoys import (
    DecoyPool, generate_all_decoys, DECOY_TYPES, _parse_metal_token,
)
from coordrep_tools.coordrep_score import (
    compute_pseudo_log_likelihood, compute_fieldwise_scores,
)
from coordrep_tools.compatibility_benchmark import (
    run_benchmark, RandomBaseline, FrequencyBaseline,
    LigandFreqBaseline, RuleBaseline, _metal_row,
)


def _extract_unique_complexes(tasks_path: str) -> List[dict]:
    """Load one token-list per unique complex from the tasks file."""
    seen = {}
    with open(tasks_path) as f:
        for line in f:
            rec = json.loads(line)
            m = re.match(r'^(tmqm_\d+)', rec["id"])
            if not m:
                continue
            cid = m.group(1)
            if cid not in seen:
                seen[cid] = {
                    "id": cid,
                    "tokens": rec["tokens"],
                }
    return list(seen.values())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint",
                        default="inputs/checkpoints/pretrain_v3/best_model.pt")
    parser.add_argument("--test",
                        default="inputs/property_benchmarks/fig5_tasks.jsonl")
    parser.add_argument("--out",
                        default="outputs/coordrep_score")
    parser.add_argument("--decoys-per-type", type=int, default=10)
    parser.add_argument("--pll-mode", default="mask_chunked",
                        choices=["mask_one_token", "mask_field", "mask_all", "mask_chunked"])
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--max-complexes", type=int, default=0,
                        help="0 = use all")
    args = parser.parse_args()

    os.makedirs(args.out, exist_ok=True)

    # ── 1. Load data ─────────────────────────────────────────────
    print("Loading unique complexes …")
    complexes = _extract_unique_complexes(args.test)
    print(f"  {len(complexes)} unique complexes")

    all_token_lists = [c["tokens"] for c in complexes]
    all_ids = [c["id"] for c in complexes]

    # ── 2. Load model ────────────────────────────────────────────
    print("Loading model …")
    model, tokenizer = load_model_and_tokenizer(args.checkpoint, device=args.device)
    print(f"  vocab={len(tokenizer.token2id)}, device={args.device}")

    # ── 3. Build decoy pool ──────────────────────────────────────
    print("Building decoy pool …")
    t0 = time.time()
    pool = DecoyPool(all_token_lists, seed=42)
    print(f"  {len(pool.metals)} metals, {len(pool.ligand_blocks)} ligand blocks, "
          f"{len(pool.stereo_blocks)} stereo blocks ({time.time()-t0:.1f}s)")

    # ── 4. Build baselines ───────────────────────────────────────
    print("Fitting baselines …")
    freq_bl = FrequencyBaseline()
    freq_bl.fit(all_token_lists)
    lig_bl = LigandFreqBaseline()
    lig_bl.fit(all_token_lists)
    rand_bl = RandomBaseline(seed=42)
    rule_bl = RuleBaseline()

    baselines = {
        "random": rand_bl,
        "frequency": freq_bl,
        "ligand_freq": lig_bl,
        "rule": rule_bl,
    }

    # ── 5. Run benchmark ─────────────────────────────────────────
    print(f"\nRunning benchmark (decoys_per_type={args.decoys_per_type}, "
          f"pll_mode={args.pll_mode}) …")
    t0 = time.time()
    results = run_benchmark(
        model=model,
        tokenizer=tokenizer,
        test_token_lists=all_token_lists,
        test_ids=all_ids,
        pool=pool,
        device=args.device,
        n_decoys_per_type=args.decoys_per_type,
        pll_mode=args.pll_mode,
        baselines=baselines,
        max_complexes=args.max_complexes,
        verbose=True,
    )
    elapsed = time.time() - t0
    print(f"  Done in {elapsed:.1f}s")

    summary = results["summary"]
    per_complex = results["per_complex"]
    by_decoy_type = results["by_decoy_type"]
    fieldwise = results["fieldwise_examples"]
    anomalies = results["anomalies"]

    # ── 6. Write outputs ─────────────────────────────────────────
    print("\nWriting outputs …")

    # 6a. compatibility_summary.csv
    path = os.path.join(args.out, "compatibility_summary.csv")
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["metric", "value"])
        for k, v in sorted(summary.items()):
            w.writerow([k, v])
    print(f"  {path}")

    # 6b. compatibility_by_decoy_type.csv
    path = os.path.join(args.out, "compatibility_by_decoy_type.csv")
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=[
            "decoy_type", "n_real", "n_decoy", "detection_auroc",
            "mean_real", "mean_decoy",
        ])
        w.writeheader()
        for dtype, d in sorted(by_decoy_type.items()):
            w.writerow({"decoy_type": dtype, **d})
    print(f"  {path}")

    # 6c. compatibility_by_cn.csv
    cn_groups = defaultdict(lambda: {"top1": 0, "mrr": 0.0, "n": 0, "delta": 0.0})
    for r in per_complex:
        cn = r["CN"]
        cn_groups[cn]["top1"] += r["top1_real"]
        cn_groups[cn]["mrr"] += r["mrr"]
        cn_groups[cn]["delta"] += r["delta_score"]
        cn_groups[cn]["n"] += 1

    path = os.path.join(args.out, "compatibility_by_cn.csv")
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["CN", "n", "top1_acc", "mrr", "mean_delta"])
        w.writeheader()
        for cn in sorted(cn_groups.keys()):
            d = cn_groups[cn]
            n = d["n"]
            w.writerow({
                "CN": cn,
                "n": n,
                "top1_acc": round(d["top1"] / n, 4),
                "mrr": round(d["mrr"] / n, 4),
                "mean_delta": round(d["delta"] / n, 4),
            })
    print(f"  {path}")

    # 6d. compatibility_by_shape_boundary.csv (stereo presence)
    stereo_groups = defaultdict(lambda: {"top1": 0, "mrr": 0.0, "n": 0})
    for r in per_complex:
        key = "has_stereo" if r.get("has_stereo") else "no_stereo"
        stereo_groups[key]["top1"] += r["top1_real"]
        stereo_groups[key]["mrr"] += r["mrr"]
        stereo_groups[key]["n"] += 1

    path = os.path.join(args.out, "compatibility_by_shape_boundary.csv")
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["group", "n", "top1_acc", "mrr"])
        w.writeheader()
        for grp in sorted(stereo_groups.keys()):
            d = stereo_groups[grp]
            n = d["n"]
            w.writerow({
                "group": grp,
                "n": n,
                "top1_acc": round(d["top1"] / n, 4),
                "mrr": round(d["mrr"] / n, 4),
            })
    print(f"  {path}")

    # 6e. top_false_positive_decoys.csv (decoys that scored higher than real)
    fp_rows = []
    for r in per_complex:
        if r["real_rank"] > 1:
            fp_rows.append(r)
    fp_rows.sort(key=lambda r: r["delta_score"])  # most negative delta first
    path = os.path.join(args.out, "top_false_positive_decoys.csv")
    with open(path, "w", newline="") as f:
        cols = ["complex_id", "metal", "CN", "metal_row", "real_score",
                "best_decoy_score", "delta_score", "real_rank", "n_decoys"]
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in fp_rows[:100]:
            w.writerow(r)
    print(f"  {path}")

    # 6f. top_anomaly_real_records.csv
    path = os.path.join(args.out, "top_anomaly_real_records.csv")
    with open(path, "w", newline="") as f:
        cols = ["complex_id", "metal", "CN", "metal_row", "n_ligands",
                "has_stereo", "real_score", "delta_score", "real_rank"]
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in anomalies[:50]:
            w.writerow(r)
    print(f"  {path}")

    # 6g. fieldwise_score_examples.jsonl
    path = os.path.join(args.out, "fieldwise_score_examples.jsonl")
    with open(path, "w") as f:
        for fw in fieldwise:
            f.write(json.dumps(fw) + "\n")
    print(f"  {path}")

    # 6h. summary.json
    summary_json = dict(summary)
    # Add key comparisons
    mlm_top1 = summary.get("real_vs_decoys_top1", 0)
    for bl_name in baselines:
        bl_top1 = summary.get(f"bl_{bl_name}_top1", 0)
        summary_json[f"score_gain_vs_{bl_name}_baseline"] = round(mlm_top1 - bl_top1, 4)

    path = os.path.join(args.out, "summary.json")
    with open(path, "w") as f:
        json.dump(summary_json, f, indent=2)
    print(f"  {path}")

    # ── 7. Print key numbers ─────────────────────────────────────
    print(f"\n{'='*60}")
    print("KEY NUMBERS")
    print(f"{'='*60}")
    print(f"  Complexes evaluated:        {summary.get('n_complexes', 0)}")
    print(f"  Total decoys:               {summary.get('n_decoys_total', 0)}")
    print(f"  MLM Top-1 real recovery:    {summary.get('real_vs_decoys_top1', 0):.1%}")
    print(f"  MLM MRR:                    {summary.get('real_vs_decoys_mrr', 0):.3f}")
    print(f"  MLM AUROC (real vs decoy):  {summary.get('real_vs_decoy_auroc', 0):.3f}")
    print(f"  Mean ΔScore (real-best):    {summary.get('mean_delta_score', 0):.2f}")
    print()
    print("  Per decoy-type detection AUROC:")
    for dtype in DECOY_TYPES:
        auc = summary.get(f"{dtype}_detection_auc", 0)
        print(f"    {dtype:<30s} {auc:.3f}")
    print()
    print("  Baselines (Top-1):")
    for bl_name in baselines:
        bl_top1 = summary.get(f"bl_{bl_name}_top1", 0)
        gain = mlm_top1 - bl_top1
        print(f"    {bl_name:<20s} {bl_top1:.1%}  (MLM gain: {gain:+.1%})")
    print()


if __name__ == "__main__":
    main()
