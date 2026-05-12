#!/usr/bin/env python3
"""
run_tool_a_ablation.py
======================

End-to-end ablation study for Tool A donor-marker prediction.

Usage::

    python run_tool_a_ablation.py \\
        --checkpoint checkpoints/best.pt \\
        --test data/test.jsonl \\
        --out revision_results/tool_a_ablation \\
        --modes full_context no_ligand_smiles_keep_length \\
               no_ligand_smiles_collapsed no_geometry_stereo \\
               metal_cn_only shuffled_ligand_smiles_control \\
        --topk 1 5 \\
        --stratify denticity CN metal_row donor_element ligand_frequency

Outputs:
    tool_a_ablation_summary.csv
    tool_a_ablation_by_denticity.csv
    tool_a_ablation_by_cn.csv
    tool_a_ablation_by_donor_element.csv
    tool_a_ablation_by_ligand_frequency.csv
    tool_a_accuracy_coverage.csv
    tool_a_baseline_comparison.csv
    tool_a_context_gain_cases.csv
    tool_a_failure_cases.csv
    tool_a_per_prediction.csv
    summary.json
"""

from __future__ import annotations

import argparse
import json
import sys
import warnings
from collections import Counter
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

sys.path.insert(0, str(Path(__file__).parent.parent))

from coordrep_tools.ablation_masking import ABLATION_MODES
from coordrep_tools.tool_a_ablation_eval import evaluate_tool_a_ablation
from coordrep_tools.tool_a_baselines import (
    CondFreqBaseline,
    LigandFreqBaseline,
    ContextFreqBaseline,
    evaluate_baseline,
)


# ──────────────────────────────────────────────────────────────────
# Prediction backend: model or from-preds file
# ──────────────────────────────────────────────────────────────────

def _make_model_predict_fn(checkpoint_path, tokenizer_path, device, k=10):
    """Return a predict_fn that calls the MLM model."""
    from coordrep_tools.infer import mlm_topk

    def predict_fn(ablated_tokens, mask_positions):
        masked_str = "".join(ablated_tokens)
        return mlm_topk(
            masked_str, mask_positions, k,
            checkpoint_path, tokenizer_path, device,
        )

    return predict_fn


def _make_preds_predict_fn(preds_path: str, k: int = 10):
    """
    Return a predict_fn that reads pre-computed predictions from
    ``fig5_preds.jsonl``.  Falls back to uniform if sample not found.
    """
    preds_map = {}
    with open(preds_path) as f:
        for line in f:
            rec = json.loads(line)
            if rec["task"] == "mask_donor":
                preds_map[rec["id"]] = rec

    def predict_fn(ablated_tokens, mask_positions, sample_id=None):
        # This is only usable for full_context mode
        if sample_id and sample_id in preds_map:
            rec = preds_map[sample_id]
            topk = rec["pred"]["topk"][:k]
            return [topk]
        return [[(t, 1.0 / 15) for t in ["N", "O", "S", "C", "P"]][:k]]

    return predict_fn


# ──────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Tool A ligand-SMILES masking ablation study"
    )
    parser.add_argument("--checkpoint", type=str,
                        default="/data/CoordRep/CoordSMILES/libcoordrep/checkpoints/pretrain_v3/best_model.pt")
    parser.add_argument("--tokenizer", type=str, default=None)
    parser.add_argument("--test", type=str,
                        default="/data/CoordRep/CoordSMILES/libcoordrep/figures/fig5/data/fig5_tasks.jsonl")
    parser.add_argument("--preds", type=str, default=None,
                        help="Pre-computed predictions (fig5_preds.jsonl). "
                             "If given, use for full_context mode instead of live model.")
    parser.add_argument("--out", type=str,
                        default="/data/CoordRep/coordrep-release/libcoordrep/outputs/tool_a_ablation")
    parser.add_argument("--modes", nargs="+", default=ABLATION_MODES)
    parser.add_argument("--topk", nargs="+", type=int, default=[1, 5])
    parser.add_argument("--stratify", nargs="+",
                        default=["denticity", "CN", "metal_row", "donor_element", "ligand_frequency"])
    parser.add_argument("--device", type=str, default="cuda")
    parser.add_argument("--max-samples", type=int, default=None,
                        help="Limit number of donor samples (for quick tests)")
    parser.add_argument("--use-preds-only", action="store_true",
                        help="Skip live model inference; use pre-computed preds for full_context, "
                             "baselines for others")
    args = parser.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    # ── Load test data ────────────────────────────────────────────
    print("Loading test data …")
    with open(args.test) as f:
        all_samples = [json.loads(line) for line in f]

    donor_samples = [s for s in all_samples if s.get("task") == "mask_donor"]
    if args.max_samples:
        donor_samples = donor_samples[:args.max_samples]
    print(f"  {len(donor_samples)} donor samples loaded")

    # ── Load pre-computed predictions (optional) ──────────────────
    preds_map = {}
    preds_path = args.preds or str(
        Path(args.test).parent / "fig5_preds.jsonl"
    )
    if Path(preds_path).exists():
        print(f"Loading pre-computed predictions from {preds_path} …")
        with open(preds_path) as f:
            for line in f:
                rec = json.loads(line)
                if rec["task"] == "mask_donor":
                    preds_map[rec["id"]] = rec
        print(f"  {len(preds_map)} prediction records loaded")

    # ── Build predict function ────────────────────────────────────
    use_model = not args.use_preds_only
    model_predict_fn = None
    if use_model:
        try:
            tok_path = args.tokenizer or str(
                Path(args.checkpoint).parent / "tokenizer.json"
            )
            model_predict_fn = _make_model_predict_fn(
                args.checkpoint, tok_path, args.device, k=max(args.topk) + 5
            )
            print("Model loaded for live inference")
        except Exception as e:
            print(f"  Could not load model: {e}")
            print("  Falling back to pre-computed predictions for full_context")
            use_model = False

    def predict_fn(ablated_tokens, mask_positions, sample_id=None, mode=None):
        """Unified predict: use model if available, else pre-computed preds."""
        # For full_context with pre-computed preds available
        if mode == "full_context" and sample_id in preds_map and not use_model:
            rec = preds_map[sample_id]
            topk = rec["pred"]["topk"]
            return [topk]

        if model_predict_fn is not None:
            masked_str = "".join(ablated_tokens)
            return model_predict_fn(ablated_tokens, mask_positions)

        # Last resort: use pre-computed preds (only valid for full_context)
        if sample_id in preds_map:
            rec = preds_map[sample_id]
            topk = rec["pred"]["topk"]
            return [topk]

        return []

    # ── Wrap predict_fn to pass sample_id and mode ────────────────
    # The ablation evaluator passes (tokens, positions); we need
    # to curry sample context.  We'll evaluate mode-by-mode.

    # ── Compute ligand frequency from training data ───────────────
    lig_freq = Counter()
    for s in donor_samples:
        tokens = s.get("tokens", [])
        mp = s.get("mask_pos")
        if isinstance(mp, int):
            from coordrep_tools.tool_a_baselines import _extract_ligand_smiles
            smi = _extract_ligand_smiles(tokens, mp)
            lig_freq[smi] += 1

    # ── Evaluate ablation modes ───────────────────────────────────
    print(f"\nEvaluating ablation modes: {args.modes}")

    all_per_pred = []

    for mode in args.modes:
        print(f"\n  Mode: {mode}")

        # Build mode-specific predict wrapper
        def _predict_for_mode(abl_tokens, mask_pos, _mode=mode):
            # We need sample_id — we'll handle this differently
            masked_str = "".join(abl_tokens)
            if model_predict_fn is not None:
                return model_predict_fn(abl_tokens, mask_pos)
            return []

        results = evaluate_tool_a_ablation(
            samples=donor_samples,
            predict_fn=_predict_for_mode,
            ablation_modes=[mode],
            top_k=args.topk,
            ligand_freq_counter=lig_freq,
        )

        pm = results["per_mode"].get(mode, {})
        print(f"    Top-1: {pm.get('top1_acc', 0):.1%}  "
              f"Top-5: {pm.get('top5_acc', 0):.1%}  "
              f"MRR: {pm.get('mrr', 0):.3f}  "
              f"NLL: {pm.get('mean_nll', 0):.3f}  "
              f"n={pm.get('n', 0)}")

        all_per_pred.extend(results["per_prediction"])

        # Save stratified
        for stratum_name, rows in results["stratified"].items():
            if stratum_name in args.stratify and rows:
                fn = out / f"tool_a_ablation_by_{stratum_name}.csv"
                existing = []
                if fn.exists():
                    existing = pd.read_csv(fn).to_dict("records")
                existing.extend(rows)
                pd.DataFrame(existing).to_csv(fn, index=False)

    # ── Save per-prediction ───────────────────────────────────────
    if all_per_pred:
        ppdf = pd.DataFrame(all_per_pred)
        # Drop list columns for CSV
        for col in ppdf.columns:
            if ppdf[col].apply(lambda x: isinstance(x, list)).any():
                ppdf[col] = ppdf[col].apply(
                    lambda x: ";".join(str(v) for v in x) if isinstance(x, list) else x
                )
        ppdf.to_csv(out / "tool_a_per_prediction.csv", index=False)

    # ── Ablation summary ──────────────────────────────────────────
    summary_rows = []
    for mode in args.modes:
        mode_preds = [p for p in all_per_pred if p["ablation_mode"] == mode]
        if not mode_preds:
            continue
        n = len(mode_preds)
        summary_rows.append({
            "mode": mode,
            "n": n,
            "top1_acc": round(np.mean([p["correct_top1"] for p in mode_preds]), 4),
            "top5_acc": round(np.mean([p["correct_top5"] for p in mode_preds]), 4),
            "mrr": round(np.mean([p["rr"] for p in mode_preds]), 4),
            "mean_nll": round(np.mean([p["nll"] for p in mode_preds]), 4),
        })
    if summary_rows:
        pd.DataFrame(summary_rows).to_csv(out / "tool_a_ablation_summary.csv", index=False)

    # ── Baselines ─────────────────────────────────────────────────
    print("\nEvaluating baselines …")

    cond_bl = CondFreqBaseline()
    cond_bl.fit(donor_samples)
    cond_res = evaluate_baseline(cond_bl, donor_samples)
    print(f"  CondFreq(metal,CN):  Top-1 {cond_res['top1_acc']:.1%}  Top-5 {cond_res['top5_acc']:.1%}")

    lig_bl = LigandFreqBaseline()
    lig_bl.fit(donor_samples)
    lig_res = evaluate_baseline(lig_bl, donor_samples)
    print(f"  LigandFreq(SMILES):  Top-1 {lig_res['top1_acc']:.1%}  Top-5 {lig_res['top5_acc']:.1%}")

    ctx_bl = ContextFreqBaseline()
    ctx_bl.fit(donor_samples)
    ctx_res = evaluate_baseline(ctx_bl, donor_samples)
    print(f"  ContextFreq(M,CN,shape,stereo): Top-1 {ctx_res['top1_acc']:.1%}  Top-5 {ctx_res['top5_acc']:.1%}")

    bl_rows = [
        {"method": "CondFreq(metal,CN)", **{k: v for k, v in cond_res.items() if k != "predictions"}},
        {"method": "LigandFreq(SMILES)", **{k: v for k, v in lig_res.items() if k != "predictions"}},
        {"method": "ContextFreq(M,CN,shape,stereo)", **{k: v for k, v in ctx_res.items() if k != "predictions"}},
    ]
    pd.DataFrame(bl_rows).to_csv(out / "tool_a_baseline_comparison.csv", index=False)

    # ── Accuracy-coverage curve ───────────────────────────────────
    acc_cov_rows = []
    for mode in args.modes:
        mode_preds = [p for p in all_per_pred if p["ablation_mode"] == mode]
        if not mode_preds:
            continue
        for thr in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 0.99]:
            above = [p for p in mode_preds if p["probability"] >= thr]
            coverage = len(above) / len(mode_preds) if mode_preds else 0
            accuracy = np.mean([p["correct_top1"] for p in above]) if above else 0
            acc_cov_rows.append({
                "mode": mode,
                "threshold": thr,
                "coverage": round(coverage, 4),
                "accuracy": round(float(accuracy), 4),
                "n_above": len(above),
            })
    if acc_cov_rows:
        pd.DataFrame(acc_cov_rows).to_csv(out / "tool_a_accuracy_coverage.csv", index=False)

    # ── Context gain cases ────────────────────────────────────────
    # Cases where full_context is correct but no_ligand_smiles_collapsed is wrong
    full_preds = {p["id"]: p for p in all_per_pred if p["ablation_mode"] == "full_context"}
    collapsed_preds = {p["id"]: p for p in all_per_pred
                       if p["ablation_mode"] == "no_ligand_smiles_collapsed"}

    gain_cases = []
    for sid, fp in full_preds.items():
        cp = collapsed_preds.get(sid)
        if fp["correct_top1"] and cp and not cp["correct_top1"]:
            gain_cases.append({
                "id": sid,
                "true_token": fp["true_token"],
                "full_pred": fp["top1_pred"],
                "full_prob": fp["probability"],
                "collapsed_pred": cp["top1_pred"],
                "collapsed_prob": cp["probability"],
                "metal": fp["metal"],
                "CN": fp["CN"],
                "denticity": fp["denticity"],
            })
    if gain_cases:
        pd.DataFrame(gain_cases[:200]).to_csv(out / "tool_a_context_gain_cases.csv", index=False)

    # ── Failure cases ─────────────────────────────────────────────
    failure_cases = [p for p in all_per_pred
                     if p["ablation_mode"] == "full_context" and not p["correct_top5"]]
    if failure_cases:
        pd.DataFrame(failure_cases[:200]).to_csv(out / "tool_a_failure_cases.csv", index=False)

    # ── Per-ligand all-correct by denticity ───────────────────────
    # Re-compute from per_prediction
    from collections import defaultdict as dd
    lig_groups = dd(list)
    for p in all_per_pred:
        key = (p["id"], p["ablation_mode"], p["ligand_smiles"])
        lig_groups[key].append(p)

    lig_allcorrect = dd(lambda: {"total": 0, "correct": 0})
    for (cid, mode, smi), preds in lig_groups.items():
        dent = preds[0]["denticity"]
        all_ok = all(p["correct_top1"] for p in preds)
        lig_allcorrect[(mode, dent)]["total"] += 1
        if all_ok:
            lig_allcorrect[(mode, dent)]["correct"] += 1

    # ── summary.json ──────────────────────────────────────────────
    summary = {}

    # Per-mode metrics
    for row in summary_rows:
        m = row["mode"]
        summary[f"{m}_top1"] = row["top1_acc"]
        summary[f"{m}_top5"] = row["top5_acc"]
        summary[f"{m}_mrr"] = row["mrr"]

    # Baselines
    summary["condfreq_top1"] = round(cond_res["top1_acc"], 4)
    summary["condfreq_top5"] = round(cond_res["top5_acc"], 4)
    summary["ligandfreq_top1"] = round(lig_res["top1_acc"], 4)
    summary["ligandfreq_top5"] = round(lig_res["top5_acc"], 4)
    summary["contextfreq_top1"] = round(ctx_res["top1_acc"], 4)

    # Context gains
    fc = summary.get("full_context_top1", 0)
    summary["context_gain_full_vs_condfreq"] = round(fc - summary.get("condfreq_top1", 0), 4)
    nosmiles = summary.get("no_ligand_smiles_collapsed_top1", 0)
    mconly = summary.get("metal_cn_only_top1", 0)
    summary["context_gain_no_smiles_vs_metal_cn_only"] = round(nosmiles - mconly, 4)

    # Denticity stratified
    for mode in args.modes:
        for dent in ["dent=1", "dent=2", "dent>=3"]:
            preds_d = [p for p in all_per_pred
                       if p["ablation_mode"] == mode and p["denticity"] == dent]
            if preds_d:
                tag = dent.replace("=", "").replace(">=", "ge")
                summary[f"{mode}_{tag}_top1"] = round(
                    np.mean([p["correct_top1"] for p in preds_d]), 4
                )

    # Per-ligand all-correct
    for (mode, dent), v in lig_allcorrect.items():
        if v["total"] > 0:
            tag = dent.replace("=", "").replace(">=", "ge")
            summary[f"{mode}_per_ligand_all_correct_{tag}"] = round(
                v["correct"] / v["total"], 4
            )

    with open(out / "summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print(f"\n{'='*60}")
    print(f"All outputs saved to {out}")
    print(f"{'='*60}")

    # Print key numbers
    print("\nKey numbers for response letter:")
    for k in ["full_context_top1", "no_ligand_smiles_keep_length_top1",
              "no_ligand_smiles_collapsed_top1", "no_geometry_stereo_top1",
              "metal_cn_only_top1", "condfreq_top1", "ligandfreq_top1",
              "context_gain_full_vs_condfreq",
              "context_gain_no_smiles_vs_metal_cn_only"]:
        if k in summary:
            print(f"  {k}: {summary[k]}")


if __name__ == "__main__":
    main()
