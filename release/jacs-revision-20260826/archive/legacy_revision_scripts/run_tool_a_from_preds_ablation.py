#!/usr/bin/env python3
"""
run_tool_a_from_preds_ablation.py
=================================

Read existing fig5_preds.jsonl and generate per-donor prediction rows
compatible with the ablation framework.

For each masked donor outputs:
    complex_id, ligand_id, donor_position, true_token, top1_pred,
    top5_pred, probability, ablation_mode, denticity, CN, metal,
    ligand_smiles_available, geometry_available, stereo_available

This is the lightweight "no-GPU" path: it reads the pre-computed
full_context predictions and computes baselines for other modes.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent))

from coordrep_tools.tool_a_baselines import (
    CondFreqBaseline,
    LigandFreqBaseline,
    ContextFreqBaseline,
    evaluate_baseline,
    _extract_ligand_smiles,
    DONOR_ATOMS,
)


_METAL_RE = re.compile(r'\[Metal:(\w+)\|.*?CN:(\d+)\]')
_3D = {"Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn"}
_4D = {"Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd"}
_5D = {"La", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg"}


def _metal_row(m):
    if m in _3D: return "3d"
    if m in _4D: return "4d"
    if m in _5D: return "5d"
    return "other"


def _get_complex_id(sid):
    m = re.match(r'^(tmqm_\d+)_donor_\d+$', sid)
    return m.group(1) if m else sid


def _get_ligand_id(tokens, mask_pos):
    if not isinstance(mask_pos, int):
        return "?"
    for i in range(mask_pos, -1, -1):
        if re.match(r'^L\d+$', tokens[i]):
            return tokens[i]
        if tokens[i] == "L" and i + 1 < len(tokens) and tokens[i + 1].isdigit():
            return f"L{tokens[i + 1]}"
    return "?"


def _build_denticity_map(tasks_by_id):
    """Build a map from sample_id -> denticity by grouping by complex+ligand."""
    from collections import defaultdict
    complex_ligs = defaultdict(lambda: defaultdict(list))
    for sid, t in tasks_by_id.items():
        cid = _get_complex_id(sid)
        lid = _get_ligand_id(t.get("tokens", []), t.get("mask_pos"))
        complex_ligs[cid][lid].append(sid)

    dent_map = {}
    for cid, ligs in complex_ligs.items():
        for lid, sids in ligs.items():
            d = len(sids)
            label = f"dent={d}" if d <= 2 else "dent>=3"
            for sid in sids:
                dent_map[sid] = label
    return dent_map


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--preds",
                        default="/data/CoordRep/CoordSMILES/libcoordrep/figures/fig5/data/fig5_preds.jsonl")
    parser.add_argument("--tasks",
                        default="/data/CoordRep/CoordSMILES/libcoordrep/figures/fig5/data/fig5_tasks.jsonl")
    parser.add_argument("--out",
                        default="/data/CoordRep/coordrep-release/libcoordrep/outputs/tool_a_ablation")
    args = parser.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    # Load tasks (for tokens)
    print("Loading tasks …")
    tasks_by_id = {}
    with open(args.tasks) as f:
        for line in f:
            rec = json.loads(line)
            if rec["task"] == "mask_donor":
                tasks_by_id[rec["id"]] = rec

    # Load preds
    print("Loading predictions …")
    preds = []
    with open(args.preds) as f:
        for line in f:
            rec = json.loads(line)
            if rec["task"] == "mask_donor":
                preds.append(rec)

    print(f"  {len(preds)} donor predictions, {len(tasks_by_id)} tasks")

    # Build denticity map
    dent_map = _build_denticity_map(tasks_by_id)

    # Build per-prediction rows
    rows = []
    for p in preds:
        sid = p["id"]
        task = tasks_by_id.get(sid, {})
        tokens = task.get("tokens", [])
        mask_pos = task.get("mask_pos")

        donor = p["y_true"].get("donor_atom", "?")
        g = p.get("group", {})
        metal = g.get("metal_element", "?")
        cn = g.get("cn", 0)

        topk = p["pred"]["topk"]
        top1_tok = topk[0][0] if topk else "?"
        top1_prob = topk[0][1] if topk else 0
        top5_toks = [t for t, _ in topk[:5]]

        dent = dent_map.get(sid, "?")
        mrow = _metal_row(metal)

        lig_smi = _extract_ligand_smiles(tokens, mask_pos) if isinstance(mask_pos, int) else "?"

        # Ligand ID
        lig_id = _get_ligand_id(tokens, mask_pos)

        rows.append({
            "complex_id": sid,
            "ligand_id": lig_id,
            "donor_position": mask_pos,
            "true_token": donor,
            "top1_pred": top1_tok,
            "top5_pred": ";".join(top5_toks),
            "probability": round(top1_prob, 6),
            "ablation_mode": "full_context",
            "denticity": dent,
            "CN": cn,
            "metal": metal,
            "metal_row": mrow,
            "ligand_smiles": lig_smi,
            "ligand_smiles_available": True,
            "geometry_available": True,
            "stereo_available": True,
            "correct_top1": p.get("correct_top1", top1_tok == donor),
            "correct_top5": p.get("correct_top5", donor in top5_toks),
        })

    df = pd.DataFrame(rows)
    df.to_csv(out / "tool_a_per_prediction.csv", index=False)
    print(f"  Saved {len(rows)} rows to tool_a_per_prediction.csv")

    # ── Summary stats ─────────────────────────────────────────────
    n = len(df)
    top1 = df["correct_top1"].mean()
    top5 = df["correct_top5"].mean()
    print(f"\n  Full-context: Top-1 {top1:.1%}  Top-5 {top5:.1%}  n={n}")

    # ── Baselines ─────────────────────────────────────────────────
    task_samples = list(tasks_by_id.values())

    cond_bl = CondFreqBaseline()
    cond_bl.fit(task_samples)
    cond_res = evaluate_baseline(cond_bl, task_samples)

    lig_bl = LigandFreqBaseline()
    lig_bl.fit(task_samples)
    lig_res = evaluate_baseline(lig_bl, task_samples)

    ctx_bl = ContextFreqBaseline()
    ctx_bl.fit(task_samples)
    ctx_res = evaluate_baseline(ctx_bl, task_samples)

    print(f"  CondFreq:    Top-1 {cond_res['top1_acc']:.1%}")
    print(f"  LigandFreq:  Top-1 {lig_res['top1_acc']:.1%}")
    print(f"  ContextFreq: Top-1 {ctx_res['top1_acc']:.1%}")

    bl_df = pd.DataFrame([
        {"method": "CondFreq(metal,CN)", "top1_acc": round(cond_res["top1_acc"], 4),
         "top5_acc": round(cond_res["top5_acc"], 4), "mrr": round(cond_res["mrr"], 4)},
        {"method": "LigandFreq(SMILES)", "top1_acc": round(lig_res["top1_acc"], 4),
         "top5_acc": round(lig_res["top5_acc"], 4), "mrr": round(lig_res["mrr"], 4)},
        {"method": "ContextFreq(M,CN,shape,stereo)", "top1_acc": round(ctx_res["top1_acc"], 4),
         "top5_acc": round(ctx_res["top5_acc"], 4), "mrr": round(ctx_res["mrr"], 4)},
    ])
    bl_df.to_csv(out / "tool_a_baseline_comparison.csv", index=False)

    # ── Stratified ────────────────────────────────────────────────
    for col, label in [("denticity", "denticity"), ("CN", "cn"),
                       ("metal_row", "metal_row"),
                       ("true_token", "donor_element")]:
        if col in df.columns:
            g = df.groupby(col).agg(
                n=("correct_top1", "count"),
                top1_acc=("correct_top1", "mean"),
                top5_acc=("correct_top5", "mean"),
            ).reset_index()
            g.to_csv(out / f"tool_a_ablation_by_{label}.csv", index=False)

    # ── summary.json ──────────────────────────────────────────────
    summary = {
        "full_context_top1": round(top1, 4),
        "full_context_top5": round(top5, 4),
        "condfreq_top1": round(cond_res["top1_acc"], 4),
        "ligandfreq_top1": round(lig_res["top1_acc"], 4),
        "contextfreq_top1": round(ctx_res["top1_acc"], 4),
        "context_gain_full_vs_condfreq": round(top1 - cond_res["top1_acc"], 4),
        "ligand_freq_coverage": round(lig_bl.get_ligand_coverage(task_samples), 4),
    }

    # Denticity breakdown
    for dent in ["dent=1", "dent=2", "dent>=3"]:
        sub = df[df["denticity"] == dent]
        if len(sub) > 0:
            tag = dent.replace("=", "").replace(">=", "ge")
            summary[f"denticity_{tag}_top1"] = round(sub["correct_top1"].mean(), 4)

    with open(out / "summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print(f"\n  All outputs → {out}")


if __name__ == "__main__":
    main()
