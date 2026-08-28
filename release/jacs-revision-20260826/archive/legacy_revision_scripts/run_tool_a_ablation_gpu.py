#!/usr/bin/env python3
"""
run_tool_a_ablation_gpu.py
==========================
Run Tool A donor-prediction ablation study on GPU.

For every ablation mode the script:
  1. Constructs the masked/ablated token list from fig5_tasks.jsonl.
  2. Converts tokens → IDs, runs the MLM checkpoint, collects top-k.
  3. Aggregates per-mode metrics (Top-1, Top-5, MRR, NLL) and
     stratified breakdowns (denticity, CN, donor element, ligand freq).
  4. Writes CSVs + summary.json to the output directory.

Usage (inside the 1701 conda env):
  python scripts/run_tool_a_ablation_gpu.py \
      --checkpoint /data/CoordRep/CoordSMILES/libcoordrep/checkpoints/pretrain_v3/best_model.pt \
      --out revision_results/tool_a_ablation
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import sys
import time
from collections import Counter, defaultdict
from typing import Dict, List, Tuple

import numpy as np
import torch

# ── project imports ──────────────────────────────────────────────
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from coordrep_tools.infer import load_model_and_tokenizer
from coordrep_tools.ablation_masking import make_ablation_input, ABLATION_MODES
from coordrep_tools.tool_a_baselines import (
    CondFreqBaseline, LigandFreqBaseline, ContextFreqBaseline,
    _extract_ligand_smiles, _extract_shape_label, _extract_stereo_class,
)

# ── constants ────────────────────────────────────────────────────
DONOR_ATOMS = [
    "C", "N", "O", "S", "P", "F", "Cl", "Br", "I",
    "Se", "Te", "As", "Si", "B", "H",
]
TOP_K = 20

_3D = {"Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn"}
_4D = {"Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd"}
_5D = {"La", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg"}


def _metal_row(m):
    if m in _3D: return "3d"
    if m in _4D: return "4d"
    if m in _5D: return "5d"
    return "other"


def _donor_elem_group(d):
    if d in ("N", "O", "S", "P", "C"):
        return d
    if d in ("F", "Cl", "Br", "I"):
        return "halide"
    return "other"


def _cn_group(cn):
    if cn <= 4: return "CN<=4"
    if cn <= 6: return "CN=5-6"
    if cn <= 8: return "CN=7-8"
    return "CN>=9"


# ── ligand-id helpers (handle split "L","2" tokens) ─────────────
def _get_ligand_id(tokens, mask_pos):
    if not isinstance(mask_pos, int):
        return "?"
    for i in range(mask_pos, -1, -1):
        if re.match(r'^L\d+$', tokens[i]):
            return tokens[i]
        if tokens[i] == "L" and i + 1 < len(tokens) and tokens[i + 1].isdigit():
            return f"L{tokens[i + 1]}"
    return "?"


def _get_complex_id(sid):
    m = re.match(r'^(tmqm_\d+)_donor_\d+$', sid)
    return m.group(1) if m else sid


def _build_denticity_map(samples):
    complex_ligs = defaultdict(lambda: defaultdict(list))
    for s in samples:
        sid = s["id"]
        cid = _get_complex_id(sid)
        lid = _get_ligand_id(s["tokens"], s["mask_pos"])
        complex_ligs[cid][lid].append(sid)
    dent_map = {}
    for cid, ligs in complex_ligs.items():
        for lid, sids in ligs.items():
            d = len(sids)
            label = f"dent={d}" if d <= 2 else "dent>=3"
            for sid in sids:
                dent_map[sid] = label
    return dent_map


# ── token-list → model inference ─────────────────────────────────

def _tokens_to_ids(tokens, tokenizer):
    """Convert a token list to ID tensor.  Unknown tokens → [UNK]."""
    mask_id = tokenizer.token2id.get("[MASK]", 3)
    unk_id = tokenizer.token2id.get("[UNK]", 1)
    ids = []
    for t in tokens:
        if t == "[MASK]":
            ids.append(mask_id)
        elif t == "[LIGMASK]":
            # Not in vocab – map to [MASK] so the model still sees
            # a "something missing" signal at reduced sequence length.
            ids.append(mask_id)
        elif t in tokenizer.token2id:
            ids.append(tokenizer.token2id[t])
        else:
            ids.append(unk_id)
    return ids


@torch.no_grad()
def _infer_batch(
    batch_ids: List[List[int]],
    batch_positions: List[List[int]],
    model,
    device: str,
    top_k: int = TOP_K,
) -> List[List[List[Tuple[str, float]]]]:
    """
    Run batched MLM inference.
    Returns: for each sample, for each mask-pos, a top-k list [(tok, prob), …].
    """
    # Pad to same length
    max_len = min(max(len(ids) for ids in batch_ids), 512)
    padded = []
    masks = []
    for ids in batch_ids:
        truncated = ids[:max_len]
        pad_len = max_len - len(truncated)
        padded.append(truncated + [0] * pad_len)
        masks.append([1] * len(truncated) + [0] * pad_len)

    input_t = torch.tensor(padded, device=device)
    mask_t = torch.tensor(masks, device=device)
    logits, _ = model(input_t, attention_mask=mask_t)

    # extract id2token once
    id2token = {v: k for k, v in model._tokenizer_ref.token2id.items()} \
        if hasattr(model, '_tokenizer_ref') else None

    results = []
    for i, positions in enumerate(batch_positions):
        sample_results = []
        for pos in positions:
            if pos >= logits.shape[1]:
                sample_results.append([])
                continue
            probs = torch.softmax(logits[i, pos, :], dim=-1)
            topk_probs, topk_ids = torch.topk(probs, min(top_k, probs.shape[0]))
            sample_results.append([
                (int(tid), float(p))
                for p, tid in zip(topk_probs.cpu().numpy(), topk_ids.cpu().numpy())
            ])
        results.append(sample_results)
    return results


# ── main ─────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint",
                        default="/data/CoordRep/CoordSMILES/libcoordrep/checkpoints/pretrain_v3/best_model.pt")
    parser.add_argument("--tasks",
                        default="/data/CoordRep/CoordSMILES/libcoordrep/figures/fig5/data/fig5_tasks.jsonl")
    parser.add_argument("--out",
                        default="/data/CoordRep/coordrep-release/libcoordrep/revision_results/tool_a_ablation")
    parser.add_argument("--batch-size", type=int, default=64)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--max-samples", type=int, default=0,
                        help="0 = use all donor samples")
    args = parser.parse_args()

    os.makedirs(args.out, exist_ok=True)

    # ── 1. load data ─────────────────────────────────────────────
    print("Loading tasks …")
    samples = []
    with open(args.tasks) as f:
        for line in f:
            rec = json.loads(line)
            if rec["task"] == "mask_donor":
                samples.append(rec)
    if args.max_samples > 0:
        samples = samples[:args.max_samples]
    print(f"  {len(samples)} donor samples")

    # ── 2. load model + tokenizer ────────────────────────────────
    print("Loading model …")
    model, tokenizer = load_model_and_tokenizer(args.checkpoint, device=args.device)
    id2token = {}
    for tok_str, tid in tokenizer.token2id.items():
        id2token[tid] = tok_str
    print(f"  vocab={len(tokenizer.token2id)}, device={args.device}")

    # ── 3. build ancillary maps ──────────────────────────────────
    dent_map = _build_denticity_map(samples)

    # Ligand-frequency counter (for stratification)
    lig_counter = Counter()
    for s in samples:
        smi = _extract_ligand_smiles(s["tokens"], s["mask_pos"])
        lig_counter[smi] += 1

    # ── 4. baselines (fit on full data, using raw sample dicts) ──
    print("Fitting baselines …")
    cond_bl = CondFreqBaseline()
    cond_bl.fit(samples)
    lig_bl = LigandFreqBaseline()
    lig_bl.fit(samples)
    ctx_bl = ContextFreqBaseline()
    ctx_bl.fit(samples)

    # Evaluate baselines
    def _bl_acc(bl, top=1):
        correct = 0
        for s in samples:
            preds = bl.predict(s, k=max(5, top))
            tops = [t for t, _ in preds[:top]]
            target = s["y_true"].get("donor_atom", "?")
            if target in tops:
                correct += 1
        return correct / len(samples)

    def _bl_mrr(bl):
        total = 0.0
        for s in samples:
            preds = bl.predict(s, k=20)
            pred_toks = [t for t, _ in preds]
            target = s["y_true"].get("donor_atom", "?")
            if target in pred_toks:
                total += 1.0 / (pred_toks.index(target) + 1)
        return total / len(samples)

    baselines = {
        "CondFreq(metal,CN)":           {"top1": _bl_acc(cond_bl, 1), "top5": _bl_acc(cond_bl, 5), "mrr": _bl_mrr(cond_bl)},
        "LigandFreq(SMILES)":           {"top1": _bl_acc(lig_bl, 1),  "top5": _bl_acc(lig_bl, 5),  "mrr": _bl_mrr(lig_bl)},
        "ContextFreq(M,CN,shape,stereo)": {"top1": _bl_acc(ctx_bl, 1), "top5": _bl_acc(ctx_bl, 5), "mrr": _bl_mrr(ctx_bl)},
    }
    for name, m in baselines.items():
        print(f"  {name}: Top-1={m['top1']:.1%}  Top-5={m['top5']:.1%}  MRR={m['mrr']:.3f}")

    # ── 5. run each ablation mode ────────────────────────────────
    mode_rows = defaultdict(list)          # mode → list of row-dicts
    mode_metrics = {}                       # mode → aggregate dict

    for mode in ABLATION_MODES:
        t0 = time.time()
        print(f"\n{'='*60}")
        print(f"Ablation mode: {mode}")

        # Build ablated inputs
        ablated_data = []   # (sample, abl_tokens, abl_positions)
        skipped = 0
        for s in samples:
            tokens = s["tokens"]
            mp = s["mask_pos"]
            donor_positions = [mp] if isinstance(mp, int) else mp
            try:
                abl_tokens, abl_pos = make_ablation_input(
                    tokens, donor_positions, mode, shuffle_seed=42,
                )
                ablated_data.append((s, abl_tokens, abl_pos))
            except Exception as e:
                skipped += 1
                continue

        if skipped:
            print(f"  Skipped {skipped} samples due to masking errors")

        # Convert to IDs
        all_ids = [_tokens_to_ids(abl_toks, tokenizer) for _, abl_toks, _ in ablated_data]
        all_positions = [abl_pos for _, _, abl_pos in ablated_data]

        # Batched inference
        n = len(ablated_data)
        all_topk = [None] * n
        bs = args.batch_size
        for start in range(0, n, bs):
            end = min(start + bs, n)
            batch_results = _infer_batch(
                all_ids[start:end],
                all_positions[start:end],
                model, args.device, TOP_K,
            )
            for j, res in enumerate(batch_results):
                all_topk[start + j] = res

        # Collect per-prediction metrics
        correct1 = correct5 = 0
        rr_sum = 0.0
        nll_sum = 0.0
        rows = []

        for idx, (s, abl_toks, abl_pos) in enumerate(ablated_data):
            g = s.get("group", {})
            metal = g.get("metal_element", "?")
            cn = g.get("cn", 0)
            y_true = s["y_true"].get("donor_atom", "?")
            dent = dent_map.get(s["id"], "?")
            mrow = _metal_row(metal)
            delem = _donor_elem_group(y_true)
            cn_grp = _cn_group(cn)
            lig_smi = _extract_ligand_smiles(s["tokens"], s["mask_pos"])
            lig_bucket = "common" if lig_counter.get(lig_smi, 0) >= 5 else "rare"

            # Only one mask position per sample for donor prediction
            topk_raw = all_topk[idx][0] if all_topk[idx] else []
            topk_tokens = [(id2token.get(tid, f"[?{tid}]"), p) for tid, p in topk_raw]

            top1_tok = topk_tokens[0][0] if topk_tokens else "?"
            top1_prob = topk_tokens[0][1] if topk_tokens else 0.0
            top5_toks = [t for t, _ in topk_tokens[:5]]

            hit1 = int(top1_tok == y_true)
            hit5 = int(y_true in top5_toks)
            correct1 += hit1
            correct5 += hit5

            # MRR
            all_pred_toks = [t for t, _ in topk_tokens]
            if y_true in all_pred_toks:
                rank = all_pred_toks.index(y_true) + 1
                rr_sum += 1.0 / rank
            else:
                rank = -1

            # NLL: find prob of true token
            true_prob = 0.0
            for t, p in topk_tokens:
                if t == y_true:
                    true_prob = p
                    break
            nll = -math.log(max(true_prob, 1e-10))
            nll_sum += nll

            row = {
                "id": s["id"],
                "ablation_mode": mode,
                "ligand_id": _get_ligand_id(s["tokens"], s["mask_pos"]),
                "true_token": y_true,
                "top1_pred": top1_tok,
                "top5_pred": ";".join(top5_toks),
                "probability": round(top1_prob, 6),
                "true_prob": round(true_prob, 6),
                "nll": round(nll, 4),
                "rank": rank,
                "correct_top1": hit1,
                "correct_top5": hit5,
                "denticity": dent,
                "CN": cn,
                "cn_group": cn_grp,
                "metal": metal,
                "metal_row": mrow,
                "donor_element": delem,
                "ligand_smiles": lig_smi,
                "ligand_frequency_bucket": lig_bucket,
            }
            rows.append(row)

        mode_rows[mode] = rows
        n_valid = len(rows)
        acc1 = correct1 / n_valid if n_valid else 0
        acc5 = correct5 / n_valid if n_valid else 0
        mrr = rr_sum / n_valid if n_valid else 0
        mean_nll = nll_sum / n_valid if n_valid else 0

        mode_metrics[mode] = {
            "n": n_valid,
            "top1_acc": round(acc1, 4),
            "top5_acc": round(acc5, 4),
            "mrr": round(mrr, 4),
            "mean_nll": round(mean_nll, 4),
        }
        elapsed = time.time() - t0
        print(f"  n={n_valid}  Top-1={acc1:.1%}  Top-5={acc5:.1%}  "
              f"MRR={mrr:.3f}  NLL={mean_nll:.3f}  ({elapsed:.1f}s)")

    # ── 6. write outputs ─────────────────────────────────────────
    print(f"\n{'='*60}")
    print("Writing outputs …")

    # 6a. Summary CSV  (one row per mode)
    import csv
    summary_csv = os.path.join(args.out, "tool_a_ablation_summary.csv")
    with open(summary_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["ablation_mode", "n", "top1_acc", "top5_acc", "mrr", "mean_nll"])
        w.writeheader()
        for mode in ABLATION_MODES:
            w.writerow({"ablation_mode": mode, **mode_metrics[mode]})
    print(f"  {summary_csv}")

    # 6b. Baseline comparison CSV
    bl_csv = os.path.join(args.out, "tool_a_baseline_comparison.csv")
    with open(bl_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["method", "top1_acc", "top5_acc", "mrr"])
        w.writeheader()
        for name, m in baselines.items():
            w.writerow({"method": name, "top1_acc": round(m["top1"], 4),
                         "top5_acc": round(m["top5"], 4), "mrr": round(m["mrr"], 4)})
    print(f"  {bl_csv}")

    # 6c. Stratified CSVs
    def _write_stratified(filename, group_key, display_key=None):
        if display_key is None:
            display_key = group_key
        path = os.path.join(args.out, filename)
        accum = defaultdict(lambda: defaultdict(lambda: {"c1": 0, "c5": 0, "rr": 0.0, "nll": 0.0, "n": 0}))
        for mode in ABLATION_MODES:
            for r in mode_rows[mode]:
                grp = r[group_key]
                d = accum[mode][grp]
                d["c1"] += r["correct_top1"]
                d["c5"] += r["correct_top5"]
                d["rr"] += (1.0 / r["rank"]) if r["rank"] > 0 else 0.0
                d["nll"] += r["nll"]
                d["n"] += 1

        rows_out = []
        for mode in ABLATION_MODES:
            for grp in sorted(accum[mode].keys(), key=str):
                d = accum[mode][grp]
                n = d["n"]
                rows_out.append({
                    "ablation_mode": mode,
                    display_key: grp,
                    "n": n,
                    "top1_acc": round(d["c1"] / n, 4) if n else 0,
                    "top5_acc": round(d["c5"] / n, 4) if n else 0,
                    "mrr": round(d["rr"] / n, 4) if n else 0,
                    "mean_nll": round(d["nll"] / n, 4) if n else 0,
                })

        with open(path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["ablation_mode", display_key, "n",
                                               "top1_acc", "top5_acc", "mrr", "mean_nll"])
            w.writeheader()
            w.writerows(rows_out)
        print(f"  {path}")

    _write_stratified("tool_a_ablation_by_denticity.csv", "denticity")
    _write_stratified("tool_a_ablation_by_cn.csv", "CN")
    _write_stratified("tool_a_ablation_by_donor_element.csv", "donor_element")
    _write_stratified("tool_a_ablation_by_ligand_freq.csv", "ligand_frequency_bucket", "ligand_frequency")

    # 6d. Per-ligand all-correct
    lig_correct = os.path.join(args.out, "tool_a_ablation_per_ligand.csv")
    lig_rows = []
    for mode in ABLATION_MODES:
        # group by (complex_id, ligand_id)
        lig_groups = defaultdict(list)
        for r in mode_rows[mode]:
            cid = _get_complex_id(r["id"])
            lid = r["ligand_id"]
            lig_groups[(cid, lid)].append(r)
        for (cid, lid), rs in lig_groups.items():
            all_ok = all(r["correct_top1"] for r in rs)
            lig_rows.append({
                "ablation_mode": mode,
                "complex_id": cid,
                "ligand_id": lid,
                "denticity": rs[0]["denticity"],
                "n_donors": len(rs),
                "all_correct_top1": int(all_ok),
            })
    with open(lig_correct, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["ablation_mode", "complex_id", "ligand_id",
                                           "denticity", "n_donors", "all_correct_top1"])
        w.writeheader()
        w.writerows(lig_rows)
    print(f"  {lig_correct}")

    # 6e. summary.json
    summary = {}
    for mode in ABLATION_MODES:
        m = mode_metrics[mode]
        summary[f"{mode}_top1"] = m["top1_acc"]
        summary[f"{mode}_top5"] = m["top5_acc"]
        summary[f"{mode}_mrr"]  = m["mrr"]
        summary[f"{mode}_nll"]  = m["mean_nll"]
    for name, m in baselines.items():
        safe = name.replace("(", "_").replace(")", "").replace(",", "_").replace(" ", "")
        summary[f"baseline_{safe}_top1"] = round(m["top1"], 4)

    # Key comparisons
    fc = mode_metrics.get("full_context", {})
    nl = mode_metrics.get("no_ligand_smiles_keep_length", {})
    nc = mode_metrics.get("no_ligand_smiles_collapsed", {})
    ng = mode_metrics.get("no_geometry_stereo", {})
    mc = mode_metrics.get("metal_cn_only", {})
    sh = mode_metrics.get("shuffled_ligand_smiles_control", {})

    summary["delta_full_vs_no_lig_keep"]      = round(fc.get("top1_acc", 0) - nl.get("top1_acc", 0), 4)
    summary["delta_full_vs_no_lig_collapsed"]  = round(fc.get("top1_acc", 0) - nc.get("top1_acc", 0), 4)
    summary["delta_full_vs_no_geom"]           = round(fc.get("top1_acc", 0) - ng.get("top1_acc", 0), 4)
    summary["delta_full_vs_shuffled"]          = round(fc.get("top1_acc", 0) - sh.get("top1_acc", 0), 4)
    summary["delta_no_lig_keep_vs_metal_only"] = round(nl.get("top1_acc", 0) - mc.get("top1_acc", 0), 4)
    summary["delta_no_lig_coll_vs_metal_only"] = round(nc.get("top1_acc", 0) - mc.get("top1_acc", 0), 4)
    summary["baseline_LigandFreq_top1"]        = round(baselines["LigandFreq(SMILES)"]["top1"], 4)

    summary_path = os.path.join(args.out, "summary.json")
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"  {summary_path}")

    # ── 7. print key numbers ─────────────────────────────────────
    print(f"\n{'='*60}")
    print("KEY NUMBERS FOR RESPONSE LETTER")
    print(f"{'='*60}")
    print(f"{'Mode':<40s} {'Top-1':>7s} {'Top-5':>7s} {'MRR':>7s} {'NLL':>7s}")
    print("-" * 70)
    for mode in ABLATION_MODES:
        m = mode_metrics[mode]
        print(f"{mode:<40s} {m['top1_acc']:7.1%} {m['top5_acc']:7.1%} {m['mrr']:7.3f} {m['mean_nll']:7.3f}")
    print("-" * 70)
    for name, m in baselines.items():
        print(f"{name:<40s} {m['top1']:7.1%} {m['top5']:7.1%} {m['mrr']:7.3f}     -")
    print()
    print("Key deltas (Top-1):")
    print(f"  A. full_context vs LigandFreq(SMILES):      {fc.get('top1_acc',0) - baselines['LigandFreq(SMILES)']['top1']:+.1%}")
    print(f"  B. no_lig_keep_length vs metal_cn_only:      {nl.get('top1_acc',0) - mc.get('top1_acc',0):+.1%}")
    print(f"  C. no_lig_collapsed vs metal_cn_only:        {nc.get('top1_acc',0) - mc.get('top1_acc',0):+.1%}")
    print(f"  D. full_context vs no_geometry_stereo:        {fc.get('top1_acc',0) - ng.get('top1_acc',0):+.1%}")
    print(f"  E. shuffled_control vs full_context:          {sh.get('top1_acc',0) - fc.get('top1_acc',0):+.1%}")
    print()


if __name__ == "__main__":
    main()
