#!/usr/bin/env python3
"""
eval_ranker_minphys.py
======================
Evaluate 4 ranker variants and produce MinPhys diagnostic comparison.

Variants:
  1. ranker_base          — existing finetuned ranker (CoordRep only)
  2. minphys_only         — only MinPhys aux tokens
  3. ranker_minphys       — CoordRep + MinPhys aux tokens
  4. shuffled_minphys     — CoordRep + shuffled MinPhys aux tokens

Outputs stratified metrics + co_ligand focus analysis.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import random
import re
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import torch

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from coordrep_tools.coordrep_ranker import load_ranker, CoordRepRanker
from coordrep_tools.ranker_dataset import (
    _extract_unique_complexes, split_complexes, _tokens_to_ids,
)
from coordrep_tools.hard_decoys import (
    HardDecoyPool, generate_all_hard_decoys, _parse_metal_token,
)
from coordrep_tools.ablation_masking import _detect_blocks
from coordrep_tools.minphys_features import MinPhysBinner
from coordrep_tools.minphys_tokens import (
    make_minphys_tokens, register_minphys_tokens_in_tokenizer,
)
from brain.tokenizer import CoordRepTokenizer, TokenizerConfig


# ── AUROC ─────────────────────────────────────────────────────────

def _auroc(labels, scores):
    if len(set(labels)) < 2:
        return 0.5
    pairs = list(zip(scores, labels))
    pairs.sort(key=lambda x: -x[0])
    tp, fp, prev_tp, prev_fp = 0, 0, 0, 0
    auc = 0.0
    prev_score = None
    for s, l in pairs:
        if prev_score is not None and s != prev_score:
            auc += (fp - prev_fp) * (tp + prev_tp) / 2
            prev_tp, prev_fp = tp, fp
        prev_score = s
        if l == 1:
            tp += 1
        else:
            fp += 1
    auc += (fp - prev_fp) * (tp + prev_tp) / 2
    return auc / max(tp * fp, 1)


def _metal_row(m):
    if m in {"Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn"}: return "3d"
    if m in {"Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd"}: return "4d"
    if m in {"La","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg"}: return "5d"
    return "other"


# ── main ──────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ranker_base",
                        default="inputs/checkpoints/coordrep_ranker/best_finetuned.pt")
    parser.add_argument("--minphys_only",
                        default="inputs/checkpoints/ranker_minphys/best_minphys_only.pt")
    parser.add_argument("--ranker_minphys",
                        default="inputs/checkpoints/ranker_minphys/best_ranker_minphys.pt")
    parser.add_argument("--shuffled_minphys",
                        default="inputs/checkpoints/ranker_minphys/best_shuffled_minphys.pt")
    parser.add_argument("--mlm_checkpoint",
                        default="inputs/checkpoints/pretrain_v3/best_model.pt")
    parser.add_argument("--data",
                        default="inputs/property_benchmarks/fig5_tasks.jsonl")
    parser.add_argument("--split_ids",
                        default="inputs/checkpoints/coordrep_ranker/split_ids.json")
    parser.add_argument("--binner_json",
                        default="inputs/checkpoints/ranker_minphys/binner.json")
    parser.add_argument("--out",
                        default="outputs/ranker_minphys")
    parser.add_argument("--decoys_per_real", type=int, default=20)
    parser.add_argument("--device", default="cuda:0")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    os.makedirs(args.out, exist_ok=True)

    print(f"{'='*60}")
    print("MinPhys Diagnostic Evaluation")
    print(f"{'='*60}")

    # ── 1. Load data & split ──────────────────────────────────
    print("\nLoading data …")
    complexes = _extract_unique_complexes(args.data)
    with open(args.split_ids) as f:
        split = json.load(f)
    id_to_cx = {c["id"]: c for c in complexes}
    train_cx = [id_to_cx[i] for i in split["train"] if i in id_to_cx]
    test_cx = [id_to_cx[i] for i in split["test"] if i in id_to_cx]
    print(f"  train={len(train_cx)}, test={len(test_cx)}")

    # ── 2. Build pool + binner ────────────────────────────────
    print("Building hard decoy pool (train split) …")
    train_toks = [c["tokens"] for c in train_cx]
    pool = HardDecoyPool(train_toks, seed=args.seed)

    binner = MinPhysBinner()
    if os.path.exists(args.binner_json):
        with open(args.binner_json) as f:
            bi = json.load(f)
        binner.ha_q33 = bi["ha_q33"]
        binner.ha_q66 = bi["ha_q66"]
        binner._fitted = True
        print(f"  Loaded binner: q33={binner.ha_q33}, q66={binner.ha_q66}")
    else:
        binner.fit(train_toks)
        print(f"  Fitted binner: q33={binner.ha_q33}, q66={binner.ha_q66}")

    # ── 3. Load tokenizer ─────────────────────────────────────
    tok_path = Path(args.mlm_checkpoint).parent / "tokenizer.json"
    with open(tok_path) as f:
        tok_data = json.load(f)
    tok_config = TokenizerConfig(**tok_data.get("config", {}))
    tokenizer = CoordRepTokenizer(tok_config)
    tokenizer.token2id = tok_data["token2id"]
    tokenizer.id2token = {v: k for k, v in tokenizer.token2id.items()}
    tokenizer.next_id = max(tokenizer.id2token.keys()) + 1
    register_minphys_tokens_in_tokenizer(tokenizer, train_toks, binner)

    pad_id = tokenizer.token2id.get("[PAD]", 0)

    def _pad(ids, max_len=512):
        ids = ids[:max_len]
        att = [1] * len(ids) + [0] * (max_len - len(ids))
        ids = ids + [pad_id] * (max_len - len(ids))
        return ids, att

    # ── 4. Load models ────────────────────────────────────────
    models = {}
    model_paths = {
        "ranker_base": args.ranker_base,
        "minphys_only": args.minphys_only,
        "ranker_minphys": args.ranker_minphys,
        "shuffled_minphys": args.shuffled_minphys,
    }
    for name, path in model_paths.items():
        if path and os.path.exists(path):
            r, _ = load_ranker(path, str(tok_path), device=args.device)
            models[name] = r
            print(f"  Loaded {name} from {path}")
        else:
            print(f"  SKIP {name} (not found: {path})")

    # ── 5. Generate test decoys ───────────────────────────────
    print("\nGenerating test decoys …")
    t0 = time.time()
    rng = random.Random(args.seed)
    test_data = []
    for c in test_cx:
        toks = c["tokens"]
        decoys = generate_all_hard_decoys(toks, pool, n_per_type=max(1, args.decoys_per_real // 5))
        if not decoys:
            continue
        if len(decoys) > args.decoys_per_real:
            rng.shuffle(decoys)
            decoys = decoys[:args.decoys_per_real]
        test_data.append({"id": c["id"], "tokens": toks, "decoys": decoys})
    print(f"  {len(test_data)} test complexes ({time.time()-t0:.1f}s)")

    # ── 6. Score all models ───────────────────────────────────
    print("\nScoring …")

    # token modes per model
    def _make_ids(toks, model_name, shuffle_aux=None):
        """Produce token IDs for a given model variant."""
        if model_name == "minphys_only":
            aux = make_minphys_tokens(toks, binner)
            return _tokens_to_ids(aux, tokenizer)
        elif model_name == "ranker_minphys":
            aux = make_minphys_tokens(toks, binner)
            return _tokens_to_ids(list(toks) + aux, tokenizer)
        elif model_name == "shuffled_minphys":
            aux = shuffle_aux if shuffle_aux else make_minphys_tokens(toks, binner)
            return _tokens_to_ids(list(toks) + aux, tokenizer)
        else:  # ranker_base
            return _tokens_to_ids(toks, tokenizer)

    dtype_scores = defaultdict(lambda: defaultdict(lambda: {"real": [], "decoy": []}))
    all_records = []

    for ti, td in enumerate(test_data):
        if ti % 200 == 0:
            print(f"  [{ti}/{len(test_data)}]", flush=True)

        toks = td["tokens"]
        decoys = td["decoys"]
        cid = td["id"]

        blocks = _detect_blocks(toks)
        metal, cn = "?", 0
        if blocks["metal_range"] is not None:
            p = _parse_metal_token(toks[blocks["metal_range"][0]])
            if p:
                metal, cn = p

        rec = {"complex_id": cid, "metal": metal, "CN": cn,
               "metal_row": _metal_row(metal), "n_decoys": len(decoys)}

        # For shuffled: use a random other test complex's aux
        shuffle_aux = None
        if "shuffled_minphys" in models:
            other = test_data[rng.randint(0, len(test_data) - 1)]
            shuffle_aux = make_minphys_tokens(other["tokens"], binner)

        for mname, ranker in models.items():
            ranker.eval()

            real_ids = _make_ids(toks, mname, shuffle_aux)
            rp, ra = _pad(real_ids)

            all_ids = [rp]
            all_att = [ra]
            for dtoks, reason in decoys:
                if mname == "shuffled_minphys":
                    d_ids = _make_ids(dtoks, mname, shuffle_aux)
                else:
                    d_ids = _make_ids(dtoks, mname)
                dp, da = _pad(d_ids)
                all_ids.append(dp)
                all_att.append(da)

            with torch.no_grad():
                ids_t = torch.tensor(all_ids, device=args.device).unsqueeze(0)
                att_t = torch.tensor(all_att, device=args.device).unsqueeze(0)
                scores = ranker.score_list(ids_t, att_t).squeeze(0)

            real_s = scores[0].item()
            decoy_s = scores[1:].cpu().numpy()

            rank = 1 + sum(1 for d in decoy_s if d >= real_s)
            wins = sum(1 for d in decoy_s if d < real_s)
            rec[f"{mname}_top1"] = int(rank == 1)
            rec[f"{mname}_mrr"] = round(1.0 / rank, 4)
            rec[f"{mname}_win_rate"] = round(wins / max(len(decoy_s), 1), 4)
            rec[f"{mname}_margin"] = round(real_s - float(max(decoy_s)) if len(decoy_s) > 0 else 0, 4)

            for di, (dtoks, reason) in enumerate(decoys):
                dtype = reason.split(":")[0]
                dtype_scores[dtype][mname]["real"].append(real_s)
                dtype_scores[dtype][mname]["decoy"].append(float(decoy_s[di]))

        all_records.append(rec)

    # ── 7. Aggregate ──────────────────────────────────────────
    print(f"\nAggregating {len(all_records)} complexes …")
    n = len(all_records)
    all_methods = list(models.keys())

    summary = {"n_test": n, "n_decoys_per": args.decoys_per_real}

    for m in all_methods:
        t1 = np.mean([r.get(f"{m}_top1", 0) for r in all_records])
        mrr = np.mean([r.get(f"{m}_mrr", 0) for r in all_records])
        wr = np.mean([r.get(f"{m}_win_rate", 0) for r in all_records])
        summary[f"{m}_top1"] = round(t1, 4)
        summary[f"{m}_mrr"] = round(mrr, 4)
        summary[f"{m}_win_rate"] = round(wr, 4)

    # Per-decoy-type AUROC
    by_decoy_type = []
    for dtype in sorted(dtype_scores.keys()):
        row = {"decoy_type": dtype}
        for m in all_methods:
            if m in dtype_scores[dtype]:
                d = dtype_scores[dtype][m]
                labels = [1] * len(d["real"]) + [0] * len(d["decoy"])
                scores = d["real"] + d["decoy"]
                auc = _auroc(labels, scores)
                row[f"{m}_auroc"] = round(auc, 4)
                summary[f"{dtype}_{m}_auroc"] = round(auc, 4)
        by_decoy_type.append(row)

    # Overall AUROC
    for m in all_methods:
        all_r, all_d = [], []
        for dtype in dtype_scores:
            if m in dtype_scores[dtype]:
                all_r.extend(dtype_scores[dtype][m]["real"])
                all_d.extend(dtype_scores[dtype][m]["decoy"])
        if all_r and all_d:
            labels = [1]*len(all_r) + [0]*len(all_d)
            scores = all_r + all_d
            summary[f"{m}_auroc"] = round(_auroc(labels, scores), 4)

    # ── 8. Write outputs ──────────────────────────────────────
    print("\nWriting outputs …")

    # summary.json
    p = os.path.join(args.out, "summary.json")
    with open(p, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"  {p}")

    # minphys_summary.csv
    p = os.path.join(args.out, "minphys_summary.csv")
    with open(p, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["metric", "value"])
        for k, v in sorted(summary.items()):
            w.writerow([k, v])
    print(f"  {p}")

    # minphys_by_decoy_type.csv
    if by_decoy_type:
        p = os.path.join(args.out, "minphys_by_decoy_type.csv")
        cols = sorted(by_decoy_type[0].keys())
        with open(p, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=cols)
            w.writeheader()
            for row in by_decoy_type:
                w.writerow(row)
        print(f"  {p}")

    # minphys_ablation.csv — side-by-side comparison
    p = os.path.join(args.out, "minphys_ablation.csv")
    with open(p, "w", newline="") as f:
        w = csv.writer(f)
        header = ["method", "top1", "mrr", "win_rate", "auroc"]
        for dtype in sorted(dtype_scores.keys()):
            header.append(f"{dtype}_auroc")
        w.writerow(header)
        for m in all_methods:
            row = [m,
                   summary.get(f"{m}_top1", ""),
                   summary.get(f"{m}_mrr", ""),
                   summary.get(f"{m}_win_rate", ""),
                   summary.get(f"{m}_auroc", "")]
            for dtype in sorted(dtype_scores.keys()):
                row.append(summary.get(f"{dtype}_{m}_auroc", ""))
            w.writerow(row)
    print(f"  {p}")

    # co_ligand_error_cases.csv
    p = os.path.join(args.out, "co_ligand_error_cases.csv")
    # Focus on co_ligand: compare ranker_base vs ranker_minphys
    co_lig_errors = []
    best_m = "ranker_minphys" if "ranker_minphys" in models else "ranker_base"
    for r in all_records:
        if r.get(f"{best_m}_top1", 1) == 0:
            co_lig_errors.append(r)
    co_lig_errors.sort(key=lambda r: r.get(f"{best_m}_margin", 0))
    with open(p, "w", newline="") as f:
        cols = ["complex_id", "metal", "CN", "metal_row"]
        for m in all_methods:
            cols.extend([f"{m}_top1", f"{m}_mrr", f"{m}_win_rate"])
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in co_lig_errors[:200]:
            w.writerow(r)
    print(f"  {p} ({len(co_lig_errors)} errors)")

    # ── 9. Print key results ──────────────────────────────────
    print(f"\n{'='*60}")
    print("MINPHYS DIAGNOSTIC — KEY NUMBERS")
    print(f"{'='*60}")
    print(f"  Test complexes: {n}\n")

    for m in all_methods:
        t1 = summary.get(f"{m}_top1", 0)
        mrr = summary.get(f"{m}_mrr", 0)
        wr = summary.get(f"{m}_win_rate", 0)
        auc = summary.get(f"{m}_auroc", "N/A")
        auc_s = f"{auc:.3f}" if isinstance(auc, float) else auc
        print(f"  {m:<25s} Top-1={t1:.1%}  MRR={mrr:.3f}  "
              f"WinRate={wr:.1%}  AUROC={auc_s}")

    print("\n  Per decoy-type AUROC:")
    for row in by_decoy_type:
        line = f"    {row['decoy_type']:<22s}"
        for m in all_methods:
            k = f"{m}_auroc"
            if k in row:
                line += f"  {m}={row[k]:.3f}"
        print(line)

    # Focus comparison
    print(f"\n  Co-ligand AUROC comparison:")
    for m in all_methods:
        k = f"co_ligand_hard_{m}_auroc"
        v = summary.get(k, "N/A")
        vs = f"{v:.3f}" if isinstance(v, float) else v
        print(f"    {m}: {vs}")

    base_co = summary.get("co_ligand_hard_ranker_base_auroc", 0.553)
    mp_co = summary.get("co_ligand_hard_ranker_minphys_auroc", 0)
    delta = mp_co - base_co if isinstance(mp_co, float) and isinstance(base_co, float) else 0
    print(f"\n  Δ(co_ligand AUROC) = {delta:+.3f}")
    if delta >= 0.05:
        print("  → MinPhys补足了CoordRep的ligand-size/charge盲区")
    elif delta > 0:
        print("  → 小幅提升，MinPhys有边际贡献")
    else:
        print("  → 无提升，co_ligand_hard可能是plausible alternatives而非错误负例")

    print()


if __name__ == "__main__":
    main()
