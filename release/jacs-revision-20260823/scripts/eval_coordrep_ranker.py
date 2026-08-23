#!/usr/bin/env python3
"""
eval_coordrep_ranker.py
=======================
Evaluate CoordRep-Ranker on held-out test complexes and compare
against all baselines.

Outputs:
  - ranker_summary.csv / summary.json
  - ranker_by_decoy_type.csv
  - ranker_by_cn.csv
  - ranker_by_boundary.csv
  - ranker_error_cases.csv
  - leakage_diagnostics (ligand set overlap check)
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
from typing import Dict, List, Optional, Tuple

import numpy as np
import torch

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from coordrep_tools.coordrep_ranker import (
    load_ranker, load_ranker_from_mlm, CoordRepRanker,
)
from coordrep_tools.infer import load_model_and_tokenizer
from coordrep_tools.ranker_dataset import (
    _extract_unique_complexes, split_complexes, _tokens_to_ids,
    RankerEvalDataset,
)
from coordrep_tools.ranker_losses import get_loss_fn
from coordrep_tools.hard_decoys import (
    HardDecoyPool, generate_all_hard_decoys, HARD_DECOY_TYPES,
    _parse_metal_token,
)
from coordrep_tools.ablation_masking import _detect_blocks
from coordrep_tools.validate import is_valid_coordrep
from coordrep_tools.coordrep_score import _SPECIAL


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


# ── Baselines ─────────────────────────────────────────────────────

class RandomBaseline:
    def __init__(self, seed=42):
        import random
        self.rng = random.Random(seed)
    def score(self, tokens):
        return self.rng.random()

class FrequencyBaseline:
    def __init__(self):
        self.freq = Counter()
    def fit(self, all_token_lists):
        for toks in all_token_lists:
            self.freq[tuple(toks)] += 1
    def score(self, tokens):
        return self.freq.get(tuple(tokens), 0)

class RuleValidatorBaseline:
    def score(self, tokens):
        s = "".join(tokens)
        return 1.0 if is_valid_coordrep(s, strict=True) else 0.0


# ── MLM PLL scorer ────────────────────────────────────────────────

class MLMFieldScorer:
    """Score using raw MLM field-targeted PLL (same as hard benchmark)."""

    def __init__(self, model, tokenizer, device, gpu_batch=64):
        self.model = model
        self.tokenizer = tokenizer
        self.device = device
        self.gpu_batch = gpu_batch
        self.pad_id = tokenizer.token2id.get("[PAD]", 0)
        self.mask_id = tokenizer.token2id.get("[MASK]", 3)

    def _forward_batch(self, id_seqs):
        if not id_seqs:
            return []
        max_l = max(len(s) for s in id_seqs)
        padded = [s + [self.pad_id] * (max_l - len(s)) for s in id_seqs]
        masks = [[1] * len(s) + [0] * (max_l - len(s)) for s in id_seqs]
        inp = torch.tensor(padded, device=self.device)
        att = torch.tensor(masks, device=self.device)
        with torch.no_grad():
            logits, _ = self.model(inp, attention_mask=att)
        return logits

    def score_pair(self, real_toks, decoy_toks):
        """Return (real_score, decoy_score) using field-targeted masking."""
        real_ids = _tokens_to_ids(real_toks, self.tokenizer)
        decoy_ids = _tokens_to_ids(decoy_toks, self.tokenizer)

        min_len = min(len(real_toks), len(decoy_toks), 512)
        diff_pos = [i for i in range(min_len)
                    if real_toks[i] != decoy_toks[i] and real_toks[i] not in _SPECIAL]
        if not diff_pos:
            return 0.0, 0.0

        if len(real_toks) == len(decoy_toks):
            # share forward pass
            ctx = list(real_ids[:512])
            for p in diff_pos:
                if p < len(ctx):
                    ctx[p] = self.mask_id
            logits = self._forward_batch([ctx])
            lp = logits[0]

            real_ll, decoy_ll = 0.0, 0.0
            for p in diff_pos:
                if p >= lp.shape[0]:
                    continue
                pr = torch.softmax(lp[p, :], dim=-1)
                real_ll += math.log(max(pr[real_ids[p]].item(), 1e-20))
                decoy_ll += math.log(max(pr[decoy_ids[p]].item(), 1e-20))
            n = max(len(diff_pos), 1)
            return real_ll / n, decoy_ll / n
        else:
            # separate forward passes
            rm = list(real_ids[:512])
            r_pos = [p for p in diff_pos if p < len(rm)]
            for p in r_pos:
                rm[p] = self.mask_id

            dm = list(decoy_ids[:512])
            d_pos = [p for p in diff_pos if p < len(dm)]
            for p in d_pos:
                dm[p] = self.mask_id

            logits = self._forward_batch([rm, dm])
            lr, ld = logits[0], logits[1]

            real_ll = sum(math.log(max(torch.softmax(lr[p, :], -1)[real_ids[p]].item(), 1e-20))
                          for p in r_pos if p < lr.shape[0])
            decoy_ll = sum(math.log(max(torch.softmax(ld[p, :], -1)[decoy_ids[p]].item(), 1e-20))
                           for p in d_pos if p < ld.shape[0])
            return real_ll / max(len(r_pos), 1), decoy_ll / max(len(d_pos), 1)


# ── metal row helper ──────────────────────────────────────────────

def _metal_row(m):
    if m in {"Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn"}: return "3d"
    if m in {"Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd"}: return "4d"
    if m in {"La","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg"}: return "5d"
    return "other"


# ── main ──────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ranker_frozen",
                        default="")
    parser.add_argument("--ranker_finetuned",
                        default="")
    parser.add_argument("--mlm_checkpoint",
                        default="/data/CoordRep/CoordSMILES/libcoordrep/checkpoints/pretrain_v3/best_model.pt")
    parser.add_argument("--data",
                        default="/data/CoordRep/CoordSMILES/libcoordrep/figures/fig5/data/fig5_tasks.jsonl")
    parser.add_argument("--split_ids", default="",
                        help="Path to split_ids.json from training")
    parser.add_argument("--out",
                        default="/data/CoordRep/coordrep-release/libcoordrep/revision_results/coordrep_ranker")
    parser.add_argument("--decoys_per_real", type=int, default=20)
    parser.add_argument("--device", default="cuda:0")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    os.makedirs(args.out, exist_ok=True)

    print(f"{'='*60}")
    print("CoordRep-Ranker Evaluation")
    print(f"{'='*60}")

    # ── 1. Load data & split ──────────────────────────────────
    print("\nLoading data …")
    complexes = _extract_unique_complexes(args.data)

    if args.split_ids:
        with open(args.split_ids) as f:
            split = json.load(f)
        id_to_cx = {c["id"]: c for c in complexes}
        train_cx = [id_to_cx[i] for i in split["train"] if i in id_to_cx]
        val_cx = [id_to_cx[i] for i in split["val"] if i in id_to_cx]
        test_cx = [id_to_cx[i] for i in split["test"] if i in id_to_cx]
    else:
        train_cx, val_cx, test_cx = split_complexes(complexes, seed=args.seed)

    print(f"  train={len(train_cx)}, val={len(val_cx)}, test={len(test_cx)}")

    # ── 2. Build hard decoy pool (train only) ─────────────────
    print("Building hard decoy pool (train split) …")
    t0 = time.time()
    train_toks = [c["tokens"] for c in train_cx]
    pool = HardDecoyPool(train_toks, seed=args.seed)
    print(f"  ({time.time()-t0:.1f}s)")

    # ── 3. Load tokenizer ─────────────────────────────────────
    from pathlib import Path
    tok_path = Path(args.mlm_checkpoint).parent / "tokenizer.json"
    with open(tok_path) as f:
        tok_data = json.load(f)
    from brain.tokenizer import CoordRepTokenizer, TokenizerConfig
    tok_config = TokenizerConfig(**tok_data.get("config", {}))
    tokenizer = CoordRepTokenizer(tok_config)
    tokenizer.token2id = tok_data["token2id"]
    tokenizer.id2token = {v: k for k, v in tokenizer.token2id.items()}

    pad_id = tokenizer.token2id.get("[PAD]", 0)

    # ── 4. Generate test decoys ───────────────────────────────
    print("Generating test decoys …")
    t0 = time.time()
    test_data = []
    for c in test_cx:
        toks = c["tokens"]
        decoys = generate_all_hard_decoys(toks, pool, n_per_type=max(1, args.decoys_per_real // 5))
        if not decoys:
            continue
        import random as _rand
        rng = _rand.Random(args.seed)
        if len(decoys) > args.decoys_per_real:
            rng.shuffle(decoys)
            decoys = decoys[:args.decoys_per_real]
        test_data.append({"id": c["id"], "tokens": toks, "decoys": decoys})
    print(f"  {len(test_data)} test complexes with decoys ({time.time()-t0:.1f}s)")

    # ── 5. Leakage diagnostics ────────────────────────────────
    print("Checking leakage …")
    train_lig_set = set()
    for c in train_cx:
        blocks = _detect_blocks(c["tokens"])
        for lb in blocks.get("ligand_blocks", []):
            if lb.get("smiles_start") is not None and lb.get("smiles_end") is not None:
                lig = tuple(c["tokens"][lb["smiles_start"]:lb["smiles_end"]])
                train_lig_set.add(lig)

    leak_count = 0
    for td in test_data:
        blocks = _detect_blocks(td["tokens"])
        for lb in blocks.get("ligand_blocks", []):
            if lb.get("smiles_start") is not None and lb.get("smiles_end") is not None:
                lig = tuple(td["tokens"][lb["smiles_start"]:lb["smiles_end"]])
                if lig in train_lig_set:
                    leak_count += 1
                    break
    print(f"  Ligand overlap: {leak_count}/{len(test_data)} test complexes "
          f"share ≥1 ligand with train ({leak_count/max(len(test_data),1):.1%})")

    # ── 6. Score all methods ──────────────────────────────────
    print("\nScoring all methods …")

    # Baselines
    freq_bl = FrequencyBaseline()
    freq_bl.fit(train_toks)
    baselines_simple = {
        "random": RandomBaseline(seed=args.seed),
        "frequency": freq_bl,
        "rule_validator": RuleValidatorBaseline(),
    }

    # MLM scorer
    mlm_model, _ = load_model_and_tokenizer(args.mlm_checkpoint, device=args.device)
    mlm_scorer = MLMFieldScorer(mlm_model, tokenizer, args.device)

    # Ranker models
    ranker_models = {}
    if args.ranker_frozen and os.path.exists(args.ranker_frozen):
        r, _ = load_ranker(args.ranker_frozen, str(tok_path), device=args.device)
        ranker_models["ranker_frozen"] = r
    if args.ranker_finetuned and os.path.exists(args.ranker_finetuned):
        r, _ = load_ranker(args.ranker_finetuned, str(tok_path), device=args.device)
        ranker_models["ranker_finetuned"] = r

    def _pad(ids, max_len=512):
        ids = ids[:max_len]
        att = [1] * len(ids) + [0] * (max_len - len(ids))
        ids = ids + [pad_id] * (max_len - len(ids))
        return ids, att

    # ── Score every (complex, decoy) pair ─────────────────────
    all_records = []  # per-complex records
    dtype_scores = defaultdict(lambda: defaultdict(lambda: {"real": [], "decoy": []}))
    cn_scores = defaultdict(lambda: defaultdict(lambda: {"top1": 0, "n": 0, "mrr": 0.0}))

    for ti, td in enumerate(test_data):
        if ti % 200 == 0:
            print(f"  [{ti}/{len(test_data)}]", flush=True)

        toks = td["tokens"]
        cid = td["id"]
        decoys = td["decoys"]

        blocks = _detect_blocks(toks)
        metal, cn = "?", 0
        if blocks["metal_range"] is not None:
            p = _parse_metal_token(toks[blocks["metal_range"][0]])
            if p:
                metal, cn = p

        rec = {"complex_id": cid, "metal": metal, "CN": cn,
               "metal_row": _metal_row(metal), "n_decoys": len(decoys)}

        # Simple baselines
        import random as _rand
        _rng = _rand.Random(ti)
        for bl_name, bl in baselines_simple.items():
            bl_real = bl.score(toks)
            bl_scores = [(bl_real + _rng.uniform(-1e-9, 1e-9), "real")]
            for dtoks, _ in decoys:
                bl_scores.append((bl.score(dtoks) + _rng.uniform(-1e-9, 1e-9), "d"))
            bl_scores.sort(key=lambda x: -x[0])
            rank = next(i+1 for i,(_, r) in enumerate(bl_scores) if r == "real")
            rec[f"bl_{bl_name}_top1"] = int(rank == 1)
            rec[f"bl_{bl_name}_mrr"] = round(1.0 / rank, 4)

        # MLM field scorer
        mlm_wins = 0
        mlm_real_scores = []
        mlm_decoy_scores = []
        for dtoks, reason in decoys:
            rs, ds = mlm_scorer.score_pair(toks, dtoks)
            mlm_real_scores.append(rs)
            mlm_decoy_scores.append(ds)
            if rs > ds:
                mlm_wins += 1
            dtype = reason.split(":")[0]
            dtype_scores[dtype]["mlm_pll"]["real"].append(rs)
            dtype_scores[dtype]["mlm_pll"]["decoy"].append(ds)

        # For MLM PLL ranking: use mean real_score as the "real" score
        # (since each pair masks different tokens, real_score varies per decoy)
        mean_real = np.mean(mlm_real_scores) if mlm_real_scores else 0
        mlm_rank = 1 + sum(1 for d in mlm_decoy_scores if d >= mean_real)
        rec["mlm_pll_top1"] = int(mlm_rank == 1)
        rec["mlm_pll_win_rate"] = round(mlm_wins / max(len(decoys), 1), 4)
        rec["mlm_pll_mrr"] = round(1.0 / mlm_rank, 4) if mlm_real_scores else 0

        # Ranker models
        real_ids = _tokens_to_ids(toks, tokenizer)
        r_padded, r_att = _pad(real_ids)

        for rname, ranker in ranker_models.items():
            ranker.eval()
            with torch.no_grad():
                all_ids = [r_padded]
                all_att = [r_att]
                for dtoks, _ in decoys:
                    d_ids = _tokens_to_ids(dtoks, tokenizer)
                    dp, da = _pad(d_ids)
                    all_ids.append(dp)
                    all_att.append(da)

                ids_t = torch.tensor(all_ids, device=args.device).unsqueeze(0)  # (1, 1+K, L)
                att_t = torch.tensor(all_att, device=args.device).unsqueeze(0)
                scores = ranker.score_list(ids_t, att_t).squeeze(0)  # (1+K,)

                real_s = scores[0].item()
                decoy_s = scores[1:].cpu().numpy()

            rank = 1 + sum(1 for d in decoy_s if d >= real_s)
            wins = sum(1 for d in decoy_s if d < real_s)
            rec[f"{rname}_top1"] = int(rank == 1)
            rec[f"{rname}_mrr"] = round(1.0 / rank, 4)
            rec[f"{rname}_win_rate"] = round(wins / max(len(decoy_s), 1), 4)
            rec[f"{rname}_margin"] = round(real_s - float(max(decoy_s)) if len(decoy_s) > 0 else 0, 4)

            # per-decoy-type
            for di, (dtoks, reason) in enumerate(decoys):
                dtype = reason.split(":")[0]
                dtype_scores[dtype][rname]["real"].append(real_s)
                dtype_scores[dtype][rname]["decoy"].append(float(decoy_s[di]))

        # CN breakdown
        for method in list(ranker_models.keys()) + ["mlm_pll"] + [f"bl_{b}" for b in baselines_simple]:
            cn_scores[cn][method]["n"] += 1
            cn_scores[cn][method]["top1"] += rec.get(f"{method}_top1", 0)
            cn_scores[cn][method]["mrr"] += rec.get(f"{method}_mrr", 0)

        all_records.append(rec)

    # ── 7. Aggregate & write outputs ──────────────────────────
    print(f"\nAggregating {len(all_records)} test complexes …")
    n = len(all_records)

    # Summary
    summary = {"n_test": n, "n_decoys_per": args.decoys_per_real}

    all_methods = (
        [f"bl_{b}" for b in baselines_simple]
        + ["mlm_pll"]
        + list(ranker_models.keys())
    )

    for method in all_methods:
        t1_key = f"{method}_top1"
        mrr_key = f"{method}_mrr"
        wr_key = f"{method}_win_rate"
        t1_vals = [r.get(t1_key, 0) for r in all_records]
        mrr_vals = [r.get(mrr_key, 0) for r in all_records]
        wr_vals = [r.get(wr_key, None) for r in all_records]
        summary[f"{method}_top1"] = round(np.mean(t1_vals), 4)
        summary[f"{method}_mrr"] = round(np.mean(mrr_vals), 4)
        wr_valid = [w for w in wr_vals if w is not None]
        if wr_valid:
            summary[f"{method}_win_rate"] = round(np.mean(wr_valid), 4)

    # AUROC per method per decoy type
    by_decoy_type = []
    for dtype in sorted(dtype_scores.keys()):
        row = {"decoy_type": dtype}
        for method in all_methods:
            if method in dtype_scores[dtype]:
                d = dtype_scores[dtype][method]
                labels = [1] * len(d["real"]) + [0] * len(d["decoy"])
                scores = d["real"] + d["decoy"]
                auc = _auroc(labels, scores)
                row[f"{method}_auroc"] = round(auc, 4)
                summary[f"{dtype}_{method}_auroc"] = round(auc, 4)
        by_decoy_type.append(row)

    # Overall AUROC
    for method in all_methods:
        all_r, all_d = [], []
        for dtype in dtype_scores:
            if method in dtype_scores[dtype]:
                all_r.extend(dtype_scores[dtype][method]["real"])
                all_d.extend(dtype_scores[dtype][method]["decoy"])
        if all_r and all_d:
            labels = [1]*len(all_r) + [0]*len(all_d)
            scores = all_r + all_d
            summary[f"{method}_auroc"] = round(_auroc(labels, scores), 4)

    # ── Write outputs ─────────────────────────────────────────
    print("\nWriting outputs …")

    # summary.json
    p = os.path.join(args.out, "summary.json")
    with open(p, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"  {p}")

    # ranker_summary.csv
    p = os.path.join(args.out, "ranker_summary.csv")
    with open(p, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["metric", "value"])
        for k, v in sorted(summary.items()):
            w.writerow([k, v])
    print(f"  {p}")

    # ranker_by_decoy_type.csv
    if by_decoy_type:
        p = os.path.join(args.out, "ranker_by_decoy_type.csv")
        cols = sorted(by_decoy_type[0].keys())
        with open(p, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=cols)
            w.writeheader()
            for row in by_decoy_type:
                w.writerow(row)
        print(f"  {p}")

    # ranker_by_cn.csv
    p = os.path.join(args.out, "ranker_by_cn.csv")
    with open(p, "w", newline="") as f:
        w = csv.writer(f)
        header = ["CN"] + [f"{m}_top1" for m in all_methods] + [f"{m}_mrr" for m in all_methods] + ["n"]
        w.writerow(header)
        for cn in sorted(cn_scores.keys()):
            row = [cn]
            for m in all_methods:
                d = cn_scores[cn][m]
                row.append(round(d["top1"] / max(d["n"], 1), 4))
            for m in all_methods:
                d = cn_scores[cn][m]
                row.append(round(d["mrr"] / max(d["n"], 1), 4))
            row.append(cn_scores[cn][all_methods[0]]["n"])
            w.writerow(row)
    print(f"  {p}")

    # ranker_error_cases.csv — cases where best ranker gets it wrong
    best_ranker = None
    for rm in ["ranker_finetuned", "ranker_frozen"]:
        if rm in ranker_models:
            best_ranker = rm
            break

    if best_ranker:
        p = os.path.join(args.out, "ranker_error_cases.csv")
        errors = [r for r in all_records if r.get(f"{best_ranker}_top1", 0) == 0]
        errors.sort(key=lambda r: r.get(f"{best_ranker}_margin", 0))
        with open(p, "w", newline="") as f:
            cols = ["complex_id", "metal", "CN", "metal_row", "n_decoys",
                    f"{best_ranker}_mrr", f"{best_ranker}_margin",
                    "mlm_pll_top1", "bl_rule_validator_top1"]
            w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
            w.writeheader()
            for r in errors[:200]:
                w.writerow(r)
        print(f"  {p} ({len(errors)} errors)")

    # ranker_by_boundary.csv — boundary metals only
    p = os.path.join(args.out, "ranker_by_boundary.csv")
    boundary = [r for r in all_records if r.get("metal_row") in ("3d", "4d", "5d")]
    with open(p, "w", newline="") as f:
        cols = ["metal_row"] + [f"{m}_top1" for m in all_methods] + ["n"]
        w = csv.writer(f)
        w.writerow(cols)
        for mrow in ["3d", "4d", "5d"]:
            sub = [r for r in boundary if r["metal_row"] == mrow]
            if not sub:
                continue
            row = [mrow]
            for m in all_methods:
                row.append(round(np.mean([r.get(f"{m}_top1", 0) for r in sub]), 4))
            row.append(len(sub))
            w.writerow(row)
    print(f"  {p}")

    # ── Print key numbers ─────────────────────────────────────
    print(f"\n{'='*60}")
    print("RANKER EVALUATION — KEY NUMBERS")
    print(f"{'='*60}")
    print(f"  Test complexes: {n}")
    print()

    for method in all_methods:
        t1 = summary.get(f"{method}_top1", 0)
        mrr = summary.get(f"{method}_mrr", 0)
        wr = summary.get(f"{method}_win_rate", "N/A")
        auc = summary.get(f"{method}_auroc", "N/A")
        wr_s = f"{wr:.1%}" if isinstance(wr, float) else wr
        auc_s = f"{auc:.3f}" if isinstance(auc, float) else auc
        print(f"  {method:<25s} Top-1={t1:.1%}  MRR={mrr:.3f}  "
              f"WinRate={wr_s}  AUROC={auc_s}")

    print()
    print("  Per decoy-type AUROC:")
    for row in by_decoy_type:
        line = f"    {row['decoy_type']:<22s}"
        for method in all_methods:
            k = f"{method}_auroc"
            if k in row:
                line += f"  {method}={row[k]:.3f}"
        print(line)

    print()


if __name__ == "__main__":
    main()
