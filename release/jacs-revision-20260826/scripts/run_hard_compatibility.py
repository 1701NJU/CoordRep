#!/usr/bin/env python3
"""
run_hard_compatibility.py
=========================
Hard-negative Counterfactual Coordination Compatibility Benchmark.

All decoys pass the CoordRep validator.  The rule-validator baseline
should drop to near-random, proving these are genuine hard negatives.

Scoring modes:
  1. Batched field-targeted (same-mask grouping)
  2. Conditional field scoring: P(field | context)
  3. PMI: log P(field|context) - log P(field)

Output → revision_results/coordrep_score_hard/
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

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from coordrep_tools.infer import load_model_and_tokenizer
from coordrep_tools.ablation_masking import _detect_blocks
from coordrep_tools.hard_decoys import (
    HardDecoyPool, generate_all_hard_decoys, HARD_DECOY_TYPES,
    _parse_metal_token,
)
from coordrep_tools.coordrep_score import (
    _tokens_to_ids, _SPECIAL, _classify_token_fields,
    score_field, compute_pmi_score,
)
from coordrep_tools.validate import is_valid_coordrep

import torch


# ── baselines ────────────────────────────────────────────────────

class RandomBaseline:
    def __init__(self, seed=42):
        self.rng = __import__("random").Random(seed)
    def score(self, tokens):
        return self.rng.random()


class FrequencyBaseline:
    def __init__(self):
        self.metal_cn_counts = Counter()
        self.total = 0
    def fit(self, all_token_lists):
        for toks in all_token_lists:
            blocks = _detect_blocks(toks)
            if blocks["metal_range"] is not None:
                parsed = _parse_metal_token(toks[blocks["metal_range"][0]])
                if parsed:
                    self.metal_cn_counts[(parsed[0], parsed[1])] += 1
                    self.total += 1
    def score(self, tokens):
        blocks = _detect_blocks(tokens)
        if blocks["metal_range"] is None:
            return -20.0
        parsed = _parse_metal_token(tokens[blocks["metal_range"][0]])
        if not parsed:
            return -20.0
        cnt = self.metal_cn_counts.get((parsed[0], parsed[1]), 0)
        return math.log(max(cnt, 1) / max(self.total, 1))


class RuleValidatorBaseline:
    """CoordRep rule validator — renamed from 'rule' baseline."""
    def score(self, tokens):
        s = "".join(tokens)
        return 1.0 if is_valid_coordrep(s, strict=True) else 0.0


# ── AUROC ────────────────────────────────────────────────────────

def _auroc(y_true, y_score):
    pairs = sorted(zip(y_score, y_true), reverse=True)
    tp = fp = tp_prev = fp_prev = 0
    auc = 0.0
    n_pos = sum(y_true)
    n_neg = len(y_true) - n_pos
    if n_pos == 0 or n_neg == 0:
        return 0.5
    prev_score = None
    for score, label in pairs:
        if prev_score is not None and score != prev_score:
            auc += (fp - fp_prev) * (tp + tp_prev) / 2.0
            tp_prev, fp_prev = tp, fp
        if label == 1:
            tp += 1
        else:
            fp += 1
        prev_score = score
    auc += (fp - fp_prev) * (tp + tp_prev) / 2.0
    return auc / (n_pos * n_neg)


# ── field mapping for decoy types ────────────────────────────────

DECOY_FIELD_MAP = {
    "metal_hard": "metal",
    "ligand_hard": "ligand_smiles",
    "stereo_hard": "stereo",
    "co_ligand_hard": "ligand_smiles",
    "boundary_hard": "metal",
}


# ── data loading ─────────────────────────────────────────────────

def _extract_unique_complexes(tasks_path):
    seen = {}
    with open(tasks_path) as f:
        for line in f:
            rec = json.loads(line)
            m = re.match(r'^(tmqm_\d+)', rec["id"])
            if not m:
                continue
            cid = m.group(1)
            if cid not in seen:
                seen[cid] = {"id": cid, "tokens": rec["tokens"]}
    return list(seen.values())


def _metal_row(m):
    if m in {"Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn"}: return "3d"
    if m in {"Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd"}: return "4d"
    if m in {"La","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg"}: return "5d"
    return "other"


# ── main ─────────────────────────────────────────────────────────

def _mega_batch_score(model_path, tokenizer_obj,
                      unique_seqs, unique_seq_positions,
                      lookups,
                      device, gpu_batch):
    """
    Score many (sequence, position, token_id) triples on one GPU.

    Parameters
    ----------
    unique_seqs : list[list[int]]
        Deduplicated masked sequences to forward.
    unique_seq_positions : list[list[int]]
        For each unique seq, the positions where we need softmax.
    lookups : list[tuple(int, list[int])]
        Each entry = (unique_seq_idx, [token_id_at_pos_0, token_id_at_pos_1, ...]).
        Returns one LL float per lookup.

    Returns
    -------
    list[float] – log-likelihood for each lookup entry.
    """
    model, _ = load_model_and_tokenizer(model_path, device=device)
    model.eval()
    pad_id = tokenizer_obj.token2id.get("[PAD]", 0)

    n_unique = len(unique_seqs)
    # Sort unique seqs by length for padding efficiency
    order = sorted(range(n_unique), key=lambda i: len(unique_seqs[i]))
    # For each unique seq, store per-position log-probs (sparse: only at requested positions)
    # pos_logprobs[uid] = {position: Tensor(vocab_size)} — but too large.
    # Instead, process batches and immediately extract all lookups.

    # Build reverse index: uid -> list of (lookup_idx, token_ids)
    uid_to_lookups = defaultdict(list)
    for li, (uid, tids) in enumerate(lookups):
        uid_to_lookups[uid].append((li, tids))

    results = [None] * len(lookups)
    n_batches = 0

    with torch.no_grad():
        for b_start in range(0, len(order), gpu_batch):
            b_end = min(b_start + gpu_batch, len(order))
            uids = order[b_start:b_end]
            seqs = [unique_seqs[u] for u in uids]
            max_l = max(len(s) for s in seqs)
            padded = [s + [pad_id] * (max_l - len(s)) for s in seqs]
            masks_a = [[1] * len(s) + [0] * (max_l - len(s)) for s in seqs]
            inp = torch.tensor(padded, device=device)
            att = torch.tensor(masks_a, device=device)
            logits, _ = model(inp, attention_mask=att)
            n_batches += 1

            for b_i, uid in enumerate(uids):
                positions = unique_seq_positions[uid]
                # Pre-compute softmax at all needed positions
                pos_probs = {}
                for p in positions:
                    if p < logits.shape[1]:
                        pos_probs[p] = torch.softmax(logits[b_i, p, :], dim=-1)

                for li, tids in uid_to_lookups[uid]:
                    ll = 0.0
                    for p, tid in zip(positions, tids):
                        if p in pos_probs:
                            ll += math.log(max(pos_probs[p][tid].item(), 1e-20))
                    results[li] = ll

    del model
    torch.cuda.empty_cache()
    return results, n_batches


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint",
                        default="inputs/checkpoints/pretrain_v3/best_model.pt")
    parser.add_argument("--test",
                        default="inputs/property_benchmarks/fig5_tasks.jsonl")
    parser.add_argument("--out",
                        default="outputs/coordrep_score_hard")
    parser.add_argument("--decoys-per-type", type=int, default=8)
    parser.add_argument("--max-complexes", type=int, default=0)
    parser.add_argument("--gpu-batch", type=int, default=128)
    parser.add_argument("--n-gpus", type=int, default=0,
                        help="0 = auto-detect")
    parser.add_argument("--cpu-workers", type=int, default=16)
    args = parser.parse_args()

    os.makedirs(args.out, exist_ok=True)

    if args.n_gpus <= 0:
        args.n_gpus = torch.cuda.device_count()
    print(f"Using {args.n_gpus} GPUs, {args.cpu_workers} CPU workers, "
          f"gpu_batch={args.gpu_batch}")

    # ── 1. Load data ─────────────────────────────────────────────
    print("Loading unique complexes …")
    complexes = _extract_unique_complexes(args.test)
    print(f"  {len(complexes)} unique complexes")
    all_token_lists = [c["tokens"] for c in complexes]
    all_ids = [c["id"] for c in complexes]

    # ── 2. Load tokenizer (model loaded per-GPU later) ───────────
    print("Loading tokenizer …")
    _, tokenizer = load_model_and_tokenizer(args.checkpoint, device="cpu")
    mask_id = tokenizer.token2id.get("[MASK]", 3)

    # ── 3. Build hard decoy pool ─────────────────────────────────
    print("Building hard decoy pool …")
    t0 = time.time()
    pool = HardDecoyPool(all_token_lists, seed=42)
    print(f"  {len(pool.metal_cns)} metals, "
          f"{sum(len(v) for v in pool.lig_by_len.values())} unique ligands, "
          f"({time.time()-t0:.1f}s)")

    # ── 4. Baselines ─────────────────────────────────────────────
    print("Fitting baselines …")
    freq_bl = FrequencyBaseline()
    freq_bl.fit(all_token_lists)
    baselines = {
        "random": RandomBaseline(seed=42),
        "frequency": freq_bl,
        "rule_validator": RuleValidatorBaseline(),
    }

    if args.max_complexes > 0:
        all_token_lists = all_token_lists[:args.max_complexes]
        all_ids = all_ids[:args.max_complexes]

    n_total = len(all_token_lists)

    # ══════════════════════════════════════════════════════════════
    # PHASE 1: CPU-parallel decoy generation
    # ══════════════════════════════════════════════════════════════
    from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
    import multiprocessing as mp

    print(f"\nPhase 1: Generating decoys ({n_total} complexes, {args.cpu_workers} workers) …")
    t1 = time.time()

    # Use ThreadPool (shares pool object) for decoy generation
    def _gen_one(args_tuple):
        idx, toks = args_tuple
        return idx, generate_all_hard_decoys(toks, pool, n_per_type=args.decoys_per_type)

    decoys_by_idx = {}
    with ThreadPoolExecutor(max_workers=args.cpu_workers) as ex:
        for result in ex.map(_gen_one, [(i, t) for i, t in enumerate(all_token_lists)]):
            idx, decoys = result
            if decoys:
                decoys_by_idx[idx] = decoys
            if len(decoys_by_idx) % 2000 == 0 and len(decoys_by_idx) > 0:
                print(f"    generated {len(decoys_by_idx)} …")

    total_decoys = sum(len(d) for d in decoys_by_idx.values())
    print(f"  Phase 1 done: {len(decoys_by_idx)} complexes, "
          f"{total_decoys} decoys ({time.time()-t1:.1f}s)")

    # ══════════════════════════════════════════════════════════════
    # PHASE 2: Build deduplicated sequence list + lookup table
    # ══════════════════════════════════════════════════════════════
    print(f"\nPhase 2: Building mega-batch …")
    t2 = time.time()

    unique_seqs = []       # deduplicated masked ID sequences
    unique_positions = []  # positions to softmax for each unique seq
    seq_to_uid = {}        # tuple(seq) -> uid for dedup
    lookups = []           # (uid, [token_ids]) — one per LL we need
    # scoring_plan[complex_idx] = list of {di, reason, real_lookup_idx, decoy_lookup_idx, n_pos_r, n_pos_d}
    scoring_plan = {}
    job_metadata = {}

    def _get_uid(seq_ids, positions):
        key = tuple(seq_ids)
        if key not in seq_to_uid:
            uid = len(unique_seqs)
            seq_to_uid[key] = uid
            unique_seqs.append(seq_ids)
            unique_positions.append(positions)
        return seq_to_uid[key]

    for idx in sorted(decoys_by_idx.keys()):
        toks = all_token_lists[idx]
        cid = all_ids[idx]
        decoys = decoys_by_idx[idx]
        real_ids = _tokens_to_ids(toks, tokenizer)
        job_metadata[idx] = (cid, toks, decoys)

        same_len_groups = defaultdict(list)
        diff_len_list = []

        for di, (dtoks, reason) in enumerate(decoys):
            if len(toks) == len(dtoks):
                diff_pos = tuple(i for i in range(min(len(toks), 512))
                                 if toks[i] != dtoks[i] and toks[i] not in _SPECIAL)
                if not diff_pos:
                    continue
                d_ids = _tokens_to_ids(dtoks, tokenizer)
                same_len_groups[diff_pos].append((di, d_ids, reason))
            else:
                diff_len_list.append((di, dtoks, reason))

        entries = []

        # Same-length: one sequence, multiple token lookups
        for mask_key, items in same_len_groups.items():
            ctx = list(real_ids[:512])
            for p in mask_key:
                if p < len(ctx):
                    ctx[p] = mask_id
            positions = list(mask_key)
            uid = _get_uid(ctx, positions)

            # Real lookup
            real_tids = [real_ids[p] if p < len(real_ids) else 1 for p in mask_key]
            real_li = len(lookups)
            lookups.append((uid, real_tids))

            for di, d_ids, reason in items:
                decoy_tids = [d_ids[p] if p < len(d_ids) else 1 for p in mask_key]
                decoy_li = len(lookups)
                lookups.append((uid, decoy_tids))
                entries.append({
                    "di": di, "reason": reason,
                    "real_li": real_li, "decoy_li": decoy_li,
                    "n_pos": len(mask_key),
                })

        # Diff-length: two sequences per decoy
        for di, dtoks, reason in diff_len_list:
            d_ids = _tokens_to_ids(dtoks, tokenizer)
            min_len = min(len(toks), len(dtoks), 512)
            diff_pos = [i for i in range(min_len)
                        if toks[i] != dtoks[i] and toks[i] not in _SPECIAL]
            if len(dtoks) > min_len:
                diff_pos.extend(range(min_len, min(len(dtoks), 512)))
            if len(toks) > min_len:
                diff_pos.extend(range(min_len, min(len(toks), 512)))
            if not diff_pos:
                continue

            r_pos = [p for p in diff_pos if p < min(len(real_ids), 512)]
            d_pos = [p for p in diff_pos if p < min(len(d_ids), 512)]

            rm = list(real_ids[:512])
            for p in r_pos: rm[p] = mask_id
            uid_r = _get_uid(rm, r_pos)
            real_li = len(lookups)
            lookups.append((uid_r, [real_ids[p] for p in r_pos]))

            dm = list(d_ids[:512])
            for p in d_pos: dm[p] = mask_id
            uid_d = _get_uid(dm, d_pos)
            decoy_li = len(lookups)
            lookups.append((uid_d, [d_ids[p] for p in d_pos]))

            entries.append({
                "di": di, "reason": reason,
                "real_li": real_li, "decoy_li": decoy_li,
                "n_pos": max(len(r_pos), len(d_pos)),
                "n_r": len(r_pos), "n_d": len(d_pos),
            })

        if entries:
            scoring_plan[idx] = entries

    print(f"  {len(unique_seqs)} unique seqs, {len(lookups)} lookups, "
          f"{len(scoring_plan)} complexes ({time.time()-t2:.1f}s)")

    # ══════════════════════════════════════════════════════════════
    # PHASE 3: Multi-GPU mega-batch scoring
    # ══════════════════════════════════════════════════════════════
    n_gpus = min(args.n_gpus, max(1, len(unique_seqs) // 50))
    n_gpus = max(1, min(n_gpus, torch.cuda.device_count()))
    print(f"\nPhase 3: Scoring {len(unique_seqs)} unique seqs / "
          f"{len(lookups)} lookups on {n_gpus} GPUs (batch={args.gpu_batch}) …")
    t3 = time.time()

    if n_gpus == 1:
        flat_results, n_bat = _mega_batch_score(
            args.checkpoint, tokenizer,
            unique_seqs, unique_positions, lookups,
            "cuda:0", args.gpu_batch)
    else:
        import threading
        # Shard unique seqs across GPUs; remap lookups to per-shard UIDs
        shard_size = (len(unique_seqs) + n_gpus - 1) // n_gpus
        flat_results = [None] * len(lookups)
        lock = threading.Lock()
        total_batches = [0]

        def _worker(rank):
            uid_lo = rank * shard_size
            uid_hi = min(uid_lo + shard_size, len(unique_seqs))
            if uid_lo >= uid_hi:
                return
            uid_set = set(range(uid_lo, uid_hi))
            # Collect lookups for this shard's UIDs
            local_lookups = []
            global_li_map = []
            for li, (uid, tids) in enumerate(lookups):
                if uid in uid_set:
                    local_lookups.append((uid - uid_lo, tids))
                    global_li_map.append(li)
            if not local_lookups:
                return
            local_seqs = unique_seqs[uid_lo:uid_hi]
            local_pos = unique_positions[uid_lo:uid_hi]
            local_results, nb = _mega_batch_score(
                args.checkpoint, tokenizer,
                local_seqs, local_pos, local_lookups,
                f"cuda:{rank}", args.gpu_batch)
            with lock:
                total_batches[0] += nb
                for i, v in enumerate(local_results):
                    flat_results[global_li_map[i]] = v

        threads = []
        for rank in range(n_gpus):
            t = threading.Thread(target=_worker, args=(rank,))
            t.start()
            threads.append(t)
        for t in threads:
            t.join()
        n_bat = total_batches[0]

    print(f"  Phase 3 done: {n_bat} batches ({time.time()-t3:.1f}s)")

    # ══════════════════════════════════════════════════════════════
    # PHASE 3b: Reconstruct per-complex decoy results
    # ══════════════════════════════════════════════════════════════
    result_dict = {}
    for idx, entries in scoring_plan.items():
        decoy_results = {}
        for e in entries:
            real_ll = flat_results[e["real_li"]] or 0.0
            decoy_ll = flat_results[e["decoy_li"]] or 0.0
            n_r = e.get("n_r", e["n_pos"])
            n_d = e.get("n_d", e["n_pos"])
            rfs = real_ll / max(n_r, 1)
            dfs = decoy_ll / max(n_d, 1)
            decoy_results[e["di"]] = (dfs, rfs, e["reason"])
        if decoy_results:
            result_dict[idx] = decoy_results

    # ══════════════════════════════════════════════════════════════
    # PHASE 4: Aggregate results
    # ══════════════════════════════════════════════════════════════
    print(f"\nPhase 4: Aggregating results …")
    import random as _rand

    per_complex = []
    all_labels = []
    all_scores_flat = []
    decoy_type_scores = defaultdict(lambda: {"real": [], "decoy": []})
    leakage_rows = []

    for idx in sorted(result_dict.keys()):
        decoy_results = result_dict[idx]
        if not decoy_results:
            continue
        cid, toks, decoys = job_metadata[idx]

        pairwise_deltas = []
        for di in sorted(decoy_results.keys()):
            dfs, rfs, reason = decoy_results[di]
            pairwise_deltas.append(rfs - dfs)

        n_wins = sum(1 for d in pairwise_deltas if d > 0)
        n_losses = sum(1 for d in pairwise_deltas if d <= 0)
        real_rank = 1 + n_losses
        win_rate = n_wins / len(pairwise_deltas)
        mean_delta = float(np.mean(pairwise_deltas))

        blocks = _detect_blocks(toks)
        metal, cn = "?", 0
        if blocks["metal_range"] is not None:
            p = _parse_metal_token(toks[blocks["metal_range"][0]])
            if p: metal, cn = p

        rec = {
            "complex_id": cid, "metal": metal, "CN": cn,
            "metal_row": _metal_row(metal),
            "n_ligs": len(blocks["ligand_blocks"]),
            "has_stereo": len(blocks.get("constraint_ranges", [])) > 0,
            "n_decoys": len(decoy_results),
            "n_wins": n_wins, "win_rate": round(win_rate, 4),
            "top1_real": int(n_wins == len(pairwise_deltas)),
            "real_rank": real_rank,
            "mrr": round(1.0 / real_rank, 4),
            "mean_delta": round(mean_delta, 4),
        }

        _rng = _rand.Random(idx)
        for bl_name, bl in baselines.items():
            bl_real = bl.score(toks)
            bl_decoy = [bl.score(dt) for dt, _ in decoys if dt is not None]
            bl_all = [(bl_real + _rng.uniform(-1e-9, 1e-9), "real")]
            for s in bl_decoy:
                bl_all.append((s + _rng.uniform(-1e-9, 1e-9), "d"))
            bl_all.sort(key=lambda x: -x[0])
            bl_rank = next(i+1 for i,(_, r) in enumerate(bl_all) if r == "real")
            rec[f"bl_{bl_name}_rank"] = bl_rank
            rec[f"bl_{bl_name}_top1"] = int(bl_rank == 1)
            rec[f"bl_{bl_name}_mrr"] = round(1.0 / bl_rank, 4)

        rule_bl = baselines["rule_validator"]
        rule_pass_count = sum(1 for dt, _ in decoys if rule_bl.score(dt) >= 1.0)
        rule_pass_rate = rule_pass_count / len(decoys) if decoys else 0
        rec["rule_pass_rate"] = round(rule_pass_rate, 4)
        if rule_pass_rate < 1.0:
            leakage_rows.append({
                "complex_id": cid, "n_decoys": len(decoys),
                "rule_pass_count": rule_pass_count,
                "rule_pass_rate": round(rule_pass_rate, 4),
            })

        per_complex.append(rec)

        for dfs, rfs, reason in decoy_results.values():
            all_labels.extend([1, 0])
            all_scores_flat.extend([rfs, dfs])
            dtype = reason.split(":")[0]
            decoy_type_scores[dtype]["real"].append(rfs)
            decoy_type_scores[dtype]["decoy"].append(dfs)

    elapsed = time.time() - t1
    print(f"  Total time: {elapsed:.1f}s")

    # ── 7. Aggregate metrics ─────────────────────────────────────
    n = len(per_complex)
    if n == 0:
        print("  No complexes scored!")
        return

    top1 = sum(r["top1_real"] for r in per_complex) / n
    mrr = sum(r["mrr"] for r in per_complex) / n
    mean_wr = sum(r["win_rate"] for r in per_complex) / n
    mean_d = sum(r["mean_delta"] for r in per_complex) / n
    auroc = _auroc(all_labels, all_scores_flat)

    summary = {
        "n_complexes": n,
        "n_decoys_total": sum(r["n_decoys"] for r in per_complex),
        "mlm_field_top1": round(top1, 4),
        "mlm_field_mrr": round(mrr, 4),
        "mlm_field_auroc": round(auroc, 4),
        "mlm_field_win_rate": round(mean_wr, 4),
        "mean_delta": round(mean_d, 4),
    }

    for bl_name in baselines:
        bt = sum(r.get(f"bl_{bl_name}_top1", 0) for r in per_complex) / n
        bm = sum(r.get(f"bl_{bl_name}_mrr", 0) for r in per_complex) / n
        summary[f"bl_{bl_name}_top1"] = round(bt, 4)
        summary[f"bl_{bl_name}_mrr"] = round(bm, 4)

    by_decoy_type = {}
    for dtype in sorted(decoy_type_scores.keys()):
        d = decoy_type_scores[dtype]
        labels = [1]*len(d["real"]) + [0]*len(d["decoy"])
        scores = d["real"] + d["decoy"]
        auc = _auroc(labels, scores)
        by_decoy_type[dtype] = {
            "n_real": len(d["real"]), "n_decoy": len(d["decoy"]),
            "detection_auroc": round(auc, 4),
            "mean_real": round(float(np.mean(d["real"])), 4) if d["real"] else 0,
            "mean_decoy": round(float(np.mean(d["decoy"])), 4) if d["decoy"] else 0,
        }
        summary[f"{dtype}_auroc"] = round(auc, 4)

    # Rule validator pass rate
    mean_rule_pass = np.mean([r["rule_pass_rate"] for r in per_complex])
    summary["mean_rule_pass_rate"] = round(float(mean_rule_pass), 4)

    # ── 8. Write outputs ─────────────────────────────────────────
    print("\nWriting outputs …")

    # 8a. summary.json
    path = os.path.join(args.out, "summary.json")
    with open(path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"  {path}")

    # 8b. hard_compatibility_summary.csv
    path = os.path.join(args.out, "hard_compatibility_summary.csv")
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["metric", "value"])
        for k, v in sorted(summary.items()):
            w.writerow([k, v])
    print(f"  {path}")

    # 8c. hard_by_decoy_type.csv
    path = os.path.join(args.out, "hard_by_decoy_type.csv")
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["decoy_type","n_real","n_decoy",
                                          "detection_auroc","mean_real","mean_decoy"])
        w.writeheader()
        for dtype, d in sorted(by_decoy_type.items()):
            w.writerow({"decoy_type": dtype, **d})
    print(f"  {path}")

    # 8d. hard_rule_validator_performance.csv
    path = os.path.join(args.out, "hard_rule_validator_performance.csv")
    with open(path, "w", newline="") as f:
        cols = ["complex_id","n_decoys","rule_pass_rate",
                "bl_rule_validator_top1","bl_rule_validator_mrr"]
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in per_complex[:200]:
            w.writerow(r)
    print(f"  {path}")

    # 8e. leakage_diagnostics.csv
    path = os.path.join(args.out, "leakage_diagnostics.csv")
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["complex_id","n_decoys",
                                          "rule_pass_count","rule_pass_rate"])
        w.writeheader()
        for r in leakage_rows:
            w.writerow(r)
    print(f"  {path}")

    # 8f. hard_mlm_fieldwise_score.csv
    path = os.path.join(args.out, "hard_mlm_fieldwise_score.csv")
    with open(path, "w", newline="") as f:
        cols = ["complex_id","metal","CN","n_decoys","n_wins","win_rate",
                "top1_real","mrr","mean_delta","rule_pass_rate",
                "bl_random_top1","bl_frequency_top1","bl_rule_validator_top1"]
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in per_complex:
            w.writerow(r)
    print(f"  {path}")

    # ── 9. Print key numbers ─────────────────────────────────────
    print(f"\n{'='*60}")
    print("HARD BENCHMARK — KEY NUMBERS")
    print(f"{'='*60}")
    print(f"  Complexes:                  {n}")
    print(f"  Total hard decoys:          {summary['n_decoys_total']}")
    print(f"  Mean rule pass rate:        {summary['mean_rule_pass_rate']:.1%}")
    print()
    print(f"  MLM field-targeted Top-1:   {summary['mlm_field_top1']:.1%}")
    print(f"  MLM field-targeted MRR:     {summary['mlm_field_mrr']:.3f}")
    print(f"  MLM field-targeted AUROC:   {summary['mlm_field_auroc']:.3f}")
    print(f"  MLM pairwise win rate:      {summary['mlm_field_win_rate']:.1%}")
    print()
    print("  Per decoy-type AUROC:")
    for dtype in HARD_DECOY_TYPES:
        auc = summary.get(f"{dtype}_auroc", 0)
        print(f"    {dtype:<25s} {auc:.3f}")
    print()
    print("  Baselines (Top-1):")
    mlm_t1 = summary["mlm_field_top1"]
    for bl_name in baselines:
        bt = summary.get(f"bl_{bl_name}_top1", 0)
        print(f"    {bl_name:<25s} {bt:.1%}  (MLM gain: {mlm_t1 - bt:+.1%})")
    print()


if __name__ == "__main__":
    main()
