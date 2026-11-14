#!/usr/bin/env python3
"""
csd_toolb_transfer_eval.py
==========================
Task 3: CSD Tool B (spellchecker/repair) transfer evaluation.
Task 4: CSD hard-negative ranker transfer evaluation.

Reads retained CSD CoordRep strings from Task 1 and:
  - Applies the same corruption suite as build_synth_corruption
  - Attempts repair with the existing MLM checkpoint (no retraining)
  - Reports parse-valid recovery, exact recovery, by corruption type
  - Generates hard negatives and scores with the trained ranker

Usage:
    python scripts/csd_toolb_transfer_eval.py [--retained PATH]
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import random
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from coordrep_tools.validate import is_valid_coordrep
from coordrep_tools.tool_b_repair import CoordRepRepairer


# ── Corruption suite (same as build_synth_corruption) ─────────

def _apply_missing_bracket(chars: List[str]) -> Tuple[List[str], bool]:
    positions = [i for i, c in enumerate(chars) if c in '()[]{}']
    if not positions:
        return chars[:], False
    pos = random.choice(positions)
    out = chars[:]
    out.pop(pos)
    return out, True


def _apply_missing_charge(chars: List[str]) -> Tuple[List[str], bool]:
    positions = [i for i, c in enumerate(chars)
                 if c in '+-' and i + 1 < len(chars) and chars[i + 1].isdigit()]
    if not positions:
        return chars[:], False
    pos = random.choice(positions)
    out = chars[:]
    out.pop(pos)
    return out, True


def _apply_delimiter_mismatch(chars: List[str]) -> Tuple[List[str], bool]:
    dmap = {';': ',', ',': ';', '|': ';'}
    positions = [i for i, c in enumerate(chars) if c in dmap]
    if not positions:
        return chars[:], False
    pos = random.choice(positions)
    out = chars[:]
    out[pos] = dmap[out[pos]]
    return out, True


def _apply_truncation(chars: List[str]) -> Tuple[List[str], bool]:
    if len(chars) <= 15:
        return chars[:], False
    n_remove = random.randint(1, max(1, len(chars) // 5))
    return chars[:-n_remove], True


def _apply_token_swap(chars: List[str]) -> Tuple[List[str], bool]:
    positions = [i for i in range(len(chars) - 1)
                 if chars[i] in '()[]{}|;,']
    if not positions:
        return chars[:], False
    pos = random.choice(positions)
    out = chars[:]
    out[pos], out[pos + 1] = out[pos + 1], out[pos]
    return out, True


CORRUPTION_FNS = {
    'missing_bracket': _apply_missing_bracket,
    'missing_charge': _apply_missing_charge,
    'delimiter_mismatch': _apply_delimiter_mismatch,
    'truncation': _apply_truncation,
    'token_swap': _apply_token_swap,
}

CORRUPTION_PROBS = {
    'missing_bracket': 0.25,
    'missing_charge': 0.20,
    'delimiter_mismatch': 0.25,
    'truncation': 0.20,
    'token_swap': 0.10,
}


def generate_corrupted(coordrep_str: str, rng: random.Random) -> Tuple[str, str]:
    """Apply one random corruption. Returns (corrupted_str, corruption_type)."""
    types = list(CORRUPTION_PROBS.keys())
    weights = list(CORRUPTION_PROBS.values())
    ctype = rng.choices(types, weights=weights)[0]
    chars = list(coordrep_str)
    fn = CORRUPTION_FNS[ctype]
    corrupted, applied = fn(chars)
    if not applied:
        return coordrep_str, "none"
    return ''.join(corrupted), ctype


# ── Task 3: Tool B transfer ──────────────────────────────────

def run_tool_b_transfer(retained: List[dict], args) -> dict:
    print(f"\n{'='*60}")
    print("Task 3: Tool B Transfer Evaluation")
    print(f"{'='*60}")

    rng = random.Random(args.seed)

    # Sample entries for repair eval
    sample_size = min(args.n_repair, len(retained))
    sampled = rng.sample(retained, sample_size)
    print(f"  Sampled {sample_size} entries for corruption+repair")

    # Generate corrupted samples
    corrupted_samples = []
    for entry in sampled:
        cstr = entry["coordrep"]
        corrupted, ctype = generate_corrupted(cstr, rng)
        if ctype == "none":
            continue
        corrupted_samples.append({
            "id": entry["refcode"],
            "coordrep_clean": cstr,
            "coordrep_corrupted": corrupted,
            "corruption_type": ctype,
        })

    print(f"  Generated {len(corrupted_samples)} corrupted samples")
    type_counts = Counter(s["corruption_type"] for s in corrupted_samples)
    for t, c in type_counts.most_common():
        print(f"    {t}: {c}")

    # Repair
    print(f"  Loading repairer (checkpoint: {args.checkpoint}) …")
    repairer = CoordRepRepairer(args.checkpoint, args.tokenizer, device=args.device)

    results = {
        "total": 0,
        "valid_before": 0,
        "valid_after": 0,
        "exact_match": 0,
        "by_type": defaultdict(lambda: {
            "total": 0, "valid_before": 0, "valid_after": 0, "exact_match": 0
        }),
        "failures": [],
    }

    t0 = time.time()
    for i, sample in enumerate(corrupted_samples):
        if i % 100 == 0 and i > 0:
            print(f"    [{i}/{len(corrupted_samples)}] "
                  f"valid_after={results['valid_after']}/{results['total']}", flush=True)

        ctype = sample["corruption_type"]
        repair_result = repairer.repair(
            sample["coordrep_corrupted"], max_iters=5, k=5)

        results["total"] += 1
        results["by_type"][ctype]["total"] += 1

        if repair_result["was_valid_before"]:
            results["valid_before"] += 1
            results["by_type"][ctype]["valid_before"] += 1

        if repair_result["is_valid_after"]:
            results["valid_after"] += 1
            results["by_type"][ctype]["valid_after"] += 1

        if repair_result["repaired"] == sample["coordrep_clean"]:
            results["exact_match"] += 1
            results["by_type"][ctype]["exact_match"] += 1
        elif not repair_result["is_valid_after"]:
            results["failures"].append({
                "refcode": sample["id"],
                "corruption_type": ctype,
                "reason": repair_result.get("failure_reason", "unknown"),
            })

    elapsed = time.time() - t0
    n = results["total"]
    print(f"  Done ({elapsed:.0f}s): valid_after={results['valid_after']}/{n} "
          f"({results['valid_after']/max(n,1):.1%}), "
          f"exact={results['exact_match']}/{n} ({results['exact_match']/max(n,1):.1%})")

    # ── Write outputs ─────────────────────────────────────
    # Summary CSV
    p = os.path.join(args.out, "csd_toolb_transfer_summary.csv")
    with open(p, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["metric", "value"])
        w.writerow(["total", n])
        w.writerow(["valid_before", results["valid_before"]])
        w.writerow(["valid_after", results["valid_after"]])
        w.writerow(["exact_match", results["exact_match"]])
        w.writerow(["valid_rate_before", round(results["valid_before"] / max(n, 1), 4)])
        w.writerow(["valid_rate_after", round(results["valid_after"] / max(n, 1), 4)])
        w.writerow(["exact_match_rate", round(results["exact_match"] / max(n, 1), 4)])
        w.writerow(["uplift", round(
            (results["valid_after"] - results["valid_before"]) / max(n, 1), 4)])
    print(f"  {p}")

    # By corruption type
    p = os.path.join(args.out, "csd_toolb_by_corruption_type.csv")
    with open(p, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["corruption_type", "total", "valid_before", "valid_after",
                     "exact_match", "valid_rate_after", "exact_rate"])
        for ctype in sorted(results["by_type"].keys()):
            d = results["by_type"][ctype]
            t = d["total"]
            w.writerow([
                ctype, t, d["valid_before"], d["valid_after"], d["exact_match"],
                round(d["valid_after"] / max(t, 1), 4),
                round(d["exact_match"] / max(t, 1), 4),
            ])
    print(f"  {p}")

    # Failure examples
    p = os.path.join(args.out, "csd_toolb_failure_examples.csv")
    with open(p, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["refcode", "corruption_type", "reason"])
        w.writeheader()
        for r in results["failures"][:200]:
            w.writerow(r)
    print(f"  {p} ({len(results['failures'])} failures)")

    return results


# ── Task 4: Ranker transfer ──────────────────────────────────

def _auroc(labels, scores):
    if len(set(labels)) < 2:
        return 0.5
    pairs = list(zip(scores, labels))
    pairs.sort(key=lambda x: -x[0])
    tp, fp, prev_tp, prev_fp = 0, 0, 0, 0
    auc, prev_score = 0.0, None
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


def run_ranker_transfer(retained: List[dict], args) -> dict:
    print(f"\n{'='*60}")
    print("Task 4: Ranker Transfer Evaluation")
    print(f"{'='*60}")

    if not args.ranker_checkpoint or not os.path.exists(args.ranker_checkpoint):
        print(f"  SKIP: ranker checkpoint not found ({args.ranker_checkpoint})")
        return {}

    rng = random.Random(args.seed)

    # Import ranker and hard decoy tools
    from coordrep_tools.coordrep_ranker import load_ranker
    from coordrep_tools.ranker_dataset import _tokens_to_ids
    from coordrep_tools.hard_decoys import HardDecoyPool, generate_all_hard_decoys
    from brain.tokenizer import CoordRepTokenizer, TokenizerConfig
    import torch

    # Load tokenizer
    tok_path = args.tokenizer
    with open(tok_path) as f:
        tok_data = json.load(f)
    tok_config = TokenizerConfig(**tok_data.get("config", {}))
    tokenizer = CoordRepTokenizer(tok_config)
    tokenizer.token2id = tok_data["token2id"]
    tokenizer.id2token = {v: k for k, v in tokenizer.token2id.items()}

    pad_id = tokenizer.token2id.get("[PAD]", 0)

    def _pad(ids, max_len=512):
        ids = ids[:max_len]
        att = [1] * len(ids) + [0] * (max_len - len(ids))
        ids = ids + [pad_id] * (max_len - len(ids))
        return ids, att

    # Load ranker
    print(f"  Loading ranker from {args.ranker_checkpoint} …")
    ranker, _ = load_ranker(args.ranker_checkpoint, tok_path, device=args.device)
    ranker.eval()

    # Build hard decoy pool from CSD retained entries
    # We need token lists — tokenize the CoordRep strings
    print("  Tokenizing CSD CoordRep strings …")
    csd_token_lists = []
    for e in retained:
        toks = tokenizer.tokenize(e["coordrep"])
        csd_token_lists.append(toks)

    print(f"  Building hard decoy pool from {len(csd_token_lists)} CSD entries …")
    pool = HardDecoyPool(csd_token_lists, seed=args.seed)

    # Sample test entries
    sample_size = min(args.n_ranker, len(retained))
    indices = rng.sample(range(len(retained)), sample_size)
    print(f"  Evaluating ranker on {sample_size} CSD complexes …")

    decoy_types_to_test = [
        "stereo_hard", "metal_hard", "boundary_hard",
        "ligand_hard", "co_ligand_hard",
    ]

    dtype_scores = defaultdict(lambda: {"real": [], "decoy": []})
    records = []
    n_decoys_per = 20

    t0 = time.time()
    for ii, idx in enumerate(indices):
        if ii % 200 == 0 and ii > 0:
            print(f"    [{ii}/{sample_size}]", flush=True)

        toks = csd_token_lists[idx]
        entry = retained[idx]

        decoys = generate_all_hard_decoys(
            toks, pool, n_per_type=max(1, n_decoys_per // 5),
            types=decoy_types_to_test)
        if not decoys:
            continue
        if len(decoys) > n_decoys_per:
            rng.shuffle(decoys)
            decoys = decoys[:n_decoys_per]

        # Score
        real_ids = _tokens_to_ids(toks, tokenizer)
        rp, ra = _pad(real_ids)
        all_ids = [rp]
        all_att = [ra]
        for dtoks, reason in decoys:
            d_ids = _tokens_to_ids(dtoks, tokenizer)
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

        rec = {
            "refcode": entry["refcode"],
            "metal": entry["metal"],
            "cn": entry["cn"],
            "top1": int(rank == 1),
            "mrr": round(1.0 / rank, 4),
            "win_rate": round(wins / max(len(decoy_s), 1), 4),
            "n_decoys": len(decoys),
        }
        records.append(rec)

        for di, (dtoks, reason) in enumerate(decoys):
            dtype = reason.split(":")[0]
            dtype_scores[dtype]["real"].append(real_s)
            if di < len(decoy_s):
                dtype_scores[dtype]["decoy"].append(float(decoy_s[di]))

    elapsed = time.time() - t0
    n = len(records)
    print(f"  Done ({elapsed:.0f}s): {n} complexes scored")

    if n == 0:
        print("  No complexes scored — skipping output")
        return {}

    # Aggregate
    mean_top1 = sum(r["top1"] for r in records) / n
    mean_mrr = sum(r["mrr"] for r in records) / n
    mean_wr = sum(r["win_rate"] for r in records) / n

    overall_r = []
    overall_d = []
    by_dtype = {}
    for dtype in sorted(dtype_scores.keys()):
        d = dtype_scores[dtype]
        labels = [1] * len(d["real"]) + [0] * len(d["decoy"])
        scores_list = d["real"] + d["decoy"]
        auc = _auroc(labels, scores_list)
        by_dtype[dtype] = round(auc, 4)
        overall_r.extend(d["real"])
        overall_d.extend(d["decoy"])

    overall_auc = _auroc(
        [1] * len(overall_r) + [0] * len(overall_d),
        overall_r + overall_d)

    # ── Write outputs ─────────────────────────────────────
    p = os.path.join(args.out, "csd_ranker_transfer_summary.csv")
    with open(p, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["metric", "value"])
        w.writerow(["n_test", n])
        w.writerow(["top1", round(mean_top1, 4)])
        w.writerow(["mrr", round(mean_mrr, 4)])
        w.writerow(["win_rate", round(mean_wr, 4)])
        w.writerow(["overall_auroc", round(overall_auc, 4)])
        for dtype, auc in by_dtype.items():
            w.writerow([f"{dtype}_auroc", auc])
    print(f"  {p}")

    p = os.path.join(args.out, "csd_ranker_by_decoy_type.csv")
    with open(p, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["decoy_type", "auroc", "n_pairs"])
        for dtype in sorted(dtype_scores.keys()):
            d = dtype_scores[dtype]
            w.writerow([dtype, by_dtype[dtype], len(d["real"])])
    print(f"  {p}")

    print(f"\n  Ranker transfer results:")
    print(f"    Top-1={mean_top1:.1%}  MRR={mean_mrr:.3f}  "
          f"WinRate={mean_wr:.1%}  AUROC={overall_auc:.3f}")
    print(f"    Per decoy-type AUROC:")
    for dtype, auc in sorted(by_dtype.items()):
        print(f"      {dtype}: {auc:.3f}")

    return {"top1": mean_top1, "mrr": mean_mrr, "auroc": overall_auc,
            "by_dtype": by_dtype}


# ── main ──────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--retained",
                        default="revision_results/csd_external/csd_retained_entries.jsonl")
    parser.add_argument("--checkpoint",
                        default="/data/CoordRep/CoordSMILES/libcoordrep/checkpoints/pretrain_v3/best_model.pt")
    parser.add_argument("--tokenizer",
                        default="/data/CoordRep/CoordSMILES/libcoordrep/checkpoints/pretrain_v3/tokenizer.json")
    parser.add_argument("--ranker_checkpoint",
                        default="/data/CoordRep/coordrep-release/libcoordrep/checkpoints/coordrep_ranker/best_finetuned.pt")
    parser.add_argument("--out", default="revision_results/csd_external")
    parser.add_argument("--n_repair", type=int, default=2000,
                        help="Number of entries for repair eval")
    parser.add_argument("--n_ranker", type=int, default=1000,
                        help="Number of entries for ranker eval")
    parser.add_argument("--device", default="cuda:0")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--skip_repair", action="store_true")
    parser.add_argument("--skip_ranker", action="store_true")
    args = parser.parse_args()

    os.makedirs(args.out, exist_ok=True)

    # Load retained entries
    print("Loading retained CSD entries …")
    retained = []
    with open(args.retained) as f:
        for line in f:
            retained.append(json.loads(line))
    print(f"  {len(retained)} entries")

    if not args.skip_repair:
        run_tool_b_transfer(retained, args)

    if not args.skip_ranker:
        run_ranker_transfer(retained, args)

    print("\nAll CSD transfer evaluations complete.")


if __name__ == "__main__":
    main()
