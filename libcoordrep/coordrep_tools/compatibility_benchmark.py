"""
compatibility_benchmark.py
==========================
Counterfactual Coordination Compatibility Benchmark.

For each real CoordRep complex, generate N decoys, score all with
MLM pseudo-log-likelihood, and evaluate whether the real complex
ranks above its decoys.

Metrics
-------
- **Top-1 recovery**: fraction where real has the highest score.
- **MRR**: mean reciprocal rank of the real among decoys.
- **AUROC**: binary classification (real=1 vs decoy=0) from scores.
- **ΔScore**: mean score gap (real − best decoy).
- Per decoy-type detection AUC.

Baselines
---------
A. Random ranking
B. Frequency score: P(metal, CN) + P(donor_pattern | metal, CN)
C. LigandFreq-only score
D. CoordRep-ID rule score (validate + heuristic)
"""

from __future__ import annotations

import math
import random
import re
from collections import Counter, defaultdict
from typing import Dict, List, Optional, Tuple

import numpy as np

from .coordrep_score import (
    compute_pseudo_log_likelihood,
    compute_normalized_pll,
    compute_fieldwise_scores,
    _classify_token_fields,
    _tokens_to_ids,
    _SPECIAL,
)
from .counterfactual_decoys import (
    DecoyPool,
    generate_all_decoys,
    DECOY_TYPES,
    _parse_metal_token,
)
from .ablation_masking import _detect_blocks
from .validate import is_valid_coordrep

# ── helpers ──────────────────────────────────────────────────────

_3D = {"Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn"}
_4D = {"Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd"}
_5D = {"La", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg"}


def _metal_row(m):
    if m in _3D: return "3d"
    if m in _4D: return "4d"
    if m in _5D: return "5d"
    return "other"


# ── baselines ────────────────────────────────────────────────────

class RandomBaseline:
    """Assign random scores."""
    def __init__(self, seed=42):
        self.rng = random.Random(seed)

    def score(self, tokens):
        return self.rng.random()


class FrequencyBaseline:
    """Score = log P(metal,CN) from training distribution."""
    def __init__(self):
        self.metal_cn_counts = Counter()
        self.total = 0

    def fit(self, all_token_lists):
        for toks in all_token_lists:
            blocks = _detect_blocks(toks)
            if blocks["metal_range"] is not None:
                mt = toks[blocks["metal_range"][0]]
                parsed = _parse_metal_token(mt)
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
        if cnt == 0:
            return -20.0
        return math.log(cnt / max(self.total, 1))


class LigandFreqBaseline:
    """Score = sum of log P(ligand_smiles) over ligand blocks."""
    def __init__(self):
        self.lig_counts = Counter()
        self.total = 0

    def fit(self, all_token_lists):
        for toks in all_token_lists:
            blocks = _detect_blocks(toks)
            for lb in blocks["ligand_blocks"]:
                smi = "".join(toks[lb["smiles_start"]:lb["smiles_end"]])
                self.lig_counts[smi] += 1
                self.total += 1

    def score(self, tokens):
        blocks = _detect_blocks(tokens)
        s = 0.0
        for lb in blocks["ligand_blocks"]:
            smi = "".join(tokens[lb["smiles_start"]:lb["smiles_end"]])
            cnt = self.lig_counts.get(smi, 0)
            s += math.log(max(cnt, 1) / max(self.total, 1))
        return s


class RuleBaseline:
    """Score = simple heuristic: valid string gets 1.0, else 0.0."""
    def score(self, tokens):
        s = "".join(tokens)
        return 1.0 if is_valid_coordrep(s, strict=True) else 0.0


# ── core benchmark ───────────────────────────────────────────────

def _auroc(y_true, y_score):
    """Compute AUROC. y_true: binary labels, y_score: scores."""
    pairs = sorted(zip(y_score, y_true), reverse=True)
    tp = fp = 0
    tp_prev = fp_prev = 0
    auc = 0.0
    n_pos = sum(y_true)
    n_neg = len(y_true) - n_pos
    if n_pos == 0 or n_neg == 0:
        return 0.5
    prev_score = None
    for score, label in pairs:
        if prev_score is not None and score != prev_score:
            auc += (fp - fp_prev) * (tp + tp_prev) / 2.0
            tp_prev = tp
            fp_prev = fp
        if label == 1:
            tp += 1
        else:
            fp += 1
        prev_score = score
    auc += (fp - fp_prev) * (tp + tp_prev) / 2.0
    return auc / (n_pos * n_neg)


def run_benchmark(
    model,
    tokenizer,
    test_token_lists: List[List[str]],
    test_ids: List[str],
    pool: DecoyPool,
    device: str = "cuda",
    n_decoys_per_type: int = 10,
    pll_mode: str = "mask_one_token",
    baselines: Optional[Dict[str, object]] = None,
    max_complexes: int = 0,
    verbose: bool = True,
) -> Dict:
    """
    Run the full compatibility benchmark.

    Returns
    -------
    dict with keys: per_complex, summary, by_decoy_type, by_cn,
    by_metal_row, fieldwise_examples, anomalies.
    """
    if max_complexes > 0:
        test_token_lists = test_token_lists[:max_complexes]
        test_ids = test_ids[:max_complexes]

    if baselines is None:
        baselines = {}

    per_complex = []
    all_real_scores = []
    all_decoy_scores = []
    all_labels = []  # 1=real, 0=decoy
    all_scores_flat = []
    decoy_type_scores = defaultdict(lambda: {"real": [], "decoy": []})
    fieldwise_examples = []

    import torch, math
    mask_id = tokenizer.token2id.get("[MASK]", 3)
    unk_id = tokenizer.token2id.get("[UNK]", 1)
    pad_id = tokenizer.token2id.get("[PAD]", 0)
    GPU_BATCH = 64  # forward-pass batch size

    def _batched_forward(id_seqs, device):
        """Run model on a batch of id-sequences, return logits list."""
        if not id_seqs:
            return []
        max_l = max(len(s) for s in id_seqs)
        padded = [s + [pad_id] * (max_l - len(s)) for s in id_seqs]
        masks = [[1] * len(s) + [0] * (max_l - len(s)) for s in id_seqs]
        inp = torch.tensor(padded, device=device)
        att = torch.tensor(masks, device=device)
        logits, _ = model(inp, attention_mask=att)
        return logits

    n_total = len(test_token_lists)
    n_fwd = 0  # forward-pass counter

    for idx, (toks, cid) in enumerate(zip(test_token_lists, test_ids)):
        if verbose and idx % 200 == 0:
            print(f"  [{idx}/{n_total}] scoring {cid} … (fwd={n_fwd})")

        # Generate decoys
        decoys = generate_all_decoys(toks, pool, n_per_type=n_decoys_per_type)
        if not decoys:
            continue

        real_ids = _tokens_to_ids(toks, tokenizer)

        # ── Classify decoys by mask-pattern ──────────────────
        # Group same-length decoys by their diff-position tuple so
        # decoys sharing the same mask pattern need only ONE forward pass.
        same_len_groups = defaultdict(list)  # mask_key -> [(decoy_idx, decoy_ids)]
        diff_len_decoys = []  # (decoy_idx, dtoks, reason)
        decoy_id_cache = {}

        for di, (dtoks, reason) in enumerate(decoys):
            if len(toks) == len(dtoks):
                diff_pos = tuple(i for i in range(min(len(toks), 512))
                                 if toks[i] != dtoks[i] and toks[i] not in _SPECIAL)
                if not diff_pos:
                    continue
                d_ids = _tokens_to_ids(dtoks, tokenizer)
                decoy_id_cache[di] = d_ids
                same_len_groups[diff_pos].append((di, d_ids, reason))
            else:
                diff_len_decoys.append((di, dtoks, reason))

        # ── Score same-length groups (batched) ───────────────
        decoy_scores = [None] * len(decoys)
        decoy_reasons_out = [None] * len(decoys)
        pairwise_deltas = []
        pairwise_real_scores = []

        # Collect all unique masked sequences for batched forward
        batch_seqs = []
        batch_meta = []  # (mask_key, group_items)

        for mask_key, items in same_len_groups.items():
            context_ids = list(real_ids[:512])
            for p in mask_key:
                if p < len(context_ids):
                    context_ids[p] = mask_id
            batch_seqs.append(context_ids)
            batch_meta.append((mask_key, items))

        # Run batched forward passes
        with torch.no_grad():
            for b_start in range(0, len(batch_seqs), GPU_BATCH):
                b_end = min(b_start + GPU_BATCH, len(batch_seqs))
                chunk = batch_seqs[b_start:b_end]
                logits = _batched_forward(chunk, device)
                n_fwd += 1

                for b_i in range(len(chunk)):
                    mask_key, items = batch_meta[b_start + b_i]
                    log_probs = logits[b_i]  # (seq_len, vocab)
                    n_pos = max(len(mask_key), 1)

                    # Compute real field score once per mask pattern
                    real_ll = 0.0
                    for p in mask_key:
                        if p >= log_probs.shape[0]:
                            continue
                        probs = torch.softmax(log_probs[p, :], dim=-1)
                        r_id = real_ids[p] if p < len(real_ids) else unk_id
                        real_ll += math.log(max(probs[r_id].item(), 1e-20))
                    real_field_score = real_ll / n_pos

                    # Score each decoy in this group (just token lookups, no extra forward)
                    for di, d_ids, reason in items:
                        decoy_ll = 0.0
                        for p in mask_key:
                            if p >= log_probs.shape[0]:
                                continue
                            probs = torch.softmax(log_probs[p, :], dim=-1)
                            did = d_ids[p] if p < len(d_ids) else unk_id
                            decoy_ll += math.log(max(probs[did].item(), 1e-20))
                        decoy_field_score = decoy_ll / n_pos

                        decoy_scores[di] = decoy_field_score
                        decoy_reasons_out[di] = reason
                        pairwise_deltas.append(real_field_score - decoy_field_score)
                        pairwise_real_scores.append(real_field_score)

        # ── Score different-length decoys (batched, 2 passes each) ──
        if diff_len_decoys:
            real_batch = []
            decoy_batch = []
            dl_meta = []
            for di, dtoks, reason in diff_len_decoys:
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
                for p in r_pos:
                    rm[p] = mask_id
                dm = list(d_ids[:512])
                for p in d_pos:
                    dm[p] = mask_id

                real_batch.append(rm)
                decoy_batch.append(dm)
                dl_meta.append((di, d_ids, reason, r_pos, d_pos))

            with torch.no_grad():
                # Batch real + decoy sequences together
                all_seqs = real_batch + decoy_batch
                all_logits_list = []
                for b_start in range(0, len(all_seqs), GPU_BATCH):
                    b_end = min(b_start + GPU_BATCH, len(all_seqs))
                    logits = _batched_forward(all_seqs[b_start:b_end], device)
                    n_fwd += 1
                    for j in range(logits.shape[0]):
                        all_logits_list.append(logits[j])

                n_dl = len(dl_meta)
                for k, (di, d_ids, reason, r_pos, d_pos) in enumerate(dl_meta):
                    log_r = all_logits_list[k]
                    log_d = all_logits_list[n_dl + k]

                    real_ll = 0.0
                    for p in r_pos:
                        if p >= log_r.shape[0]: continue
                        probs = torch.softmax(log_r[p, :], dim=-1)
                        real_ll += math.log(max(probs[real_ids[p]].item(), 1e-20))

                    decoy_ll = 0.0
                    for p in d_pos:
                        if p >= log_d.shape[0]: continue
                        probs = torch.softmax(log_d[p, :], dim=-1)
                        decoy_ll += math.log(max(probs[d_ids[p]].item(), 1e-20))

                    rfs = real_ll / max(len(r_pos), 1)
                    dfs = decoy_ll / max(len(d_pos), 1)
                    decoy_scores[di] = dfs
                    decoy_reasons_out[di] = reason
                    pairwise_deltas.append(rfs - dfs)
                    pairwise_real_scores.append(rfs)

        # Filter out None entries (decoys that were skipped)
        valid = [(s, r) for s, r in zip(decoy_scores, decoy_reasons_out) if s is not None]
        if not valid:
            continue
        decoy_scores_clean = [v[0] for v in valid]
        decoy_reasons_clean = [v[1] for v in valid]

        # Pairwise ranking
        n_wins = sum(1 for d in pairwise_deltas if d > 0)
        all_beat = int(n_wins == len(pairwise_deltas))
        n_losses = sum(1 for d in pairwise_deltas if d <= 0)
        real_rank = 1 + n_losses
        mean_delta = np.mean(pairwise_deltas) if pairwise_deltas else 0.0
        real_score_full = np.mean(pairwise_real_scores) if pairwise_real_scores else 0.0

        # Extract metadata
        blocks = _detect_blocks(toks)
        metal = "?"
        cn = 0
        if blocks["metal_range"] is not None:
            parsed = _parse_metal_token(toks[blocks["metal_range"][0]])
            if parsed:
                metal, cn = parsed

        mrow = _metal_row(metal)
        n_ligs = len(blocks["ligand_blocks"])
        has_stereo = len(blocks.get("constraint_ranges", [])) > 0

        best_decoy_score = max(decoy_scores) if decoy_scores else float("-inf")
        pairwise_win_rate = n_wins / len(pairwise_deltas) if pairwise_deltas else 0.0

        rec = {
            "complex_id": cid,
            "metal": metal,
            "CN": cn,
            "metal_row": mrow,
            "n_ligands": n_ligs,
            "has_stereo": has_stereo,
            "real_score": round(real_score_full, 4),
            "best_decoy_score": round(best_decoy_score, 4),
            "delta_score": round(mean_delta, 4),
            "real_rank": real_rank,
            "n_decoys": len(decoys),
            "n_wins": n_wins,
            "win_rate": round(pairwise_win_rate, 4),
            "top1_real": all_beat,
            "mrr": round(1.0 / real_rank, 4),
        }

        # Baseline scores
        for bl_name, bl in baselines.items():
            bl_real = bl.score(toks)
            bl_decoy_scores = [bl.score(dt) for dt, _ in decoys]
            bl_all = [(bl_real, "real")] + [(s, "decoy") for s in bl_decoy_scores]
            bl_all.sort(key=lambda x: -x[0])
            bl_rank = next(i + 1 for i, (_, r) in enumerate(bl_all) if r == "real")
            rec[f"bl_{bl_name}_rank"] = bl_rank
            rec[f"bl_{bl_name}_top1"] = int(bl_rank == 1)
            rec[f"bl_{bl_name}_mrr"] = round(1.0 / bl_rank, 4)

        per_complex.append(rec)

        # Collect for AUROC: use pairwise deltas as scores for binary classification
        # label=1 → real field, label=0 → decoy field
        for i_d, (ds, reason, pw_real, pw_delta) in enumerate(
            zip(decoy_scores, decoy_reasons, pairwise_real_scores, pairwise_deltas)
        ):
            all_labels.append(1)
            all_scores_flat.append(pw_real)
            all_labels.append(0)
            all_scores_flat.append(ds)

            dtype = reason.split(":")[0]
            decoy_type_scores[dtype]["real"].append(pw_real)
            decoy_type_scores[dtype]["decoy"].append(ds)

        # Fieldwise for first 20 complexes
        if idx < 20:
            fw = compute_fieldwise_scores(model, tokenizer, toks, device)
            fw["complex_id"] = cid
            fieldwise_examples.append(fw)

    # ── aggregate metrics ────────────────────────────────────────
    n = len(per_complex)
    if n == 0:
        return {"per_complex": [], "summary": {}, "by_decoy_type": {},
                "fieldwise_examples": [], "anomalies": []}

    top1 = sum(r["top1_real"] for r in per_complex) / n
    mrr = sum(r["mrr"] for r in per_complex) / n
    mean_delta = sum(r["delta_score"] for r in per_complex) / n
    mean_win_rate = sum(r["win_rate"] for r in per_complex) / n
    auroc = _auroc(all_labels, all_scores_flat)

    summary = {
        "n_complexes": n,
        "n_decoys_total": sum(r["n_decoys"] for r in per_complex),
        "real_vs_decoys_top1": round(top1, 4),
        "real_vs_decoys_mrr": round(mrr, 4),
        "real_vs_decoy_auroc": round(auroc, 4),
        "mean_pairwise_win_rate": round(mean_win_rate, 4),
        "mean_delta_score": round(mean_delta, 4),
    }

    # Baseline summaries
    for bl_name in baselines:
        bl_top1 = sum(r.get(f"bl_{bl_name}_top1", 0) for r in per_complex) / n
        bl_mrr = sum(r.get(f"bl_{bl_name}_mrr", 0) for r in per_complex) / n
        summary[f"bl_{bl_name}_top1"] = round(bl_top1, 4)
        summary[f"bl_{bl_name}_mrr"] = round(bl_mrr, 4)

    # By decoy type AUROC
    by_decoy_type = {}
    for dtype in sorted(decoy_type_scores.keys()):
        d = decoy_type_scores[dtype]
        labels = [1] * len(d["real"]) + [0] * len(d["decoy"])
        scores = d["real"] + d["decoy"]
        auc = _auroc(labels, scores)
        by_decoy_type[dtype] = {
            "n_real": len(d["real"]),
            "n_decoy": len(d["decoy"]),
            "detection_auroc": round(auc, 4),
            "mean_real": round(np.mean(d["real"]), 4) if d["real"] else 0,
            "mean_decoy": round(np.mean(d["decoy"]), 4) if d["decoy"] else 0,
        }
        summary[f"{dtype}_detection_auc"] = round(auc, 4)

    # Anomaly mining: real records with lowest scores
    sorted_by_score = sorted(per_complex, key=lambda r: r["real_score"])
    n_anomaly = max(1, n // 100)  # bottom 1%
    anomalies = sorted_by_score[:max(n_anomaly, 50)]

    return {
        "per_complex": per_complex,
        "summary": summary,
        "by_decoy_type": by_decoy_type,
        "fieldwise_examples": fieldwise_examples,
        "anomalies": anomalies,
    }
