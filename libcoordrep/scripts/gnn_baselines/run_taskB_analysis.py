#!/usr/bin/env python3
"""
run_taskB_analysis.py
=====================
Three focused Task B analyses:
  1. Strict vs plausible decoy regrouping
  2. CoordRep-Ranker field ablation
  3. EGNN + metadata control experiment

No new large model training — reuses existing checkpoints.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import random as stdlib_random
import re
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from coordrep_tools.ranker_dataset import (
    _extract_unique_complexes, split_complexes, _tokens_to_ids,
)
from coordrep_tools.hard_decoys import (
    HardDecoyPool, generate_all_hard_decoys, _parse_metal_token,
)
from coordrep_tools.ablation_masking import _detect_blocks


# ── Constants ─────────────────────────────────────────────

FIG5_PATH = "/data/CoordRep/CoordSMILES/libcoordrep/figures/fig5/data/fig5_tasks.jsonl"
MLM_CKPT = "/data/CoordRep/CoordSMILES/libcoordrep/checkpoints/pretrain_v3/best_model.pt"
RANKER_CKPT = "/data/CoordRep/coordrep-release/libcoordrep/checkpoints/coordrep_ranker/best_finetuned.pt"
TOK_PATH = "/data/CoordRep/CoordSMILES/libcoordrep/checkpoints/pretrain_v3/tokenizer.json"

STRICT_TYPES = {'stereo_hard', 'metal_hard', 'boundary_hard'}
PLAUSIBLE_TYPES = {'ligand_hard', 'co_ligand_hard'}


# ── AUROC ─────────────────────────────────────────────────

def _auroc(labels, scores):
    if len(set(labels)) < 2:
        return 0.5
    pairs = sorted(zip(scores, labels), key=lambda x: -x[0])
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


# ── Shared data loading ──────────────────────────────────

def load_test_data(seed=42, n_decoys=20):
    """Load test complexes with pre-generated decoys."""
    complexes = _extract_unique_complexes(FIG5_PATH)
    train_cx, val_cx, test_cx = split_complexes(complexes, seed=seed)

    train_toks = [c['tokens'] for c in train_cx]
    pool = HardDecoyPool(train_toks, seed=seed)

    rng = stdlib_random.Random(seed)
    test_data = []
    for c in test_cx:
        toks = c['tokens']
        decoys = generate_all_hard_decoys(toks, pool, n_per_type=max(1, n_decoys // 5))
        if not decoys:
            continue
        if len(decoys) > n_decoys:
            rng.shuffle(decoys)
            decoys = decoys[:n_decoys]
        test_data.append({
            'id': c['id'], 'tokens': toks, 'decoys': decoys
        })
    return test_data, train_cx, pool


# ══════════════════════════════════════════════════════════
# ANALYSIS 1: Strict vs Plausible Regrouping
# ══════════════════════════════════════════════════════════

def run_analysis_1(test_data, device, out_dir):
    """Evaluate all methods on strict_semantic vs plausible_alternative subsets."""
    print("\n" + "="*70)
    print("ANALYSIS 1: Strict vs Plausible Decoy Regrouping")
    print("="*70)

    # Load tokenizer
    with open(TOK_PATH) as f:
        tok_data = json.load(f)
    from brain.tokenizer import CoordRepTokenizer, TokenizerConfig
    tok_config = TokenizerConfig(**tok_data.get("config", {}))
    tokenizer = CoordRepTokenizer(tok_config)
    tokenizer.token2id = tok_data["token2id"]
    tokenizer.id2token = {v: k for k, v in tokenizer.token2id.items()}
    pad_id = tokenizer.token2id.get("[PAD]", 0)

    # Load CoordRep-Ranker
    from coordrep_tools.coordrep_ranker import load_ranker
    ranker, _ = load_ranker(RANKER_CKPT, TOK_PATH, device=device)
    ranker.eval()

    # Load MLM scorer
    from coordrep_tools.infer import load_model_and_tokenizer
    from coordrep_tools.coordrep_score import _SPECIAL
    mlm_model, _ = load_model_and_tokenizer(MLM_CKPT, device=device)

    mask_id = tokenizer.token2id.get("[MASK]", 3)

    def mlm_pll_score(toks_a, toks_b):
        """Field-targeted PLL: mask differing positions, score both."""
        ids_a = _tokens_to_ids(toks_a, tokenizer)
        ids_b = _tokens_to_ids(toks_b, tokenizer)
        min_len = min(len(toks_a), len(toks_b), 512)
        diff_pos = [i for i in range(min_len)
                    if toks_a[i] != toks_b[i] and toks_a[i] not in _SPECIAL]
        if not diff_pos:
            return 0.0, 0.0
        ctx = list(ids_a[:512])
        for p in diff_pos:
            if p < len(ctx):
                ctx[p] = mask_id
        padded = ctx + [pad_id] * (512 - len(ctx)) if len(ctx) < 512 else ctx[:512]
        att = [1]*min(len(ctx), 512) + [0]*(512 - min(len(ctx), 512))
        inp = torch.tensor([padded], device=device)
        att_t = torch.tensor([att], device=device)
        with torch.no_grad():
            logits, _ = mlm_model(inp, attention_mask=att_t)
        lp = logits[0]
        real_ll, decoy_ll = 0.0, 0.0
        for p in diff_pos:
            if p >= lp.shape[0]:
                continue
            pr = torch.softmax(lp[p, :], dim=-1)
            real_ll += math.log(max(pr[ids_a[p]].item(), 1e-20))
            decoy_ll += math.log(max(pr[ids_b[p]].item(), 1e-20))
        n = max(len(diff_pos), 1)
        return real_ll / n, decoy_ll / n

    def _pad(ids, max_len=512):
        ids = ids[:max_len]
        att = [1]*len(ids) + [0]*(max_len - len(ids))
        ids = ids + [pad_id]*(max_len - len(ids))
        return ids, att

    def ranker_score_group(toks, decoy_list):
        """Score real + all decoys with CoordRep-Ranker."""
        real_ids = _tokens_to_ids(toks, tokenizer)
        r_pad, r_att = _pad(real_ids)
        all_ids = [r_pad]
        all_att = [r_att]
        for dtoks, _ in decoy_list:
            d_ids = _tokens_to_ids(dtoks, tokenizer)
            dp, da = _pad(d_ids)
            all_ids.append(dp)
            all_att.append(da)
        ids_t = torch.tensor(all_ids, device=device).unsqueeze(0)
        att_t = torch.tensor(all_att, device=device).unsqueeze(0)
        with torch.no_grad():
            scores = ranker.score_list(ids_t, att_t).squeeze(0)
        return scores[0].item(), scores[1:].cpu().numpy()

    # Evaluate per subset
    rng = stdlib_random.Random(42)

    # Accumulate per-method per-subset
    # methods: random, mlm_pll, coordrep_ranker
    subsets = {'strict': STRICT_TYPES, 'plausible': PLAUSIBLE_TYPES, 'all': STRICT_TYPES | PLAUSIBLE_TYPES}

    # Per (method, subset): lists
    metrics = defaultdict(lambda: defaultdict(lambda: {
        'top1': [], 'mrr': [], 'win_rate': [],
        'real_scores': [], 'decoy_scores': [],
        'type_real': defaultdict(list), 'type_decoy': defaultdict(list),
    }))

    print(f"Scoring {len(test_data)} test groups …")
    for ti, td in enumerate(test_data):
        if ti % 200 == 0:
            print(f"  [{ti}/{len(test_data)}]", flush=True)

        toks = td['tokens']
        all_decoys = td['decoys']

        # Classify decoys
        decoys_by_subset = {'strict': [], 'plausible': [], 'all': all_decoys}
        for dtoks, reason in all_decoys:
            dtype = reason.split(':')[0]
            if dtype in STRICT_TYPES:
                decoys_by_subset['strict'].append((dtoks, reason))
            elif dtype in PLAUSIBLE_TYPES:
                decoys_by_subset['plausible'].append((dtoks, reason))

        # CoordRep-Ranker: score all at once
        if all_decoys:
            real_s, all_decoy_s = ranker_score_group(toks, all_decoys)
        else:
            continue

        # MLM PLL: score each pair
        mlm_real_scores = []
        mlm_decoy_scores = []
        for dtoks, reason in all_decoys:
            rs, ds = mlm_pll_score(toks, dtoks)
            mlm_real_scores.append(rs)
            mlm_decoy_scores.append(ds)

        # Now compute metrics per subset
        for subset_name, subset_decoys in decoys_by_subset.items():
            if not subset_decoys:
                continue

            # Get indices into all_decoys for this subset
            subset_indices = []
            for di, (dtoks, reason) in enumerate(all_decoys):
                dtype = reason.split(':')[0]
                if subset_name == 'all' or dtype in subsets[subset_name]:
                    subset_indices.append(di)

            if not subset_indices:
                continue

            # Random
            rnd_real = rng.random()
            rnd_decoy = [rng.random() for _ in subset_indices]
            rnd_wins = sum(1 for d in rnd_decoy if rnd_real > d)
            rnd_rank = 1 + sum(1 for d in rnd_decoy if d >= rnd_real)
            m = metrics['Random'][subset_name]
            m['top1'].append(int(rnd_rank == 1))
            m['mrr'].append(1.0 / rnd_rank)
            m['win_rate'].append(rnd_wins / len(subset_indices))

            # CoordRep-Ranker
            sub_decoy_s = all_decoy_s[subset_indices]
            cr_wins = int(np.sum(sub_decoy_s < real_s))
            cr_rank = 1 + int(np.sum(sub_decoy_s >= real_s))
            m = metrics['CoordRep_Ranker'][subset_name]
            m['top1'].append(int(cr_rank == 1))
            m['mrr'].append(1.0 / cr_rank)
            m['win_rate'].append(cr_wins / len(subset_indices))
            m['real_scores'].append(real_s)
            for di in subset_indices:
                m['decoy_scores'].append(float(all_decoy_s[di]))
                dtype = all_decoys[di][1].split(':')[0]
                m['type_real'][dtype].append(real_s)
                m['type_decoy'][dtype].append(float(all_decoy_s[di]))

            # MLM PLL
            sub_mlm_real = [mlm_real_scores[di] for di in subset_indices]
            sub_mlm_decoy = [mlm_decoy_scores[di] for di in subset_indices]
            mean_real_mlm = np.mean(sub_mlm_real)
            mlm_wins = sum(1 for d in sub_mlm_decoy if mean_real_mlm > d)
            mlm_rank = 1 + sum(1 for d in sub_mlm_decoy if d >= mean_real_mlm)
            m = metrics['CoordRep_MLM_PLL'][subset_name]
            m['top1'].append(int(mlm_rank == 1))
            m['mrr'].append(1.0 / mlm_rank)
            m['win_rate'].append(mlm_wins / len(subset_indices))
            m['real_scores'].extend(sub_mlm_real)
            m['decoy_scores'].extend(sub_mlm_decoy)
            for di in subset_indices:
                dtype = all_decoys[di][1].split(':')[0]
                m['type_real'][dtype].append(mlm_real_scores[di])
                m['type_decoy'][dtype].append(mlm_decoy_scores[di])

    # Compute summary
    rows = []
    summary_json = {}
    for method in ['Random', 'CoordRep_MLM_PLL', 'CoordRep_Ranker']:
        for subset in ['strict', 'plausible', 'all']:
            m = metrics[method][subset]
            if not m['top1']:
                continue
            overall_auroc = _auroc(
                [1]*len(m['real_scores']) + [0]*len(m['decoy_scores']),
                list(m['real_scores']) + list(m['decoy_scores'])
            ) if m['real_scores'] and m['decoy_scores'] else 0.5

            row = {
                'method': method,
                'subset': subset,
                'n_groups': len(m['top1']),
                'top1': round(float(np.mean(m['top1'])), 4),
                'mrr': round(float(np.mean(m['mrr'])), 4),
                'win_rate': round(float(np.mean(m['win_rate'])), 4),
                'auroc': round(overall_auroc, 4),
            }

            # Per-type AUROC
            for dtype in sorted(STRICT_TYPES | PLAUSIBLE_TYPES):
                tr = m['type_real'].get(dtype, [])
                td_s = m['type_decoy'].get(dtype, [])
                if tr and td_s:
                    ta = _auroc([1]*len(tr) + [0]*len(td_s), tr + td_s)
                    row[f'auroc_{dtype}'] = round(ta, 4)
                else:
                    row[f'auroc_{dtype}'] = ''

            rows.append(row)
            summary_json[f'{method}_{subset}'] = row

    # Print
    print(f"\n{'Method':<22} {'Subset':<12} {'Top1':>7} {'MRR':>7} {'WinRate':>8} {'AUROC':>7}")
    print("-" * 70)
    for r in rows:
        print(f"{r['method']:<22} {r['subset']:<12} {r['top1']:>7.4f} "
              f"{r['mrr']:>7.4f} {r['win_rate']:>8.4f} {r['auroc']:>7.4f}")

    # Write CSV
    p = os.path.join(out_dir, "taskB_strict_vs_plausible.csv")
    fields = ['method', 'subset', 'n_groups', 'top1', 'mrr', 'win_rate', 'auroc'] + \
             [f'auroc_{dt}' for dt in sorted(STRICT_TYPES | PLAUSIBLE_TYPES)]
    with open(p, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"\n  {p}")

    # Write JSON
    p = os.path.join(out_dir, "taskB_strict_summary.json")
    with open(p, 'w') as f:
        json.dump(summary_json, f, indent=2, default=str)
    print(f"  {p}")

    return summary_json


# ══════════════════════════════════════════════════════════
# ANALYSIS 2: CoordRep-Ranker Field Ablation
# ══════════════════════════════════════════════════════════

def _mask_metal_field(tokens):
    """Replace metal token with [MASK]."""
    out = list(tokens)
    for i, t in enumerate(out):
        if t.startswith('[Metal:'):
            out[i] = '[MASK]'
    return out

def _mask_stereo_tokens(tokens):
    """Remove/mask constraint blocks ({trans:...}, {cis:...})."""
    blocks = _detect_blocks(tokens)
    out = list(tokens)
    for cs, ce in blocks['constraint_ranges']:
        for j in range(cs, ce):
            out[j] = '[MASK]'
    return out

def _mask_shape_boundary(tokens):
    """Remove/mask shape block (<Shape: ... >)."""
    blocks = _detect_blocks(tokens)
    out = list(tokens)
    if blocks['shape_range']:
        s, e = blocks['shape_range']
        for j in range(s, e):
            out[j] = '[MASK]'
    return out

def _mask_ligand_blocks(tokens):
    """Mask all ligand SMILES content."""
    blocks = _detect_blocks(tokens)
    out = list(tokens)
    for lb in blocks['ligand_blocks']:
        if lb['smiles_start'] is not None and lb['smiles_end'] is not None:
            for j in range(lb['smiles_start'], lb['smiles_end']):
                out[j] = '[MASK]'
    return out


def run_analysis_2(test_data, device, out_dir):
    """CoordRep-Ranker field ablation."""
    print("\n" + "="*70)
    print("ANALYSIS 2: CoordRep-Ranker Field Ablation")
    print("="*70)

    with open(TOK_PATH) as f:
        tok_data = json.load(f)
    from brain.tokenizer import CoordRepTokenizer, TokenizerConfig
    tok_config = TokenizerConfig(**tok_data.get("config", {}))
    tokenizer = CoordRepTokenizer(tok_config)
    tokenizer.token2id = tok_data["token2id"]
    tokenizer.id2token = {v: k for k, v in tokenizer.token2id.items()}
    pad_id = tokenizer.token2id.get("[PAD]", 0)

    from coordrep_tools.coordrep_ranker import load_ranker
    ranker, _ = load_ranker(RANKER_CKPT, TOK_PATH, device=device)
    ranker.eval()

    def _pad(ids, max_len=512):
        ids = ids[:max_len]
        att = [1]*len(ids) + [0]*(max_len - len(ids))
        ids = ids + [pad_id]*(max_len - len(ids))
        return ids, att

    ablation_configs = {
        'full': lambda t: t,
        'no_stereo': _mask_stereo_tokens,
        'no_metal': _mask_metal_field,
        'no_shape': _mask_shape_boundary,
        'no_ligand': _mask_ligand_blocks,
    }

    # Per ablation: accumulate per-type scores
    abl_metrics = {}
    for abl_name in ablation_configs:
        abl_metrics[abl_name] = {
            'top1': [], 'mrr': [], 'win_rate': [],
            'real_scores': [], 'decoy_scores': [],
            'type_real': defaultdict(list), 'type_decoy': defaultdict(list),
        }

    print(f"Scoring {len(test_data)} groups × {len(ablation_configs)} ablations …")

    for ti, td in enumerate(test_data):
        if ti % 200 == 0:
            print(f"  [{ti}/{len(test_data)}]", flush=True)

        toks = td['tokens']
        decoys = td['decoys']
        if not decoys:
            continue

        for abl_name, abl_fn in ablation_configs.items():
            # Ablate real and all decoys
            abl_real = abl_fn(toks)
            real_ids = _tokens_to_ids(abl_real, tokenizer)
            r_pad, r_att = _pad(real_ids)

            all_ids = [r_pad]
            all_att = [r_att]
            for dtoks, reason in decoys:
                abl_d = abl_fn(dtoks)
                d_ids = _tokens_to_ids(abl_d, tokenizer)
                dp, da = _pad(d_ids)
                all_ids.append(dp)
                all_att.append(da)

            ids_t = torch.tensor(all_ids, device=device).unsqueeze(0)
            att_t = torch.tensor(all_att, device=device).unsqueeze(0)
            with torch.no_grad():
                scores = ranker.score_list(ids_t, att_t).squeeze(0)

            real_s = scores[0].item()
            decoy_s = scores[1:].cpu().numpy()

            wins = int(np.sum(decoy_s < real_s))
            rank = 1 + int(np.sum(decoy_s >= real_s))

            m = abl_metrics[abl_name]
            m['top1'].append(int(rank == 1))
            m['mrr'].append(1.0 / rank)
            m['win_rate'].append(wins / len(decoy_s))
            m['real_scores'].append(real_s)

            for di, (dtoks, reason) in enumerate(decoys):
                dtype = reason.split(':')[0]
                m['decoy_scores'].append(float(decoy_s[di]))
                m['type_real'][dtype].append(real_s)
                m['type_decoy'][dtype].append(float(decoy_s[di]))

    # Aggregate
    rows = []
    summary_json = {}
    all_types = sorted(STRICT_TYPES | PLAUSIBLE_TYPES)

    for abl_name in ablation_configs:
        m = abl_metrics[abl_name]
        if not m['top1']:
            continue

        overall_auroc = _auroc(
            [1]*len(m['real_scores']) + [0]*len(m['decoy_scores']),
            list(m['real_scores']) + list(m['decoy_scores'])
        )

        row = {
            'ablation': abl_name,
            'n_groups': len(m['top1']),
            'top1': round(float(np.mean(m['top1'])), 4),
            'mrr': round(float(np.mean(m['mrr'])), 4),
            'win_rate': round(float(np.mean(m['win_rate'])), 4),
            'auroc': round(overall_auroc, 4),
        }

        for dtype in all_types:
            tr = m['type_real'].get(dtype, [])
            td_s = m['type_decoy'].get(dtype, [])
            if tr and td_s:
                row[f'auroc_{dtype}'] = round(_auroc(
                    [1]*len(tr) + [0]*len(td_s), tr + td_s), 4)
            else:
                row[f'auroc_{dtype}'] = ''

        rows.append(row)
        summary_json[abl_name] = row

    # Print
    print(f"\n{'Ablation':<15} {'Top1':>7} {'MRR':>7} {'WinRate':>8} {'AUROC':>7} "
          f"{'stereo':>8} {'metal':>8} {'boundary':>8}")
    print("-" * 85)
    for r in rows:
        st = r.get('auroc_stereo_hard', '')
        mt = r.get('auroc_metal_hard', '')
        bd = r.get('auroc_boundary_hard', '')
        print(f"{r['ablation']:<15} {r['top1']:>7.4f} {r['mrr']:>7.4f} "
              f"{r['win_rate']:>8.4f} {r['auroc']:>7.4f} "
              f"{st!s:>8} {mt!s:>8} {bd!s:>8}")

    # Write CSV
    p = os.path.join(out_dir, "taskB_coordrep_field_ablation.csv")
    fields = ['ablation', 'n_groups', 'top1', 'mrr', 'win_rate', 'auroc'] + \
             [f'auroc_{dt}' for dt in all_types]
    with open(p, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"\n  {p}")

    p = os.path.join(out_dir, "taskB_field_ablation_summary.json")
    with open(p, 'w') as f:
        json.dump(summary_json, f, indent=2, default=str)
    print(f"  {p}")

    return summary_json


# ══════════════════════════════════════════════════════════
# ANALYSIS 3: EGNN + Metadata Control
# ══════════════════════════════════════════════════════════

def run_analysis_3(test_data, train_cx, pool, device, out_dir):
    """EGNN embedding + metadata MLP control experiment."""
    print("\n" + "="*70)
    print("ANALYSIS 3: EGNN + Metadata Control")
    print("="*70)

    from coordrep_tools.gnn_baselines.datasets import parse_tmqm_xyz, TMQM_XYZ_DIR
    from coordrep_tools.gnn_baselines.datasets import one_hot
    from coordrep_tools.gnn_baselines.models import EGNNLayer
    from torch_geometric.data import Data, Batch
    from torch_geometric.nn import global_mean_pool
    from scipy.spatial.distance import cdist

    ELEMENTS_3D = ['H', 'C', 'N', 'O', 'S', 'P', 'F', 'Cl', 'Br', 'Se',
                   'Si', 'B', 'I', 'Fe', 'Cu', 'Ni', 'Co', 'Zn', 'Mn',
                   'Pd', 'Pt', 'Ru', 'Rh', 'Ir', 'Cr', 'V', 'Ti', 'other']
    TRANSITION_METALS = {
        'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
        'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
        'La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
    }
    METAL_LIST = ['Fe', 'Cu', 'Ni', 'Co', 'Zn', 'Mn', 'Cr', 'V', 'Ti',
                  'Pd', 'Pt', 'Ru', 'Rh', 'Ir', 'Os', 'Ag', 'Au', 'Cd',
                  'Hg', 'Sc', 'Mo', 'W', 'Re', 'Ta', 'Nb', 'Zr', 'Hf',
                  'Y', 'La', 'other']
    SHAPE_LABELS = ['L', 'Td', 'SP', 'Oh', 'TP', 'TBP', 'SPY', 'TPr', 'other']
    STEREO_CLASSES = ['none', 'trans', 'cis', 'fac', 'mer', 'other']

    def get_metal_from_tokens(tokens):
        if tokens and tokens[0].startswith('[Metal:'):
            m = re.search(r'\[Metal:([A-Z][a-z]?)', tokens[0])
            if m:
                return m.group(1)
        return None

    def get_cn_from_tokens(tokens):
        if tokens and tokens[0].startswith('[Metal:'):
            m = re.search(r'CN:(\d+)', tokens[0])
            if m:
                return int(m.group(1))
        return 0

    def get_shape_from_tokens(tokens):
        blocks = _detect_blocks(tokens)
        if blocks['shape_range']:
            s, e = blocks['shape_range']
            for t in tokens[s:e]:
                if t in {'Td', 'SP', 'Oh', 'TP', 'TBP', 'SPY', 'TPr', 'L'}:
                    return t
        return 'other'

    def get_stereo_from_tokens(tokens):
        blocks = _detect_blocks(tokens)
        for cs, ce in blocks['constraint_ranges']:
            for t in tokens[cs:ce]:
                if t.startswith('{trans:'):
                    return 'trans'
                if t.startswith('{cis:'):
                    return 'cis'
        # Check fac/mer in tokens
        for t in tokens:
            if 'fac' in t.lower():
                return 'fac'
            if 'mer' in t.lower():
                return 'mer'
        return 'none'

    def get_donor_pattern(tokens):
        """Extract donor elements as sorted tuple."""
        blocks = _detect_blocks(tokens)
        donors = []
        for i, t in enumerate(tokens):
            if re.match(r'^:[A-Z][a-z]?:\d+$', t):
                elem = t.split(':')[1]
                donors.append(elem)
        return tuple(sorted(donors)) if donors else ('unknown',)

    def encode_metadata(tokens):
        """Encode metadata as a fixed-length vector."""
        metal = get_metal_from_tokens(tokens) or 'other'
        cn = get_cn_from_tokens(tokens)
        shape = get_shape_from_tokens(tokens)
        stereo = get_stereo_from_tokens(tokens)
        donors = get_donor_pattern(tokens)

        feat = []
        # Metal one-hot
        feat.extend(one_hot(metal, METAL_LIST))
        # CN (normalised)
        feat.append(cn / 10.0)
        # Shape one-hot
        feat.extend(one_hot(shape, SHAPE_LABELS))
        # Stereo one-hot
        feat.extend(one_hot(stereo, STEREO_CLASSES))
        # Donor pattern: count of common donor elements
        for elem in ['N', 'O', 'S', 'P', 'C', 'Cl']:
            feat.append(donors.count(elem) / 6.0)
        # Total donors
        feat.append(len(donors) / 10.0)

        return torch.tensor(feat, dtype=torch.float)

    meta_dim = len(encode_metadata(['[Metal:Fe|CN:6]']))  # probe dim

    # Build tmqm_N -> CSD mapping (lightweight)
    idx_to_csd = {}
    _idx = 0
    for _xf in sorted(['tmQM_X1.xyz', 'tmQM_X2.xyz', 'tmQM_X3.xyz']):
        _xp = os.path.join(TMQM_XYZ_DIR, _xf)
        if not os.path.exists(_xp):
            continue
        with open(_xp) as _f:
            while True:
                _line = _f.readline()
                if not _line:
                    break
                try:
                    _na = int(_line.strip())
                except ValueError:
                    continue
                _comment = _f.readline().strip()
                _m = re.search(r'CSD_code\s*=\s*(\S+)', _comment)
                idx_to_csd[_idx] = _m.group(1) if _m else f'unk_{_idx}'
                for _ in range(_na):
                    _f.readline()
                _idx += 1

    def tmqm_to_csd(tmqm_id):
        m = re.match(r'tmqm_(\d+)', tmqm_id)
        return idx_to_csd.get(int(m.group(1))) if m else None

    # Collect needed CSD codes
    all_ids = set()
    for c in train_cx[:3000]:
        csd = tmqm_to_csd(c['id'])
        if csd:
            all_ids.add(csd)
    for td in test_data:
        csd = tmqm_to_csd(td['id'])
        if csd:
            all_ids.add(csd)

    print(f"Parsing XYZ for {len(all_ids)} molecules …")
    xyz_data = {}
    for xf in ['tmQM_X1.xyz', 'tmQM_X2.xyz', 'tmQM_X3.xyz']:
        xp = os.path.join(TMQM_XYZ_DIR, xf)
        if os.path.exists(xp):
            xyz_data.update(parse_tmqm_xyz(xp, all_ids))
    print(f"  Loaded {len(xyz_data)}")

    def _radius_edges(coords, cutoff=4.0):
        D = cdist(coords, coords)
        mask = (D < cutoff) & (D > 0)
        src, dst = np.where(mask)
        return torch.tensor(np.stack([src, dst]), dtype=torch.long) if len(src) else torch.zeros((2, 0), dtype=torch.long)

    def build_graph(mol_data, metal_elem):
        atoms = mol_data['atoms']
        coords = np.array(mol_data['coords'], dtype=np.float32)
        n = len(atoms)
        if n == 0 or n > 300:
            return None
        metal_idx = None
        for i, a in enumerate(atoms):
            if a == metal_elem or a in TRANSITION_METALS:
                metal_idx = i
                break
        if metal_idx is None:
            return None
        x_list = []
        for i, a in enumerate(atoms):
            feat = one_hot(a, ELEMENTS_3D)
            feat.append(1.0 if i == metal_idx else 0.0)
            feat.append(1.0 if a in TRANSITION_METALS else 0.0)
            x_list.append(feat)
        x = torch.tensor(x_list, dtype=torch.float)
        edge_index = _radius_edges(coords)
        pos = torch.tensor(coords, dtype=torch.float)
        return Data(x=x, edge_index=edge_index, pos=pos, num_nodes=n)

    # Model: EGNN encoder + metadata MLP → scalar score
    class EGNNMetaRanker(nn.Module):
        def __init__(self, atom_dim, meta_dim, hidden=128, n_layers=4):
            super().__init__()
            self.embed = nn.Linear(atom_dim, hidden)
            self.layers = nn.ModuleList([EGNNLayer(hidden) for _ in range(n_layers)])
            self.meta_proj = nn.Linear(meta_dim, hidden)
            self.head = nn.Sequential(
                nn.Linear(hidden * 2, hidden),
                nn.SiLU(),
                nn.Linear(hidden, 1),
            )
            self.dropout = 0.1

        def forward(self, x, edge_index, pos, batch, meta):
            h = self.embed(x)
            for layer in self.layers:
                h, pos = layer(h, pos, edge_index)
                h = F.dropout(h, p=self.dropout, training=self.training)
            pooled = global_mean_pool(h, batch)
            meta_h = F.silu(self.meta_proj(meta))
            combined = torch.cat([pooled, meta_h], dim=-1)
            return self.head(combined).squeeze(-1)

    # Also build EGNN-only (no metadata) as control
    class EGNNOnlyRanker(nn.Module):
        def __init__(self, atom_dim, hidden=128, n_layers=4):
            super().__init__()
            self.embed = nn.Linear(atom_dim, hidden)
            self.layers = nn.ModuleList([EGNNLayer(hidden) for _ in range(n_layers)])
            self.head = nn.Sequential(
                nn.Linear(hidden, hidden),
                nn.SiLU(),
                nn.Linear(hidden, 1),
            )
            self.dropout = 0.1

        def forward(self, x, edge_index, pos, batch, meta=None):
            h = self.embed(x)
            for layer in self.layers:
                h, pos = layer(h, pos, edge_index)
                h = F.dropout(h, p=self.dropout, training=self.training)
            pooled = global_mean_pool(h, batch)
            return self.head(pooled).squeeze(-1)

    # Build training data
    print("Building training groups …")
    COORD_CONSISTENT_TYPES = {'stereo_hard', 'metal_hard', 'boundary_hard'}
    train_groups = []
    for c in train_cx[:2000]:
        csd = tmqm_to_csd(c['id'])
        if csd is None or csd not in xyz_data:
            continue
        metal = get_metal_from_tokens(c['tokens'])
        if metal is None:
            continue
        g = build_graph(xyz_data[csd], metal)
        if g is None:
            continue
        meta_real = encode_metadata(c['tokens'])
        decoys = generate_all_hard_decoys(c['tokens'], pool, n_per_type=2)
        filtered = [(dt, r) for dt, r in decoys if r.split(':')[0] in COORD_CONSISTENT_TYPES]
        if not filtered:
            continue
        decoy_metas = []
        decoy_graphs = []
        for dt, reason in filtered[:8]:
            dtype = reason.split(':')[0]
            dm = get_metal_from_tokens(dt) if dtype == 'metal_hard' else metal
            dg = build_graph(xyz_data[csd], dm or metal)
            if dg is not None:
                decoy_graphs.append(dg)
                decoy_metas.append(encode_metadata(dt))
        if decoy_graphs:
            train_groups.append((g, meta_real, decoy_graphs, decoy_metas))

    print(f"  Train groups: {len(train_groups)}")
    if not train_groups:
        print("  SKIP: no training data")
        return None

    atom_dim = train_groups[0][0].x.size(1)

    # Train both models
    results = {}
    for model_name, ModelClass in [('EGNN_only', EGNNOnlyRanker), ('EGNN_meta', EGNNMetaRanker)]:
        print(f"\nTraining {model_name} …")
        if model_name == 'EGNN_meta':
            model = ModelClass(atom_dim, meta_dim, hidden=128, n_layers=4).to(device)
        else:
            model = ModelClass(atom_dim, hidden=128, n_layers=4).to(device)

        optimizer = torch.optim.Adam(model.parameters(), lr=1e-3, weight_decay=1e-5)
        best_loss = float('inf')
        best_state = None

        for epoch in range(10):
            model.train()
            total_loss = 0
            nb = 0
            stdlib_random.shuffle(train_groups)
            for real_g, meta_real, decoy_gs, decoy_metas in train_groups:
                all_gs = [real_g] + decoy_gs
                all_metas = [meta_real] + decoy_metas
                batch = Batch.from_data_list(all_gs).to(device)
                meta_t = torch.stack(all_metas).to(device)

                scores = model(batch.x, batch.edge_index, batch.pos, batch.batch, meta_t)
                target = torch.zeros(1, dtype=torch.long, device=device)
                loss = F.cross_entropy(scores.unsqueeze(0), target)
                optimizer.zero_grad()
                loss.backward()
                torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
                optimizer.step()
                total_loss += loss.item()
                nb += 1

            avg = total_loss / max(nb, 1)
            if avg < best_loss:
                best_loss = avg
                best_state = {k: v.clone() for k, v in model.state_dict().items()}
            if (epoch + 1) % 5 == 0:
                print(f"    Epoch {epoch+1}: loss={avg:.4f}")

        if best_state:
            model.load_state_dict(best_state)

        # Evaluate
        model.eval()
        m_results = {
            'top1': [], 'mrr': [], 'win_rate': [],
            'real_scores': [], 'decoy_scores': [],
            'type_real': defaultdict(list), 'type_decoy': defaultdict(list),
        }

        with torch.no_grad():
            for td in test_data[:1501]:
                csd = tmqm_to_csd(td['id'])
                if csd is None or csd not in xyz_data:
                    continue
                metal = get_metal_from_tokens(td['tokens'])
                if metal is None:
                    continue
                real_g = build_graph(xyz_data[csd], metal)
                if real_g is None:
                    continue
                meta_real = encode_metadata(td['tokens'])

                # Score real
                batch_r = Batch.from_data_list([real_g]).to(device)
                meta_r = meta_real.unsqueeze(0).to(device)
                real_s = model(batch_r.x, batch_r.edge_index, batch_r.pos, batch_r.batch, meta_r).item()

                decoy_scores = []
                for dtoks, reason in td['decoys']:
                    dtype = reason.split(':')[0]
                    if dtype not in COORD_CONSISTENT_TYPES:
                        continue
                    dm = get_metal_from_tokens(dtoks) if dtype == 'metal_hard' else metal
                    dg = build_graph(xyz_data[csd], dm or metal)
                    if dg is None:
                        continue
                    meta_d = encode_metadata(dtoks)
                    batch_d = Batch.from_data_list([dg]).to(device)
                    meta_dt = meta_d.unsqueeze(0).to(device)
                    ds = model(batch_d.x, batch_d.edge_index, batch_d.pos, batch_d.batch, meta_dt).item()
                    decoy_scores.append((ds, reason))

                if not decoy_scores:
                    continue

                wins = sum(1 for ds, _ in decoy_scores if real_s > ds)
                rank = 1 + sum(1 for ds, _ in decoy_scores if ds >= real_s)

                m_results['top1'].append(int(rank == 1))
                m_results['mrr'].append(1.0 / rank)
                m_results['win_rate'].append(wins / len(decoy_scores))
                m_results['real_scores'].append(real_s)
                for ds, reason in decoy_scores:
                    m_results['decoy_scores'].append(ds)
                    dtype = reason.split(':')[0]
                    m_results['type_real'][dtype].append(real_s)
                    m_results['type_decoy'][dtype].append(ds)

        overall_auroc = _auroc(
            [1]*len(m_results['real_scores']) + [0]*len(m_results['decoy_scores']),
            m_results['real_scores'] + m_results['decoy_scores']
        ) if m_results['real_scores'] else 0.5

        row = {
            'model': model_name,
            'n': len(m_results['top1']),
            'top1': round(float(np.mean(m_results['top1'])), 4) if m_results['top1'] else 0,
            'mrr': round(float(np.mean(m_results['mrr'])), 4) if m_results['mrr'] else 0,
            'win_rate': round(float(np.mean(m_results['win_rate'])), 4) if m_results['win_rate'] else 0,
            'auroc': round(overall_auroc, 4),
        }
        for dtype in COORD_CONSISTENT_TYPES:
            tr = m_results['type_real'].get(dtype, [])
            td_s = m_results['type_decoy'].get(dtype, [])
            if tr and td_s:
                row[f'auroc_{dtype}'] = round(_auroc([1]*len(tr) + [0]*len(td_s), tr + td_s), 4)
            else:
                row[f'auroc_{dtype}'] = ''

        results[model_name] = row
        print(f"  {model_name}: Top1={row['top1']} MRR={row['mrr']} "
              f"WinRate={row['win_rate']} AUROC={row['auroc']}")

        del model
        torch.cuda.empty_cache()

    # Write
    p = os.path.join(out_dir, "taskB_egnn_metadata_control.csv")
    all_rows = list(results.values())
    if all_rows:
        fields = list(all_rows[0].keys())
        with open(p, 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=fields)
            w.writeheader()
            for r in all_rows:
                w.writerow(r)
    print(f"\n  {p}")

    return results


# ── Main ─────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--device", default="cuda:0")
    parser.add_argument("--out", default="revision_results/gnn_baselines")
    parser.add_argument("--skip_1", action="store_true")
    parser.add_argument("--skip_2", action="store_true")
    parser.add_argument("--skip_3", action="store_true")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    os.makedirs(args.out, exist_ok=True)
    device = torch.device(args.device if torch.cuda.is_available() else 'cpu')

    print("Loading test data …")
    test_data, train_cx, pool = load_test_data(seed=args.seed)
    print(f"  {len(test_data)} test groups")

    a1_results = None
    a2_results = None
    a3_results = None

    if not args.skip_1:
        a1_results = run_analysis_1(test_data, device, args.out)

    if not args.skip_2:
        a2_results = run_analysis_2(test_data, device, args.out)

    if not args.skip_3:
        a3_results = run_analysis_3(test_data, train_cx, pool, device, args.out)

    print(f"\n{'='*70}")
    print("All analyses complete.")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
