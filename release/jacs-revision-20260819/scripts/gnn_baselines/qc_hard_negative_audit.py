#!/usr/bin/env python3
"""
qc_hard_negative_audit.py
=========================
QC audit of hard-negative benchmark metrics.

Answers:
1. Are GIN and CoordRep-Ranker on the same test groups? → YES (same split_complexes seed=42)
2. Candidates per group: 1 real + N decoys
3. Top-1 = listwise real-ranked-first (rank==1)
4. 78% in README was WRONG — it was pairwise win_rate (74.3%), not top-1 (13.1%)
5. Only stereo_hard? No — all 5 decoy types
6. Train/val/test leakage? No — same split fn
7. Output 20 groups with real_score, decoy_scores, real_rank
8. Generate hard_negative_summary_fixed.csv with all required fields

Evaluates GIN Ranker on SAME test groups as CoordRep-Ranker (full 1501).
"""

from __future__ import annotations

import csv
import json
import math
import os
import random as stdlib_random
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import List

import numpy as np
import torch
import torch.nn.functional as F
from torch_geometric.data import Batch

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from coordrep_tools.hard_decoys import (
    HardDecoyPool, generate_all_hard_decoys, HARD_DECOY_TYPES,
    _parse_metal_token,
)
from coordrep_tools.ablation_masking import _detect_blocks
from coordrep_tools.ranker_dataset import (
    _extract_unique_complexes, split_complexes,
)
from coordrep_tools.gnn_baselines.models import GINRanker

# Re-use the graph builder from the ranker script
sys.path.insert(0, str(Path(__file__).parent))
from run_hard_negative_ranker import (
    tokens_to_simple_graph, RandomBaseline, FrequencyBaseline,
    WLHashBaseline, RuleValidatorBaseline,
)


def compute_auroc(real_scores, decoy_scores):
    """Compute AUROC from lists of real and decoy scores."""
    labels = [1] * len(real_scores) + [0] * len(decoy_scores)
    scores = list(real_scores) + list(decoy_scores)
    if len(set(labels)) < 2 or not scores:
        return 0.5
    # Sort by score descending
    pairs = sorted(zip(scores, labels), key=lambda x: -x[0])
    n_pos = sum(labels)
    n_neg = len(labels) - n_pos
    if n_pos == 0 or n_neg == 0:
        return 0.5
    tp = 0
    auc = 0.0
    for s, l in pairs:
        if l == 1:
            tp += 1
        else:
            auc += tp
    return auc / (n_pos * n_neg)


def main():
    out_dir = "revision_results/gnn_baselines"
    os.makedirs(out_dir, exist_ok=True)

    fig5_path = "/data/CoordRep/CoordSMILES/libcoordrep/figures/fig5/data/fig5_tasks.jsonl"
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    seed = 42

    # ── Load & split (identical to CoordRep-Ranker) ───────
    print("Loading complexes (same split as CoordRep-Ranker) …")
    complexes = _extract_unique_complexes(fig5_path)
    train_cx, val_cx, test_cx = split_complexes(complexes, seed=seed)
    print(f"  train={len(train_cx)} val={len(val_cx)} test={len(test_cx)}")

    # Verify against saved split
    cr_split_path = "/data/CoordRep/coordrep-release/libcoordrep/checkpoints/coordrep_ranker/split_ids.json"
    with open(cr_split_path) as f:
        cr_split = json.load(f)
    assert set(c['id'] for c in test_cx) == set(cr_split['test']), "SPLIT MISMATCH!"
    print("  ✓ Test split matches CoordRep-Ranker checkpoint split_ids.json")

    # ── Build pool from train only ────────────────────────
    print("Building hard decoy pool (train only) …")
    t0 = time.time()
    train_toks = [c['tokens'] for c in train_cx]
    pool = HardDecoyPool(train_toks, seed=seed)
    print(f"  Done ({time.time()-t0:.1f}s)")

    # ── Generate decoys for ALL 1501 test complexes ───────
    print("Generating decoys for all test complexes (1 real + 20 decoys) …")
    n_decoys_target = 20
    test_data = []
    rng = stdlib_random.Random(seed)
    for c in test_cx:
        decoys = generate_all_hard_decoys(c['tokens'], pool,
                                          n_per_type=max(1, n_decoys_target // 5))
        if not decoys:
            continue
        if len(decoys) > n_decoys_target:
            rng.shuffle(decoys)
            decoys = decoys[:n_decoys_target]
        test_data.append((c['id'], c['tokens'], decoys))
    print(f"  {len(test_data)} test groups, target {n_decoys_target} decoys each")

    # Stats on actual group sizes
    sizes = [len(d) for _, _, d in test_data]
    print(f"  Decoys per group: min={min(sizes)} max={max(sizes)} "
          f"mean={np.mean(sizes):.1f} median={np.median(sizes):.0f}")

    # ── Train GIN Ranker (same as before but on full data) ─
    print("\nTraining GIN Ranker (3000 train, full 1501 test) …")
    from run_hard_negative_ranker import train_gin_ranker

    train_graph_data = []
    for c in train_cx[:3000]:
        decoys = generate_all_hard_decoys(c['tokens'], pool, n_per_type=2)
        if not decoys:
            continue
        real_g = tokens_to_simple_graph(c['tokens'])
        decoy_gs = [(tokens_to_simple_graph(dt), reason) for dt, reason in decoys]
        train_graph_data.append((real_g, [dg for dg, _ in decoy_gs],
                                 [r for _, r in decoy_gs]))

    val_graph_data = []
    for c in val_cx[:500]:
        decoys = generate_all_hard_decoys(c['tokens'], pool, n_per_type=2)
        if not decoys:
            continue
        real_g = tokens_to_simple_graph(c['tokens'])
        decoy_gs = [(tokens_to_simple_graph(dt), reason) for dt, reason in decoys]
        val_graph_data.append((real_g, [dg for dg, _ in decoy_gs],
                                [r for _, r in decoy_gs]))

    print(f"  Train graphs: {len(train_graph_data)}, Val: {len(val_graph_data)}")

    in_dim = train_graph_data[0][0].x.size(1) if train_graph_data else 11
    gin_ranker = GINRanker(in_dim=in_dim, hidden_dim=128, n_layers=4).to(device)
    gin_ranker = train_gin_ranker(
        gin_ranker, train_graph_data, val_graph_data, device,
        epochs=15, lr=1e-3, patience=5)

    # ── Evaluate all methods ──────────────────────────────
    print("\nEvaluating all methods on full test set …")

    # Non-neural baselines
    random_bl = RandomBaseline(seed)
    freq_bl = FrequencyBaseline()
    freq_bl.fit(train_toks)
    wl_bl = WLHashBaseline()
    wl_bl.fit(train_toks)

    methods = {
        'Random': random_bl.score,
        'Frequency': freq_bl.score,
        'WLHash': wl_bl.score,
    }

    # Results accumulator: method -> {top1, mrr, win_rate, real_scores, decoy_scores, by_type}
    all_metrics = {}
    for name in list(methods.keys()) + ['GIN_Ranker']:
        all_metrics[name] = {
            'top1_list': [],
            'mrr_list': [],
            'win_rate_list': [],
            'real_scores': [],
            'decoy_scores': [],
            'by_type_real': defaultdict(list),
            'by_type_decoy': defaultdict(list),
        }

    # Evaluate non-neural
    for name, score_fn in methods.items():
        for cid, real_toks, decoys in test_data:
            real_s = score_fn(real_toks)
            dec_scores = [(score_fn(dt), reason) for dt, reason in decoys]
            n_wins = sum(1 for ds, _ in dec_scores if real_s > ds)
            rank = 1 + sum(1 for ds, _ in dec_scores if ds >= real_s)

            all_metrics[name]['top1_list'].append(int(rank == 1))
            all_metrics[name]['mrr_list'].append(1.0 / rank)
            all_metrics[name]['win_rate_list'].append(n_wins / len(dec_scores))
            all_metrics[name]['real_scores'].append(real_s)
            for ds, reason in dec_scores:
                all_metrics[name]['decoy_scores'].append(ds)
                dtype = reason.split(':')[0]
                all_metrics[name]['by_type_real'][dtype].append(real_s)
                all_metrics[name]['by_type_decoy'][dtype].append(ds)

    # Evaluate GIN Ranker
    print("  Scoring GIN Ranker …")
    gin_ranker.eval()
    with torch.no_grad():
        for cid, real_toks, decoys in test_data:
            real_g = tokens_to_simple_graph(real_toks)
            all_graphs = [real_g] + [tokens_to_simple_graph(dt) for dt, _ in decoys]
            batch = Batch.from_data_list(all_graphs).to(device)
            scores = gin_ranker(batch.x, batch.edge_index, batch.batch).cpu().numpy()

            real_s = float(scores[0])
            dec_s = scores[1:]
            n_wins = sum(1 for ds in dec_s if real_s > ds)
            rank = 1 + sum(1 for ds in dec_s if ds >= real_s)

            all_metrics['GIN_Ranker']['top1_list'].append(int(rank == 1))
            all_metrics['GIN_Ranker']['mrr_list'].append(1.0 / rank)
            all_metrics['GIN_Ranker']['win_rate_list'].append(n_wins / len(dec_s))
            all_metrics['GIN_Ranker']['real_scores'].append(real_s)
            for i, (dt, reason) in enumerate(decoys):
                all_metrics['GIN_Ranker']['decoy_scores'].append(float(dec_s[i]))
                dtype = reason.split(':')[0]
                all_metrics['GIN_Ranker']['by_type_real'][dtype].append(real_s)
                all_metrics['GIN_Ranker']['by_type_decoy'][dtype].append(float(dec_s[i]))

    # ── Add CoordRep-Ranker (finetuned) from saved results ─
    # These are the authoritative numbers from the ranker evaluation
    all_metrics['CoordRep_Ranker_finetuned'] = {
        'source': 'revision_results/coordrep_ranker/ranker_summary.csv',
        'n_test': 1501,
        'n_decoys_per': 20,
    }
    cr_by_type_auroc = {
        'boundary_hard': 0.6979,
        'co_ligand_hard': 0.5527,
        'ligand_hard': 0.6539,
        'metal_hard': 0.7913,
        'stereo_hard': 0.9475,
    }

    # ── Compute summary ───────────────────────────────────
    print("\n" + "="*70)
    print("QC AUDIT: HARD-NEGATIVE BENCHMARK")
    print("="*70)

    print(f"\n1. Same test groups? YES — verified split_ids.json match")
    print(f"2. Candidates per group: 1 real + {n_decoys_target} decoys (actual mean={np.mean(sizes):.1f})")
    print(f"3. Top-1 = listwise (rank==1 among 1+K candidates)")
    print(f"4. README '78%' was WRONG — should be win_rate=74.3% or top1=13.1%")
    print(f"5. All 5 decoy types used (not just stereo_hard)")
    print(f"6. No train/test mixing — same split function, seed=42")

    # ── Answer 7: Output 20 groups ────────────────────────
    print(f"\n7. First 20 groups (GIN Ranker scores):")
    print(f"   {'CID':<12} {'RealScore':>10} {'BestDecoy':>10} {'Rank':>5} {'WinRate':>8}")
    gin_ranker.eval()
    sample_groups = []
    with torch.no_grad():
        for i, (cid, real_toks, decoys) in enumerate(test_data[:20]):
            real_g = tokens_to_simple_graph(real_toks)
            all_graphs = [real_g] + [tokens_to_simple_graph(dt) for dt, _ in decoys]
            batch = Batch.from_data_list(all_graphs).to(device)
            scores = gin_ranker(batch.x, batch.edge_index, batch.batch).cpu().numpy()
            real_s = float(scores[0])
            dec_s = [float(s) for s in scores[1:]]
            rank = 1 + sum(1 for ds in dec_s if ds >= real_s)
            wr = sum(1 for ds in dec_s if real_s > ds) / len(dec_s)
            print(f"   {cid:<12} {real_s:>10.4f} {max(dec_s):>10.4f} {rank:>5d} {wr:>8.3f}")
            sample_groups.append({
                'complex_id': cid,
                'n_decoys': len(dec_s),
                'real_score': round(real_s, 6),
                'decoy_scores': [round(ds, 6) for ds in dec_s],
                'real_rank': rank,
                'win_rate': round(wr, 4),
            })

    # Save sample groups
    p = os.path.join(out_dir, "hard_negative_20_groups_debug.json")
    with open(p, 'w') as f:
        json.dump(sample_groups, f, indent=2)
    print(f"   → Saved to {p}")

    # ── Answer 8: Generate hard_negative_summary_fixed.csv ─
    print(f"\n8. Generating hard_negative_summary_fixed.csv …")

    rows = []
    for name in ['Random', 'Frequency', 'WLHash', 'GIN_Ranker']:
        m = all_metrics[name]
        n = len(m['top1_list'])
        listwise_top1 = float(np.mean(m['top1_list']))
        mrr = float(np.mean(m['mrr_list']))
        pairwise_win_rate = float(np.mean(m['win_rate_list']))
        auroc = compute_auroc(m['real_scores'], m['decoy_scores'])

        # Per decoy type AUROC
        type_aurocs = {}
        for dtype in HARD_DECOY_TYPES:
            r_scores = m['by_type_real'].get(dtype, [])
            d_scores = m['by_type_decoy'].get(dtype, [])
            if r_scores and d_scores:
                type_aurocs[dtype] = compute_auroc(r_scores, d_scores)
            else:
                type_aurocs[dtype] = None

        row = {
            'method': name,
            'n_test_groups': n,
            'n_decoys_per_group': n_decoys_target,
            'listwise_top1': round(listwise_top1, 4),
            'mrr': round(mrr, 4),
            'pairwise_win_rate': round(pairwise_win_rate, 4),
            'auroc': round(auroc, 4),
        }
        for dtype in HARD_DECOY_TYPES:
            row[f'auroc_{dtype}'] = round(type_aurocs[dtype], 4) if type_aurocs[dtype] else ''
        rows.append(row)

    # Add CoordRep-Ranker from saved results
    rows.append({
        'method': 'CoordRep_Ranker_finetuned',
        'n_test_groups': 1501,
        'n_decoys_per_group': 20,
        'listwise_top1': 0.1306,
        'mrr': 0.3337,
        'pairwise_win_rate': 0.7431,
        'auroc': 0.6953,
        'auroc_metal_hard': 0.7913,
        'auroc_ligand_hard': 0.6539,
        'auroc_stereo_hard': 0.9475,
        'auroc_co_ligand_hard': 0.5527,
        'auroc_boundary_hard': 0.6979,
    })

    # Add MLM-PLL from saved results
    rows.append({
        'method': 'CoordRep_MLM_PLL',
        'n_test_groups': 1501,
        'n_decoys_per_group': 20,
        'listwise_top1': 0.0040,
        'mrr': 0.1139,
        'pairwise_win_rate': 0.5770,
        'auroc': 0.5535,
        'auroc_metal_hard': 0.6773,
        'auroc_ligand_hard': 0.6513,
        'auroc_stereo_hard': 0.7448,
        'auroc_co_ligand_hard': 0.4838,
        'auroc_boundary_hard': 0.6079,
    })

    # Write CSV
    fieldnames = ['method', 'n_test_groups', 'n_decoys_per_group',
                  'listwise_top1', 'mrr', 'pairwise_win_rate', 'auroc']
    for dtype in HARD_DECOY_TYPES:
        fieldnames.append(f'auroc_{dtype}')

    p = os.path.join(out_dir, "hard_negative_summary_fixed.csv")
    with open(p, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for row in rows:
            w.writerow(row)
    print(f"   → {p}")

    # Also JSON
    p = os.path.join(out_dir, "hard_negative_summary_fixed.json")
    with open(p, 'w') as f:
        json.dump(rows, f, indent=2)
    print(f"   → {p}")

    # ── Final summary table ───────────────────────────────
    print(f"\n{'='*70}")
    print("CORRECTED HARD-NEGATIVE SUMMARY")
    print(f"{'='*70}")
    print(f"  {'Method':<25} {'N':>5} {'Top1':>7} {'MRR':>7} {'WinRate':>8} {'AUROC':>7}")
    print(f"  {'-'*62}")
    for row in rows:
        print(f"  {row['method']:<25} {row['n_test_groups']:>5} "
              f"{row['listwise_top1']:>7.4f} {row['mrr']:>7.4f} "
              f"{row['pairwise_win_rate']:>8.4f} {row['auroc']:>7.4f}")

    print(f"\n  Per-decoy-type AUROC:")
    print(f"  {'Method':<25} {'metal':>7} {'ligand':>7} {'stereo':>7} {'co_lig':>7} {'bound':>7}")
    print(f"  {'-'*62}")
    for row in rows:
        vals = []
        for dtype in HARD_DECOY_TYPES:
            v = row.get(f'auroc_{dtype}', '')
            vals.append(f"{v:>7.4f}" if isinstance(v, (int, float)) else f"{'':>7}")
        print(f"  {row['method']:<25} {'  '.join(vals)}")


if __name__ == "__main__":
    main()
