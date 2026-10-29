#!/usr/bin/env python3
"""
run_hard_negative_ranker.py
===========================
Part 2: Hard-negative compatibility GNN ranker baseline.

Compares:
- Rule validator baseline (should be ~random)
- Random baseline
- Frequency baseline
- WL graph hash baseline
- GIN graph ranker (2D)
- CoordRep-Ranker (from existing results, loaded if checkpoint exists)
- CoordRep-MLM field scorer

Uses the same hard decoy generation as the existing eval pipeline.

Outputs: revision_results/gnn_baselines/hard_negative_*.csv
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
from typing import List, Tuple

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.data import Data, Batch
from torch_geometric.loader import DataLoader

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from coordrep_tools.hard_decoys import (
    HardDecoyPool, generate_all_hard_decoys, HARD_DECOY_TYPES,
    _parse_metal_token,
)
from coordrep_tools.ablation_masking import _detect_blocks
from coordrep_tools.ranker_dataset import (
    _extract_unique_complexes, split_complexes, _tokens_to_ids,
)
from coordrep_tools.gnn_baselines.models import GINRanker
from coordrep_tools.gnn_baselines.non_neural import wl_hash, canonical_smiles_multiset_key

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError:
    Chem = None


# ── Utility: tokens → graph ───────────────────────────────

ATOM_CHARS = set("BCNOPSFIHcnospbfli")
BOND_CHARS = set("=#\\/-.")

def tokens_to_simple_graph(tokens: List[str]) -> Data:
    """
    Convert CoordRep token list to a simplified graph for GNN ranker.
    Nodes = token embeddings (one-hot over token type categories).
    Edges = sequential (token i ↔ token i+1) + structural.
    """
    # Categorize each token
    categories = {
        'metal': 0, 'shape': 1, 'constraint': 2,
        'lig_id': 3, 'donor_marker': 4, 'smiles_atom': 5,
        'smiles_bond': 6, 'delimiter': 7, 'special': 8, 'other': 9,
    }
    n_cat = len(categories) + 1

    blocks = _detect_blocks(tokens)
    n = len(tokens)
    if n == 0:
        return Data(
            x=torch.zeros((1, n_cat), dtype=torch.float),
            edge_index=torch.zeros((2, 0), dtype=torch.long),
        )

    # Node features
    x = torch.zeros((n, n_cat), dtype=torch.float)
    for i, tok in enumerate(tokens):
        if tok.startswith('[Metal:') or tok.startswith('[') and not tok.startswith('[MASK'):
            x[i, categories['metal']] = 1
        elif tok.startswith('<Shape:') or tok.startswith('V_'):
            x[i, categories['shape']] = 1
        elif tok.startswith('{') or tok == '}':
            x[i, categories['constraint']] = 1
        elif re.match(r'^L\d+$', tok):
            x[i, categories['lig_id']] = 1
        elif re.match(r'^:[A-Z][a-z]?:\d+$', tok):
            x[i, categories['donor_marker']] = 1
        elif len(tok) == 1 and tok in ATOM_CHARS:
            x[i, categories['smiles_atom']] = 1
        elif len(tok) == 1 and tok in BOND_CHARS:
            x[i, categories['smiles_bond']] = 1
        elif tok in ('|', '=', '--', ','):
            x[i, categories['delimiter']] = 1
        elif tok in ('[CLS]', '[SEP]', '[PAD]', '[MASK]', '[UNK]'):
            x[i, categories['special']] = 1
        else:
            x[i, categories['other']] = 1

    # Sequential edges
    edge_index = []
    for i in range(n - 1):
        edge_index.append([i, i + 1])
        edge_index.append([i + 1, i])

    if edge_index:
        edge_index = torch.tensor(edge_index, dtype=torch.long).t().contiguous()
    else:
        edge_index = torch.zeros((2, 0), dtype=torch.long)

    return Data(x=x, edge_index=edge_index, num_nodes=n)


# ── WL Hash Baseline ─────────────────────────────────────

class WLHashBaseline:
    """Score by counting WL-hash matches between ligands and training set."""

    def __init__(self):
        self.train_hashes = Counter()

    def fit(self, train_token_lists):
        for token_list in train_token_lists:
            blocks = _detect_blocks(token_list)
            for lb in blocks.get('ligand_blocks', []):
                s, e = lb.get('smiles_start'), lb.get('smiles_end')
                if s is not None and e is not None:
                    smi = ''.join(token_list[s:e])
                    h = wl_hash(smi)
                    self.train_hashes[h] += 1

    def score(self, tokens):
        blocks = _detect_blocks(tokens)
        total = 0
        for lb in blocks.get('ligand_blocks', []):
            s, e = lb.get('smiles_start'), lb.get('smiles_end')
            if s is not None and e is not None:
                smi = ''.join(tokens[s:e])
                h = wl_hash(smi)
                total += math.log(self.train_hashes.get(h, 1))
        return total


class RuleValidatorBaseline:
    """Score = 1 if valid, 0 otherwise. Should be ~random on rule-passing decoys."""

    def score(self, tokens):
        from coordrep_tools.validate import is_valid_coordrep
        s = ''.join(tokens)
        return 1.0 if is_valid_coordrep(s, strict=True) else 0.0


class RandomBaseline:
    def __init__(self, seed=42):
        self.rng = stdlib_random.Random(seed)

    def score(self, tokens):
        return self.rng.random()


class FrequencyBaseline:
    """Score = sum of log-frequencies of tokens."""

    def __init__(self):
        self.token_freq = Counter()

    def fit(self, train_token_lists):
        for toks in train_token_lists:
            self.token_freq.update(toks)

    def score(self, tokens):
        return sum(math.log(self.token_freq.get(t, 1)) for t in tokens)


# ── GIN Ranker training ──────────────────────────────────

def train_gin_ranker(model, train_data, val_data, device,
                     epochs=20, lr=1e-3, patience=5, batch_size=16):
    """
    Train GIN ranker with listwise softmax loss.
    train_data: list of (real_graph, [decoy_graphs], decoy_types)
    """
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    best_val_loss = float('inf')
    best_state = None
    patience_counter = 0

    for epoch in range(epochs):
        model.train()
        total_loss = 0
        n_batches = 0
        stdlib_random.shuffle(train_data)

        for i in range(0, len(train_data), batch_size):
            batch_items = train_data[i:i + batch_size]

            # Build batch: score all graphs, label 0 = real
            all_graphs = []
            labels = []
            group_sizes = []
            for real_g, decoy_gs, _ in batch_items:
                all_graphs.append(real_g)
                for dg in decoy_gs:
                    all_graphs.append(dg)
                group_sizes.append(1 + len(decoy_gs))

            if not all_graphs:
                continue

            batch = Batch.from_data_list(all_graphs).to(device)
            scores = model(batch.x, batch.edge_index, batch.batch)

            # Listwise softmax loss: real should have highest score
            loss = 0
            offset = 0
            for gs in group_sizes:
                group_scores = scores[offset:offset + gs]
                # Label: first one is real
                target = torch.zeros(gs, device=device)
                target[0] = 1.0
                loss += F.cross_entropy(group_scores.unsqueeze(0), torch.tensor([0], device=device))
                offset += gs

            loss = loss / len(group_sizes)
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            total_loss += loss.item()
            n_batches += 1

        # Validation
        model.eval()
        val_loss = 0
        n_val = 0
        with torch.no_grad():
            for i in range(0, len(val_data), batch_size):
                batch_items = val_data[i:i + batch_size]
                all_graphs = []
                group_sizes = []
                for real_g, decoy_gs, _ in batch_items:
                    all_graphs.append(real_g)
                    for dg in decoy_gs:
                        all_graphs.append(dg)
                    group_sizes.append(1 + len(decoy_gs))

                if not all_graphs:
                    continue
                batch = Batch.from_data_list(all_graphs).to(device)
                scores = model(batch.x, batch.edge_index, batch.batch)

                offset = 0
                for gs in group_sizes:
                    group_scores = scores[offset:offset + gs]
                    val_loss += F.cross_entropy(
                        group_scores.unsqueeze(0),
                        torch.tensor([0], device=device)).item()
                    n_val += 1
                    offset += gs

        val_loss /= max(n_val, 1)
        if val_loss < best_val_loss:
            best_val_loss = val_loss
            best_state = {k: v.clone() for k, v in model.state_dict().items()}
            patience_counter = 0
        else:
            patience_counter += 1
            if patience_counter >= patience:
                break

        if (epoch + 1) % 5 == 0:
            print(f"    Epoch {epoch+1}: train={total_loss/max(n_batches,1):.4f} val={val_loss:.4f}")

    if best_state:
        model.load_state_dict(best_state)
    return model


# ── Evaluation ────────────────────────────────────────────

def evaluate_ranker(score_fn, test_data, method_name):
    """
    Evaluate a scoring function on test hard-negative data.
    score_fn(tokens) → float
    test_data: list of (complex_id, real_tokens, [(decoy_tokens, reason), ...])
    """
    results = []
    by_type = defaultdict(lambda: {'n': 0, 'wins': 0})

    for cid, real_toks, decoys in test_data:
        if not decoys:
            continue
        real_score = score_fn(real_toks)
        decoy_scores = [(score_fn(dtoks), reason) for dtoks, reason in decoys]

        n_wins = sum(1 for ds, _ in decoy_scores if real_score > ds)
        rank = 1 + sum(1 for ds, _ in decoy_scores if ds >= real_score)
        mrr = 1.0 / rank
        top1 = int(rank == 1)
        win_rate = n_wins / len(decoy_scores) if decoy_scores else 0

        results.append({
            'complex_id': cid,
            'method': method_name,
            'n_decoys': len(decoys),
            'real_score': real_score,
            'top1': top1,
            'mrr': mrr,
            'win_rate': win_rate,
            'rank': rank,
        })

        for ds, reason in decoy_scores:
            dtype = reason.split(':')[0]
            by_type[dtype]['n'] += 1
            if real_score > ds:
                by_type[dtype]['wins'] += 1

    return results, dict(by_type)


def evaluate_gin_ranker(model, test_items, device, method_name):
    """Evaluate GIN ranker on graph-level scoring."""
    model.eval()
    results = []
    by_type = defaultdict(lambda: {'n': 0, 'wins': 0})

    with torch.no_grad():
        for cid, real_g, decoy_items in test_items:
            all_graphs = [real_g] + [dg for dg, _ in decoy_items]
            batch = Batch.from_data_list(all_graphs).to(device)
            scores = model(batch.x, batch.edge_index, batch.batch).cpu().numpy()

            real_score = scores[0]
            decoy_scores = scores[1:]

            n_wins = sum(1 for ds in decoy_scores if real_score > ds)
            rank = 1 + sum(1 for ds in decoy_scores if ds >= real_score)
            mrr = 1.0 / rank
            top1 = int(rank == 1)

            results.append({
                'complex_id': cid,
                'method': method_name,
                'n_decoys': len(decoy_scores),
                'top1': top1,
                'mrr': mrr,
                'win_rate': n_wins / max(len(decoy_scores), 1),
                'rank': rank,
            })

            for i, (_, reason) in enumerate(decoy_items):
                dtype = reason.split(':')[0]
                by_type[dtype]['n'] += 1
                if real_score > decoy_scores[i]:
                    by_type[dtype]['wins'] += 1

    return results, dict(by_type)


# ── Main ──────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fig5",
                        default="inputs/property_benchmarks/fig5_tasks.jsonl")
    parser.add_argument("--max_test", type=int, default=1000)
    parser.add_argument("--decoys_per_real", type=int, default=20)
    parser.add_argument("--gin_epochs", type=int, default=15)
    parser.add_argument("--device", default="cuda:0")
    parser.add_argument("--out", default="revision_results/gnn_baselines")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    os.makedirs(args.out, exist_ok=True)
    device = torch.device(args.device if torch.cuda.is_available() else 'cpu')

    # ── Load complexes ────────────────────────────────────
    print("Loading unique complexes …")
    complexes = _extract_unique_complexes(args.fig5)
    train_cx, val_cx, test_cx = split_complexes(complexes, seed=args.seed)
    print(f"  train={len(train_cx)} val={len(val_cx)} test={len(test_cx)}")

    test_cx = test_cx[:args.max_test]

    # ── Build hard decoy pool ─────────────────────────────
    print("Building hard decoy pool (train only) …")
    t0 = time.time()
    train_toks = [c['tokens'] for c in train_cx]
    pool = HardDecoyPool(train_toks, seed=args.seed)
    print(f"  Done ({time.time()-t0:.1f}s)")

    # ── Generate test decoys ──────────────────────────────
    print("Generating test decoys …")
    test_data = []
    for c in test_cx:
        decoys = generate_all_hard_decoys(
            c['tokens'], pool,
            n_per_type=max(1, args.decoys_per_real // 5))
        if decoys:
            if len(decoys) > args.decoys_per_real:
                stdlib_random.Random(args.seed).shuffle(decoys)
                decoys = decoys[:args.decoys_per_real]
            test_data.append((c['id'], c['tokens'], decoys))
    print(f"  {len(test_data)} test complexes with decoys")

    # ── Non-neural baselines ──────────────────────────────
    print(f"\n{'='*60}")
    print("Non-neural baselines")
    print(f"{'='*60}")

    random_bl = RandomBaseline(args.seed)
    rule_bl = RuleValidatorBaseline()
    freq_bl = FrequencyBaseline()
    freq_bl.fit(train_toks)
    wl_bl = WLHashBaseline()
    wl_bl.fit(train_toks)

    all_results = {}
    all_by_type = {}

    for name, bl in [('Random', random_bl), ('RuleValidator', rule_bl),
                      ('Frequency', freq_bl), ('WLHash', wl_bl)]:
        res, bt = evaluate_ranker(bl.score, test_data, name)
        all_results[name] = res
        all_by_type[name] = bt
        if res:
            avg_top1 = np.mean([r['top1'] for r in res])
            avg_mrr = np.mean([r['mrr'] for r in res])
            print(f"  {name}: Top-1={avg_top1:.3f} MRR={avg_mrr:.3f}")

    # ── GIN Graph Ranker ──────────────────────────────────
    print(f"\n{'='*60}")
    print("Training GIN Graph Ranker …")
    print(f"{'='*60}")

    # Build graph training data
    train_graph_data = []
    for c in train_cx[:3000]:  # limit for speed
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

    test_graph_items = []
    for cid, toks, decoys in test_data:
        real_g = tokens_to_simple_graph(toks)
        decoy_items = [(tokens_to_simple_graph(dt), reason) for dt, reason in decoys]
        test_graph_items.append((cid, real_g, decoy_items))

    print(f"  Train graphs: {len(train_graph_data)}, Val: {len(val_graph_data)}")

    if train_graph_data:
        in_dim = train_graph_data[0][0].x.size(1)
        gin_ranker = GINRanker(in_dim=in_dim, hidden_dim=128, n_layers=4).to(device)

        gin_ranker = train_gin_ranker(
            gin_ranker, train_graph_data, val_graph_data, device,
            epochs=args.gin_epochs, lr=1e-3, patience=5)

        gin_res, gin_bt = evaluate_gin_ranker(gin_ranker, test_graph_items, device, 'GIN_Ranker')
        all_results['GIN_Ranker'] = gin_res
        all_by_type['GIN_Ranker'] = gin_bt

        if gin_res:
            avg_top1 = np.mean([r['top1'] for r in gin_res])
            avg_mrr = np.mean([r['mrr'] for r in gin_res])
            print(f"  GIN_Ranker: Top-1={avg_top1:.3f} MRR={avg_mrr:.3f}")

        del gin_ranker
        torch.cuda.empty_cache()

    # ── Write outputs ─────────────────────────────────────
    print(f"\n{'='*60}")
    print("Writing outputs")
    print(f"{'='*60}")

    # Summary CSV
    p = os.path.join(args.out, "hard_negative_summary.csv")
    with open(p, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["method", "n", "top1", "mrr", "win_rate"])
        for name, res in all_results.items():
            if res:
                w.writerow([
                    name, len(res),
                    f"{np.mean([r['top1'] for r in res]):.4f}",
                    f"{np.mean([r['mrr'] for r in res]):.4f}",
                    f"{np.mean([r['win_rate'] for r in res]):.4f}",
                ])
    print(f"  {p}")

    # By decoy type
    p = os.path.join(args.out, "hard_negative_by_decoy_type.csv")
    with open(p, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["method", "decoy_type", "n", "win_rate"])
        for name, bt in all_by_type.items():
            for dtype, vals in sorted(bt.items()):
                wr = vals['wins'] / max(vals['n'], 1)
                w.writerow([name, dtype, vals['n'], f"{wr:.4f}"])
    print(f"  {p}")

    # Summary JSON
    summary = {}
    for name, res in all_results.items():
        if res:
            summary[name] = {
                'n': len(res),
                'top1': float(np.mean([r['top1'] for r in res])),
                'mrr': float(np.mean([r['mrr'] for r in res])),
                'win_rate': float(np.mean([r['win_rate'] for r in res])),
            }
    p = os.path.join(args.out, "hard_negative_summary.json")
    with open(p, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"  {p}")

    # Final table
    print(f"\n{'='*60}")
    print("HARD NEGATIVE RANKER SUMMARY")
    print(f"{'='*60}")
    print(f"  {'Method':<20s} {'N':>6s} {'Top-1':>7s} {'MRR':>7s} {'WinRate':>8s}")
    print(f"  {'-'*50}")
    for name, res in all_results.items():
        if res:
            print(f"  {name:<20s} {len(res):>6d} "
                  f"{np.mean([r['top1'] for r in res]):>7.3f} "
                  f"{np.mean([r['mrr'] for r in res]):>7.3f} "
                  f"{np.mean([r['win_rate'] for r in res]):>8.3f}")


if __name__ == "__main__":
    main()
