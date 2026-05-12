#!/usr/bin/env python3
"""
run_donor_annotation.py
=======================
Part 1: Donor annotation GNN baseline comparison.

Compares:
- Non-neural: Random, GlobalFreq, CondFreq(metal,CN), LigandFreq(SMILES)
- 2D GNN: GCN, GIN (ligand_only and ligand_plus_context)
- 3D GNN: EGNN upper bound (geometry-aware)
- CoordRep-MLM reference (from existing results)

Outputs to revision_results/gnn_baselines/donor_annotation_*.csv
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.data import Data, Batch
from torch_geometric.loader import DataLoader

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from coordrep_tools.gnn_baselines.datasets import (
    load_pipeline_data, build_unified_split, leakage_diagnostic,
    LigandGraphDataset, mol_to_graph, smiles_to_graph_with_donors,
    atom_features, parse_tmqm_xyz, build_3d_graph, TMQM_XYZ_DIR,
)
from coordrep_tools.gnn_baselines.models import (
    GCNDonorPredictor, GINDonorPredictor, GATDonorPredictor, EGNNDonorPredictor,
)
from coordrep_tools.gnn_baselines.non_neural import (
    RandomDonorBaseline, GlobalFreqDonorBaseline,
    CondFreqDonorBaseline, LigandFreqDonorBaseline,
)


# ── Metrics ───────────────────────────────────────────────

def donor_metrics(pred_set: set, true_set: set, n_atoms: int):
    """Compute donor annotation metrics for one ligand."""
    tp = len(pred_set & true_set)
    fp = len(pred_set - true_set)
    fn = len(true_set - pred_set)

    prec = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    rec = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1 = 2 * prec * rec / (prec + rec) if (prec + rec) > 0 else 0.0

    # Jaccard
    union = len(pred_set | true_set)
    jaccard = len(pred_set & true_set) / union if union > 0 else 0.0

    all_correct = pred_set == true_set
    return {
        'precision': prec, 'recall': rec, 'f1': f1,
        'jaccard': jaccard, 'all_correct': int(all_correct),
    }


def element_topk(pred_elements: list, true_elements: list, k: int = 5):
    """Element-level top-k accuracy."""
    true_set = set(true_elements)
    top1 = int(pred_elements[0] in true_set) if pred_elements else 0
    topk = int(bool(set(pred_elements[:k]) & true_set)) if pred_elements else 0
    return top1, topk


# ── GNN Training & Evaluation ────────────────────────────

def build_pyg_dataset(samples, include_context=False):
    """Build PyG Data list from LigandGraphDataset samples."""
    data_list = []
    meta_list = []
    for s in samples:
        data = smiles_to_graph_with_donors(
            s['smiles'], s['donor_elements'], s['denticity'])
        if data is None:
            continue
        if include_context:
            # Encode metal + CN as context vector
            METALS = ['Fe', 'Cu', 'Ni', 'Co', 'Zn', 'Mn', 'Cr', 'V',
                      'Ti', 'Pd', 'Pt', 'Ru', 'Rh', 'Ir', 'Os', 'Ag',
                      'Au', 'Cd', 'Hg', 'Sc', 'Mo', 'W', 'Re', 'Ta',
                      'Nb', 'Zr', 'Hf', 'Y', 'La', 'other']
            metal = s.get('metal', 'other')
            metal_oh = [0] * (len(METALS) + 1)
            try:
                metal_oh[METALS.index(metal)] = 1
            except ValueError:
                metal_oh[-1] = 1
            cn_feat = [s.get('cn', 0) / 8.0]
            dent_feat = [s.get('denticity', 1) / 6.0]
            ctx = torch.tensor(metal_oh + cn_feat + dent_feat, dtype=torch.float)
            data.context = ctx

        data_list.append(data)
        meta_list.append(s)
    return data_list, meta_list


def train_gnn_donor(model, train_loader, val_loader, device,
                    epochs=30, lr=1e-3, patience=5):
    """Train a GNN donor predictor with early stopping."""
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    best_val_loss = float('inf')
    best_state = None
    patience_counter = 0

    for epoch in range(epochs):
        model.train()
        total_loss = 0
        n_batches = 0
        for batch in train_loader:
            batch = batch.to(device)
            # Context: per-graph vector → per-node via batch index
            ctx = None
            if hasattr(batch, 'context') and batch.context is not None:
                ctx_raw = batch.context
                if ctx_raw.dim() == 1:
                    # Infer context dim from model
                    n_graphs = batch.batch.max().item() + 1
                    ctx_dim = ctx_raw.numel() // n_graphs
                    ctx_raw = ctx_raw.view(n_graphs, ctx_dim)
                ctx = ctx_raw  # shape (n_graphs, ctx_dim)

            if hasattr(model, 'forward') and 'pos' in model.forward.__code__.co_varnames:
                logits = model(batch.x, batch.edge_index, batch.pos, batch.batch, ctx)
            else:
                logits = model(batch.x, batch.edge_index, batch.batch, ctx)

            loss = F.binary_cross_entropy_with_logits(logits, batch.y)
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
            for batch in val_loader:
                batch = batch.to(device)
                ctx = None
                if hasattr(batch, 'context') and batch.context is not None:
                    ctx_raw = batch.context
                    if ctx_raw.dim() == 1:
                        n_graphs = batch.batch.max().item() + 1
                        ctx_dim = ctx_raw.numel() // n_graphs
                        ctx_raw = ctx_raw.view(n_graphs, ctx_dim)
                    ctx = ctx_raw
                if hasattr(model, 'forward') and 'pos' in model.forward.__code__.co_varnames:
                    logits = model(batch.x, batch.edge_index, batch.pos, batch.batch, ctx)
                else:
                    logits = model(batch.x, batch.edge_index, batch.batch, ctx)
                loss = F.binary_cross_entropy_with_logits(logits, batch.y)
                val_loss += loss.item()
                n_val += 1

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
            print(f"      Epoch {epoch+1}: train={total_loss/max(n_batches,1):.4f} val={val_loss:.4f}")

    if best_state:
        model.load_state_dict(best_state)
    return model


def evaluate_gnn_donor(model, test_loader, test_meta, device, is_3d=False):
    """Evaluate GNN donor predictor on test set."""
    model.eval()
    all_results = []

    offset = 0
    with torch.no_grad():
        for batch in test_loader:
            batch = batch.to(device)
            ctx = None
            if hasattr(batch, 'context') and batch.context is not None:
                ctx_raw = batch.context
                if ctx_raw.dim() == 1:
                    n_graphs = batch.batch.max().item() + 1
                    ctx_dim = ctx_raw.numel() // n_graphs
                    ctx_raw = ctx_raw.view(n_graphs, ctx_dim)
                ctx = ctx_raw

            if is_3d:
                logits = model(batch.x, batch.edge_index, batch.pos, batch.batch, ctx)
            else:
                logits = model(batch.x, batch.edge_index, batch.batch, ctx)

            probs = torch.sigmoid(logits)

            # Per-graph evaluation
            ptr = batch.ptr if hasattr(batch, 'ptr') else None
            if ptr is not None:
                for gi in range(len(ptr) - 1):
                    start, end = ptr[gi].item(), ptr[gi + 1].item()
                    graph_probs = probs[start:end].cpu().numpy()
                    graph_labels = batch.y[start:end].cpu().numpy()

                    # Predict top-k donors
                    dent = int(graph_labels.sum())
                    if dent == 0:
                        dent = 1
                    pred_indices = set(np.argsort(-graph_probs)[:dent].tolist())
                    true_indices = set(np.where(graph_labels > 0.5)[0].tolist())

                    meta = test_meta[offset + gi] if (offset + gi) < len(test_meta) else {}
                    m = donor_metrics(pred_indices, true_indices, len(graph_probs))
                    m.update({
                        'complex_id': meta.get('complex_id', ''),
                        'denticity': meta.get('denticity', 0),
                        'metal': meta.get('metal', '?'),
                        'cn': meta.get('cn', 0),
                        'smiles': meta.get('smiles', ''),
                        'source': meta.get('source', ''),
                    })
                    all_results.append(m)
                offset += len(ptr) - 1

    return all_results


# ── Non-neural evaluation ─────────────────────────────────

def evaluate_non_neural_baselines(train_samples, test_samples):
    """Evaluate all non-neural baselines on test samples."""
    results = {}

    # Random
    rnd = RandomDonorBaseline()
    # GlobalFreq
    gfreq = GlobalFreqDonorBaseline().fit(train_samples)
    # CondFreq
    cfreq = CondFreqDonorBaseline().fit(train_samples)
    # LigandFreq
    lfreq = LigandFreqDonorBaseline().fit(train_samples)

    for name, baseline in [('Random', rnd), ('GlobalFreq', gfreq),
                            ('CondFreq', cfreq), ('LigandFreq', lfreq)]:
        preds = []
        for s in test_samples:
            true_donors = s.get('donor_elements', [])
            dent = s.get('denticity', 1)
            metal = s.get('metal', '?')
            cn = s.get('cn', 0)
            smi = s.get('smiles', '')

            if name == 'Random':
                # Element prediction: random from known elements
                pred_elems = list(np.random.choice(
                    ['C', 'N', 'O', 'S', 'P'], size=min(dent, 5), replace=True))
            elif name == 'GlobalFreq':
                pred_elems = baseline.predict_elements(dent)
            elif name == 'CondFreq':
                pred_elems = baseline.predict_elements(metal, cn, dent)
            elif name == 'LigandFreq':
                pred_elems = baseline.predict_elements(smi, dent)

            # Element-level metrics
            true_set = set(true_donors)
            pred_set = set(pred_elems[:len(true_donors)])

            m = donor_metrics(pred_set, true_set, 0)
            top1, top5 = element_topk(pred_elems, true_donors)
            m['top1'] = top1
            m['top5'] = top5
            m['denticity'] = dent
            m['metal'] = metal
            m['cn'] = cn
            m['smiles'] = smi
            m['source'] = s.get('source', '')
            m['is_rare'] = not lfreq.has_smiles(smi)
            preds.append(m)

        results[name] = preds

    return results


# ── Aggregate ─────────────────────────────────────────────

def aggregate_results(results_list, model_name):
    """Aggregate per-sample results into summary metrics."""
    n = len(results_list)
    if n == 0:
        return {}

    summary = {
        'model': model_name,
        'n': n,
        'precision': np.mean([r['precision'] for r in results_list]),
        'recall': np.mean([r['recall'] for r in results_list]),
        'f1': np.mean([r['f1'] for r in results_list]),
        'jaccard': np.mean([r['jaccard'] for r in results_list]),
        'all_correct': np.mean([r['all_correct'] for r in results_list]),
    }
    if 'top1' in results_list[0]:
        summary['top1'] = np.mean([r['top1'] for r in results_list])
        summary['top5'] = np.mean([r.get('top5', 0) for r in results_list])

    return summary


def aggregate_by_field(results_list, field, model_name):
    """Aggregate by a categorical field (denticity, source, etc.)."""
    groups = defaultdict(list)
    for r in results_list:
        val = r.get(field, '?')
        if field == 'denticity':
            if val >= 3:
                val = '>=3'
            else:
                val = str(val)
        elif field == 'is_rare':
            val = 'rare' if val else 'common'
        else:
            val = str(val)
        groups[val].append(r)

    rows = []
    for gval, gresults in sorted(groups.items()):
        row = aggregate_results(gresults, model_name)
        row[field] = gval
        rows.append(row)
    return rows


# ── Main ──────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--max_data", type=int, default=20000)
    parser.add_argument("--max_test", type=int, default=3000)
    parser.add_argument("--epochs", type=int, default=30)
    parser.add_argument("--batch_size", type=int, default=64)
    parser.add_argument("--device", type=str, default="cuda:0")
    parser.add_argument("--out", type=str,
                        default="revision_results/gnn_baselines")
    parser.add_argument("--skip_3d", action="store_true",
                        help="Skip 3D EGNN (slow XYZ parsing)")
    args = parser.parse_args()

    os.makedirs(args.out, exist_ok=True)
    device = torch.device(args.device if torch.cuda.is_available() else 'cpu')

    # ── Load data ─────────────────────────────────────────
    print("Loading pipeline data …")
    entries = load_pipeline_data(args.max_data)
    splits = build_unified_split(entries)
    diag = leakage_diagnostic(splits)
    print(f"  train={diag['train_size']} val={diag['val_size']} "
          f"test={diag['test_size']}")
    print(f"  mol_id leakage: train-test={diag['train_test_mol_id_leak']}")
    print(f"  SMILES overlap train-test: {diag['ligand_smiles_train_test_overlap']}"
          f"/{diag['ligand_smiles_test_total']}")

    with open(os.path.join(args.out, "leakage_diagnostics.json"), "w") as f:
        json.dump(diag, f, indent=2)

    # Build ligand datasets
    print("Building ligand graph datasets …")
    train_ds = LigandGraphDataset(splits['train'])
    val_ds = LigandGraphDataset(splits['val'])
    test_ds = LigandGraphDataset(splits['test'])
    print(f"  Ligand samples: train={len(train_ds)} val={len(val_ds)} test={len(test_ds)}")

    test_samples = test_ds.samples[:args.max_test]

    # ── Part A: Non-neural baselines ──────────────────────
    print(f"\n{'='*60}")
    print("Non-neural baselines")
    print(f"{'='*60}")
    nn_results = evaluate_non_neural_baselines(train_ds.samples, test_samples)

    all_summaries = []
    all_by_dent = []
    all_by_freq = []

    for name, preds in nn_results.items():
        s = aggregate_results(preds, name)
        all_summaries.append(s)
        all_by_dent.extend(aggregate_by_field(preds, 'denticity', name))
        all_by_freq.extend(aggregate_by_field(preds, 'is_rare', name))
        print(f"  {name}: F1={s['f1']:.3f} Jaccard={s['jaccard']:.3f} "
              f"AllCorrect={s['all_correct']:.3f}")

    # ── Part B: GNN baselines ─────────────────────────────
    gnn_configs = [
        ('GCN_ligand_only', GCNDonorPredictor, False),
        ('GIN_ligand_only', GINDonorPredictor, False),
        ('GIN_ligand_context', GINDonorPredictor, True),
    ]

    for gnn_name, ModelClass, use_context in gnn_configs:
        print(f"\n{'='*60}")
        print(f"Training {gnn_name} …")
        print(f"{'='*60}")

        train_pyg, train_meta = build_pyg_dataset(
            train_ds.samples, include_context=use_context)
        val_pyg, val_meta = build_pyg_dataset(
            val_ds.samples, include_context=use_context)
        test_pyg, test_meta = build_pyg_dataset(
            test_samples, include_context=use_context)

        if not train_pyg or not test_pyg:
            print(f"  SKIP: no valid graphs")
            continue

        in_dim = train_pyg[0].x.size(1)
        ctx_dim = train_pyg[0].context.size(0) if use_context and hasattr(train_pyg[0], 'context') else 0

        model = ModelClass(in_dim=in_dim, hidden_dim=128, n_layers=4,
                           context_dim=ctx_dim).to(device)

        train_loader = DataLoader(train_pyg, batch_size=args.batch_size, shuffle=True)
        val_loader = DataLoader(val_pyg, batch_size=args.batch_size)
        test_loader = DataLoader(test_pyg, batch_size=args.batch_size)

        t0 = time.time()
        model = train_gnn_donor(model, train_loader, val_loader, device,
                                epochs=args.epochs, lr=1e-3, patience=7)
        dt = time.time() - t0
        print(f"  Training done ({dt:.0f}s)")

        gnn_results = evaluate_gnn_donor(model, test_loader, test_meta, device)

        if gnn_results:
            s = aggregate_results(gnn_results, gnn_name)
            all_summaries.append(s)
            all_by_dent.extend(aggregate_by_field(gnn_results, 'denticity', gnn_name))
            print(f"  {gnn_name}: F1={s['f1']:.3f} Jaccard={s['jaccard']:.3f} "
                  f"AllCorrect={s['all_correct']:.3f}")

        del model
        torch.cuda.empty_cache()

    # ── Part C: 3D EGNN upper bound ───────────────────────
    if not args.skip_3d:
        print(f"\n{'='*60}")
        print("3D EGNN upper bound (parsing tmQM XYZ …)")
        print(f"{'='*60}")

        # Get mol_ids from test set
        test_mol_ids = {s.get('complex_id', s.get('mol_id', ''))
                        for s in test_samples}
        # Also need train mol_ids
        train_mol_ids = {s.get('complex_id', s.get('mol_id', ''))
                         for s in train_ds.samples}
        all_mol_ids = test_mol_ids | train_mol_ids

        # Parse XYZ files
        xyz_data = {}
        for xyz_file in ['tmQM_X1.xyz', 'tmQM_X2.xyz']:
            xyz_path = os.path.join(TMQM_XYZ_DIR, xyz_file)
            if os.path.exists(xyz_path):
                print(f"  Parsing {xyz_file} …")
                parsed = parse_tmqm_xyz(xyz_path, all_mol_ids)
                xyz_data.update(parsed)
                print(f"    Found {len(parsed)} molecules")

        if xyz_data:
            # Build 3D graphs for train + test
            train_3d = []
            for s in train_ds.samples[:5000]:  # limit for speed
                mid = s.get('complex_id', '')
                if mid in xyz_data:
                    g = build_3d_graph(xyz_data[mid], s['metal'], s.get('donor_elements'))
                    if g is not None:
                        g.complex_id = mid
                        train_3d.append(g)

            test_3d = []
            test_3d_meta = []
            for s in test_samples:
                mid = s.get('complex_id', '')
                if mid in xyz_data:
                    g = build_3d_graph(xyz_data[mid], s['metal'], s.get('donor_elements'))
                    if g is not None:
                        g.complex_id = mid
                        test_3d.append(g)
                        test_3d_meta.append(s)

            if train_3d and test_3d:
                in_dim_3d = train_3d[0].x.size(1)
                egnn = EGNNDonorPredictor(in_dim=in_dim_3d, hidden_dim=128,
                                          n_layers=4).to(device)

                train_loader_3d = DataLoader(train_3d, batch_size=32, shuffle=True)
                val_3d = train_3d[:len(train_3d)//5]  # quick val
                val_loader_3d = DataLoader(val_3d, batch_size=32)
                test_loader_3d = DataLoader(test_3d, batch_size=32)

                egnn = train_gnn_donor(egnn, train_loader_3d, val_loader_3d, device,
                                       epochs=args.epochs, lr=1e-3)
                egnn_results = evaluate_gnn_donor(egnn, test_loader_3d,
                                                  test_3d_meta, device, is_3d=True)
                if egnn_results:
                    s = aggregate_results(egnn_results, 'EGNN_3D_upper')
                    all_summaries.append(s)
                    all_by_dent.extend(aggregate_by_field(
                        egnn_results, 'denticity', 'EGNN_3D_upper'))
                    print(f"  EGNN_3D: F1={s['f1']:.3f} Jaccard={s['jaccard']:.3f}")

                del egnn
                torch.cuda.empty_cache()
            else:
                print("  SKIP: not enough 3D data")
        else:
            print("  SKIP: no XYZ data found")
    else:
        print("\n  Skipping 3D EGNN (--skip_3d)")

    # ── Write outputs ─────────────────────────────────────
    print(f"\n{'='*60}")
    print("Writing outputs")
    print(f"{'='*60}")

    # Summary CSV
    p = os.path.join(args.out, "donor_annotation_summary.csv")
    if all_summaries:
        keys = list(all_summaries[0].keys())
        with open(p, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            for row in all_summaries:
                w.writerow({k: f"{v:.4f}" if isinstance(v, float) else v
                            for k, v in row.items()})
    print(f"  {p}")

    # By denticity
    p = os.path.join(args.out, "donor_annotation_by_denticity.csv")
    if all_by_dent:
        keys = list(all_by_dent[0].keys())
        with open(p, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            for row in all_by_dent:
                w.writerow({k: f"{v:.4f}" if isinstance(v, float) else v
                            for k, v in row.items()})
    print(f"  {p}")

    # By ligand frequency
    p = os.path.join(args.out, "donor_annotation_by_ligand_freq.csv")
    if all_by_freq:
        keys = list(all_by_freq[0].keys())
        with open(p, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            for row in all_by_freq:
                w.writerow({k: f"{v:.4f}" if isinstance(v, float) else v
                            for k, v in row.items()})
    print(f"  {p}")

    # Summary JSON
    summary_json = {
        'models': {s['model']: s for s in all_summaries},
        'leakage': diag,
    }
    p = os.path.join(args.out, "donor_annotation_summary.json")
    with open(p, "w") as f:
        json.dump(summary_json, f, indent=2, default=str)
    print(f"  {p}")

    # Final table
    print(f"\n{'='*60}")
    print("DONOR ANNOTATION SUMMARY")
    print(f"{'='*60}")
    print(f"  {'Model':<25s} {'N':>6s} {'Prec':>7s} {'Rec':>7s} "
          f"{'F1':>7s} {'Jacc':>7s} {'AllCorr':>8s}")
    print(f"  {'-'*70}")
    for s in all_summaries:
        print(f"  {s['model']:<25s} {s['n']:>6d} "
              f"{s['precision']:>7.3f} {s['recall']:>7.3f} "
              f"{s['f1']:>7.3f} {s['jaccard']:>7.3f} "
              f"{s['all_correct']:>8.3f}")


if __name__ == "__main__":
    main()
