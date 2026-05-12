#!/usr/bin/env python3
"""
factorized_token_ablation.py
============================
Tasks 2-4: Build factorized tokenizer, train minimal MLM,
and evaluate field-masking for both composite and factorized variants.

Strategy:
- composite: use existing pretrained checkpoint (production)
- factorized: re-tokenize data with factorized metal block,
  build new tokenizer, and do continued MLM training from
  scratch with same architecture/hyperparams

Evaluation:
A. Donor-marker completion (Top-1/5, MRR, by denticity)
B. Metal masking (Top-1/5)
C. CN masking (Top-1/5)
D. Joint metal+CN masking (exact joint, individual)
E. Shape field masking (optional)
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
import random
import time
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader

from brain.tokenizer import CoordRepTokenizer, TokenizerConfig, SPECIAL_TOKENS
from brain.model import CoordRepForMLM, CoordRepModelConfig
from coordrep_tools.tokenizer_variants import retokenize_metal_factorized
from coordrep_tools.infer import load_model_and_tokenizer


# ── Constants ─────────────────────────────────────────────

SEED = 42
METAL_ROW = {}
for m in ['Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn']:
    METAL_ROW[m] = '3d'
for m in ['Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd']:
    METAL_ROW[m] = '4d'
for m in ['La','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg']:
    METAL_ROW[m] = '5d'


# ── Data Loading ──────────────────────────────────────────

def load_coordrep_strings(path, max_entries=0):
    """Load CoordRep strings from pipeline output."""
    records = []
    with open(path) as f:
        for i, line in enumerate(f):
            if max_entries > 0 and i >= max_entries:
                break
            d = json.loads(line)
            cr = d.get('coordrep', '')
            if not cr:
                continue
            records.append({
                'mol_id': d.get('mol_id', ''),
                'coordrep': cr,
                'metal': d.get('metal', '?'),
                'cn': d.get('cn', 0),
                'ligand_dents': d.get('ligand_dents', []),
            })
    return records


def build_factorized_tokenizer(records, vocab_size=10000):
    """
    Build a new tokenizer with factorized metal tokens.
    Re-tokenize all strings to factorized format first.
    """
    config = TokenizerConfig(vocab_size=vocab_size)
    tok = CoordRepTokenizer(config)

    # Tokenize all records to build vocab
    for rec in records:
        factorized_str = retokenize_metal_factorized(rec['coordrep'])
        tok.encode(factorized_str, add_special=True)

    return tok


# ── MLM Dataset ───────────────────────────────────────────

class MLMDataset(Dataset):
    """Simple MLM dataset for token ablation."""

    def __init__(self, records, tokenizer, max_len=512, mask_rate=0.15,
                 factorized=False):
        self.tokenizer = tokenizer
        self.max_len = max_len
        self.mask_rate = mask_rate
        self.factorized = factorized

        self.samples = []
        mask_id = tokenizer.token2id.get("[MASK]", 3)
        pad_id = tokenizer.token2id.get("[PAD]", 2)

        for rec in records:
            cr = rec['coordrep']
            if factorized:
                cr = retokenize_metal_factorized(cr)
            ids = tokenizer.encode(cr, add_special=True)
            if len(ids) > max_len:
                ids = ids[:max_len]
            self.samples.append(ids)

        self.mask_id = mask_id
        self.pad_id = pad_id

    def __len__(self):
        return len(self.samples)

    def __getitem__(self, idx):
        ids = list(self.samples[idx])
        labels = [-100] * len(ids)

        for i in range(1, len(ids) - 1):  # skip CLS/SEP
            if random.random() < self.mask_rate:
                labels[i] = ids[i]
                ids[i] = self.mask_id

        # Pad
        pad_len = self.max_len - len(ids)
        attention_mask = [1] * len(ids) + [0] * pad_len
        ids = ids + [self.pad_id] * pad_len
        labels = labels + [-100] * pad_len

        return {
            'input_ids': torch.tensor(ids, dtype=torch.long),
            'attention_mask': torch.tensor(attention_mask, dtype=torch.long),
            'labels': torch.tensor(labels, dtype=torch.long),
        }


# ── Training ──────────────────────────────────────────────

def train_mlm(model, train_dataset, val_dataset, device, epochs=5,
              batch_size=32, lr=1e-4, save_dir=None):
    """Train MLM with early stopping."""
    train_loader = DataLoader(train_dataset, batch_size=batch_size,
                              shuffle=True, num_workers=0)
    val_loader = DataLoader(val_dataset, batch_size=batch_size,
                            shuffle=False, num_workers=0)

    optimizer = torch.optim.AdamW(model.parameters(), lr=lr)
    loss_fn = nn.CrossEntropyLoss(ignore_index=-100)

    best_val_loss = float('inf')
    patience = 2
    patience_counter = 0

    for epoch in range(epochs):
        model.train()
        total_loss = 0
        n_batches = 0

        for batch in train_loader:
            input_ids = batch['input_ids'].to(device)
            attention_mask = batch['attention_mask'].to(device)
            labels = batch['labels'].to(device)

            output = model(input_ids, attention_mask)
            logits = output[0] if isinstance(output, tuple) else output
            loss = loss_fn(logits.view(-1, logits.size(-1)), labels.view(-1))

            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

            total_loss += loss.item()
            n_batches += 1

        train_loss = total_loss / max(n_batches, 1)

        # Validation
        model.eval()
        val_loss = 0
        n_val = 0
        with torch.no_grad():
            for batch in val_loader:
                input_ids = batch['input_ids'].to(device)
                attention_mask = batch['attention_mask'].to(device)
                labels = batch['labels'].to(device)
                output = model(input_ids, attention_mask)
                logits = output[0] if isinstance(output, tuple) else output
                loss = loss_fn(logits.view(-1, logits.size(-1)), labels.view(-1))
                val_loss += loss.item()
                n_val += 1
        val_loss /= max(n_val, 1)

        print(f"    Epoch {epoch+1}: train_loss={train_loss:.4f} val_loss={val_loss:.4f}")

        if val_loss < best_val_loss:
            best_val_loss = val_loss
            patience_counter = 0
            if save_dir:
                os.makedirs(save_dir, exist_ok=True)
                torch.save({
                    'model_state_dict': model.state_dict(),
                    'epoch': epoch,
                    'val_loss': val_loss,
                }, os.path.join(save_dir, "best_model.pt"))
        else:
            patience_counter += 1
            if patience_counter >= patience:
                print(f"    Early stopping at epoch {epoch+1}")
                break

    return best_val_loss


# ── Evaluation Helpers ────────────────────────────────────

def identify_field_positions(tokens, tokenizer):
    """
    Identify which token positions correspond to which fields.
    Returns dict of field_name -> list of positions.
    """
    fields = {
        'metal': [],
        'cn': [],
        'ox': [],
        'd_count': [],
        'donor': [],
        'shape': [],
        'ligand_smiles': [],
    }

    for i, tok in enumerate(tokens):
        t = tokenizer.id2token.get(tok, '')
        if t.startswith('[Metal:') or re.match(r'^\[[A-Z][a-z]?\]$', t):
            fields['metal'].append(i)
        elif t.startswith(';CN=') or 'CN:' in t:
            fields['cn'].append(i)
        elif t.startswith(';ox=') or 'ox:' in t:
            fields['ox'].append(i)
        elif t.startswith(';d=') or 'd:d' in t:
            fields['d_count'].append(i)
        elif ':N:' in t or ':O:' in t or ':S:' in t or ':P:' in t or ':Cl:' in t:
            fields['donor'].append(i)
        elif t in ('Td', 'SP', 'Oh', 'TP', 'TBP', 'SPY', 'TPr', 'L'):
            fields['shape'].append(i)

    return fields


def evaluate_field_masking(model, tokenizer, records, device,
                           factorized=False, max_eval=2000):
    """
    Evaluate field-specific masking.
    Returns dict of task -> metrics.
    """
    model.eval()
    mask_id = tokenizer.token2id.get("[MASK]", 3)
    pad_id = tokenizer.token2id.get("[PAD]", 2)
    max_len = 512

    results = defaultdict(lambda: {'n': 0, 'top1': 0, 'top5': 0, 'mrr': 0.0})
    dent_results = defaultdict(lambda: {'n': 0, 'top1': 0, 'top5': 0})

    # For joint metal+CN
    joint_results = {'n': 0, 'joint_top1': 0, 'metal_top1': 0, 'cn_top1': 0}

    for rec_idx, rec in enumerate(records[:max_eval]):
        cr = rec['coordrep']
        if factorized:
            cr = retokenize_metal_factorized(cr)

        ids = tokenizer.encode(cr, add_special=True)
        if len(ids) > max_len:
            ids = ids[:max_len]

        tokens_str = [tokenizer.id2token.get(t, '') for t in ids]

        # Identify fields
        fields = {}
        metal_pos = []
        cn_pos = []
        donor_pos = []

        for i, t in enumerate(tokens_str):
            if factorized:
                if re.match(r'^\[[A-Z][a-z]?\]$', t):
                    metal_pos.append(i)
                elif t.startswith(';CN='):
                    cn_pos.append(i)
            else:
                if t.startswith('[Metal:'):
                    metal_pos.append(i)
                    cn_pos.append(i)  # same token in composite

            if ':N:' in t or ':O:' in t or ':S:' in t or ':P:' in t or ':Cl:' in t:
                donor_pos.append(i)

        # ── Task A: Donor masking ─────────────────────────
        for dp in donor_pos:
            masked_ids = list(ids)
            true_id = masked_ids[dp]
            masked_ids[dp] = mask_id

            padded = masked_ids + [pad_id] * (max_len - len(masked_ids))
            att = [1] * len(masked_ids) + [0] * (max_len - len(masked_ids))

            with torch.no_grad():
                inp = torch.tensor([padded], dtype=torch.long, device=device)
                att_t = torch.tensor([att], dtype=torch.long, device=device)
                output = model(inp, att_t)
                logits = output[0] if isinstance(output, tuple) else output

            probs = logits[0, dp]
            topk = probs.topk(5).indices.tolist()

            results['donor']['n'] += 1
            if topk[0] == true_id:
                results['donor']['top1'] += 1
            if true_id in topk:
                results['donor']['top5'] += 1
            rank = (probs.argsort(descending=True) == true_id).nonzero(as_tuple=True)[0]
            if len(rank) > 0:
                results['donor']['mrr'] += 1.0 / (rank[0].item() + 1)

        # ── Task B: Metal masking ─────────────────────────
        if metal_pos:
            for mp in metal_pos[:1]:  # just first
                masked_ids = list(ids)
                true_id = masked_ids[mp]
                masked_ids[mp] = mask_id

                padded = masked_ids + [pad_id] * (max_len - len(masked_ids))
                att = [1] * len(masked_ids) + [0] * (max_len - len(masked_ids))

                with torch.no_grad():
                    inp = torch.tensor([padded], dtype=torch.long, device=device)
                    att_t = torch.tensor([att], dtype=torch.long, device=device)
                    output = model(inp, att_t)
                    logits = output[0] if isinstance(output, tuple) else output

                probs = logits[0, mp]
                topk = probs.topk(5).indices.tolist()

                results['metal']['n'] += 1
                if topk[0] == true_id:
                    results['metal']['top1'] += 1
                if true_id in topk:
                    results['metal']['top5'] += 1

        # ── Task C: CN masking ────────────────────────────
        if cn_pos and factorized:
            for cp in cn_pos[:1]:
                masked_ids = list(ids)
                true_id = masked_ids[cp]
                masked_ids[cp] = mask_id

                padded = masked_ids + [pad_id] * (max_len - len(masked_ids))
                att = [1] * len(masked_ids) + [0] * (max_len - len(masked_ids))

                with torch.no_grad():
                    inp = torch.tensor([padded], dtype=torch.long, device=device)
                    att_t = torch.tensor([att], dtype=torch.long, device=device)
                    output = model(inp, att_t)
                    logits = output[0] if isinstance(output, tuple) else output

                probs = logits[0, cp]
                topk = probs.topk(5).indices.tolist()

                results['cn']['n'] += 1
                if topk[0] == true_id:
                    results['cn']['top1'] += 1
                if true_id in topk:
                    results['cn']['top5'] += 1
        elif cn_pos and not factorized:
            # Composite: metal and CN are same token, so CN masking = metal masking
            results['cn']['n'] += 1
            # Copy metal result
            if results['metal']['top1'] > results['cn'].get('_prev_top1', 0):
                results['cn']['top1'] = results['metal']['top1']
                results['cn']['top5'] = results['metal']['top5']
                results['cn']['_prev_top1'] = results['metal']['top1']

        # ── Task D: Joint metal+CN ────────────────────────
        if factorized and metal_pos and cn_pos:
            masked_ids = list(ids)
            metal_true = masked_ids[metal_pos[0]]
            cn_true = masked_ids[cn_pos[0]]
            masked_ids[metal_pos[0]] = mask_id
            masked_ids[cn_pos[0]] = mask_id

            padded = masked_ids + [pad_id] * (max_len - len(masked_ids))
            att = [1] * len(masked_ids) + [0] * (max_len - len(masked_ids))

            with torch.no_grad():
                inp = torch.tensor([padded], dtype=torch.long, device=device)
                att_t = torch.tensor([att], dtype=torch.long, device=device)
                output = model(inp, att_t)
                logits = output[0] if isinstance(output, tuple) else output

            metal_pred = logits[0, metal_pos[0]].argmax().item()
            cn_pred = logits[0, cn_pos[0]].argmax().item()

            joint_results['n'] += 1
            if metal_pred == metal_true:
                joint_results['metal_top1'] += 1
            if cn_pred == cn_true:
                joint_results['cn_top1'] += 1
            if metal_pred == metal_true and cn_pred == cn_true:
                joint_results['joint_top1'] += 1
        elif not factorized and metal_pos:
            # Composite: masking the single token recovers both
            joint_results['n'] += 1
            # Same as metal masking
            if results['metal']['n'] > 0:
                last_hit = results['metal']['top1']
                joint_results['metal_top1'] = last_hit
                joint_results['cn_top1'] = last_hit
                joint_results['joint_top1'] = last_hit

        if (rec_idx + 1) % 500 == 0:
            print(f"    [{rec_idx+1}/{min(len(records), max_eval)}]")

    # Compute rates
    metrics = {}
    for task in ['donor', 'metal', 'cn']:
        n = results[task]['n']
        if n > 0:
            metrics[f'{task}_top1'] = round(results[task]['top1'] / n, 4)
            metrics[f'{task}_top5'] = round(results[task]['top5'] / n, 4)
            if 'mrr' in results[task]:
                metrics[f'{task}_mrr'] = round(results[task]['mrr'] / n, 4)
            metrics[f'{task}_n'] = n

    n_j = joint_results['n']
    if n_j > 0:
        metrics['joint_metal_cn_top1'] = round(joint_results['joint_top1'] / n_j, 4)
        metrics['joint_metal_only_top1'] = round(joint_results['metal_top1'] / n_j, 4)
        metrics['joint_cn_only_top1'] = round(joint_results['cn_top1'] / n_j, 4)

    return metrics


# ── Main ──────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data",
                        default="/data/CoordRep/CoordSMILES/pipeline_full_output/results.jsonl")
    parser.add_argument("--checkpoint",
                        default="/data/CoordRep/CoordSMILES/libcoordrep/checkpoints/pretrain_v3/best_model.pt")
    parser.add_argument("--tokenizer_path",
                        default="/data/CoordRep/CoordSMILES/libcoordrep/checkpoints/pretrain_v3/tokenizer.json")
    parser.add_argument("--out", default="revision_results/factorized_token")
    parser.add_argument("--device", default="cuda:0")
    parser.add_argument("--n_train", type=int, default=10000,
                        help="Training samples for factorized variant")
    parser.add_argument("--n_eval", type=int, default=2000)
    parser.add_argument("--epochs", type=int, default=5)
    parser.add_argument("--batch_size", type=int, default=32)
    parser.add_argument("--lr", type=float, default=5e-5)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    os.makedirs(args.out, exist_ok=True)
    random.seed(args.seed)
    np.random.seed(args.seed)
    torch.manual_seed(args.seed)

    device = args.device if torch.cuda.is_available() else "cpu"

    # ── Load data ─────────────────────────────────────────
    print("Loading data …")
    all_records = load_coordrep_strings(args.data)
    random.shuffle(all_records)

    n_total = len(all_records)
    n_train = min(args.n_train, int(n_total * 0.8))
    n_val = min(2000, int(n_total * 0.1))
    n_test = min(args.n_eval, n_total - n_train - n_val)

    train_records = all_records[:n_train]
    val_records = all_records[n_train:n_train + n_val]
    test_records = all_records[n_train + n_val:n_train + n_val + n_test]

    print(f"  train={len(train_records)}, val={len(val_records)}, test={len(test_records)}")

    # ══════════════════════════════════════════════════════
    # COMPOSITE VARIANT (existing pretrained model)
    # ══════════════════════════════════════════════════════
    print(f"\n{'='*60}")
    print("COMPOSITE VARIANT (existing pretrained)")
    print(f"{'='*60}")

    comp_model, comp_tokenizer = load_model_and_tokenizer(
        args.checkpoint, args.tokenizer_path, device)

    print("  Evaluating composite variant …")
    comp_metrics = evaluate_field_masking(
        comp_model, comp_tokenizer, test_records, device,
        factorized=False, max_eval=args.n_eval)

    print("  Composite metrics:")
    for k, v in sorted(comp_metrics.items()):
        print(f"    {k}: {v}")

    # ══════════════════════════════════════════════════════
    # FACTORIZED VARIANT (new tokenizer + continued training)
    # ══════════════════════════════════════════════════════
    print(f"\n{'='*60}")
    print("FACTORIZED VARIANT (new tokenizer + training)")
    print(f"{'='*60}")

    print("  Building factorized tokenizer …")
    fact_tokenizer = build_factorized_tokenizer(all_records[:20000])
    print(f"    Vocab size: {len(fact_tokenizer.token2id)}")

    # Save tokenizer
    fact_tok_dir = os.path.join(args.out, "factorized_tokenizer")
    os.makedirs(fact_tok_dir, exist_ok=True)
    fact_tokenizer.save(os.path.join(fact_tok_dir, "tokenizer.json"))

    # Build datasets
    print("  Building datasets …")
    train_ds = MLMDataset(train_records, fact_tokenizer, factorized=True)
    val_ds = MLMDataset(val_records, fact_tokenizer, factorized=True)

    # Create model
    model_config = CoordRepModelConfig.small(max_length=768)
    model_config.vocab_size = 10000
    fact_model = CoordRepForMLM(model_config).to(device)

    fact_save_dir = "checkpoints/token_ablation/factorized"

    print(f"  Training factorized MLM ({args.epochs} epochs, {len(train_ds)} samples) …")
    t0 = time.time()
    best_loss = train_mlm(
        fact_model, train_ds, val_ds, device,
        epochs=args.epochs, batch_size=args.batch_size,
        lr=args.lr, save_dir=fact_save_dir)
    print(f"  Done ({time.time()-t0:.0f}s), best_val_loss={best_loss:.4f}")

    # Load best checkpoint
    ckpt = torch.load(os.path.join(fact_save_dir, "best_model.pt"),
                      map_location='cpu', weights_only=False)
    fact_model.load_state_dict(ckpt['model_state_dict'])
    fact_model.eval().to(device)

    # Save tokenizer alongside checkpoint
    fact_tokenizer.save(os.path.join(fact_save_dir, "tokenizer.json"))

    print("  Evaluating factorized variant …")
    fact_metrics = evaluate_field_masking(
        fact_model, fact_tokenizer, test_records, device,
        factorized=True, max_eval=args.n_eval)

    print("  Factorized metrics:")
    for k, v in sorted(fact_metrics.items()):
        print(f"    {k}: {v}")

    # Also train composite from scratch for fair comparison
    print(f"\n{'='*60}")
    print("COMPOSITE VARIANT (matched training for fair comparison)")
    print(f"{'='*60}")

    comp_train_ds = MLMDataset(train_records, comp_tokenizer, factorized=False)
    comp_val_ds = MLMDataset(val_records, comp_tokenizer, factorized=False)

    comp_scratch_model = CoordRepForMLM(model_config).to(device)
    comp_save_dir = "checkpoints/token_ablation/composite"

    print(f"  Training composite MLM ({args.epochs} epochs, {len(comp_train_ds)} samples) …")
    t0 = time.time()
    best_loss_comp = train_mlm(
        comp_scratch_model, comp_train_ds, comp_val_ds, device,
        epochs=args.epochs, batch_size=args.batch_size,
        lr=args.lr, save_dir=comp_save_dir)
    print(f"  Done ({time.time()-t0:.0f}s), best_val_loss={best_loss_comp:.4f}")

    ckpt_c = torch.load(os.path.join(comp_save_dir, "best_model.pt"),
                        map_location='cpu', weights_only=False)
    comp_scratch_model.load_state_dict(ckpt_c['model_state_dict'])
    comp_scratch_model.eval().to(device)
    comp_tokenizer.save(os.path.join(comp_save_dir, "tokenizer.json"))

    print("  Evaluating matched-training composite …")
    comp_matched_metrics = evaluate_field_masking(
        comp_scratch_model, comp_tokenizer, test_records, device,
        factorized=False, max_eval=args.n_eval)

    # ══════════════════════════════════════════════════════
    # WRITE OUTPUTS
    # ══════════════════════════════════════════════════════
    print(f"\n{'='*60}")
    print("WRITING OUTPUTS")
    print(f"{'='*60}")

    # ── factorized_ablation_summary.csv ───────────────────
    p = os.path.join(args.out, "factorized_ablation_summary.csv")
    with open(p, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["variant", "metric", "value"])
        for k, v in sorted(comp_metrics.items()):
            w.writerow(["composite_pretrained", k, v])
        for k, v in sorted(comp_matched_metrics.items()):
            w.writerow(["composite_matched", k, v])
        for k, v in sorted(fact_metrics.items()):
            w.writerow(["factorized", k, v])
    print(f"  {p}")

    # ── factorized_by_task.csv ────────────────────────────
    tasks = ['donor', 'metal', 'cn']
    p = os.path.join(args.out, "factorized_by_task.csv")
    with open(p, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["task", "composite_pretrained_top1", "composite_matched_top1",
                     "factorized_top1", "composite_pretrained_top5",
                     "composite_matched_top5", "factorized_top5"])
        for task in tasks:
            w.writerow([
                task,
                comp_metrics.get(f'{task}_top1', ''),
                comp_matched_metrics.get(f'{task}_top1', ''),
                fact_metrics.get(f'{task}_top1', ''),
                comp_metrics.get(f'{task}_top5', ''),
                comp_matched_metrics.get(f'{task}_top5', ''),
                fact_metrics.get(f'{task}_top5', ''),
            ])
    print(f"  {p}")

    # ── summary.json ──────────────────────────────────────
    summary = {
        "composite_donor_top1": comp_metrics.get('donor_top1', 0),
        "factorized_donor_top1": fact_metrics.get('donor_top1', 0),
        "composite_donor_top5": comp_metrics.get('donor_top5', 0),
        "factorized_donor_top5": fact_metrics.get('donor_top5', 0),
        "composite_metal_top1": comp_metrics.get('metal_top1', 0),
        "factorized_metal_top1": fact_metrics.get('metal_top1', 0),
        "composite_cn_top1": comp_metrics.get('cn_top1', 0),
        "factorized_cn_top1": fact_metrics.get('cn_top1', 0),
        "composite_joint_metal_cn_top1": comp_metrics.get('joint_metal_cn_top1', 0),
        "factorized_joint_metal_cn_top1": fact_metrics.get('joint_metal_cn_top1', 0),
        "composite_matched_donor_top1": comp_matched_metrics.get('donor_top1', 0),
        "composite_matched_metal_top1": comp_matched_metrics.get('metal_top1', 0),
        "n_train": len(train_records),
        "n_test": len(test_records),
        "epochs": args.epochs,
        "seed": args.seed,
    }

    # Determine recommendation
    comp_d = comp_matched_metrics.get('donor_top1', 0)
    fact_d = fact_metrics.get('donor_top1', 0)
    fact_cn = fact_metrics.get('cn_top1', 0)
    fact_m = fact_metrics.get('metal_top1', 0)

    if fact_d >= comp_d - 0.02 and fact_cn > 0:
        summary['recommendation'] = 'use_factorized'
        summary['rationale'] = ('Factorized maintains donor accuracy while enabling '
                               'independent metal/CN masking and evaluation.')
    elif comp_d > fact_d + 0.05:
        summary['recommendation'] = 'keep_composite'
        summary['rationale'] = ('Composite has substantially better donor accuracy; '
                               'the metal-CN prior aids prediction.')
    else:
        summary['recommendation'] = 'both_supported'
        summary['rationale'] = ('Performance is comparable; factorized offers '
                               'interpretability benefits.')

    p = os.path.join(args.out, "summary.json")
    with open(p, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"  {p}")

    # ── Print final comparison ────────────────────────────
    print(f"\n{'='*60}")
    print("FACTORIZED TOKEN ABLATION — RESULTS")
    print(f"{'='*60}")
    print(f"  {'Task':<20s} {'Comp(pretrained)':<18s} {'Comp(matched)':<18s} {'Factorized':<18s}")
    print(f"  {'-'*74}")
    for task in ['donor', 'metal', 'cn']:
        cp = comp_metrics.get(f'{task}_top1', '-')
        cm = comp_matched_metrics.get(f'{task}_top1', '-')
        ft = fact_metrics.get(f'{task}_top1', '-')
        print(f"  {task+'_top1':<20s} {str(cp):<18s} {str(cm):<18s} {str(ft):<18s}")
    jc = comp_metrics.get('joint_metal_cn_top1', '-')
    jf = fact_metrics.get('joint_metal_cn_top1', '-')
    jcm = comp_matched_metrics.get('joint_metal_cn_top1', '-')
    print(f"  {'joint_top1':<20s} {str(jc):<18s} {str(jcm):<18s} {str(jf):<18s}")
    print(f"\n  Recommendation: {summary['recommendation']}")
    print(f"  Rationale: {summary['rationale']}")
    print()


if __name__ == "__main__":
    main()
