#!/usr/bin/env python3
"""
factorized_token_qc_and_standard_eval.py
========================================
1. QC: diagnose why factorized ablation donor Top-1 (14–27%) differs
   from standard Tool A full-context (85%).
2. Re-evaluate composite and factorized models under the STANDARD
   Tool A protocol using fig5_tasks.jsonl.

Standard Tool A protocol:
- Pre-tokenized samples from fig5_tasks.jsonl
- One donor atom masked per sample (single SMILES character)
- y_true is a donor element (C, N, O, S, P, F, I, Cl)
- Prediction space: ~7 donor atom candidates
- Ablation modes: full_context, no_ligand_smiles_keep_length, metal_cn_only
- Stratified by denticity, CN, metal_row
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
import random
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
import torch

from brain.tokenizer import CoordRepTokenizer, TokenizerConfig
from brain.model import CoordRepForMLM, CoordRepModelConfig
from coordrep_tools.infer import load_model_and_tokenizer
from coordrep_tools.ablation_masking import (
    make_ablation_input,
    mask_donor_markers,
    ABLATION_MODES,
)
from coordrep_tools.tokenizer_variants import retokenize_metal_factorized


# ── Constants ─────────────────────────────────────────────

DONOR_ELEMENTS = {'C', 'N', 'O', 'S', 'P', 'F', 'I', 'Cl', 'Br', 'Se', 'As', 'Te'}

METAL_ROW = {}
for m in ['Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn']:
    METAL_ROW[m] = '3d'
for m in ['Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd']:
    METAL_ROW[m] = '4d'
for m in ['La','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg']:
    METAL_ROW[m] = '5d'


# ── QC Diagnostic ─────────────────────────────────────────

def generate_qc_report(fig5_path, out_dir):
    """Generate diagnostic QC comparing old vs new evaluation protocols."""
    with open(fig5_path) as f:
        all_samples = [json.loads(l) for l in f]
    donors = [s for s in all_samples if s.get('task') == 'mask_donor']

    y_counts = Counter(s['y_true']['donor_atom'] for s in donors)

    # Check what token type is at mask_pos
    token_types = Counter()
    for s in donors:
        tok = s['tokens'][s['mask_pos']]
        if re.match(r'^:[A-Z][a-z]?:\d+$', tok):
            token_types['donor_marker'] += 1
        elif len(tok) <= 2 and tok.upper() in DONOR_ELEMENTS:
            token_types['smiles_element'] += 1
        else:
            token_types[f'other:{tok}'] += 1

    report = []
    report.append("# Evaluator QC Report")
    report.append("")
    report.append("## Root Cause of Accuracy Discrepancy")
    report.append("")
    report.append("The factorized ablation script (`factorized_token_ablation.py`)")
    report.append("and the standard Tool A evaluator use **fundamentally different")
    report.append("masking protocols**:")
    report.append("")
    report.append("| Aspect | Standard Tool A | Factorized Ablation |")
    report.append("|--------|----------------|---------------------|")
    report.append("| Test set | fig5_tasks.jsonl (19,992 samples) | Random 2,000 from pipeline |")
    report.append("| Masked token | Single SMILES character (donor atom) | Donor marker (`:N:1` etc.) in constraint block |")
    report.append(f"| y_true type | Element char: C, N, O (from ~{len(y_counts)} classes) | Full token (from {657}+ vocab) |")
    report.append("| Candidates | ~7 donor elements | Entire vocabulary |")
    report.append("| Masks per sample | 1 | All donors simultaneously |")
    report.append("| Tokenization | Pre-tokenized with composite metal block | Re-tokenized at runtime |")
    report.append("")
    report.append("## Impact on Reported Accuracy")
    report.append("")
    report.append("- **Standard Tool A Top-1 = 85.7%**: Predicts which of ~7 elements")
    report.append("  (C, N, O, S, P, F, I) fills the masked SMILES position.")
    report.append("- **Factorized ablation Top-1 = 14.6–27.5%**: Predicts which of 657+")
    report.append("  vocabulary tokens fills the masked donor-marker position.")
    report.append("")
    report.append("These numbers are **not comparable**. The factorized ablation metric")
    report.append("is valid for relative comparison (composite vs factorized), but the")
    report.append("absolute values cannot be compared with Tool A.")
    report.append("")
    report.append("## Token Type at Masked Positions (fig5_tasks)")
    report.append("")
    for ttype, cnt in token_types.most_common():
        report.append(f"- **{ttype}**: {cnt} ({cnt/len(donors)*100:.1f}%)")
    report.append("")
    report.append("## y_true Distribution (fig5_tasks)")
    report.append("")
    for elem, cnt in y_counts.most_common():
        report.append(f"- **{elem}**: {cnt} ({cnt/len(donors)*100:.1f}%)")
    report.append("")
    report.append("## Resolution")
    report.append("")
    report.append("The standard Tool A evaluation must be re-run using:")
    report.append("1. The same fig5_tasks.jsonl test set")
    report.append("2. The same single-donor masking per sample")
    report.append("3. For factorized: re-tokenize the string, remap mask position")
    report.append("4. Restrict Top-k prediction to donor element tokens only")
    report.append("5. Stratify by denticity, CN, and ablation mode")
    report.append("")

    p = os.path.join(out_dir, "evaluator_qc.md")
    with open(p, "w") as f:
        f.write("\n".join(report))
    print(f"  {p}")

    return donors


# ── Standard Evaluation ───────────────────────────────────

def build_factorized_token_map(composite_tokens):
    """
    Given a composite token list, produce factorized token list
    and a position mapping: composite_pos -> factorized_pos.
    """
    fact_tokens = []
    pos_map = {}  # composite_pos -> factorized_pos

    for ci, tok in enumerate(composite_tokens):
        if tok.startswith('[Metal:'):
            # Factorize: [Metal:Fe|ox:+2|d:d6|CN:6] → [Fe], ;ox=+2, ;d=6, ;CN=6
            inner = tok[1:-1]  # Metal:Fe|ox:+2|d:d6|CN:6
            parts = inner.split('|')
            metal = parts[0].replace('Metal:', '')
            pos_map[ci] = len(fact_tokens)
            fact_tokens.append(f'[{metal}]')
            for part in parts[1:]:
                if part.startswith('ox:'):
                    fact_tokens.append(f';ox={part[3:]}')
                elif part.startswith('d:d'):
                    fact_tokens.append(f';d={part[3:]}')
                elif part.startswith('d:'):
                    fact_tokens.append(f';d={part[2:]}')
                elif part.startswith('CN:'):
                    fact_tokens.append(f';CN={part[3:]}')
        else:
            pos_map[ci] = len(fact_tokens)
            fact_tokens.append(tok)

    return fact_tokens, pos_map


def evaluate_standard_protocol(
    model, tokenizer, donor_samples, device,
    ablation_modes, is_factorized=False,
    max_eval=None, label=""):
    """
    Evaluate using standard Tool A protocol:
    - One donor masked per sample
    - Prediction restricted to donor element tokens
    - y_true is donor element character
    """
    mask_id = tokenizer.token2id.get("[MASK]", 3)
    pad_id = tokenizer.token2id.get("[PAD]", 2)
    max_len = 512

    # Build donor element token IDs
    donor_token_ids = {}
    for elem in DONOR_ELEMENTS:
        if elem in tokenizer.token2id:
            donor_token_ids[elem] = tokenizer.token2id[elem]

    results_by_mode = defaultdict(lambda: {
        'n': 0, 'top1': 0, 'top5': 0, 'mrr': 0.0,
        'by_dent': defaultdict(lambda: {'n': 0, 'top1': 0, 'top5': 0}),
        'by_cn': defaultdict(lambda: {'n': 0, 'top1': 0, 'top5': 0}),
        'by_mrow': defaultdict(lambda: {'n': 0, 'top1': 0, 'top5': 0}),
        'per_ligand': defaultdict(list),  # (sample_id, lig_smi) -> list of bools
    })

    examples = []  # For donor_prediction_examples CSV
    n_eval = len(donor_samples) if max_eval is None else min(max_eval, len(donor_samples))

    for si, sample in enumerate(donor_samples[:n_eval]):
        tokens = sample['tokens']
        mask_pos = sample['mask_pos']
        y_true = sample['y_true']['donor_atom']
        group = sample.get('group', {})
        sample_id = sample.get('id', '')

        cn = group.get('cn', 0)
        metal = group.get('metal_element', '?')
        mrow = METAL_ROW.get(metal, '?')

        # Determine denticity from sample group
        dent = group.get('denticity', '?')

        for mode in ablation_modes:
            try:
                abl_tokens, abl_pos_list = make_ablation_input(
                    tokens, [mask_pos], mode, shuffle_seed=42)
            except Exception:
                continue

            if not abl_pos_list:
                continue

            abl_pos = abl_pos_list[0]

            # Factorize if needed
            if is_factorized:
                fact_tokens, pos_map = build_factorized_token_map(abl_tokens)
                new_pos = pos_map.get(abl_pos, abl_pos)
                abl_tokens = fact_tokens
                abl_pos = new_pos

            # Encode
            token_ids = []
            for tok in abl_tokens:
                if tok == '[MASK]':
                    token_ids.append(mask_id)
                elif tok in tokenizer.token2id:
                    token_ids.append(tokenizer.token2id[tok])
                else:
                    token_ids.append(tokenizer.token2id.get('[UNK]', 4))

            if len(token_ids) > max_len:
                if abl_pos >= max_len:
                    continue
                token_ids = token_ids[:max_len]

            # Forward
            padded = token_ids + [pad_id] * (max_len - len(token_ids))
            att = [1] * len(token_ids) + [0] * (max_len - len(token_ids))

            with torch.no_grad():
                inp = torch.tensor([padded], dtype=torch.long, device=device)
                att_t = torch.tensor([att], dtype=torch.long, device=device)
                output = model(inp, att_t)
                logits = output[0] if isinstance(output, tuple) else output

            # Extract donor-element logits only
            pos_logits = logits[0, abl_pos]
            donor_logits = {}
            for elem, tid in donor_token_ids.items():
                donor_logits[elem] = pos_logits[tid].item()

            # Rank by logit
            ranked = sorted(donor_logits.items(), key=lambda x: -x[1])
            top5_elems = [e for e, _ in ranked[:5]]
            top1 = top5_elems[0] if top5_elems else ''

            is_top1 = top1 == y_true
            is_top5 = y_true in top5_elems
            rank = (top5_elems.index(y_true) + 1) if y_true in top5_elems else 0
            rr = 1.0 / rank if rank > 0 else 0.0

            res = results_by_mode[mode]
            res['n'] += 1
            if is_top1:
                res['top1'] += 1
            if is_top5:
                res['top5'] += 1
            res['mrr'] += rr

            # Stratified
            d_key = str(dent)
            res['by_dent'][d_key]['n'] += 1
            if is_top1: res['by_dent'][d_key]['top1'] += 1
            if is_top5: res['by_dent'][d_key]['top5'] += 1

            cn_key = str(cn)
            res['by_cn'][cn_key]['n'] += 1
            if is_top1: res['by_cn'][cn_key]['top1'] += 1
            if is_top5: res['by_cn'][cn_key]['top5'] += 1

            res['by_mrow'][mrow]['n'] += 1
            if is_top1: res['by_mrow'][mrow]['top1'] += 1
            if is_top5: res['by_mrow'][mrow]['top5'] += 1

            # Per-ligand all-correct (only for full_context)
            if mode == 'full_context':
                res['per_ligand'][(sample_id,)].append(is_top1)

            # Examples
            if len(examples) < 200 and mode == 'full_context':
                examples.append({
                    'sample_id': sample_id,
                    'model': label,
                    'y_true': y_true,
                    'top1_pred': top1,
                    'top5_pred': ','.join(top5_elems),
                    'correct_top1': is_top1,
                    'correct_top5': is_top5,
                    'cn': cn,
                    'metal': metal,
                    'denticity': dent,
                })

        if (si + 1) % 2000 == 0:
            print(f"    [{si+1}/{n_eval}] {label}")

    # Aggregate
    metrics = {}
    for mode, res in results_by_mode.items():
        n = res['n']
        if n == 0:
            continue
        metrics[mode] = {
            'n': n,
            'top1': round(res['top1'] / n, 4),
            'top5': round(res['top5'] / n, 4),
            'mrr': round(res['mrr'] / n, 4),
            'by_dent': {},
            'by_cn': {},
            'by_mrow': {},
        }
        for k, v in res['by_dent'].items():
            if v['n'] > 0:
                metrics[mode]['by_dent'][k] = {
                    'n': v['n'],
                    'top1': round(v['top1'] / v['n'], 4),
                    'top5': round(v['top5'] / v['n'], 4),
                }
        for k, v in res['by_cn'].items():
            if v['n'] > 0:
                metrics[mode]['by_cn'][k] = {
                    'n': v['n'],
                    'top1': round(v['top1'] / v['n'], 4),
                    'top5': round(v['top5'] / v['n'], 4),
                }
        for k, v in res['by_mrow'].items():
            if v['n'] > 0:
                metrics[mode]['by_mrow'][k] = {
                    'n': v['n'],
                    'top1': round(v['top1'] / v['n'], 4),
                    'top5': round(v['top5'] / v['n'], 4),
                }
        # Per-ligand all-correct
        all_correct = [all(v) for v in res['per_ligand'].values()]
        metrics[mode]['ligand_all_correct'] = (
            round(sum(all_correct) / len(all_correct), 4) if all_correct else 0)

    return metrics, examples


# ── Main ──────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fig5",
                        default="/data/CoordRep/CoordSMILES/libcoordrep/figures/fig5/data/fig5_tasks.jsonl")
    parser.add_argument("--composite_checkpoint",
                        default="/data/CoordRep/CoordSMILES/libcoordrep/checkpoints/pretrain_v3/best_model.pt")
    parser.add_argument("--composite_tokenizer",
                        default="/data/CoordRep/CoordSMILES/libcoordrep/checkpoints/pretrain_v3/tokenizer.json")
    parser.add_argument("--matched_checkpoint",
                        default="checkpoints/token_ablation/composite/best_model.pt")
    parser.add_argument("--matched_tokenizer",
                        default="checkpoints/token_ablation/composite/tokenizer.json")
    parser.add_argument("--factorized_checkpoint",
                        default="checkpoints/token_ablation/factorized/best_model.pt")
    parser.add_argument("--factorized_tokenizer",
                        default="checkpoints/token_ablation/factorized/tokenizer.json")
    parser.add_argument("--out", default="revision_results/factorized_token")
    parser.add_argument("--device", default="cuda:0")
    parser.add_argument("--max_eval", type=int, default=5000,
                        help="Max samples per model (full set=19992)")
    parser.add_argument("--modes", nargs="+",
                        default=["full_context", "no_ligand_smiles_keep_length", "metal_cn_only"])
    args = parser.parse_args()

    os.makedirs(args.out, exist_ok=True)

    # ── Step 1: QC Report ─────────────────────────────────
    print("Generating QC report …")
    donors = generate_qc_report(args.fig5, args.out)
    print(f"  {len(donors)} donor samples")

    # ── Step 2: Write old vs new example comparison ───────
    # (We generate this during eval below)

    # ── Step 3: Evaluate all three models ─────────────────
    all_examples = []
    all_metrics = {}

    # 3a. Composite pretrained
    print(f"\n{'='*60}")
    print("Evaluating: COMPOSITE PRETRAINED (original)")
    print(f"{'='*60}")
    comp_model, comp_tok = load_model_and_tokenizer(
        args.composite_checkpoint, args.composite_tokenizer, args.device)
    comp_metrics, comp_examples = evaluate_standard_protocol(
        comp_model, comp_tok, donors, args.device,
        args.modes, is_factorized=False, max_eval=args.max_eval,
        label="composite_pretrained")
    all_metrics['composite_pretrained'] = comp_metrics
    all_examples.extend(comp_examples)
    del comp_model
    torch.cuda.empty_cache()

    for mode, m in comp_metrics.items():
        print(f"  {mode}: Top-1={m['top1']}, Top-5={m['top5']}, MRR={m['mrr']}")

    # 3b. Composite matched training
    print(f"\n{'='*60}")
    print("Evaluating: COMPOSITE MATCHED TRAINING")
    print(f"{'='*60}")
    if os.path.exists(args.matched_checkpoint):
        matched_tok = CoordRepTokenizer.load(args.matched_tokenizer)
        model_config = CoordRepModelConfig.small(max_length=768)
        model_config.vocab_size = 10000
        matched_model = CoordRepForMLM(model_config)
        ckpt = torch.load(args.matched_checkpoint, map_location='cpu', weights_only=False)
        matched_model.load_state_dict(ckpt['model_state_dict'])
        matched_model.eval().to(args.device)

        matched_metrics, matched_examples = evaluate_standard_protocol(
            matched_model, matched_tok, donors, args.device,
            args.modes, is_factorized=False, max_eval=args.max_eval,
            label="composite_matched")
        all_metrics['composite_matched'] = matched_metrics
        all_examples.extend(matched_examples)
        del matched_model
        torch.cuda.empty_cache()

        for mode, m in matched_metrics.items():
            print(f"  {mode}: Top-1={m['top1']}, Top-5={m['top5']}, MRR={m['mrr']}")
    else:
        print(f"  SKIP: {args.matched_checkpoint} not found")
        all_metrics['composite_matched'] = {}

    # 3c. Factorized
    print(f"\n{'='*60}")
    print("Evaluating: FACTORIZED")
    print(f"{'='*60}")
    if os.path.exists(args.factorized_checkpoint):
        fact_tok = CoordRepTokenizer.load(args.factorized_tokenizer)
        model_config = CoordRepModelConfig.small(max_length=768)
        model_config.vocab_size = 10000
        fact_model = CoordRepForMLM(model_config)
        ckpt = torch.load(args.factorized_checkpoint, map_location='cpu', weights_only=False)
        fact_model.load_state_dict(ckpt['model_state_dict'])
        fact_model.eval().to(args.device)

        fact_metrics, fact_examples = evaluate_standard_protocol(
            fact_model, fact_tok, donors, args.device,
            args.modes, is_factorized=True, max_eval=args.max_eval,
            label="factorized")
        all_metrics['factorized'] = fact_metrics
        all_examples.extend(fact_examples)
        del fact_model
        torch.cuda.empty_cache()

        for mode, m in fact_metrics.items():
            print(f"  {mode}: Top-1={m['top1']}, Top-5={m['top5']}, MRR={m['mrr']}")
    else:
        print(f"  SKIP: {args.factorized_checkpoint} not found")
        all_metrics['factorized'] = {}

    # ── Step 4: Write outputs ─────────────────────────────
    print(f"\n{'='*60}")
    print("Writing outputs")
    print(f"{'='*60}")

    # 4a. tool_a_standard_eval_comparison.csv
    p = os.path.join(args.out, "tool_a_standard_eval_comparison.csv")
    with open(p, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["model", "ablation_mode", "n", "top1", "top5", "mrr",
                     "ligand_all_correct"])
        for model_name, mode_dict in all_metrics.items():
            for mode, m in mode_dict.items():
                w.writerow([
                    model_name, mode,
                    m.get('n', 0),
                    m.get('top1', 0),
                    m.get('top5', 0),
                    m.get('mrr', 0),
                    m.get('ligand_all_correct', 0),
                ])
    print(f"  {p}")

    # 4b. By denticity
    p = os.path.join(args.out, "factorized_by_denticity.csv")
    with open(p, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["model", "ablation_mode", "denticity", "n", "top1", "top5"])
        for model_name, mode_dict in all_metrics.items():
            for mode, m in mode_dict.items():
                for dent, dv in sorted(m.get('by_dent', {}).items()):
                    w.writerow([model_name, mode, dent, dv['n'], dv['top1'], dv['top5']])
    print(f"  {p}")

    # 4c. By CN
    p = os.path.join(args.out, "factorized_by_cn.csv")
    with open(p, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["model", "ablation_mode", "cn", "n", "top1", "top5"])
        for model_name, mode_dict in all_metrics.items():
            for mode, m in mode_dict.items():
                for cn_k, cv in sorted(m.get('by_cn', {}).items()):
                    w.writerow([model_name, mode, cn_k, cv['n'], cv['top1'], cv['top5']])
    print(f"  {p}")

    # 4d. By metal row
    p = os.path.join(args.out, "factorized_by_metal_row.csv")
    with open(p, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["model", "ablation_mode", "metal_row", "n", "top1", "top5"])
        for model_name, mode_dict in all_metrics.items():
            for mode, m in mode_dict.items():
                for mrow, mv in sorted(m.get('by_mrow', {}).items()):
                    w.writerow([model_name, mode, mrow, mv['n'], mv['top1'], mv['top5']])
    print(f"  {p}")

    # 4e. Donor prediction examples (old vs new)
    p = os.path.join(args.out, "donor_prediction_examples_old_vs_new.csv")
    with open(p, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["sample_id", "model", "y_true", "top1_pred", "top5_pred",
                     "correct_top1", "correct_top5", "cn", "metal", "denticity"])
        for ex in all_examples:
            w.writerow([
                ex['sample_id'], ex['model'], ex['y_true'],
                ex['top1_pred'], ex['top5_pred'],
                ex['correct_top1'], ex['correct_top5'],
                ex['cn'], ex['metal'], ex['denticity'],
            ])
    print(f"  {p}")

    # 4f. Summary JSON
    summary = {}
    for model_name, mode_dict in all_metrics.items():
        fc = mode_dict.get('full_context', {})
        nls = mode_dict.get('no_ligand_smiles_keep_length', {})
        mco = mode_dict.get('metal_cn_only', {})
        summary[f'{model_name}_full_context_top1'] = fc.get('top1', 0)
        summary[f'{model_name}_full_context_top5'] = fc.get('top5', 0)
        summary[f'{model_name}_full_context_mrr'] = fc.get('mrr', 0)
        summary[f'{model_name}_no_lig_smiles_top1'] = nls.get('top1', 0)
        summary[f'{model_name}_metal_cn_only_top1'] = mco.get('top1', 0)
        summary[f'{model_name}_ligand_all_correct'] = fc.get('ligand_all_correct', 0)

    p = os.path.join(args.out, "summary.json")
    with open(p, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"  {p}")

    # ── Print final comparison ────────────────────────────
    print(f"\n{'='*60}")
    print("STANDARD TOOL A PROTOCOL COMPARISON")
    print(f"{'='*60}")
    print(f"  {'Model':<25s} {'Mode':<32s} {'N':>6s} {'Top1':>7s} {'Top5':>7s} {'MRR':>7s}")
    print(f"  {'-'*85}")
    for model_name, mode_dict in all_metrics.items():
        for mode, m in sorted(mode_dict.items()):
            print(f"  {model_name:<25s} {mode:<32s} {m['n']:>6d} "
                  f"{m['top1']:>7.4f} {m['top5']:>7.4f} {m['mrr']:>7.4f}")
    print()


if __name__ == "__main__":
    main()
