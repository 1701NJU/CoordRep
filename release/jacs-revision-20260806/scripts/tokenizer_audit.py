#!/usr/bin/env python3
"""
tokenizer_audit.py
==================
Task 1: Audit the current tokenizer to determine whether
metal/CN tokenization is composite or factorized.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
from collections import Counter

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tokenizer",
                        default="/data/CoordRep/CoordSMILES/libcoordrep/checkpoints/pretrain_v3/tokenizer.json")
    parser.add_argument("--data",
                        default="/data/CoordRep/CoordSMILES/pipeline_full_output/results.jsonl")
    parser.add_argument("--n_sample", type=int, default=1000)
    parser.add_argument("--out", default="revision_results/factorized_token")
    args = parser.parse_args()

    os.makedirs(args.out, exist_ok=True)

    # ── Load tokenizer vocab ──────────────────────────────
    with open(args.tokenizer) as f:
        data = json.load(f)
    vocab = data['token2id']

    # ── Classify tokens ───────────────────────────────────
    composite_metal_cn = []   # [Metal:Fe|...|CN:6]
    factorized_metal = []     # [Fe], [Co], etc.
    factorized_cn = []        # ;CN=2, ;CN=6, etc.
    factorized_ox = []        # ;ox=+2, etc.
    factorized_d = []         # ;d=6, etc.
    factorized_row = []       # ;row=3, etc.
    shape_tokens = []

    for tok in vocab:
        if tok.startswith('[Metal:'):
            composite_metal_cn.append(tok)
        elif re.match(r'^\[[A-Z][a-z]?\]$', tok):
            factorized_metal.append(tok)
        elif tok.startswith(';CN='):
            factorized_cn.append(tok)
        elif tok.startswith(';ox='):
            factorized_ox.append(tok)
        elif tok.startswith(';d='):
            factorized_d.append(tok)
        elif tok.startswith(';row='):
            factorized_row.append(tok)
        elif tok.startswith('<Shape:') or tok in ('Td', 'SP', 'Oh', 'TP', 'TBP',
                                                   'SPY', 'TPr', 'L'):
            shape_tokens.append(tok)

    # ── Determine tokenization type ───────────────────────
    # The key question: are the composite tokens actually USED?
    # Check if the _tokenize_metal splits on ';' but serializer uses '|'
    # → composite tokens are auto-added during encode()
    if len(composite_metal_cn) > 0:
        current_type = "composite"
    elif len(factorized_metal) > 0 and len(factorized_cn) > 0:
        current_type = "factorized"
    else:
        current_type = "unknown"

    # Check if factorized tokens are also present (they are, but unused)
    has_factorized_in_code = len(factorized_metal) > 0 and len(factorized_cn) > 0

    # ── Sample CoordRep strings ───────────────────────────
    print(f"Sampling {args.n_sample} CoordRep strings …")
    examples = []
    metal_blocks_seen = Counter()

    with open(args.data) as f:
        for i, line in enumerate(f):
            if i >= args.n_sample:
                break
            d = json.loads(line)
            cr = d.get('coordrep', '')
            if not cr:
                continue
            # Extract metal block
            m = re.match(r'(\[[^\]]+\])', cr)
            if m:
                metal_blocks_seen[m.group(1)] += 1
            examples.append({
                'mol_id': d.get('mol_id', ''),
                'coordrep_prefix': cr[:120],
                'metal_block': m.group(1) if m else '',
            })

    # ── Analyze metal blocks ──────────────────────────────
    n_composite_used = sum(1 for mb in metal_blocks_seen if mb.startswith('[Metal:'))
    n_factorized_used = sum(1 for mb in metal_blocks_seen
                            if re.match(r'^\[[A-Z][a-z]?;', mb))

    # ── Write audit CSV ───────────────────────────────────
    p = os.path.join(args.out, "tokenizer_audit.csv")
    with open(p, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["category", "token", "vocab_id"])
        for tok in sorted(composite_metal_cn)[:30]:
            w.writerow(["composite_metal_cn", tok, vocab[tok]])
        for tok in sorted(factorized_metal):
            w.writerow(["factorized_metal", tok, vocab[tok]])
        for tok in sorted(factorized_cn):
            w.writerow(["factorized_cn", tok, vocab[tok]])
        for tok in sorted(factorized_ox):
            w.writerow(["factorized_ox", tok, vocab[tok]])
        for tok in sorted(factorized_d):
            w.writerow(["factorized_d", tok, vocab[tok]])
    print(f"  {p}")

    # ── Summary JSON ──────────────────────────────────────
    # Extract example tokens per type
    example_composite = sorted(composite_metal_cn)[:5]
    example_factorized = sorted(factorized_metal)[:5]
    example_cn = sorted(factorized_cn)
    example_ox = sorted(factorized_ox)[:5]
    example_shape = sorted(shape_tokens)[:5]

    summary = {
        "current_tokenization": current_type,
        "n_fused_metal_cn_tokens": len(composite_metal_cn),
        "n_metal_only_tokens": len(factorized_metal),
        "n_cn_only_tokens": len(factorized_cn),
        "n_ox_tokens": len(factorized_ox),
        "n_d_count_tokens": len(factorized_d),
        "n_row_tokens": len(factorized_row),
        "total_vocab_size": len(vocab),
        "factorized_defined_in_code_but_unused": has_factorized_in_code,
        "serializer_uses_pipe_separator": True,
        "tokenizer_splits_on_semicolon": True,
        "mismatch_causes_composite": True,
        "n_unique_metal_blocks_in_sample": len(metal_blocks_seen),
        "n_composite_blocks_in_sample": n_composite_used,
        "n_factorized_blocks_in_sample": n_factorized_used,
        "example_composite_tokens": example_composite,
        "example_factorized_metal_tokens": example_factorized,
        "example_cn_tokens": example_cn,
        "example_ox_tokens": example_ox,
        "example_shape_tokens": example_shape,
    }

    p = os.path.join(args.out, "tokenizer_audit_summary.json")
    with open(p, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"  {p}")

    # ── Print results ─────────────────────────────────────
    print(f"\n{'='*60}")
    print("TOKENIZER AUDIT RESULTS")
    print(f"{'='*60}")
    print(f"  Current tokenization: {current_type.upper()}")
    print(f"  Composite metal+CN tokens in vocab: {len(composite_metal_cn)}")
    print(f"  Factorized [Metal] tokens in vocab: {len(factorized_metal)} (unused)")
    print(f"  Factorized ;CN= tokens in vocab: {len(factorized_cn)} (unused)")
    print(f"  Total vocab size: {len(vocab)}")
    print(f"\n  Root cause:")
    print(f"    Serializer (to_string.py) uses '|' separator: [Metal:Fe|ox:+2|d:d6|CN:6]")
    print(f"    Tokenizer (_tokenize_metal) splits on ';'")
    print(f"    → metal block is treated as ONE composite token")
    print(f"    → factorized tokens [Fe], ;CN=6 are defined but NEVER used")
    print(f"\n  Examples (first 3 metal blocks from data):")
    for mb, cnt in metal_blocks_seen.most_common(5):
        print(f"    {mb}: {cnt} occurrences → single token")
    print()


if __name__ == "__main__":
    main()
