#!/usr/bin/env python3
"""
run_toolb_baselines.py
======================
Tool B sequence-repair baselines comparison.

Methods:
  A. corrupted_only   — no repair, measure valid parse rate
  B. rule_only        — deterministic bracket/delimiter/charge fixes
  C. edit_distance    — nearest-legal-token heuristic repair
  D. char_lm          — character-level bigram LM repair (no CoordRep tokenizer)
  E. coordrep_mlm     — CoordRep-specific MLM repair (Tool B)

Datasets:
  1. Original synthetic corruption suite
  2. CSD-derived corruption transfer set

Outputs:
  revision_results/toolb_sequence_baselines/toolb_baseline_summary.csv
  revision_results/toolb_sequence_baselines/toolb_by_corruption_type.csv
  revision_results/toolb_sequence_baselines/toolb_csd_transfer_summary.csv
  revision_results/toolb_sequence_baselines/toolb_failure_examples.csv
  revision_results/toolb_sequence_baselines/summary.json
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import random
import re
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

sys.path.insert(0, str(Path(__file__).parent.parent))

from coordrep_tools.validate import is_valid_coordrep, validate_brackets, validate_structure


# ══════════════════════════════════════════════════════════
# Field-level semantic parsing for comparison
# ══════════════════════════════════════════════════════════

def _parse_fields(s):
    """Extract key semantic fields from a CoordRep string."""
    fields = {}
    # Metal
    m = re.search(r'\[Metal:([A-Z][a-z]?)\|', s)
    fields['metal'] = m.group(1) if m else ''
    # CN
    m = re.search(r'CN:(\d+)', s)
    fields['cn'] = m.group(1) if m else ''
    # Oxidation state
    m = re.search(r'ox:([^|\]]+)', s)
    fields['ox'] = m.group(1) if m else ''
    # Ligand SMILES (sorted set)
    ligs = re.findall(r'\|L\d+=([^|]+)', s)
    fields['ligands'] = tuple(sorted(ligs))
    # Stereo relations
    stereo = re.findall(r'\{((?:trans|cis|mer|fac):[^}]+)\}', s)
    fields['stereo'] = tuple(sorted(stereo))
    # Shape
    m = re.search(r'<ShapeBest:([^|>]+)', s)
    fields['shape'] = m.group(1) if m else ''
    return fields


def _field_match_score(clean_str, repaired_str):
    """Compute fraction of semantic fields that match between clean and repaired.
    Returns (score, detail_dict)."""
    cf = _parse_fields(clean_str)
    rf = _parse_fields(repaired_str)
    checks = ['metal', 'cn', 'ox', 'ligands', 'stereo', 'shape']
    matches = sum(1 for k in checks if cf[k] == rf[k])
    return matches / len(checks)


# ══════════════════════════════════════════════════════════
# Corruption suite (mirrored from csd_toolb_transfer_eval)
# ══════════════════════════════════════════════════════════

def _apply_missing_bracket(chars):
    positions = [i for i, c in enumerate(chars) if c in '()[]{}']
    if not positions:
        return chars[:], False
    pos = random.choice(positions)
    out = chars[:]
    out.pop(pos)
    return out, True

def _apply_missing_charge(chars):
    positions = [i for i, c in enumerate(chars)
                 if c in '+-' and i + 1 < len(chars) and chars[i + 1].isdigit()]
    if not positions:
        return chars[:], False
    pos = random.choice(positions)
    out = chars[:]
    out.pop(pos)
    return out, True

def _apply_delimiter_mismatch(chars):
    dmap = {';': ',', ',': ';', '|': ';'}
    positions = [i for i, c in enumerate(chars) if c in dmap]
    if not positions:
        return chars[:], False
    pos = random.choice(positions)
    out = chars[:]
    out[pos] = dmap[out[pos]]
    return out, True

def _apply_truncation(chars):
    if len(chars) <= 15:
        return chars[:], False
    n_remove = random.randint(1, max(1, len(chars) // 5))
    return chars[:-n_remove], True

def _apply_token_swap(chars):
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


def generate_corrupted(coordrep_str, rng):
    types = list(CORRUPTION_PROBS.keys())
    weights = list(CORRUPTION_PROBS.values())
    ctype = rng.choices(types, weights=weights)[0]
    chars = list(coordrep_str)
    fn = CORRUPTION_FNS[ctype]
    corrupted, applied = fn(chars)
    if not applied:
        return coordrep_str, "none"
    return ''.join(corrupted), ctype


# ══════════════════════════════════════════════════════════
# Normalise corruption type names across datasets
# ══════════════════════════════════════════════════════════

_CTYPE_MAP = {
    'delimiter_swap': 'delimiter_mismatch',
    'bracket_swap': 'token_swap',
}

def _normalise_ctype(ct):
    return _CTYPE_MAP.get(ct, ct)


# ══════════════════════════════════════════════════════════
# Method A: Corrupted-only (no repair)
# ══════════════════════════════════════════════════════════

def repair_none(corrupted_str, clean_str=None):
    return corrupted_str


# ══════════════════════════════════════════════════════════
# Method B: Rule-only repair
# ══════════════════════════════════════════════════════════

def repair_rule_only(corrupted_str, clean_str=None):
    """Deterministic rule-based fixes: brackets, delimiters, charges."""
    s = corrupted_str

    # Fix 1: Balance brackets greedily
    s = _fix_brackets(s)

    # Fix 2: Fix common delimiter issues
    s = _fix_delimiters(s)

    # Fix 3: Fix missing charge sign before digit in Metal token
    s = _fix_charge(s)

    return s


def _fix_brackets(s):
    """Attempt to fix unbalanced brackets by insertion."""
    pairs = {'(': ')', '[': ']', '{': '}', '<': '>'}
    close_to_open = {v: k for k, v in pairs.items()}

    max_iters = 5
    for _ in range(max_iters):
        valid, err = validate_brackets(s)
        if valid:
            break

        # Find first error
        stack = []
        fixed = False
        for i, c in enumerate(s):
            if c in pairs:
                stack.append((c, i))
            elif c in close_to_open:
                if not stack:
                    # Unmatched closer → insert opener before it
                    opener = close_to_open[c]
                    s = s[:i] + opener + s[i:]
                    fixed = True
                    break
                top_c, _ = stack[-1]
                if pairs[top_c] != c:
                    # Mismatched → replace closer with expected
                    s = s[:i] + pairs[top_c] + s[i+1:]
                    stack.pop()
                    fixed = True
                    break
                else:
                    stack.pop()

        if not fixed and stack:
            # Unmatched opener → append closer at end
            top_c, _ = stack[-1]
            s = s + pairs[top_c]

    return s


def _fix_delimiters(s):
    """Fix obvious delimiter issues: ;; → ;, || in ligand dict → |."""
    s = re.sub(r';;+', ';', s)
    s = re.sub(r'\|\|+', '|', s)
    # Fix ';' where '|' expected in ligand dictionary positions
    # Pattern: |L\d+= ... ; L\d+= → should be | not ;
    s = re.sub(r';(L\d+=)', r'|\1', s)
    return s


def _fix_charge(s):
    """If metal token has ox field missing sign, insert +."""
    # Pattern: ox:2| → ox:+2|
    s = re.sub(r'ox:(\d)', r'ox:+\1', s)
    return s


# ══════════════════════════════════════════════════════════
# Method C: Edit-distance / heuristic repair
# ══════════════════════════════════════════════════════════

# Common legal tokens for nearest-match
_LEGAL_DELIMS = set('|;,')
_LEGAL_BRACKETS = set('()[]{}')
_BRACKET_PAIRS = {'(': ')', ')': '(', '[': ']', ']': '[',
                   '{': '}', '}': '{', '<': '>', '>': '<'}

def repair_edit_distance(corrupted_str, clean_str=None):
    """Heuristic: fix brackets by local context, fix delimiters by position."""
    s = corrupted_str

    # Step 1: Fix brackets (same greedy approach as rule, but also tries
    # swapping adjacent chars at error position)
    s = _fix_brackets(s)

    # Step 2: If still invalid, try swapping each adjacent pair of structural chars
    if not is_valid_coordrep(s):
        best = s
        for i in range(len(s) - 1):
            if s[i] in '()[]{}|;,' or s[i+1] in '()[]{}|;,':
                trial = list(s)
                trial[i], trial[i+1] = trial[i+1], trial[i]
                trial_s = ''.join(trial)
                if is_valid_coordrep(trial_s):
                    best = trial_s
                    break
        s = best

    # Step 3: If still invalid and looks truncated, try appending closing brackets
    if not is_valid_coordrep(s):
        # Count open brackets
        pairs = {'(': ')', '[': ']', '{': '}', '<': '>'}
        stack = []
        for c in s:
            if c in pairs:
                stack.append(c)
            elif c in pairs.values():
                if stack:
                    stack.pop()
        # Close any remaining open brackets
        while stack:
            s = s + pairs[stack.pop()]

    # Step 4: Fix delimiters
    s = _fix_delimiters(s)
    s = _fix_charge(s)

    return s


# ══════════════════════════════════════════════════════════
# Method D: Character-level bigram LM repair
# ══════════════════════════════════════════════════════════

class CharBigramRepairer:
    """Char bigram LM: learns P(c_i | c_{i-1}) from training data,
    then at each invalid position tries all printable chars and picks
    the one that maximises local likelihood while making the sequence valid."""

    def __init__(self):
        self.bigram_counts = defaultdict(lambda: defaultdict(int))
        self.unigram_counts = defaultdict(int)
        self._fitted = False

    def fit(self, train_strings):
        for s in train_strings:
            prev = '<S>'
            for c in s:
                self.bigram_counts[prev][c] += 1
                self.unigram_counts[c] += 1
                prev = c
            self.bigram_counts[prev]['<E>'] += 1
        self._fitted = True

    def _score_char(self, prev_char, c):
        total = sum(self.bigram_counts[prev_char].values())
        if total == 0:
            return 1e-10
        return (self.bigram_counts[prev_char].get(c, 0) + 1) / (total + 256)

    def repair(self, corrupted_str, clean_str=None):
        s = corrupted_str

        # Strategy: apply rule-based bracket fix first (char LM can't reason
        # about long-range bracket matching), then use bigram for local fixes
        s = _fix_brackets(s)
        s = _fix_delimiters(s)
        s = _fix_charge(s)

        if is_valid_coordrep(s):
            return s

        # Try single-char substitutions at structural positions
        chars_to_try = list('()[]{}|;,:<>+-=') + [c for c in set(s) if c.isalpha()]
        best = s
        for i in range(len(s)):
            if s[i] not in '()[]{}|;,:<>+-=':
                continue
            prev = s[i-1] if i > 0 else '<S>'
            best_score = self._score_char(prev, s[i])
            for c in chars_to_try:
                if c == s[i]:
                    continue
                score = self._score_char(prev, c)
                if score > best_score:
                    trial = s[:i] + c + s[i+1:]
                    if is_valid_coordrep(trial):
                        best = trial
                        best_score = score
                        break
            if is_valid_coordrep(best):
                break

        return best


# ══════════════════════════════════════════════════════════
# Method E: CoordRep-MLM repair (existing Tool B)
# ══════════════════════════════════════════════════════════

class CoordRepMLMRepairWrapper:
    """Wraps CoordRepRepairer for the same interface."""

    def __init__(self, checkpoint, tokenizer, device):
        from coordrep_tools.tool_b_repair import CoordRepRepairer
        self.repairer = CoordRepRepairer(checkpoint, tokenizer, device)

    def repair(self, corrupted_str, clean_str=None):
        result = self.repairer.repair(corrupted_str, max_iters=5, k=5)
        return result['repaired']


# ══════════════════════════════════════════════════════════
# Evaluation engine
# ══════════════════════════════════════════════════════════

def evaluate_method(method_name, repair_fn, samples):
    """Evaluate a repair method on a list of samples.

    Each sample: {coordrep_clean, coordrep_corrupted, corruption_type}
    Returns per-sample records and aggregate stats.
    """
    records = []
    agg = {
        'total': 0,
        'valid_before': 0,
        'valid_after': 0,
        'strict_valid_after': 0,
        'exact_match': 0,
        'field_match_sum': 0.0,
        'by_type': defaultdict(lambda: {
            'total': 0, 'valid_before': 0, 'valid_after': 0,
            'strict_valid_after': 0, 'exact_match': 0, 'field_match_sum': 0.0
        }),
    }

    for sample in samples:
        clean = sample['coordrep_clean']
        corrupted = sample['coordrep_corrupted']
        ctype = _normalise_ctype(sample.get('corruption_type', 'unknown'))

        vb = is_valid_coordrep(corrupted)
        repaired = repair_fn(corrupted, clean)
        va = is_valid_coordrep(repaired)
        sv = is_valid_coordrep(repaired, strict=True)
        exact = (repaired == clean)
        fm = _field_match_score(clean, repaired)

        agg['total'] += 1
        if vb: agg['valid_before'] += 1
        if va: agg['valid_after'] += 1
        if sv: agg['strict_valid_after'] += 1
        if exact: agg['exact_match'] += 1
        agg['field_match_sum'] += fm

        agg['by_type'][ctype]['total'] += 1
        if vb: agg['by_type'][ctype]['valid_before'] += 1
        if va: agg['by_type'][ctype]['valid_after'] += 1
        if sv: agg['by_type'][ctype]['strict_valid_after'] += 1
        if exact: agg['by_type'][ctype]['exact_match'] += 1
        agg['by_type'][ctype]['field_match_sum'] += fm

        records.append({
            'method': method_name,
            'corruption_type': ctype,
            'valid_before': vb,
            'valid_after': va,
            'strict_valid': sv,
            'exact_match': exact,
            'field_match': round(fm, 4),
            'clean_snippet': clean[:60],
            'corrupted_snippet': corrupted[:60],
            'repaired_snippet': repaired[:60],
        })

    n = agg['total']
    agg['valid_rate_before'] = agg['valid_before'] / max(n, 1)
    agg['valid_rate_after'] = agg['valid_after'] / max(n, 1)
    agg['strict_valid_rate'] = agg['strict_valid_after'] / max(n, 1)
    agg['exact_rate'] = agg['exact_match'] / max(n, 1)
    agg['field_match_rate'] = agg['field_match_sum'] / max(n, 1)
    agg['recovery_rate'] = ((agg['valid_after'] - agg['valid_before'])
                            / max(n - agg['valid_before'], 1))

    return records, agg


# ══════════════════════════════════════════════════════════
# Sampling helper
# ══════════════════════════════════════════════════════════

def _stratified_sample(samples, max_n, seed=42):
    """Stratified sample: equal number per corruption type."""
    by_type = defaultdict(list)
    for s in samples:
        ct = _normalise_ctype(s.get('corruption_type', 'unknown'))
        by_type[ct].append(s)
    n_types = len(by_type)
    per_type = max(1, max_n // max(n_types, 1))
    rng = random.Random(seed)
    out = []
    for ct in sorted(by_type.keys()):
        pool = by_type[ct]
        rng.shuffle(pool)
        out.extend(pool[:per_type])
    return out[:max_n]


# ══════════════════════════════════════════════════════════
# Data loading helpers
# ══════════════════════════════════════════════════════════

def load_synth_corrupted(path):
    samples = []
    with open(path) as f:
        for line in f:
            d = json.loads(line)
            samples.append({
                'coordrep_clean': d['coordrep_clean'],
                'coordrep_corrupted': d['coordrep_corrupted'],
                'corruption_type': d.get('corruption_type', 'unknown'),
                'id': d.get('id', ''),
            })
    return samples


def generate_csd_corrupted(retained_path, n=2000, seed=42):
    rng = random.Random(seed)
    entries = []
    with open(retained_path) as f:
        for line in f:
            entries.append(json.loads(line))

    sample_size = min(n, len(entries))
    sampled = rng.sample(entries, sample_size)

    samples = []
    for entry in sampled:
        cstr = entry['coordrep']
        corrupted, ctype = generate_corrupted(cstr, rng)
        if ctype == 'none':
            continue
        samples.append({
            'coordrep_clean': cstr,
            'coordrep_corrupted': corrupted,
            'corruption_type': ctype,
            'id': entry.get('refcode', ''),
        })
    return samples


# ══════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--synth_data',
                        default='/data/CoordRep/CoordSMILES/libcoordrep/outputs/fig5/synth_corrupted.jsonl')
    parser.add_argument('--csd_retained',
                        default='revision_results/csd_external/csd_retained_entries.jsonl')
    parser.add_argument('--checkpoint',
                        default='/data/CoordRep/CoordSMILES/libcoordrep/checkpoints/pretrain_v3/best_model.pt')
    parser.add_argument('--tokenizer',
                        default='/data/CoordRep/CoordSMILES/libcoordrep/checkpoints/pretrain_v3/tokenizer.json')
    parser.add_argument('--out', default='revision_results/toolb_sequence_baselines')
    parser.add_argument('--device', default='cuda:0')
    parser.add_argument('--n_csd', type=int, default=2000)
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--max_mlm', type=int, default=500,
                        help='Max samples for MLM repair (slow)')
    args = parser.parse_args()

    random.seed(args.seed)
    os.makedirs(args.out, exist_ok=True)

    print("="*70)
    print("Tool B Sequence Baselines Comparison")
    print("="*70)

    # ── Load data ─────────────────────────────────────────
    print("\n1. Loading data …")

    synth_samples = load_synth_corrupted(args.synth_data)
    print(f"   Synthetic corruption suite: {len(synth_samples)}")
    synth_types = Counter(_normalise_ctype(s['corruption_type']) for s in synth_samples)
    for t, c in synth_types.most_common():
        print(f"     {t}: {c}")

    csd_samples = generate_csd_corrupted(args.csd_retained, n=args.n_csd, seed=args.seed)
    print(f"   CSD transfer set: {len(csd_samples)}")
    csd_types = Counter(s['corruption_type'] for s in csd_samples)
    for t, c in csd_types.most_common():
        print(f"     {t}: {c}")

    # ── Fit char bigram from clean training data ──────────
    print("\n2. Fitting char bigram LM …")
    clean_strs = [s['coordrep_clean'] for s in synth_samples]
    char_lm = CharBigramRepairer()
    char_lm.fit(clean_strs)
    # Also fit on CSD clean strings
    csd_clean = [s['coordrep_clean'] for s in csd_samples]
    char_lm.fit(csd_clean)
    print(f"   Fitted on {len(clean_strs) + len(csd_clean)} strings")

    # ── Build methods ─────────────────────────────────────
    methods = {
        'corrupted_only': repair_none,
        'rule_only': repair_rule_only,
        'edit_distance': repair_edit_distance,
        'char_bigram': char_lm.repair,
    }

    # MLM wrapper (expensive — limit sample count)
    print("\n3. Loading CoordRep-MLM repairer …")
    mlm_wrapper = CoordRepMLMRepairWrapper(args.checkpoint, args.tokenizer, args.device)
    methods['coordrep_mlm'] = mlm_wrapper.repair

    # ── Evaluate on synthetic corruption ──────────────────
    print("\n4. Evaluating on synthetic corruption suite …")
    all_records_synth = []
    all_agg_synth = {}

    for mname in ['corrupted_only', 'rule_only', 'edit_distance', 'char_bigram']:
        t0 = time.time()
        recs, agg = evaluate_method(mname, methods[mname], synth_samples)
        elapsed = time.time() - t0
        all_records_synth.extend(recs)
        all_agg_synth[mname] = agg
        print(f"   {mname:<20s}: valid_after={agg['valid_rate_after']:.3f}  "
              f"exact={agg['exact_rate']:.3f}  recovery={agg['recovery_rate']:.3f}  "
              f"({elapsed:.1f}s)")

    # MLM on stratified subset (ensure all corruption types represented)
    mlm_subset = _stratified_sample(synth_samples, args.max_mlm)
    print(f"   CoordRep-MLM (on {len(mlm_subset)} stratified samples) …")
    t0 = time.time()
    recs_mlm, agg_mlm = evaluate_method('coordrep_mlm', methods['coordrep_mlm'], mlm_subset)
    elapsed = time.time() - t0
    all_records_synth.extend(recs_mlm)
    all_agg_synth['coordrep_mlm'] = agg_mlm
    print(f"   {'coordrep_mlm':<20s}: valid_after={agg_mlm['valid_rate_after']:.3f}  "
          f"exact={agg_mlm['exact_rate']:.3f}  recovery={agg_mlm['recovery_rate']:.3f}  "
          f"({elapsed:.1f}s)")

    # ── Evaluate on CSD transfer ──────────────────────────
    print("\n5. Evaluating on CSD transfer set …")
    all_records_csd = []
    all_agg_csd = {}

    for mname in ['corrupted_only', 'rule_only', 'edit_distance', 'char_bigram']:
        t0 = time.time()
        recs, agg = evaluate_method(mname, methods[mname], csd_samples)
        elapsed = time.time() - t0
        all_records_csd.extend(recs)
        all_agg_csd[mname] = agg
        print(f"   {mname:<20s}: valid_after={agg['valid_rate_after']:.3f}  "
              f"exact={agg['exact_rate']:.3f}  recovery={agg['recovery_rate']:.3f}  "
              f"({elapsed:.1f}s)")

    mlm_csd_subset = _stratified_sample(csd_samples, args.max_mlm)
    print(f"   CoordRep-MLM (on {len(mlm_csd_subset)} stratified samples) …")
    t0 = time.time()
    recs_csd_mlm, agg_csd_mlm = evaluate_method('coordrep_mlm', methods['coordrep_mlm'], mlm_csd_subset)
    elapsed = time.time() - t0
    all_records_csd.extend(recs_csd_mlm)
    all_agg_csd['coordrep_mlm'] = agg_csd_mlm
    print(f"   {'coordrep_mlm':<20s}: valid_after={agg_csd_mlm['valid_rate_after']:.3f}  "
          f"exact={agg_csd_mlm['exact_rate']:.3f}  recovery={agg_csd_mlm['recovery_rate']:.3f}  "
          f"({elapsed:.1f}s)")

    # ── Same-subset QC: re-evaluate all baselines on MLM's 500 ──
    print("\n5b. Same-subset QC (all methods on MLM's 500) …")
    same500_rows = []
    # Reuse already-computed MLM agg; only re-run cheap baselines on same subset
    mlm_agg_cache = {'synthetic': agg_mlm, 'csd_transfer': agg_csd_mlm}
    for ds_label, subset in [('synthetic', mlm_subset), ('csd_transfer', mlm_csd_subset)]:
        for mname in ['corrupted_only', 'rule_only', 'edit_distance', 'char_bigram', 'coordrep_mlm']:
            if mname == 'coordrep_mlm':
                agg = mlm_agg_cache[ds_label]
            else:
                _, agg = evaluate_method(mname, methods[mname], subset)
            same500_rows.append({
                'dataset': ds_label, 'method': mname, 'n': agg['total'],
                'valid_rate_before': round(agg['valid_rate_before'], 4),
                'valid_rate_after': round(agg['valid_rate_after'], 4),
                'strict_valid_rate': round(agg['strict_valid_rate'], 4),
                'exact_rate': round(agg['exact_rate'], 4),
                'field_match_rate': round(agg['field_match_rate'], 4),
                'recovery_rate': round(agg['recovery_rate'], 4),
            })
            print(f"   [{ds_label}] {mname:<20s}: valid={agg['valid_rate_after']:.3f} "
                  f"strict={agg['strict_valid_rate']:.3f} exact={agg['exact_rate']:.3f} "
                  f"field={agg['field_match_rate']:.3f}")

    p500 = os.path.join(args.out, 'toolb_same500_subset.csv')
    with open(p500, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(same500_rows[0].keys()))
        w.writeheader()
        for r in same500_rows:
            w.writerow(r)
    print(f"   → {p500}")

    # ══════════════════════════════════════════════════════
    # Write outputs
    # ══════════════════════════════════════════════════════
    print("\n6. Writing outputs …")

    # ── toolb_baseline_summary.csv ────────────────────────
    p = os.path.join(args.out, 'toolb_baseline_summary.csv')
    with open(p, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['dataset', 'method', 'n', 'valid_before', 'valid_after',
                     'strict_valid_after', 'exact_match',
                     'valid_rate_before', 'valid_rate_after',
                     'strict_valid_rate', 'exact_rate',
                     'field_match_rate', 'recovery_rate'])
        for mname in ['corrupted_only', 'rule_only', 'edit_distance', 'char_bigram', 'coordrep_mlm']:
            a = all_agg_synth.get(mname)
            if a:
                w.writerow(['synthetic', mname, a['total'],
                            a['valid_before'], a['valid_after'],
                            a['strict_valid_after'], a['exact_match'],
                            round(a['valid_rate_before'], 4),
                            round(a['valid_rate_after'], 4),
                            round(a['strict_valid_rate'], 4),
                            round(a['exact_rate'], 4),
                            round(a['field_match_rate'], 4),
                            round(a['recovery_rate'], 4)])
        for mname in ['corrupted_only', 'rule_only', 'edit_distance', 'char_bigram', 'coordrep_mlm']:
            a = all_agg_csd.get(mname)
            if a:
                w.writerow(['csd_transfer', mname, a['total'],
                            a['valid_before'], a['valid_after'],
                            a['strict_valid_after'], a['exact_match'],
                            round(a['valid_rate_before'], 4),
                            round(a['valid_rate_after'], 4),
                            round(a['strict_valid_rate'], 4),
                            round(a['exact_rate'], 4),
                            round(a['field_match_rate'], 4),
                            round(a['recovery_rate'], 4)])
    print(f"   {p}")

    # ── toolb_by_corruption_type.csv ──────────────────────
    p = os.path.join(args.out, 'toolb_by_corruption_type.csv')
    with open(p, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['dataset', 'method', 'corruption_type', 'n',
                     'valid_before', 'valid_after', 'strict_valid_after',
                     'exact_match', 'valid_rate_after', 'strict_valid_rate',
                     'exact_rate', 'field_match_rate'])
        for ds_name, agg_dict in [('synthetic', all_agg_synth), ('csd_transfer', all_agg_csd)]:
            for mname in ['corrupted_only', 'rule_only', 'edit_distance', 'char_bigram', 'coordrep_mlm']:
                a = agg_dict.get(mname)
                if not a:
                    continue
                for ctype in sorted(a['by_type'].keys()):
                    d = a['by_type'][ctype]
                    t = d['total']
                    w.writerow([ds_name, mname, ctype, t,
                                d['valid_before'], d['valid_after'],
                                d['strict_valid_after'], d['exact_match'],
                                round(d['valid_after'] / max(t, 1), 4),
                                round(d['strict_valid_after'] / max(t, 1), 4),
                                round(d['exact_match'] / max(t, 1), 4),
                                round(d['field_match_sum'] / max(t, 1), 4)])
    print(f"   {p}")

    # ── toolb_csd_transfer_summary.csv ────────────────────
    p = os.path.join(args.out, 'toolb_csd_transfer_summary.csv')
    with open(p, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['method', 'n', 'valid_rate_before', 'valid_rate_after',
                     'strict_valid_rate', 'exact_rate',
                     'field_match_rate', 'recovery_rate'])
        for mname in ['corrupted_only', 'rule_only', 'edit_distance', 'char_bigram', 'coordrep_mlm']:
            a = all_agg_csd.get(mname)
            if a:
                w.writerow([mname, a['total'],
                            round(a['valid_rate_before'], 4),
                            round(a['valid_rate_after'], 4),
                            round(a['strict_valid_rate'], 4),
                            round(a['exact_rate'], 4),
                            round(a['field_match_rate'], 4),
                            round(a['recovery_rate'], 4)])
    print(f"   {p}")

    # ── toolb_failure_examples.csv ────────────────────────
    p = os.path.join(args.out, 'toolb_failure_examples.csv')
    failures = []
    for rec in all_records_synth:
        if rec['method'] == 'coordrep_mlm' and not rec['valid_after'] and not rec['valid_before']:
            failures.append(rec)
    for rec in all_records_csd:
        if rec['method'] == 'coordrep_mlm' and not rec['valid_after'] and not rec['valid_before']:
            failures.append(rec)
    with open(p, 'w', newline='') as f:
        if failures:
            w = csv.DictWriter(f, fieldnames=list(failures[0].keys()))
            w.writeheader()
            for rec in failures[:200]:
                w.writerow(rec)
    print(f"   {p} ({len(failures)} failure examples)")

    # ── summary.json ──────────────────────────────────────
    summary = {
        'synthetic': {},
        'csd_transfer': {},
    }
    for mname in ['corrupted_only', 'rule_only', 'edit_distance', 'char_bigram', 'coordrep_mlm']:
        for ds_name, agg_dict in [('synthetic', all_agg_synth), ('csd_transfer', all_agg_csd)]:
            a = agg_dict.get(mname)
            if not a:
                continue
            entry = {
                'n': a['total'],
                'valid_rate_before': round(a['valid_rate_before'], 4),
                'valid_rate_after': round(a['valid_rate_after'], 4),
                'strict_valid_rate': round(a['strict_valid_rate'], 4),
                'exact_rate': round(a['exact_rate'], 4),
                'field_match_rate': round(a['field_match_rate'], 4),
                'recovery_rate': round(a['recovery_rate'], 4),
                'by_type': {},
            }
            for ctype in sorted(a['by_type'].keys()):
                d = a['by_type'][ctype]
                t = d['total']
                entry['by_type'][ctype] = {
                    'n': t,
                    'valid_rate_after': round(d['valid_after'] / max(t, 1), 4),
                    'strict_valid_rate': round(d['strict_valid_after'] / max(t, 1), 4),
                    'exact_rate': round(d['exact_match'] / max(t, 1), 4),
                    'field_match_rate': round(d['field_match_sum'] / max(t, 1), 4),
                }
            summary[ds_name][mname] = entry

    p = os.path.join(args.out, 'summary.json')
    with open(p, 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"   {p}")

    # ── Print summary table ───────────────────────────────
    print(f"\n{'='*70}")
    print("TOOL B BASELINE SUMMARY")
    print(f"{'='*70}")
    for ds_name, agg_dict in [('Synthetic', all_agg_synth), ('CSD Transfer', all_agg_csd)]:
        print(f"\n  {ds_name}:")
        print(f"  {'Method':<20s} {'N':>6s} {'ValidBef':>9s} {'ValidAft':>9s} "
              f"{'StrictV':>8s} {'Exact':>7s} {'FieldM':>7s} {'Recovery':>9s}")
        print(f"  {'-'*80}")
        for mname in ['corrupted_only', 'rule_only', 'edit_distance', 'char_bigram', 'coordrep_mlm']:
            a = agg_dict.get(mname)
            if not a:
                continue
            print(f"  {mname:<20s} {a['total']:>6d} {a['valid_rate_before']:>9.3f} "
                  f"{a['valid_rate_after']:>9.3f} {a['strict_valid_rate']:>8.3f} "
                  f"{a['exact_rate']:>7.3f} {a['field_match_rate']:>7.3f} "
                  f"{a['recovery_rate']:>9.3f}")

    print(f"\n{'='*70}")
    print("Done.")
    print(f"{'='*70}")


if __name__ == '__main__':
    main()
