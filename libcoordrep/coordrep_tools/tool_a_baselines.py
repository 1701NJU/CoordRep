#!/usr/bin/env python3
"""
tool_a_baselines.py
===================

Three statistical baselines for Tool A donor-marker prediction:

A.  CondFreq(metal, CN)
    – Predict most frequent donor element for (metal, CN) pair.
    – Already exists in baselines.py; re-exported here with matching API.

B.  LigandFreq(SMILES)
    – Predict donor markers from ligand canonical SMILES alone.
    – Answers: "how much can you get from just knowing the ligand?"

C.  ContextFreq(metal, CN, shape_label | stereo_class)
    – Low-dimensional coordination-context lookup.
    – Answers: "does metal + geometry context add anything over
      just metal + CN?"
"""

from __future__ import annotations

import re
from collections import Counter, defaultdict
from typing import Dict, List, Optional, Tuple


# Standard donor atoms (same as baselines.py)
DONOR_ATOMS = ['C', 'N', 'O', 'S', 'P', 'F', 'Cl', 'Br', 'I',
               'Se', 'Te', 'As', 'Si', 'B', 'H']


# ──────────────────────────────────────────────────────────────────
# Helpers for extracting context from token lists or group dicts
# ──────────────────────────────────────────────────────────────────

_METAL_RE = re.compile(r'\[Metal:(\w+)\|.*?CN:(\d+)\]')
_SHAPE_BEST_RE = re.compile(r'^(Td|SP|Oh|TP|TBP|SPY|TPr|L)$')


def _extract_metal_cn(tokens: List[str]) -> Tuple[str, int]:
    for tok in tokens:
        m = _METAL_RE.match(tok)
        if m:
            return m.group(1), int(m.group(2))
    return "?", 0


def _extract_shape_label(tokens: List[str]) -> str:
    in_shape = False
    for tok in tokens:
        if tok == "<Shape:":
            in_shape = True
            continue
        if in_shape and _SHAPE_BEST_RE.match(tok):
            return tok
        if tok == ">":
            in_shape = False
    return "?"


def _extract_stereo_class(tokens: List[str]) -> str:
    has_trans = any(tok == "{trans:" for tok in tokens)
    has_cis = any(tok == "{cis:" for tok in tokens)
    if has_trans and has_cis:
        return "trans+cis"
    if has_trans:
        return "trans"
    if has_cis:
        return "cis"
    return "none"


def _extract_ligand_smiles(tokens: List[str], mask_pos: int) -> str:
    """
    Extract the canonical SMILES of the ligand block containing *mask_pos*.
    Returns the SMILES string (without donor markers) or "?" if not found.

    Handles both single-token Lid (``L1``) and split Lid (``L``, ``2``)
    tokenisation variants.
    """
    # Walk backwards to find = preceded by Lid
    lig_start = None
    for i in range(mask_pos, -1, -1):
        if tokens[i] == "=":
            # Case 1: single-token Lid  e.g. tokens[i-1] == "L1"
            if i > 0 and re.match(r'^L\d+$', tokens[i - 1]):
                lig_start = i + 1
                break
            # Case 2: split Lid  e.g. tokens[i-2] == "L", tokens[i-1] == "2"
            if (i > 1 and tokens[i - 2] == "L"
                    and tokens[i - 1].isdigit()):
                lig_start = i + 1
                break
    if lig_start is None:
        return "?"

    lig_end = mask_pos
    for j in range(mask_pos + 1, len(tokens)):
        if tokens[j] == "|":
            lig_end = j
            break

    smiles_toks = tokens[lig_start:lig_end]
    return "".join(smiles_toks)


# ──────────────────────────────────────────────────────────────────
# Baseline A: CondFreq(metal, CN)
# ──────────────────────────────────────────────────────────────────

class CondFreqBaseline:
    """
    Predict most frequent donor element for a (metal, CN) pair.
    Falls back to global frequency when the pair is unseen.
    """

    def __init__(self):
        self.cond_freq: Dict[Tuple, Counter] = defaultdict(Counter)
        self.global_freq: Counter = Counter()
        self._fitted = False

    def fit(self, samples: List[dict]):
        self.cond_freq = defaultdict(Counter)
        self.global_freq = Counter()
        for s in samples:
            donor = s.get("y_true", {}).get("donor_atom")
            if not donor:
                continue
            g = s.get("group", {})
            metal = g.get("metal_element", "?")
            cn = g.get("cn", 0)
            self.cond_freq[(metal, cn)][donor] += 1
            self.global_freq[donor] += 1
        self._fitted = True

    def predict(self, sample: dict, k: int = 5) -> List[Tuple[str, float]]:
        g = sample.get("group", {})
        metal = g.get("metal_element", "?")
        cn = g.get("cn", 0)
        freq = self.cond_freq.get((metal, cn), self.global_freq)
        total = sum(freq.values()) or 1
        ranked = freq.most_common(k)
        return [(tok, cnt / total) for tok, cnt in ranked]


# ──────────────────────────────────────────────────────────────────
# Baseline B: LigandFreq(SMILES)
# ──────────────────────────────────────────────────────────────────

class LigandFreqBaseline:
    """
    For each ligand SMILES seen in training, record the donor-atom
    distribution.  At test time, predict the most likely donor purely
    from the ligand SMILES.
    """

    def __init__(self):
        self.lig_freq: Dict[str, Counter] = defaultdict(Counter)
        self.global_freq: Counter = Counter()
        self._fitted = False

    def fit(self, samples: List[dict]):
        self.lig_freq = defaultdict(Counter)
        self.global_freq = Counter()
        for s in samples:
            donor = s.get("y_true", {}).get("donor_atom")
            if not donor:
                continue
            tokens = s.get("tokens", [])
            mask_pos = s.get("mask_pos")
            if isinstance(mask_pos, int):
                smiles = _extract_ligand_smiles(tokens, mask_pos)
            else:
                smiles = "?"
            self.lig_freq[smiles][donor] += 1
            self.global_freq[donor] += 1
        self._fitted = True

    def predict(self, sample: dict, k: int = 5) -> List[Tuple[str, float]]:
        tokens = sample.get("tokens", [])
        mask_pos = sample.get("mask_pos")
        if isinstance(mask_pos, int):
            smiles = _extract_ligand_smiles(tokens, mask_pos)
        else:
            smiles = "?"
        freq = self.lig_freq.get(smiles, self.global_freq)
        total = sum(freq.values()) or 1
        ranked = freq.most_common(k)
        return [(tok, cnt / total) for tok, cnt in ranked]

    def get_ligand_coverage(self, samples: List[dict]) -> float:
        """Fraction of test samples whose ligand SMILES was seen in training."""
        seen = 0
        for s in samples:
            tokens = s.get("tokens", [])
            mask_pos = s.get("mask_pos")
            if isinstance(mask_pos, int):
                smiles = _extract_ligand_smiles(tokens, mask_pos)
                if smiles in self.lig_freq:
                    seen += 1
        return seen / len(samples) if samples else 0


# ──────────────────────────────────────────────────────────────────
# Baseline C: ContextFreq(metal, CN, shape, stereo)
# ──────────────────────────────────────────────────────────────────

class ContextFreqBaseline:
    """
    Condition on (metal, CN, best_shape_label, stereo_class).
    Falls back to (metal, CN) then global.
    """

    def __init__(self):
        self.full_freq: Dict[Tuple, Counter] = defaultdict(Counter)
        self.metal_cn_freq: Dict[Tuple, Counter] = defaultdict(Counter)
        self.global_freq: Counter = Counter()
        self._fitted = False

    def fit(self, samples: List[dict]):
        self.full_freq = defaultdict(Counter)
        self.metal_cn_freq = defaultdict(Counter)
        self.global_freq = Counter()
        for s in samples:
            donor = s.get("y_true", {}).get("donor_atom")
            if not donor:
                continue
            g = s.get("group", {})
            tokens = s.get("tokens", [])
            metal = g.get("metal_element", "?")
            cn = g.get("cn", 0)
            shape = _extract_shape_label(tokens)
            stereo = _extract_stereo_class(tokens)
            self.full_freq[(metal, cn, shape, stereo)][donor] += 1
            self.metal_cn_freq[(metal, cn)][donor] += 1
            self.global_freq[donor] += 1
        self._fitted = True

    def predict(self, sample: dict, k: int = 5) -> List[Tuple[str, float]]:
        g = sample.get("group", {})
        tokens = sample.get("tokens", [])
        metal = g.get("metal_element", "?")
        cn = g.get("cn", 0)
        shape = _extract_shape_label(tokens)
        stereo = _extract_stereo_class(tokens)

        key_full = (metal, cn, shape, stereo)
        if key_full in self.full_freq:
            freq = self.full_freq[key_full]
        elif (metal, cn) in self.metal_cn_freq:
            freq = self.metal_cn_freq[(metal, cn)]
        else:
            freq = self.global_freq

        total = sum(freq.values()) or 1
        ranked = freq.most_common(k)
        return [(tok, cnt / total) for tok, cnt in ranked]


# ──────────────────────────────────────────────────────────────────
# Unified evaluation
# ──────────────────────────────────────────────────────────────────

def evaluate_baseline(
    baseline,
    test_samples: List[dict],
    k: int = 5,
) -> Dict:
    """
    Evaluate a baseline on test samples.

    Returns dict with top1_acc, top5_acc, mrr, n, per-sample predictions.
    """
    correct_top1 = 0
    correct_top5 = 0
    rr_sum = 0.0
    n = 0
    predictions = []

    for s in test_samples:
        donor = s.get("y_true", {}).get("donor_atom")
        if not donor:
            continue
        n += 1
        preds = baseline.predict(s, k=k)
        pred_tokens = [t for t, p in preds]

        is_top1 = len(pred_tokens) > 0 and pred_tokens[0] == donor
        is_top5 = donor in pred_tokens
        rank = (pred_tokens.index(donor) + 1) if donor in pred_tokens else 0

        if is_top1:
            correct_top1 += 1
        if is_top5:
            correct_top5 += 1
        if rank > 0:
            rr_sum += 1.0 / rank

        predictions.append({
            "id": s.get("id", ""),
            "target": donor,
            "top1_pred": pred_tokens[0] if pred_tokens else "",
            "top5_pred": pred_tokens[:5],
            "correct_top1": is_top1,
            "correct_top5": is_top5,
            "rank": rank,
        })

    return {
        "top1_acc": correct_top1 / n if n else 0,
        "top5_acc": correct_top5 / n if n else 0,
        "mrr": rr_sum / n if n else 0,
        "n": n,
        "predictions": predictions,
    }
