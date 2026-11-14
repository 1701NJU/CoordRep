"""
ranker_dataset.py
=================
Datasets for CoordRep-Ranker training and evaluation.

Provides listwise examples: 1 real + K decoys per sample.
Decoy generation is done on-the-fly using HardDecoyPool.
Supports decoy-type labels for stratified evaluation.
"""

from __future__ import annotations

import json
import random
import re
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

import torch
from torch.utils.data import Dataset

from .hard_decoys import HardDecoyPool, generate_all_hard_decoys, HARD_DECOY_TYPES
from .ablation_masking import _detect_blocks


# ── helpers ───────────────────────────────────────────────────────

def _tokens_to_ids(tokens: List[str], tokenizer) -> List[int]:
    unk = tokenizer.token2id.get("[UNK]", 1)
    return [tokenizer.token2id.get(t, unk) for t in tokens]


def _extract_unique_complexes(path: str) -> List[dict]:
    """Load unique complexes from fig5_tasks.jsonl."""
    seen = {}
    with open(path) as f:
        for line in f:
            rec = json.loads(line)
            m = re.match(r'^(tmqm_\d+)', rec["id"])
            if not m:
                continue
            cid = m.group(1)
            if cid not in seen:
                seen[cid] = {"id": cid, "tokens": rec["tokens"]}
    return list(seen.values())


def split_complexes(
    complexes: List[dict],
    train_frac: float = 0.7,
    val_frac: float = 0.15,
    seed: int = 42,
) -> Tuple[List[dict], List[dict], List[dict]]:
    """Deterministic train/val/test split by complex ID."""
    rng = random.Random(seed)
    ids_sorted = sorted(complexes, key=lambda c: c["id"])
    rng.shuffle(ids_sorted)
    n = len(ids_sorted)
    n_train = int(n * train_frac)
    n_val = int(n * val_frac)
    return ids_sorted[:n_train], ids_sorted[n_train:n_train+n_val], ids_sorted[n_train+n_val:]


# ── listwise dataset ─────────────────────────────────────────────

class RankerListwiseDataset(Dataset):
    """
    Each __getitem__ returns:
      input_ids:      (1+K, max_len)   — [real, decoy_1, …, decoy_K]
      attention_mask:  (1+K, max_len)
      decoy_types:    list of K strings

    Decoys are generated on-the-fly from the pool,
    so the pool must be built from the TRAIN split only.
    """

    def __init__(
        self,
        complexes: List[dict],
        pool: HardDecoyPool,
        tokenizer,
        n_decoys: int = 20,
        max_len: int = 512,
        seed: int = 42,
    ):
        self.complexes = complexes
        self.pool = pool
        self.tokenizer = tokenizer
        self.n_decoys = n_decoys
        self.max_len = max_len
        self.pad_id = tokenizer.token2id.get("[PAD]", 0)
        self.rng = random.Random(seed)

    def __len__(self):
        return len(self.complexes)

    def _pad(self, ids: List[int]) -> Tuple[List[int], List[int]]:
        ids = ids[:self.max_len]
        att = [1] * len(ids) + [0] * (self.max_len - len(ids))
        ids = ids + [self.pad_id] * (self.max_len - len(ids))
        return ids, att

    def __getitem__(self, idx: int) -> dict:
        c = self.complexes[idx]
        toks = c["tokens"]
        real_ids = _tokens_to_ids(toks, self.tokenizer)

        # generate hard decoys
        decoys = generate_all_hard_decoys(toks, self.pool, n_per_type=max(1, self.n_decoys // 5))
        if len(decoys) > self.n_decoys:
            self.rng.shuffle(decoys)
            decoys = decoys[:self.n_decoys]

        # pad to exactly n_decoys (repeat if needed)
        while len(decoys) < self.n_decoys:
            if decoys:
                decoys.append(decoys[self.rng.randint(0, len(decoys) - 1)])
            else:
                # fallback: use real as decoy (loss = 0)
                decoys.append((list(toks), "fallback"))

        all_ids = []
        all_att = []
        decoy_types = []

        # real first
        r_ids, r_att = self._pad(real_ids)
        all_ids.append(r_ids)
        all_att.append(r_att)

        for dtoks, reason in decoys:
            d_ids = _tokens_to_ids(dtoks, self.tokenizer)
            d_padded, d_att = self._pad(d_ids)
            all_ids.append(d_padded)
            all_att.append(d_att)
            decoy_types.append(reason.split(":")[0] if reason else "unknown")

        return {
            "input_ids": torch.tensor(all_ids, dtype=torch.long),
            "attention_mask": torch.tensor(all_att, dtype=torch.long),
            "decoy_types": decoy_types,
            "complex_id": c["id"],
        }


# ── eval dataset (pre-generated decoys) ──────────────────────────

class RankerEvalDataset(Dataset):
    """
    Pre-generates all decoys once for deterministic evaluation.
    Returns same format as RankerListwiseDataset.
    """

    def __init__(
        self,
        complexes: List[dict],
        pool: HardDecoyPool,
        tokenizer,
        n_decoys: int = 20,
        max_len: int = 512,
        seed: int = 42,
    ):
        self.tokenizer = tokenizer
        self.max_len = max_len
        self.pad_id = tokenizer.token2id.get("[PAD]", 0)
        self.samples = []

        rng = random.Random(seed)
        for c in complexes:
            toks = c["tokens"]
            decoys = generate_all_hard_decoys(toks, pool, n_per_type=max(1, n_decoys // 5))
            if len(decoys) > n_decoys:
                rng.shuffle(decoys)
                decoys = decoys[:n_decoys]
            while len(decoys) < n_decoys and decoys:
                decoys.append(decoys[rng.randint(0, len(decoys) - 1)])
            if not decoys:
                continue
            self.samples.append({
                "id": c["id"],
                "tokens": toks,
                "decoys": decoys,
            })

    def __len__(self):
        return len(self.samples)

    def _pad(self, ids: List[int]) -> Tuple[List[int], List[int]]:
        ids = ids[:self.max_len]
        att = [1] * len(ids) + [0] * (self.max_len - len(ids))
        ids = ids + [self.pad_id] * (self.max_len - len(ids))
        return ids, att

    def __getitem__(self, idx: int) -> dict:
        s = self.samples[idx]
        real_ids = _tokens_to_ids(s["tokens"], self.tokenizer)
        all_ids, all_att, dtypes = [], [], []

        r_ids, r_att = self._pad(real_ids)
        all_ids.append(r_ids)
        all_att.append(r_att)

        for dtoks, reason in s["decoys"]:
            d_ids = _tokens_to_ids(dtoks, self.tokenizer)
            d_padded, d_att = self._pad(d_ids)
            all_ids.append(d_padded)
            all_att.append(d_att)
            dtypes.append(reason.split(":")[0] if reason else "unknown")

        return {
            "input_ids": torch.tensor(all_ids, dtype=torch.long),
            "attention_mask": torch.tensor(all_att, dtype=torch.long),
            "decoy_types": dtypes,
            "complex_id": s["id"],
        }
