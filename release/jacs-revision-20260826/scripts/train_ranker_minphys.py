#!/usr/bin/env python3
"""
train_ranker_minphys.py
=======================
Train three MinPhys diagnostic variants of CoordRep-Ranker:

  1. minphys_only:          only MinPhys aux tokens  (no CoordRep)
  2. ranker_minphys:        CoordRep tokens + MinPhys aux tokens
  3. ranker_shuffled_minphys: CoordRep tokens + batch-shuffled MinPhys tokens

Uses the same train/val/test split and decoy pool as the base ranker.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import random
import sys
import time
from pathlib import Path
from typing import List, Tuple

import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from coordrep_tools.coordrep_ranker import (
    load_ranker_from_mlm, save_ranker, CoordRepRanker,
)
from coordrep_tools.ranker_dataset import (
    _extract_unique_complexes, split_complexes, _tokens_to_ids,
)
from coordrep_tools.ranker_losses import get_loss_fn
from coordrep_tools.hard_decoys import HardDecoyPool, generate_all_hard_decoys
from coordrep_tools.minphys_features import MinPhysBinner
from coordrep_tools.minphys_tokens import (
    make_minphys_tokens, append_minphys, register_minphys_tokens_in_tokenizer,
)
from brain.model import CoordRepModelConfig
from brain.tokenizer import CoordRepTokenizer, TokenizerConfig


# ── dataset ───────────────────────────────────────────────────────

class MinPhysListwiseDataset(Dataset):
    """
    Listwise dataset with MinPhys token modes:
      mode="minphys_only"   → only Aux tokens (no CoordRep)
      mode="coordrep+minphys" → CoordRep + Aux tokens appended
      mode="shuffled_minphys" → CoordRep + Aux tokens from random other sample
    """

    def __init__(
        self,
        complexes: List[dict],
        pool: HardDecoyPool,
        tokenizer,
        binner: MinPhysBinner,
        mode: str = "coordrep+minphys",
        n_decoys: int = 20,
        max_len: int = 512,
        seed: int = 42,
    ):
        self.complexes = complexes
        self.pool = pool
        self.tokenizer = tokenizer
        self.binner = binner
        self.mode = mode
        self.n_decoys = n_decoys
        self.max_len = max_len
        self.pad_id = tokenizer.token2id.get("[PAD]", 0)
        self.rng = random.Random(seed)

        # Pre-compute MinPhys tokens for all complexes (for shuffled mode)
        if mode == "shuffled_minphys":
            self._all_aux = [make_minphys_tokens(c["tokens"], binner) for c in complexes]

    def __len__(self):
        return len(self.complexes)

    def _tokenise(self, toks: List[str], shuffle_aux: List[str] = None) -> List[int]:
        """Convert tokens to IDs based on mode."""
        if self.mode == "minphys_only":
            aux = make_minphys_tokens(toks, self.binner)
            return _tokens_to_ids(aux, self.tokenizer)
        elif self.mode == "shuffled_minphys":
            aux = shuffle_aux if shuffle_aux else make_minphys_tokens(toks, self.binner)
            combined = list(toks) + aux
            return _tokens_to_ids(combined, self.tokenizer)
        else:  # coordrep+minphys
            aux = make_minphys_tokens(toks, self.binner)
            combined = list(toks) + aux
            return _tokens_to_ids(combined, self.tokenizer)

    def _pad(self, ids: List[int]) -> Tuple[List[int], List[int]]:
        ids = ids[:self.max_len]
        att = [1] * len(ids) + [0] * (self.max_len - len(ids))
        ids = ids + [self.pad_id] * (self.max_len - len(ids))
        return ids, att

    def __getitem__(self, idx: int) -> dict:
        c = self.complexes[idx]
        toks = c["tokens"]

        # generate hard decoys
        decoys = generate_all_hard_decoys(toks, self.pool, n_per_type=max(1, self.n_decoys // 5))
        if len(decoys) > self.n_decoys:
            self.rng.shuffle(decoys)
            decoys = decoys[:self.n_decoys]
        while len(decoys) < self.n_decoys:
            if decoys:
                decoys.append(decoys[self.rng.randint(0, len(decoys) - 1)])
            else:
                decoys.append((list(toks), "fallback"))

        # For shuffled mode: pick a random other complex's aux tokens
        shuffle_aux = None
        if self.mode == "shuffled_minphys":
            other_idx = self.rng.randint(0, len(self.complexes) - 1)
            while other_idx == idx and len(self.complexes) > 1:
                other_idx = self.rng.randint(0, len(self.complexes) - 1)
            shuffle_aux = self._all_aux[other_idx]

        all_ids, all_att, dtypes = [], [], []

        # real
        r_ids = self._tokenise(toks, shuffle_aux)
        rp, ra = self._pad(r_ids)
        all_ids.append(rp)
        all_att.append(ra)

        for dtoks, reason in decoys:
            # IMPORTANT: for each decoy, recompute MinPhys tokens from DECOY tokens
            if self.mode == "shuffled_minphys":
                d_ids = self._tokenise(dtoks, shuffle_aux)
            else:
                d_ids = self._tokenise(dtoks)
            dp, da = self._pad(d_ids)
            all_ids.append(dp)
            all_att.append(da)
            dtypes.append(reason.split(":")[0] if reason else "unknown")

        return {
            "input_ids": torch.tensor(all_ids, dtype=torch.long),
            "attention_mask": torch.tensor(all_att, dtype=torch.long),
            "decoy_types": dtypes,
            "complex_id": c["id"],
        }


def listwise_collate(batch):
    return {
        "input_ids": torch.stack([b["input_ids"] for b in batch]),
        "attention_mask": torch.stack([b["attention_mask"] for b in batch]),
        "decoy_types": [b["decoy_types"] for b in batch],
        "complex_id": [b["complex_id"] for b in batch],
    }


# ── eval helper ───────────────────────────────────────────────────

@torch.no_grad()
def evaluate(model, loader, loss_fn, device):
    model.eval()
    total_loss, total_top1, total_mrr, n = 0.0, 0, 0.0, 0
    for batch in loader:
        ids = batch["input_ids"].to(device)
        att = batch["attention_mask"].to(device)
        scores = model.score_list(ids, att)
        loss = loss_fn(scores)
        total_loss += loss.item() * ids.size(0)
        ranks = (scores[:, 1:] >= scores[:, :1]).sum(dim=1) + 1
        total_top1 += (ranks == 1).sum().item()
        total_mrr += (1.0 / ranks.float()).sum().item()
        n += ids.size(0)
    model.train()
    return {
        "loss": total_loss / max(n, 1),
        "top1": total_top1 / max(n, 1),
        "mrr": total_mrr / max(n, 1),
        "n": n,
    }


# ── train one variant ─────────────────────────────────────────────

def train_variant(
    variant_name: str,
    mode: str,
    train_cx, val_cx, pool, tokenizer, binner,
    args, device,
):
    print(f"\n{'='*60}")
    print(f"Training variant: {variant_name} (mode={mode})")
    print(f"{'='*60}")

    # model
    ranker, _ = load_ranker_from_mlm(
        args.mlm_checkpoint, device=device,
        pooling=args.pooling, freeze_encoder=False,
    )
    print(f"  Trainable: {ranker.trainable_parameters():,}")

    # datasets
    train_ds = MinPhysListwiseDataset(
        train_cx, pool, tokenizer, binner,
        mode=mode, n_decoys=args.decoys_per_real,
        max_len=args.max_len, seed=args.seed,
    )
    val_ds = MinPhysListwiseDataset(
        val_cx, pool, tokenizer, binner,
        mode=mode, n_decoys=args.decoys_per_real,
        max_len=args.max_len, seed=args.seed + 1,
    )

    train_loader = DataLoader(
        train_ds, batch_size=args.batch_size, shuffle=True,
        collate_fn=listwise_collate, num_workers=args.num_workers,
        pin_memory=True, drop_last=True,
    )
    val_loader = DataLoader(
        val_ds, batch_size=args.batch_size * 2, shuffle=False,
        collate_fn=listwise_collate, num_workers=args.num_workers,
        pin_memory=True,
    )

    loss_fn = get_loss_fn("listwise").to(device)
    optimizer = torch.optim.AdamW([
        {"params": ranker.head.parameters(), "lr": args.lr},
        {"params": ranker.encoder.parameters(), "lr": args.encoder_lr},
    ], weight_decay=args.weight_decay)

    total_steps = len(train_loader) * args.epochs
    warmup = args.warmup_steps

    def lr_schedule(step):
        if step < warmup:
            return step / max(warmup, 1)
        progress = (step - warmup) / max(total_steps - warmup, 1)
        return 0.5 * (1 + math.cos(math.pi * progress))

    scheduler = torch.optim.lr_scheduler.LambdaLR(optimizer, lr_schedule)

    best_val_mrr = 0.0
    patience_counter = 0
    history = []

    for epoch in range(1, args.epochs + 1):
        ranker.train()
        epoch_loss, epoch_top1, epoch_n = 0.0, 0, 0
        t0 = time.time()

        for step, batch in enumerate(train_loader):
            ids = batch["input_ids"].to(device)
            att = batch["attention_mask"].to(device)
            scores = ranker.score_list(ids, att)
            loss = loss_fn(scores)

            optimizer.zero_grad()
            loss.backward()
            nn.utils.clip_grad_norm_(ranker.parameters(), args.grad_clip)
            optimizer.step()
            scheduler.step()

            B = ids.size(0)
            epoch_loss += loss.item() * B
            ranks = (scores[:, 1:] >= scores[:, :1]).sum(dim=1) + 1
            epoch_top1 += (ranks == 1).sum().item()
            epoch_n += B

            if step % 100 == 0:
                print(f"  E{epoch} step {step}/{len(train_loader)} "
                      f"loss={loss.item():.4f}", flush=True)

        train_loss = epoch_loss / max(epoch_n, 1)
        train_top1 = epoch_top1 / max(epoch_n, 1)
        val_m = evaluate(ranker, val_loader, loss_fn, device)
        elapsed = time.time() - t0

        rec = {"epoch": epoch,
               "train_loss": round(train_loss, 4), "train_top1": round(train_top1, 4),
               "val_loss": round(val_m["loss"], 4), "val_top1": round(val_m["top1"], 4),
               "val_mrr": round(val_m["mrr"], 4), "time": round(elapsed, 1)}
        history.append(rec)
        print(f"  E{epoch}: train_loss={train_loss:.4f} train_top1={train_top1:.1%} | "
              f"val_loss={val_m['loss']:.4f} val_top1={val_m['top1']:.1%} "
              f"val_mrr={val_m['mrr']:.3f} ({elapsed:.0f}s)")

        if val_m["mrr"] > best_val_mrr:
            best_val_mrr = val_m["mrr"]
            patience_counter = 0
            ckpt_path = os.path.join(args.out, f"best_{variant_name}.pt")
            save_ranker(ranker, ckpt_path, extra={"epoch": epoch, "val_mrr": best_val_mrr})
            print(f"  ★ New best val_mrr={best_val_mrr:.4f}")
        else:
            patience_counter += 1
            if patience_counter >= args.patience:
                print(f"  Early stopping at epoch {epoch}")
                break

    hist_path = os.path.join(args.out, f"history_{variant_name}.json")
    with open(hist_path, "w") as f:
        json.dump(history, f, indent=2)
    print(f"  Best val_mrr={best_val_mrr:.4f}")

    return best_val_mrr


# ── main ──────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mlm_checkpoint",
                        default="inputs/checkpoints/pretrain_v3/best_model.pt")
    parser.add_argument("--data",
                        default="inputs/property_benchmarks/fig5_tasks.jsonl")
    parser.add_argument("--split_ids",
                        default="inputs/checkpoints/coordrep_ranker/split_ids.json")
    parser.add_argument("--out",
                        default="outputs/checkpoints/ranker_minphys")
    parser.add_argument("--pooling", default="mean")
    parser.add_argument("--decoys_per_real", type=int, default=20)
    parser.add_argument("--epochs", type=int, default=10)
    parser.add_argument("--batch_size", type=int, default=4)
    parser.add_argument("--lr", type=float, default=5e-4)
    parser.add_argument("--encoder_lr", type=float, default=2e-5)
    parser.add_argument("--weight_decay", type=float, default=1e-4)
    parser.add_argument("--warmup_steps", type=int, default=300)
    parser.add_argument("--grad_clip", type=float, default=1.0)
    parser.add_argument("--patience", type=int, default=4)
    parser.add_argument("--max_len", type=int, default=512)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--device", default="cuda:0")
    parser.add_argument("--num_workers", type=int, default=4)
    # Select which variants to train
    parser.add_argument("--variants", nargs="+",
                        default=["minphys_only", "ranker_minphys", "shuffled_minphys"])
    args = parser.parse_args()

    os.makedirs(args.out, exist_ok=True)
    torch.manual_seed(args.seed)
    np.random.seed(args.seed)

    # ── 1. Load data & split ──────────────────────────────────
    print("Loading data …")
    complexes = _extract_unique_complexes(args.data)
    with open(args.split_ids) as f:
        split = json.load(f)
    id_to_cx = {c["id"]: c for c in complexes}
    train_cx = [id_to_cx[i] for i in split["train"] if i in id_to_cx]
    val_cx = [id_to_cx[i] for i in split["val"] if i in id_to_cx]
    test_cx = [id_to_cx[i] for i in split["test"] if i in id_to_cx]
    print(f"  train={len(train_cx)}, val={len(val_cx)}, test={len(test_cx)}")

    # ── 2. Build hard decoy pool (train only) ─────────────────
    print("Building hard decoy pool …")
    t0 = time.time()
    train_toks = [c["tokens"] for c in train_cx]
    pool = HardDecoyPool(train_toks, seed=args.seed)
    print(f"  ({time.time()-t0:.1f}s)")

    # ── 3. Fit MinPhys binner on train ────────────────────────
    print("Fitting MinPhys binner …")
    binner = MinPhysBinner()
    binner.fit(train_toks)
    print(f"  HA quantiles: q33={binner.ha_q33}, q66={binner.ha_q66}")

    # ── 4. Build tokenizer with aux tokens ────────────────────
    print("Loading tokenizer + registering MinPhys tokens …")
    tok_path = Path(args.mlm_checkpoint).parent / "tokenizer.json"
    with open(tok_path) as f:
        tok_data = json.load(f)
    tok_config = TokenizerConfig(**tok_data.get("config", {}))
    tokenizer = CoordRepTokenizer(tok_config)
    tokenizer.token2id = tok_data["token2id"]
    tokenizer.id2token = {v: k for k, v in tokenizer.token2id.items()}
    tokenizer.next_id = max(tokenizer.id2token.keys()) + 1

    n_new = register_minphys_tokens_in_tokenizer(tokenizer, train_toks, binner)
    print(f"  Registered {n_new} MinPhys token types, vocab now {tokenizer.next_id}")

    # Save binner for eval
    binner_info = {"ha_q33": binner.ha_q33, "ha_q66": binner.ha_q66}
    with open(os.path.join(args.out, "binner.json"), "w") as f:
        json.dump(binner_info, f)

    # ── 5. Train variants ─────────────────────────────────────
    variant_modes = {
        "minphys_only": "minphys_only",
        "ranker_minphys": "coordrep+minphys",
        "shuffled_minphys": "shuffled_minphys",
    }

    for vname in args.variants:
        mode = variant_modes[vname]
        train_variant(vname, mode, train_cx, val_cx, pool, tokenizer, binner, args, args.device)

    print(f"\nAll variants trained. Checkpoints in {args.out}")


if __name__ == "__main__":
    main()
