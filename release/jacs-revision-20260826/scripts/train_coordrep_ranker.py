#!/usr/bin/env python3
"""
train_coordrep_ranker.py
========================
Train CoordRep-Ranker: contrastive compatibility model.

Initialises encoder from pretrained MLM checkpoint, trains a scoring
head with listwise cross-entropy loss on hard-negative decoys.

Usage:
  python scripts/train_coordrep_ranker.py --freeze_encoder
  python scripts/train_coordrep_ranker.py --no-freeze_encoder --lr 5e-5
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import DataLoader

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from coordrep_tools.coordrep_ranker import (
    load_ranker_from_mlm, save_ranker,
)
from coordrep_tools.ranker_dataset import (
    _extract_unique_complexes, split_complexes,
    RankerListwiseDataset, RankerEvalDataset,
)
from coordrep_tools.ranker_losses import get_loss_fn
from coordrep_tools.hard_decoys import HardDecoyPool


# ── collate for listwise ─────────────────────────────────────────

def listwise_collate(batch):
    """Stack listwise samples into a batch."""
    input_ids = torch.stack([b["input_ids"] for b in batch])       # (B, 1+K, L)
    attention_mask = torch.stack([b["attention_mask"] for b in batch])
    return {
        "input_ids": input_ids,
        "attention_mask": attention_mask,
        "decoy_types": [b["decoy_types"] for b in batch],
        "complex_id": [b["complex_id"] for b in batch],
    }


# ── eval helper ───────────────────────────────────────────────────

@torch.no_grad()
def evaluate(model, loader, loss_fn, device):
    model.eval()
    total_loss, total_top1, total_mrr, n = 0.0, 0, 0.0, 0
    for batch in loader:
        ids = batch["input_ids"].to(device)      # (B, 1+K, L)
        att = batch["attention_mask"].to(device)
        scores = model.score_list(ids, att)       # (B, 1+K)
        loss = loss_fn(scores)
        total_loss += loss.item() * ids.size(0)

        ranks = (scores[:, 1:] >= scores[:, :1]).sum(dim=1) + 1  # rank of real
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


# ── main ──────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mlm_checkpoint",
                        default="inputs/checkpoints/pretrain_v3/best_model.pt")
    parser.add_argument("--data",
                        default="inputs/property_benchmarks/fig5_tasks.jsonl")
    parser.add_argument("--out",
                        default="outputs/checkpoints/coordrep_ranker")
    parser.add_argument("--loss", choices=["margin", "bce", "listwise"], default="listwise")
    parser.add_argument("--freeze_encoder", action="store_true")
    parser.add_argument("--no-freeze_encoder", dest="freeze_encoder", action="store_false")
    parser.set_defaults(freeze_encoder=True)
    parser.add_argument("--pooling", default="mean", choices=["mean", "cls"])
    parser.add_argument("--decoys_per_real", type=int, default=20)
    parser.add_argument("--epochs", type=int, default=10)
    parser.add_argument("--batch_size", type=int, default=16)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--encoder_lr", type=float, default=5e-5)
    parser.add_argument("--weight_decay", type=float, default=1e-4)
    parser.add_argument("--warmup_steps", type=int, default=200)
    parser.add_argument("--grad_clip", type=float, default=1.0)
    parser.add_argument("--patience", type=int, default=3)
    parser.add_argument("--max_len", type=int, default=512)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--device", default="cuda:0")
    parser.add_argument("--num_workers", type=int, default=4)
    args = parser.parse_args()

    os.makedirs(args.out, exist_ok=True)
    torch.manual_seed(args.seed)
    np.random.seed(args.seed)

    tag = "frozen" if args.freeze_encoder else "finetuned"
    print(f"{'='*60}")
    print(f"CoordRep-Ranker Training ({tag})")
    print(f"{'='*60}")
    print(f"  loss={args.loss}, pooling={args.pooling}")
    print(f"  freeze_encoder={args.freeze_encoder}")
    print(f"  decoys_per_real={args.decoys_per_real}, batch_size={args.batch_size}")
    print(f"  lr={args.lr}, encoder_lr={args.encoder_lr}")

    # ── 1. Load data & split ──────────────────────────────────
    print("\nLoading data …")
    complexes = _extract_unique_complexes(args.data)
    train_cx, val_cx, test_cx = split_complexes(complexes, seed=args.seed)
    print(f"  train={len(train_cx)}, val={len(val_cx)}, test={len(test_cx)}")

    # Save split IDs for reproducibility
    split_path = os.path.join(args.out, "split_ids.json")
    with open(split_path, "w") as f:
        json.dump({
            "train": [c["id"] for c in train_cx],
            "val": [c["id"] for c in val_cx],
            "test": [c["id"] for c in test_cx],
        }, f)
    print(f"  Saved split IDs → {split_path}")

    # ── 2. Build hard decoy pool (train only!) ────────────────
    print("Building hard decoy pool (train split only) …")
    t0 = time.time()
    train_toks = [c["tokens"] for c in train_cx]
    pool = HardDecoyPool(train_toks, seed=args.seed)
    print(f"  {len(pool.metal_cns)} metals, "
          f"{sum(len(v) for v in pool.lig_by_len.values())} ligands "
          f"({time.time()-t0:.1f}s)")

    # ── 3. Build model ────────────────────────────────────────
    print("Loading ranker from MLM checkpoint …")
    ranker, tokenizer = load_ranker_from_mlm(
        args.mlm_checkpoint,
        device=args.device,
        pooling=args.pooling,
        freeze_encoder=args.freeze_encoder,
    )
    print(f"  Total params: {ranker.total_parameters():,}")
    print(f"  Trainable:    {ranker.trainable_parameters():,}")

    # ── 4. Datasets & loaders ─────────────────────────────────
    print("Building datasets …")
    train_ds = RankerListwiseDataset(
        train_cx, pool, tokenizer,
        n_decoys=args.decoys_per_real, max_len=args.max_len, seed=args.seed,
    )
    val_ds = RankerEvalDataset(
        val_cx, pool, tokenizer,
        n_decoys=args.decoys_per_real, max_len=args.max_len, seed=args.seed + 1,
    )
    print(f"  train: {len(train_ds)}, val: {len(val_ds)}")

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

    # ── 5. Optimiser & scheduler ──────────────────────────────
    loss_fn = get_loss_fn(args.loss)
    loss_fn = loss_fn.to(args.device)

    if args.freeze_encoder:
        optimizer = torch.optim.AdamW(
            ranker.head.parameters(), lr=args.lr, weight_decay=args.weight_decay,
        )
    else:
        # differential LR: head gets higher LR, encoder gets lower
        optimizer = torch.optim.AdamW([
            {"params": ranker.head.parameters(), "lr": args.lr},
            {"params": ranker.encoder.parameters(), "lr": args.encoder_lr},
        ], weight_decay=args.weight_decay)

    total_steps = len(train_loader) * args.epochs

    def lr_schedule(step):
        if step < args.warmup_steps:
            return step / max(args.warmup_steps, 1)
        progress = (step - args.warmup_steps) / max(total_steps - args.warmup_steps, 1)
        return 0.5 * (1 + math.cos(math.pi * progress))

    scheduler = torch.optim.lr_scheduler.LambdaLR(optimizer, lr_schedule)

    # ── 6. Training loop ──────────────────────────────────────
    best_val_mrr = 0.0
    patience_counter = 0
    history = []

    for epoch in range(1, args.epochs + 1):
        ranker.train()
        epoch_loss = 0.0
        epoch_top1 = 0
        epoch_n = 0
        t_epoch = time.time()

        for step, batch in enumerate(train_loader):
            ids = batch["input_ids"].to(args.device)      # (B, 1+K, L)
            att = batch["attention_mask"].to(args.device)

            scores = ranker.score_list(ids, att)  # (B, 1+K)
            loss = loss_fn(scores)

            optimizer.zero_grad()
            loss.backward()
            if args.grad_clip > 0:
                nn.utils.clip_grad_norm_(ranker.parameters(), args.grad_clip)
            optimizer.step()
            scheduler.step()

            B = ids.size(0)
            epoch_loss += loss.item() * B
            ranks = (scores[:, 1:] >= scores[:, :1]).sum(dim=1) + 1
            epoch_top1 += (ranks == 1).sum().item()
            epoch_n += B

            if step % 50 == 0:
                print(f"  Epoch {epoch} step {step}/{len(train_loader)} "
                      f"loss={loss.item():.4f} "
                      f"lr={scheduler.get_last_lr()[0]:.2e}", flush=True)

        train_loss = epoch_loss / max(epoch_n, 1)
        train_top1 = epoch_top1 / max(epoch_n, 1)

        # Validate
        val_metrics = evaluate(ranker, val_loader, loss_fn, args.device)
        elapsed = time.time() - t_epoch

        rec = {
            "epoch": epoch,
            "train_loss": round(train_loss, 4),
            "train_top1": round(train_top1, 4),
            "val_loss": round(val_metrics["loss"], 4),
            "val_top1": round(val_metrics["top1"], 4),
            "val_mrr": round(val_metrics["mrr"], 4),
            "time": round(elapsed, 1),
        }
        history.append(rec)

        print(f"\n  Epoch {epoch}: "
              f"train_loss={train_loss:.4f} train_top1={train_top1:.1%} | "
              f"val_loss={val_metrics['loss']:.4f} val_top1={val_metrics['top1']:.1%} "
              f"val_mrr={val_metrics['mrr']:.3f} ({elapsed:.0f}s)")

        # Early stopping on val MRR
        if val_metrics["mrr"] > best_val_mrr:
            best_val_mrr = val_metrics["mrr"]
            patience_counter = 0
            save_ranker(ranker, os.path.join(args.out, f"best_{tag}.pt"),
                        extra={"epoch": epoch, "val_mrr": best_val_mrr})
            print(f"  ★ New best val_mrr={best_val_mrr:.4f}, saved.")
        else:
            patience_counter += 1
            if patience_counter >= args.patience:
                print(f"  Early stopping at epoch {epoch}")
                break

    # Save final
    save_ranker(ranker, os.path.join(args.out, f"final_{tag}.pt"),
                extra={"epoch": epoch})

    # Save history
    hist_path = os.path.join(args.out, f"history_{tag}.json")
    with open(hist_path, "w") as f:
        json.dump(history, f, indent=2)
    print(f"\nTraining complete. Best val_mrr={best_val_mrr:.4f}")
    print(f"  Checkpoints: {args.out}/best_{tag}.pt, final_{tag}.pt")
    print(f"  History: {hist_path}")


if __name__ == "__main__":
    main()
