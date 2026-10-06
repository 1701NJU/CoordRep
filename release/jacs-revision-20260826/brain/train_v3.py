"""
CoordRep Mini Brain v3

1. Ligand-only
2. 4
3. Constraint exact match (T2)
4. Valid syntax/chem @top-k
5. Epoch 1
"""

import os
import sys
import json
import argparse
import hashlib
import re
import warnings
from pathlib import Path
from datetime import datetime
from collections import defaultdict
import numpy as np

warnings.filterwarnings('ignore')
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

import torch
import torch.nn as nn
from torch.optim import AdamW
from torch.optim.lr_scheduler import CosineAnnealingWarmRestarts
from torch.utils.data import DataLoader
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).parent.parent))

from brain.tokenizer import CoordRepTokenizer
from brain.masking import CoordinatedMasker, MaskingConfig
from brain.dataloader import CoordRepDataset, DataConfig
from brain.model import CoordRepModelConfig, CoordRepForMLM, count_parameters
from brain.smart_truncation import TruncationMonitor


class EnhancedMetricsLogger:
    """"""

    def __init__(self, output_dir: Path):
        self.output_dir = output_dir
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.metrics_file = output_dir / "metrics.jsonl"
        self.step = 0

    def log(self, metrics: dict, step: int = None):
        if step is not None:
            self.step = step

        entry = {
            "step": self.step,
            "timestamp": datetime.now().isoformat(),
            **metrics
        }

        with open(self.metrics_file, 'a') as f:
            f.write(json.dumps(entry) + '\n')

        self.step += 1

    def log_epoch(self, epoch: int, metrics: dict):
        """ epoch """
        epoch_file = self.output_dir / "epoch_metrics.jsonl"
        entry = {
            "epoch": epoch,
            "timestamp": datetime.now().isoformat(),
            **metrics
        }
        with open(epoch_file, 'a') as f:
            f.write(json.dumps(entry) + '\n')


def syntax_valid(s: str) -> bool:
    """"""
    try:
        if not s.strip():
            return False
        refs = set(re.findall(r'L(\d+)(?=[:@\)])', s.split('|')[0]))
        defs = set(re.findall(r'\|L(\d+)=', s))
        if refs and not refs.issubset(defs):
            return False
        return True
    except:
        return False


def chem_valid(s: str) -> bool:
    """"""
    if not syntax_valid(s):
        return False
    return True


def constraint_exact_match(pred_str: str, target_str: str) -> float:
    """ exact match"""
    pred_trans = set(re.findall(r'\{trans:([^}]+)\}', pred_str))
    target_trans = set(re.findall(r'\{trans:([^}]+)\}', target_str))

    if not target_trans:
        return 1.0 if not pred_trans else 0.0

    return len(pred_trans & target_trans) / len(target_trans)


def validate_generation_enhanced(
    model,
    tokenizer,
    device,
    n_samples: int = 50,
    topk: int = 5
):
    """"""
    model.eval()

    results = {
        "syntax_valid_top1": 0,
        "syntax_valid_top5": 0,
        "chem_valid_top1": 0,
        "constraint_exact_match": 0,
    }
    total = 0

    templates = [
        "[Fe;ox=+2;CN=6]<Shape:[MASK]>[MASK]|L1=NCCN|",
        "[Zn;ox=+2;CN=4]<Shape:[MASK]>|L1=Cl|",
        "[Pt;ox=+2;CN=4]<Shape:[MASK]>{trans:[MASK]}|L1=Cl||L2=N|",
    ]

    for template in templates * (n_samples // len(templates)):
        try:
            input_ids = tokenizer.encode(template)
            input_tensor = torch.tensor([input_ids]).to(device)

            with torch.no_grad():
                logits, _ = model(input_tensor)

            # Top-1
            top1_preds = logits.argmax(-1)[0]
            top1_output = tokenizer.decode(top1_preds.tolist())

            if syntax_valid(top1_output):
                results["syntax_valid_top1"] += 1
            if chem_valid(top1_output):
                results["chem_valid_top1"] += 1

            # Top-k
            topk_valid = False
            for k in range(min(topk, logits.size(-1))):
                topk_preds = logits.topk(k+1, dim=-1).indices[0, :, k]
                topk_output = tokenizer.decode(topk_preds.tolist())
                if syntax_valid(topk_output):
                    topk_valid = True
                    break

            if topk_valid:
                results["syntax_valid_top5"] += 1

            # Constraint match
            match = constraint_exact_match(top1_output, template)
            results["constraint_exact_match"] += match

            total += 1
        except:
            continue

    if total > 0:
        for k in results:
            results[k] /= total

    return results


def train_epoch(
    model,
    dataloader,
    optimizer,
    scheduler,
    device,
    tokenizer,
    logger: EnhancedMetricsLogger,
    trunc_monitor: TruncationMonitor,
    log_every: int = 50
):
    """ epoch"""
    model.train()
    total_loss = 0
    num_batches = 0

    field_stats = defaultdict(lambda: {"correct": 0, "total": 0})

    pbar = tqdm(dataloader, desc="Training")
    for batch_idx, batch in enumerate(pbar):
        input_ids = batch["input_ids"].to(device)
        attention_mask = batch["attention_mask"].to(device)
        labels = batch["labels"].to(device)

        optimizer.zero_grad()
        logits, loss = model(input_ids, attention_mask, labels)

        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        optimizer.step()
        scheduler.step()

        total_loss += loss.item()
        num_batches += 1

        preds = logits.argmax(-1)
        mask = labels != -100

        for i in range(labels.size(0)):
            for j in range(labels.size(1)):
                if labels[i, j] != -100:
                    tok_id = labels[i, j].item()
                    tok = tokenizer.id2token.get(tok_id, "")

                    if "Shape" in tok or tok.startswith("V_"):
                        field = "Shape"
                    elif "trans" in tok or "cis" in tok:
                        field = "Constraint"
                    elif tok.startswith("[") or tok.startswith(";"):
                        field = "Metal"
                    else:
                        field = "Ligand"

                    field_stats[field]["total"] += 1
                    if preds[i, j] == labels[i, j]:
                        field_stats[field]["correct"] += 1

        pbar.set_postfix({
            "loss": f"{loss.item():.4f}",
            "lr": f"{scheduler.get_last_lr()[0]:.2e}"
        })

        if (batch_idx + 1) % log_every == 0:
            field_accs = {}
            for field in ["Metal", "Shape", "Constraint", "Ligand"]:
                stats = field_stats[field]
                if stats["total"] > 0:
                    field_accs[f"acc_{field.lower()}"] = stats["correct"] / stats["total"]

            trunc_rates = trunc_monitor.get_rates()

            logger.log({
                "train_loss": loss.item(),
                "learning_rate": scheduler.get_last_lr()[0],
                "perplexity": np.exp(loss.item()),
                **field_accs,
                **trunc_rates
            })

    field_accs = {}
    for field in ["Metal", "Shape", "Constraint", "Ligand"]:
        stats = field_stats[field]
        if stats["total"] > 0:
            field_accs[f"acc_{field.lower()}"] = stats["correct"] / stats["total"]

    return total_loss / num_batches, field_accs


def quick_eval(model, tokenizer, device, epoch: int):
    """ epoch """
    print(f"\n--- Quick Eval (Epoch {epoch}) ---")

    metrics = validate_generation_enhanced(model, tokenizer, device, n_samples=50)

    print(f"  syntax_valid@top1: {metrics['syntax_valid_top1']*100:.1f}%")
    print(f"  syntax_valid@top5: {metrics['syntax_valid_top5']*100:.1f}%")
    print(f"  chem_valid@top1: {metrics['chem_valid_top1']*100:.1f}%")
    print(f"  constraint_exact_match: {metrics['constraint_exact_match']*100:.1f}%")

    return metrics


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tmqm-dir", default="../tmQM-master/tmQM")
    parser.add_argument("--cod-dir", default=None)
    parser.add_argument("--output-dir", default="./checkpoints/train_v3")
    parser.add_argument("--epochs", type=int, default=10)
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--lr", type=float, default=2e-4)
    parser.add_argument("--warmup-ratio", type=float, default=0.05)
    parser.add_argument("--max-length", type=int, default=768)
    parser.add_argument("--config", choices=["small", "tiny"], default="small")
    parser.add_argument("--log-every", type=int, default=50)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--curriculum", action="store_true")
    args = parser.parse_args()

    torch.manual_seed(args.seed)
    np.random.seed(args.seed)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Device: {device}")

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    config_dict = vars(args)
    config_dict["device"] = str(device)
    config_dict["timestamp"] = datetime.now().isoformat()
    config_dict["git_hash"] = "unknown"

    with open(output_dir / "config.json", 'w') as f:
        json.dump(config_dict, f, indent=2)

    print("\n=== Initializing ===")

    tokenizer = CoordRepTokenizer()
    masker = CoordinatedMasker()
    logger = EnhancedMetricsLogger(output_dir)
    trunc_monitor = TruncationMonitor()

    data_config = DataConfig(
        tmqm_dir=args.tmqm_dir,
        cod_dir=None, # Curriculum: tmQM
        max_length=args.max_length,
        batch_size=args.batch_size,
        shuffle_ligands=True,
        n_augments=2
    )

    model_config = (
        CoordRepModelConfig.small(max_length=args.max_length) if args.config == "small"
        else CoordRepModelConfig.tiny(max_length=args.max_length)
    )
    model_config.vocab_size = max(tokenizer.config.vocab_size, tokenizer.vocab_size + 100)

    print("\n=== Loading Data ===")
    train_dataset = CoordRepDataset(data_config, tokenizer, masker, split="train")
    train_loader = DataLoader(
        train_dataset, batch_size=args.batch_size, shuffle=True, num_workers=0
    )

    print("\n=== Creating Model ===")
    model = CoordRepForMLM(model_config).to(device)
    print(f"Parameters: {count_parameters(model):,}")

    optimizer = AdamW(model.parameters(), lr=args.lr, weight_decay=0.01)
    scheduler = CosineAnnealingWarmRestarts(optimizer, T_0=len(train_loader), T_mult=1)

    print("\n=== Training ===")
    best_loss = float("inf")

    for epoch in range(args.epochs):
        print(f"\n{'='*60}")
        print(f"Epoch {epoch + 1}/{args.epochs}")
        print(f"{'='*60}")

        if args.curriculum and epoch == 2 and args.cod_dir:
            print("  Adding COD data (curriculum)")
            data_config.cod_dir = args.cod_dir
            train_dataset = CoordRepDataset(data_config, tokenizer, masker, split="train")
            train_loader = DataLoader(
                train_dataset, batch_size=args.batch_size, shuffle=True, num_workers=0
            )

        train_loss, field_accs = train_epoch(
            model, train_loader, optimizer, scheduler, device,
            tokenizer, logger, trunc_monitor, args.log_every
        )

        print(f"\nTrain Loss: {train_loss:.4f}, PPL: {np.exp(train_loss):.2f}")
        print(f"Field Accuracies: {field_accs}")

        eval_metrics = quick_eval(model, tokenizer, device, epoch + 1)

        trunc_rates = trunc_monitor.get_rates()
        logger.log_epoch(epoch + 1, {
            "train_loss": train_loss,
            "perplexity": np.exp(train_loss),
            **field_accs,
            **eval_metrics,
            **trunc_rates
        })

        # Save
        checkpoint = {
            "epoch": epoch + 1,
            "model_state_dict": model.state_dict(),
            "optimizer_state_dict": optimizer.state_dict(),
            "loss": train_loss,
            "config": model_config.__dict__,
            "metrics": {**field_accs, **eval_metrics}
        }

        torch.save(checkpoint, output_dir / f"checkpoint_epoch_{epoch+1}.pt")

        if train_loss < best_loss:
            best_loss = train_loss
            torch.save(checkpoint, output_dir / "best_model.pt")
            print(f"  ✓ Saved best model")

        if epoch == 0:
            if eval_metrics['syntax_valid_top1'] < 0.5:
                print("\n  ⚠️ Warning: syntax_valid@top1 < 50% after epoch 1")
                print("  ⚠️ Consider checking masking/truncation settings")

    tokenizer.save(str(output_dir / "tokenizer.json"))

    print(f"\n{'='*60}")
    print("=== Training Complete ===")
    print(f"{'='*60}")
    print(f"Best loss: {best_loss:.4f}")
    print(f"Models saved to {output_dir}")


if __name__ == "__main__":
    main()
