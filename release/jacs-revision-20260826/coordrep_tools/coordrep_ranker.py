"""
coordrep_ranker.py
==================
CoordRep-Ranker: contrastive compatibility scorer.

Initialises the encoder from a pretrained CoordRepForMLM checkpoint,
adds a lightweight MLP scoring head that maps pooled embeddings → scalar
compatibility score.
"""

from __future__ import annotations

import copy
import json
import math
from pathlib import Path
from typing import Optional, Tuple

import torch
import torch.nn as nn
import torch.nn.functional as F

from brain.model import CoordRepModelConfig, CoordRepEncoder, CoordRepForMLM
from brain.tokenizer import CoordRepTokenizer, TokenizerConfig


# ── scoring head ──────────────────────────────────────────────────

class ScoringHead(nn.Module):
    """MLP: hidden_size → scalar score."""

    def __init__(self, hidden_size: int, dropout: float = 0.1):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(hidden_size, hidden_size),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_size, hidden_size // 2),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_size // 2, 1),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """x: (batch, hidden_size) → (batch,)"""
        return self.net(x).squeeze(-1)


# ── ranker model ──────────────────────────────────────────────────

class CoordRepRanker(nn.Module):
    """
    Encoder + pooling + scoring head.

    Supports two pooling modes:
      - "cls"  : use the [CLS] token embedding
      - "mean" : attention-masked mean pooling
    """

    def __init__(
        self,
        config: CoordRepModelConfig,
        pooling: str = "mean",
        head_dropout: float = 0.1,
    ):
        super().__init__()
        self.config = config
        self.pooling = pooling
        self.encoder = CoordRepEncoder(config)
        self.head = ScoringHead(config.hidden_size, dropout=head_dropout)

    # ── forward ───────────────────────────────────────────────

    def encode(
        self,
        input_ids: torch.Tensor,
        attention_mask: Optional[torch.Tensor] = None,
    ) -> torch.Tensor:
        """Return pooled embedding (batch, hidden_size)."""
        hidden = self.encoder(input_ids, attention_mask)  # (B, L, H)
        if self.pooling == "cls":
            return hidden[:, 0]
        # mean pooling
        if attention_mask is not None:
            mask = attention_mask.unsqueeze(-1).float()  # (B, L, 1)
            return (hidden * mask).sum(1) / mask.sum(1).clamp(min=1e-9)
        return hidden.mean(1)

    def forward(
        self,
        input_ids: torch.Tensor,
        attention_mask: Optional[torch.Tensor] = None,
    ) -> torch.Tensor:
        """Return scalar compatibility score (batch,)."""
        emb = self.encode(input_ids, attention_mask)
        return self.head(emb)

    def score_pair(
        self,
        real_ids: torch.Tensor,
        real_mask: torch.Tensor,
        decoy_ids: torch.Tensor,
        decoy_mask: torch.Tensor,
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        """Score a (real, decoy) pair. Returns (real_score, decoy_score)."""
        ids = torch.cat([real_ids, decoy_ids], dim=0)
        masks = torch.cat([real_mask, decoy_mask], dim=0)
        scores = self.forward(ids, masks)
        B = real_ids.size(0)
        return scores[:B], scores[B:]

    def score_list(
        self,
        input_ids: torch.Tensor,
        attention_mask: torch.Tensor,
    ) -> torch.Tensor:
        """
        Score a listwise batch.
        input_ids:      (batch, 1+K, seq_len)  — first is real, rest are decoys
        attention_mask:  (batch, 1+K, seq_len)

        Returns: (batch, 1+K) scores.
        """
        B, N, L = input_ids.shape
        flat_ids = input_ids.view(B * N, L)
        flat_mask = attention_mask.view(B * N, L)
        scores = self.forward(flat_ids, flat_mask)
        return scores.view(B, N)

    # ── freeze / unfreeze ─────────────────────────────────────

    def freeze_encoder(self):
        for p in self.encoder.parameters():
            p.requires_grad = False

    def unfreeze_encoder(self):
        for p in self.encoder.parameters():
            p.requires_grad = True

    def trainable_parameters(self):
        return sum(p.numel() for p in self.parameters() if p.requires_grad)

    def total_parameters(self):
        return sum(p.numel() for p in self.parameters())


# ── initialisation helpers ────────────────────────────────────────

def load_ranker_from_mlm(
    checkpoint_path: str,
    device: str = "cuda",
    pooling: str = "mean",
    head_dropout: float = 0.1,
    freeze_encoder: bool = False,
) -> Tuple[CoordRepRanker, CoordRepTokenizer]:
    """
    Create a CoordRepRanker, initialise the encoder from a pretrained
    CoordRepForMLM checkpoint, and return (ranker, tokenizer).
    """
    ckpt_dir = Path(checkpoint_path).parent
    tok_path = ckpt_dir / "tokenizer.json"

    # tokenizer
    with open(tok_path) as f:
        data = json.load(f)
    config_dict = data.get("config", {})
    tok_config = TokenizerConfig(**config_dict)
    tokenizer = CoordRepTokenizer(tok_config)
    tokenizer.token2id = data["token2id"]
    tokenizer.id2token = {v: k for k, v in tokenizer.token2id.items()}

    # model config (must match checkpoint)
    model_config = CoordRepModelConfig.small(max_length=768)
    model_config.vocab_size = 10000

    # load MLM checkpoint
    mlm = CoordRepForMLM(model_config)
    ckpt = torch.load(checkpoint_path, map_location="cpu", weights_only=False)
    mlm.load_state_dict(ckpt["model_state_dict"])

    # build ranker, copy encoder weights
    ranker = CoordRepRanker(model_config, pooling=pooling, head_dropout=head_dropout)
    ranker.encoder.load_state_dict(mlm.encoder.state_dict())

    if freeze_encoder:
        ranker.freeze_encoder()

    ranker = ranker.to(device)
    return ranker, tokenizer


def save_ranker(ranker: CoordRepRanker, path: str, extra: dict = None):
    """Save ranker checkpoint."""
    state = {
        "model_state_dict": ranker.state_dict(),
        "config": {
            "hidden_size": ranker.config.hidden_size,
            "num_hidden_layers": ranker.config.num_hidden_layers,
            "num_attention_heads": ranker.config.num_attention_heads,
            "intermediate_size": ranker.config.intermediate_size,
            "vocab_size": ranker.config.vocab_size,
            "max_position_embeddings": ranker.config.max_position_embeddings,
        },
        "pooling": ranker.pooling,
    }
    if extra:
        state.update(extra)
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    torch.save(state, path)


def load_ranker(path: str, tokenizer_path: str, device: str = "cuda") -> Tuple[CoordRepRanker, CoordRepTokenizer]:
    """Load a saved ranker checkpoint."""
    ckpt = torch.load(path, map_location="cpu", weights_only=False)
    cfg = ckpt["config"]
    model_config = CoordRepModelConfig(
        vocab_size=cfg["vocab_size"],
        hidden_size=cfg["hidden_size"],
        num_hidden_layers=cfg["num_hidden_layers"],
        num_attention_heads=cfg["num_attention_heads"],
        intermediate_size=cfg["intermediate_size"],
        max_position_embeddings=cfg["max_position_embeddings"],
    )
    pooling = ckpt.get("pooling", "mean")
    ranker = CoordRepRanker(model_config, pooling=pooling)
    ranker.load_state_dict(ckpt["model_state_dict"])
    ranker.eval().to(device)

    with open(tokenizer_path) as f:
        data = json.load(f)
    tok_config = TokenizerConfig(**data.get("config", {}))
    tokenizer = CoordRepTokenizer(tok_config)
    tokenizer.token2id = data["token2id"]
    tokenizer.id2token = {v: k for k, v in tokenizer.token2id.items()}

    return ranker, tokenizer
