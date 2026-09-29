"""
ranker_losses.py
================
Three loss functions for CoordRep-Ranker training.

A. PairwiseMarginLoss  — max(0, margin - s_real + s_decoy)
B. PairwiseBCELoss     — BCE with real→1, decoy→0
C. ListwiseCELoss      — cross-entropy over [real, decoy_1, …, decoy_K]
"""

from __future__ import annotations

import torch
import torch.nn as nn
import torch.nn.functional as F


class PairwiseMarginLoss(nn.Module):
    """max(0, margin - score_real + score_decoy)"""

    def __init__(self, margin: float = 1.0):
        super().__init__()
        self.margin = margin

    def forward(self, score_real: torch.Tensor, score_decoy: torch.Tensor) -> torch.Tensor:
        """
        score_real:  (batch,)
        score_decoy: (batch,)
        """
        return F.relu(self.margin - score_real + score_decoy).mean()


class PairwiseBCELoss(nn.Module):
    """BCE loss: real → 1, decoy → 0."""

    def forward(self, score_real: torch.Tensor, score_decoy: torch.Tensor) -> torch.Tensor:
        scores = torch.cat([score_real, score_decoy], dim=0)
        labels = torch.cat([
            torch.ones_like(score_real),
            torch.zeros_like(score_decoy),
        ], dim=0)
        return F.binary_cross_entropy_with_logits(scores, labels)


class ListwiseCELoss(nn.Module):
    """
    Cross-entropy over a list of [real, decoy_1, …, decoy_K].

    Input: scores of shape (batch, 1+K).
    Target: index 0 (real) is the correct class.
    """

    def __init__(self, temperature: float = 1.0):
        super().__init__()
        self.temperature = temperature

    def forward(self, scores: torch.Tensor) -> torch.Tensor:
        """scores: (batch, 1+K) where column 0 is real."""
        targets = torch.zeros(scores.size(0), dtype=torch.long, device=scores.device)
        return F.cross_entropy(scores / self.temperature, targets)


def get_loss_fn(name: str, **kwargs) -> nn.Module:
    """Factory helper."""
    if name == "margin":
        return PairwiseMarginLoss(margin=kwargs.get("margin", 1.0))
    if name == "bce":
        return PairwiseBCELoss()
    if name == "listwise":
        return ListwiseCELoss(temperature=kwargs.get("temperature", 1.0))
    raise ValueError(f"Unknown loss: {name}")
