#!/usr/bin/env python
"""Official PyG ViSNet trunk with CoordStatePairs two-property readouts."""

from __future__ import annotations

import numpy as np
import torch
import torch.nn as nn
from rdkit import Chem
from torch_geometric.nn import global_add_pool, global_mean_pool
from torch_geometric.nn.models import ViSNet


class MultiTargetViSNet(nn.Module):
    """Official PyG ViSNet vector-scalar trunk with invariant molecular heads.

    The constructor mirrors ``MultiTargetSchNet`` so the exact frozen OOF
    runner, cohort, folds, validation policy, and metrics can be reused.
    ``interactions`` maps to the number of ViSNet layers; the frozen production
    configuration uses four layers and eight attention heads.
    """

    def __init__(
        self,
        hidden: int,
        interactions: int,
        num_gaussians: int,
        cutoff: float,
        max_neighbors: int,
        n_targets: int,
        target_mean: np.ndarray,
        target_scale: np.ndarray,
    ) -> None:
        super().__init__()
        if n_targets != 2:
            raise ValueError("The CoordStatePairs benchmark expects exactly two targets")
        model = ViSNet(
            lmax=1,
            vecnorm_type=None,
            trainable_vecnorm=False,
            num_heads=8,
            num_layers=int(interactions),
            hidden_channels=int(hidden),
            num_rbf=int(num_gaussians),
            trainable_rbf=False,
            max_z=100,
            cutoff=float(cutoff),
            max_num_neighbors=int(max_neighbors),
            vertex=False,
            reduce_op="sum",
            derivative=False,
        )
        # Register only the official vector-scalar representation trunk. The
        # stock scalar energy head is intentionally replaced by the same two
        # target-appropriate heads used for the common-cohort 3D comparison.
        self.representation_model = model.representation_model
        head_hidden = max(int(hidden) // 2, 16)
        self.readout_norm = nn.LayerNorm(int(hidden))
        self.gap_head = nn.Sequential(
            nn.Linear(int(hidden), head_hidden),
            nn.SiLU(),
            nn.Linear(head_hidden, 1),
        )
        self.charge_head = nn.Sequential(
            nn.Linear(int(hidden), head_hidden),
            nn.SiLU(),
            nn.Linear(head_hidden, 1),
        )

        periodic = Chem.GetPeriodicTable()
        atomic_mass = torch.zeros(120, dtype=torch.float32)
        for atomic_number in range(1, 119):
            atomic_mass[atomic_number] = float(periodic.GetAtomicWeight(atomic_number))
        self.register_buffer("atomic_mass", atomic_mass)
        self.register_buffer(
            "target_mean", torch.as_tensor(target_mean, dtype=torch.float32).view(1, 2)
        )
        self.register_buffer(
            "target_scale", torch.as_tensor(target_scale, dtype=torch.float32).view(1, 2)
        )

    def forward(self, batch: object) -> torch.Tensor:
        scalar, _vector = self.representation_model(batch.z, batch.pos, batch.batch)
        scalar = self.readout_norm(scalar)
        gap_scaled = global_mean_pool(self.gap_head(scalar), batch.batch)

        mass = self.atomic_mass[batch.z].view(-1, 1)
        total_mass = global_add_pool(mass, batch.batch).clamp_min(1.0e-12)
        center_of_mass = global_add_pool(mass * batch.pos, batch.batch) / total_mass
        # These unconstrained latent scalar weights are used only to construct
        # an invariant vector-norm readout; they are not physical charges.
        latent_weight = self.charge_head(scalar)
        dipole_vector = global_add_pool(
            latent_weight * (batch.pos - center_of_mass.index_select(0, batch.batch)),
            batch.batch,
        )
        dipole_raw = torch.linalg.vector_norm(dipole_vector, dim=-1, keepdim=True)
        dipole_scaled = (dipole_raw - self.target_mean[:, 1:2]) / self.target_scale[:, 1:2]
        return torch.cat([gap_scaled, dipole_scaled], dim=-1)
