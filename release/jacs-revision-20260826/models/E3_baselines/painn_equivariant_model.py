#!/usr/bin/env python
"""PaiNN-style E(3)-equivariant baseline for the CoordStatePairs audit.

The implementation follows the scalar/vector message-passing construction of
Schuett et al., ICML 2021. Scalar channels are rotation invariant and vector
channels transform equivariantly. The model consumes only atomic numbers and
Cartesian coordinates, matching the information supplied to the SchNet
comparator. The two readouts are an invariant scalar gap head and an
equivariant latent-charge dipole vector whose norm is invariant.

This module intentionally exposes the same constructor signature as
``MultiTargetSchNet`` so that the frozen CoordStatePairs OOF runner can be
reused without changing folds, targets, validation policy, or output schema.
"""

from __future__ import annotations

import math

import numpy as np
import torch
import torch.nn as nn
from rdkit import Chem
from torch_geometric.nn import global_add_pool, global_mean_pool, radius_graph


class GaussianRBF(nn.Module):
    def __init__(self, n_rbf: int, cutoff: float) -> None:
        super().__init__()
        if n_rbf < 2:
            raise ValueError("n_rbf must be at least two")
        centers = torch.linspace(0.0, float(cutoff), int(n_rbf))
        spacing = float(centers[1] - centers[0])
        self.register_buffer("centers", centers)
        self.gamma = 0.5 / max(spacing * spacing, 1.0e-12)
        self.cutoff = float(cutoff)

    def forward(self, distance: torch.Tensor) -> torch.Tensor:
        rbf = torch.exp(-self.gamma * (distance.unsqueeze(-1) - self.centers) ** 2)
        envelope = 0.5 * (
            torch.cos(math.pi * distance.clamp(max=self.cutoff) / self.cutoff) + 1.0
        )
        envelope = envelope * (distance < self.cutoff).to(distance.dtype)
        return rbf * envelope.unsqueeze(-1)


class PaiNNInteraction(nn.Module):
    """One equivariant message and invariant/equivariant update block."""

    def __init__(self, hidden: int, n_rbf: int) -> None:
        super().__init__()
        self.hidden = int(hidden)
        self.scalar_message = nn.Sequential(
            nn.Linear(hidden, hidden),
            nn.SiLU(),
            nn.Linear(hidden, 3 * hidden),
        )
        self.radial_filter = nn.Sequential(
            nn.Linear(n_rbf, hidden),
            nn.SiLU(),
            nn.Linear(hidden, 3 * hidden),
        )
        self.vector_u = nn.Linear(hidden, hidden, bias=False)
        self.vector_v = nn.Linear(hidden, hidden, bias=False)
        self.update_mlp = nn.Sequential(
            nn.Linear(2 * hidden, hidden),
            nn.SiLU(),
            nn.Linear(hidden, 3 * hidden),
        )
        self.scalar_norm = nn.LayerNorm(hidden)

    def forward(
        self,
        scalar: torch.Tensor,
        vector: torch.Tensor,
        edge_index: torch.Tensor,
        direction: torch.Tensor,
        rbf: torch.Tensor,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        source, target = edge_index[0], edge_index[1]
        message = self.scalar_message(scalar.index_select(0, source))
        message = message * self.radial_filter(rbf)
        message_scalar, message_vector_scale, message_direction_scale = message.chunk(3, -1)

        vector_source = vector.index_select(0, source)
        vector_message = (
            vector_source * message_vector_scale.unsqueeze(1)
            + direction.unsqueeze(-1) * message_direction_scale.unsqueeze(1)
        )
        # ``radius_graph`` aggregation targets are kept in the module state
        # dtype under AMP, so cast autocast messages before ``index_add_``.
        message_scalar = message_scalar.to(scalar.dtype)
        vector_message = vector_message.to(vector.dtype)
        scalar_aggregate = torch.zeros_like(scalar)
        vector_aggregate = torch.zeros_like(vector)
        scalar_aggregate.index_add_(0, target, message_scalar)
        vector_aggregate.index_add_(0, target, vector_message)
        scalar = scalar + scalar_aggregate
        vector = vector + vector_aggregate

        vector_u = self.vector_u(vector)
        vector_v = self.vector_v(vector)
        vector_norm = torch.sqrt(torch.sum(vector_v * vector_v, dim=1) + 1.0e-8)
        update = self.update_mlp(torch.cat([self.scalar_norm(scalar), vector_norm], dim=-1))
        scalar_shift, scalar_vector_shift, vector_scale = update.chunk(3, -1)
        vector_dot = torch.sum(vector_u * vector_v, dim=1)
        scalar = scalar + scalar_shift + scalar_vector_shift * vector_dot
        vector = vector + vector_scale.unsqueeze(1) * vector_u
        return scalar, vector


class MultiTargetPaiNN(nn.Module):
    """PaiNN scalar/vector trunk with invariant gap and dipole-magnitude heads."""

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
        self.hidden = int(hidden)
        self.cutoff = float(cutoff)
        self.max_neighbors = int(max_neighbors)
        self.embedding = nn.Embedding(120, hidden, padding_idx=0)
        self.rbf = GaussianRBF(num_gaussians, cutoff)
        self.interactions = nn.ModuleList(
            [PaiNNInteraction(hidden, num_gaussians) for _ in range(int(interactions))]
        )
        self.readout_norm = nn.LayerNorm(hidden)
        head_hidden = max(hidden // 2, 16)
        self.gap_head = nn.Sequential(
            nn.Linear(hidden, head_hidden), nn.SiLU(), nn.Linear(head_hidden, 1)
        )
        self.charge_head = nn.Sequential(
            nn.Linear(hidden, head_hidden), nn.SiLU(), nn.Linear(head_hidden, 1)
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
        z = batch.z
        pos = batch.pos
        batch_index = batch.batch
        edge_index = radius_graph(
            pos,
            r=self.cutoff,
            batch=batch_index,
            loop=False,
            max_num_neighbors=self.max_neighbors,
        )
        source, target = edge_index[0], edge_index[1]
        displacement = pos.index_select(0, source) - pos.index_select(0, target)
        distance = torch.linalg.vector_norm(displacement, dim=-1).clamp_min(1.0e-8)
        direction = displacement / distance.unsqueeze(-1)
        rbf = self.rbf(distance)

        scalar = self.embedding(z)
        vector = scalar.new_zeros((scalar.shape[0], 3, self.hidden))
        for interaction in self.interactions:
            scalar, vector = interaction(scalar, vector, edge_index, direction, rbf)
        scalar = self.readout_norm(scalar)

        gap_scaled = global_mean_pool(self.gap_head(scalar), batch_index)
        mass = self.atomic_mass[z].view(-1, 1)
        total_mass = global_add_pool(mass, batch_index).clamp_min(1.0e-12)
        center_of_mass = global_add_pool(mass * pos, batch_index) / total_mass
        latent_charge = self.charge_head(scalar)
        dipole_vector = global_add_pool(
            latent_charge * (pos - center_of_mass.index_select(0, batch_index)),
            batch_index,
        )
        dipole_raw = torch.linalg.vector_norm(dipole_vector, dim=-1, keepdim=True)
        dipole_scaled = (dipole_raw - self.target_mean[:, 1:2]) / self.target_scale[:, 1:2]
        return torch.cat([gap_scaled, dipole_scaled], dim=-1)
