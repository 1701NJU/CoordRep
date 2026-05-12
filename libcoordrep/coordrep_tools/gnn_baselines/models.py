"""
models.py
=========
GNN architectures for donor annotation and compatibility ranking.

Models:
1. GCN — Graph Convolutional Network (node classification)
2. GIN — Graph Isomorphism Network (node classification)
3. GAT — Graph Attention Network (node classification)
4. EGNN — E(n) Equivariant GNN (3D, node classification)
5. GINRanker — GIN-based graph-level ranker (for hard negatives)
"""

from __future__ import annotations

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import (
    GCNConv, GINConv, GATConv,
    global_mean_pool, global_add_pool,
)


# ── Donor Annotation Models (node classification) ────────

class GCNDonorPredictor(nn.Module):
    """GCN for atom-level donor annotation."""

    def __init__(self, in_dim: int, hidden_dim: int = 128, n_layers: int = 4,
                 dropout: float = 0.1, context_dim: int = 0):
        super().__init__()
        self.convs = nn.ModuleList()
        self.bns = nn.ModuleList()
        self.convs.append(GCNConv(in_dim + context_dim, hidden_dim))
        self.bns.append(nn.BatchNorm1d(hidden_dim))
        for _ in range(n_layers - 1):
            self.convs.append(GCNConv(hidden_dim, hidden_dim))
            self.bns.append(nn.BatchNorm1d(hidden_dim))
        self.head = nn.Linear(hidden_dim, 1)
        self.dropout = dropout

    def forward(self, x, edge_index, batch=None, context=None):
        if context is not None:
            # Expand context (graph-level) to each node
            if batch is not None:
                ctx = context[batch]
            else:
                ctx = context.expand(x.size(0), -1)
            x = torch.cat([x, ctx], dim=-1)

        for conv, bn in zip(self.convs, self.bns):
            x = conv(x, edge_index)
            x = bn(x)
            x = F.relu(x)
            x = F.dropout(x, p=self.dropout, training=self.training)
        return self.head(x).squeeze(-1)


class GINDonorPredictor(nn.Module):
    """GIN for atom-level donor annotation."""

    def __init__(self, in_dim: int, hidden_dim: int = 128, n_layers: int = 4,
                 dropout: float = 0.1, context_dim: int = 0):
        super().__init__()
        self.convs = nn.ModuleList()
        self.bns = nn.ModuleList()
        # First layer
        mlp0 = nn.Sequential(
            nn.Linear(in_dim + context_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim),
        )
        self.convs.append(GINConv(mlp0))
        self.bns.append(nn.BatchNorm1d(hidden_dim))
        for _ in range(n_layers - 1):
            mlp = nn.Sequential(
                nn.Linear(hidden_dim, hidden_dim),
                nn.ReLU(),
                nn.Linear(hidden_dim, hidden_dim),
            )
            self.convs.append(GINConv(mlp))
            self.bns.append(nn.BatchNorm1d(hidden_dim))
        self.head = nn.Linear(hidden_dim, 1)
        self.dropout = dropout

    def forward(self, x, edge_index, batch=None, context=None):
        if context is not None:
            if batch is not None:
                ctx = context[batch]
            else:
                ctx = context.expand(x.size(0), -1)
            x = torch.cat([x, ctx], dim=-1)

        for conv, bn in zip(self.convs, self.bns):
            x = conv(x, edge_index)
            x = bn(x)
            x = F.relu(x)
            x = F.dropout(x, p=self.dropout, training=self.training)
        return self.head(x).squeeze(-1)


class GATDonorPredictor(nn.Module):
    """GAT for atom-level donor annotation."""

    def __init__(self, in_dim: int, hidden_dim: int = 128, n_layers: int = 4,
                 heads: int = 4, dropout: float = 0.1, context_dim: int = 0):
        super().__init__()
        self.convs = nn.ModuleList()
        self.bns = nn.ModuleList()
        self.convs.append(GATConv(in_dim + context_dim, hidden_dim // heads, heads=heads))
        self.bns.append(nn.BatchNorm1d(hidden_dim))
        for _ in range(n_layers - 1):
            self.convs.append(GATConv(hidden_dim, hidden_dim // heads, heads=heads))
            self.bns.append(nn.BatchNorm1d(hidden_dim))
        self.head = nn.Linear(hidden_dim, 1)
        self.dropout = dropout

    def forward(self, x, edge_index, batch=None, context=None):
        if context is not None:
            if batch is not None:
                ctx = context[batch]
            else:
                ctx = context.expand(x.size(0), -1)
            x = torch.cat([x, ctx], dim=-1)

        for conv, bn in zip(self.convs, self.bns):
            x = conv(x, edge_index)
            x = bn(x)
            x = F.relu(x)
            x = F.dropout(x, p=self.dropout, training=self.training)
        return self.head(x).squeeze(-1)


# ── EGNN (simplified E(n)-equivariant) ───────────────────

class EGNNLayer(nn.Module):
    """Simplified EGNN layer: message passing with coordinate updates."""

    def __init__(self, hidden_dim: int):
        super().__init__()
        self.edge_mlp = nn.Sequential(
            nn.Linear(2 * hidden_dim + 1, hidden_dim),
            nn.SiLU(),
            nn.Linear(hidden_dim, hidden_dim),
        )
        self.node_mlp = nn.Sequential(
            nn.Linear(2 * hidden_dim, hidden_dim),
            nn.SiLU(),
            nn.Linear(hidden_dim, hidden_dim),
        )
        self.coord_mlp = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim),
            nn.SiLU(),
            nn.Linear(hidden_dim, 1, bias=False),
        )

    def forward(self, h, pos, edge_index):
        row, col = edge_index
        diff = pos[row] - pos[col]
        dist = (diff ** 2).sum(dim=-1, keepdim=True).clamp(min=1e-6)

        # Edge messages
        edge_feat = torch.cat([h[row], h[col], dist], dim=-1)
        m_ij = self.edge_mlp(edge_feat)

        # Aggregate
        agg = torch.zeros_like(h)
        agg.index_add_(0, row, m_ij)

        # Update nodes
        h_new = self.node_mlp(torch.cat([h, agg], dim=-1))
        h = h + h_new

        # Update coordinates
        coord_weights = self.coord_mlp(m_ij)
        coord_update = torch.zeros_like(pos)
        coord_update.index_add_(0, row, diff * coord_weights)
        pos = pos + coord_update

        return h, pos


class EGNNDonorPredictor(nn.Module):
    """EGNN for 3D-aware donor annotation (upper bound)."""

    def __init__(self, in_dim: int, hidden_dim: int = 128, n_layers: int = 4,
                 dropout: float = 0.1):
        super().__init__()
        self.embed = nn.Linear(in_dim, hidden_dim)
        self.layers = nn.ModuleList([EGNNLayer(hidden_dim) for _ in range(n_layers)])
        self.head = nn.Linear(hidden_dim, 1)
        self.dropout = dropout

    def forward(self, x, edge_index, pos, batch=None, context=None):
        h = self.embed(x)
        for layer in self.layers:
            h, pos = layer(h, pos, edge_index)
            h = F.dropout(h, p=self.dropout, training=self.training)
        return self.head(h).squeeze(-1)


# ── Graph-level Ranker (for hard-negative compatibility) ──

class GINRanker(nn.Module):
    """GIN-based graph-level ranker producing a scalar score."""

    def __init__(self, in_dim: int, hidden_dim: int = 128, n_layers: int = 4,
                 dropout: float = 0.1):
        super().__init__()
        self.convs = nn.ModuleList()
        self.bns = nn.ModuleList()
        mlp0 = nn.Sequential(
            nn.Linear(in_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim),
        )
        self.convs.append(GINConv(mlp0))
        self.bns.append(nn.BatchNorm1d(hidden_dim))
        for _ in range(n_layers - 1):
            mlp = nn.Sequential(
                nn.Linear(hidden_dim, hidden_dim),
                nn.ReLU(),
                nn.Linear(hidden_dim, hidden_dim),
            )
            self.convs.append(GINConv(mlp))
            self.bns.append(nn.BatchNorm1d(hidden_dim))

        self.pool_head = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, 1),
        )
        self.dropout = dropout

    def forward(self, x, edge_index, batch):
        for conv, bn in zip(self.convs, self.bns):
            x = conv(x, edge_index)
            x = bn(x)
            x = F.relu(x)
            x = F.dropout(x, p=self.dropout, training=self.training)
        pooled = global_add_pool(x, batch)
        return self.pool_head(pooled).squeeze(-1)

    def embed(self, x, edge_index, batch):
        """Return graph-level embedding without scoring head."""
        for conv, bn in zip(self.convs, self.bns):
            x = conv(x, edge_index)
            x = bn(x)
            x = F.relu(x)
        return global_add_pool(x, batch)
