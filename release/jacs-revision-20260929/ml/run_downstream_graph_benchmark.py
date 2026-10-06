#!/usr/bin/env python
"""Matched CoordRep-field, molecular-graph, and hybrid property benchmark.

This script is intentionally independent of the legacy CoordRep tokenizer.  It
uses the typed fields emitted by the current pipeline and compares them with a
standard attributed molecular-graph GIN on exactly the same tmQM records,
split, targets, training budget, and evaluation metrics.

The benchmark is designed for the JACS revision and writes auditable JSON/CSV
outputs rather than editing manuscript files.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import random
import re
import time
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from rdkit import Chem, RDLogger
from sklearn.feature_extraction import DictVectorizer
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import GroupShuffleSplit
from sklearn.preprocessing import StandardScaler
from torch.utils.data import DataLoader as TensorDataLoader
from torch.utils.data import TensorDataset
from torch_geometric.data import Data
from torch_geometric.loader import DataLoader as GraphDataLoader
from torch_geometric.nn import GINEConv, global_mean_pool


RDLogger.DisableLog("rdApp.*")


ATOM_DIMS = (119, 16, 13, 13, 10, 16, 2, 2)
BOND_DIMS = (7, 8, 2, 2)
FORMULA_RE = re.compile(r"([A-Z][a-z]?)(\d*)")
TARGET_SPECS = {
    "hl_gap_ev": ("hl_gap", 27.211386245988),
    "dipole_moment": ("dipole_moment", 1.0),
    "polarizability": ("polarizability", 1.0),
}


@dataclass
class PreparedData:
    ids: List[str]
    smiles: List[str]
    graphs: List[Data]
    fields: np.ndarray
    targets: np.ndarray
    target_names: List[str]
    field_names: List[str]
    train_idx: np.ndarray
    val_idx: np.ndarray
    test_idx: np.ndarray
    target_mean: np.ndarray
    target_scale: np.ndarray
    field_scaler_mean: np.ndarray
    field_scaler_scale: np.ndarray
    audit: Dict[str, object]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path("tmp/model_audit/data/CoordRep/CoordSMILES"),
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("revision_experiments/results/downstream_graph"),
    )
    parser.add_argument("--targets", nargs="+", default=["hl_gap_ev", "dipole_moment"])
    parser.add_argument("--split-seed", type=int, default=20260723)
    parser.add_argument("--seeds", nargs="+", type=int, default=[11, 22, 33])
    parser.add_argument("--epochs", type=int, default=40)
    parser.add_argument("--patience", type=int, default=8)
    parser.add_argument("--batch-size", type=int, default=192)
    parser.add_argument("--hidden", type=int, default=128)
    parser.add_argument("--layers", type=int, default=4)
    parser.add_argument("--lr", type=float, default=1.0e-3)
    parser.add_argument("--weight-decay", type=float, default=1.0e-5)
    parser.add_argument("--max-records", type=int, default=0)
    parser.add_argument("--models", nargs="+", default=["coordrep_mlp", "gin", "hybrid"])
    parser.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    parser.add_argument("--num-workers", type=int, default=0)
    parser.add_argument("--smoke", action="store_true")
    return parser.parse_args()


def set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_current_records(path: Path) -> Tuple[Dict[str, dict], Dict[str, int]]:
    records: Dict[str, dict] = {}
    counters = Counter()
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            counters["pipeline_rows"] += 1
            row = json.loads(line)
            mol_id = row.get("mol_id")
            if not mol_id or not row.get("coordrep"):
                counters["pipeline_missing_id_or_record"] += 1
                continue
            if mol_id in records:
                counters["pipeline_duplicate_ids"] += 1
                continue
            records[mol_id] = row
    counters["pipeline_unique_usable"] = len(records)
    return records, dict(counters)


def formula_counts(formula: str) -> Counter:
    counts: Counter = Counter()
    if not formula:
        return counts
    for element, raw_count in FORMULA_RE.findall(str(formula)):
        counts[element] += int(raw_count) if raw_count else 1
    return counts


def record_to_field_dict(record: dict) -> Dict[str, float]:
    """Extract only fields encoded in the current CoordRep structural record."""
    features: Dict[str, float] = {}
    metal = str(record.get("metal") or "UNK")
    features[f"metal={metal}"] = 1.0

    for key in ("oxidation", "dcount", "cn", "shape_cn", "n_ligands"):
        value = record.get(key)
        features[f"{key}_missing"] = float(value is None)
        if value is not None:
            features[key] = float(value)

    refs = record.get("shape_ref") or []
    values = record.get("shape_values") or []
    if refs and values:
        pairs = [(str(ref), float(value)) for ref, value in zip(refs, values)]
        for ref, value in pairs:
            features[f"cshm_{ref}"] = value
            features[f"shape_candidate={ref}"] = 1.0
        best_ref, best_value = min(pairs, key=lambda item: item[1])
        features[f"shape_best={best_ref}"] = 1.0
        features["cshm_best"] = best_value
        if len(pairs) > 1:
            ordered = sorted(value for _, value in pairs)
            features["cshm_margin"] = ordered[1] - ordered[0]

    dent_counts = Counter(int(x) for x in (record.get("ligand_dents") or []))
    for denticity, count in dent_counts.items():
        features[f"ligand_denticity_{denticity}"] = float(count)
    features["ligand_denticity_sum"] = float(sum((record.get("ligand_dents") or [])))
    features["ligand_denticity_max"] = float(max(record.get("ligand_dents") or [0]))

    donor_counts = Counter(str(x) for x in (record.get("donor_elements") or []))
    for element, count in donor_counts.items():
        features[f"donor_{element}"] = float(count)

    ligand_formula_total: Counter = Counter()
    for formula in record.get("ligand_smiles") or []:
        ligand_formula_total.update(formula_counts(str(formula)))
    for element, count in ligand_formula_total.items():
        features[f"ligand_formula_{element}"] = float(count)

    # Relation counts are explicitly serialized in braces in the current record.
    text = str(record.get("coordrep") or "")
    features["n_trans_relations"] = float(text.count("{trans:"))
    features["n_cis_relations"] = float(text.count("{cis:"))
    return features


def bond_type_index(bond: Chem.Bond) -> int:
    bond_type = bond.GetBondType()
    mapping = {
        Chem.BondType.SINGLE: 0,
        Chem.BondType.DOUBLE: 1,
        Chem.BondType.TRIPLE: 2,
        Chem.BondType.AROMATIC: 3,
        Chem.BondType.DATIVE: 4,
        Chem.BondType.ZERO: 5,
    }
    return mapping.get(bond_type, 6)


def mol_to_graph(mol: Chem.Mol) -> Data:
    atoms: List[List[int]] = []
    for atom in mol.GetAtoms():
        atomic_number = min(max(atom.GetAtomicNum(), 0), 118)
        chirality = min(max(int(atom.GetChiralTag()), 0), 15)
        degree = min(max(atom.GetTotalDegree(), 0), 12)
        formal_charge = min(max(atom.GetFormalCharge(), -6), 6) + 6
        total_h = min(max(atom.GetTotalNumHs(includeNeighbors=True), 0), 9)
        hybrid = min(max(int(atom.GetHybridization()), 0), 15)
        atoms.append(
            [
                atomic_number,
                chirality,
                degree,
                formal_charge,
                total_h,
                hybrid,
                int(atom.GetIsAromatic()),
                int(atom.IsInRing()),
            ]
        )

    edges: List[List[int]] = []
    edge_features: List[List[int]] = []
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        feature = [
            bond_type_index(bond),
            min(max(int(bond.GetStereo()), 0), 7),
            int(bond.GetIsConjugated()),
            int(bond.IsInRing()),
        ]
        edges.extend([[i, j], [j, i]])
        edge_features.extend([feature, feature])

    x = torch.tensor(atoms, dtype=torch.long)
    if edges:
        edge_index = torch.tensor(edges, dtype=torch.long).t().contiguous()
        edge_attr = torch.tensor(edge_features, dtype=torch.long)
    else:
        edge_index = torch.empty((2, 0), dtype=torch.long)
        edge_attr = torch.empty((0, len(BOND_DIMS)), dtype=torch.long)
    return Data(x=x, edge_index=edge_index, edge_attr=edge_attr)


def make_strata(
    metals: Sequence[str], targets: np.ndarray, min_count: int = 10
) -> np.ndarray:
    # Joint stratification by metal and coarse gap quantile, with rare strata folded.
    q = pd.qcut(targets[:, 0], q=5, labels=False, duplicates="drop")
    raw = np.asarray([f"{m}|{int(b)}" for m, b in zip(metals, q)], dtype=object)
    counts = Counter(raw.tolist())
    folded = np.asarray(
        [x if counts[x] >= min_count else str(x).split("|")[0] for x in raw],
        dtype=object,
    )
    folded_counts = Counter(folded.tolist())
    return np.asarray(
        [x if folded_counts[x] >= min_count else "OTHER" for x in folded],
        dtype=object,
    )


def make_group_split(smiles: Sequence[str], split_seed: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    indices = np.arange(len(smiles))
    first = GroupShuffleSplit(n_splits=1, test_size=0.20, random_state=split_seed)
    train_idx, held_idx = next(first.split(indices, groups=np.asarray(smiles, dtype=object)))
    held_smiles = np.asarray(smiles, dtype=object)[held_idx]
    second = GroupShuffleSplit(n_splits=1, test_size=0.50, random_state=split_seed + 1)
    val_local, test_local = next(
        second.split(np.arange(len(held_idx)), groups=held_smiles)
    )
    return (
        np.asarray(train_idx),
        np.asarray(held_idx[val_local]),
        np.asarray(held_idx[test_local]),
    )


def prepare_data(args: argparse.Namespace) -> PreparedData:
    project_root = args.project_root.resolve()
    property_path = project_root / "data" / "processed_full" / "dataset_coordsmiles.csv"
    record_path = project_root / "pipeline_full_output" / "results.jsonl"
    if not property_path.exists() or not record_path.exists():
        raise FileNotFoundError(f"Missing inputs: {property_path} or {record_path}")

    records, pipeline_audit = load_current_records(record_path)
    frame = pd.read_csv(property_path)
    raw_property_rows = len(frame)
    exact_duplicate_rows = int(frame.duplicated().sum())
    duplicate_id_rows = int(frame.duplicated("csd_code").sum())
    frame = frame.drop_duplicates("csd_code", keep="first").copy()
    frame = frame[frame["csd_code"].isin(records)].copy()
    frame = frame[frame["original_smiles"].notna()].copy()
    unknown_targets = sorted(set(args.targets) - set(TARGET_SPECS))
    if unknown_targets:
        raise ValueError(f"Unknown targets: {unknown_targets}; choose from {sorted(TARGET_SPECS)}")
    for target in args.targets:
        source_column, _ = TARGET_SPECS[target]
        frame = frame[np.isfinite(pd.to_numeric(frame[source_column], errors="coerce"))]

    ids: List[str] = []
    smiles: List[str] = []
    graphs: List[Data] = []
    field_dicts: List[Dict[str, float]] = []
    target_rows: List[List[float]] = []
    metals: List[str] = []
    parse_failures: List[str] = []

    for row in frame.itertuples(index=False):
        mol_id = str(row.csd_code)
        smi = str(row.original_smiles)
        mol = Chem.MolFromSmiles(smi)
        if mol is None or mol.GetNumAtoms() == 0:
            parse_failures.append(mol_id)
            continue
        record = records[mol_id]
        ids.append(mol_id)
        smiles.append(Chem.MolToSmiles(mol, canonical=True))
        graphs.append(mol_to_graph(mol))
        field_dicts.append(record_to_field_dict(record))
        target_rows.append(
            [
                float(getattr(row, TARGET_SPECS[target][0])) * TARGET_SPECS[target][1]
                for target in args.targets
            ]
        )
        metals.append(str(record.get("metal") or "UNK"))

    targets = np.asarray(target_rows, dtype=np.float32)
    train_idx, val_idx, test_idx = make_group_split(smiles, args.split_seed)

    if args.max_records and args.max_records < len(ids):
        # Development-only balanced subsample.  Re-split after selection.
        rng = np.random.default_rng(args.split_seed)
        chosen = np.sort(rng.choice(len(ids), size=args.max_records, replace=False))
        ids = [ids[i] for i in chosen]
        smiles = [smiles[i] for i in chosen]
        graphs = [graphs[i] for i in chosen]
        field_dicts = [field_dicts[i] for i in chosen]
        targets = targets[chosen]
        metals = [metals[i] for i in chosen]
        train_idx, val_idx, test_idx = make_group_split(smiles, args.split_seed)

    vectorizer = DictVectorizer(sparse=False, sort=True)
    vectorizer.fit([field_dicts[i] for i in train_idx])
    field_raw = vectorizer.transform(field_dicts).astype(np.float32)
    field_scaler = StandardScaler()
    field_scaler.fit(field_raw[train_idx])
    fields = field_scaler.transform(field_raw).astype(np.float32)

    target_mean = targets[train_idx].mean(axis=0)
    target_scale = targets[train_idx].std(axis=0)
    target_scale[target_scale < 1.0e-12] = 1.0

    y_scaled = (targets - target_mean) / target_scale
    for i, graph in enumerate(graphs):
        graph.coord_x = torch.from_numpy(fields[i : i + 1])
        graph.y = torch.from_numpy(y_scaled[i : i + 1].astype(np.float32))
        graph.sample_index = torch.tensor([i], dtype=torch.long)

    split_sets = {
        "train": set(ids[i] for i in train_idx),
        "validation": set(ids[i] for i in val_idx),
        "test": set(ids[i] for i in test_idx),
    }
    intersections = {
        "train_validation": len(split_sets["train"] & split_sets["validation"]),
        "train_test": len(split_sets["train"] & split_sets["test"]),
        "validation_test": len(split_sets["validation"] & split_sets["test"]),
    }
    canonical_smiles_split_overlap = {
        "train_validation": len(
            set(smiles[i] for i in train_idx) & set(smiles[i] for i in val_idx)
        ),
        "train_test": len(set(smiles[i] for i in train_idx) & set(smiles[i] for i in test_idx)),
        "validation_test": len(
            set(smiles[i] for i in val_idx) & set(smiles[i] for i in test_idx)
        ),
    }
    audit: Dict[str, object] = {
        **pipeline_audit,
        "property_rows_raw": raw_property_rows,
        "property_exact_duplicate_rows": exact_duplicate_rows,
        "property_duplicate_id_rows": duplicate_id_rows,
        "property_unique_ids": int(raw_property_rows - duplicate_id_rows),
        "joined_with_current_coordrep": int(len(frame)),
        "rdkit_parse_failures": len(parse_failures),
        "final_records": len(ids),
        "train_records": len(train_idx),
        "validation_records": len(val_idx),
        "test_records": len(test_idx),
        "split_id_intersections": intersections,
        "canonical_smiles_overlap": canonical_smiles_split_overlap,
        "field_dimension": int(fields.shape[1]),
        "target_names": list(args.targets),
        "target_train_mean": target_mean.tolist(),
        "target_train_scale": target_scale.tolist(),
        "metal_counts": dict(Counter(metals)),
        "split_seed": args.split_seed,
        "split_strategy": "80/10/10 GroupShuffleSplit by canonical RDKit SMILES; exact graph duplicates cannot cross splits",
        "data_sha256": {
            "properties": file_sha256(property_path),
            "coordrep_records": file_sha256(record_path),
        },
    }
    return PreparedData(
        ids=ids,
        smiles=smiles,
        graphs=graphs,
        fields=fields,
        targets=targets,
        target_names=list(args.targets),
        field_names=list(vectorizer.get_feature_names_out()),
        train_idx=np.asarray(train_idx),
        val_idx=np.asarray(val_idx),
        test_idx=np.asarray(test_idx),
        target_mean=target_mean,
        target_scale=target_scale,
        field_scaler_mean=field_scaler.mean_,
        field_scaler_scale=field_scaler.scale_,
        audit=audit,
    )


class CategoricalEncoder(nn.Module):
    def __init__(self, dimensions: Sequence[int], hidden: int):
        super().__init__()
        self.embeddings = nn.ModuleList([nn.Embedding(size, hidden) for size in dimensions])
        for embedding in self.embeddings:
            nn.init.xavier_uniform_(embedding.weight)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        out = self.embeddings[0](x[:, 0])
        for column, embedding in enumerate(self.embeddings[1:], start=1):
            out = out + embedding(x[:, column])
        return out


class CoordRepMLP(nn.Module):
    def __init__(self, field_dim: int, hidden: int, n_targets: int):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(field_dim, hidden),
            nn.LayerNorm(hidden),
            nn.SiLU(),
            nn.Dropout(0.10),
            nn.Linear(hidden, hidden),
            nn.SiLU(),
            nn.Dropout(0.10),
            nn.Linear(hidden, n_targets),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.net(x)


class GINRegressor(nn.Module):
    def __init__(
        self,
        field_dim: int,
        hidden: int,
        layers: int,
        n_targets: int,
        use_coordrep: bool,
    ):
        super().__init__()
        self.use_coordrep = use_coordrep
        self.atom_encoder = CategoricalEncoder(ATOM_DIMS, hidden)
        self.bond_encoder = CategoricalEncoder(BOND_DIMS, hidden)
        self.convs = nn.ModuleList()
        self.norms = nn.ModuleList()
        for _ in range(layers):
            mlp = nn.Sequential(
                nn.Linear(hidden, hidden * 2),
                nn.BatchNorm1d(hidden * 2),
                nn.SiLU(),
                nn.Linear(hidden * 2, hidden),
            )
            self.convs.append(GINEConv(mlp, train_eps=True))
            self.norms.append(nn.BatchNorm1d(hidden))
        if use_coordrep:
            self.coord_encoder = nn.Sequential(
                nn.Linear(field_dim, hidden),
                nn.LayerNorm(hidden),
                nn.SiLU(),
                nn.Dropout(0.10),
            )
            head_in = hidden * 2
        else:
            self.coord_encoder = None
            head_in = hidden
        self.head = nn.Sequential(
            nn.Linear(head_in, hidden),
            nn.SiLU(),
            nn.Dropout(0.10),
            nn.Linear(hidden, n_targets),
        )

    def forward(self, batch: Data) -> torch.Tensor:
        h = self.atom_encoder(batch.x)
        edge_h = self.bond_encoder(batch.edge_attr)
        for conv, norm in zip(self.convs, self.norms):
            residual = h
            h = conv(h, batch.edge_index, edge_h)
            h = norm(h)
            h = F.silu(h)
            h = F.dropout(h, p=0.10, training=self.training)
            h = h + residual
        pooled = global_mean_pool(h, batch.batch)
        if self.use_coordrep:
            pooled = torch.cat([pooled, self.coord_encoder(batch.coord_x)], dim=-1)
        return self.head(pooled)


def parameter_count(model: nn.Module) -> int:
    return sum(parameter.numel() for parameter in model.parameters() if parameter.requires_grad)


@torch.no_grad()
def predict_coordrep(
    model: nn.Module, loader: TensorDataLoader, device: torch.device
) -> Tuple[np.ndarray, np.ndarray]:
    model.eval()
    predictions: List[np.ndarray] = []
    labels: List[np.ndarray] = []
    for x, y in loader:
        predictions.append(model(x.to(device)).cpu().numpy())
        labels.append(y.numpy())
    return np.concatenate(predictions), np.concatenate(labels)


@torch.no_grad()
def predict_graph(
    model: nn.Module, loader: GraphDataLoader, device: torch.device
) -> Tuple[np.ndarray, np.ndarray]:
    model.eval()
    predictions: List[np.ndarray] = []
    labels: List[np.ndarray] = []
    for batch in loader:
        batch = batch.to(device)
        predictions.append(model(batch).cpu().numpy())
        labels.append(batch.y.cpu().numpy())
    return np.concatenate(predictions), np.concatenate(labels)


def normalized_validation_mae(pred: np.ndarray, truth: np.ndarray) -> float:
    return float(np.mean(np.abs(pred - truth)))


def train_coordrep_model(
    data: PreparedData, args: argparse.Namespace, seed: int, device: torch.device
) -> Tuple[nn.Module, Dict[str, object]]:
    set_seed(seed)
    x = torch.from_numpy(data.fields)
    y = torch.from_numpy(((data.targets - data.target_mean) / data.target_scale).astype(np.float32))
    generator = torch.Generator().manual_seed(seed)
    train_loader = TensorDataLoader(
        TensorDataset(x[data.train_idx], y[data.train_idx]),
        batch_size=args.batch_size,
        shuffle=True,
        generator=generator,
    )
    val_loader = TensorDataLoader(
        TensorDataset(x[data.val_idx], y[data.val_idx]),
        batch_size=args.batch_size * 2,
        shuffle=False,
    )
    model = CoordRepMLP(data.fields.shape[1], args.hidden, len(data.target_names)).to(device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    best_state = None
    best_epoch = 0
    best_val = math.inf
    remaining = args.patience
    history = []
    for epoch in range(1, args.epochs + 1):
        model.train()
        losses = []
        for batch_x, batch_y in train_loader:
            optimizer.zero_grad(set_to_none=True)
            loss = F.mse_loss(model(batch_x.to(device)), batch_y.to(device))
            loss.backward()
            optimizer.step()
            losses.append(float(loss.detach().cpu()))
        val_pred, val_truth = predict_coordrep(model, val_loader, device)
        val_mae = normalized_validation_mae(val_pred, val_truth)
        history.append({"epoch": epoch, "loss": float(np.mean(losses)), "val_norm_mae": val_mae})
        if val_mae < best_val - 1.0e-5:
            best_val = val_mae
            best_epoch = epoch
            best_state = {key: value.detach().cpu().clone() for key, value in model.state_dict().items()}
            remaining = args.patience
        else:
            remaining -= 1
            if remaining <= 0:
                break
    if best_state is None:
        raise RuntimeError("No CoordRep MLP checkpoint was selected")
    model.load_state_dict(best_state)
    return model, {
        "best_epoch": best_epoch,
        "best_validation_normalized_mae": best_val,
        "history": history,
        "parameters": parameter_count(model),
    }


def graph_loader(
    graphs: Sequence[Data],
    indices: np.ndarray,
    batch_size: int,
    shuffle: bool,
    seed: int,
    num_workers: int,
) -> GraphDataLoader:
    subset = [graphs[int(i)] for i in indices]
    generator = torch.Generator().manual_seed(seed)
    return GraphDataLoader(
        subset,
        batch_size=batch_size,
        shuffle=shuffle,
        generator=generator,
        num_workers=num_workers,
        pin_memory=torch.cuda.is_available(),
    )


def train_graph_model(
    model_name: str,
    data: PreparedData,
    args: argparse.Namespace,
    seed: int,
    device: torch.device,
) -> Tuple[nn.Module, Dict[str, object]]:
    set_seed(seed)
    train_loader = graph_loader(
        data.graphs, data.train_idx, args.batch_size, True, seed, args.num_workers
    )
    val_loader = graph_loader(
        data.graphs, data.val_idx, args.batch_size * 2, False, seed, args.num_workers
    )
    model_hidden = args.hidden + 5 if model_name == "gin_wide" else args.hidden
    model = GINRegressor(
        field_dim=data.fields.shape[1],
        hidden=model_hidden,
        layers=args.layers,
        n_targets=len(data.target_names),
        use_coordrep=model_name == "hybrid",
    ).to(device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    best_state = None
    best_epoch = 0
    best_val = math.inf
    remaining = args.patience
    history = []
    for epoch in range(1, args.epochs + 1):
        model.train()
        losses = []
        for batch in train_loader:
            batch = batch.to(device)
            optimizer.zero_grad(set_to_none=True)
            loss = F.mse_loss(model(batch), batch.y)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            optimizer.step()
            losses.append(float(loss.detach().cpu()))
        val_pred, val_truth = predict_graph(model, val_loader, device)
        val_mae = normalized_validation_mae(val_pred, val_truth)
        history.append({"epoch": epoch, "loss": float(np.mean(losses)), "val_norm_mae": val_mae})
        if val_mae < best_val - 1.0e-5:
            best_val = val_mae
            best_epoch = epoch
            best_state = {key: value.detach().cpu().clone() for key, value in model.state_dict().items()}
            remaining = args.patience
        else:
            remaining -= 1
            if remaining <= 0:
                break
    if best_state is None:
        raise RuntimeError(f"No {model_name} checkpoint was selected")
    model.load_state_dict(best_state)
    return model, {
        "best_epoch": best_epoch,
        "best_validation_normalized_mae": best_val,
        "history": history,
        "parameters": parameter_count(model),
    }


def metric_rows(
    model_name: str,
    seed: int,
    prediction_scaled: np.ndarray,
    truth_scaled: np.ndarray,
    data: PreparedData,
    metadata: Dict[str, object],
) -> Tuple[List[Dict[str, object]], np.ndarray, np.ndarray]:
    prediction = prediction_scaled * data.target_scale + data.target_mean
    truth = truth_scaled * data.target_scale + data.target_mean
    rows = []
    for column, target in enumerate(data.target_names):
        y_true = truth[:, column]
        y_pred = prediction[:, column]
        rows.append(
            {
                "model": model_name,
                "seed": seed,
                "target": target,
                "n_test": len(y_true),
                "mae": float(mean_absolute_error(y_true, y_pred)),
                "rmse": float(mean_squared_error(y_true, y_pred) ** 0.5),
                "r2": float(r2_score(y_true, y_pred)),
                "best_epoch": metadata["best_epoch"],
                "parameters": metadata["parameters"],
            }
        )
    return rows, prediction, truth


def aggregate(rows: Sequence[Dict[str, object]]) -> List[Dict[str, object]]:
    frame = pd.DataFrame(rows)
    result: List[Dict[str, object]] = []
    for (model, target), group in frame.groupby(["model", "target"], sort=True):
        entry: Dict[str, object] = {
            "model": model,
            "target": target,
            "n_test": int(group["n_test"].iloc[0]),
            "n_seeds": int(group["seed"].nunique()),
            "parameters": int(group["parameters"].iloc[0]),
        }
        for metric in ("mae", "rmse", "r2"):
            values = group[metric].astype(float).to_numpy()
            entry[f"{metric}_mean"] = float(values.mean())
            entry[f"{metric}_sd"] = float(values.std(ddof=1)) if len(values) > 1 else 0.0
        result.append(entry)
    return result


def write_csv(path: Path, rows: Sequence[Dict[str, object]]) -> None:
    if not rows:
        return
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    if args.smoke:
        args.max_records = args.max_records or 1800
        args.epochs = min(args.epochs, 2)
        args.patience = min(args.patience, 2)
        args.seeds = args.seeds[:1]
    args.out.mkdir(parents=True, exist_ok=True)
    started = time.time()
    device = torch.device(args.device)
    print(f"Preparing matched dataset on {device} ...", flush=True)
    data = prepare_data(args)
    with (args.out / "dataset_audit.json").open("w", encoding="utf-8") as handle:
        json.dump(data.audit, handle, indent=2)
    with (args.out / "field_names.json").open("w", encoding="utf-8") as handle:
        json.dump(data.field_names, handle, indent=2)
    split_rows = []
    split_lookup = {int(i): "train" for i in data.train_idx}
    split_lookup.update({int(i): "validation" for i in data.val_idx})
    split_lookup.update({int(i): "test" for i in data.test_idx})
    for i, mol_id in enumerate(data.ids):
        split_rows.append({"mol_id": mol_id, "split": split_lookup[i]})
    write_csv(args.out / "split_ids.csv", split_rows)
    print(
        f"Prepared {len(data.ids):,} records; train/val/test = "
        f"{len(data.train_idx):,}/{len(data.val_idx):,}/{len(data.test_idx):,}; "
        f"CoordRep fields = {data.fields.shape[1]}",
        flush=True,
    )

    all_rows: List[Dict[str, object]] = []
    prediction_rows: List[Dict[str, object]] = []
    run_metadata: Dict[str, object] = {}
    for model_name in args.models:
        if model_name not in {"coordrep_mlp", "gin", "gin_wide", "hybrid"}:
            raise ValueError(f"Unsupported model: {model_name}")
        for seed in args.seeds:
            print(f"Training {model_name}, seed={seed} ...", flush=True)
            run_start = time.time()
            if model_name == "coordrep_mlp":
                model, metadata = train_coordrep_model(data, args, seed, device)
                x = torch.from_numpy(data.fields[data.test_idx])
                y = torch.from_numpy(
                    ((data.targets[data.test_idx] - data.target_mean) / data.target_scale).astype(np.float32)
                )
                loader = TensorDataLoader(
                    TensorDataset(x, y), batch_size=args.batch_size * 2, shuffle=False
                )
                pred_scaled, truth_scaled = predict_coordrep(model, loader, device)
            else:
                model, metadata = train_graph_model(model_name, data, args, seed, device)
                loader = graph_loader(
                    data.graphs,
                    data.test_idx,
                    args.batch_size * 2,
                    False,
                    seed,
                    args.num_workers,
                )
                pred_scaled, truth_scaled = predict_graph(model, loader, device)
            rows, prediction, truth = metric_rows(
                model_name, seed, pred_scaled, truth_scaled, data, metadata
            )
            all_rows.extend(rows)
            for local_idx, global_idx in enumerate(data.test_idx):
                row: Dict[str, object] = {
                    "model": model_name,
                    "seed": seed,
                    "mol_id": data.ids[int(global_idx)],
                }
                for column, target in enumerate(data.target_names):
                    row[f"{target}_true"] = float(truth[local_idx, column])
                    row[f"{target}_pred"] = float(prediction[local_idx, column])
                prediction_rows.append(row)
            metadata["elapsed_seconds"] = time.time() - run_start
            run_metadata[f"{model_name}_seed{seed}"] = metadata
            print(
                f"Finished {model_name}, seed={seed}, best_epoch={metadata['best_epoch']}, "
                f"elapsed={metadata['elapsed_seconds']:.1f}s",
                flush=True,
            )
            del model
            if torch.cuda.is_available():
                torch.cuda.empty_cache()

    summary = aggregate(all_rows)
    write_csv(args.out / "metrics_by_seed.csv", all_rows)
    write_csv(args.out / "metrics_summary.csv", summary)
    write_csv(args.out / "test_predictions.csv", prediction_rows)
    with (args.out / "run_metadata.json").open("w", encoding="utf-8") as handle:
        json.dump(
            {
                "arguments": vars(args) | {"project_root": str(args.project_root), "out": str(args.out)},
                "device": str(device),
                "torch_version": torch.__version__,
                "cuda_available": torch.cuda.is_available(),
                "total_elapsed_seconds": time.time() - started,
                "runs": run_metadata,
            },
            handle,
            indent=2,
        )
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
