#!/usr/bin/env python
"""Five-fold graph-grouped SchNet OOF runner for CoordStatePairs.

This runner is the 3D counterpart to the CoordStatePairs 2D/record models.  It
uses an optional shared manifest with ``mol_id``, ``graph_group``, and ``fold``
columns.  When no manifest is supplied, the exact current 48,057-record cohort
is assigned to five folds by a deterministic hash of canonical attributed
RDKit SMILES; identical graphs therefore never cross folds.

The primary output is a long-form prediction table with merge keys compatible
with the planned 2D runner:

    model, fold, seed, mol_id, graph_group, canonical_graph,
    target, y_true, y_pred

The model receives atomic numbers and optimized Cartesian coordinates, while
the chemically attributed 2D graph/record inputs contain bond-order and charge
fields absent from this Cartesian interface. The inputs are therefore
non-nested, and this script must not be used to describe an equal-information
representation ablation or a state-of-the-art 3D leaderboard.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import time
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
from rdkit import Chem, RDLogger
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score

from analyze_graph_degenerate_states import (
    discrete_coordination_signature,
    load_records,
    shape_vectors_differ,
)
from run_schnet_3d_benchmark import (
    MultiTargetSchNet,
    make_loader,
    parameter_count,
    parse_selected_xyz,
    predict,
    set_seed,
    sha256,
    xyz_paths,
)


RDLogger.DisableLog("rdApp.*")

HARTREE_TO_EV = 27.211386245988
SCHEMA_VERSION = "coordstatepairs-oof-v1"
FROZEN_RECORDS_MANIFEST_SHA256 = (
    "c5f8f0a4dfb7cb9d06acff264774ad502732131345c4db6822eb110b324beb4e"
)
FROZEN_PAIRS_MANIFEST_SHA256 = (
    "2da2258128d1843ce40f085b54100cff0fc9b79c28b2318719df0398c33fe951"
)
MODEL_LABEL = "schnet_3d"
MODEL_ARCHITECTURE = "PyG SchNet distance-message-passing trunk"
MODEL_GUARDRAIL = (
    "SchNet receives atomic numbers and optimized Cartesian coordinates but omits bond "
    "order and formal/total charge attributes available to the chemically attributed 2D "
    "graph or CoordRep fields. This is a non-nested common-cohort comparison, not an "
    "equal-information ablation or a uniformly retuned 3D leaderboard."
)
TARGETS = {
    "hl_gap_ev": {
        "source": "hl_gap",
        "scale": HARTREE_TO_EV,
        "pair_threshold": 0.1,
        "pair_left": "left_hl_gap_ev",
        "pair_right": "right_hl_gap_ev",
        "pair_delta": "delta_hl_gap_ev",
        "pair_abs_delta": "abs_delta_hl_gap_ev",
    },
    "dipole_moment": {
        "source": "dipole_moment",
        "scale": 1.0,
        "pair_threshold": 1.0,
        "pair_left": "left_dipole_moment_d",
        "pair_right": "right_dipole_moment_d",
        "pair_delta": "delta_dipole_moment_d",
        "pair_abs_delta": "abs_delta_dipole_moment_d",
    },
}
PAIR_LEVELS = {
    "G0": "is_g0",
    "S1": "is_s1",
    "S2": "is_s2",
    "CONTEXT_EXACT": "is_context_exact",
    "S3": "is_s3",
}
PREDICTION_COLUMNS = [
    "schema_version",
    "model",
    "fold",
    "seed",
    "mol_id",
    "graph_group",
    "canonical_graph",
    "target",
    "y_true",
    "y_pred",
    "is_graph_degenerate",
    "continuous_only",
]
PAIR_COLUMNS = [
    "schema_version",
    "subset",
    "model",
    "fold",
    "seed",
    "target",
    "pair_id",
    "graph_group",
    "canonical_graph",
    "graph_group_size",
    "left_id",
    "right_id",
    "y_true_left",
    "y_true_right",
    "y_pred_left",
    "y_pred_right",
    "true_delta",
    "pred_delta",
    "true_abs_delta",
    "material_threshold",
    "abs_delta_error",
    "concordance",
    "prediction_tie",
    "prediction_tie_tolerance",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path("tmp/model_audit/data/CoordRep/CoordSMILES"),
    )
    parser.add_argument(
        "--cohort-csv",
        type=Path,
        default=Path("revision_experiments/results/downstream_graph_grouped/split_ids.csv"),
        help="ID source used only when --manifest is omitted.",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("revision_experiments/results/coordstate_pairs_manifest_v1/records.csv"),
        help=(
            "Optional shared CoordStatePairs manifest. Required logical fields are mol_id, "
            "graph_group, and fold; common aliases are accepted."
        ),
    )
    parser.add_argument(
        "--pairs-manifest",
        type=Path,
        default=Path("revision_experiments/results/coordstate_pairs_manifest_v1/pairs.csv"),
        help=(
            "Optional frozen CoordStatePairs pair table. If --manifest points to records.csv "
            "and a sibling pairs.csv exists, that table is selected automatically. The frozen "
            "G0/S1/S2/CONTEXT_EXACT/S3 flags define all primary pair subsets."
        ),
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("revision_experiments/results/coordstatepairs_schnet_3d_oof"),
    )
    parser.add_argument("--n-folds", type=int, default=5)
    parser.add_argument("--folds", nargs="+", type=int, default=None)
    parser.add_argument("--seeds", nargs="+", type=int, default=[11, 22, 33])
    parser.add_argument("--fold-salt", default="coordstatepairs-v1")
    parser.add_argument("--validation-fraction", type=float, default=0.125)
    parser.add_argument("--epochs", type=int, default=30)
    parser.add_argument("--patience", type=int, default=6)
    parser.add_argument("--batch-size", type=int, default=48)
    parser.add_argument("--hidden", type=int, default=128)
    parser.add_argument("--interactions", type=int, default=4)
    parser.add_argument("--num-gaussians", type=int, default=32)
    parser.add_argument("--cutoff", type=float, default=5.0)
    parser.add_argument(
        "--max-neighbors",
        type=int,
        default=96,
        help=(
            "Maximum neighbors retained inside the cutoff sphere. Production uses 96, "
            "which exceeds the audited maximum 5 A neighborhood in the frozen cohort."
        ),
    )
    parser.add_argument("--lr", type=float, default=1.0e-3)
    parser.add_argument("--weight-decay", type=float, default=1.0e-5)
    parser.add_argument("--num-workers", type=int, default=0)
    parser.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    parser.add_argument(
        "--precision",
        choices=["float16", "bfloat16", "float32"],
        default="bfloat16",
        help=(
            "CUDA training precision. Complete-radius production uses bfloat16 after an "
            "FP16 engineering smoke test triggered the non-finite-gradient guard."
        ),
    )
    parser.add_argument("--prediction-tie-tolerance", type=float, default=1.0e-5)
    parser.add_argument("--cshm-difference-tolerance", type=float, default=1.0e-3)
    parser.add_argument("--bootstrap-replicates", type=int, default=20_000)
    parser.add_argument("--bootstrap-seed", type=int, default=20260723)
    parser.add_argument(
        "--allow-nonfrozen-inputs",
        action="store_true",
        help=(
            "Permit a manifest/pairs hash other than the frozen CoordStatePairs v1 inputs. "
            "Never use this flag for the manuscript production runs."
        ),
    )
    parser.add_argument(
        "--smoke",
        action="store_true",
        help="Run one requested/default fold, one seed, and one epoch on the full fold cohort.",
    )
    return parser.parse_args()


def stable_hash(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def deterministic_fold(canonical_graph: str, n_folds: int, salt: str) -> int:
    value = int(stable_hash(f"{salt}|outer|{canonical_graph}")[:16], 16)
    return value % n_folds


def canonical_smiles(text: str) -> str | None:
    mol = Chem.MolFromSmiles(str(text))
    if mol is None or mol.GetNumAtoms() == 0:
        return None
    return Chem.MolToSmiles(mol, canonical=True)


def find_column(frame: pd.DataFrame, aliases: Sequence[str], required: bool = True) -> str | None:
    lower = {str(column).lower(): str(column) for column in frame.columns}
    for alias in aliases:
        if alias.lower() in lower:
            return lower[alias.lower()]
    if required:
        raise ValueError(f"Manifest is missing one of the required columns: {list(aliases)}")
    return None


def parse_fold_value(value: object) -> int:
    text = str(value).strip().lower()
    if text.startswith("fold_"):
        text = text[5:]
    elif text.startswith("fold"):
        text = text[4:]
    return int(float(text))


def load_property_frame(project_root: Path) -> Tuple[pd.DataFrame, Path]:
    path = project_root / "data" / "processed_full" / "dataset_coordsmiles.csv"
    frame = pd.read_csv(
        path,
        usecols=["csd_code", "original_smiles", "hl_gap", "dipole_moment"],
        dtype={"csd_code": str},
    )
    frame = frame.drop_duplicates("csd_code", keep="first").copy()
    return frame, path


def cohort_ids_from_csv(path: Path) -> List[str]:
    frame = pd.read_csv(path, dtype=str)
    id_column = find_column(frame, ["mol_id", "csd_code", "id"])
    ids = frame[id_column].astype(str).tolist()
    if len(ids) != len(set(ids)):
        raise ValueError("Cohort CSV contains duplicate IDs")
    return ids


def prepare_fold_manifest(args: argparse.Namespace) -> Tuple[pd.DataFrame, Path]:
    properties, property_path = load_property_frame(args.project_root)
    if args.manifest is not None:
        raw = pd.read_csv(args.manifest, dtype=str)
        id_column = find_column(raw, ["mol_id", "csd_code", "id", "record_id"])
        fold_column = find_column(raw, ["fold", "outer_fold", "oof_fold"])
        group_column = find_column(
            raw,
            [
                "graph_group",
                "graph_sha256",
                "g0_group",
                "canonical_graph_group",
                "group_id",
            ],
            required=False,
        )
        manifest = pd.DataFrame(
            {
                "mol_id": raw[id_column].astype(str),
                "fold": raw[fold_column].map(parse_fold_value).astype(int),
            }
        )
        if group_column is not None:
            manifest["provided_graph_group"] = raw[group_column].astype(str)
        fold_source = f"shared manifest: {args.manifest.resolve()}"
    else:
        ids = cohort_ids_from_csv(args.cohort_csv)
        manifest = pd.DataFrame({"mol_id": ids})
        fold_source = "deterministic SHA-256 hash of canonical attributed RDKit SMILES"

    if manifest["mol_id"].duplicated().any():
        raise ValueError("Fold manifest contains duplicate molecule IDs")
    merged = manifest.merge(
        properties,
        left_on="mol_id",
        right_on="csd_code",
        how="left",
        validate="one_to_one",
    )
    missing_property = merged["original_smiles"].isna()
    if missing_property.any():
        raise ValueError(f"Missing property/SMILES rows for {int(missing_property.sum())} manifest IDs")
    merged["canonical_graph"] = merged["original_smiles"].map(canonical_smiles)
    if merged["canonical_graph"].isna().any():
        raise ValueError(
            f"RDKit failed on {int(merged['canonical_graph'].isna().sum())} manifest structures"
        )
    computed_group = merged["canonical_graph"].map(lambda x: f"g0_{stable_hash(x)}")
    if "provided_graph_group" in merged:
        merged["graph_group"] = merged["provided_graph_group"].astype(str)
        canonical_to_group = merged.groupby("canonical_graph")["graph_group"].nunique()
        group_to_canonical = merged.groupby("graph_group")["canonical_graph"].nunique()
        if int((canonical_to_group != 1).sum()) or int((group_to_canonical != 1).sum()):
            raise ValueError("Provided graph_group is not one-to-one with canonical attributed graph")
    else:
        merged["graph_group"] = computed_group

    if "fold" not in merged:
        merged["fold"] = merged["canonical_graph"].map(
            lambda x: deterministic_fold(x, args.n_folds, args.fold_salt)
        )
    invalid_folds = sorted(set(merged["fold"].astype(int)) - set(range(args.n_folds)))
    if invalid_folds:
        raise ValueError(f"Manifest folds outside 0..{args.n_folds - 1}: {invalid_folds}")
    merged["fold"] = merged["fold"].astype(int)
    split_graphs = merged.groupby("graph_group")["fold"].nunique()
    if int((split_graphs != 1).sum()):
        raise ValueError("At least one exact graph_group is split across outer folds")

    merged["hl_gap_ev"] = pd.to_numeric(merged["hl_gap"], errors="coerce") * HARTREE_TO_EV
    merged["dipole_moment"] = pd.to_numeric(merged["dipole_moment"], errors="coerce")
    if not np.isfinite(merged[["hl_gap_ev", "dipole_moment"]].to_numpy(dtype=float)).all():
        raise ValueError("Non-finite target in OOF cohort")
    merged["fold_source"] = fold_source
    return merged.reset_index(drop=True), property_path


def deterministic_validation_indices(
    manifest: pd.DataFrame,
    test_fold: int,
    fraction: float,
    salt: str,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    if not 0.0 < fraction < 0.5:
        raise ValueError("--validation-fraction must be between 0 and 0.5")
    test_mask = manifest["fold"].to_numpy() == test_fold
    group_score: Dict[str, float] = {}
    for group in manifest.loc[~test_mask, "graph_group"].unique():
        raw = int(stable_hash(f"{salt}|validation|fold={test_fold}|{group}")[:16], 16)
        group_score[str(group)] = raw / float(16**16 - 1)
    val_groups = {group for group, score in group_score.items() if score < fraction}
    val_mask = (~test_mask) & manifest["graph_group"].isin(val_groups).to_numpy()
    train_mask = ~(test_mask | val_mask)
    train_idx = np.flatnonzero(train_mask)
    val_idx = np.flatnonzero(val_mask)
    test_idx = np.flatnonzero(test_mask)
    if not len(train_idx) or not len(val_idx) or not len(test_idx):
        raise ValueError(
            f"Empty split for fold {test_fold}: train={len(train_idx)}, val={len(val_idx)}, "
            f"test={len(test_idx)}"
        )
    return train_idx, val_idx, test_idx


def annotate_state_subsets(
    manifest: pd.DataFrame,
    records: Dict[str, dict],
    cshm_tolerance: float,
) -> Tuple[pd.DataFrame, Dict[str, object]]:
    annotated = manifest.copy()
    group_sizes = annotated.groupby("graph_group")["mol_id"].transform("size")
    annotated["is_graph_degenerate"] = group_sizes >= 2
    continuous_groups: set[str] = set()
    discrete_groups = 0
    for group_id, group in annotated.loc[annotated["is_graph_degenerate"]].groupby(
        "graph_group", sort=False
    ):
        group_records = [records[mol_id] for mol_id in group["mol_id"]]
        signatures = {json.dumps(discrete_coordination_signature(row), default=str) for row in group_records}
        if len(signatures) > 1:
            discrete_groups += 1
        if len(signatures) == 1 and shape_vectors_differ(group_records, cshm_tolerance):
            continuous_groups.add(str(group_id))
    annotated["continuous_only"] = annotated["graph_group"].isin(continuous_groups)
    audit = {
        "graph_degenerate_groups": int(
            annotated.loc[annotated["is_graph_degenerate"], "graph_group"].nunique()
        ),
        "graph_degenerate_members": int(annotated["is_graph_degenerate"].sum()),
        "groups_with_multiple_discrete_coordrep_state": discrete_groups,
        "continuous_only_groups": len(continuous_groups),
        "continuous_only_members": int(annotated["continuous_only"].sum()),
    }
    return annotated, audit


def resolve_pairs_manifest_path(args: argparse.Namespace) -> Path | None:
    """Resolve an explicit pair table or the conventional sibling of records.csv."""
    if args.pairs_manifest is not None:
        return args.pairs_manifest.resolve()
    if args.manifest is not None:
        sibling = args.manifest.resolve().with_name("pairs.csv")
        if sibling.exists():
            return sibling
    return None


def fallback_g0_pairs(manifest: pd.DataFrame) -> pd.DataFrame:
    """Construct G0 pairs only when a frozen hierarchy table is unavailable."""
    rows: List[dict] = []
    for graph_group, group in manifest.groupby("graph_group", sort=True):
        members = group.sort_values("mol_id").to_dict("records")
        if len(members) < 2:
            continue
        for left_index in range(len(members) - 1):
            for right_index in range(left_index + 1, len(members)):
                left = members[left_index]
                right = members[right_index]
                gap_delta = float(left["hl_gap_ev"] - right["hl_gap_ev"])
                dipole_delta = float(left["dipole_moment"] - right["dipole_moment"])
                rows.append(
                    {
                        "pair_id": stable_hash(
                            f"{SCHEMA_VERSION}|{graph_group}|{left['mol_id']}|{right['mol_id']}"
                        ),
                        "graph_sha256": str(graph_group),
                        "graph_group_size": len(members),
                        "fold": int(left["fold"]),
                        "left_mol_id": str(left["mol_id"]),
                        "right_mol_id": str(right["mol_id"]),
                        "is_g0": 1,
                        "is_s1": 0,
                        "is_s2": 0,
                        "is_context_exact": 0,
                        "is_s3": 0,
                        "left_hl_gap_ev": float(left["hl_gap_ev"]),
                        "right_hl_gap_ev": float(right["hl_gap_ev"]),
                        "delta_hl_gap_ev": gap_delta,
                        "abs_delta_hl_gap_ev": abs(gap_delta),
                        "left_dipole_moment_d": float(left["dipole_moment"]),
                        "right_dipole_moment_d": float(right["dipole_moment"]),
                        "delta_dipole_moment_d": dipole_delta,
                        "abs_delta_dipole_moment_d": abs(dipole_delta),
                    }
                )
    return pd.DataFrame(rows)


def load_pair_manifest(
    path: Path | None,
    manifest: pd.DataFrame,
) -> Tuple[pd.DataFrame, Dict[str, object]]:
    """Load and defensively validate the frozen pair hierarchy."""
    source = "derived G0-only pairs from the current canonical graph groups"
    if path is None:
        pairs = fallback_g0_pairs(manifest)
    else:
        pairs = pd.read_csv(path, dtype=str)
        source = f"frozen pair manifest: {path.resolve()}"

    required = {
        "pair_id",
        "graph_sha256",
        "graph_group_size",
        "fold",
        "left_mol_id",
        "right_mol_id",
        *PAIR_LEVELS.values(),
        *[
            str(TARGETS[target][field])
            for target in TARGETS
            for field in ("pair_left", "pair_right", "pair_delta", "pair_abs_delta")
        ],
    }
    missing_columns = sorted(required - set(pairs.columns))
    if missing_columns:
        raise ValueError(f"Pair manifest is missing required columns: {missing_columns}")
    pairs = pairs.copy()
    for column in ("pair_id", "graph_sha256", "left_mol_id", "right_mol_id"):
        pairs[column] = pairs[column].astype(str)
    pairs["fold"] = pairs["fold"].map(parse_fold_value).astype(int)
    pairs["graph_group_size"] = pd.to_numeric(
        pairs["graph_group_size"], errors="raise"
    ).astype(int)
    for flag in PAIR_LEVELS.values():
        pairs[flag] = pd.to_numeric(pairs[flag], errors="raise").astype(int)
        invalid = sorted(set(pairs[flag]) - {0, 1})
        if invalid:
            raise ValueError(f"Pair manifest flag {flag} contains values outside 0/1: {invalid}")
    numeric_target_columns = sorted(
        {
            str(TARGETS[target][field])
            for target in TARGETS
            for field in ("pair_left", "pair_right", "pair_delta", "pair_abs_delta")
        }
    )
    for column in numeric_target_columns:
        pairs[column] = pd.to_numeric(pairs[column], errors="raise")
    if not np.isfinite(pairs[numeric_target_columns].to_numpy(dtype=float)).all():
        raise ValueError("Pair manifest contains non-finite target values")
    if pairs["pair_id"].duplicated().any():
        raise ValueError("Pair manifest contains duplicate pair_id values")
    if (pairs["left_mol_id"] >= pairs["right_mol_id"]).any():
        raise ValueError("Pair manifest endpoints must use unique lexicographic left/right order")
    if (pairs["is_g0"] != 1).any():
        raise ValueError("Every pair-manifest row must be a G0 same-graph pair")
    for child, parent in (
        ("is_s1", "is_g0"),
        ("is_s2", "is_s1"),
        ("is_context_exact", "is_s2"),
        ("is_s3", "is_context_exact"),
    ):
        if ((pairs[child] == 1) & (pairs[parent] != 1)).any():
            raise ValueError(f"Pair hierarchy is not nested: {child} is not a subset of {parent}")

    record_index = manifest.set_index("mol_id", verify_integrity=True)
    missing_left = sorted(set(pairs["left_mol_id"]) - set(record_index.index))
    missing_right = sorted(set(pairs["right_mol_id"]) - set(record_index.index))
    if missing_left or missing_right:
        raise ValueError(
            f"Pair endpoints missing from records: left={len(missing_left)}, right={len(missing_right)}"
        )
    left_records = record_index.loc[pairs["left_mol_id"]].reset_index(drop=True)
    right_records = record_index.loc[pairs["right_mol_id"]].reset_index(drop=True)
    graph_mismatch = (
        (left_records["graph_group"].astype(str).to_numpy() != pairs["graph_sha256"].to_numpy())
        | (right_records["graph_group"].astype(str).to_numpy() != pairs["graph_sha256"].to_numpy())
    )
    fold_mismatch = (
        (left_records["fold"].to_numpy(dtype=int) != pairs["fold"].to_numpy(dtype=int))
        | (right_records["fold"].to_numpy(dtype=int) != pairs["fold"].to_numpy(dtype=int))
    )
    if graph_mismatch.any() or fold_mismatch.any():
        raise ValueError(
            "Pair manifest does not match record graph/fold assignments: "
            f"graph_mismatches={int(graph_mismatch.sum())}, "
            f"fold_mismatches={int(fold_mismatch.sum())}"
        )
    actual_group_sizes = manifest.groupby("graph_group")["mol_id"].size()
    expected_group_sizes = pairs["graph_sha256"].map(actual_group_sizes).to_numpy(dtype=int)
    if not np.array_equal(expected_group_sizes, pairs["graph_group_size"].to_numpy(dtype=int)):
        raise ValueError("Pair manifest graph_group_size disagrees with the records manifest")

    target_consistency_failures: Dict[str, int] = {}
    for target, details in TARGETS.items():
        left_values = left_records[target].to_numpy(dtype=float)
        right_values = right_records[target].to_numpy(dtype=float)
        pair_left = pairs[str(details["pair_left"])].to_numpy(dtype=float)
        pair_right = pairs[str(details["pair_right"])].to_numpy(dtype=float)
        pair_delta = pairs[str(details["pair_delta"])].to_numpy(dtype=float)
        pair_abs_delta = pairs[str(details["pair_abs_delta"])].to_numpy(dtype=float)
        failures = ~(
            np.isclose(left_values, pair_left, rtol=0.0, atol=1.0e-9)
            & np.isclose(right_values, pair_right, rtol=0.0, atol=1.0e-9)
            & np.isclose(pair_left - pair_right, pair_delta, rtol=0.0, atol=1.0e-9)
            & np.isclose(np.abs(pair_delta), pair_abs_delta, rtol=0.0, atol=1.0e-9)
        )
        target_consistency_failures[target] = int(failures.sum())
    if sum(target_consistency_failures.values()):
        raise ValueError(
            f"Pair target values disagree with records/orientation: {target_consistency_failures}"
        )

    audit = {
        "source": source,
        "rows": int(len(pairs)),
        "level_pair_counts": {
            level: int(pairs[flag].sum()) for level, flag in PAIR_LEVELS.items()
        },
        "material_pair_counts": {
            level: {
                target: int(
                    (
                        (pairs[flag] == 1)
                        & (
                            pairs[str(details["pair_abs_delta"])]
                            >= float(details["pair_threshold"])
                        )
                    ).sum()
                )
                for target, details in TARGETS.items()
            }
            for level, flag in PAIR_LEVELS.items()
        },
        "target_consistency_failures": target_consistency_failures,
        "graph_assignment_mismatches": int(graph_mismatch.sum()),
        "fold_assignment_mismatches": int(fold_mismatch.sum()),
    }
    return pairs.reset_index(drop=True), audit


def make_base_dataset(
    manifest: pd.DataFrame,
    structures: Dict[str, Tuple[torch.Tensor, torch.Tensor]],
) -> List[object]:
    from torch_geometric.data import Data

    dataset = []
    for row_id, mol_id in enumerate(manifest["mol_id"]):
        z, pos = structures[str(mol_id)]
        dataset.append(
            Data(
                z=z,
                pos=pos,
                y=torch.zeros((1, len(TARGETS)), dtype=torch.float32),
                row_id=torch.tensor([row_id], dtype=torch.long),
            )
        )
    return dataset


def complete_radius_graph_audit(
    dataset: Sequence[object], max_neighbors: int, cutoff: float
) -> Dict[str, object]:
    """Census every cutoff neighborhood and fail if the cap can truncate it."""

    if not dataset:
        raise RuntimeError("Cannot audit an empty Cartesian dataset")
    margin = 1.0e-6
    threshold = float(cutoff) + margin
    legacy_cap = 32
    maximum_atoms = 0
    maximum_inclusive_neighbors = 0
    records_above_legacy_cap = 0
    atoms_above_legacy_cap = 0
    total_inclusive_edges = 0
    omitted_inclusive_edges_at_legacy_cap = 0
    total_atoms = 0
    for data in dataset:
        positions = data.pos.detach().cpu().to(dtype=torch.float64)
        n_atoms = int(positions.shape[0])
        maximum_atoms = max(maximum_atoms, n_atoms)
        total_atoms += n_atoms
        counts = (torch.cdist(positions, positions) <= threshold).sum(dim=1)
        maximum_inclusive_neighbors = max(
            maximum_inclusive_neighbors, int(counts.max().item())
        )
        above = counts > legacy_cap
        if bool(above.any()):
            records_above_legacy_cap += 1
        atoms_above_legacy_cap += int(above.sum().item())
        total_inclusive_edges += int(counts.sum().item())
        omitted_inclusive_edges_at_legacy_cap += int(
            torch.clamp(counts - legacy_cap, min=0).sum().item()
        )
    minimum_complete_cap = maximum_inclusive_neighbors
    if int(max_neighbors) < minimum_complete_cap:
        raise RuntimeError(
            "Neighbor cap can truncate the cohort radius graph: "
            f"max_neighbors={max_neighbors}, but the conservative cutoff census "
            f"requires at least {minimum_complete_cap}."
        )
    return {
        "cutoff_angstrom": float(cutoff),
        "census_distance_margin_angstrom": margin,
        "cohort_records": int(len(dataset)),
        "cohort_atoms": int(total_atoms),
        "maximum_atoms_per_record": int(maximum_atoms),
        "maximum_neighbors_within_cutoff_including_self": int(
            maximum_inclusive_neighbors
        ),
        "minimum_complete_radius_max_neighbors": int(minimum_complete_cap),
        "configured_max_neighbors": int(max_neighbors),
        "untruncated_radius_graph_guaranteed": True,
        "legacy_cap": legacy_cap,
        "records_exceeding_legacy_cap": int(records_above_legacy_cap),
        "atoms_exceeding_legacy_cap": int(atoms_above_legacy_cap),
        "total_inclusive_directed_edges": int(total_inclusive_edges),
        "inclusive_directed_edges_omitted_at_legacy_cap": int(
            omitted_inclusive_edges_at_legacy_cap
        ),
        "proof": (
            "Every atom-level neighborhood was counted at cutoff + 1e-6 A, including "
            "the central atom as a conservative margin; configured max_neighbors is no "
            "smaller than the observed maximum."
        ),
    }


def _random_orthogonal_matrix(seed: int) -> torch.Tensor:
    generator = torch.Generator().manual_seed(seed)
    matrix = torch.randn(3, 3, generator=generator)
    q, _ = torch.linalg.qr(matrix)
    if torch.det(q) < 0:
        q[:, 0] *= -1
    return q


@torch.no_grad()
def invariance_preflight(
    out: Path,
    device: torch.device,
    args: argparse.Namespace,
    radius_audit: Dict[str, object],
) -> Dict[str, object]:
    """Regression test invariant scalar outputs on a cap-triggering dense graph."""

    set_seed(20260803)
    n_atoms = 48
    generator = torch.Generator().manual_seed(20260803)
    z_pattern = torch.tensor([1, 6, 7, 8, 15, 16, 17, 26], dtype=torch.long)
    z = z_pattern.repeat(math.ceil(n_atoms / len(z_pattern)))[:n_atoms].clone()
    # All pair distances are safely below 5 A.  With a cap of 32 this fixture
    # necessarily truncates; with the production cap it contains all 47
    # neighbors per atom and detects the historical ordering failure mode.
    pos = 0.35 * torch.randn(n_atoms, 3, generator=generator)
    maximum_pair_distance = float(torch.cdist(pos, pos).max())
    if maximum_pair_distance >= float(args.cutoff):
        raise RuntimeError("Dense preflight fixture unexpectedly exceeds the cutoff")
    if int(args.max_neighbors) < n_atoms - 1:
        raise RuntimeError(
            "Dense permutation preflight requires at least 47 neighbors per atom"
        )

    model = MultiTargetSchNet(
        hidden=int(args.hidden),
        interactions=int(args.interactions),
        num_gaussians=int(args.num_gaussians),
        cutoff=float(args.cutoff),
        max_neighbors=int(args.max_neighbors),
        n_targets=len(TARGETS),
        target_mean=np.asarray([2.9, 5.7], dtype=np.float32),
        target_scale=np.asarray([0.93, 4.14], dtype=np.float32),
    ).to(device)
    model.eval()

    from torch_geometric.data import Data

    def evaluate(z_value: torch.Tensor, pos_value: torch.Tensor) -> np.ndarray:
        data = Data(z=z_value, pos=pos_value)
        data.batch = torch.zeros(len(z_value), dtype=torch.long)
        return model(data.to(device)).detach().cpu().numpy().reshape(-1)

    base = evaluate(z, pos)
    rotation = _random_orthogonal_matrix(20260804)
    reflected = pos.clone()
    reflected[:, 0] *= -1
    transformed = {
        "proper_rotation": evaluate(z, pos @ rotation.T),
        "reflection": evaluate(z, reflected),
        "translation": evaluate(z, pos + torch.tensor([3.1, -2.7, 5.4])),
    }
    for permutation_index in range(10):
        permutation = torch.randperm(
            n_atoms,
            generator=torch.Generator().manual_seed(20260805 + permutation_index),
        )
        transformed[f"atom_permutation_dense_48_{permutation_index + 1:02d}"] = evaluate(
            z[permutation], pos[permutation]
        )
    differences = {
        name: float(np.max(np.abs(value - base))) for name, value in transformed.items()
    }
    tolerance = 1.0e-5
    audit = {
        "model": MODEL_LABEL,
        "architecture": MODEL_ARCHITECTURE,
        "device": str(device),
        "torch_version": torch.__version__,
        "torch_geometric_version": __import__("torch_geometric").__version__,
        "dense_fixture_atoms": n_atoms,
        "dense_fixture_other_atoms_within_cutoff": n_atoms - 1,
        "dense_fixture_conservative_slots_including_self": n_atoms,
        "dense_fixture_maximum_pair_distance_angstrom": maximum_pair_distance,
        "configured_max_neighbors": int(args.max_neighbors),
        "cohort_radius_graph_audit": radius_audit,
        "base_prediction": base.tolist(),
        "maximum_absolute_differences": differences,
        "tolerance": tolerance,
        "passed": bool(max(differences.values()) <= tolerance),
        "scope": (
            "Implementation regression under proper rotation, reflection, translation, "
            "and ten atom permutations on a dense graph with more than 32 neighbors per "
            "atom; this verifies invariant scalar outputs, "
            "not an empirical property-performance result."
        ),
    }
    (out / "INVARIANCE_PREFLIGHT.json").write_text(
        json.dumps(audit, indent=2, allow_nan=False), encoding="utf-8"
    )
    if not audit["passed"]:
        raise RuntimeError(f"{MODEL_LABEL} invariance preflight failed: {differences}")
    return audit


def write_model_card(
    out: Path,
    args: argparse.Namespace,
    preflight: Dict[str, object],
) -> None:
    source_paths = [Path(__file__).resolve(), Path(__file__).with_name("run_schnet_3d_benchmark.py").resolve()]
    card = {
        "model_label": MODEL_LABEL,
        "architecture": MODEL_ARCHITECTURE,
        "production_configuration": {
            "layers": int(args.interactions),
            "hidden_channels": int(args.hidden),
            "radial_basis_functions": int(args.num_gaussians),
            "cutoff_angstrom": float(args.cutoff),
            "max_neighbors": int(args.max_neighbors),
            "requested_training_precision": str(args.precision),
            "effective_training_precision": (
                str(args.precision)
                if str(preflight.get("device", "")).startswith("cuda")
                else "float32"
            ),
        },
        "inputs": ["atomic number", "optimized Cartesian coordinates in angstrom"],
        "comparison_role": MODEL_GUARDRAIL,
        "invariance_preflight": preflight,
        "source_sha256": {path.name: sha256(path) for path in source_paths},
    }
    (out / "MODEL_CARD.json").write_text(
        json.dumps(card, indent=2, allow_nan=False), encoding="utf-8"
    )


def set_scaled_targets(
    dataset: Sequence[object],
    targets: np.ndarray,
    mean: np.ndarray,
    scale: np.ndarray,
) -> None:
    scaled = ((targets - mean) / scale).astype(np.float32)
    for row_id, data in enumerate(dataset):
        data.y = torch.tensor(scaled[row_id : row_id + 1], dtype=torch.float32)


def train_fold_seed(
    dataset: Sequence[object],
    targets: np.ndarray,
    train_idx: np.ndarray,
    val_idx: np.ndarray,
    test_idx: np.ndarray,
    fold: int,
    seed: int,
    args: argparse.Namespace,
    device: torch.device,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, Dict[str, object]]:
    set_seed(seed)
    target_mean = targets[train_idx].mean(axis=0).astype(np.float32)
    target_scale = targets[train_idx].std(axis=0).astype(np.float32)
    target_scale[target_scale < 1.0e-12] = 1.0
    set_scaled_targets(dataset, targets, target_mean, target_scale)
    train_loader = make_loader(
        dataset, train_idx, args.batch_size, True, seed, args.num_workers
    )
    val_loader = make_loader(
        dataset, val_idx, args.batch_size * 2, False, seed, args.num_workers
    )
    test_loader = make_loader(
        dataset, test_idx, args.batch_size * 2, False, seed, args.num_workers
    )
    model = MultiTargetSchNet(
        hidden=args.hidden,
        interactions=args.interactions,
        num_gaussians=args.num_gaussians,
        cutoff=args.cutoff,
        max_neighbors=args.max_neighbors,
        n_targets=len(TARGETS),
        target_mean=target_mean,
        target_scale=target_scale,
    ).to(device)
    optimizer = torch.optim.AdamW(
        model.parameters(), lr=args.lr, weight_decay=args.weight_decay
    )
    use_amp = device.type == "cuda" and args.precision != "float32"
    amp_dtype = torch.bfloat16 if args.precision == "bfloat16" else torch.float16
    use_grad_scaler = use_amp and args.precision == "float16"
    scaler = torch.amp.GradScaler("cuda", enabled=use_grad_scaler)
    best_state = None
    best_epoch = 0
    best_val = math.inf
    patience_left = args.patience
    history: List[dict] = []
    started = time.time()
    for epoch in range(1, args.epochs + 1):
        model.train()
        batch_losses: List[float] = []
        for batch_number, batch in enumerate(train_loader, start=1):
            batch = batch.to(device, non_blocking=True)
            optimizer.zero_grad(set_to_none=True)
            with torch.amp.autocast("cuda", enabled=use_amp, dtype=amp_dtype):
                batch_prediction = model(batch)
                loss = F.mse_loss(batch_prediction, batch.y)
            if not bool(torch.isfinite(batch_prediction).all()):
                row_ids = batch.row_id.detach().cpu().view(-1).tolist()
                raise FloatingPointError(
                    f"Non-finite training prediction for fold={fold}, seed={seed}, "
                    f"epoch={epoch}, batch={batch_number}, row_ids={row_ids[:12]}"
                )
            if not bool(torch.isfinite(loss)):
                raise FloatingPointError(
                    f"Non-finite training loss for fold={fold}, seed={seed}, epoch={epoch}, "
                    f"batch={batch_number}"
                )
            scaler.scale(loss).backward()
            scaler.unscale_(optimizer)
            torch.nn.utils.clip_grad_norm_(
                model.parameters(), 5.0, error_if_nonfinite=True
            )
            scaler.step(optimizer)
            scaler.update()
            batch_losses.append(float(loss.detach().cpu()))
        val_pred, val_truth, _ = predict(model, val_loader, device)
        if not np.isfinite(val_pred).all() or not np.isfinite(val_truth).all():
            raise FloatingPointError(
                f"Non-finite validation values for fold={fold}, seed={seed}, epoch={epoch}"
            )
        val_norm_mae = float(np.mean(np.abs(val_pred - val_truth)))
        history.append(
            {
                "epoch": epoch,
                "train_loss": float(np.mean(batch_losses)),
                "validation_normalized_mae": val_norm_mae,
            }
        )
        print(
            f"fold={fold} seed={seed} epoch={epoch} loss={np.mean(batch_losses):.5f} "
            f"val_norm_mae={val_norm_mae:.5f}",
            flush=True,
        )
        if val_norm_mae < best_val - 1.0e-5:
            best_val = val_norm_mae
            best_epoch = epoch
            best_state = {
                key: value.detach().cpu().clone() for key, value in model.state_dict().items()
            }
            patience_left = args.patience
        else:
            patience_left -= 1
            if patience_left <= 0:
                break
    if best_state is None:
        raise RuntimeError(f"No checkpoint selected for fold={fold}, seed={seed}")
    model.load_state_dict(best_state)
    pred_scaled, truth_scaled, row_ids = predict(model, test_loader, device)
    if not np.isfinite(pred_scaled).all() or not np.isfinite(truth_scaled).all():
        raise FloatingPointError(f"Non-finite test values for fold={fold}, seed={seed}")
    predictions = pred_scaled * target_scale + target_mean
    truths = truth_scaled * target_scale + target_mean
    if not np.isfinite(predictions).all() or not np.isfinite(truths).all():
        raise FloatingPointError(
            f"Non-finite unscaled test values for fold={fold}, seed={seed}"
        )
    effective_precision = args.precision if device.type == "cuda" else "float32"
    metadata = {
        "fold": fold,
        "seed": seed,
        "train_records": int(len(train_idx)),
        "validation_records": int(len(val_idx)),
        "test_records": int(len(test_idx)),
        "best_epoch": best_epoch,
        "best_validation_normalized_mae": best_val,
        "parameters": parameter_count(model),
        "target_train_mean": target_mean.tolist(),
        "target_train_scale": target_scale.tolist(),
        "elapsed_seconds": time.time() - started,
        "requested_training_precision": args.precision,
        "effective_training_precision": effective_precision,
        "history": history,
    }
    del model
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
    return predictions, truths, row_ids.astype(int), metadata


def prediction_rows(
    manifest: pd.DataFrame,
    fold: int,
    seed: int,
    row_ids: np.ndarray,
    predictions: np.ndarray,
    truths: np.ndarray,
) -> List[dict]:
    rows: List[dict] = []
    target_names = list(TARGETS)
    for local_row, row_id in enumerate(row_ids):
        meta = manifest.iloc[int(row_id)]
        for target_column, target in enumerate(target_names):
            rows.append(
                {
                    "schema_version": SCHEMA_VERSION,
                    "model": MODEL_LABEL,
                    "fold": fold,
                    "seed": seed,
                    "mol_id": str(meta["mol_id"]),
                    "graph_group": str(meta["graph_group"]),
                    "canonical_graph": str(meta["canonical_graph"]),
                    "target": target,
                    "y_true": float(truths[local_row, target_column]),
                    "y_pred": float(predictions[local_row, target_column]),
                    "is_graph_degenerate": bool(meta["is_graph_degenerate"]),
                    "continuous_only": bool(meta["continuous_only"]),
                }
            )
    return rows


def append_csv(path: Path, rows: Sequence[dict], columns: Sequence[str]) -> None:
    first = not path.exists()
    with path.open("a", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(columns))
        if first:
            writer.writeheader()
        writer.writerows(rows)


def ordinary_metric_rows(
    prediction_frame: pd.DataFrame,
    run_metadata: Dict[Tuple[int, int], dict],
) -> pd.DataFrame:
    rows = []
    for (fold, seed, target), group in prediction_frame.groupby(
        ["fold", "seed", "target"], sort=True
    ):
        y_true = group["y_true"].to_numpy(dtype=float)
        y_pred = group["y_pred"].to_numpy(dtype=float)
        metadata = run_metadata[(int(fold), int(seed))]
        rows.append(
            {
                "model": MODEL_LABEL,
                "fold": int(fold),
                "seed": int(seed),
                "target": target,
                "n_test": len(group),
                "mae": float(mean_absolute_error(y_true, y_pred)),
                "rmse": float(mean_squared_error(y_true, y_pred) ** 0.5),
                "r2": float(r2_score(y_true, y_pred)),
                "spearman": float(pd.Series(y_true).corr(pd.Series(y_pred), method="spearman")),
                "best_epoch": metadata["best_epoch"],
                "parameters": metadata["parameters"],
            }
        )
    return pd.DataFrame(rows)


def pair_rows_for_run(
    run: pd.DataFrame,
    pair_manifest: pd.DataFrame,
    tie_tolerance: float,
) -> List[dict]:
    """Join one fold/seed/target prediction block to frozen oriented pairs."""
    fold = int(run["fold"].iloc[0])
    seed = int(run["seed"].iloc[0])
    target = str(run["target"].iloc[0])
    model = str(run["model"].iloc[0])
    details = TARGETS[target]
    threshold = float(details["pair_threshold"])
    eligible = pair_manifest[
        (pair_manifest["fold"] == fold)
        & (pair_manifest[str(details["pair_abs_delta"])] >= threshold)
    ].copy()
    if eligible.empty:
        return []

    prediction = run[
        ["mol_id", "graph_group", "canonical_graph", "y_true", "y_pred"]
    ].copy()
    left_prediction = prediction.rename(
        columns={
            "mol_id": "left_mol_id",
            "graph_group": "left_graph_group",
            "canonical_graph": "canonical_graph",
            "y_true": "y_true_left",
            "y_pred": "y_pred_left",
        }
    )
    right_prediction = prediction.rename(
        columns={
            "mol_id": "right_mol_id",
            "graph_group": "right_graph_group",
            "canonical_graph": "right_canonical_graph",
            "y_true": "y_true_right",
            "y_pred": "y_pred_right",
        }
    )
    joined = eligible.merge(
        left_prediction,
        on="left_mol_id",
        how="left",
        validate="many_to_one",
    ).merge(
        right_prediction,
        on="right_mol_id",
        how="left",
        validate="many_to_one",
    )
    missing_prediction = joined[["y_pred_left", "y_pred_right"]].isna().any(axis=1)
    if missing_prediction.any():
        raise RuntimeError(
            f"Frozen pairs miss OOF endpoint predictions for fold={fold}, seed={seed}, "
            f"target={target}: {int(missing_prediction.sum())} rows"
        )
    graph_mismatch = (
        (joined["left_graph_group"] != joined["graph_sha256"])
        | (joined["right_graph_group"] != joined["graph_sha256"])
    )
    if graph_mismatch.any():
        raise RuntimeError(f"Prediction/pair graph mismatch in {int(graph_mismatch.sum())} rows")
    frozen_left = joined[str(details["pair_left"])].to_numpy(dtype=float)
    frozen_right = joined[str(details["pair_right"])].to_numpy(dtype=float)
    predicted_truth_left = joined["y_true_left"].to_numpy(dtype=float)
    predicted_truth_right = joined["y_true_right"].to_numpy(dtype=float)
    truth_mismatch = ~(
        np.isclose(predicted_truth_left, frozen_left, rtol=1.0e-6, atol=2.0e-6)
        & np.isclose(predicted_truth_right, frozen_right, rtol=1.0e-6, atol=2.0e-6)
    )
    if truth_mismatch.any():
        raise RuntimeError(
            f"OOF y_true values disagree with frozen pair targets in {int(truth_mismatch.sum())} rows"
        )

    rows: List[dict] = []
    for item in joined.itertuples(index=False):
        frozen_true_delta = float(getattr(item, str(details["pair_delta"])))
        pred_delta = float(item.y_pred_left - item.y_pred_right)
        tied = abs(pred_delta) <= tie_tolerance
        concordance = (
            0.5
            if tied
            else float(math.copysign(1.0, pred_delta) == math.copysign(1.0, frozen_true_delta))
        )
        for subset, flag in PAIR_LEVELS.items():
            if int(getattr(item, flag)) != 1:
                continue
            rows.append(
                {
                    "schema_version": SCHEMA_VERSION,
                    "subset": subset,
                    "model": model,
                    "fold": fold,
                    "seed": seed,
                    "target": target,
                    "pair_id": str(item.pair_id),
                    "graph_group": str(item.graph_sha256),
                    "canonical_graph": str(item.canonical_graph),
                    "graph_group_size": int(item.graph_group_size),
                    "left_id": str(item.left_mol_id),
                    "right_id": str(item.right_mol_id),
                    "y_true_left": float(item.y_true_left),
                    "y_true_right": float(item.y_true_right),
                    "y_pred_left": float(item.y_pred_left),
                    "y_pred_right": float(item.y_pred_right),
                    "true_delta": frozen_true_delta,
                    "pred_delta": pred_delta,
                    "true_abs_delta": float(getattr(item, str(details["pair_abs_delta"]))),
                    "material_threshold": threshold,
                    "abs_delta_error": abs(pred_delta - frozen_true_delta),
                    "concordance": concordance,
                    "prediction_tie": tied,
                    "prediction_tie_tolerance": tie_tolerance,
                }
            )
    return rows


def pair_metric_block(block: pd.DataFrame) -> Dict[str, object]:
    return {
        "n_graph_groups": int(block["graph_group"].nunique()),
        "n_pairs": int(block["pair_id"].nunique()),
        "n_member_records": int(
            len(set(block["left_id"].astype(str)) | set(block["right_id"].astype(str)))
        ),
        "true_abs_delta_mean": float(block["true_abs_delta"].mean()),
        "pair_delta_mae": float(block["abs_delta_error"].mean()),
        "pair_concordance_ties_half": float(block["concordance"].mean()),
        "pair_tie_fraction": float(block["prediction_tie"].mean()),
    }


def pair_outputs(
    prediction_frame: pd.DataFrame,
    pair_manifest: pd.DataFrame,
    tie_tolerance: float,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    all_pair_rows: List[dict] = []
    for (_, _, _), run in prediction_frame.groupby(
        ["fold", "seed", "target"], sort=True
    ):
        all_pair_rows.extend(pair_rows_for_run(run, pair_manifest, tie_tolerance))
    pair_frame = pd.DataFrame(all_pair_rows, columns=PAIR_COLUMNS)

    fold_metric_rows: List[dict] = []
    for (subset, fold, seed, target), block in pair_frame.groupby(
        ["subset", "fold", "seed", "target"], sort=True
    ):
        fold_metric_rows.append(
            {
                "subset": subset,
                "model": MODEL_LABEL,
                "fold": int(fold),
                "seed": int(seed),
                "target": target,
                **pair_metric_block(block),
            }
        )

    seed_metric_rows: List[dict] = []
    for (subset, seed, target), block in pair_frame.groupby(
        ["subset", "seed", "target"], sort=True
    ):
        seed_metric_rows.append(
            {
                "subset": subset,
                "model": MODEL_LABEL,
                "seed": int(seed),
                "target": target,
                "n_folds": int(block["fold"].nunique()),
                **pair_metric_block(block),
            }
        )
    return pair_frame, pd.DataFrame(fold_metric_rows), pd.DataFrame(seed_metric_rows)


def bootstrap_pair_metrics(
    pair_frame: pd.DataFrame,
    replicates: int,
    seed: int,
) -> pd.DataFrame:
    rows: List[dict] = []
    if pair_frame.empty:
        return pd.DataFrame()
    for (subset, target), block in pair_frame.groupby(["subset", "target"], sort=True):
        averaged = (
            block.groupby(
                [
                    "graph_group",
                    "canonical_graph",
                    "pair_id",
                    "left_id",
                    "right_id",
                ],
                as_index=False,
            )[["abs_delta_error", "concordance", "prediction_tie"]]
            .mean()
        )
        metric_map = {
            "pair_delta_mae": "abs_delta_error",
            "pair_concordance_ties_half": "concordance",
            "pair_tie_fraction": "prediction_tie",
        }
        for metric_name, column in metric_map.items():
            per_group = averaged.groupby("graph_group")[column].mean().to_numpy(dtype=float)
            if not len(per_group):
                continue
            rng = np.random.default_rng(seed + len(rows))
            draws = np.empty(replicates, dtype=float)
            for i in range(replicates):
                draws[i] = float(rng.choice(per_group, size=len(per_group), replace=True).mean())
            rows.append(
                {
                    "model": MODEL_LABEL,
                    "subset": subset,
                    "target": target,
                    "metric": metric_name,
                    "value": float(per_group.mean()),
                    "ci95_low": float(np.quantile(draws, 0.025)),
                    "ci95_high": float(np.quantile(draws, 0.975)),
                    "n_graph_clusters": int(len(per_group)),
                    "n_pairs": int(len(averaged)),
                    "n_seeds": int(block["seed"].nunique()),
                    "n_folds": int(block["fold"].nunique()),
                    "bootstrap_replicates": replicates,
                    "interval_scope": (
                        "graph-group cluster bootstrap over the selected OOF folds, conditional "
                        "on the fixed trained seeds"
                    ),
                }
            )
    return pd.DataFrame(rows)


def aggregate_ordinary_metrics(metrics: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for target, group in metrics.groupby("target", sort=True):
        row = {
            "model": MODEL_LABEL,
            "target": target,
            "n_fold_seed_runs": int(len(group)),
            "n_folds": int(group["fold"].nunique()),
            "n_seeds": int(group["seed"].nunique()),
            "parameters": int(group["parameters"].iloc[0]),
        }
        for metric in ("mae", "rmse", "r2", "spearman"):
            values = group[metric].to_numpy(dtype=float)
            row[f"{metric}_mean"] = float(values.mean())
            row[f"{metric}_sd"] = float(values.std(ddof=1)) if len(values) > 1 else 0.0
        rows.append(row)
    return pd.DataFrame(rows)


def write_schema_files(out: Path) -> None:
    prediction_schema = {
        "schema_version": SCHEMA_VERSION,
        "primary_file": "oof_predictions_long.csv",
        "columns": PREDICTION_COLUMNS,
        "merge_keys_with_2d_runner": [
            "schema_version",
            "fold",
            "seed",
            "mol_id",
            "graph_group",
            "target",
        ],
        "value_columns": ["y_true", "y_pred"],
        "model_label": MODEL_LABEL,
        "guardrail": MODEL_GUARDRAIL,
    }
    pair_schema = {
        "schema_version": SCHEMA_VERSION,
        "pair_prediction_file": "pair_predictions_long.csv",
        "pair_prediction_columns": PAIR_COLUMNS,
        "pair_metric_file": "pair_metrics.csv",
        "pair_comparison_merge_keys": ["subset", "target", "metric"],
        "frozen_pair_levels": list(PAIR_LEVELS),
        "pair_level_flag_columns": PAIR_LEVELS,
        "thresholds": {
            target: details["pair_threshold"] for target, details in TARGETS.items()
        },
        "material_pair_rule": "absolute frozen true pair delta >= the target threshold",
        "orientation": "true_delta and pred_delta are left_mol_id minus right_mol_id",
        "concordance_policy": (
            "sign match scores 1, sign mismatch scores 0, and abs(pred_delta) <= "
            "prediction_tie_tolerance scores 0.5"
        ),
        "fallback_note": (
            "Without a frozen pairs manifest only G0 can be reconstructed from canonical "
            "graphs; S1/S2/CONTEXT_EXACT/S3 require pairs.csv."
        ),
    }
    (out / "prediction_schema.json").write_text(
        json.dumps(prediction_schema, indent=2), encoding="utf-8"
    )
    (out / "pair_metric_schema.json").write_text(
        json.dumps(pair_schema, indent=2), encoding="utf-8"
    )


def main() -> None:
    args = parse_args()
    args.project_root = args.project_root.resolve()
    args.cohort_csv = args.cohort_csv.resolve()
    if args.manifest is not None:
        args.manifest = args.manifest.resolve()
    args.pairs_manifest = resolve_pairs_manifest_path(args)
    if args.manifest is None or args.pairs_manifest is None:
        raise ValueError("Production requires explicit records.csv and pairs.csv manifests")
    observed_manifest_hash = sha256(args.manifest)
    observed_pairs_hash = sha256(args.pairs_manifest)
    if not args.allow_nonfrozen_inputs:
        if observed_manifest_hash != FROZEN_RECORDS_MANIFEST_SHA256:
            raise RuntimeError(
                "Frozen records manifest hash mismatch: "
                f"{observed_manifest_hash} != {FROZEN_RECORDS_MANIFEST_SHA256}"
            )
        if observed_pairs_hash != FROZEN_PAIRS_MANIFEST_SHA256:
            raise RuntimeError(
                "Frozen pairs manifest hash mismatch: "
                f"{observed_pairs_hash} != {FROZEN_PAIRS_MANIFEST_SHA256}"
            )
    args.out = args.out.resolve()
    if args.out.exists() and any(args.out.iterdir()):
        raise FileExistsError(
            f"Refusing to mix or overwrite a nonempty output directory: {args.out}"
        )
    args.out.mkdir(parents=True, exist_ok=True)
    selected_folds = list(range(args.n_folds)) if args.folds is None else list(args.folds)
    if args.smoke:
        selected_folds = selected_folds[:1]
        args.seeds = args.seeds[:1]
        args.epochs = 1
        args.patience = 1
    invalid_selected = sorted(set(selected_folds) - set(range(args.n_folds)))
    if invalid_selected:
        raise ValueError(f"Requested folds outside 0..{args.n_folds - 1}: {invalid_selected}")
    if not selected_folds or not args.seeds:
        raise ValueError("At least one fold and seed are required")

    started = time.time()
    device = torch.device(args.device)
    print(f"Preparing CoordStatePairs OOF cohort on {device} ...", flush=True)
    manifest, property_path = prepare_fold_manifest(args)
    record_path = args.project_root / "pipeline_full_output" / "results.jsonl"
    records = load_records(record_path)
    missing_records = sorted(set(manifest["mol_id"]) - set(records))
    if missing_records:
        raise ValueError(f"Current CoordRep record file misses {len(missing_records)} manifest IDs")
    manifest, state_audit = annotate_state_subsets(
        manifest, records, args.cshm_difference_tolerance
    )
    pair_manifest, pair_manifest_audit = load_pair_manifest(args.pairs_manifest, manifest)
    coordinate_paths = xyz_paths(args.project_root)
    structures, coordinate_audit = parse_selected_xyz(coordinate_paths, manifest["mol_id"].tolist())
    dataset = make_base_dataset(manifest, structures)
    radius_graph_audit = complete_radius_graph_audit(
        dataset, args.max_neighbors, args.cutoff
    )
    preflight = invariance_preflight(args.out, device, args, radius_graph_audit)
    targets = manifest[list(TARGETS)].to_numpy(dtype=np.float32)

    fold_assignment_path = args.out / "fold_assignments.csv"
    manifest[
        [
            "mol_id",
            "graph_group",
            "canonical_graph",
            "fold",
            "fold_source",
            "is_graph_degenerate",
            "continuous_only",
        ]
    ].to_csv(fold_assignment_path, index=False)
    write_schema_files(args.out)

    prediction_path = args.out / "oof_predictions_long.csv"
    if prediction_path.exists():
        prediction_path.unlink()
    run_metadata: Dict[Tuple[int, int], dict] = {}
    for fold in selected_folds:
        train_idx, val_idx, test_idx = deterministic_validation_indices(
            manifest, fold, args.validation_fraction, args.fold_salt
        )
        # Defensive graph-group disjointness audit.
        train_groups = set(manifest.iloc[train_idx]["graph_group"])
        val_groups = set(manifest.iloc[val_idx]["graph_group"])
        test_groups = set(manifest.iloc[test_idx]["graph_group"])
        if train_groups & val_groups or train_groups & test_groups or val_groups & test_groups:
            raise RuntimeError(f"Graph-group leakage detected for fold {fold}")
        for seed in args.seeds:
            print(
                f"Training {MODEL_LABEL} OOF fold={fold}, seed={seed}; "
                f"train/val/test={len(train_idx):,}/{len(val_idx):,}/{len(test_idx):,}",
                flush=True,
            )
            predictions, truths, row_ids, metadata = train_fold_seed(
                dataset,
                targets,
                train_idx,
                val_idx,
                test_idx,
                fold,
                seed,
                args,
                device,
            )
            rows = prediction_rows(
                manifest, fold, seed, row_ids, predictions, truths
            )
            append_csv(prediction_path, rows, PREDICTION_COLUMNS)
            run_metadata[(fold, seed)] = metadata
            print(
                f"Finished fold={fold}, seed={seed}, best_epoch={metadata['best_epoch']}, "
                f"elapsed={metadata['elapsed_seconds']:.1f}s, prediction_rows={len(rows):,}",
                flush=True,
            )

    predictions = pd.read_csv(prediction_path)
    predictions["mol_id"] = predictions["mol_id"].astype(str)
    merge_key_columns = ["model", "fold", "seed", "mol_id", "graph_group", "target"]
    duplicate_prediction_keys = int(predictions.duplicated(merge_key_columns).sum())
    if duplicate_prediction_keys:
        raise RuntimeError(f"Duplicate OOF prediction merge keys: {duplicate_prediction_keys}")
    if not np.isfinite(predictions[["y_true", "y_pred"]].to_numpy(dtype=float)).all():
        raise RuntimeError("Non-finite OOF prediction output")

    ordinary_metrics = ordinary_metric_rows(predictions, run_metadata)
    ordinary_summary = aggregate_ordinary_metrics(ordinary_metrics)
    pair_predictions, pair_by_fold_seed, pair_by_seed = pair_outputs(
        predictions, pair_manifest, args.prediction_tie_tolerance
    )
    pair_summary = bootstrap_pair_metrics(
        pair_predictions, args.bootstrap_replicates, args.bootstrap_seed
    )
    ordinary_metrics.to_csv(args.out / "oof_metrics_by_fold_seed.csv", index=False)
    ordinary_summary.to_csv(args.out / "oof_metrics_summary.csv", index=False)
    pair_predictions.to_csv(args.out / "pair_predictions_long.csv", index=False)
    pair_by_fold_seed.to_csv(args.out / "pair_metrics_by_fold_seed.csv", index=False)
    pair_by_seed.to_csv(args.out / "pair_metrics_by_seed.csv", index=False)
    pair_summary.to_csv(args.out / "pair_metrics.csv", index=False)

    selected_ids = set(manifest.loc[manifest["fold"].isin(selected_folds), "mol_id"])
    expected_prediction_rows = len(selected_ids) * len(args.seeds) * len(TARGETS)
    observed_test_ids = set(predictions["mol_id"])
    audit = {
        "schema_version": SCHEMA_VERSION,
        "run_mode": "smoke" if args.smoke else "requested_oof_run",
        "model": MODEL_LABEL,
        "cohort_records": int(len(manifest)),
        "cohort_graph_groups": int(manifest["graph_group"].nunique()),
        "n_folds_defined": args.n_folds,
        "selected_folds": selected_folds,
        "selected_seeds": args.seeds,
        "fold_record_counts": {
            str(key): int(value) for key, value in manifest["fold"].value_counts().sort_index().items()
        },
        "graph_groups_split_across_folds": int(
            (manifest.groupby("graph_group")["fold"].nunique() != 1).sum()
        ),
        "prediction_rows": int(len(predictions)),
        "expected_prediction_rows": int(expected_prediction_rows),
        "unique_predicted_molecules": int(len(observed_test_ids)),
        "expected_predicted_molecules": int(len(selected_ids)),
        "missing_selected_test_ids": int(len(selected_ids - observed_test_ids)),
        "unexpected_test_ids": int(len(observed_test_ids - selected_ids)),
        "duplicate_prediction_merge_keys": duplicate_prediction_keys,
        "finite_prediction_values": bool(
            np.isfinite(predictions[["y_true", "y_pred"]].to_numpy(dtype=float)).all()
        ),
        "pair_prediction_rows": int(len(pair_predictions)),
        "pair_metric_rows": int(len(pair_summary)),
        "pair_manifest_audit": pair_manifest_audit,
        "state_subset_audit": state_audit,
        "coordinate_audit": coordinate_audit,
        "radius_graph_audit": radius_graph_audit,
        "invariance_preflight": preflight,
        "fold_assignment_source": str(manifest["fold_source"].iloc[0]),
        "validation_policy": (
            f"deterministic graph-group hash within each outer-training partition; "
            f"fraction={args.validation_fraction} of non-test graph groups"
        ),
        "target_leakage_controls": [
            "outer test membership determined only by graph group/manifest",
            "identical graph groups cannot cross train/validation/test",
            "target scaling fitted on the training subset of each fold only",
            "checkpoint selected on validation normalized MAE only",
            "test targets used only after checkpoint selection",
            "XYZ comment metadata were excluded from model inputs",
        ],
        "cartesian_input_guardrail": MODEL_GUARDRAIL,
        "requested_training_precision": args.precision,
        "effective_training_precision": (
            args.precision if device.type == "cuda" else "float32"
        ),
        "raw_runner_summary_scope": (
            "Diagnostic production-run summaries. Manuscript estimands are computed by "
            "the separately frozen postprocessor after recordwise seed averaging."
        ),
        "pair_interval_scope": (
            "graph-group cluster bootstrap over selected OOF folds, conditional on fixed trained seeds"
        ),
        "input_sha256": {
            "properties": sha256(property_path),
            "coordrep_records": sha256(record_path),
            "cohort_csv": sha256(args.cohort_csv) if args.manifest is None else None,
            "shared_manifest": sha256(args.manifest) if args.manifest is not None else None,
            "shared_pairs_manifest": (
                sha256(args.pairs_manifest) if args.pairs_manifest is not None else None
            ),
            **{path.name: sha256(path) for path in coordinate_paths},
        },
        "total_elapsed_seconds": time.time() - started,
    }
    with (args.out / "dataset_audit.json").open("w", encoding="utf-8") as handle:
        json.dump(audit, handle, indent=2, ensure_ascii=False, allow_nan=False)
    with (args.out / "run_metadata.json").open("w", encoding="utf-8") as handle:
        json.dump(
            {
                "arguments": {
                    **vars(args),
                    "project_root": str(args.project_root),
                    "cohort_csv": str(args.cohort_csv),
                    "manifest": str(args.manifest) if args.manifest else None,
                    "pairs_manifest": (
                        str(args.pairs_manifest) if args.pairs_manifest else None
                    ),
                    "out": str(args.out),
                },
                "device": str(device),
                "torch_version": torch.__version__,
                "runs": {f"fold{fold}_seed{seed}": value for (fold, seed), value in run_metadata.items()},
            },
            handle,
            indent=2,
            allow_nan=False,
        )
    write_model_card(args.out, args, preflight)
    print(json.dumps(audit, indent=2, ensure_ascii=False), flush=True)
    print("\nOOF metric smoke/full summary:", flush=True)
    print(ordinary_summary.to_string(index=False), flush=True)
    print("\nPair metric summary:", flush=True)
    print(pair_summary.to_string(index=False) if not pair_summary.empty else "No threshold pairs", flush=True)


if __name__ == "__main__":
    main()
