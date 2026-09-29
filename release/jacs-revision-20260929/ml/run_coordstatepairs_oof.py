#!/usr/bin/env python
"""Five-fold OOF benchmark for coordination-state pairs.

The benchmark targets a precise information gap: two tmQM records can have
the same canonical attributed 2D graph but different coordination states.
An invariant 2D GINE must return the same prediction for those graph members,
whereas CoordRep fields can distinguish them.  This runner compares the
capacity-matched GINE-wide baseline with the current GINE+CoordRep hybrid.

Outer folds are read from the released SHA256 graph-fold manifest
``coordstate_pairs_manifest_v1/records.csv``.  The released ``pairs.csv`` is
also used by default so that G0/S1/S2/CONTEXT_EXACT/S3 membership is frozen.
In every outer fold, the validation set is drawn by graph hash from the outer
training pool, and DictVectorizer, StandardScaler, and target normalization
are fitted on the final training partition only.  A legacy SGKF generator is
retained solely for explicitly requested development runs.

The script writes auditable member-level OOF predictions and long-format
same-graph pair predictions/metrics.  For material property differences,
predicted graph ties receive 0.5 in the primary concordance metric; non-tie
direction accuracy is reported separately.
"""

from __future__ import annotations

import argparse
import csv
import fnmatch
import hashlib
import json
import time
from collections import Counter
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np
import pandas as pd
import torch
from rdkit import Chem, RDLogger
from sklearn.feature_extraction import DictVectorizer
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import GroupShuffleSplit, StratifiedGroupKFold
from sklearn.preprocessing import StandardScaler

from run_downstream_graph_benchmark import (
    PreparedData,
    TARGET_SPECS,
    file_sha256,
    graph_loader,
    load_current_records,
    make_strata,
    mol_to_graph,
    predict_graph,
    record_to_field_dict,
    train_graph_model,
)


RDLogger.DisableLog("rdApp.*")


FORMAL_MANIFEST_DIR = Path(
    "revision_experiments/results/coordstate_pairs_manifest_v1"
)
FORMAL_RECORDS_MANIFEST = FORMAL_MANIFEST_DIR / "records.csv"
FORMAL_PAIRS_MANIFEST = FORMAL_MANIFEST_DIR / "pairs.csv"
FORMAL_FOLD_SIZES = {0: 9573, 1: 9706, 2: 9631, 3: 9576, 4: 9571}
FORMAL_PAIR_FOLD_SIZES = {0: 511, 1: 572, 2: 721, 3: 668, 4: 505}
TREX_MATCHED_COHORT = Path(
    "revision_experiments/results/"
    "trex_common_arch_coordstatepairs_oof_smoke_v3/matched_cohort.csv"
)
TREX_MATCHED_FOLD_SIZES = {0: 6236, 1: 6348, 2: 6233, 3: 6139, 4: 6163}

# This named variant is deliberately built in so its scientific meaning cannot
# drift between formal runs.  A selector without a glob also masks its
# ``*_missing`` companion, preventing missingness from leaking the held-out
# field back to the model.
BUILTIN_HYBRID_FIELD_MASKS = {
    "hybrid_mask_cn_oxidation": ("cn", "oxidation", "dcount"),
}


@dataclass
class RawCohort:
    ids: List[str]
    smiles: List[str]
    graphs: list
    field_dicts: List[Dict[str, float]]
    targets: np.ndarray
    target_names: List[str]
    metals: List[str]
    coordination_numbers: List[int]
    shapes: List[str]
    coordrep_hashes: List[str]
    audit: Dict[str, object]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path("tmp/model_audit/data/CoordRep/CoordSMILES"),
    )
    parser.add_argument(
        "--record-jsonl",
        type=Path,
        default=None,
        help=(
            "Optional non-overwriting CoordRep record JSONL. When omitted, "
            "<project-root>/pipeline_full_output/results.jsonl is used."
        ),
    )
    parser.add_argument(
        "--property-csv",
        type=Path,
        default=None,
        help=(
            "Optional property table. When omitted, "
            "<project-root>/data/processed_full/dataset_coordsmiles.csv is used."
        ),
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("revision_experiments/results/coordstate_pairs_oof_v1"),
    )
    parser.add_argument(
        "--fold-manifest",
        type=Path,
        default=FORMAL_RECORDS_MANIFEST,
        help="Frozen record/fold manifest (formal default: coordstatepairs-v1).",
    )
    parser.add_argument(
        "--pair-manifest",
        type=Path,
        default=FORMAL_PAIRS_MANIFEST,
        help="Frozen G0/S1/S2/CONTEXT_EXACT/S3 pair manifest.",
    )
    cohort_group = parser.add_mutually_exclusive_group()
    cohort_group.add_argument(
        "--include-ids",
        nargs="+",
        default=None,
        help="Restrict evaluation to these mol_id values while retaining frozen folds.",
    )
    cohort_group.add_argument(
        "--cohort-manifest",
        type=Path,
        default=None,
        help="CSV containing a mol_id column used to restrict the frozen cohort.",
    )
    parser.add_argument(
        "--development-derive-pairs-from-members",
        action="store_true",
        help="DEVELOPMENT ONLY: ignore pairs.csv and derive G0 pairs at run time.",
    )
    parser.add_argument(
        "--development-create-sgkf-manifest-if-missing",
        action="store_true",
        help=(
            "DEVELOPMENT ONLY: create a legacy SGKF manifest at a non-formal path "
            "when --fold-manifest is missing."
        ),
    )
    parser.add_argument("--n-folds", type=int, default=5)
    parser.add_argument("--folds", nargs="+", type=int, default=None)
    parser.add_argument("--targets", nargs="+", default=["hl_gap_ev", "dipole_moment"])
    parser.add_argument("--split-seed", type=int, default=20260723)
    parser.add_argument("--seeds", nargs="+", type=int, default=[11, 22, 33])
    parser.add_argument("--models", nargs="+", default=["gine_wide", "hybrid"])
    parser.add_argument(
        "--field-mode",
        action="append",
        default=[],
        metavar="MODEL=SELECTOR[,SELECTOR...]",
        help=(
            "Register an additional masked-hybrid model. Exact selectors also "
            "mask SELECTOR_missing; shell-style wildcards are supported. The "
            "frozen built-in mode hybrid_mask_cn_oxidation masks cn, oxidation, "
            "dcount, and their missingness indicators."
        ),
    )
    parser.add_argument("--validation-fraction", type=float, default=0.10)
    parser.add_argument("--epochs", type=int, default=40)
    parser.add_argument("--patience", type=int, default=8)
    parser.add_argument("--batch-size", type=int, default=256)
    parser.add_argument("--hidden", type=int, default=128)
    parser.add_argument("--layers", type=int, default=4)
    parser.add_argument("--lr", type=float, default=1.0e-3)
    parser.add_argument("--weight-decay", type=float, default=1.0e-5)
    parser.add_argument("--num-workers", type=int, default=0)
    parser.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    # Same canonical graphs can differ at ~1e-6 solely from GPU reduction
    # order after RDKit atom reordering.  A 1e-5 tolerance safely absorbs that
    # numerical noise while remaining negligible for both reported targets.
    parser.add_argument("--prediction-tie-tolerance", type=float, default=1.0e-5)
    parser.add_argument("--gap-min-delta", type=float, default=0.1)
    parser.add_argument("--dipole-min-delta", type=float, default=1.0)
    parser.add_argument(
        "--smoke",
        action="store_true",
        help="Run fold 0, first seed, and one epoch on the full cohort.",
    )
    return parser.parse_args()


def parse_field_modes(specs: Sequence[str]) -> Dict[str, Tuple[str, ...]]:
    """Return immutable built-in plus user-named hybrid field masks."""
    modes = dict(BUILTIN_HYBRID_FIELD_MASKS)
    for raw_spec in specs:
        if "=" not in raw_spec:
            raise ValueError(
                f"Invalid --field-mode {raw_spec!r}; expected MODEL=SELECTOR[,SELECTOR...]"
            )
        model_label, raw_selectors = raw_spec.split("=", 1)
        model_label = model_label.strip()
        selectors = tuple(
            dict.fromkeys(
                selector.strip()
                for selector in raw_selectors.split(",")
                if selector.strip()
            )
        )
        if not model_label.startswith("hybrid_mask_"):
            raise ValueError(
                f"Field-mode model label must start with 'hybrid_mask_': {model_label!r}"
            )
        if not selectors:
            raise ValueError(f"Field mode {model_label!r} has no selectors")
        if model_label in modes:
            raise ValueError(
                f"Field mode {model_label!r} is already defined and may not be overridden"
            )
        modes[model_label] = selectors
    return modes


def resolve_field_mask(
    field_names: Sequence[str], selectors: Sequence[str]
) -> Tuple[List[int], Dict[str, List[str]]]:
    """Resolve auditable selectors against the train-fitted field vocabulary."""
    names = [str(name) for name in field_names]
    matched_by_selector: Dict[str, List[str]] = {}
    selected_indices = set()
    for selector in selectors:
        has_glob = any(token in selector for token in ("*", "?", "["))
        if has_glob:
            matched = [
                (index, name)
                for index, name in enumerate(names)
                if fnmatch.fnmatchcase(name, selector)
            ]
        else:
            allowed = {selector, f"{selector}_missing"}
            matched = [
                (index, name)
                for index, name in enumerate(names)
                if name in allowed
            ]
        if not matched:
            raise ValueError(
                f"Field selector {selector!r} matched no train-fitted columns"
            )
        matched_by_selector[selector] = [name for _, name in matched]
        selected_indices.update(index for index, _ in matched)
    return sorted(selected_indices), matched_by_selector


def apply_model_field_view(
    data: PreparedData,
    model_label: str,
    field_modes: Dict[str, Tuple[str, ...]],
) -> Dict[str, object]:
    """Attach either full or zero-masked fields to every graph.

    Zero is the training mean after train-only StandardScaler fitting.  Keeping
    the full input width makes the current and masked hybrid exactly
    parameter-matched; only field information changes.
    """
    selectors = field_modes.get(model_label, ())
    masked_indices, matched_by_selector = resolve_field_mask(
        data.field_names, selectors
    )
    model_fields = data.fields if not masked_indices else data.fields.copy()
    if masked_indices:
        model_fields[:, masked_indices] = 0.0

    for index, graph in enumerate(data.graphs):
        graph.coord_x = torch.from_numpy(model_fields[index : index + 1])

    partition_max_abs: Dict[str, float] = {}
    if masked_indices:
        for partition, indices in (
            ("train", data.train_idx),
            ("validation", data.val_idx),
            ("test", data.test_idx),
        ):
            values = model_fields[np.asarray(indices)][:, masked_indices]
            partition_max_abs[partition] = (
                float(np.max(np.abs(values))) if values.size else 0.0
            )
        if any(value != 0.0 for value in partition_max_abs.values()):
            raise RuntimeError(
                f"Nonzero masked field survived in {model_label}: {partition_max_abs}"
            )

    return {
        "model_label": model_label,
        "internal_architecture": (
            "gine_wide" if model_label == "gine_wide" else "hybrid"
        ),
        "uses_coordrep_fields": model_label != "gine_wide",
        "field_view": "zero_masked" if masked_indices else "full",
        "requested_selectors": list(selectors),
        "matched_columns_by_selector": matched_by_selector,
        "masked_columns": [str(data.field_names[index]) for index in masked_indices],
        "n_input_columns": int(data.fields.shape[1]),
        "n_masked_columns": len(masked_indices),
        "input_dimension_preserved": True,
        "mask_value_after_train_only_standardization": 0.0,
        "mask_value_interpretation": (
            "training-partition mean; no selected-field variation remains"
            if masked_indices
            else "not_applicable"
        ),
        "masked_partition_max_abs": partition_max_abs,
    }


def audit_hybrid_parameter_control(
    run_details: Dict[str, object],
    selected_folds: Sequence[int],
    seeds: Sequence[int],
    model_labels: Sequence[str],
    field_modes: Dict[str, Tuple[str, ...]],
) -> Dict[str, object]:
    """Assert exact parameter matching for full versus masked hybrid runs."""
    masked_labels = [label for label in model_labels if label in field_modes]
    comparisons: List[Dict[str, object]] = []
    if "hybrid" in model_labels:
        for fold in selected_folds:
            for seed in seeds:
                reference_key = f"fold{fold}_hybrid_seed{seed}"
                reference_parameters = int(run_details[reference_key]["parameters"])
                for model_label in masked_labels:
                    masked_key = f"fold{fold}_{model_label}_seed{seed}"
                    masked_parameters = int(run_details[masked_key]["parameters"])
                    equal = reference_parameters == masked_parameters
                    comparisons.append(
                        {
                            "fold": int(fold),
                            "seed": int(seed),
                            "reference_model": "hybrid",
                            "masked_model": model_label,
                            "reference_parameters": reference_parameters,
                            "masked_parameters": masked_parameters,
                            "exactly_equal": equal,
                        }
                    )
                    if not equal:
                        raise RuntimeError(
                            f"Hybrid parameter-control failure: {reference_key}="
                            f"{reference_parameters}, {masked_key}={masked_parameters}"
                        )
    return {
        "policy": (
            "Keep the train-fitted field dimension and hybrid architecture fixed; "
            "replace selected standardized columns by zero in train, validation, "
            "and test. Same fold/seed calls reset the RNG before construction."
        ),
        "reference_model_requested": "hybrid" in model_labels,
        "masked_models_requested": masked_labels,
        "n_paired_comparisons": len(comparisons),
        "all_paired_parameter_counts_exactly_equal": (
            all(row["exactly_equal"] for row in comparisons)
            if comparisons
            else None
        ),
        "comparisons": comparisons,
    }


def sha256_text(text: str, length: int = 20) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()[:length]


def full_sha256_text(text: str) -> str:
    return sha256_text(text, length=64)


def best_shape(record: dict) -> str:
    refs = record.get("shape_ref") or []
    values = record.get("shape_values") or []
    if not refs or not values:
        return "NONE"
    return str(min(zip(refs, values), key=lambda item: float(item[1]))[0])


def write_csv(path: Path, rows: Sequence[Dict[str, object]]) -> None:
    if not rows:
        return
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def load_raw_cohort(args: argparse.Namespace) -> RawCohort:
    root = args.project_root.resolve()
    property_path = (
        args.property_csv.resolve()
        if args.property_csv is not None
        else root / "data" / "processed_full" / "dataset_coordsmiles.csv"
    )
    record_path = (
        args.record_jsonl.resolve()
        if args.record_jsonl is not None
        else root / "pipeline_full_output" / "results.jsonl"
    )
    if not property_path.exists() or not record_path.exists():
        raise FileNotFoundError(f"Missing {property_path} or {record_path}")

    unknown_targets = sorted(set(args.targets) - set(TARGET_SPECS))
    if unknown_targets:
        raise ValueError(f"Unknown targets: {unknown_targets}")

    records, pipeline_audit = load_current_records(record_path)
    needed_columns = ["csd_code", "original_smiles"] + [
        TARGET_SPECS[target][0] for target in args.targets
    ]
    frame = pd.read_csv(property_path, usecols=needed_columns)
    raw_rows = len(frame)
    duplicate_id_rows = int(frame.duplicated("csd_code").sum())
    frame = frame.drop_duplicates("csd_code", keep="first").copy()
    frame["csd_code"] = frame["csd_code"].astype(str)
    frame = frame[frame["csd_code"].isin(records) & frame["original_smiles"].notna()].copy()
    for target in args.targets:
        source, _ = TARGET_SPECS[target]
        frame = frame[np.isfinite(pd.to_numeric(frame[source], errors="coerce"))]

    ids: List[str] = []
    smiles: List[str] = []
    graphs = []
    field_dicts: List[Dict[str, float]] = []
    target_rows: List[List[float]] = []
    metals: List[str] = []
    coordination_numbers: List[int] = []
    shapes: List[str] = []
    coordrep_hashes: List[str] = []
    parse_failures: List[str] = []

    for row in frame.itertuples(index=False):
        mol_id = str(row.csd_code)
        molecule = Chem.MolFromSmiles(str(row.original_smiles))
        if molecule is None or molecule.GetNumAtoms() == 0:
            parse_failures.append(mol_id)
            continue
        record = records[mol_id]
        canonical = Chem.MolToSmiles(
            molecule, canonical=True, isomericSmiles=True
        )
        ids.append(mol_id)
        smiles.append(canonical)
        graphs.append(mol_to_graph(molecule))
        field_dicts.append(record_to_field_dict(record))
        target_rows.append(
            [
                float(getattr(row, TARGET_SPECS[target][0])) * TARGET_SPECS[target][1]
                for target in args.targets
            ]
        )
        metals.append(str(record.get("metal") or "UNK"))
        coordination_numbers.append(int(record.get("cn") or -1))
        shapes.append(best_shape(record))
        coordrep_hashes.append(full_sha256_text(str(record.get("coordrep") or "")))

    targets = np.asarray(target_rows, dtype=np.float32)
    graph_counts = Counter(smiles)
    same_graph_groups = sum(count >= 2 for count in graph_counts.values())
    same_graph_records = sum(count for count in graph_counts.values() if count >= 2)
    audit: Dict[str, object] = {
        **pipeline_audit,
        "property_rows_raw": raw_rows,
        "property_duplicate_id_rows": duplicate_id_rows,
        "joined_rows_before_rdkit": len(frame),
        "rdkit_parse_failures": len(parse_failures),
        "final_records": len(ids),
        "unique_canonical_graphs": len(graph_counts),
        "same_graph_groups": same_graph_groups,
        "same_graph_records": same_graph_records,
        "target_names": list(args.targets),
        "metal_counts": dict(Counter(metals)),
        "cn_counts": {str(k): v for k, v in sorted(Counter(coordination_numbers).items())},
        "shape_counts": dict(Counter(shapes)),
        "source_sha256": {
            "properties": file_sha256(property_path),
            "coordrep_records": file_sha256(record_path),
        },
    }
    return RawCohort(
        ids=ids,
        smiles=smiles,
        graphs=graphs,
        field_dicts=field_dicts,
        targets=targets,
        target_names=list(args.targets),
        metals=metals,
        coordination_numbers=coordination_numbers,
        shapes=shapes,
        coordrep_hashes=coordrep_hashes,
        audit=audit,
    )


def create_development_fold_manifest(
    cohort: RawCohort, args: argparse.Namespace
) -> pd.DataFrame:
    if args.n_folds < 2:
        raise ValueError("--n-folds must be at least 2")
    strata = make_strata(cohort.metals, cohort.targets)
    splitter = StratifiedGroupKFold(
        n_splits=args.n_folds, shuffle=True, random_state=args.split_seed
    )
    fold_by_row = np.full(len(cohort.ids), -1, dtype=int)
    dummy = np.zeros(len(cohort.ids), dtype=np.int8)
    for fold, (_, test_idx) in enumerate(
        splitter.split(dummy, y=strata, groups=np.asarray(cohort.smiles, dtype=object))
    ):
        fold_by_row[test_idx] = fold
    if np.any(fold_by_row < 0):
        raise RuntimeError("Manifest generation did not assign every record")

    manifest = pd.DataFrame({"mol_id": cohort.ids, "fold": fold_by_row})
    graph_fold_count = (
        pd.DataFrame({"graph": cohort.smiles, "fold": fold_by_row})
        .groupby("graph")["fold"]
        .nunique()
    )
    if int((graph_fold_count > 1).sum()) != 0:
        raise RuntimeError("Canonical graph leakage detected while generating manifest")
    args.fold_manifest.parent.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(args.fold_manifest, index=False)
    metadata = {
        "status": "DEVELOPMENT_ONLY_NOT_FOR_FORMAL_ANALYSIS",
        "n_records": len(manifest),
        "n_folds": args.n_folds,
        "split_seed": args.split_seed,
        "strategy": (
            "DEVELOPMENT ONLY: StratifiedGroupKFold by metal and gap quintile; "
            "canonical RDKit SMILES groups"
        ),
        "fold_counts": {
            str(int(k)): int(v) for k, v in manifest["fold"].value_counts().sort_index().items()
        },
        "cohort_source_sha256": cohort.audit["source_sha256"],
        "manifest_sha256": file_sha256(args.fold_manifest),
    }
    with (args.fold_manifest.parent / "manifest_metadata.json").open(
        "w", encoding="utf-8"
    ) as handle:
        json.dump(metadata, handle, indent=2)
    return manifest


def load_fold_manifest(cohort: RawCohort, args: argparse.Namespace) -> pd.DataFrame:
    formal_path = args.fold_manifest.resolve() == FORMAL_RECORDS_MANIFEST.resolve()
    if args.fold_manifest.exists():
        manifest = pd.read_csv(args.fold_manifest, dtype={"mol_id": str, "fold": int})
    elif args.development_create_sgkf_manifest_if_missing:
        if formal_path:
            raise ValueError(
                "The formal records.csv may not be replaced by an SGKF manifest. "
                "Choose a separate --fold-manifest path for a development run."
            )
        manifest = create_development_fold_manifest(cohort, args)
    else:
        raise FileNotFoundError(
            f"Frozen manifest not found: {args.fold_manifest}. "
            "Formal runs require the released records.csv."
        )
    required = {"mol_id", "fold"}
    if not required.issubset(manifest.columns):
        raise ValueError(f"Manifest needs columns {sorted(required)}")
    if manifest["mol_id"].duplicated().any():
        raise ValueError("Manifest has duplicate molecule IDs")
    expected, actual = set(cohort.ids), set(manifest["mol_id"])
    if expected != actual:
        raise ValueError(
            f"Manifest/cohort mismatch: missing={len(expected-actual)}, extra={len(actual-expected)}"
        )
    invalid = sorted(set(manifest["fold"].astype(int)) - set(range(args.n_folds)))
    if invalid:
        raise ValueError(f"Manifest has invalid folds: {invalid}")
    order = {mol_id: i for i, mol_id in enumerate(cohort.ids)}
    manifest = manifest.assign(_order=manifest["mol_id"].map(order)).sort_values("_order")
    manifest = manifest.drop(columns="_order").reset_index(drop=True)

    computed_graph_hashes = np.asarray(
        [full_sha256_text(smiles) for smiles in cohort.smiles], dtype=object
    )
    if "graph_sha256" in manifest:
        frozen_hashes = manifest["graph_sha256"].astype(str).to_numpy()
        mismatches = int(np.count_nonzero(frozen_hashes != computed_graph_hashes))
        if mismatches:
            raise ValueError(
                f"Frozen graph SHA256 disagrees with the loaded cohort for {mismatches} rows"
            )
    else:
        manifest["graph_sha256"] = computed_graph_hashes

    fallback_columns = {
        "metal": cohort.metals,
        "cn": cohort.coordination_numbers,
        "best_shape": cohort.shapes,
        "context_sha256": [""] * len(cohort.ids),
        "coordrep_no_v_sha256": cohort.coordrep_hashes,
    }
    for column, values in fallback_columns.items():
        if column not in manifest:
            manifest[column] = values

    target_manifest_columns = {
        "hl_gap_ev": "hl_gap_ev",
        "dipole_moment": "dipole_moment_d",
    }
    target_checks: Dict[str, object] = {}
    for target_index, target in enumerate(cohort.target_names):
        manifest_column = target_manifest_columns.get(target)
        if manifest_column is None:
            continue
        loaded = cohort.targets[:, target_index].astype(float)
        if manifest_column not in manifest:
            manifest[manifest_column] = loaded
            continue
        frozen = pd.to_numeric(manifest[manifest_column], errors="raise").to_numpy(float)
        max_difference = float(np.max(np.abs(frozen - loaded)))
        if not np.allclose(frozen, loaded, rtol=0.0, atol=2.0e-5):
            raise ValueError(
                f"Frozen {manifest_column} disagrees with the loaded cohort; "
                f"max absolute difference={max_difference:.6g}"
            )
        target_checks[manifest_column] = {"max_absolute_difference": max_difference}

    graph_folds = (
        manifest[["graph_sha256", "fold"]]
        .rename(columns={"graph_sha256": "graph"})
        .groupby("graph")["fold"]
        .nunique()
    )
    if int((graph_folds > 1).sum()) != 0:
        raise ValueError("Frozen manifest leaks graph hashes across outer folds")

    fold_counts = {
        int(k): int(v)
        for k, v in manifest["fold"].value_counts().sort_index().items()
    }
    if formal_path and fold_counts != FORMAL_FOLD_SIZES:
        raise ValueError(
            f"Formal manifest fold sizes changed: {fold_counts}; "
            f"expected {FORMAL_FOLD_SIZES}"
        )
    cohort.audit["fold_manifest_audit"] = {
        "status": "FORMAL_FROZEN_HASH_FOLDS" if formal_path else "DEVELOPMENT_CUSTOM",
        "path": str(args.fold_manifest.resolve()),
        "sha256": file_sha256(args.fold_manifest),
        "fold_counts": {str(k): v for k, v in fold_counts.items()},
        "graph_sha256_mismatch_count": 0,
        "graph_assigned_to_multiple_folds": int((graph_folds > 1).sum()),
        "target_checks": target_checks,
    }
    return manifest


def restrict_analysis_cohort(
    cohort: RawCohort, manifest: pd.DataFrame, args: argparse.Namespace
) -> Tuple[RawCohort, pd.DataFrame]:
    if args.include_ids is None and args.cohort_manifest is None:
        fold_counts = {
            int(k): int(v)
            for k, v in manifest["fold"].value_counts().sort_index().items()
        }
        cohort.audit["analysis_cohort"] = {
            "status": "FULL_FORMAL_COHORT",
            "n_records": len(cohort.ids),
            "fold_counts": {str(k): v for k, v in fold_counts.items()},
        }
        return cohort, manifest

    cohort_manifest_fold_mismatches = 0
    if args.cohort_manifest is not None:
        if not args.cohort_manifest.exists():
            raise FileNotFoundError(f"Cohort manifest not found: {args.cohort_manifest}")
        header = pd.read_csv(args.cohort_manifest, nrows=0).columns
        if "mol_id" not in header:
            raise ValueError("--cohort-manifest must contain a mol_id column")
        usecols = ["mol_id"] + (["fold"] if "fold" in header else [])
        cohort_frame = pd.read_csv(
            args.cohort_manifest,
            usecols=usecols,
            dtype={"mol_id": str},
        )
        duplicate_ids = int(cohort_frame["mol_id"].duplicated().sum())
        if duplicate_ids:
            raise ValueError(
                f"Cohort manifest has {duplicate_ids} duplicate mol_id rows"
            )
        requested_ids = cohort_frame["mol_id"].astype(str).tolist()
        source = "COHORT_MANIFEST"
        source_path = str(args.cohort_manifest.resolve())
        source_sha256 = file_sha256(args.cohort_manifest)
    else:
        requested_ids = [str(value) for value in args.include_ids]
        duplicate_ids = len(requested_ids) - len(set(requested_ids))
        if duplicate_ids:
            raise ValueError(f"--include-ids contains {duplicate_ids} duplicate IDs")
        cohort_frame = None
        source = "INCLUDE_IDS_ARGUMENT"
        source_path = None
        source_sha256 = None

    full_id_set = set(cohort.ids)
    requested_set = set(requested_ids)
    missing_ids = sorted(requested_set - full_id_set)
    if missing_ids:
        preview = ", ".join(missing_ids[:10])
        raise ValueError(
            f"Requested cohort has {len(missing_ids)} IDs outside the formal cohort: {preview}"
        )
    if not requested_set:
        raise ValueError("Requested analysis cohort is empty")

    indices = np.asarray(
        [index for index, mol_id in enumerate(cohort.ids) if mol_id in requested_set],
        dtype=int,
    )
    subset_manifest = manifest.iloc[indices].reset_index(drop=True).copy()
    if cohort_frame is not None and "fold" in cohort_frame:
        frozen_fold_by_id = subset_manifest.set_index("mol_id")["fold"]
        supplied_folds = pd.to_numeric(cohort_frame["fold"], errors="raise").astype(int)
        expected_folds = cohort_frame["mol_id"].map(frozen_fold_by_id).astype(int)
        cohort_manifest_fold_mismatches = int((supplied_folds != expected_folds).sum())
        if cohort_manifest_fold_mismatches:
            raise ValueError(
                "Cohort manifest attempts to change "
                f"{cohort_manifest_fold_mismatches} frozen fold assignments"
            )

    fold_counts = {
        int(k): int(v)
        for k, v in subset_manifest["fold"].value_counts().sort_index().items()
    }
    is_trex_matched = (
        args.cohort_manifest is not None
        and args.cohort_manifest.resolve() == TREX_MATCHED_COHORT.resolve()
    )
    known_excluded = {"ZEMNAJ", "ZEMNEN"}
    known_excluded_present = sorted(known_excluded & requested_set)
    if is_trex_matched:
        if len(requested_set) != 31119:
            raise ValueError(
                f"Frozen T-REX matched cohort has {len(requested_set)} IDs; expected 31119"
            )
        if fold_counts != TREX_MATCHED_FOLD_SIZES:
            raise ValueError(
                f"T-REX matched fold sizes changed: {fold_counts}; "
                f"expected {TREX_MATCHED_FOLD_SIZES}"
            )
        if known_excluded_present:
            raise ValueError(
                "T-REX matched cohort must exclude ZEMNAJ and ZEMNEN; present="
                f"{known_excluded_present}"
            )

    def select_list(values: Sequence[object]) -> list:
        return [values[index] for index in indices]

    audit = dict(cohort.audit)
    audit["analysis_cohort"] = {
        "status": "FROZEN_TREX_MATCHED_COHORT" if is_trex_matched else source,
        "source_path": source_path,
        "source_sha256": source_sha256,
        "n_requested_ids": len(requested_ids),
        "n_unique_requested_ids": len(requested_set),
        "n_selected_records": len(indices),
        "missing_id_count": 0,
        "duplicate_id_count": duplicate_ids,
        "cohort_manifest_fold_mismatch_count": cohort_manifest_fold_mismatches,
        "fold_counts": {str(k): v for k, v in fold_counts.items()},
        "known_jointly_excluded_ids": sorted(known_excluded),
        "known_jointly_excluded_ids_present": known_excluded_present,
    }
    subset = RawCohort(
        ids=select_list(cohort.ids),
        smiles=select_list(cohort.smiles),
        graphs=select_list(cohort.graphs),
        field_dicts=select_list(cohort.field_dicts),
        targets=cohort.targets[indices].copy(),
        target_names=list(cohort.target_names),
        metals=select_list(cohort.metals),
        coordination_numbers=select_list(cohort.coordination_numbers),
        shapes=select_list(cohort.shapes),
        coordrep_hashes=select_list(cohort.coordrep_hashes),
        audit=audit,
    )
    return subset, subset_manifest


def load_pair_manifest(
    manifest: pd.DataFrame,
    args: argparse.Namespace,
    source_manifest: pd.DataFrame | None = None,
) -> Tuple[pd.DataFrame | None, Dict[str, object]]:
    if args.development_derive_pairs_from_members:
        return None, {
            "status": "DEVELOPMENT_DERIVED_FROM_MEMBERS",
            "warning": "Pair levels are not frozen; do not use for formal analysis.",
        }
    if not args.pair_manifest.exists():
        raise FileNotFoundError(
            f"Frozen pair manifest not found: {args.pair_manifest}. "
            "Formal runs require pairs.csv."
        )
    pairs = pd.read_csv(
        args.pair_manifest,
        dtype={
            "pair_id": str,
            "graph_sha256": str,
            "left_mol_id": str,
            "right_mol_id": str,
            "fold": int,
        },
    )
    required = {
        "pair_id",
        "graph_sha256",
        "fold",
        "left_mol_id",
        "right_mol_id",
        "is_g0",
        "is_s1",
        "is_s2",
        "is_context_exact",
        "is_s3",
    }
    if not required.issubset(pairs.columns):
        raise ValueError(f"Pair manifest needs columns {sorted(required)}")
    if pairs["pair_id"].duplicated().any():
        raise ValueError("Pair manifest has duplicate pair IDs")
    unordered_keys = pairs.apply(
        lambda row: "\0".join(sorted((row["left_mol_id"], row["right_mol_id"]))),
        axis=1,
    )
    if unordered_keys.duplicated().any():
        raise ValueError("Pair manifest has duplicate unordered molecule pairs")

    flag_columns = ["is_g0", "is_s1", "is_s2", "is_context_exact", "is_s3"]
    for column in flag_columns:
        pairs[column] = pairs[column].astype(bool)
    nesting_violations = int(
        (
            (~pairs["is_g0"])
            | (pairs["is_s1"] & ~pairs["is_g0"])
            | (pairs["is_s2"] & ~pairs["is_s1"])
            | (pairs["is_context_exact"] & ~pairs["is_s2"])
            | (pairs["is_s3"] & ~pairs["is_context_exact"])
        ).sum()
    )
    if nesting_violations:
        raise ValueError(f"Pair-level nesting violations: {nesting_violations}")

    source_records = (
        source_manifest if source_manifest is not None else manifest
    ).set_index("mol_id", drop=False)
    endpoint_ids = set(pairs["left_mol_id"]) | set(pairs["right_mol_id"])
    missing_source_endpoints = endpoint_ids - set(source_records.index)
    if missing_source_endpoints:
        raise ValueError(
            f"Pair manifest has {len(missing_source_endpoints)} endpoints outside "
            "the source record manifest"
        )
    source_left_fold = pairs["left_mol_id"].map(source_records["fold"])
    source_right_fold = pairs["right_mol_id"].map(source_records["fold"])
    source_left_graph = pairs["left_mol_id"].map(source_records["graph_sha256"])
    source_right_graph = pairs["right_mol_id"].map(source_records["graph_sha256"])
    source_fold_mismatches = int(
        (
            (source_left_fold != pairs["fold"])
            | (source_right_fold != pairs["fold"])
        ).sum()
    )
    source_graph_mismatches = int(
        (
            (source_left_graph != pairs["graph_sha256"])
            | (source_right_graph != pairs["graph_sha256"])
        ).sum()
    )
    if source_fold_mismatches or source_graph_mismatches:
        raise ValueError(
            "Pair/source-record manifest mismatch: "
            f"fold={source_fold_mismatches}, "
            f"graph_sha256={source_graph_mismatches}"
        )

    source_fold_counts = {
        int(k): int(v) for k, v in pairs["fold"].value_counts().sort_index().items()
    }
    source_level_counts = {
        column.removeprefix("is_").upper(): int(pairs[column].sum())
        for column in flag_columns
    }
    formal_path = args.pair_manifest.resolve() == FORMAL_PAIRS_MANIFEST.resolve()
    if formal_path and source_fold_counts != FORMAL_PAIR_FOLD_SIZES:
        raise ValueError(
            f"Formal pair counts changed: {source_fold_counts}; "
            f"expected {FORMAL_PAIR_FOLD_SIZES}"
        )

    selected_ids = set(manifest["mol_id"].astype(str))
    endpoint_intersection = pairs["left_mol_id"].isin(selected_ids) & pairs[
        "right_mol_id"
    ].isin(selected_ids)
    source_pair_count = len(pairs)
    pairs = pairs.loc[endpoint_intersection].reset_index(drop=True).copy()
    records = manifest.set_index("mol_id", drop=False)
    left_fold = pairs["left_mol_id"].map(records["fold"])
    right_fold = pairs["right_mol_id"].map(records["fold"])
    left_graph = pairs["left_mol_id"].map(records["graph_sha256"])
    right_graph = pairs["right_mol_id"].map(records["graph_sha256"])
    fold_mismatches = int(
        ((left_fold != pairs["fold"]) | (right_fold != pairs["fold"])).sum()
    )
    graph_mismatches = int(
        (
            (left_graph != pairs["graph_sha256"])
            | (right_graph != pairs["graph_sha256"])
        ).sum()
    )
    if fold_mismatches or graph_mismatches:
        raise ValueError(
            "Pair/record manifest mismatch: "
            f"fold={fold_mismatches}, graph_sha256={graph_mismatches}"
        )

    fold_counts = {
        int(k): int(v) for k, v in pairs["fold"].value_counts().sort_index().items()
    }
    audit = {
        "status": "FORMAL_FROZEN_PAIR_LEVELS" if formal_path else "CUSTOM_PAIR_MANIFEST",
        "path": str(args.pair_manifest.resolve()),
        "sha256": file_sha256(args.pair_manifest),
        "pair_endpoint_policy": "retain only pairs with both endpoints in analysis cohort",
        "n_source_pairs": source_pair_count,
        "n_pairs": len(pairs),
        "n_pairs_removed_by_cohort_intersection": source_pair_count - len(pairs),
        "source_fold_counts": {
            str(k): v for k, v in source_fold_counts.items()
        },
        "fold_counts": {str(k): v for k, v in fold_counts.items()},
        "source_level_counts": source_level_counts,
        "level_counts": {
            column.removeprefix("is_").upper(): int(pairs[column].sum())
            for column in flag_columns
        },
        "duplicate_pair_id_count": 0,
        "duplicate_unordered_pair_count": 0,
        "missing_source_endpoint_count": 0,
        "source_fold_mismatch_count": source_fold_mismatches,
        "source_graph_sha256_mismatch_count": source_graph_mismatches,
        "fold_mismatch_count": fold_mismatches,
        "graph_sha256_mismatch_count": graph_mismatches,
        "nesting_violation_count": nesting_violations,
    }
    return pairs, audit


def make_fold_data(
    cohort: RawCohort, manifest: pd.DataFrame, fold: int, args: argparse.Namespace
) -> Tuple[PreparedData, Dict[str, object]]:
    outer_test = np.flatnonzero(manifest["fold"].to_numpy() == fold)
    outer_train = np.flatnonzero(manifest["fold"].to_numpy() != fold)
    graph_hashes = manifest["graph_sha256"].astype(str).to_numpy()
    outer_groups = graph_hashes[outer_train]
    inner = GroupShuffleSplit(
        n_splits=1,
        test_size=args.validation_fraction,
        random_state=args.split_seed + 1000 + fold,
    )
    train_local, val_local = next(inner.split(outer_train, groups=outer_groups))
    train_idx = outer_train[train_local]
    val_idx = outer_train[val_local]
    test_idx = outer_test

    vectorizer = DictVectorizer(sparse=False, sort=True)
    vectorizer.fit([cohort.field_dicts[i] for i in train_idx])
    field_raw = vectorizer.transform(cohort.field_dicts).astype(np.float32)
    scaler = StandardScaler()
    scaler.fit(field_raw[train_idx])
    fields = scaler.transform(field_raw).astype(np.float32)

    target_mean = cohort.targets[train_idx].mean(axis=0)
    target_scale = cohort.targets[train_idx].std(axis=0)
    target_scale[target_scale < 1.0e-12] = 1.0
    y_scaled = (cohort.targets - target_mean) / target_scale
    for i, graph in enumerate(cohort.graphs):
        graph.coord_x = torch.from_numpy(fields[i : i + 1])
        graph.y = torch.from_numpy(y_scaled[i : i + 1].astype(np.float32))

    split_graphs = {
        "train": set(graph_hashes[train_idx]),
        "validation": set(graph_hashes[val_idx]),
        "test": set(graph_hashes[test_idx]),
    }
    overlaps = {
        "train_validation": len(split_graphs["train"] & split_graphs["validation"]),
        "train_test": len(split_graphs["train"] & split_graphs["test"]),
        "validation_test": len(split_graphs["validation"] & split_graphs["test"]),
    }
    if any(overlaps.values()):
        raise RuntimeError(f"Fold {fold} graph leakage: {overlaps}")

    fold_audit: Dict[str, object] = {
        "fold": fold,
        "train_records": len(train_idx),
        "validation_records": len(val_idx),
        "test_records": len(test_idx),
        "field_dimension": int(fields.shape[1]),
        "field_names": [str(name) for name in vectorizer.get_feature_names_out()],
        "model_field_views": {},
        "graph_sha256_overlaps": overlaps,
        "vectorizer_fit_partition": "train_only",
        "field_scaler_fit_partition": "train_only",
        "target_scaler_fit_partition": "train_only",
        "target_train_mean": target_mean.tolist(),
        "target_train_scale": target_scale.tolist(),
    }
    prepared = PreparedData(
        ids=cohort.ids,
        smiles=cohort.smiles,
        graphs=cohort.graphs,
        fields=fields,
        targets=cohort.targets,
        target_names=cohort.target_names,
        field_names=list(vectorizer.get_feature_names_out()),
        train_idx=np.asarray(train_idx),
        val_idx=np.asarray(val_idx),
        test_idx=np.asarray(test_idx),
        target_mean=target_mean,
        target_scale=target_scale,
        field_scaler_mean=scaler.mean_,
        field_scaler_scale=scaler.scale_,
        audit=fold_audit,
    )
    return prepared, fold_audit


def member_rows(
    cohort: RawCohort,
    manifest: pd.DataFrame,
    data: PreparedData,
    prediction_scaled: np.ndarray,
    truth_scaled: np.ndarray,
    fold: int,
    seed: int,
    model_label: str,
) -> List[Dict[str, object]]:
    prediction = prediction_scaled * data.target_scale + data.target_mean
    truth = truth_scaled * data.target_scale + data.target_mean
    rows: List[Dict[str, object]] = []
    for local, index in enumerate(data.test_idx):
        frozen = manifest.iloc[index]
        row: Dict[str, object] = {
            "model": model_label,
            "seed": seed,
            "fold": fold,
            "mol_id": cohort.ids[index],
            "graph_sha256": str(frozen["graph_sha256"]),
            "metal": str(frozen["metal"]),
            "cn": int(frozen["cn"]),
            "shape": str(frozen["best_shape"]),
            "context_sha256": str(frozen["context_sha256"]),
            "coordrep_no_v_sha256": str(frozen["coordrep_no_v_sha256"]),
            "coordrep_full_sha256": cohort.coordrep_hashes[index],
        }
        for column, target in enumerate(cohort.target_names):
            row[f"{target}_true"] = float(truth[local, column])
            row[f"{target}_pred"] = float(prediction[local, column])
        rows.append(row)
    return rows


def pair_detail_rows(
    rows: Sequence[Dict[str, object]],
    pair_manifest: pd.DataFrame | None,
    args: argparse.Namespace,
) -> List[Dict[str, object]]:
    frame = pd.DataFrame(rows)
    if frame.empty:
        return []
    details: List[Dict[str, object]] = []
    thresholds = {
        "hl_gap_ev": args.gap_min_delta,
        "dipole_moment": args.dipole_min_delta,
    }
    members_by_id = {str(row["mol_id"]): row for row in rows}
    fold = int(frame["fold"].iloc[0])
    if pair_manifest is not None:
        pair_records = pair_manifest.loc[pair_manifest["fold"] == fold].to_dict("records")
        pair_source = "FROZEN_PAIR_MANIFEST"
    else:
        pair_records = []
        for graph_sha256, group in frame.groupby("graph_sha256", sort=True):
            members = group.sort_values("mol_id").to_dict("records")
            for left, right in combinations(members, 2):
                is_s1 = (
                    left["metal"] == right["metal"]
                    and int(left["cn"]) == int(right["cn"])
                )
                pair_records.append(
                    {
                        "pair_id": full_sha256_text(
                            f"{graph_sha256}|{left['mol_id']}|{right['mol_id']}"
                        ),
                        "graph_sha256": graph_sha256,
                        "fold": fold,
                        "left_mol_id": left["mol_id"],
                        "right_mol_id": right["mol_id"],
                        "is_g0": True,
                        "is_s1": is_s1,
                        # Exact S2/context/S3 definitions need the frozen manifest.
                        "is_s2": False,
                        "is_context_exact": False,
                        "is_s3": False,
                    }
                )
        pair_source = "DEVELOPMENT_DERIVED_FROM_MEMBERS"

    frozen_delta_columns = {
        "hl_gap_ev": "delta_hl_gap_ev",
        "dipole_moment": "delta_dipole_moment_d",
    }
    for pair in pair_records:
        left = members_by_id.get(str(pair["left_mol_id"]))
        right = members_by_id.get(str(pair["right_mol_id"]))
        if left is None or right is None:
            raise RuntimeError(
                f"Fold {fold} pair {pair['pair_id']} has an endpoint outside test members"
            )
        graph_sha256 = str(pair["graph_sha256"])
        if (
            left["graph_sha256"] != graph_sha256
            or right["graph_sha256"] != graph_sha256
        ):
            raise RuntimeError(f"Pair {pair['pair_id']} graph SHA256 mismatch")
        coordrep_diff = left["coordrep_full_sha256"] != right["coordrep_full_sha256"]
        shape_diff = left["shape"] != right["shape"]
        cn_diff = int(left["cn"]) != int(right["cn"])
        sp_td = {str(left["shape"]), str(right["shape"])} == {"SP", "Td"}
        for target_name, threshold in thresholds.items():
            true_column = f"{target_name}_true"
            if true_column not in frame:
                continue
            # The released pair contract orients signed deltas as left minus right.
            true_delta = float(left[true_column] - right[true_column])
            pred_delta = float(
                left[f"{target_name}_pred"] - right[f"{target_name}_pred"]
            )
            frozen_delta_column = frozen_delta_columns[target_name]
            if pair_manifest is not None and frozen_delta_column in pair:
                frozen_delta = float(pair[frozen_delta_column])
                if not np.isclose(true_delta, frozen_delta, rtol=0.0, atol=2.0e-5):
                    raise RuntimeError(
                        f"Pair {pair['pair_id']} {target_name} true delta disagrees "
                        f"with pairs.csv: loaded={true_delta}, frozen={frozen_delta}"
                    )
            else:
                frozen_delta = np.nan
            pred_tie = abs(pred_delta) <= args.prediction_tie_tolerance
            material = abs(true_delta) > threshold
            if not material:
                primary_concordance = np.nan
                non_tie_direction_correct = np.nan
            elif pred_tie:
                primary_concordance = 0.5
                non_tie_direction_correct = np.nan
            else:
                correct = bool(np.sign(true_delta) == np.sign(pred_delta))
                primary_concordance = 1.0 if correct else 0.0
                non_tie_direction_correct = correct
            details.append(
                {
                    "model": left["model"],
                    "seed": left["seed"],
                    "fold": fold,
                    "pair_id": str(pair["pair_id"]),
                    "pair_source": pair_source,
                    "graph_sha256": graph_sha256,
                    "left_mol_id": left["mol_id"],
                    "right_mol_id": right["mol_id"],
                    "left_metal": left["metal"],
                    "right_metal": right["metal"],
                    "left_cn": left["cn"],
                    "right_cn": right["cn"],
                    "left_shape": left["shape"],
                    "right_shape": right["shape"],
                    "is_g0": bool(pair["is_g0"]),
                    "is_s1": bool(pair["is_s1"]),
                    "is_s2": bool(pair["is_s2"]),
                    "is_context_exact": bool(pair["is_context_exact"]),
                    "is_s3": bool(pair["is_s3"]),
                    "coordrep_diff": coordrep_diff,
                    "shape_diff": shape_diff,
                    "cn_diff": cn_diff,
                    "sp_td": sp_td,
                    "target": target_name,
                    "true_left": left[true_column],
                    "true_right": right[true_column],
                    "pred_left": left[f"{target_name}_pred"],
                    "pred_right": right[f"{target_name}_pred"],
                    "true_delta_left_minus_right": true_delta,
                    "frozen_true_delta_left_minus_right": frozen_delta,
                    "pred_delta_left_minus_right": pred_delta,
                    "delta_abs_error": abs(true_delta - pred_delta),
                    "prediction_is_tie": pred_tie,
                    "prediction_tie_tolerance": args.prediction_tie_tolerance,
                    "material_difference_threshold": threshold,
                    "material_difference_strictly_greater": material,
                    "primary_concordance_ties_half": primary_concordance,
                    "non_tie_direction_correct": non_tie_direction_correct,
                }
            )
    return details


def pair_metric_rows(
    details: Sequence[Dict[str, object]], members: Sequence[Dict[str, object]]
) -> List[Dict[str, object]]:
    pair_frame = pd.DataFrame(details)
    member_frame = pd.DataFrame(members)
    if pair_frame.empty:
        return []
    subsets = {
        "G0": pair_frame["is_g0"].astype(bool).to_numpy(),
        "S1": pair_frame["is_s1"].astype(bool).to_numpy(),
        "S2": pair_frame["is_s2"].astype(bool).to_numpy(),
        "CONTEXT_EXACT": pair_frame["is_context_exact"].astype(bool).to_numpy(),
        "S3": pair_frame["is_s3"].astype(bool).to_numpy(),
        "coordrep_diff": pair_frame["coordrep_diff"].astype(bool).to_numpy(),
        "shape_diff": pair_frame["shape_diff"].astype(bool).to_numpy(),
        "sp_td": pair_frame["sp_td"].astype(bool).to_numpy(),
        "cn_diff": pair_frame["cn_diff"].astype(bool).to_numpy(),
    }
    metrics: List[Dict[str, object]] = []
    for subset, mask in subsets.items():
        selected = pair_frame.loc[mask]
        if selected.empty:
            continue
        for target, group in selected.groupby("target", sort=True):
            material = group[group["material_difference_strictly_greater"].astype(bool)]
            non_tie_material = material[~material["prediction_is_tie"].astype(bool)]
            graph_hashes = set(group["graph_sha256"])
            selected_members = member_frame[
                member_frame["graph_sha256"].isin(graph_hashes)
            ].copy()
            true_column = f"{target}_true"
            pred_column = f"{target}_pred"
            selected_members["true_centered"] = selected_members[true_column] - selected_members.groupby(
                "graph_sha256"
            )[true_column].transform("mean")
            selected_members["pred_centered"] = selected_members[pred_column] - selected_members.groupby(
                "graph_sha256"
            )[pred_column].transform("mean")
            centered_mae = float(
                np.mean(
                    np.abs(
                        selected_members["true_centered"]
                        - selected_members["pred_centered"]
                    )
                )
            )
            metrics.append(
                {
                    "model": str(group["model"].iloc[0]),
                    "seed": int(group["seed"].iloc[0]),
                    "fold": int(group["fold"].iloc[0]),
                    "subset": subset,
                    "target": target,
                    "n_pairs": len(group),
                    "n_graph_groups": int(group["graph_sha256"].nunique()),
                    "n_member_records": int(selected_members["mol_id"].nunique()),
                    "delta_mae_all_pairs": float(group["delta_abs_error"].mean()),
                    "true_abs_delta_mean": float(
                        group["true_delta_left_minus_right"].abs().mean()
                    ),
                    "true_abs_delta_median": float(
                        group["true_delta_left_minus_right"].abs().median()
                    ),
                    "within_graph_centered_mae": centered_mae,
                    "prediction_tie_rate": float(group["prediction_is_tie"].mean()),
                    "n_material_pairs": len(material),
                    "material_prediction_tie_rate": (
                        float(material["prediction_is_tie"].mean())
                        if len(material)
                        else np.nan
                    ),
                    "primary_concordance_ties_half": (
                        float(material["primary_concordance_ties_half"].mean())
                        if len(material)
                        else np.nan
                    ),
                    "n_non_tie_material_pairs": len(non_tie_material),
                    "non_tie_direction_accuracy": (
                        float(non_tie_material["non_tie_direction_correct"].mean())
                        if len(non_tie_material)
                        else np.nan
                    ),
                }
            )
    return metrics


def member_metric_rows(rows: Sequence[Dict[str, object]]) -> List[Dict[str, object]]:
    frame = pd.DataFrame(rows)
    metrics: List[Dict[str, object]] = []
    for target in ("hl_gap_ev", "dipole_moment"):
        if f"{target}_true" not in frame:
            continue
        truth = frame[f"{target}_true"].to_numpy(float)
        prediction = frame[f"{target}_pred"].to_numpy(float)
        metrics.append(
            {
                "model": str(frame["model"].iloc[0]),
                "seed": int(frame["seed"].iloc[0]),
                "fold": int(frame["fold"].iloc[0]),
                "target": target,
                "n_test": len(frame),
                "mae": float(mean_absolute_error(truth, prediction)),
                "rmse": float(mean_squared_error(truth, prediction) ** 0.5),
                "r2": float(r2_score(truth, prediction)),
            }
        )
    return metrics


def aggregate_metrics(rows: Sequence[Dict[str, object]], keys: Sequence[str]) -> List[Dict[str, object]]:
    frame = pd.DataFrame(rows)
    if frame.empty:
        return []
    numeric = [
        column
        for column in frame.select_dtypes(include=[np.number]).columns
        if column not in {"seed", "fold"}
    ]
    result: List[Dict[str, object]] = []
    for labels, group in frame.groupby(list(keys), sort=True, dropna=False):
        labels = labels if isinstance(labels, tuple) else (labels,)
        entry = {key: value for key, value in zip(keys, labels)}
        entry["n_fold_seed_runs"] = len(group)
        for column in numeric:
            values = group[column].to_numpy(float)
            finite = values[np.isfinite(values)]
            entry[f"{column}_mean"] = float(finite.mean()) if len(finite) else np.nan
            entry[f"{column}_sd"] = (
                float(finite.std(ddof=1)) if len(finite) > 1 else 0.0
            )
        result.append(entry)
    return result


def main() -> None:
    args = parse_args()
    if args.smoke:
        args.folds = [0]
        args.seeds = [args.seeds[0]]
        args.epochs = 1
        args.patience = 1
    field_modes = parse_field_modes(args.field_mode)
    allowed_models = {"gine_wide", "hybrid"} | set(field_modes)
    unknown_models = sorted(set(args.models) - allowed_models)
    if unknown_models:
        raise ValueError(f"Unknown models: {unknown_models}; choose from {sorted(allowed_models)}")
    selected_folds = args.folds if args.folds is not None else list(range(args.n_folds))
    invalid_folds = sorted(set(selected_folds) - set(range(args.n_folds)))
    if invalid_folds:
        raise ValueError(f"Invalid --folds: {invalid_folds}")

    args.out.mkdir(parents=True, exist_ok=True)
    started = time.time()
    print("Loading matched tmQM/CoordRep cohort...", flush=True)
    cohort = load_raw_cohort(args)
    full_manifest = load_fold_manifest(cohort, args)
    cohort, manifest = restrict_analysis_cohort(cohort, full_manifest, args)
    pair_manifest, pair_manifest_audit = load_pair_manifest(
        manifest, args, source_manifest=full_manifest
    )
    cohort.audit["pair_manifest_audit"] = pair_manifest_audit
    manifest_columns = [
        column
        for column in (
            "mol_id",
            "graph_sha256",
            "fold",
            "metal",
            "cn",
            "best_shape",
            "context_sha256",
            "coordrep_no_v_sha256",
            "hl_gap_ev",
            "dipole_moment_d",
        )
        if column in manifest
    ]
    manifest[manifest_columns].to_csv(
        args.out / "fold_manifest_used.csv", index=False
    )
    if pair_manifest is not None:
        pair_manifest.to_csv(args.out / "pair_manifest_used.csv", index=False)
    with (args.out / "dataset_audit.json").open("w", encoding="utf-8") as handle:
        json.dump(cohort.audit, handle, indent=2)

    device = torch.device(args.device)
    all_member_rows: List[Dict[str, object]] = []
    all_pair_rows: List[Dict[str, object]] = []
    all_pair_metrics: List[Dict[str, object]] = []
    all_member_metrics: List[Dict[str, object]] = []
    fold_audits: List[Dict[str, object]] = []
    run_details: Dict[str, object] = {}

    for fold in selected_folds:
        print(f"Preparing outer fold {fold}...", flush=True)
        data, fold_audit = make_fold_data(cohort, manifest, fold, args)
        fold_audits.append(fold_audit)
        for seed in args.seeds:
            for model_label in args.models:
                internal_name = "gin_wide" if model_label == "gine_wide" else "hybrid"
                key = f"fold{fold}_{model_label}_seed{seed}"
                field_view_audit = apply_model_field_view(
                    data, model_label, field_modes
                )
                prior_field_view = fold_audit["model_field_views"].get(model_label)
                if prior_field_view is not None and prior_field_view != field_view_audit:
                    raise RuntimeError(
                        f"Field view changed across seeds for fold {fold}, {model_label}"
                    )
                fold_audit["model_field_views"][model_label] = field_view_audit
                print(f"Training {key}...", flush=True)
                run_start = time.time()
                model, training_metadata = train_graph_model(
                    internal_name, data, args, seed, device
                )
                test_loader = graph_loader(
                    data.graphs,
                    data.test_idx,
                    args.batch_size * 2,
                    False,
                    seed,
                    args.num_workers,
                )
                pred_scaled, truth_scaled = predict_graph(model, test_loader, device)
                members = member_rows(
                    cohort,
                    manifest,
                    data,
                    pred_scaled,
                    truth_scaled,
                    fold,
                    seed,
                    model_label,
                )
                pairs = pair_detail_rows(members, pair_manifest, args)
                pair_metrics = pair_metric_rows(pairs, members)
                all_member_rows.extend(members)
                all_pair_rows.extend(pairs)
                all_pair_metrics.extend(pair_metrics)
                all_member_metrics.extend(member_metric_rows(members))
                run_details[key] = {
                    **training_metadata,
                    "model_label": model_label,
                    "internal_architecture": internal_name,
                    "field_view": field_view_audit,
                    "paired_initialization_seed": seed,
                    "elapsed_seconds": time.time() - run_start,
                    "n_member_predictions": len(members),
                    "n_pair_target_rows": len(pairs),
                }
                del model
                if torch.cuda.is_available():
                    torch.cuda.empty_cache()

    write_csv(args.out / "member_predictions.csv", all_member_rows)
    write_csv(args.out / "same_graph_pair_details.csv", all_pair_rows)
    write_csv(args.out / "pair_metrics_by_fold_seed.csv", all_pair_metrics)
    write_csv(args.out / "member_metrics_by_fold_seed.csv", all_member_metrics)
    write_csv(
        args.out / "pair_metrics_summary.csv",
        aggregate_metrics(all_pair_metrics, ["model", "subset", "target"]),
    )
    write_csv(
        args.out / "member_metrics_summary.csv",
        aggregate_metrics(all_member_metrics, ["model", "target"]),
    )
    parameter_control = audit_hybrid_parameter_control(
        run_details,
        selected_folds,
        args.seeds,
        args.models,
        field_modes,
    )
    with (args.out / "fold_audits.json").open("w", encoding="utf-8") as handle:
        json.dump(fold_audits, handle, indent=2)
    metadata = {
        "arguments": {
            key: str(value) if isinstance(value, Path) else value
            for key, value in vars(args).items()
        },
        "device": str(device),
        "torch_version": torch.__version__,
        "cuda_available": torch.cuda.is_available(),
        "selected_folds": selected_folds,
        "fold_manifest_sha256": file_sha256(args.fold_manifest),
        "pair_manifest_sha256": (
            file_sha256(args.pair_manifest) if pair_manifest is not None else None
        ),
        "formal_fold_sizes": {str(k): v for k, v in FORMAL_FOLD_SIZES.items()},
        "analysis_fold_sizes": {
            str(int(k)): int(v)
            for k, v in manifest["fold"].value_counts().sort_index().items()
        },
        "primary_evidence_unit": (
            "molecule-level out-of-fold metrics on the identical analysis cohort; "
            "same-graph pair metrics are secondary diagnostics"
        ),
        "prediction_tie_policy": (
            "For abs(true delta) strictly greater than the material threshold, "
            "abs(predicted delta) <= tolerance is a tie and receives 0.5 in the "
            "primary concordance. Non-tie direction accuracy is reported separately."
        ),
        "material_difference_thresholds": {
            "hl_gap_ev": args.gap_min_delta,
            "dipole_moment_d": args.dipole_min_delta,
            "comparison": "strictly_greater_than",
        },
        "field_quality_ablation": {
            "built_in_modes": {
                label: list(selectors)
                for label, selectors in BUILTIN_HYBRID_FIELD_MASKS.items()
            },
            "resolved_modes": {
                label: list(selectors) for label, selectors in field_modes.items()
            },
            "selected_models": list(args.models),
            "masking_stage": (
                "after DictVectorizer and train-only StandardScaler; before graph "
                "attachment and model training"
            ),
            "masked_partitions": ["train", "validation", "test"],
            "parameter_control": parameter_control,
        },
        "total_elapsed_seconds": time.time() - started,
        "runs": run_details,
    }
    with (args.out / "run_metadata.json").open("w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2)
    print(f"Done. Outputs: {args.out}", flush=True)


if __name__ == "__main__":
    main()
