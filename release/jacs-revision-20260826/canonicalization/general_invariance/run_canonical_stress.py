#!/usr/bin/env python
"""Reproducible CoordRep rotation/permutation invariance stress test.

The publication cohort is a fixed, graph-disjoint, fold-by-CN stratified
CN4--6 sample from the tmQMg/PBE overlap.  It is *not* a prevalence sample
from the April 2025 CSD.  Four independent source ledgers are joined before
selection:

1. corrected CoordRep 1.1.2rc2 PBE records;
2. the matched-record graph/fold ledger;
3. the independent single-metal XYZ audit; and
4. the official individual tmQMg/PBE XYZ ZIP.

The script supports two canonicalizer modes:

``source``
    Use the supplied release tree without runtime modification.  This is the
    only mode permitted for a publication-scale run.

``minimal_smoke``
    Runtime diagnostic that removes raw ``attach_atoms`` from the duplicate-
    ligand equivalence key while retaining exact Cartesian-product search.
    This mode is deliberately restricted to the 100 x 20 gate.  It is not a
    publication result because a ligand-local canonical attachment/orbit key
    is required to distinguish non-equivalent binding sites safely.

The four requested comparison representations are evaluated from every
re-encoded raw molecule:

* raw ordered first-sphere Cartesian sequence;
* ordered first-sphere distance matrix;
* source-label-free connectivity-only record; and
* the full canonical CoordRep string.

An uncanonicalized CoordRep string is retained as an explicit
canonicalization-off control.  Three transformation arms are run independently:
rigid motion only, atom permutation only, and the combined transformation.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import platform
import shutil
import subprocess
import sys
import time
import zipfile
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from itertools import permutations, product
from math import factorial
from pathlib import Path
from typing import Any, Iterable

import numpy as np


SCRIPT = Path(__file__).resolve()
OUT = SCRIPT.parent
ROOT = SCRIPT.parents[3]

PUBLIC_REPO = ROOT / "tmp" / "github_audit_coordrep_20260806"
PUBLIC_RELEASE = PUBLIC_REPO / "release" / "jacs-revision-20260819"
CANDIDATE_FIX = (
    ROOT
    / "rew20260817"
    / "figure2_redesign_20260822"
    / "canonical_fix_20260822"
)
RC2_RECORDS = (
    ROOT
    / "revision_experiments"
    / "results"
    / "coordrep_1_1_2_corrected_cshm_pbe31121_rc2_v2"
    / "records.jsonl"
)
MATCHED_RECORDS = (
    ROOT
    / "revision_experiments"
    / "results"
    / "trex_coordrep_conflict_audit_v1"
    / "matched_records.csv"
)
SINGLE_METAL_AUDIT = (
    ROOT
    / "revision_experiments"
    / "results"
    / "frozen_coordrep_pbe_serializer_full31121_v1"
    / "record_audit.csv"
)
PBE_ZIP = (
    ROOT
    / "revision_experiments"
    / "source_cache"
    / "official_tmqmg_github_xyz"
    / "tmQMg_xyz.zip"
)
PBE_ZIP_MANIFEST = PBE_ZIP.with_name("SOURCE_MANIFEST.json")
FAC_MOL2 = (
    ROOT
    / "rew20260817"
    / "Fig2_Fig5_FullCSD_Redraw_20260818"
    / "structures"
    / "EBAGAR.mol2"
)
MER_MOL2 = FAC_MOL2.with_name("EBAGEV.mol2")

EXPECTED_PUBLIC_COMMIT = "8a51e6eac370fd5ff20d0b26df834e2b61f5e90b"
EXPECTED_PUBLIC_CANONICALIZER_SHA256 = (
    "eeee838333f78e0f1f16a3f2c51a009a076640dc28e83ae31f51703a23aa8f7f"
)
EXPECTED_PUBLIC_ENCODER_SHA256 = (
    "cbfb61492de3cb8853deb0eb24c7a0182fb455ca5d94f30da8f10f70a014d9a5"
)
EXPECTED_PBE_ZIP_SHA256 = (
    "e0d15a70bcba294717cd9f9792e7fac99ef0c5c61c3a6e08dcc8a8643f53660a"
)
EXPECTED_CANDIDATE_SOURCE_SHA256 = {
    "coordrep/encode.py": "cbfb61492de3cb8853deb0eb24c7a0182fb455ca5d94f30da8f10f70a014d9a5",
    "coordrep/canonical/canonicalize.py": "753f13ed93dd9b41221097f72cd8f2a037c27bd3c202097a07df5c1b91866967",
    "coordrep/core.py": "6bc9c204c0ba53f08a09003dda777d161eb82ad741e9929489ebfee4db1bd08f",
    "coordrep/graph/ligand_module.py": "4320ed2e0ca7565c35d231cb7216b7200e3f2664b86a5c2655c106783854949a",
    "coordrep/graph/donor_sites.py": "0461a81a250f3f80f4d8110773420cd98cc4f0f2df36011b225baf3178e82648",
    "coordrep/geometry/rel_config.py": "9854495bc8b836576733137f7bcbd32b83d97182dbbdfe24a14930c130a78748",
    "coordrep/serialize/to_string.py": "9bf507cad69d7de95adb1201a4117704e2bc5382bd23b616176e79de8dbfdc43",
}

SELECTION_SEED = 20260822
TRANSFORM_SEED = 930_220_260
FOLDS = tuple(range(5))
MAIN_CNS = (4, 5, 6)
# Every fold contributes exactly 200 records.  Rotating the 67/67/66
# allocation gives the locked overall CN4/CN5/CN6 counts 334/333/333.
CN_QUOTA_BY_FOLD = {
    0: {4: 67, 5: 67, 6: 66},
    1: {4: 67, 5: 66, 6: 67},
    2: {4: 66, 5: 67, 6: 67},
    3: {4: 67, 5: 67, 6: 66},
    4: {4: 67, 5: 66, 6: 67},
}
# The N=100 gate mirrors that design at one tenth of the full allocation.
GATE_QUOTA_BY_FOLD = {
    0: {4: 7, 5: 7, 6: 6},
    1: {4: 7, 5: 6, 6: 7},
    2: {4: 6, 5: 7, 6: 7},
    3: {4: 7, 5: 7, 6: 6},
    4: {4: 7, 5: 6, 6: 7},
}
EXPECTED_N = sum(sum(quota.values()) for quota in CN_QUOTA_BY_FOLD.values())
EXPECTED_GATE_N = sum(sum(quota.values()) for quota in GATE_QUOTA_BY_FOLD.values())
EXPECTED_ELIGIBLE_UNIQUE_GRAPHS = 6_427
ARRAY_SERIALIZATION_DECIMALS = 6
ARMS = ("rigid", "permutation", "combined")
ARM_CODE = {"rigid": 1, "permutation": 2, "combined": 3}
REPRESENTATIONS = (
    "raw_cartesian",
    "ordered_distance_matrix",
    "connectivity_only",
    "uncanonicalized_coordrep",
    "patched_coordrep",
)
PUBLICATION_REPRESENTATIONS = (
    "raw_cartesian",
    "ordered_distance_matrix",
    "connectivity_only",
    "patched_coordrep",
)

SUMMARY_FIELDS = [
    "run_stage",
    "selection_order",
    "mol_id",
    "fold",
    "metal",
    "cn",
    "graph_sha256",
    "arm",
    "k",
    "representation",
    "reference_sha256",
    "matched_variants",
    "mismatched_variants",
    "match_fraction",
    "unique_variant_hashes",
    "unique_hashes_with_reference",
    "all_variants_match_reference",
]

MISMATCH_FIELDS = [
    "run_stage",
    "selection_order",
    "mol_id",
    "fold",
    "metal",
    "cn",
    "arm",
    "variant_index",
    "representation",
    "reference_sha256",
    "variant_sha256",
    "reference_value",
    "variant_value",
    "failure_stage",
    "error",
]

_WORKER: dict[str, Any] = {}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release-root", type=Path, default=PUBLIC_RELEASE)
    parser.add_argument(
        "--patch-mode",
        choices=("source", "minimal_smoke"),
        default="minimal_smoke",
    )
    parser.add_argument(
        "--run-full",
        action="store_true",
        help="After the gate, run N=1000 x K=100 per arm.",
    )
    parser.add_argument("--workers", type=int, default=min(8, os.cpu_count() or 1))
    parser.add_argument("--gate-n", type=int, default=100)
    parser.add_argument("--gate-k", type=int, default=20)
    parser.add_argument("--full-k", type=int, default=100)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Defaults to this script directory for source mode and to minimal_fix_smoke/ for diagnostic mode.",
    )
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def deterministic_score(namespace: str, *values: object) -> str:
    text = "|".join([str(SELECTION_SEED), namespace, *(str(v) for v in values)])
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def truth(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes"}


def write_csv(path: Path, fieldnames: list[str], rows: Iterable[dict[str, Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def git_value(repo: Path, *args: str) -> str:
    try:
        return subprocess.check_output(
            ["git", "-C", str(repo), *args],
            text=True,
            encoding="utf-8",
            stderr=subprocess.DEVNULL,
        ).strip()
    except Exception:
        return ""


def load_joined_candidates() -> tuple[list[dict[str, Any]], dict[str, Any]]:
    matched: dict[str, dict[str, str]] = {}
    with MATCHED_RECORDS.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            matched[row["mol_id"]] = {
                "mol_id": row["mol_id"],
                "graph_sha256": row["graph_sha256"],
                "fold": row["fold"],
                "metal": row["metal"],
                "matched_cn": row["cn"],
                "canonical_smiles": row.get("canonical_smiles", ""),
                "official_status": row.get("official_status", ""),
            }

    single: dict[str, dict[str, str]] = {}
    with SINGLE_METAL_AUDIT.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            single[row["mol_id"]] = {
                "n_transition_metals_in_xyz": row["n_transition_metals_in_xyz"],
                "official_xyz_sha256": row["official_xyz_sha256"],
                "official_xyz_bytes": row["official_xyz_bytes"],
            }

    audit = Counter()
    candidates: list[dict[str, Any]] = []
    n_rc2 = 0
    with RC2_RECORDS.open(encoding="utf-8") as handle:
        for line in handle:
            n_rc2 += 1
            record = json.loads(line)
            mol_id = str(record["mol_id"])
            audit["rc2_records"] += 1
            if mol_id not in matched or mol_id not in single:
                audit["excluded_missing_four_source_join"] += 1
                continue
            audit["four_source_ledger_join"] += 1

            if not all(
                truth(record.get(field))
                for field in ("parser_ok", "encode_ok", "serialize_ok", "identity_ok")
            ):
                audit["excluded_pipeline_status"] += 1
                continue
            audit["pipeline_status_pass"] += 1

            # Locked strict-supported-state gate.  Warnings are not silently
            # treated as successes in this benchmark; any rc2 validation issue
            # (principally UNUSUAL_GEOMETRY) excludes the record.
            if list(record.get("validation_issues") or []):
                audit["excluded_nonempty_validation_issues"] += 1
                continue
            audit["empty_validation_issues_pass"] += 1

            if single[mol_id]["n_transition_metals_in_xyz"] != "1":
                audit["excluded_not_single_metal_xyz"] += 1
                continue
            audit["single_metal_xyz_pass"] += 1

            cn = int(record.get("cn", -1))
            if cn not in MAIN_CNS:
                audit["excluded_cn_outside_4_6"] += 1
                continue
            audit["cn_4_6_pass"] += 1

            provenance = list(record.get("ligand_payload_provenance") or [])
            statuses = list(record.get("ligand_connectivity_status") or [])
            issue_codes = {
                str(item.get("code", ""))
                for item in (record.get("raw_encoder_issues") or [])
                if isinstance(item, dict)
            }
            if (
                not provenance
                or any(item != "SMILES" for item in provenance)
                or any("formula" in str(item).lower() for item in statuses)
                or "FORMULA_FALLBACK" in issue_codes
            ):
                audit["excluded_formula_or_non_smiles_payload"] += 1
                continue
            audit["smiles_payload_pass"] += 1

            match = matched[mol_id]
            graph = match["graph_sha256"].strip()
            fold_text = match["fold"].strip()
            if not graph or fold_text not in {str(v) for v in FOLDS}:
                audit["excluded_missing_graph_or_fold"] += 1
                continue
            if str(record.get("metal", "")) != match["metal"]:
                audit["excluded_metal_join_mismatch"] += 1
                continue
            audit["graph_fold_metal_pass"] += 1

            expected_hash = str(record.get("source_payload_sha256", ""))
            if (
                not expected_hash
                or expected_hash != single[mol_id]["official_xyz_sha256"]
                or not truth(record.get("source_payload_matches_frozen_zenodo"))
            ):
                audit["excluded_source_hash_lineage"] += 1
                continue
            audit["source_hash_lineage_pass"] += 1

            candidates.append(
                {
                    "mol_id": mol_id,
                    "fold": int(fold_text),
                    "metal": match["metal"],
                    "cn": cn,
                    "graph_sha256": graph,
                    "canonical_smiles": match["canonical_smiles"],
                    "source_payload_sha256": expected_hash,
                    "source_payload_bytes": int(record["source_payload_bytes"]),
                    "n_atoms": int(record["n_atoms"]),
                    "source_member": f"xyz/{mol_id}.xyz",
                    "ligand_payloads": list(record.get("ligand_smiles") or []),
                    "ligand_payload_provenance": provenance,
                    "connectivity_status": statuses,
                }
            )

    if n_rc2 != 31_121:
        raise ValueError(f"Expected 31,121 rc2 PBE rows, found {n_rc2:,}")
    if len(matched) != 31_121 or len(single) != 31_121:
        raise ValueError(
            f"Unexpected join ledger sizes: matched={len(matched):,}, single={len(single):,}"
        )
    audit["strict_eligible_rows"] = len(candidates)
    audit["strict_eligible_unique_graphs"] = len(
        {row["graph_sha256"] for row in candidates}
    )
    if audit["strict_eligible_unique_graphs"] != EXPECTED_ELIGIBLE_UNIQUE_GRAPHS:
        raise ValueError(
            "Strict CN4--6 eligible unique-graph denominator changed: "
            f"expected {EXPECTED_ELIGIBLE_UNIQUE_GRAPHS:,}, "
            f"found {audit['strict_eligible_unique_graphs']:,}"
        )
    return candidates, dict(audit)


def select_graph_disjoint_cohort(
    candidates: list[dict[str, Any]],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    by_stratum: dict[tuple[int, int], list[dict[str, Any]]] = defaultdict(list)
    for row in candidates:
        by_stratum[(int(row["fold"]), int(row["cn"]))].append(row)

    selected: list[dict[str, Any]] = []
    used_graphs: set[str] = set()
    stratum_audit: list[dict[str, Any]] = []
    # Rarer strata first prevents a duplicate graph from consuming scarce
    # capacity.  Ties use the fixed fold/CN order.
    strata = [(fold, cn) for fold in FOLDS for cn in MAIN_CNS]
    strata.sort(
        key=lambda key: (
            len({row["graph_sha256"] for row in by_stratum[key]})
            / CN_QUOTA_BY_FOLD[key[0]][key[1]],
            key[0],
            key[1],
        )
    )
    for fold, cn in strata:
        quota = CN_QUOTA_BY_FOLD[fold][cn]
        ranked = sorted(
            by_stratum[(fold, cn)],
            key=lambda row: deterministic_score(
                "within_stratum", fold, cn, row["graph_sha256"], row["mol_id"]
            ),
        )
        picked: list[dict[str, Any]] = []
        local_graphs: set[str] = set()
        for row in ranked:
            graph = row["graph_sha256"]
            if graph in used_graphs or graph in local_graphs:
                continue
            local_graphs.add(graph)
            used_graphs.add(graph)
            picked.append(dict(row))
            if len(picked) == quota:
                break
        if len(picked) != quota:
            raise ValueError(
                f"Insufficient graph-disjoint rows for fold={fold}, CN={cn}: "
                f"needed {quota}, selected {len(picked)}"
            )
        gate_quota = GATE_QUOTA_BY_FOLD[fold][cn]
        for stratum_rank, row in enumerate(picked, start=1):
            row["stratum_selection_rank"] = stratum_rank
            row["gate_member"] = stratum_rank <= gate_quota
        selected.extend(picked)
        stratum_audit.append(
            {
                "fold": fold,
                "cn": cn,
                "eligible_rows": len(ranked),
                "eligible_unique_graphs": len({row["graph_sha256"] for row in ranked}),
                "quota": quota,
                "selected": len(picked),
            }
        )

    if len(selected) != EXPECTED_N:
        raise AssertionError((len(selected), EXPECTED_N))
    if len({row["graph_sha256"] for row in selected}) != EXPECTED_N:
        raise AssertionError("Selected graph hashes are not unique")
    if len({row["mol_id"] for row in selected}) != EXPECTED_N:
        raise AssertionError("Selected molecule IDs are not unique")

    if sum(bool(row["gate_member"]) for row in selected) != EXPECTED_GATE_N:
        raise AssertionError("Gate allocation does not contain exactly 100 records")
    selected.sort(
        key=lambda row: (
            not bool(row["gate_member"]),
            deterministic_score("global_order", row["graph_sha256"], row["mol_id"]),
        )
    )
    for selection_order, row in enumerate(selected, start=1):
        row["selection_order"] = selection_order
        row["selection_score_sha256"] = deterministic_score(
            "global_order", row["graph_sha256"], row["mol_id"]
        )
    return selected, {"strata": stratum_audit}


def load_and_verify_payloads(selected: list[dict[str, Any]]) -> dict[str, bytes]:
    payloads: dict[str, bytes] = {}
    with zipfile.ZipFile(PBE_ZIP) as archive:
        names = set(archive.namelist())
        for row in selected:
            member = row["source_member"]
            if member not in names:
                raise FileNotFoundError(member)
            payload = archive.read(member)
            digest = sha256_bytes(payload)
            if digest != row["source_payload_sha256"]:
                raise ValueError(
                    f"ZIP payload hash mismatch for {row['mol_id']}: {digest}"
                )
            if len(payload) != int(row["source_payload_bytes"]):
                raise ValueError(f"ZIP payload byte-size mismatch for {row['mol_id']}")
            payloads[row["mol_id"]] = payload
    return payloads


def _install_minimal_smoke_patch(module: Any) -> None:
    """Install the diagnostic-only source-label-free duplicate-group patch."""

    def exact_source_label_free_relabeling(cc: Any) -> Any:
        groups: dict[tuple[Any, ...], list[str]] = defaultdict(list)
        for ligand in cc.ligands:
            key = (
                ligand.payload_provenance,
                ligand.smiles,
                tuple(sorted(str(value) for value in ligand.donor_elements)),
                ligand.dent,
                ligand.eta,
                ligand.charge,
                ligand.connectivity_status,
            )
            groups[key].append(ligand.lig_id)
        duplicates = [ids for ids in groups.values() if len(ids) > 1]
        if not duplicates:
            return cc
        candidate_count = math.prod(factorial(len(group)) for group in duplicates)
        if candidate_count > 100_000:
            raise ValueError(
                "diagnostic exact equivalent-ligand search exceeds 100,000 candidates"
            )
        best = cc
        best_string = cc.to_string()
        families = [tuple(permutations(group)) for group in duplicates]
        for choices in product(*families):
            mapping: dict[str, str] = {}
            for group, choice in zip(duplicates, choices):
                mapping.update(dict(zip(group, choice)))
            candidate = module._apply_permutation(cc, mapping)
            candidate_string = candidate.to_string()
            if candidate_string < best_string:
                best = candidate
                best_string = candidate_string
        return best

    module._apply_equivalent_relabeling = exact_source_label_free_relabeling


def initialize_worker(release_root: str, patch_mode: str) -> None:
    release = Path(release_root).resolve()
    release_text = str(release)
    if release_text not in sys.path:
        sys.path.insert(0, release_text)
    from rdkit import RDLogger

    RDLogger.DisableLog("rdApp.*")
    import coordrep
    import coordrep.canonical.canonicalize as canonical_module
    from coordrep import encode_molecule
    from coordrep.io.tmqm_reader import Atom, RawMolecule

    package_path = Path(coordrep.__file__).resolve()
    if release not in package_path.parents:
        raise ImportError(f"CoordRep import escaped release root: {package_path}")
    if patch_mode == "minimal_smoke":
        _install_minimal_smoke_patch(canonical_module)
    _WORKER.update(
        {
            "Atom": Atom,
            "RawMolecule": RawMolecule,
            "encode_molecule": encode_molecule,
            "release_root": release,
            "patch_mode": patch_mode,
        }
    )


def parse_xyz_payload(payload: bytes, mol_id: str) -> Any:
    text = payload.decode("utf-8")
    lines = text.splitlines()
    n_atoms = int(lines[0].strip())
    if len(lines) < n_atoms + 2:
        raise ValueError(f"Truncated XYZ payload for {mol_id}")
    atoms = []
    for index, line in enumerate(lines[2 : n_atoms + 2]):
        fields = line.split()
        atoms.append(
            _WORKER["Atom"](
                index=index,
                element=fields[0],
                x=float(fields[1]),
                y=float(fields[2]),
                z=float(fields[3]),
            )
        )
    return _WORKER["RawMolecule"](
        mol_id=mol_id,
        atoms=atoms,
        bond_orders=None,
        properties={},
    )


def proper_rotation_matrix(rng: np.random.Generator) -> np.ndarray:
    """Shoemake Haar-uniform proper rotation, deterministic from NumPy RNG."""

    u1, u2, u3 = rng.random(3)
    x = math.sqrt(1.0 - u1) * math.sin(2.0 * math.pi * u2)
    y = math.sqrt(1.0 - u1) * math.cos(2.0 * math.pi * u2)
    z = math.sqrt(u1) * math.sin(2.0 * math.pi * u3)
    w = math.sqrt(u1) * math.cos(2.0 * math.pi * u3)
    return np.array(
        [
            [1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
            [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
            [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)],
        ],
        dtype=np.float64,
    )


def make_variant(
    molecule: Any,
    selection_order: int,
    arm: str,
    variant_index: int,
) -> Any:
    rng = np.random.default_rng(
        np.random.SeedSequence(
            [TRANSFORM_SEED, int(selection_order), ARM_CODE[arm], int(variant_index)]
        )
    )
    coords = molecule.get_coords().astype(np.float64, copy=True)
    if arm in {"rigid", "combined"}:
        coords = coords @ proper_rotation_matrix(rng).T
        coords = coords + rng.uniform(-25.0, 25.0, size=3)
    order = np.arange(len(molecule.atoms), dtype=np.int64)
    if arm in {"permutation", "combined"}:
        order = rng.permutation(order)

    atoms = []
    for new_index, old_index_raw in enumerate(order):
        old_index = int(old_index_raw)
        source_atom = molecule.atoms[old_index]
        xyz = coords[old_index]
        atoms.append(
            _WORKER["Atom"](
                index=new_index,
                element=source_atom.element,
                x=float(xyz[0]),
                y=float(xyz[1]),
                z=float(xyz[2]),
            )
        )

    bond_orders = None
    if molecule.bond_orders is not None:
        matrix = np.asarray(molecule.bond_orders)
        bond_orders = matrix[np.ix_(order, order)].copy()
    return _WORKER["RawMolecule"](
        mol_id=molecule.mol_id,
        atoms=atoms,
        bond_orders=bond_orders,
        properties=dict(molecule.properties or {}),
    )


def digest_array(label: str, elements: list[str], array: np.ndarray) -> str:
    digest = hashlib.sha256()
    digest.update(label.encode("ascii"))
    digest.update("\0".join(elements).encode("ascii"))
    rounded = np.ascontiguousarray(
        np.round(array, decimals=ARRAY_SERIALIZATION_DECIMALS), dtype="<f8"
    )
    digest.update(rounded.tobytes(order="C"))
    return digest.hexdigest()


def encode_representation_keys(
    molecule: Any,
    *,
    return_strings: bool = False,
) -> tuple[dict[str, str], dict[str, str], dict[str, Any]]:
    record = _WORKER["encode_molecule"](molecule)
    if record.metal.element == "?" or record.graph is None:
        raise ValueError(f"Encoder placeholder for {molecule.mol_id}")
    metal_index = int(record.graph.metal_idx)
    donors = [int(value) for value in record.graph.donor_indices]
    # The diagnostic controls deliberately retain the source atom order of the
    # first sphere.  The tested CoordRep key below uses the same encoded object
    # and then invokes canonicalize().
    indices = [metal_index, *sorted(donors)]
    elements = [molecule.atoms[index].element for index in indices]
    coords = np.asarray([molecule.atoms[index].coords for index in indices], dtype=float)
    distances = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=2)

    raw_hash = digest_array("raw-first-sphere-cartesian-v1", elements, coords)
    distance_hash = digest_array("ordered-first-sphere-distance-v1", elements, distances)

    ligand_records = []
    for ligand in record.ligands:
        ligand_records.append(
            {
                "payload_provenance": ligand.payload_provenance,
                "smiles": ligand.smiles,
                "donor_elements": sorted(str(v) for v in ligand.donor_elements),
                "dent": int(ligand.dent),
                "eta": ligand.eta,
                "charge": ligand.charge,
                "connectivity_status": ligand.connectivity_status,
            }
        )
    ligand_records.sort(
        key=lambda value: json.dumps(value, sort_keys=True, separators=(",", ":"))
    )
    connectivity_object = {
        "metal": record.metal.element,
        "cn": len(donors),
        "ligands": ligand_records,
    }
    connectivity = json.dumps(
        connectivity_object, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    )
    uncanonicalized = record.to_string()
    canonicalized = record.canonicalize().to_string()

    strings = {
        "connectivity_only": connectivity,
        "uncanonicalized_coordrep": uncanonicalized,
        "patched_coordrep": canonicalized,
    }
    hashes = {
        "raw_cartesian": raw_hash,
        "ordered_distance_matrix": distance_hash,
        **{key: sha256_bytes(value.encode("utf-8")) for key, value in strings.items()},
    }
    meta = {
        "metal": record.metal.element,
        "cn": len(donors),
        "n_ligands": len(record.ligands),
        "issues": [
            {"code": issue.code, "severity": issue.severity}
            for issue in record.issues
        ],
    }
    return hashes, (strings if return_strings else {}), meta


def evaluate_molecule(
    molecule: Any,
    metadata: dict[str, Any],
    *,
    k: int,
    run_stage: str,
    capture_coordrep_mismatches: bool = True,
) -> dict[str, Any]:
    reference, reference_strings, encoded_meta = encode_representation_keys(
        molecule, return_strings=True
    )
    if encoded_meta["metal"] != metadata["metal"] or encoded_meta["cn"] != int(metadata["cn"]):
        raise ValueError(
            f"Baseline re-encoding changed metal/CN for {metadata['mol_id']}: "
            f"{encoded_meta['metal']}/CN{encoded_meta['cn']}"
        )

    state: dict[tuple[str, str], dict[str, Any]] = {}
    for arm in ARMS:
        for representation in REPRESENTATIONS:
            state[(arm, representation)] = {
                "matched": 0,
                "hashes": set(),
            }
    mismatches: list[dict[str, Any]] = []

    for arm in ARMS:
        for variant_index in range(k):
            variant = make_variant(
                molecule,
                int(metadata["selection_order"]),
                arm,
                variant_index,
            )
            try:
                observed, observed_strings, _ = encode_representation_keys(
                    variant, return_strings=capture_coordrep_mismatches
                )
            except Exception as error:
                mismatches.append(
                    {
                        "run_stage": run_stage,
                        "selection_order": metadata["selection_order"],
                        "mol_id": metadata["mol_id"],
                        "fold": metadata["fold"],
                        "metal": metadata["metal"],
                        "cn": metadata["cn"],
                        "arm": arm,
                        "variant_index": variant_index,
                        "representation": "patched_coordrep",
                        "reference_sha256": reference["patched_coordrep"],
                        "variant_sha256": "",
                        "reference_value": reference_strings["patched_coordrep"],
                        "variant_value": "",
                        "failure_stage": "encode_variant",
                        "error": f"{type(error).__name__}: {error}",
                    }
                )
                continue
            for representation in REPRESENTATIONS:
                value = observed[representation]
                slot = state[(arm, representation)]
                slot["hashes"].add(value)
                if value == reference[representation]:
                    slot["matched"] += 1
                elif capture_coordrep_mismatches and representation == "patched_coordrep":
                    mismatches.append(
                        {
                            "run_stage": run_stage,
                            "selection_order": metadata["selection_order"],
                            "mol_id": metadata["mol_id"],
                            "fold": metadata["fold"],
                            "metal": metadata["metal"],
                            "cn": metadata["cn"],
                            "arm": arm,
                            "variant_index": variant_index,
                            "representation": representation,
                            "reference_sha256": reference[representation],
                            "variant_sha256": value,
                            "reference_value": reference_strings[representation],
                            "variant_value": observed_strings[representation],
                            "failure_stage": "invariance_mismatch",
                            "error": "",
                        }
                    )

    summaries: list[dict[str, Any]] = []
    for arm in ARMS:
        for representation in REPRESENTATIONS:
            slot = state[(arm, representation)]
            matched = int(slot["matched"])
            hashes = set(slot["hashes"])
            summaries.append(
                {
                    "run_stage": run_stage,
                    "selection_order": metadata["selection_order"],
                    "mol_id": metadata["mol_id"],
                    "fold": metadata["fold"],
                    "metal": metadata["metal"],
                    "cn": metadata["cn"],
                    "graph_sha256": metadata.get("graph_sha256", ""),
                    "arm": arm,
                    "k": k,
                    "representation": representation,
                    "reference_sha256": reference[representation],
                    "matched_variants": matched,
                    "mismatched_variants": k - matched,
                    "match_fraction": matched / k,
                    "unique_variant_hashes": len(hashes),
                    "unique_hashes_with_reference": len(hashes | {reference[representation]}),
                    "all_variants_match_reference": matched == k,
                }
            )
    return {
        "summaries": summaries,
        "mismatches": mismatches,
        "reference_hashes": reference,
        "reference_strings": reference_strings,
        "encoded_meta": encoded_meta,
    }


def worker_task(task: tuple[dict[str, Any], bytes, int, str]) -> dict[str, Any]:
    metadata, payload, k, run_stage = task
    molecule = parse_xyz_payload(payload, metadata["mol_id"])
    try:
        result = evaluate_molecule(
            molecule,
            metadata,
            k=k,
            run_stage=run_stage,
        )
        result["fatal"] = None
        return result
    except Exception as error:
        return {
            "summaries": [],
            "mismatches": [
                {
                    "run_stage": run_stage,
                    "selection_order": metadata["selection_order"],
                    "mol_id": metadata["mol_id"],
                    "fold": metadata["fold"],
                    "metal": metadata["metal"],
                    "cn": metadata["cn"],
                    "arm": "",
                    "variant_index": "",
                    "representation": "patched_coordrep",
                    "reference_sha256": "",
                    "variant_sha256": "",
                    "failure_stage": "baseline_or_worker",
                    "error": f"{type(error).__name__}: {error}",
                }
            ],
            "reference_hashes": {},
            "reference_strings": {},
            "encoded_meta": {},
            "fatal": f"{type(error).__name__}: {error}",
        }


def run_parallel(
    rows: list[dict[str, Any]],
    payloads: dict[str, bytes],
    *,
    k: int,
    run_stage: str,
    release_root: Path,
    patch_mode: str,
    workers: int,
    log: Any,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[str]]:
    summaries: list[dict[str, Any]] = []
    mismatches: list[dict[str, Any]] = []
    fatals: list[str] = []
    tasks = [(row, payloads[row["mol_id"]], k, run_stage) for row in rows]
    started = time.perf_counter()
    with ProcessPoolExecutor(
        max_workers=workers,
        initializer=initialize_worker,
        initargs=(str(release_root.resolve()), patch_mode),
    ) as executor:
        futures = [executor.submit(worker_task, task) for task in tasks]
        for completed, future in enumerate(as_completed(futures), start=1):
            result = future.result()
            summaries.extend(result["summaries"])
            mismatches.extend(result["mismatches"])
            if result.get("fatal"):
                fatals.append(result["fatal"])
            if completed % max(1, min(25, len(tasks))) == 0 or completed == len(tasks):
                elapsed = time.perf_counter() - started
                log(
                    f"{run_stage}: {completed:,}/{len(tasks):,} structures; "
                    f"{completed/elapsed:.2f} structures/s"
                )
    summaries.sort(
        key=lambda row: (
            int(row["selection_order"]),
            ARMS.index(row["arm"]),
            REPRESENTATIONS.index(row["representation"]),
        )
    )
    mismatches.sort(
        key=lambda row: (
            int(row["selection_order"]),
            str(row["arm"]),
            int(row["variant_index"]) if str(row["variant_index"]) else -1,
        )
    )
    return summaries, mismatches, fatals


def aggregate_summaries(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    groups: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        groups[(row["arm"], row["representation"])].append(row)
    output = []
    for arm in ARMS:
        for representation in REPRESENTATIONS:
            values = groups.get((arm, representation), [])
            k_total = sum(int(row["k"]) for row in values)
            matches = sum(int(row["matched_variants"]) for row in values)
            output.append(
                {
                    "arm": arm,
                    "representation": representation,
                    "n_structures": len(values),
                    "variant_encodings": k_total,
                    "matched_encodings": matches,
                    "mismatched_encodings": k_total - matches,
                    "invariant_match_fraction": matches / k_total if k_total else None,
                    "fully_collapsed_structures": sum(
                        truth(row["all_variants_match_reference"]) for row in values
                    ),
                    "structure_collapse_fraction": (
                        sum(truth(row["all_variants_match_reference"]) for row in values)
                        / len(values)
                        if values
                        else None
                    ),
                }
            )
    return output


def parse_mol2(path: Path, mol_id: str) -> Any:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    atom_start = lines.index("@<TRIPOS>ATOM") + 1
    bond_start = lines.index("@<TRIPOS>BOND") + 1
    next_sections = [
        index
        for index in range(bond_start, len(lines))
        if lines[index].startswith("@<TRIPOS>")
    ]
    bond_end = next_sections[0] if next_sections else len(lines)
    atom_lines = lines[atom_start : bond_start - 1]
    atoms = []
    id_to_index: dict[int, int] = {}
    for new_index, line in enumerate(atom_lines):
        fields = line.split()
        atom_id = int(fields[0])
        atom_type = fields[5].split(".")[0]
        element = "".join(character for character in atom_type if character.isalpha())
        if len(element) >= 2 and element[:2].title() in {"Cl", "Br", "Ru"}:
            element = element[:2].title()
        else:
            element = element[:1].upper()
        id_to_index[atom_id] = new_index
        atoms.append(
            _WORKER["Atom"](
                index=new_index,
                element=element,
                x=float(fields[2]),
                y=float(fields[3]),
                z=float(fields[4]),
            )
        )
    bond_orders = np.zeros((len(atoms), len(atoms)), dtype=np.float32)
    order_map = {"1": 1.0, "2": 2.0, "3": 3.0, "ar": 1.5, "am": 1.0}
    for line in lines[bond_start:bond_end]:
        fields = line.split()
        if len(fields) < 4:
            continue
        left = id_to_index[int(fields[1])]
        right = id_to_index[int(fields[2])]
        value = order_map.get(fields[3].lower(), 1.0)
        bond_orders[left, right] = value
        bond_orders[right, left] = value
    return _WORKER["RawMolecule"](
        mol_id=mol_id,
        atoms=atoms,
        bond_orders=bond_orders,
        properties={},
    )


def run_fac_mer_specificity(
    release_root: Path,
    patch_mode: str,
    k: int,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    initialize_worker(str(release_root.resolve()), patch_mode)
    molecules = {
        "fac": parse_mol2(FAC_MOL2, "EBAGAR"),
        "mer": parse_mol2(MER_MOL2, "EBAGEV"),
    }
    results: dict[str, dict[str, Any]] = {}
    for index, (label, molecule) in enumerate(molecules.items(), start=1):
        results[label] = evaluate_molecule(
            molecule,
            {
                "selection_order": 20_000 + index,
                "mol_id": molecule.mol_id,
                "fold": "specificity",
                "metal": "Ru",
                "cn": 6,
                "graph_sha256": "same-mol2-topology-audited-separately",
            },
            k=k,
            run_stage="fac_mer_specificity",
            capture_coordrep_mismatches=True,
        )
    rows = []
    for representation in REPRESENTATIONS:
        row: dict[str, Any] = {
            "representation": representation,
            "fac_reference_sha256": results["fac"]["reference_hashes"][representation],
            "mer_reference_sha256": results["mer"]["reference_hashes"][representation],
            "cross_isomer_distinct": (
                results["fac"]["reference_hashes"][representation]
                != results["mer"]["reference_hashes"][representation]
            ),
        }
        for label in ("fac", "mer"):
            for arm in ARMS:
                match = next(
                    value
                    for value in results[label]["summaries"]
                    if value["representation"] == representation and value["arm"] == arm
                )
                row[f"{label}_{arm}_matched"] = match["matched_variants"]
                row[f"{label}_{arm}_k"] = match["k"]
                row[f"{label}_{arm}_collapse"] = match["all_variants_match_reference"]
        rows.append(row)

    audit = {
        "fac_id": "EBAGAR",
        "mer_id": "EBAGEV",
        "same_connectivity_only_key": (
            results["fac"]["reference_hashes"]["connectivity_only"]
            == results["mer"]["reference_hashes"]["connectivity_only"]
        ),
        "distinct_full_coordrep_keys": (
            results["fac"]["reference_hashes"]["patched_coordrep"]
            != results["mer"]["reference_hashes"]["patched_coordrep"]
        ),
        "fac_coordrep": results["fac"]["reference_strings"]["patched_coordrep"],
        "mer_coordrep": results["mer"]["reference_strings"]["patched_coordrep"],
        "coordrep_mismatch_rows": len(results["fac"]["mismatches"])
        + len(results["mer"]["mismatches"]),
    }
    return rows, audit


def prepare_output(directory: Path, force: bool) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    known = [
        "cohort_manifest.csv",
        "variant_summary.csv",
        "mismatch_variants.csv",
        "fac_mer_specificity.csv",
        "summary.json",
        "PROVENANCE.json",
        "README.md",
        "SHA256SUMS.csv",
        "run_log.txt",
    ]
    existing = [directory / name for name in known if (directory / name).exists()]
    if existing and not force:
        raise FileExistsError(
            "Refusing to overwrite existing outputs; pass --force: "
            + ", ".join(path.name for path in existing)
        )
    if force:
        resolved = directory.resolve()
        if OUT.resolve() not in (resolved, *resolved.parents):
            raise ValueError(f"Refusing cleanup outside task directory: {resolved}")
        for path in existing:
            path.unlink()


def main() -> int:
    args = parse_args()
    if args.workers <= 0 or args.gate_n <= 0 or args.gate_k <= 0 or args.full_k <= 0:
        raise ValueError("workers and sample/trial counts must be positive")
    if args.gate_n > EXPECTED_N:
        raise ValueError("gate_n exceeds the locked cohort")
    if args.patch_mode == "minimal_smoke" and args.run_full:
        raise ValueError(
            "minimal_smoke is diagnostic only; publication-scale execution requires "
            "a source tree with ligand-local canonical attachment/orbit keys"
        )
    output_dir = args.output_dir
    if output_dir is None:
        output_dir = OUT if args.patch_mode == "source" else OUT / "minimal_fix_smoke"
    output_dir = output_dir.resolve()
    prepare_output(output_dir, args.force)

    log_path = output_dir / "run_log.txt"

    def log(message: str) -> None:
        line = f"[{utc_now()}] {message}"
        print(line, flush=True)
        with log_path.open("a", encoding="utf-8", newline="\n") as handle:
            handle.write(line + "\n")

    started_utc = utc_now()
    started = time.perf_counter()
    release_root = args.release_root.resolve()
    canonicalizer = release_root / "coordrep" / "canonical" / "canonicalize.py"
    encoder = release_root / "coordrep" / "encode.py"
    if not canonicalizer.exists() or not encoder.exists():
        raise FileNotFoundError(f"Invalid CoordRep release root: {release_root}")
    if release_root == PUBLIC_RELEASE.resolve():
        if sha256_file(canonicalizer) != EXPECTED_PUBLIC_CANONICALIZER_SHA256:
            raise ValueError("Public canonicalizer source hash mismatch")
        if sha256_file(encoder) != EXPECTED_PUBLIC_ENCODER_SHA256:
            raise ValueError("Public encoder source hash mismatch")
        if git_value(PUBLIC_REPO, "rev-parse", "HEAD") != EXPECTED_PUBLIC_COMMIT:
            raise ValueError("Public checkout commit mismatch")
    task_local_candidate = release_root == CANDIDATE_FIX.resolve()
    if task_local_candidate:
        for relative_path, expected_sha256 in EXPECTED_CANDIDATE_SOURCE_SHA256.items():
            observed_sha256 = sha256_file(release_root / relative_path)
            if observed_sha256 != expected_sha256:
                raise ValueError(
                    f"Frozen task-local candidate source hash mismatch: {relative_path}; "
                    f"expected {expected_sha256}, found {observed_sha256}"
                )
    if sha256_file(PBE_ZIP) != EXPECTED_PBE_ZIP_SHA256:
        raise ValueError("Official tmQMg/PBE XYZ ZIP hash mismatch")

    log("loading and auditing the four-source tmQMg/PBE join")
    candidates, join_audit = load_joined_candidates()
    selected, selection_audit = select_graph_disjoint_cohort(candidates)
    payloads = load_and_verify_payloads(selected)
    log(
        f"selected {len(selected):,} unique graphs; verified {len(payloads):,} ZIP payload hashes"
    )

    manifest_fields = [
        "selection_order",
        "gate_member",
        "stratum_selection_rank",
        "mol_id",
        "fold",
        "metal",
        "cn",
        "graph_sha256",
        "canonical_smiles",
        "source_member",
        "source_payload_bytes",
        "source_payload_sha256",
        "n_atoms",
        "selection_score_sha256",
        "ligand_payloads_json",
        "ligand_payload_provenance_json",
        "connectivity_status_json",
    ]
    manifest_rows = []
    for row in selected:
        manifest_rows.append(
            {
                **row,
                "ligand_payloads_json": json.dumps(row["ligand_payloads"], separators=(",", ":")),
                "ligand_payload_provenance_json": json.dumps(
                    row["ligand_payload_provenance"], separators=(",", ":")
                ),
                "connectivity_status_json": json.dumps(
                    row["connectivity_status"], separators=(",", ":")
                ),
            }
        )
    write_csv(output_dir / "cohort_manifest.csv", manifest_fields, manifest_rows)

    gate_rows = selected[: args.gate_n]
    if args.gate_n == EXPECTED_GATE_N and not all(
        bool(row["gate_member"]) for row in gate_rows
    ):
        raise AssertionError("Locked N=100 gate is not the fold/CN-balanced gate")
    log(
        f"running gate: N={len(gate_rows):,}, K={args.gate_k} per arm, "
        f"patch_mode={args.patch_mode}, workers={args.workers}"
    )
    gate_summary_rows, gate_mismatches, gate_fatals = run_parallel(
        gate_rows,
        payloads,
        k=args.gate_k,
        run_stage="gate",
        release_root=release_root,
        patch_mode=args.patch_mode,
        workers=args.workers,
        log=log,
    )
    gate_aggregate = aggregate_summaries(gate_summary_rows)
    gate_coordrep = [
        row for row in gate_aggregate if row["representation"] == "patched_coordrep"
    ]
    gate_connectivity = [
        row for row in gate_aggregate if row["representation"] == "connectivity_only"
    ]
    gate_pass = (
        not gate_fatals
        and not gate_mismatches
        and all(row["mismatched_encodings"] == 0 for row in gate_coordrep)
        and all(row["mismatched_encodings"] == 0 for row in gate_connectivity)
    )
    log(f"gate status: {'PASS' if gate_pass else 'FAIL'}")

    run_stage = "gate"
    final_rows = gate_summary_rows
    final_mismatches = gate_mismatches
    full_aggregate: list[dict[str, Any]] = []
    full_fatals: list[str] = []
    if args.run_full and gate_pass:
        run_stage = "full"
        log(
            f"running locked full stress: N={len(selected):,}, K={args.full_k} per arm, "
            f"total CoordRep variant encodings={len(selected)*args.full_k*len(ARMS):,}"
        )
        final_rows, final_mismatches, full_fatals = run_parallel(
            selected,
            payloads,
            k=args.full_k,
            run_stage="full",
            release_root=release_root,
            patch_mode=args.patch_mode,
            workers=args.workers,
            log=log,
        )
        full_aggregate = aggregate_summaries(final_rows)
    elif args.run_full and not gate_pass:
        log("full run skipped because the gate failed")

    write_csv(output_dir / "variant_summary.csv", SUMMARY_FIELDS, final_rows)
    write_csv(output_dir / "mismatch_variants.csv", MISMATCH_FIELDS, final_mismatches)

    log("running EBAGAR/EBAGEV fac/mer specificity control")
    facmer_rows, facmer_audit = run_fac_mer_specificity(
        release_root, args.patch_mode, k=min(args.full_k, 100)
    )
    facmer_fields = list(facmer_rows[0].keys())
    write_csv(output_dir / "fac_mer_specificity.csv", facmer_fields, facmer_rows)

    full_coordrep = [
        row for row in full_aggregate if row["representation"] == "patched_coordrep"
    ]
    publication_pass = bool(full_aggregate) and (
        not full_fatals
        and not final_mismatches
        and all(row["mismatched_encodings"] == 0 for row in full_coordrep)
        and facmer_audit["same_connectivity_only_key"]
        and facmer_audit["distinct_full_coordrep_keys"]
        and facmer_audit["coordrep_mismatch_rows"] == 0
    )
    if args.patch_mode == "minimal_smoke":
        status = "SMOKE_ONLY_NOT_FOR_PUBLICATION"
    elif task_local_candidate and args.run_full and publication_pass:
        status = "TASK_LOCAL_CANDIDATE_STRESS_PASS"
    elif task_local_candidate and args.run_full:
        status = "TASK_LOCAL_CANDIDATE_STRESS_FAIL"
    elif task_local_candidate:
        status = "TASK_LOCAL_CANDIDATE_GATE_PASS" if gate_pass else "TASK_LOCAL_CANDIDATE_GATE_FAIL"
    elif args.run_full and publication_pass:
        status = "PUBLICATION_STRESS_PASS"
    elif args.run_full:
        status = "PUBLICATION_STRESS_FAIL"
    else:
        status = "SOURCE_GATE_PASS" if gate_pass else "SOURCE_GATE_FAIL"

    summary = {
        "status": status,
        "claim_boundary": (
            "This fixed tmQMg/PBE graph-disjoint benchmark is not an April-2025-CSD prevalence sample. "
            "minimal_smoke results are diagnostic only. canonical_fix_20260822 is a frozen task-local "
            "candidate implementation, not an already integrated or public release."
        ),
        "patch_mode": args.patch_mode,
        "release_root": str(release_root),
        "task_local_candidate": task_local_candidate,
        "selection_seed": SELECTION_SEED,
        "transform_seed": TRANSFORM_SEED,
        "cohort_n": len(selected),
        "unique_graphs": len({row["graph_sha256"] for row in selected}),
        "cn_scope": list(MAIN_CNS),
        "quota_by_fold": CN_QUOTA_BY_FOLD,
        "gate_quota_by_fold": GATE_QUOTA_BY_FOLD,
        "eligible_unique_graph_denominator": EXPECTED_ELIGIBLE_UNIQUE_GRAPHS,
        "join_audit": join_audit,
        "selection_audit": selection_audit,
        "gate": {
            "n": len(gate_rows),
            "k_per_arm": args.gate_k,
            "passed": gate_pass,
            "fatal_count": len(gate_fatals),
            "coordrep_mismatch_rows": len(gate_mismatches),
            "aggregate": gate_aggregate,
            "end_to_end_raw_molecule_permutation_regression": next(
                (
                    row
                    for row in gate_coordrep
                    if row["arm"] == "permutation"
                ),
                None,
            ),
        },
        "full": {
            "requested": bool(args.run_full),
            "executed": bool(full_aggregate),
            "n": len(selected) if full_aggregate else 0,
            "k_per_arm": args.full_k if full_aggregate else 0,
            "fatal_count": len(full_fatals),
            "coordrep_mismatch_rows": len(final_mismatches) if full_aggregate else None,
            "publication_pass": publication_pass,
            "aggregate": full_aggregate,
        },
        "fac_mer_specificity": facmer_audit,
        "representations": list(REPRESENTATIONS),
        "publication_figure_representations": list(PUBLICATION_REPRESENTATIONS),
        "arms": list(ARMS),
    }
    (output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )

    dependency_versions: dict[str, str] = {"numpy": np.__version__}
    for module_name in ("scipy", "networkx", "rdkit"):
        try:
            module = __import__(module_name)
            dependency_versions[module_name] = str(getattr(module, "__version__", "unknown"))
        except Exception as error:
            dependency_versions[module_name] = f"unavailable: {error}"

    provenance = {
        "created_utc": utc_now(),
        "started_utc": started_utc,
        "elapsed_seconds": time.perf_counter() - started,
        "status": status,
        "script": str(SCRIPT),
        "script_sha256": sha256_file(SCRIPT),
        "python": sys.version,
        "platform": platform.platform(),
        "dependencies": dependency_versions,
        "workers": args.workers,
        "release": {
            "root": str(release_root),
            "patch_mode": args.patch_mode,
            "task_local_candidate": task_local_candidate,
            "candidate_source_sha256_lock": (
                EXPECTED_CANDIDATE_SOURCE_SHA256 if task_local_candidate else None
            ),
            "canonicalizer_sha256": sha256_file(canonicalizer),
            "encoder_sha256": sha256_file(encoder),
            "public_repo_commit": git_value(PUBLIC_REPO, "rev-parse", "HEAD"),
            "public_repo_branch": git_value(PUBLIC_REPO, "branch", "--show-current"),
        },
        "minimal_smoke_patch": {
            "applied": args.patch_mode == "minimal_smoke",
            "description": (
                "Diagnostic runtime replacement of duplicate-ligand equivalence grouping: "
                "raw/global attach_atoms are omitted; exact Cartesian-product minimization is retained."
            ),
            "publication_safe": False,
            "known_limit": (
                "Removing raw attachment indices without a ligand-local canonical attachment/orbit key "
                "can merge the same payload bound through non-equivalent sites."
            ),
        },
        "sources": {
            "rc2_records": file_provenance(RC2_RECORDS),
            "matched_records": file_provenance(MATCHED_RECORDS),
            "single_metal_audit": file_provenance(SINGLE_METAL_AUDIT),
            "official_pbe_xyz_zip": file_provenance(PBE_ZIP),
            "official_pbe_xyz_manifest": file_provenance(PBE_ZIP_MANIFEST),
            "fac_mol2": file_provenance(FAC_MOL2),
            "mer_mol2": file_provenance(MER_MOL2),
        },
        "selection": {
            "seed": SELECTION_SEED,
            "folds": list(FOLDS),
            "cn_scope": list(MAIN_CNS),
            "quota_by_fold": CN_QUOTA_BY_FOLD,
            "overall_cn_quota": {4: 334, 5: 333, 6: 333},
            "gate_quota_by_fold": GATE_QUOTA_BY_FOLD,
            "eligible_unique_graph_denominator": EXPECTED_ELIGIBLE_UNIQUE_GRAPHS,
            "graph_disjoint": True,
            "formula_fallback_excluded": True,
            "single_transition_metal_xyz_required": True,
            "zip_payload_sha256_verified_per_record": True,
        },
        "transformations": {
            "seed": TRANSFORM_SEED,
            "rigid": "Shoemake Haar-uniform SO(3) rotation plus uniform translation in [-25,25]^3 Angstrom",
            "permutation": "uniform random permutation of all atom rows; bond matrix permuted when present",
            "combined": "independent rigid transform followed by all-atom permutation",
            "coordinate_rounding_decimals_for_raw_controls": ARRAY_SERIALIZATION_DECIMALS,
            "coordrep_rounding": "release default",
        },
    }
    (output_dir / "PROVENANCE.json").write_text(
        json.dumps(provenance, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )

    aggregate_to_report = full_aggregate if full_aggregate else gate_aggregate
    readme = build_readme(
        status=status,
        patch_mode=args.patch_mode,
        gate_pass=gate_pass,
        full_executed=bool(full_aggregate),
        publication_pass=publication_pass,
        aggregate=aggregate_to_report,
        facmer=facmer_audit,
        args=args,
        task_local_candidate=task_local_candidate,
    )
    (output_dir / "README.md").write_text(readme, encoding="utf-8")
    log(f"completed with status {status}; outputs: {output_dir}")
    write_checksums(output_dir)
    if args.patch_mode == "minimal_smoke":
        return 0 if gate_pass else 1
    if args.run_full:
        return 0 if publication_pass else 1
    return 0 if gate_pass else 1


def file_provenance(path: Path) -> dict[str, Any]:
    return {
        "path": str(path.resolve()),
        "bytes": path.stat().st_size,
        "sha256": sha256_file(path),
    }


def build_readme(
    *,
    status: str,
    patch_mode: str,
    gate_pass: bool,
    full_executed: bool,
    publication_pass: bool,
    aggregate: list[dict[str, Any]],
    facmer: dict[str, Any],
    args: argparse.Namespace,
    task_local_candidate: bool,
) -> str:
    table_lines = [
        "| Arm | Representation | Matched encodings | Match fraction | Fully collapsed structures |",
        "|---|---|---:|---:|---:|",
    ]
    for row in aggregate:
        table_lines.append(
            f"| {row['arm']} | {row['representation']} | "
            f"{row['matched_encodings']:,}/{row['variant_encodings']:,} | "
            f"{row['invariant_match_fraction']:.6f} | "
            f"{row['fully_collapsed_structures']:,}/{row['n_structures']:,} |"
        )
    warning = (
        "**Diagnostic only.** The minimal patch removes raw attachment indices but does not yet "
        "replace them with ligand-local canonical attachment/orbit keys. These numbers must not be "
        "used in the manuscript or response letter."
        if patch_mode == "minimal_smoke"
        else (
            (
                "**Implementation boundary.** The tested `canonical_fix_20260822` source is a "
                "frozen task-local candidate. It is not yet integrated into the public release."
            )
            if task_local_candidate
            else (
                "The full result passed the locked publication gate."
                if publication_pass
                else "The source-mode result has not passed the locked full publication gate."
            )
        )
    )
    return f"""# CoordRep canonical invariance stress test

Status: **{status}**

{warning}

This package uses a fixed **tmQMg/PBE graph-disjoint CN4--6 benchmark cohort**.
It is not an April 2025 CSD prevalence sample. The selection joins the corrected
rc2 PBE ledger, graph/fold assignments, an independent single-transition-metal
XYZ audit, and the official individual tmQMg/PBE XYZ ZIP. Formula fallbacks
are excluded. The strict eligible denominator is 6,427 unique molecular graphs.
The locked cohort contains CN4/CN5/CN6 = 334/333/333 records, with 200 records
from each fold; `cohort_manifest.csv` records the exact fold-by-CN allocation.

## Execution

- Patch mode: `{patch_mode}`
- Gate: N={args.gate_n}, K={args.gate_k} per arm; **{'PASS' if gate_pass else 'FAIL'}**
- Full N=1,000, K={args.full_k} per arm executed: **{str(full_executed).lower()}**
- Publication pass: **{str(publication_pass).lower()}**
- Arms: rigid motion only; atom permutation only; combined.
- Full CoordRep mismatch details are written to `mismatch_variants.csv`; the
  file retains its header even when no mismatches occur.

## Results

{chr(10).join(table_lines)}

## fac/mer specificity control

EBAGAR and EBAGEV have the same connectivity-only key:
**{str(facmer['same_connectivity_only_key']).lower()}**. Their full CoordRep
keys are distinct: **{str(facmer['distinct_full_coordrep_keys']).lower()}**.
This control distinguishes nuisance invariance within each isomer from the
required preservation of a genuine fac/mer relation difference.

## Files

- `run_canonical_stress.py`: executable harness.
- `cohort_manifest.csv`: locked 1,000-record selection and source hashes.
- `variant_summary.csv`: structure-level results for every arm and representation.
- `mismatch_variants.csv`: full-CoordRep mismatches/exceptions (header retained if empty).
- `fac_mer_specificity.csv`: EBAGAR/EBAGEV control.
- `summary.json`: aggregate results and gate decisions.
- `PROVENANCE.json`: source, environment, transformation, and patch lineage.
- `SHA256SUMS.csv`: checksums of this evidence directory.
"""


def write_checksums(directory: Path) -> None:
    rows = []
    for path in sorted(directory.iterdir(), key=lambda value: value.name):
        if not path.is_file() or path.name == "SHA256SUMS.csv":
            continue
        rows.append(
            {
                "path": path.name,
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        )
    write_csv(directory / "SHA256SUMS.csv", ["path", "bytes", "sha256"], rows)


if __name__ == "__main__":
    raise SystemExit(main())
