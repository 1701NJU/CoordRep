#!/usr/bin/env python
"""Scan and challenge real legal canonicalization orbits in the 6,427 pool."""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import math
import os
import platform
import re
import sys
import time
import zipfile
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from copy import deepcopy
from datetime import datetime, timezone
from math import factorial
from pathlib import Path
from typing import Any, Iterable

import numpy as np


SCRIPT = Path(__file__).resolve()
OUT = SCRIPT.parent
ROOT = SCRIPT.parents[3]
STRESS_DIR = OUT.parent / "canonical_stress_20260822"
STRESS_SCRIPT = STRESS_DIR / "run_canonical_stress.py"
GENERAL_MANIFEST = STRESS_DIR / "candidate_fix_source" / "cohort_manifest.csv"
CANDIDATE_FIX = OUT.parent / "canonical_fix_20260822"
PBE_ZIP = (
    ROOT
    / "revision_experiments"
    / "source_cache"
    / "official_tmqmg_github_xyz"
    / "tmQMg_xyz.zip"
)

SELECTION_SEED = 20260822
TRANSFORM_SEED = 151_120_260
TARGET_N = 1000
DEFAULT_K = 100
MAX_EXACT_CANDIDATES = 100_000
CN_SCOPE = (4, 5, 6)
EXPECTED_UNIQUE_GRAPHS = 6427

_WORKER: dict[str, Any] = {}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=("scan", "run", "all"), default="all")
    parser.add_argument("--workers", type=int, default=min(8, os.cpu_count() or 1))
    parser.add_argument("--k", type=int, default=DEFAULT_K)
    parser.add_argument("--target-n", type=int, default=TARGET_N)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def sha256_text(value: str) -> str:
    return sha256_bytes(value.encode("utf-8"))


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def score(namespace: str, *values: object) -> str:
    payload = "|".join([str(SELECTION_SEED), namespace, *(str(v) for v in values)])
    return sha256_text(payload)


def percentile(values: list[float], q: float) -> float | None:
    if not values:
        return None
    return float(np.percentile(np.asarray(values, dtype=float), q))


def write_csv(path: Path, fields: list[str], rows: Iterable[dict[str, Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def load_stress_module() -> Any:
    spec = importlib.util.spec_from_file_location("locked_stress_harness", STRESS_SCRIPT)
    if spec is None or spec.loader is None:
        raise ImportError(STRESS_SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def unique_graph_population() -> tuple[list[dict[str, Any]], dict[str, Any]]:
    stress = load_stress_module()
    candidates, join_audit = stress.load_joined_candidates()
    by_graph: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in candidates:
        by_graph[str(row["graph_sha256"])].append(dict(row))
    representatives = []
    for graph_sha256, rows in by_graph.items():
        chosen = min(
            rows,
            key=lambda row: score(
                "unique_graph_representative", graph_sha256, row["mol_id"]
            ),
        )
        chosen["graph_representative_score_sha256"] = score(
            "unique_graph_representative", graph_sha256, chosen["mol_id"]
        )
        representatives.append(chosen)
    representatives.sort(key=lambda row: row["graph_representative_score_sha256"])
    if len(representatives) != EXPECTED_UNIQUE_GRAPHS:
        raise ValueError(
            f"Expected {EXPECTED_UNIQUE_GRAPHS:,} unique graphs; "
            f"found {len(representatives):,}"
        )
    return representatives, join_audit


def load_payloads(rows: list[dict[str, Any]]) -> dict[str, bytes]:
    payloads: dict[str, bytes] = {}
    with zipfile.ZipFile(PBE_ZIP) as archive:
        names = set(archive.namelist())
        for row in rows:
            member = str(row["source_member"])
            if member not in names:
                raise FileNotFoundError(member)
            payload = archive.read(member)
            if sha256_bytes(payload) != str(row["source_payload_sha256"]):
                raise ValueError(f"Source hash mismatch: {row['mol_id']}")
            payloads[str(row["mol_id"])] = payload
    return payloads


def initialize_worker() -> None:
    release_text = str(CANDIDATE_FIX.resolve())
    if release_text not in sys.path:
        sys.path.insert(0, release_text)
    from rdkit import RDLogger

    RDLogger.DisableLog("rdApp.*")
    import coordrep
    import coordrep.canonical.canonicalize as canonical_module
    from coordrep import encode_molecule
    from coordrep.io.tmqm_reader import Atom, RawMolecule
    from coordrep.serialize.parse_string import parse_constraint_block

    package_path = Path(coordrep.__file__).resolve()
    if CANDIDATE_FIX.resolve() not in package_path.parents:
        raise ImportError(f"CoordRep import escaped candidate root: {package_path}")
    _WORKER.update(
        {
            "Atom": Atom,
            "RawMolecule": RawMolecule,
            "encode_molecule": encode_molecule,
            "canonical": canonical_module,
            "parse_constraint_block": parse_constraint_block,
        }
    )


def parse_xyz(payload: bytes, mol_id: str) -> Any:
    lines = payload.decode("utf-8").splitlines()
    n_atoms = int(lines[0].strip())
    if len(lines) < n_atoms + 2:
        raise ValueError(f"Truncated XYZ: {mol_id}")
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
        mol_id=mol_id, atoms=atoms, bond_orders=None, properties={}
    )


def permute_molecule(molecule: Any, selection_order: int, variant_index: int) -> Any:
    rng = np.random.default_rng(
        np.random.SeedSequence(
            [TRANSFORM_SEED, int(selection_order), int(variant_index)]
        )
    )
    order = rng.permutation(np.arange(len(molecule.atoms), dtype=np.int64))
    atoms = []
    for new_index, old_raw in enumerate(order):
        old_index = int(old_raw)
        source = molecule.atoms[old_index]
        atoms.append(
            _WORKER["Atom"](
                index=new_index,
                element=source.element,
                x=float(source.x),
                y=float(source.y),
                z=float(source.z),
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


def payload_key(ligand: Any) -> tuple[Any, ...]:
    module = _WORKER["canonical"]
    return module._ligand_payload_key(ligand)


def record_orbit_metrics(record: Any) -> dict[str, Any]:
    module = _WORKER["canonical"]
    payload_groups: dict[tuple[Any, ...], list[Any]] = defaultdict(list)
    legal_groups: dict[tuple[Any, ...], list[Any]] = defaultdict(list)
    for ligand in record.ligands:
        payload_groups[payload_key(ligand)].append(ligand)
        key = module._ligand_equivalence_key(ligand)
        if key is not None:
            legal_groups[key].append(ligand)

    repeated_groups = [group for group in payload_groups.values() if len(group) > 1]
    missing_identity_groups = [
        group
        for group in repeated_groups
        if any(not ligand.attachment_set_key for ligand in group)
    ]
    duplicate_sizes = sorted(
        (len(group) for group in legal_groups.values() if len(group) > 1),
        reverse=True,
    )
    donor_orbit_sizes = sorted(
        (
            len(orbit)
            for ligand in record.ligands
            for orbit in ligand.donor_rank_orbits
            if len(orbit) > 1
        ),
        reverse=True,
    )
    candidate_count = 1
    for size in duplicate_sizes + donor_orbit_sizes:
        candidate_count *= factorial(size)
    return {
        "candidate_count": candidate_count,
        "challenge_eligible": candidate_count > 1 and not missing_identity_groups,
        "repeated_payload": bool(repeated_groups),
        "repeated_payload_group_sizes": repeated_groups and sorted(
            (len(group) for group in repeated_groups), reverse=True
        ) or [],
        "legal_ligand_group_sizes": duplicate_sizes,
        "donor_rank_orbit_sizes": donor_orbit_sizes,
        "missing_attachment_identity_group_count": len(missing_identity_groups),
        "above_exact_limit": candidate_count > MAX_EXACT_CANDIDATES,
    }


_METAL_RE = re.compile(r"^\[Metal:[A-Z][a-z]?(?:\|[^\[\]]+)*\]$")
_SHAPE_RE = re.compile(r"^<(?:ShapeBest|ShapeStatus):[^<>]+>$")
_REL_RE = re.compile(
    r"^\{(?:trans|cis):L\d+:[A-Z][a-z]?:\d+--L\d+:[A-Z][a-z]?:\d+\}$"
)
_FM_RE = re.compile(r"^\{fm:(?:fac|mer)\}$")
_LIG_RE = re.compile(r"^(L\d+)=(SMILES|FORMULA):(.+)$")


def grammar_roundtrip(value: str) -> bool:
    """Independently parse, validate, and byte-reconstruct a mononuclear record."""
    if not value.startswith("["):
        return False
    metal_end = value.find("]")
    if metal_end < 0:
        return False
    metal = value[: metal_end + 1]
    if not _METAL_RE.fullmatch(metal):
        return False
    pos = metal_end + 1
    if pos >= len(value) or value[pos] != "<":
        return False
    shape_end = value.find(">", pos)
    if shape_end < 0:
        return False
    shape = value[pos : shape_end + 1]
    if not _SHAPE_RE.fullmatch(shape):
        return False
    pos = shape_end + 1
    relations: list[str] = []
    while pos < len(value) and value[pos] == "{":
        end = value.find("}", pos)
        if end < 0:
            return False
        relation = value[pos : end + 1]
        if not (_REL_RE.fullmatch(relation) or _FM_RE.fullmatch(relation)):
            return False
        relations.append(relation)
        pos = end + 1
    ligand_text = value[pos:]
    if not (ligand_text.startswith("|") and ligand_text.endswith("|")):
        return False
    ligand_items = ligand_text[1:-1].split("|") if ligand_text != "||" else []
    parsed_ligands: list[tuple[str, str, str]] = []
    for item in ligand_items:
        match = _LIG_RE.fullmatch(item)
        if match is None:
            return False
        parsed_ligands.append(match.groups())
    if [item[0] for item in parsed_ligands] != [
        f"L{index}" for index in range(1, len(parsed_ligands) + 1)
    ]:
        return False
    parsed_constraint = _WORKER["parse_constraint_block"]("".join(relations))
    relation_rebuilt = parsed_constraint.to_string()
    if relation_rebuilt != "".join(relations):
        return False
    ligand_rebuilt = "|" + "|".join(
        f"{lig_id}={provenance}:{payload}"
        for lig_id, provenance, payload in parsed_ligands
    ) + "|"
    rebuilt = metal + shape + relation_rebuilt + ligand_rebuilt
    return rebuilt == value


def signature_only_string(record: Any) -> str:
    """Initial invariant sort/relabel/constraint rewrite without exact search."""
    module = _WORKER["canonical"]
    cc = deepcopy(record)
    module._validate_ligand_attachment_metadata(cc.ligands)
    sorted_ligands = sorted(cc.ligands, key=lambda ligand: ligand.get_sort_key())
    id_map: dict[str, str] = {}
    for index, ligand in enumerate(sorted_ligands, start=1):
        old_id = ligand.lig_id
        ligand.lig_id = f"L{index}"
        id_map[old_id] = ligand.lig_id
    cc.ligands = sorted_ligands
    cc.constraints = module._rewrite_constraints(cc.constraints, id_map)
    if cc.shape is not None:
        cc.shape.round(decimals=2)
    cc.is_canonical = True
    return cc.to_string()


def exact_string(record: Any) -> str:
    return record.canonicalize().to_string()


def screen_worker(task: tuple[dict[str, Any], bytes]) -> dict[str, Any]:
    row, payload = task
    started = time.perf_counter()
    output = {
        "mol_id": row["mol_id"],
        "fold": row["fold"],
        "metal": row["metal"],
        "cn": row["cn"],
        "graph_sha256": row["graph_sha256"],
        "source_member": row["source_member"],
        "source_payload_sha256": row["source_payload_sha256"],
        "graph_representative_score_sha256": row["graph_representative_score_sha256"],
        "screen_error": "",
    }
    try:
        molecule = parse_xyz(payload, str(row["mol_id"]))
        encode_started = time.perf_counter()
        record = _WORKER["encode_molecule"](molecule)
        encode_seconds = time.perf_counter() - encode_started
        if record.metal.element != row["metal"] or len(record.graph.donor_indices) != int(row["cn"]):
            raise ValueError(
                f"metal/CN changed to {record.metal.element}/CN{len(record.graph.donor_indices)}"
            )
        metrics = record_orbit_metrics(record)
        exact_started = time.perf_counter()
        full = exact_string(record)
        exact_seconds = time.perf_counter() - exact_started
        output.update(metrics)
        output.update(
            {
                "n_ligands": len(record.ligands),
                "encode_seconds": encode_seconds,
                "exact_seconds": exact_seconds,
                "exact_sha256": sha256_text(full),
                "grammar_roundtrip": grammar_roundtrip(full),
                "screen_seconds": time.perf_counter() - started,
            }
        )
    except Exception as error:
        output.update(
            {
                "challenge_eligible": False,
                "screen_error": f"{type(error).__name__}: {error}",
                "screen_seconds": time.perf_counter() - started,
            }
        )
    return output


def parallel_map(
    worker: Any,
    tasks: list[Any],
    workers: int,
    label: str,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    started = time.perf_counter()
    with ProcessPoolExecutor(max_workers=workers, initializer=initialize_worker) as pool:
        futures = [pool.submit(worker, task) for task in tasks]
        for completed, future in enumerate(as_completed(futures), start=1):
            rows.append(future.result())
            if completed % max(1, min(100, len(tasks))) == 0 or completed == len(tasks):
                elapsed = time.perf_counter() - started
                print(
                    f"[{utc_now()}] {label}: {completed:,}/{len(tasks):,}; "
                    f"{completed/elapsed:.2f} structures/s",
                    flush=True,
                )
    return rows


def hamilton_quotas(counts: dict[int, int], target: int) -> dict[int, int]:
    total = sum(counts.values())
    if target >= total:
        return dict(counts)
    raw = {cn: target * counts[cn] / total for cn in CN_SCOPE}
    quotas = {cn: int(math.floor(raw[cn])) for cn in CN_SCOPE}
    remainder = target - sum(quotas.values())
    order = sorted(CN_SCOPE, key=lambda cn: (-(raw[cn] - quotas[cn]), cn))
    for cn in order[:remainder]:
        quotas[cn] += 1
    return quotas


def load_general_graphs() -> set[str]:
    with GENERAL_MANIFEST.open(newline="", encoding="utf-8") as handle:
        return {str(row["graph_sha256"]) for row in csv.DictReader(handle)}


def select_challenge(screen_rows: list[dict[str, Any]], target_n: int) -> tuple[list[dict[str, Any]], dict[int, int]]:
    eligible = [
        row
        for row in screen_rows
        if row.get("challenge_eligible")
        and not row.get("screen_error")
        and not row.get("above_exact_limit")
    ]
    counts = Counter(int(row["cn"]) for row in eligible)
    quotas = hamilton_quotas({cn: counts[cn] for cn in CN_SCOPE}, min(target_n, len(eligible)))
    general_graphs = load_general_graphs()
    selected: list[dict[str, Any]] = []
    for cn in CN_SCOPE:
        ranked = sorted(
            (row for row in eligible if int(row["cn"]) == cn),
            key=lambda row: score(
                "challenge_selection", cn, row["graph_sha256"], row["mol_id"]
            ),
        )
        for rank, row in enumerate(ranked[: quotas[cn]], start=1):
            chosen = dict(row)
            chosen["cn_selection_rank"] = rank
            chosen["selection_score_sha256"] = score(
                "challenge_selection", cn, row["graph_sha256"], row["mol_id"]
            )
            chosen["in_general_n1000"] = row["graph_sha256"] in general_graphs
            selected.append(chosen)
    selected.sort(key=lambda row: row["selection_score_sha256"])
    for index, row in enumerate(selected, start=1):
        row["selection_order"] = index
    return selected, quotas


def evaluate_challenge(task: tuple[dict[str, Any], bytes, int]) -> dict[str, Any]:
    row, payload, k = task
    output = {
        "selection_order": row["selection_order"],
        "mol_id": row["mol_id"],
        "fold": row["fold"],
        "metal": row["metal"],
        "cn": row["cn"],
        "graph_sha256": row["graph_sha256"],
        "candidate_count": row["candidate_count"],
        "k": k,
        "signature_matches": 0,
        "exact_matches": 0,
        "signature_roundtrip_passes": 0,
        "exact_roundtrip_passes": 0,
        "candidate_count_matches": 0,
        "variant_errors": 0,
        "error_examples_json": "[]",
    }
    signature_times: list[float] = []
    exact_times: list[float] = []
    encode_times: list[float] = []
    errors: list[str] = []
    try:
        molecule = parse_xyz(payload, str(row["mol_id"]))
        baseline = _WORKER["encode_molecule"](molecule)
        baseline_metrics = record_orbit_metrics(baseline)
        if int(baseline_metrics["candidate_count"]) != int(row["candidate_count"]):
            raise ValueError("baseline candidate count differs from frozen screen")
        start = time.perf_counter()
        signature_reference = signature_only_string(baseline)
        signature_times.append(time.perf_counter() - start)
        start = time.perf_counter()
        exact_reference = exact_string(baseline)
        exact_times.append(time.perf_counter() - start)
        if not grammar_roundtrip(signature_reference) or not grammar_roundtrip(exact_reference):
            raise ValueError("baseline grammar round trip failed")
        output["signature_reference_sha256"] = sha256_text(signature_reference)
        output["exact_reference_sha256"] = sha256_text(exact_reference)
        for variant_index in range(k):
            try:
                variant = permute_molecule(
                    molecule, int(row["selection_order"]), variant_index
                )
                start = time.perf_counter()
                record = _WORKER["encode_molecule"](variant)
                encode_times.append(time.perf_counter() - start)
                metrics = record_orbit_metrics(record)
                if int(metrics["candidate_count"]) == int(row["candidate_count"]):
                    output["candidate_count_matches"] += 1
                start = time.perf_counter()
                signature_value = signature_only_string(record)
                signature_times.append(time.perf_counter() - start)
                start = time.perf_counter()
                exact_value = exact_string(record)
                exact_times.append(time.perf_counter() - start)
                if signature_value == signature_reference:
                    output["signature_matches"] += 1
                if exact_value == exact_reference:
                    output["exact_matches"] += 1
                if grammar_roundtrip(signature_value):
                    output["signature_roundtrip_passes"] += 1
                if grammar_roundtrip(exact_value):
                    output["exact_roundtrip_passes"] += 1
            except Exception as error:
                output["variant_errors"] += 1
                if len(errors) < 5:
                    errors.append(
                        f"variant {variant_index}: {type(error).__name__}: {error}"
                    )
        output.update(
            {
                "signature_fully_collapsed": output["signature_matches"] == k,
                "exact_fully_collapsed": output["exact_matches"] == k,
                "candidate_count_stable": output["candidate_count_matches"] == k,
                "signature_roundtrip_all": output["signature_roundtrip_passes"] == k,
                "exact_roundtrip_all": output["exact_roundtrip_passes"] == k,
                "encode_times": encode_times,
                "signature_times": signature_times,
                "exact_times": exact_times,
                "error_examples_json": json.dumps(errors, ensure_ascii=True),
                "fatal_error": "",
            }
        )
    except Exception as error:
        output.update(
            {
                "signature_fully_collapsed": False,
                "exact_fully_collapsed": False,
                "candidate_count_stable": False,
                "signature_roundtrip_all": False,
                "exact_roundtrip_all": False,
                "encode_times": encode_times,
                "signature_times": signature_times,
                "exact_times": exact_times,
                "fatal_error": f"{type(error).__name__}: {error}",
            }
        )
    return output


def run_negative_controls() -> dict[str, Any]:
    from coordrep.canonical.canonicalize import _ligand_equivalence_key
    from coordrep.core import (
        AssemblyGraph,
        ConstraintSet,
        CoordComplex,
        LigandModule,
        MetalState,
    )
    from coordrep.graph.ligand_module import LigandExtractor
    from coordrep.io.tmqm_reader import Atom
    import networkx as nx

    atoms = [Atom(index=i, element="C", x=float(i), y=0.0, z=0.0) for i in range(6)]
    graph = nx.Graph()
    graph.add_nodes_from((i, {"element": "C"}) for i in range(6))
    for index in range(6):
        graph.add_edge(index, (index + 1) % 6, bond_order=1.0)
    extractor = LigandExtractor()

    def descriptor(donors: list[int]) -> tuple[str, list[str]]:
        key, donor_map = extractor._generate_attachment_descriptors(
            atoms, graph, set(range(6)), donors
        )
        return key, [donor_map[index] for index in donors]

    adjacent_key, adjacent_donors = descriptor([0, 1])
    opposite_key, opposite_donors = descriptor([0, 3])

    def ligand(key: str, donor_keys: list[str], attach: list[int]) -> Any:
        return LigandModule(
            lig_id="source",
            smiles="C1CCCCC1",
            attach_atoms=attach,
            donor_elements=["C", "C"],
            dent=2,
            attachment_set_key=key,
            donor_attachment_keys=donor_keys,
        )

    def complex_for(item: Any) -> Any:
        return CoordComplex(
            metal=MetalState("Ir", oxidation=3),
            ligands=[item],
            graph=AssemblyGraph(
                metal_idx=0,
                donor_indices=[1, 2],
                lig_assignments={},
                bond_orders={},
            ),
            shape=None,
            constraints=ConstraintSet(),
        )

    adjacent_ligand = ligand(adjacent_key, adjacent_donors, [101, 102])
    opposite_ligand = ligand(opposite_key, opposite_donors, [401, 404])
    adjacent_string = complex_for(adjacent_ligand).canonicalize().to_string()
    opposite_string = complex_for(opposite_ligand).canonicalize().to_string()

    unresolved = [
        LigandModule(
            lig_id=f"U{index}",
            smiles="c1ccncc1",
            attach_atoms=[10 + index],
            donor_elements=["N"],
            dent=1,
        )
        for index in range(2)
    ]
    fail_closed = False
    fail_closed_error = ""
    try:
        CoordComplex(
            metal=MetalState("Ir", oxidation=3),
            ligands=unresolved,
            graph=AssemblyGraph(
                metal_idx=0,
                donor_indices=[1, 2],
                lig_assignments={},
                bond_orders={},
            ),
            shape=None,
            constraints=ConstraintSet(),
        ).canonicalize()
    except ValueError as error:
        fail_closed = "attachment_set_key" in str(error)
        fail_closed_error = str(error)

    return {
        "non_automorphic_pairs_tested": 1,
        "attachment_keys_distinct": adjacent_key != opposite_key,
        "equivalence_keys_distinct": (
            _ligand_equivalence_key(adjacent_ligand)
            != _ligand_equivalence_key(opposite_ligand)
        ),
        "canonical_strings_distinct": adjacent_string != opposite_string,
        "false_merges": int(adjacent_string == opposite_string),
        "adjacent_attachment_key_sha256": sha256_text(adjacent_key),
        "opposite_attachment_key_sha256": sha256_text(opposite_key),
        "adjacent_record_sha256": sha256_text(adjacent_string),
        "opposite_record_sha256": sha256_text(opposite_string),
        "missing_metadata_fail_closed": fail_closed,
        "missing_metadata_error": fail_closed_error,
    }


def list_distribution(rows: list[dict[str, Any]], field: str) -> list[dict[str, int]]:
    counts = Counter(int(row[field]) for row in rows if str(row.get(field, "")) != "")
    return [{field: key, "records": counts[key]} for key in sorted(counts)]


def flatten_times(rows: list[dict[str, Any]], field: str) -> list[float]:
    return [float(value) for row in rows for value in row.get(field, [])]


def timing_summary(values: list[float]) -> dict[str, Any]:
    return {
        "calls": len(values),
        "median_seconds": percentile(values, 50),
        "p95_seconds": percentile(values, 95),
        "max_seconds": max(values) if values else None,
    }


def write_checksums() -> None:
    rows = []
    for path in sorted(OUT.iterdir(), key=lambda item: item.name.lower()):
        if not path.is_file() or path.name == "checksums.sha256":
            continue
        rows.append(f"{sha256_file(path)}  {path.name}")
    (OUT / "checksums.sha256").write_text("\n".join(rows) + "\n", encoding="utf-8")


def clean_outputs(force: bool, stage: str) -> None:
    generated = [
        "screen_manifest.csv",
        "manifest.csv",
        "source_data.csv",
        "summary.json",
        "checksums.sha256",
        "run_log.txt",
    ]
    if stage in {"scan", "all"}:
        existing = [OUT / name for name in generated if (OUT / name).exists()]
        if existing and not force:
            raise FileExistsError("Outputs exist; pass --force")
        for path in existing:
            path.unlink()


def main() -> int:
    args = parse_args()
    if args.workers < 1 or args.k < 1 or args.target_n < 1:
        raise ValueError("workers, k, and target-n must be positive")
    clean_outputs(args.force, args.stage)
    log_path = OUT / "run_log.txt"

    def log(message: str) -> None:
        line = f"[{utc_now()}] {message}"
        print(line, flush=True)
        with log_path.open("a", encoding="utf-8", newline="\n") as handle:
            handle.write(line + "\n")

    started_utc = utc_now()
    started = time.perf_counter()
    population, join_audit = unique_graph_population()
    payloads = load_payloads(population)
    log(f"verified {len(population):,} graph representatives and source payloads")

    screen_path = OUT / "screen_manifest.csv"
    if args.stage in {"scan", "all"}:
        screen_rows = parallel_map(
            screen_worker,
            [(row, payloads[row["mol_id"]]) for row in population],
            args.workers,
            "screen",
        )
        screen_rows.sort(key=lambda row: str(row["graph_representative_score_sha256"]))
        screen_fields = [
            "mol_id", "fold", "metal", "cn", "graph_sha256", "source_member",
            "source_payload_sha256", "graph_representative_score_sha256",
            "n_ligands", "repeated_payload", "repeated_payload_group_sizes",
            "legal_ligand_group_sizes", "donor_rank_orbit_sizes", "candidate_count",
            "challenge_eligible", "missing_attachment_identity_group_count",
            "above_exact_limit", "grammar_roundtrip", "exact_sha256",
            "encode_seconds", "exact_seconds", "screen_seconds", "screen_error",
        ]
        serialized_screen = []
        for row in screen_rows:
            rendered = dict(row)
            for field in (
                "repeated_payload_group_sizes",
                "legal_ligand_group_sizes",
                "donor_rank_orbit_sizes",
            ):
                rendered[field] = json.dumps(row.get(field, []), separators=(",", ":"))
            serialized_screen.append(rendered)
        write_csv(screen_path, screen_fields, serialized_screen)
        log(
            "challenge eligible: "
            f"{sum(bool(row.get('challenge_eligible')) for row in screen_rows):,}/"
            f"{len(screen_rows):,}"
        )
    else:
        if not screen_path.exists():
            raise FileNotFoundError("Run --stage scan first")
        with screen_path.open(newline="", encoding="utf-8") as handle:
            screen_rows = list(csv.DictReader(handle))
        for row in screen_rows:
            for field in (
                "repeated_payload_group_sizes",
                "legal_ligand_group_sizes",
                "donor_rank_orbit_sizes",
            ):
                row[field] = json.loads(row.get(field) or "[]")
            for field in (
                "challenge_eligible", "repeated_payload", "above_exact_limit",
                "grammar_roundtrip",
            ):
                row[field] = str(row.get(field, "")).lower() == "true"
            for field in (
                "candidate_count", "missing_attachment_identity_group_count", "cn", "fold",
            ):
                if str(row.get(field, "")):
                    row[field] = int(row[field])

    selected, quotas = select_challenge(screen_rows, args.target_n)
    manifest_fields = [
        "selection_order", "cn_selection_rank", "selection_score_sha256",
        "mol_id", "fold", "metal", "cn", "graph_sha256", "source_member",
        "source_payload_sha256", "candidate_count", "repeated_payload_group_sizes",
        "legal_ligand_group_sizes", "donor_rank_orbit_sizes", "in_general_n1000",
    ]
    manifest_rows = []
    for row in selected:
        rendered = dict(row)
        for field in (
            "repeated_payload_group_sizes",
            "legal_ligand_group_sizes",
            "donor_rank_orbit_sizes",
        ):
            rendered[field] = json.dumps(row.get(field, []), separators=(",", ":"))
        manifest_rows.append(rendered)
    write_csv(OUT / "manifest.csv", manifest_fields, manifest_rows)
    log(f"locked challenge cohort N={len(selected):,}; quotas={quotas}")

    if args.stage == "scan":
        scan_summary = {
            "status": "SCAN_COMPLETE_RUN_NOT_EXECUTED",
            "population_n": len(screen_rows),
            "eligible_n": sum(bool(row.get("challenge_eligible")) for row in screen_rows),
            "candidate_count_distribution": list_distribution(
                [row for row in screen_rows if row.get("challenge_eligible")],
                "candidate_count",
            ),
            "selected_n": len(selected),
            "selected_cn_quotas": quotas,
        }
        (OUT / "summary.json").write_text(
            json.dumps(scan_summary, indent=2) + "\n", encoding="utf-8"
        )
        write_checksums()
        return 0

    run_rows = parallel_map(
        evaluate_challenge,
        [(row, payloads[row["mol_id"]], args.k) for row in selected],
        args.workers,
        "challenge",
    )
    run_rows.sort(key=lambda row: int(row["selection_order"]))
    # The process-pool initializer runs only in child processes.  Initialize
    # the parent once before executing the local specificity controls.
    initialize_worker()
    negative_controls = run_negative_controls()
    source_fields = [
        "selection_order", "mol_id", "fold", "metal", "cn", "graph_sha256",
        "candidate_count", "k", "signature_matches", "exact_matches",
        "signature_fully_collapsed", "exact_fully_collapsed",
        "signature_roundtrip_passes", "exact_roundtrip_passes",
        "signature_roundtrip_all", "exact_roundtrip_all", "candidate_count_matches",
        "candidate_count_stable", "variant_errors", "signature_reference_sha256",
        "exact_reference_sha256", "encode_median_seconds", "signature_median_seconds",
        "exact_median_seconds", "error_examples_json", "fatal_error",
    ]
    source_rows = []
    for row in run_rows:
        rendered = dict(row)
        rendered["encode_median_seconds"] = percentile(row.get("encode_times", []), 50)
        rendered["signature_median_seconds"] = percentile(row.get("signature_times", []), 50)
        rendered["exact_median_seconds"] = percentile(row.get("exact_times", []), 50)
        source_rows.append(rendered)
    write_csv(OUT / "source_data.csv", source_fields, source_rows)

    eligible_rows = [row for row in screen_rows if row.get("challenge_eligible")]
    signature_matches = sum(int(row["signature_matches"]) for row in run_rows)
    exact_matches = sum(int(row["exact_matches"]) for row in run_rows)
    total_variants = len(run_rows) * args.k
    encode_times = flatten_times(run_rows, "encode_times")
    signature_times = flatten_times(run_rows, "signature_times")
    exact_times = flatten_times(run_rows, "exact_times")
    screen_errors = [row for row in screen_rows if row.get("screen_error")]
    fatal_rows = [row for row in run_rows if row.get("fatal_error")]
    variant_errors = sum(int(row["variant_errors"]) for row in run_rows)
    full_pass = (
        exact_matches == total_variants
        and all(bool(row["exact_fully_collapsed"]) for row in run_rows)
        and all(bool(row["exact_roundtrip_all"]) for row in run_rows)
        and all(bool(row["candidate_count_stable"]) for row in run_rows)
        and not fatal_rows
        and variant_errors == 0
        and negative_controls["false_merges"] == 0
        and negative_controls["missing_metadata_fail_closed"]
    )
    summary = {
        "status": "TASK_LOCAL_ORBIT_CHALLENGE_PASS" if full_pass else "TASK_LOCAL_ORBIT_CHALLENGE_FAIL",
        "claim_boundary": (
            "Fixed graph-disjoint tmQMg/PBE molecular CN4-6 snapshots; "
            "not April-2025-CSD prevalence. Frozen task-local candidate, not public release."
        ),
        "protocol_sha256": sha256_file(OUT / "protocol.md"),
        "source_population": {
            "strict_eligible_unique_graphs": len(screen_rows),
            "expected_unique_graphs": EXPECTED_UNIQUE_GRAPHS,
            "join_audit": join_audit,
            "screen_errors": len(screen_errors),
            "screen_grammar_roundtrip_failures": sum(
                not bool(row.get("grammar_roundtrip"))
                for row in screen_rows
                if not row.get("screen_error")
            ),
            "repeated_payload_records": sum(
                bool(row.get("repeated_payload")) for row in screen_rows
            ),
            "challenge_eligible_records": len(eligible_rows),
            "challenge_eligible_by_cn": dict(
                sorted(Counter(int(row["cn"]) for row in eligible_rows).items())
            ),
            "candidate_count_distribution": list_distribution(
                eligible_rows, "candidate_count"
            ),
            "missing_attachment_identity_records": sum(
                int(row.get("missing_attachment_identity_group_count") or 0) > 0
                for row in screen_rows
            ),
            "above_exact_limit_records": sum(
                bool(row.get("above_exact_limit")) for row in screen_rows
            ),
        },
        "executed_cohort": {
            "n": len(selected),
            "k_atom_reindexings": args.k,
            "total_variant_encodings": total_variants,
            "cn_quotas": quotas,
            "overlap_with_general_n1000": sum(bool(row["in_general_n1000"]) for row in selected),
            "unique_graphs": len({row["graph_sha256"] for row in selected}),
        },
        "signature_only_ablation": {
            "definition": (
                "candidate invariant ligand sort + L relabel + constraint rewrite/sort + "
                "shape rounding; exact legal-orbit enumeration and whole-record minimization skipped"
            ),
            "exact_variant_matches": signature_matches,
            "total_variants": total_variants,
            "variant_match_fraction": signature_matches / total_variants,
            "fully_collapsed_structures": sum(
                bool(row["signature_fully_collapsed"]) for row in run_rows
            ),
            "n_structures": len(run_rows),
            "grammar_roundtrip_passes": sum(
                int(row["signature_roundtrip_passes"]) for row in run_rows
            ),
        },
        "full_exact_coordrep": {
            "exact_variant_matches": exact_matches,
            "total_variants": total_variants,
            "variant_match_fraction": exact_matches / total_variants,
            "fully_collapsed_structures": sum(
                bool(row["exact_fully_collapsed"]) for row in run_rows
            ),
            "n_structures": len(run_rows),
            "grammar_roundtrip_passes": sum(
                int(row["exact_roundtrip_passes"]) for row in run_rows
            ),
            "candidate_count_stable_structures": sum(
                bool(row["candidate_count_stable"]) for row in run_rows
            ),
            "variant_errors": variant_errors,
            "fatal_structures": len(fatal_rows),
        },
        "specificity_controls": negative_controls,
        "timing": {
            "scope": "descriptive wall-clock seconds per call on this machine",
            "encode_variant": timing_summary(encode_times),
            "signature_only": timing_summary(signature_times),
            "full_exact": timing_summary(exact_times),
        },
        "environment": {
            "created_utc": utc_now(),
            "started_utc": started_utc,
            "elapsed_seconds": time.perf_counter() - started,
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "workers": args.workers,
            "script_sha256": sha256_file(SCRIPT),
            "candidate_fix_root": str(CANDIDATE_FIX.resolve()),
            "candidate_canonicalizer_sha256": sha256_file(
                CANDIDATE_FIX / "coordrep" / "canonical" / "canonicalize.py"
            ),
            "source_zip_sha256": sha256_file(PBE_ZIP),
            "selection_seed": SELECTION_SEED,
            "transform_seed": TRANSFORM_SEED,
        },
    }
    (OUT / "summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    log(
        f"complete: status={summary['status']}; exact={exact_matches:,}/{total_variants:,}; "
        f"signature={signature_matches:,}/{total_variants:,}"
    )
    write_checksums()
    return 0 if full_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
