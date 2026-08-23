#!/usr/bin/env python3
"""Launch, resume, validate, and merge a sharded frozen-CSD audit.

The launcher is deliberately orchestration-only.  Each worker writes an
independent restartable shard; completed shards are immutable and skipped on
subsequent invocations.  Source-fidelity validation is run only after all
scan shards complete successfully, and merging is the final gated step.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import ccdc
import numpy as np
import scipy


BASE = Path(__file__).resolve().parent.parent
SCRIPTS = BASE / "scripts"
if str(BASE) not in sys.path:
    sys.path.insert(0, str(BASE))

from coordrep.audit.csd_adapter import PROTOCOL_ID


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--database", required=True)
    parser.add_argument("--source-release", required=True)
    parser.add_argument("--database-sha256", required=True)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--num-shards", type=int, default=32)
    parser.add_argument("--max-workers", type=int, default=16)
    parser.add_argument("--expected-database-entries", type=int, required=True)
    parser.add_argument("--expected-tm-entries", type=int, required=True)
    parser.add_argument("--checkpoint-every", type=int, default=1000)
    parser.add_argument("--progress-every", type=int, default=10000)
    parser.add_argument("--max-targets", type=int)
    parser.add_argument(
        "--stage",
        choices=("scan", "validate", "merge", "all"),
        default="all",
    )
    return parser.parse_args()


def _shard_dir(root: Path, shard_id: int, num_shards: int) -> Path:
    return root / "shards" / f"shard_{shard_id:04d}_of_{num_shards:04d}"


def _run_logged(command: Sequence[str], log_path: Path) -> Tuple[int, float]:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    started = time.time()
    with log_path.open("ab") as handle:
        handle.write(("\nCOMMAND " + json.dumps(list(command)) + "\n").encode("utf-8"))
        handle.flush()
        environment = os.environ.copy()
        environment.update({
            "OMP_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
            "OPENBLAS_NUM_THREADS": "1",
            "NUMEXPR_NUM_THREADS": "1",
        })
        completed = subprocess.run(
            list(command),
            stdout=handle,
            stderr=subprocess.STDOUT,
            env=environment,
            check=False,
        )
        elapsed = time.time() - started
        handle.write(
            (f"\nEXIT {completed.returncode} ELAPSED_SECONDS {elapsed:.3f}\n").encode("utf-8")
        )
    return completed.returncode, elapsed


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _line_count(path: Path) -> int:
    with path.open("rb") as handle:
        return sum(1 for line in handle if line.strip())


def _expected_sources() -> Dict[str, str]:
    paths = [
        BASE / "coordrep" / "audit" / "records.py",
        BASE / "coordrep" / "audit" / "csd_adapter.py",
        BASE / "coordrep" / "geometry" / "shape.py",
        BASE / "coordrep" / "io" / "tmqm_reader.py",
        SCRIPTS / "run_csd_release_record_shard.py",
    ]
    return {str(path.relative_to(BASE)): _sha256(path) for path in paths}


def _expected_runtime() -> Dict[str, str]:
    return {
        "python": sys.version.replace("\n", " "),
        "implementation": platform.python_implementation(),
        "platform": platform.platform(),
        "ccdc": str(getattr(ccdc, "__version__", "unknown")),
        "numpy": str(np.__version__),
        "scipy": str(scipy.__version__),
    }


def _completed_scan_is_consistent(
    directory: Path, args: argparse.Namespace, shard_id: int
) -> bool:
    summary_path = directory / "SUMMARY.json"
    paths = (
        directory / "RECORDS_INTERNAL.jsonl",
        directory / "ENTRY_OUTCOMES.jsonl",
        directory / "FAILURES_INTERNAL.jsonl",
    )
    if not summary_path.exists() or not all(path.exists() for path in paths):
        return False
    payload = json.loads(summary_path.read_text(encoding="utf-8"))
    counters = payload["counters"]
    config = payload["config"]
    return (
        config.get("protocol_id") == PROTOCOL_ID
        and config.get("source_release") == args.source_release
        and config.get("database_sha256") == args.database_sha256.lower()
        and config.get("shard_id") == shard_id
        and config.get("num_shards") == args.num_shards
        and config.get("max_targets") == args.max_targets
        and config.get("sources") == _expected_sources()
        and config.get("runtime") == _expected_runtime()
        and _line_count(paths[0]) == int(counters.get("records_emitted", 0))
        and _line_count(paths[1]) == int(counters.get("tm_entries", 0))
        and _line_count(paths[2])
        == int(counters.get("controlled_record_failure", 0))
        + int(counters.get("untyped_pipeline_failure", 0))
        + int(counters.get("entry_read_failed", 0))
        + int(counters.get("classification_failed", 0))
    )


def _scan_one(args: argparse.Namespace, shard_id: int) -> Dict[str, object]:
    directory = _shard_dir(args.run_root, shard_id, args.num_shards)
    directory.mkdir(parents=True, exist_ok=True)
    summary = directory / "SUMMARY.json"
    if summary.exists():
        status = (
            "skipped_complete"
            if _completed_scan_is_consistent(directory, args, shard_id)
            else "corrupt_complete"
        )
        return {"shard": shard_id, "stage": "scan", "status": status}
    command: List[str] = [
        sys.executable,
        str(SCRIPTS / "run_csd_release_record_shard.py"),
        "--database",
        args.database,
        "--source-release",
        args.source_release,
        "--database-sha256",
        args.database_sha256.lower(),
        "--output-dir",
        str(directory),
        "--shard-id",
        str(shard_id),
        "--num-shards",
        str(args.num_shards),
        "--checkpoint-every",
        str(args.checkpoint_every),
        "--progress-every",
        str(args.progress_every),
    ]
    if (directory / "checkpoint.json").exists():
        command.append("--resume")
    if args.max_targets is not None:
        command.extend(("--max-targets", str(args.max_targets)))
    code, elapsed = _run_logged(command, directory / "scan.log")
    return {
        "shard": shard_id,
        "stage": "scan",
        "status": "complete" if code == 0 else "failed",
        "returncode": code,
        "elapsed_seconds": round(elapsed, 3),
    }


def _validate_one(args: argparse.Namespace, shard_id: int) -> Dict[str, object]:
    directory = _shard_dir(args.run_root, shard_id, args.num_shards)
    records = directory / "RECORDS_INTERNAL.jsonl"
    output = directory / "SOURCE_FIDELITY.jsonl"
    summary = output.with_suffix(output.suffix + ".summary.json")
    if not records.exists():
        return {"shard": shard_id, "stage": "validate", "status": "missing_records"}
    if summary.exists():
        payload = json.loads(summary.read_text(encoding="utf-8"))
        scan_summary = json.loads((directory / "SUMMARY.json").read_text(encoding="utf-8"))
        consistent = (
            output.exists()
            and int(payload.get("failed", -1)) == 0
            and int(payload.get("records", -1))
            == int(scan_summary["counters"].get("records_emitted", -2))
            and payload.get("records_sha256") == _sha256(records)
            and payload.get("validation_sha256") == _sha256(output)
            and payload.get("validator_sha256")
            == _sha256(SCRIPTS / "validate_csd_release_record_fidelity.py")
            and payload.get("runtime") == _expected_runtime()
            and int(payload.get("database", {}).get("entry_count", -1))
            == args.expected_database_entries
        )
        status = "skipped_complete" if consistent else "failed"
        return {"shard": shard_id, "stage": "validate", "status": status}
    command = [
        sys.executable,
        str(SCRIPTS / "validate_csd_release_record_fidelity.py"),
        "--database",
        args.database,
        "--records",
        str(records),
        "--output",
        str(output),
    ]
    code, elapsed = _run_logged(command, directory / "validate.log")
    status = "complete" if code == 0 and summary.exists() else "failed"
    if status == "complete":
        payload = json.loads(summary.read_text(encoding="utf-8"))
        if int(payload.get("failed", -1)) != 0:
            status = "failed"
    return {
        "shard": shard_id,
        "stage": "validate",
        "status": status,
        "returncode": code,
        "elapsed_seconds": round(elapsed, 3),
    }


def _parallel_stage(args: argparse.Namespace, function) -> List[Dict[str, object]]:
    results: List[Dict[str, object]] = []
    with ThreadPoolExecutor(max_workers=args.max_workers) as executor:
        futures = {
            executor.submit(function, args, shard_id): shard_id
            for shard_id in range(args.num_shards)
        }
        for future in as_completed(futures):
            result = future.result()
            results.append(result)
            print(json.dumps(result, sort_keys=True), flush=True)
    results.sort(key=lambda row: int(row["shard"]))
    return results


def _require_success(results: Sequence[Dict[str, object]], stage: str) -> None:
    accepted = {"complete", "skipped_complete"}
    failed = [row for row in results if row.get("status") not in accepted]
    if failed:
        raise SystemExit(f"{stage} failed for {len(failed)} shard(s): {failed}")


def _validate_gate(args: argparse.Namespace) -> Dict[str, object]:
    total_records = 0
    total_passed = 0
    total_failed = 0
    summaries = []
    reference_validator_sha256 = None
    reference_runtime = None
    expected_validator_sha256 = _sha256(
        SCRIPTS / "validate_csd_release_record_fidelity.py"
    )
    expected_runtime = _expected_runtime()
    for shard_id in range(args.num_shards):
        path = (
            _shard_dir(args.run_root, shard_id, args.num_shards)
            / "SOURCE_FIDELITY.jsonl.summary.json"
        )
        if not path.exists():
            raise SystemExit(f"missing validation summary: {path}")
        payload = json.loads(path.read_text(encoding="utf-8"))
        directory = _shard_dir(args.run_root, shard_id, args.num_shards)
        records_path = directory / "RECORDS_INTERNAL.jsonl"
        validation_path = directory / "SOURCE_FIDELITY.jsonl"
        scan_summary = json.loads((directory / "SUMMARY.json").read_text(encoding="utf-8"))
        if (
            not validation_path.exists()
            or payload.get("records_sha256") != _sha256(records_path)
            or payload.get("validation_sha256") != _sha256(validation_path)
            or int(payload["records"])
            != int(scan_summary["counters"].get("records_emitted", -1))
            or payload.get("validator_sha256") != expected_validator_sha256
            or payload.get("runtime") != expected_runtime
            or int(payload.get("database", {}).get("entry_count", -1))
            != args.expected_database_entries
        ):
            raise SystemExit(f"stale or inconsistent validation artifacts in shard {shard_id}")
        if reference_validator_sha256 is None:
            reference_validator_sha256 = payload.get("validator_sha256")
            reference_runtime = payload.get("runtime")
        if payload.get("validator_sha256") != reference_validator_sha256:
            raise SystemExit(f"validator fingerprint mismatch in shard {shard_id}")
        if payload.get("runtime") != reference_runtime:
            raise SystemExit(f"validator runtime mismatch in shard {shard_id}")
        total_records += int(payload["records"])
        total_passed += int(payload["passed"])
        total_failed += int(payload["failed"])
        summaries.append(payload)
    aggregate = {
        "records": total_records,
        "passed": total_passed,
        "failed": total_failed,
        "pass_rate": total_passed / total_records if total_records else None,
        "shards": args.num_shards,
        "validator_sha256": reference_validator_sha256,
        "runtime": reference_runtime,
        "database_sha256": args.database_sha256.lower(),
        "shard_summaries": summaries,
    }
    destination = args.run_root / "merged" / "SOURCE_FIDELITY_AGGREGATE.json"
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(aggregate, indent=2, sort_keys=True), encoding="utf-8")
    if total_failed:
        raise SystemExit(f"source-fidelity gate failed for {total_failed} record(s)")
    return aggregate


def _merge(args: argparse.Namespace) -> None:
    merged = args.run_root / "merged"
    public_summary_path = merged / "SUMMARY_PUBLIC_AGGREGATE.json"
    if public_summary_path.exists():
        payload = json.loads(public_summary_path.read_text(encoding="utf-8"))
        artifact_paths = {
            "records_internal_sha256": merged / "RECORDS_INTERNAL.jsonl",
            "outcomes_internal_sha256": merged / "ENTRY_OUTCOMES_INTERNAL.jsonl",
            "failures_internal_sha256": merged / "FAILURES_INTERNAL.jsonl",
        }
        artifacts = payload.get("artifacts", {})
        consistent = (
            payload.get("protocol_id") == PROTOCOL_ID
            and payload.get("source_release") == args.source_release
            and str(payload.get("database_sha256", "")).lower()
            == args.database_sha256.lower()
            and int(payload.get("database_entries_evaluated", -1))
            == args.expected_database_entries
            and int(payload.get("transition_metal_entries", -1))
            == args.expected_tm_entries
            and int(payload.get("shards", -1)) == args.num_shards
            and payload.get("source_fingerprints") == _expected_sources()
            and all(
                path.exists() and artifacts.get(key) == _sha256(path)
                for key, path in artifact_paths.items()
            )
        )
        if not consistent:
            raise SystemExit("existing merged output is stale or inconsistent")
        print(json.dumps({"stage": "merge", "status": "skipped_complete"}), flush=True)
        return
    command = [
        sys.executable,
        str(SCRIPTS / "merge_csd_release_record_shards.py"),
        "--root",
        str(args.run_root / "shards"),
        "--num-shards",
        str(args.num_shards),
        "--output-dir",
        str(merged),
        "--expected-database-entries",
        str(args.expected_database_entries),
        "--expected-tm-entries",
        str(args.expected_tm_entries),
    ]
    code, _ = _run_logged(command, args.run_root / "merge.log")
    if code != 0:
        raise SystemExit(f"merge failed; inspect {args.run_root / 'merge.log'}")


def main() -> None:
    args = parse_args()
    if args.num_shards < 1 or args.max_workers < 1:
        raise SystemExit("num-shards and max-workers must be positive")
    args.run_root.mkdir(parents=True, exist_ok=True)
    database_path = Path(args.database)
    if not database_path.is_file():
        raise SystemExit("the frozen audit requires an explicit database file")
    database_digest = hashlib.sha256()
    with database_path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            database_digest.update(block)
    observed_digest = database_digest.hexdigest()
    if observed_digest.lower() != args.database_sha256.lower():
        raise SystemExit(
            f"database SHA-256 mismatch: {observed_digest} != {args.database_sha256.lower()}"
        )
    config_path = args.run_root / "LAUNCH_CONFIG.json"
    config = {
        "database": args.database,
        "source_release": args.source_release,
        "database_sha256": args.database_sha256.lower(),
        "num_shards": args.num_shards,
        "max_workers": args.max_workers,
        "expected_database_entries": args.expected_database_entries,
        "expected_tm_entries": args.expected_tm_entries,
        "max_targets": args.max_targets,
    }
    if config_path.exists():
        existing = json.loads(config_path.read_text(encoding="utf-8"))
        if existing != config:
            raise SystemExit("launch configuration differs from the frozen run configuration")
    else:
        config_path.write_text(json.dumps(config, indent=2, sort_keys=True), encoding="utf-8")

    if args.stage in ("scan", "all"):
        _require_success(_parallel_stage(args, _scan_one), "scan")
    if args.stage in ("validate", "all"):
        _require_success(_parallel_stage(args, _validate_one), "validation")
        aggregate = _validate_gate(args)
        print(json.dumps({"source_fidelity": aggregate}, sort_keys=True), flush=True)
    if args.stage in ("merge", "all"):
        _validate_gate(args)
        _merge(args)
        print(json.dumps({"stage": "merge", "status": "complete"}), flush=True)


if __name__ == "__main__":
    main()
