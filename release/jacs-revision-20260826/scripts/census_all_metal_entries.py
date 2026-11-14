#!/usr/bin/env python3
"""Count the frozen metal domain independently of structural-record support.

The structural-record audit targets three-dimensional metal-containing CSD
entries.  This companion census also counts metal-containing entries without
usable 3D structure, so the publication can distinguish the complete
metal-containing census from the operational structural target.  It emits
aggregate evidence only and never redistributes CSD-derived rows.
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
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, Iterable, Sequence, Tuple

import ccdc
from ccdc.io import EntryReader


BASE = Path(__file__).resolve().parent.parent
if str(BASE) not in sys.path:
    sys.path.insert(0, str(BASE))

from coordrep.audit.metal_policy import (
    ALL_METALS,
    METAL_POLICY_ID,
    METAL_POLICY_SHA256,
    is_target_metal_symbol,
    metal_blocks,
)


PROTOCOL_ID = "csd-all-metal-entry-census-20260825-v1"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--database", required=True)
    parser.add_argument("--database-sha256", required=True)
    parser.add_argument("--source-release", required=True)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--expected-database-entries", type=int, required=True)
    parser.add_argument("--expected-3d-target-entries", type=int)
    parser.add_argument("--num-shards", type=int, default=48)
    parser.add_argument("--max-workers", type=int, default=24)
    parser.add_argument("--progress-every", type=int, default=10000)
    parser.add_argument("--shard-id", type=int, help=argparse.SUPPRESS)
    return parser.parse_args()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _atomic_json(path: Path, payload: Dict[str, Any]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8"
    )
    os.replace(temporary, path)


def _runtime() -> Dict[str, str]:
    return {
        "python": sys.version.replace("\n", " "),
        "implementation": platform.python_implementation(),
        "platform": platform.platform(),
        "ccdc": str(getattr(ccdc, "__version__", "unknown")),
    }


def _range(total: int, shard_id: int, num_shards: int) -> Tuple[int, int]:
    return total * shard_id // num_shards, total * (shard_id + 1) // num_shards


def _symbols(entry: Any) -> Tuple[str, ...]:
    molecule = entry.molecule
    if molecule is None:
        raise RuntimeError("CSD entry has no readable molecular object")
    return tuple(
        str(getattr(atom, "atomic_symbol", "?"))
        for atom in molecule.atoms
        if is_target_metal_symbol(getattr(atom, "atomic_symbol", "?"))
    )


def _increment_presence(
    counters: Counter, symbols: Iterable[str], suffix: str
) -> None:
    unique_elements = tuple(sorted(set(symbols)))
    unique_blocks = metal_blocks(unique_elements)
    for element in unique_elements:
        counters[f"element_presence:{suffix}:{element}"] += 1
    for block in unique_blocks:
        counters[f"block_presence:{suffix}:{block}"] += 1
    if len(unique_blocks) > 1:
        counters[f"mixed_blocks:{suffix}"] += 1


def _shard_path(args: argparse.Namespace, shard_id: int) -> Path:
    return (
        args.run_root
        / "shards"
        / f"shard_{shard_id:04d}_of_{args.num_shards:04d}.json"
    )


def _worker(args: argparse.Namespace) -> None:
    assert args.shard_id is not None
    reader = EntryReader(args.database)
    if len(reader) != args.expected_database_entries:
        raise SystemExit("database entry count differs from the frozen denominator")
    start, stop = _range(len(reader), args.shard_id, args.num_shards)
    counters: Counter = Counter()
    failures = []
    started = time.time()
    for index in range(start, stop):
        counters["inputs_attempted"] += 1
        try:
            entry = reader[index]
            has_3d = bool(getattr(entry, "has_3d_structure", False))
            counters["entries_with_3d" if has_3d else "entries_without_3d"] += 1
            symbols = _symbols(entry)
        except Exception as exc:
            counters["classification_failed"] += 1
            if len(failures) < 20:
                failures.append({
                    "source_index": index,
                    "exception_type": type(exc).__name__,
                    "exception_message": str(exc)[:300],
                })
            continue
        if symbols:
            suffix = "metal_3d" if has_3d else "metal_no3d"
            counters["metal_entries"] += 1
            counters[f"{suffix}_entries"] += 1
            _increment_presence(counters, symbols, "metal_all")
            _increment_presence(counters, symbols, suffix)
        else:
            counters["nonmetal_entries"] += 1
        if counters["inputs_attempted"] % args.progress_every == 0:
            print(json.dumps({
                "shard": args.shard_id,
                "attempted": counters["inputs_attempted"],
                "metal_entries": counters["metal_entries"],
                "metal_3d_entries": counters["metal_3d_entries"],
            }, sort_keys=True), flush=True)
    payload = {
        "config": {
            "protocol_id": PROTOCOL_ID,
            "source_release": args.source_release,
            "database": args.database,
            "database_sha256": args.database_sha256.lower(),
            "database_entries": len(reader),
            "metal_policy_id": METAL_POLICY_ID,
            "metal_policy_sha256": METAL_POLICY_SHA256,
            "script_sha256": _sha256(Path(__file__).resolve()),
            "runtime": _runtime(),
            "shard_id": args.shard_id,
            "num_shards": args.num_shards,
            "start": start,
            "stop": stop,
        },
        "counters": dict(counters),
        "failure_examples": failures,
        "elapsed_seconds": round(time.time() - started, 3),
    }
    destination = _shard_path(args, args.shard_id)
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        raise SystemExit(f"refusing to overwrite completed shard: {destination}")
    _atomic_json(destination, payload)


def _run_logged(command: Sequence[str], path: Path) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    environment = os.environ.copy()
    environment.update({
        "OMP_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "NUMEXPR_NUM_THREADS": "1",
    })
    with path.open("ab") as handle:
        completed = subprocess.run(
            list(command),
            stdout=handle,
            stderr=subprocess.STDOUT,
            env=environment,
            check=False,
        )
    return completed.returncode


def _launch_one(args: argparse.Namespace, shard_id: int) -> Dict[str, Any]:
    destination = _shard_path(args, shard_id)
    if destination.exists():
        return {"shard": shard_id, "status": "skipped_complete"}
    command = [
        sys.executable,
        str(Path(__file__).resolve()),
        "--database", args.database,
        "--database-sha256", args.database_sha256.lower(),
        "--source-release", args.source_release,
        "--run-root", str(args.run_root),
        "--expected-database-entries", str(args.expected_database_entries),
        "--num-shards", str(args.num_shards),
        "--max-workers", str(args.max_workers),
        "--progress-every", str(args.progress_every),
        "--shard-id", str(shard_id),
    ]
    if args.expected_3d_target_entries is not None:
        command.extend((
            "--expected-3d-target-entries",
            str(args.expected_3d_target_entries),
        ))
    code = _run_logged(command, args.run_root / "logs" / f"shard_{shard_id:04d}.log")
    return {
        "shard": shard_id,
        "status": "complete" if code == 0 else "failed",
        "returncode": code,
    }


def _selected(counter: Counter, prefix: str) -> Dict[str, int]:
    return {
        key.removeprefix(prefix): value
        for key, value in sorted(counter.items())
        if key.startswith(prefix)
    }


def _parent(args: argparse.Namespace) -> None:
    if args.num_shards < 1 or args.max_workers < 1:
        raise SystemExit("num-shards and max-workers must be positive")
    if args.progress_every < 1:
        raise SystemExit("progress-every must be positive")
    database = Path(args.database)
    if not database.is_file():
        raise SystemExit("the frozen census requires an explicit database file")
    if _sha256(database).lower() != args.database_sha256.lower():
        raise SystemExit("database SHA-256 mismatch")
    args.run_root.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(max_workers=args.max_workers) as executor:
        futures = {
            executor.submit(_launch_one, args, shard_id): shard_id
            for shard_id in range(args.num_shards)
        }
        results = []
        for future in as_completed(futures):
            result = future.result()
            results.append(result)
            print(json.dumps(result, sort_keys=True), flush=True)
    failed = [row for row in results if row["status"] == "failed"]
    if failed:
        raise SystemExit(f"metal census failed for {len(failed)} shard(s): {failed}")

    aggregate: Counter = Counter()
    reference_config = None
    previous_stop = 0
    for shard_id in range(args.num_shards):
        payload = json.loads(_shard_path(args, shard_id).read_text(encoding="utf-8"))
        config = payload["config"]
        start, stop = _range(args.expected_database_entries, shard_id, args.num_shards)
        if config["start"] != start or config["stop"] != stop or start != previous_stop:
            raise SystemExit(f"range gap, overlap, or mismatch in shard {shard_id}")
        previous_stop = stop
        invariant = {
            key: value
            for key, value in config.items()
            if key not in {"shard_id", "start", "stop"}
        }
        if reference_config is None:
            reference_config = invariant
        elif invariant != reference_config:
            raise SystemExit(f"configuration mismatch in shard {shard_id}")
        aggregate.update(payload["counters"])
    if previous_stop != args.expected_database_entries:
        raise SystemExit("shards do not close to the database denominator")
    if aggregate["classification_failed"]:
        raise SystemExit("one or more CSD entries could not be classified")
    if aggregate["inputs_attempted"] != args.expected_database_entries:
        raise SystemExit("input accounting does not close")
    if (
        aggregate["entries_with_3d"] + aggregate["entries_without_3d"]
        != aggregate["inputs_attempted"]
    ):
        raise SystemExit("3D/non-3D partition does not close")
    if (
        aggregate["metal_entries"] + aggregate["nonmetal_entries"]
        != aggregate["inputs_attempted"]
    ):
        raise SystemExit("metal/nonmetal partition does not close")
    if (
        aggregate["metal_3d_entries"] + aggregate["metal_no3d_entries"]
        != aggregate["metal_entries"]
    ):
        raise SystemExit("metal 3D/non-3D partition does not close")
    if (
        args.expected_3d_target_entries is not None
        and aggregate["metal_3d_entries"] != args.expected_3d_target_entries
    ):
        raise SystemExit(
            "independent 3D metal census disagrees with the structural-audit target"
        )

    output = {
        "protocol_id": PROTOCOL_ID,
        "source_release": args.source_release,
        "database_entries": args.expected_database_entries,
        "database_sha256": args.database_sha256.lower(),
        "metal_policy_id": METAL_POLICY_ID,
        "metal_policy_sha256": METAL_POLICY_SHA256,
        "metal_policy_symbols": ALL_METALS,
        "metal_containing_entries": aggregate["metal_entries"],
        "metal_containing_3d_entries": aggregate["metal_3d_entries"],
        "metal_containing_without_3d_entries": aggregate["metal_no3d_entries"],
        "entries_with_3d": aggregate["entries_with_3d"],
        "entries_without_3d": aggregate["entries_without_3d"],
        "metal_block_entry_counts": _selected(
            aggregate, "block_presence:metal_all:"
        ),
        "metal_3d_block_entry_counts": _selected(
            aggregate, "block_presence:metal_3d:"
        ),
        "metal_no3d_block_entry_counts": _selected(
            aggregate, "block_presence:metal_no3d:"
        ),
        "mixed_metal_block_entries": aggregate["mixed_blocks:metal_all"],
        "mixed_metal_3d_entries": aggregate["mixed_blocks:metal_3d"],
        "mixed_metal_no3d_entries": aggregate["mixed_blocks:metal_no3d"],
        "metal_element_entry_counts": _selected(
            aggregate, "element_presence:metal_all:"
        ),
        "metal_3d_element_entry_counts": _selected(
            aggregate, "element_presence:metal_3d:"
        ),
        "metal_no3d_element_entry_counts": _selected(
            aggregate, "element_presence:metal_no3d:"
        ),
        "classification_failures": aggregate["classification_failed"],
        "shards": args.num_shards,
        "runtime_fingerprint": reference_config["runtime"],
        "script_sha256": reference_config["script_sha256"],
    }
    _atomic_json(args.run_root / "METAL_DOMAIN_CENSUS_PUBLIC.json", output)
    print(json.dumps(output, indent=2, sort_keys=True))


def main() -> None:
    args = parse_args()
    if args.shard_id is None:
        _parent(args)
    else:
        if not 0 <= args.shard_id < args.num_shards:
            raise SystemExit("shard-id is outside the configured range")
        _worker(args)


if __name__ == "__main__":
    main()
