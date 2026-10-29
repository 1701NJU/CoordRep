#!/usr/bin/env python3
"""Run one restartable shard of the frozen CSD CoordRep-Record audit.

The runner emits one terminal outcome for every in-domain metal-containing 3D
entry encountered.  Schema-valid structural records and controlled failures
are written to separate streams so a failure envelope can never inflate the
record-coverage numerator.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import sys
import time
from collections import Counter
from pathlib import Path
from typing import Any, Dict, Iterable, Iterator, List, Optional, Tuple


BASE = Path(__file__).resolve().parent.parent
if str(BASE) not in sys.path:
    sys.path.insert(0, str(BASE))

from ccdc.io import EntryReader
import ccdc
import numpy as np
import scipy

from coordrep.audit.csd_adapter import CSDRecordError, PROTOCOL_ID, adapt_csd_entry
from coordrep.audit.metal_policy import (
    METAL_POLICY_ID,
    METAL_POLICY_SHA256,
    is_target_metal_symbol,
    metal_blocks,
)
from coordrep.audit.records import canonical_json


def _atomic_json(path: Path, payload: Dict[str, Any]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    os.replace(temporary, path)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _database_fingerprint(reader: Any, database_arg: str) -> Dict[str, Any]:
    path_text = str(getattr(reader, "file_name", database_arg))
    path = Path(path_text)
    payload: Dict[str, Any] = {
        "argument": database_arg,
        "resolved_file": path_text,
        "entry_count": len(reader),
    }
    if path.exists():
        stat = path.stat()
        payload.update({"file_size": stat.st_size, "mtime_ns": stat.st_mtime_ns})
    return payload


def _source_fingerprint() -> Dict[str, str]:
    paths = [
        BASE / "coordrep" / "audit" / "records.py",
        BASE / "coordrep" / "audit" / "csd_adapter.py",
        BASE / "coordrep" / "audit" / "metal_policy.py",
        BASE / "coordrep" / "geometry" / "shape.py",
        Path(__file__).resolve(),
    ]
    return {str(path.relative_to(BASE)): _sha256(path) for path in paths}


def _runtime_fingerprint() -> Dict[str, str]:
    return {
        "python": sys.version.replace("\n", " "),
        "implementation": platform.python_implementation(),
        "platform": platform.platform(),
        "ccdc": str(getattr(ccdc, "__version__", "unknown")),
        "numpy": str(np.__version__),
        "scipy": str(scipy.__version__),
    }


def _resolve_range(total: int, shard_id: int, num_shards: int) -> Tuple[int, int]:
    if num_shards < 1 or not 0 <= shard_id < num_shards:
        raise SystemExit("require num_shards >= 1 and 0 <= shard_id < num_shards")
    return total * shard_id // num_shards, total * (shard_id + 1) // num_shards


def _load_manifest(path: Optional[Path]) -> Optional[List[Dict[str, Any]]]:
    if path is None:
        return None
    rows: List[Dict[str, Any]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("{"):
            payload = json.loads(stripped)
            refcode = payload.get("refcode") or payload.get("identifier")
            if refcode:
                rows.append({
                    "refcode": str(refcode),
                    "source_index": int(payload.get("source_index", payload.get("database_index", -1))),
                })
        else:
            rows.append({"refcode": stripped.split(",", 1)[0], "source_index": -1})
    return rows


def _metal_symbols(entry: Any) -> Tuple[str, ...]:
    """Return in-domain metal symbols for a usable three-dimensional entry."""

    if not bool(getattr(entry, "has_3d_structure", False)):
        return ()
    molecule = entry.molecule
    if molecule is None:
        return ()
    return tuple(
        str(getattr(atom, "atomic_symbol", "?"))
        for atom in molecule.atoms
        if is_target_metal_symbol(getattr(atom, "atomic_symbol", "?"))
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--database", default="CSD")
    parser.add_argument("--database-sha256")
    parser.add_argument("--source-release", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--shard-id", type=int, default=0)
    parser.add_argument("--num-shards", type=int, default=1)
    parser.add_argument("--refcodes-file", type=Path)
    parser.add_argument("--checkpoint-every", type=int, default=1000)
    parser.add_argument("--progress-every", type=int, default=10000)
    parser.add_argument("--max-targets", type=int)
    parser.add_argument("--resume", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.checkpoint_every <= 0 or args.progress_every <= 0:
        raise SystemExit("checkpoint-every and progress-every must be positive")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    records_partial = args.output_dir / "RECORDS_INTERNAL.jsonl.partial"
    outcomes_partial = args.output_dir / "ENTRY_OUTCOMES.jsonl.partial"
    failures_partial = args.output_dir / "FAILURES_INTERNAL.jsonl.partial"
    records_final = args.output_dir / "RECORDS_INTERNAL.jsonl"
    outcomes_final = args.output_dir / "ENTRY_OUTCOMES.jsonl"
    failures_final = args.output_dir / "FAILURES_INTERNAL.jsonl"
    checkpoint_path = args.output_dir / "checkpoint.json"
    summary_path = args.output_dir / "SUMMARY.json"

    if any(path.exists() for path in (records_final, outcomes_final, failures_final, summary_path)):
        raise SystemExit("completed output already exists; choose a new output directory")

    reader = EntryReader(args.database)
    database_fingerprint = _database_fingerprint(reader, args.database)
    source_fingerprint = _source_fingerprint()
    manifest = _load_manifest(args.refcodes_file)
    if manifest is None:
        start, stop = _resolve_range(len(reader), args.shard_id, args.num_shards)
        work_total = stop - start
        work_kind = "database_range"
    else:
        start, stop = 0, len(manifest)
        work_total = len(manifest)
        work_kind = "refcode_manifest"

    config = {
        "protocol_id": PROTOCOL_ID,
        "source_release": args.source_release,
        "database": database_fingerprint,
        "database_sha256": args.database_sha256,
        "runtime": _runtime_fingerprint(),
        "sources": source_fingerprint,
        "metal_policy_id": METAL_POLICY_ID,
        "metal_policy_sha256": METAL_POLICY_SHA256,
        "work_kind": work_kind,
        "start": start,
        "stop": stop,
        "shard_id": args.shard_id,
        "num_shards": args.num_shards,
        "refcodes_file": str(args.refcodes_file.resolve()) if args.refcodes_file else None,
        "refcodes_sha256": _sha256(args.refcodes_file) if args.refcodes_file else None,
        "max_targets": args.max_targets,
    }

    counters: Counter = Counter()
    next_position = start
    offsets = {"records": 0, "outcomes": 0, "failures": 0}
    started_at = time.time()
    if args.resume:
        payload = json.loads(checkpoint_path.read_text(encoding="utf-8"))
        if payload["config"] != config:
            raise SystemExit("checkpoint configuration does not match this run")
        next_position = int(payload["next_position"])
        offsets = {key: int(value) for key, value in payload["offsets"].items()}
        counters.update(payload["counters"])
        started_at = float(payload["started_at_epoch"])
        for path, key in (
            (records_partial, "records"),
            (outcomes_partial, "outcomes"),
            (failures_partial, "failures"),
        ):
            with path.open("r+b") as handle:
                handle.truncate(offsets[key])
    else:
        if checkpoint_path.exists() or any(
            path.exists() for path in (records_partial, outcomes_partial, failures_partial)
        ):
            raise SystemExit("partial output exists; use --resume or a new output directory")

    def write_checkpoint(position: int, handles: Dict[str, Any]) -> None:
        for handle in handles.values():
            handle.flush()
            os.fsync(handle.fileno())
        current_offsets = {key: handle.tell() for key, handle in handles.items()}
        _atomic_json(checkpoint_path, {
            "config": config,
            "next_position": position,
            "offsets": current_offsets,
            "counters": dict(counters),
            "started_at_epoch": started_at,
            "checkpointed_at_epoch": time.time(),
        })

    modes = "ab" if args.resume else "wb"
    handles = {
        "records": records_partial.open(modes),
        "outcomes": outcomes_partial.open(modes),
        "failures": failures_partial.open(modes),
    }
    last_checkpoint = next_position
    if not args.resume:
        write_checkpoint(next_position, handles)

    def checkpoint_and_progress(position: int) -> None:
        nonlocal last_checkpoint
        if position - last_checkpoint >= args.checkpoint_every:
            write_checkpoint(position, handles)
            last_checkpoint = position
        if counters["inputs_attempted"] % args.progress_every == 0:
            elapsed = time.time() - started_at
            print(json.dumps({
                "position": position,
                "attempted": counters["inputs_attempted"],
                "target_entries": counters["target_entries"],
                "records": counters["records_emitted"],
                "structural": counters["structural_records_emitted"],
                "elapsed_s": round(elapsed, 2),
            }), flush=True)

    try:
        for position in range(next_position, stop):
            if args.max_targets is not None and counters["target_entries"] >= args.max_targets:
                break
            counters["inputs_attempted"] += 1
            source_index = position if manifest is None else int(manifest[position]["source_index"])
            source_key = position if manifest is None else str(manifest[position]["refcode"])
            try:
                entry = reader[source_key] if manifest is None else reader.entry(source_key)
            except Exception as exc:
                counters["entry_read_failed"] += 1
                failure = {
                    "source_index": source_index,
                    "refcode": str(source_key) if manifest is not None else "",
                    "status": "FAILED_INPUT",
                    "issue_code": "ENTRY_READ_FAILED",
                    "exception_type": type(exc).__name__,
                    "exception_message": str(exc)[:500],
                }
                handles["failures"].write((canonical_json(failure) + "\n").encode("utf-8"))
                next_position = position + 1
                checkpoint_and_progress(next_position)
                continue

            try:
                target_symbols = _metal_symbols(entry)
                n_metals = len(target_symbols)
            except Exception as exc:
                counters["classification_failed"] += 1
                failure = {
                    "source_index": source_index,
                    "refcode": str(getattr(entry, "identifier", source_key)),
                    "status": "FAILED_INPUT",
                    "issue_code": "CLASSIFICATION_FAILED",
                    "exception_type": type(exc).__name__,
                    "exception_message": str(exc)[:500],
                }
                handles["failures"].write((canonical_json(failure) + "\n").encode("utf-8"))
                next_position = position + 1
                checkpoint_and_progress(next_position)
                continue
            if n_metals == 0:
                counters["non_target_entries"] += 1
                next_position = position + 1
            else:
                counters["target_entries"] += 1
                unique_elements = tuple(sorted(set(target_symbols)))
                unique_blocks = metal_blocks(unique_elements)
                for element in unique_elements:
                    counters[f"element_presence:{element}"] += 1
                for block in unique_blocks:
                    counters[f"metal_block_presence:{block}"] += 1
                if len(unique_blocks) > 1:
                    counters["mixed_metal_block_entries"] += 1
                try:
                    record = adapt_csd_entry(
                        entry,
                        source_index,
                        source_release=args.source_release,
                    )
                    record_payload = record.to_dict()
                    handles["records"].write((canonical_json(record_payload) + "\n").encode("utf-8"))
                    structural_metal_sites = sum(
                        site.record_level != "audit_only" for site in record.metal_sites
                    )
                    audit_only_metal_sites = record.expected_metal_sites - structural_metal_sites
                    structural = structural_metal_sites > 0
                    aggregated_issues = sorted({
                        *record.issue_codes,
                        *(issue for site in record.metal_sites for issue in site.issues),
                        *(issue for group in record.donor_groups for issue in group.issues),
                    })
                    outcome = {
                        "source_index": source_index,
                        "refcode": record.refcode,
                        "status": record.entry_resolution,
                        "scope": record.scope,
                        "scope_flags": record.scope_flags,
                        "metal_elements": unique_elements,
                        "metal_blocks": unique_blocks,
                        "metal_sites": record.expected_metal_sites,
                        "structural_metal_sites": structural_metal_sites,
                        "audit_only_metal_sites": audit_only_metal_sites,
                        "all_metal_sites_structural": audit_only_metal_sites == 0,
                        "donor_groups": len(record.donor_groups),
                        "incidences": len(record.incidences),
                        "structural_record": structural,
                        "schema_valid": True,
                        "record_checksum_sha256": record.record_checksum_sha256,
                        "issue_codes": aggregated_issues,
                    }
                    handles["outcomes"].write((canonical_json(outcome) + "\n").encode("utf-8"))
                    counters[record.entry_resolution] += 1
                    counters[f"scope:{record.scope}"] += 1
                    counters["records_emitted"] += 1
                    if structural:
                        counters["structural_records_emitted"] += 1
                    else:
                        counters["audit_only_records_emitted"] += 1
                    counters["metal_sites_accounted"] += record.expected_metal_sites
                except CSDRecordError as exc:
                    counters["controlled_record_failure"] += 1
                    failure = {
                        "source_index": source_index,
                        "refcode": str(getattr(entry, "identifier", source_key)),
                        "status": "FAILED_INPUT",
                        "issue_code": str(exc),
                        "exception_type": type(exc).__name__,
                        "exception_message": str(exc)[:500],
                        "metal_sites": n_metals,
                        "metal_elements": unique_elements,
                        "metal_blocks": unique_blocks,
                        "structural_metal_sites": 0,
                        "audit_only_metal_sites": 0,
                        "failed_metal_sites": n_metals,
                    }
                    handles["failures"].write((canonical_json(failure) + "\n").encode("utf-8"))
                    handles["outcomes"].write((canonical_json(failure) + "\n").encode("utf-8"))
                except Exception as exc:
                    counters["untyped_pipeline_failure"] += 1
                    failure = {
                        "source_index": source_index,
                        "refcode": str(getattr(entry, "identifier", source_key)),
                        "status": "FAILED_PIPELINE",
                        "issue_code": "UNTYPED_PIPELINE_EXCEPTION",
                        "exception_type": type(exc).__name__,
                        "exception_message": str(exc)[:500],
                        "metal_sites": n_metals,
                        "metal_elements": unique_elements,
                        "metal_blocks": unique_blocks,
                        "structural_metal_sites": 0,
                        "audit_only_metal_sites": 0,
                        "failed_metal_sites": n_metals,
                    }
                    handles["failures"].write((canonical_json(failure) + "\n").encode("utf-8"))
                    handles["outcomes"].write((canonical_json(failure) + "\n").encode("utf-8"))
                next_position = position + 1

            checkpoint_and_progress(next_position)
    except BaseException:
        write_checkpoint(next_position, handles)
        raise
    finally:
        for handle in handles.values():
            handle.close()

    elapsed = time.time() - started_at
    processed_to = next_position
    completed_range = processed_to >= stop or (
        args.max_targets is not None and counters["target_entries"] >= args.max_targets
    )
    if not completed_range:
        raise SystemExit("run stopped before its configured terminal condition")

    os.replace(records_partial, records_final)
    os.replace(outcomes_partial, outcomes_final)
    os.replace(failures_partial, failures_final)
    if checkpoint_path.exists():
        checkpoint_path.unlink()
    summary = {
        "protocol_id": PROTOCOL_ID,
        "config": config,
        "counters": dict(counters),
        "processed_to": processed_to,
        "elapsed_seconds": round(elapsed, 3),
        "records_file": records_final.name,
        "outcomes_file": outcomes_final.name,
        "failures_file": failures_final.name,
    }
    _atomic_json(summary_path, summary)
    print(json.dumps(summary, indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
