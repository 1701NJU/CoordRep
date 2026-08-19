#!/usr/bin/env python3
"""Verify and merge all shards of the frozen CSD CoordRep-Record audit."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import string
import sys
from collections import Counter
from pathlib import Path
from typing import Any, Dict, Iterable, Iterator, List


BASE = Path(__file__).resolve().parent.parent
if str(BASE) not in sys.path:
    sys.path.insert(0, str(BASE))

from coordrep.audit.records import ENTRY_RESOLUTIONS, record_digest


FAILED_STATUSES = {"FAILED_INPUT", "FAILED_PIPELINE"}
OUTPUT_SCOPES = {"mono_eta1", "mono_haptic", "multi_eta1", "multi_haptic"}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _atomic_json(path: Path, payload: Dict[str, Any]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    os.replace(temporary, path)


def _copy_concat(paths: Iterable[Path], destination: Path) -> None:
    temporary = destination.with_suffix(destination.suffix + ".partial")
    with temporary.open("wb") as output:
        for path in paths:
            with path.open("rb") as source:
                for block in iter(lambda: source.read(8 * 1024 * 1024), b""):
                    output.write(block)
        output.flush()
        os.fsync(output.fileno())
    os.replace(temporary, destination)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--num-shards", type=int, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--expected-database-entries", type=int, required=True)
    parser.add_argument("--expected-tm-entries", type=int, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    protected_outputs = (
        args.output_dir / "RECORDS_INTERNAL.jsonl",
        args.output_dir / "ENTRY_OUTCOMES_INTERNAL.jsonl",
        args.output_dir / "FAILURES_INTERNAL.jsonl",
        args.output_dir / "SUMMARY_PUBLIC_AGGREGATE.json",
        args.output_dir / "SUMMARY_INTERNAL.json",
    )
    if any(path.exists() for path in protected_outputs):
        raise SystemExit("merged output already exists; choose a new output directory")
    shard_dirs = [
        args.root / f"shard_{index:04d}_of_{args.num_shards:04d}"
        for index in range(args.num_shards)
    ]
    summaries = []
    for directory in shard_dirs:
        summary_path = directory / "SUMMARY.json"
        if not summary_path.exists():
            raise SystemExit(f"missing shard summary: {summary_path}")
        summaries.append(json.loads(summary_path.read_text(encoding="utf-8")))

    reference_config = summaries[0]["config"]
    reference_database = reference_config["database"]
    reference_sources = reference_config["sources"]
    reference_release = reference_config["source_release"]
    reference_protocol = reference_config["protocol_id"]
    reference_runtime = reference_config.get("runtime")
    reference_database_sha256 = reference_config.get("database_sha256")
    if (
        not isinstance(reference_database_sha256, str)
        or len(reference_database_sha256) != 64
        or any(character not in string.hexdigits for character in reference_database_sha256)
    ):
        raise SystemExit("a complete frozen-database SHA-256 is required")
    aggregate = Counter()
    previous_stop = 0
    for shard_id, summary in enumerate(summaries):
        config = summary["config"]
        expected_start = args.expected_database_entries * shard_id // args.num_shards
        expected_stop = args.expected_database_entries * (shard_id + 1) // args.num_shards
        if config["database"] != reference_database:
            raise SystemExit(f"database fingerprint mismatch in shard {shard_id}")
        if config["sources"] != reference_sources:
            raise SystemExit(f"source fingerprint mismatch in shard {shard_id}")
        if config["source_release"] != reference_release:
            raise SystemExit(f"source-release mismatch in shard {shard_id}")
        if config["protocol_id"] != reference_protocol:
            raise SystemExit(f"protocol mismatch in shard {shard_id}")
        if config.get("runtime") != reference_runtime:
            raise SystemExit(f"runtime fingerprint mismatch in shard {shard_id}")
        if config.get("database_sha256") != reference_database_sha256:
            raise SystemExit(f"database SHA-256 mismatch in shard {shard_id}")
        if config["shard_id"] != shard_id or config["num_shards"] != args.num_shards:
            raise SystemExit(f"shard metadata mismatch in shard {shard_id}")
        if config["work_kind"] != "database_range":
            raise SystemExit(f"shard {shard_id} is not a full database range")
        if config["start"] != expected_start or config["stop"] != expected_stop:
            raise SystemExit(f"range mismatch in shard {shard_id}")
        if expected_start != previous_stop:
            raise SystemExit(f"gap or overlap before shard {shard_id}")
        if summary["processed_to"] != expected_stop:
            raise SystemExit(f"shard {shard_id} did not reach its stop")
        attempted = int(summary["counters"].get("inputs_attempted", 0))
        if attempted != expected_stop - expected_start:
            raise SystemExit(f"input accounting mismatch in shard {shard_id}")
        previous_stop = expected_stop
        aggregate.update(summary["counters"])

    if previous_stop != args.expected_database_entries:
        raise SystemExit("shard ranges do not cover the complete database")
    if int(reference_database["entry_count"]) != args.expected_database_entries:
        raise SystemExit("database entry count does not match the frozen release")
    if aggregate["inputs_attempted"] != args.expected_database_entries:
        raise SystemExit("aggregate input count does not match the frozen release")
    if aggregate["tm_entries"] != args.expected_tm_entries:
        raise SystemExit(
            f"transition-metal count is {aggregate['tm_entries']}, expected {args.expected_tm_entries}"
        )
    if aggregate["entry_read_failed"] or aggregate["classification_failed"]:
        raise SystemExit("one or more database entries could not be classified")
    if (
        aggregate["non_target_entries"] + aggregate["tm_entries"]
        != aggregate["inputs_attempted"]
    ):
        raise SystemExit("database input classification counters do not close")
    if (
        aggregate["records_emitted"]
        + aggregate["controlled_record_failure"]
        + aggregate["untyped_pipeline_failure"]
        != aggregate["tm_entries"]
    ):
        raise SystemExit("transition-metal terminal counters do not close")
    if sum(aggregate[status] for status in ENTRY_RESOLUTIONS) != aggregate[
        "records_emitted"
    ]:
        raise SystemExit("emitted record-resolution counters do not close")
    if (
        aggregate["structural_records_emitted"]
        + aggregate["audit_only_records_emitted"]
        != aggregate["records_emitted"]
    ):
        raise SystemExit("structural and audit-only record counters do not close")

    record_paths = [directory / "RECORDS_INTERNAL.jsonl" for directory in shard_dirs]
    outcome_paths = [directory / "ENTRY_OUTCOMES.jsonl" for directory in shard_dirs]
    failure_paths = [directory / "FAILURES_INTERNAL.jsonl" for directory in shard_dirs]
    for path in record_paths + outcome_paths + failure_paths:
        if not path.exists():
            raise SystemExit(f"missing shard payload: {path}")

    # Before concatenation, prove that every terminal outcome is paired with
    # exactly one matching emitted record or failure row in its own range.
    for shard_id, (record_path, outcome_path, failure_path) in enumerate(
        zip(record_paths, outcome_paths, failure_paths)
    ):
        expected_start = args.expected_database_entries * shard_id // args.num_shards
        expected_stop = args.expected_database_entries * (shard_id + 1) // args.num_shards
        with (
            record_path.open(encoding="utf-8") as record_handle,
            outcome_path.open(encoding="utf-8") as outcome_handle,
            failure_path.open(encoding="utf-8") as failure_handle,
        ):
            records = (json.loads(line) for line in record_handle if line.strip())
            failures = (json.loads(line) for line in failure_handle if line.strip())
            current_record = next(records, None)
            current_failure = next(failures, None)
            for line in outcome_handle:
                if not line.strip():
                    continue
                outcome = json.loads(line)
                source_index = int(outcome["source_index"])
                if not expected_start <= source_index < expected_stop:
                    raise SystemExit(
                        f"outcome source index outside shard {shard_id}: {source_index}"
                    )
                status = str(outcome["status"])
                if status.startswith("EMITTED_"):
                    if status not in ENTRY_RESOLUTIONS:
                        raise SystemExit(f"invalid emitted status in shard {shard_id}: {status}")
                    if outcome.get("schema_valid") is not True:
                        raise SystemExit(f"non-valid schema outcome in shard {shard_id}")
                    if current_record is None:
                        raise SystemExit(f"missing record paired to an outcome in shard {shard_id}")
                    keys = ("source_index", "refcode", "record_checksum_sha256")
                    if any(current_record.get(key) != outcome.get(key) for key in keys):
                        raise SystemExit(f"record/outcome mismatch in shard {shard_id}")
                    if current_record.get("entry_resolution") != status:
                        raise SystemExit(f"record status mismatch in shard {shard_id}")
                    if current_record.get("scope") not in OUTPUT_SCOPES:
                        raise SystemExit(f"invalid output scope in shard {shard_id}")
                    if int(current_record.get("expected_metal_sites", -1)) != len(
                        current_record.get("metal_sites", [])
                    ):
                        raise SystemExit(f"metal-site accounting mismatch in shard {shard_id}")
                    structural_body = {
                        key: current_record[key]
                        for key in (
                            "schema_version",
                            "scope",
                            "scope_flags",
                            "metal_sites",
                            "donor_groups",
                            "incidences",
                            "metal_relations",
                            "external_center_relations",
                            "external_donor_relations",
                            "donor_group_relations",
                            "donor_pair_relations",
                        )
                    }
                    if record_digest(structural_body) != current_record.get(
                        "record_checksum_sha256"
                    ):
                        raise SystemExit(f"record digest mismatch in shard {shard_id}")
                    sites = current_record["metal_sites"]
                    structural_sites = sum(
                        site.get("record_level") != "audit_only" for site in sites
                    )
                    if (status == "EMITTED_AUDIT_ONLY") != (structural_sites == 0):
                        raise SystemExit(
                            f"audit-only status/content mismatch in shard {shard_id}"
                        )
                    site_resolutions = [str(site.get("resolution", "")) for site in sites]
                    if status == "EMITTED_RESOLVED" and any(
                        resolution != "resolved" for resolution in site_resolutions
                    ):
                        raise SystemExit(
                            f"resolved status contains a non-resolved site in shard {shard_id}"
                        )
                    if status == "EMITTED_PARTIAL" and not any(
                        resolution in {"partial", "unresolved"}
                        for resolution in site_resolutions
                    ):
                        raise SystemExit(
                            f"partial status lacks a partial site in shard {shard_id}"
                        )
                    if status == "EMITTED_AMBIGUOUS" and not any(
                        resolution == "ambiguous" for resolution in site_resolutions
                    ):
                        raise SystemExit(
                            f"ambiguous status lacks an ambiguous site in shard {shard_id}"
                        )
                    expected_fields = {
                        "scope": current_record["scope"],
                        "scope_flags": current_record["scope_flags"],
                        "metal_sites": len(sites),
                        "donor_groups": len(current_record["donor_groups"]),
                        "incidences": len(current_record["incidences"]),
                        "structural_record": structural_sites > 0,
                        "structural_metal_sites": structural_sites,
                        "audit_only_metal_sites": len(sites) - structural_sites,
                        "all_metal_sites_structural": structural_sites == len(sites),
                    }
                    if any(outcome.get(key) != value for key, value in expected_fields.items()):
                        raise SystemExit(f"record/outcome statistics mismatch in shard {shard_id}")
                    expected_issues = sorted({
                        *current_record.get("issue_codes", []),
                        *(
                            issue
                            for site in current_record.get("metal_sites", [])
                            for issue in site.get("issues", [])
                        ),
                        *(
                            issue
                            for group in current_record.get("donor_groups", [])
                            for issue in group.get("issues", [])
                        ),
                    })
                    if outcome.get("issue_codes") != expected_issues:
                        raise SystemExit(f"record/outcome issue ledger mismatch in shard {shard_id}")
                    current_record = next(records, None)
                elif status.startswith("FAILED_"):
                    if status not in FAILED_STATUSES:
                        raise SystemExit(f"invalid failure status in shard {shard_id}: {status}")
                    if current_failure is None or current_failure != outcome:
                        raise SystemExit(f"failure/outcome mismatch in shard {shard_id}")
                    failed_sites = int(outcome.get("failed_metal_sites", 0))
                    metal_sites = int(outcome.get("metal_sites", 0))
                    if (
                        metal_sites < 1
                        or failed_sites != metal_sites
                        or int(outcome.get("structural_metal_sites", 0)) != 0
                        or int(outcome.get("audit_only_metal_sites", 0)) != 0
                    ):
                        raise SystemExit(
                            f"failed-outcome site fields do not close in shard {shard_id}"
                        )
                    current_failure = next(failures, None)
                else:
                    raise SystemExit(f"unknown terminal status in shard {shard_id}: {status}")
            if current_record is not None or next(records, None) is not None:
                raise SystemExit(f"unpaired record row(s) in shard {shard_id}")
            if current_failure is not None or next(failures, None) is not None:
                raise SystemExit(f"unpaired failure row(s) in shard {shard_id}")

    records_out = args.output_dir / "RECORDS_INTERNAL.jsonl"
    outcomes_out = args.output_dir / "ENTRY_OUTCOMES_INTERNAL.jsonl"
    failures_out = args.output_dir / "FAILURES_INTERNAL.jsonl"
    _copy_concat(record_paths, records_out)
    _copy_concat(outcome_paths, outcomes_out)
    _copy_concat(failure_paths, failures_out)

    outcome_counts = Counter()
    source_indices = set()
    refcodes = set()
    outcome_rows = 0
    with outcomes_out.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            row = json.loads(line)
            outcome_rows += 1
            source_index = int(row["source_index"])
            refcode = str(row["refcode"])
            if source_index in source_indices:
                raise SystemExit(f"duplicate source index: {source_index}")
            if refcode in refcodes:
                raise SystemExit(f"duplicate refcode: {refcode}")
            source_indices.add(source_index)
            refcodes.add(refcode)
            outcome_counts[f"status:{row['status']}"] += 1
            outcome_counts[f"scope:{row.get('scope', 'failed')}"] += 1
            outcome_counts["structural_records"] += int(bool(row.get("structural_record")))
            outcome_counts["schema_valid"] += int(bool(row.get("schema_valid")))
            outcome_counts["metal_sites"] += int(row.get("metal_sites", 0))
            outcome_counts["structural_metal_sites"] += int(
                row.get("structural_metal_sites", 0)
            )
            outcome_counts["audit_only_metal_sites"] += int(
                row.get("audit_only_metal_sites", 0)
            )
            outcome_counts["failed_metal_sites"] += int(row.get("failed_metal_sites", 0))
            outcome_counts["entries_all_metal_sites_structural"] += int(
                bool(row.get("all_metal_sites_structural"))
            )
            outcome_counts["incidences"] += int(row.get("incidences", 0))
            for flag in row.get("scope_flags", []):
                outcome_counts[f"flag:{flag}"] += 1
            for issue in row.get("issue_codes", []):
                outcome_counts[f"issue:{issue}"] += 1
            if str(row.get("status", "")).startswith("FAILED_"):
                outcome_counts[f"failure_issue:{row.get('issue_code', 'UNKNOWN')}"] += 1

    if outcome_rows != args.expected_tm_entries:
        raise SystemExit(f"outcome rows are {outcome_rows}, expected {args.expected_tm_entries}")
    record_rows = sum(1 for line in records_out.open(encoding="utf-8") if line.strip())
    failure_rows = sum(1 for line in failures_out.open(encoding="utf-8") if line.strip())
    if record_rows != aggregate["records_emitted"]:
        raise SystemExit("merged record row count disagrees with shard counters")
    failed_outcomes = sum(
        count for key, count in outcome_counts.items() if key.startswith("status:FAILED_")
    )
    if record_rows + failed_outcomes != outcome_rows:
        raise SystemExit("emitted and failed terminal outcomes do not close")
    if failure_rows != failed_outcomes:
        raise SystemExit("failure rows and failed terminal outcomes do not close")
    if outcome_counts["schema_valid"] != record_rows:
        raise SystemExit("schema-valid outcome count and emitted records do not close")
    if (
        outcome_counts["structural_metal_sites"]
        + outcome_counts["audit_only_metal_sites"]
        + outcome_counts["failed_metal_sites"]
        != outcome_counts["metal_sites"]
    ):
        raise SystemExit("merged metal-site fields do not close")
    if (
        outcome_counts["structural_records"]
        + outcome_counts["status:EMITTED_AUDIT_ONLY"]
        != record_rows
    ):
        raise SystemExit("structural and audit-only emitted outcomes do not close")

    public_summary = {
        "protocol_id": reference_config["protocol_id"],
        "source_release": reference_config["source_release"],
        "database_entries_evaluated": args.expected_database_entries,
        "transition_metal_entries": args.expected_tm_entries,
        "entry_outcomes": outcome_rows,
        "records_emitted": record_rows,
        "record_emission_rate": record_rows / args.expected_tm_entries,
        "structural_records_emitted": outcome_counts["structural_records"],
        "structural_record_coverage": (
            outcome_counts["structural_records"] / args.expected_tm_entries
        ),
        "entries_all_metal_sites_structural": outcome_counts[
            "entries_all_metal_sites_structural"
        ],
        "all_metal_sites_structural_entry_coverage": (
            outcome_counts["entries_all_metal_sites_structural"]
            / args.expected_tm_entries
        ),
        "schema_valid_records": outcome_counts["schema_valid"],
        "audit_only_records": outcome_counts["status:EMITTED_AUDIT_ONLY"],
        "failed_outcomes": sum(
            count for key, count in outcome_counts.items() if key.startswith("status:FAILED_")
        ),
        "metal_sites_accounted": outcome_counts["metal_sites"],
        "structural_metal_sites": outcome_counts["structural_metal_sites"],
        "audit_only_metal_sites": outcome_counts["audit_only_metal_sites"],
        "failed_metal_sites": outcome_counts["failed_metal_sites"],
        "structural_metal_site_coverage": (
            outcome_counts["structural_metal_sites"] / outcome_counts["metal_sites"]
            if outcome_counts["metal_sites"]
            else None
        ),
        "metal_site_accounting_closed": (
            outcome_counts["metal_sites"]
            == outcome_counts["structural_metal_sites"]
            + outcome_counts["audit_only_metal_sites"]
            + outcome_counts["failed_metal_sites"]
        ),
        "coordination_incidences": outcome_counts["incidences"],
        "status_counts": {
            key.removeprefix("status:"): value
            for key, value in sorted(outcome_counts.items())
            if key.startswith("status:")
        },
        "scope_counts": {
            key.removeprefix("scope:"): value
            for key, value in sorted(outcome_counts.items())
            if key.startswith("scope:")
        },
        "scope_flag_counts": {
            key.removeprefix("flag:"): value
            for key, value in sorted(outcome_counts.items())
            if key.startswith("flag:")
        },
        "issue_counts": {
            key.removeprefix("issue:"): value
            for key, value in sorted(outcome_counts.items())
            if key.startswith("issue:")
        },
        "failure_issue_counts": {
            key.removeprefix("failure_issue:"): value
            for key, value in sorted(outcome_counts.items())
            if key.startswith("failure_issue:")
        },
        "failure_rows": failure_rows,
        "artifacts": {
            "records_internal_sha256": _sha256(records_out),
            "outcomes_internal_sha256": _sha256(outcomes_out),
            "failures_internal_sha256": _sha256(failures_out),
        },
        "database_fingerprint": reference_database,
        "database_sha256": reference_database_sha256,
        "runtime_fingerprint": reference_runtime,
        "source_fingerprints": reference_sources,
        "shards": args.num_shards,
    }
    _atomic_json(args.output_dir / "SUMMARY_PUBLIC_AGGREGATE.json", public_summary)
    internal_summary = {
        **public_summary,
        "aggregate_shard_counters": dict(aggregate),
        "outcome_counters": dict(outcome_counts),
        "shard_summaries": summaries,
    }
    _atomic_json(args.output_dir / "SUMMARY_INTERNAL.json", internal_summary)
    print(json.dumps(public_summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
