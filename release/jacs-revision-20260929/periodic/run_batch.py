#!/usr/bin/env python3
"""Stream a selected prefix (or all) of the pinned public MOF ZIP for MID/SID.

Only aggregate JSON is printed. The optional --compare-frozen-entries reads a
private historical ledger for local QA and emits only pass/fail counts; it is
not required to run from the public ZIP.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import re
import sys
import time
import zipfile
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any, Dict, List, Mapping

from periodic_v3 import canonical_entry_multiset_id
from periodic_v3.adapter import extract_entry


ARCHIVE_SHA256 = "ABF9F2A7288FD83E2A59F19B717C202414E744A1B2C8C1F9EEC1F909554346FA"
EXPECTED_COLLECTION_ENTRIES = 15906
CANONICALIZER_SHA256 = "292AFD4DB89D1DD552F9A03433B207FD648BB30E9DED717BFB8DBCE17D4694AC"
QUOTIENT_SHA256 = "8C9DE81E1EBABDECE083BD6393786357BDC35CE3D6AD3142AAB205E2ECB19EEA"
_WORKER_ARCHIVE = None


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest().upper()


def truthy_marker(value: Any) -> bool:
    return str(value).strip() not in {"", "-", "0", "False", "false", "No", "no"}


def collection_index(archive: zipfile.ZipFile) -> List[Dict[str, Any]]:
    details = archive.read("CSD_MOF_Collection/Framework details.csv")
    rows = list(csv.DictReader(io.StringIO(details.decode("utf-8-sig"))))
    cif_lookup = {
        Path(name).stem.lower(): name
        for name in archive.namelist()
        if name.lower().endswith(".cif")
    }
    if len(rows) != EXPECTED_COLLECTION_ENTRIES or len(cif_lookup) != len(rows):
        raise ValueError("Framework details/CIF count differs from frozen collection")
    result = []
    seen = set()
    for index, row in enumerate(rows):
        refcode = row["CSD refcode"].strip().upper()
        if refcode in seen:
            raise ValueError("Duplicate collection refcode")
        seen.add(refcode)
        stem = row["CIF filename"].strip()
        member = cif_lookup.get(stem.lower())
        if member is None:
            raise ValueError(f"CIF listed in metadata is missing: {stem}")
        result.append(
            {
                "index": index,
                "refcode": refcode,
                "cif_member": member,
                "suspect_chemistry": truthy_marker(row["Unreliable chemistry"]),
                "sohncke_flag": row["Sohncke space group"].strip(),
            }
        )
    return result


def frozen_comparison(path: Path, selected_refcodes: set) -> Dict[str, Mapping[str, Any]]:
    result = {}
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            row = json.loads(line)
            refcode = str(row.get("refcode", "")).upper()
            if refcode in selected_refcodes:
                if refcode in result:
                    raise ValueError("Duplicate refcode in frozen ledger")
                result[refcode] = row
    if set(result) != selected_refcodes:
        raise ValueError("Frozen ledger lacks one or more selected refcodes")
    return result


def _worker_setup(zip_path: str) -> None:
    """Give each spawned process an independent ZIP handle and CCDC runtime."""
    global _WORKER_ARCHIVE
    _WORKER_ARCHIVE = zipfile.ZipFile(zip_path)


def _extract_one(record: Mapping[str, Any], archive: zipfile.ZipFile) -> Dict[str, Any]:
    try:
        from ccdc.crystal import Crystal
        cif_text = archive.read(record["cif_member"]).decode("utf-8")
    except Exception as exc:
        return {"ok": False, "error_stage": "archive_read", "error_type": type(exc).__name__}
    try:
        crystal = Crystal.from_string(cif_text, format="cif")
    except Exception as exc:
        return {"ok": False, "error_stage": "ccdc_parse", "error_type": type(exc).__name__}
    try:
        return {
            "ok": True,
            "entry": extract_entry(crystal, suspect_chemistry=record["suspect_chemistry"]),
        }
    except Exception as exc:
        return {
            "ok": False,
            "error_stage": "site_extract",
            "error_type": type(exc).__name__,
        }


def _worker_extract(record: Mapping[str, Any]) -> Dict[str, Any]:
    if _WORKER_ARCHIVE is None:
        raise RuntimeError("Worker archive was not initialized")
    return _extract_one(record, _WORKER_ARCHIVE)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mof-zip", type=Path, required=True)
    selection = parser.add_mutually_exclusive_group(required=True)
    selection.add_argument("--limit", type=int, help="Process first N metadata rows")
    selection.add_argument("--all", action="store_true", help="Explicit full-collection run")
    parser.add_argument("--start-index", type=int, default=0, help="Zero-based start for --limit (default: 0)")
    parser.add_argument("--workers", type=int, default=1, help="Isolated CCDC processes (default: 1)")
    parser.add_argument(
        "--compare-frozen-entries", type=Path,
        help="Optional private local QA ledger; never required or redistributed",
    )
    args = parser.parse_args()
    if not args.mof_zip.is_file():
        parser.error(f"MOF ZIP not found: {args.mof_zip}")
    if args.limit is not None and not 1 <= args.limit <= EXPECTED_COLLECTION_ENTRIES:
        parser.error("--limit must be between 1 and 15906")
    if args.start_index < 0 or (args.all and args.start_index != 0):
        parser.error("--start-index must be nonnegative and is only used with --limit")
    if args.limit is not None and args.start_index + args.limit > EXPECTED_COLLECTION_ENTRIES:
        parser.error("--start-index + --limit exceeds 15906")
    if not 1 <= args.workers <= 4:
        parser.error("--workers must be between 1 and 4")
    archive_hash = sha256(args.mof_zip)
    if archive_hash != ARCHIVE_SHA256:
        raise ValueError("MOF ZIP hash differs from the frozen v3 source")
    code_root = Path(__file__).resolve().parent / "periodic_v3"
    if sha256(code_root / "canonical_periodic_site.py") != CANONICALIZER_SHA256:
        raise ValueError("Frozen canonicalizer hash mismatch")
    if sha256(code_root / "canonical_csd_quotient.py") != QUOTIENT_SHA256:
        raise ValueError("Frozen quotient builder hash mismatch")
    try:
        import ccdc
    except ImportError as exc:
        raise SystemExit("Authorized CCDC Python API is required") from exc

    with zipfile.ZipFile(args.mof_zip) as archive:
        records = collection_index(archive)
        selected = records if args.all else records[args.start_index: args.start_index + args.limit]
        historical = (
            frozen_comparison(args.compare_frozen_entries, {row["refcode"] for row in selected})
            if args.compare_frozen_entries else None
        )
        started = time.perf_counter()
        totals: Counter[str] = Counter()
        errors: Counter[str] = Counter()
        gate_counts: Counter[str] = Counter()
        metric_to_stereo: Dict[str, set] = defaultdict(set)
        entry_pair_groups: Counter[str] = Counter()
        entry_state_groups: Counter[str] = Counter()
        family_groups: Dict[str, set] = defaultdict(set)
        paired_families: set = set()
        frozen_checks: Counter[str] = Counter()
        if args.workers == 1:
            extracted = (_extract_one(record, archive) for record in selected)
            executor = None
        else:
            executor = ProcessPoolExecutor(
                max_workers=args.workers,
                initializer=_worker_setup,
                initargs=(str(args.mof_zip),),
            )
            extracted = executor.map(_worker_extract, selected, chunksize=1)
        try:
            for ordinal, (record, outcome) in enumerate(zip(selected, extracted), 1):
              refcode = record["refcode"]
              try:
                  if not outcome["ok"]:
                      raise ValueError(outcome["error_stage"] + ":" + outcome["error_type"])
                  entry = outcome["entry"]
                  states = entry["states"]
                  mid_ids = [state["metric_id"] for state in states]
                  sid_ids = [state["stereo_id"] for state in states]
                  mid_hash = canonical_entry_multiset_id(mid_ids) if mid_ids else ""
                  sid_hash = (
                      canonical_entry_multiset_id(sid_ids, version="CR-PERIODIC-STEREO-ATLAS/1")
                      if sid_ids else ""
                  )
                  totals["processing_success"] += 1
                  totals["metal_sites"] += int(entry["metal_sites"])
                  totals["periodic_atom_orbits"] += int(entry["periodic_atom_orbits"])
                  totals["emitted_sites"] += len(states)
                  totals["sites_with_nonzero_translation"] += sum(
                      bool(state["has_nonzero_translation"]) for state in states
                  )
                  gate_counts.update(entry["gate_counts"])
                  entry_map: Dict[str, set] = defaultdict(set)
                  for state in states:
                      mid = state["metric_id"]
                      sid = state["stereo_id"]
                      metric_to_stereo[mid].add(sid)
                      entry_map[mid].add(sid)
                  if states:
                      group = (
                          "Sohncke" if record["sohncke_flag"].lower() == "yes"
                          else "not flagged Sohncke"
                      )
                      entry_state_groups[group] += 1
                      family_match = re.match(r"^([A-Z]+)", refcode)
                      family = family_match.group(1) if family_match else refcode
                      family_groups[family].add(group)
                      if any(len(values) > 1 for values in entry_map.values()):
                          entry_pair_groups[group] += 1
                          paired_families.add(family)
                  if historical is not None:
                      frozen = historical[refcode]
                      checks = {
                          "processing_success": bool(frozen.get("processing_success")),
                          "periodic_atom_orbits": int(frozen.get("n_periodic_atom_groups", -1))
                          == entry["periodic_atom_orbits"],
                          "metal_sites": int(frozen.get("n_in_domain_metal_sites", -1))
                          == entry["metal_sites"],
                          "emitted_sites": int(frozen.get("n_sites_representation_emitted", -1))
                          == len(states),
                          "nonzero_translation_sites": int(frozen.get("n_sites_nonzero_translation_edge", -1))
                          == entry["nonzero_translation_sites_all"],
                          "mid_multiset": frozen.get("canonical_metric_atlas_multiset_id", "") == mid_hash,
                          "sid_multiset": frozen.get("canonical_stereo_atlas_multiset_id", "") == sid_hash,
                      }
                      if not all(checks.values()):
                          raise AssertionError(
                              f"Frozen-entry comparison failed for {refcode}: "
                              + ",".join(key for key, ok in checks.items() if not ok)
                          )
                      frozen_checks["entries_exact"] += 1
              except Exception as exc:
                  if isinstance(exc, AssertionError):
                      raise
                  totals["processing_failed"] += 1
                  if not outcome["ok"]:
                      errors[f"{outcome['error_stage']}:{outcome['error_type']}"] += 1
                  else:
                      errors[f"{type(exc).__name__}: {str(exc)[:160]}"] += 1
                  if historical is not None:
                      frozen = historical[refcode]
                      if bool(frozen.get("processing_success")):
                          raise AssertionError(f"New failure for {refcode}") from exc
                      frozen_checks["matched_failure_entries"] += 1
              progress_every = 500 if args.all else 25
              if ordinal % progress_every == 0 or ordinal == len(selected):
                  print(f"[periodic-v3] {ordinal}/{len(selected)} entries", file=sys.stderr, flush=True)
        finally:
            if executor is not None:
                executor.shutdown(wait=True, cancel_futures=True)

    multiplicity = Counter(len(values) for values in metric_to_stereo.values())
    family_totals: Counter[str] = Counter()
    family_pairs: Counter[str] = Counter()
    for family, labels in family_groups.items():
        label = next(iter(labels)) if len(labels) == 1 else "mixed"
        family_totals[label] += 1
        family_pairs[label] += int(family in paired_families)
    elapsed = time.perf_counter() - started
    result = {
        "status": "partial_smoke" if not args.all else "full_zip_run_requires_claim_review",
        "source_archive_sha256": archive_hash,
        "environment": {"ccdc_api": ccdc.__version__, "numpy": __import__("numpy").__version__},
        "selected_entries": len(selected),
        "start_index": args.start_index,
        "workers": args.workers,
        "totals": dict(sorted(totals.items())),
        "gate_counts_nonexclusive": dict(sorted(gate_counts.items())),
        "error_counts": dict(sorted(errors.items())),
        "unique_metric_ids": len(metric_to_stereo),
        "unique_stereo_ids": len({sid for values in metric_to_stereo.values() for sid in values}),
        "metric_to_stereo_multiplicity": {str(key): value for key, value in sorted(multiplicity.items())},
        "state_bearing_entries_by_flag": dict(sorted(entry_state_groups.items())),
        "paired_entries_by_flag": dict(sorted(entry_pair_groups.items())),
        "state_bearing_families_by_flag": dict(sorted(family_totals.items())),
        "paired_families_by_flag": dict(sorted(family_pairs.items())),
        "private_frozen_QA": dict(sorted(frozen_checks.items())) if historical is not None else "not_requested",
        "elapsed_seconds": round(elapsed, 3),
        "seconds_per_selected_entry": round(elapsed / len(selected), 3),
        "claim_boundary": "Only a full run with exact locked counts is a CIF-level reproduction of the atlas",
    }
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
