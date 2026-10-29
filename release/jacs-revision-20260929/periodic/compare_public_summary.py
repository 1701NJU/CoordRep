#!/usr/bin/env python3
"""Check aggregate-only batch output against the public Figure 6 summary.

This does not read private ledgers, CIF/site rows, or refcode-level records.
The secondary CShM/shape-conflict metric is intentionally outside the
code-only primary MID/SID runner's scope.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Dict


def _percent(part: int, whole: int) -> float:
    return 100.0 * part / whole if whole else 0.0


def _compare(expected: Any, actual: Any, field: str) -> None:
    if isinstance(actual, dict):
        if not isinstance(expected, dict):
            raise AssertionError(f"{field}: expected a mapping")
        if set(expected) != set(actual):
            raise AssertionError(f"{field}: group labels differ")
        for key, value in actual.items():
            _compare(expected[key], value, f"{field}.{key}")
    elif isinstance(actual, float):
        if not math.isclose(float(expected), actual, rel_tol=0.0, abs_tol=1e-9):
            raise AssertionError(f"{field}: public={expected}, recomputed={actual}")
    elif expected != actual:
        raise AssertionError(f"{field}: public={expected}, recomputed={actual}")


def extract_headline(batch: Dict[str, Any]) -> Dict[str, Any]:
    if batch.get("status") != "full_zip_run_requires_claim_review":
        raise ValueError("Batch output is not from an explicit --all run")
    if batch.get("selected_entries") != 15906 or batch.get("start_index") != 0:
        raise ValueError("Batch output does not cover the full pinned collection")
    multiplicity = batch["metric_to_stereo_multiplicity"]
    state_entries = batch["state_bearing_entries_by_flag"]
    paired_entries = batch["paired_entries_by_flag"]
    state_families = batch["state_bearing_families_by_flag"]
    paired_families = batch["paired_families_by_flag"]
    n_state = sum(state_entries.values())
    n_paired = sum(paired_entries.values())
    split = sum(count for cardinality, count in multiplicity.items() if int(cardinality) > 1)
    return {
        "emitted_site_instances": batch["totals"]["emitted_sites"],
        "state_bearing_entries": n_state,
        "unique_metric_ids": batch["unique_metric_ids"],
        "unique_stereo_ids": batch["unique_stereo_ids"],
        "metric_ids_with_multiple_stereo_ids": split,
        "metric_split_fraction_percent": _percent(split, batch["unique_metric_ids"]),
        "metric_stereo_multiplicity": multiplicity,
        "entries_with_within_entry_pair": n_paired,
        "paired_entry_fraction_percent": _percent(n_paired, n_state),
        "spacegroup_blind_join": {
            label: {
                "state_bearing_entries": number,
                "paired_entries": paired_entries.get(label, 0),
                "paired_fraction_percent": _percent(paired_entries.get(label, 0), number),
            }
            for label, number in state_entries.items()
        },
        "refcode_family_control": {
            label: {
                "state_bearing_families": number,
                "paired_families": paired_families.get(label, 0),
                "paired_fraction_percent": _percent(paired_families.get(label, 0), number),
            }
            for label, number in state_families.items()
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batch-json", type=Path, required=True)
    parser.add_argument("--public-summary", type=Path, required=True)
    args = parser.parse_args()
    with args.batch_json.open("r", encoding="utf-8-sig") as handle:
        batch = json.load(handle)
    with args.public_summary.open("r", encoding="utf-8-sig") as handle:
        public = json.load(handle)
    observed = extract_headline(batch)
    for field, value in observed.items():
        _compare(public[field], value, field)
    print(json.dumps({
        "status": "PASS_public_Figure6_primary_MID_SID_aggregates_match_independent_ZIP_run",
        "checked_top_level_fields": sorted(observed),
        "excluded_secondary_field": "shape_conflicting_metric_ids (requires corrected CShM)",
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
