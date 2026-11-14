#!/usr/bin/env python3
"""Summarize the preregistered four-stratum CSD audit pilot."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path


def _wilson(successes: int, total: int, z: float = 1.959963984540054) -> dict:
    if total == 0:
        return {"lower": None, "upper": None}
    p = successes / total
    denominator = 1 + z * z / total
    centre = (p + z * z / (2 * total)) / denominator
    half = z * math.sqrt(p * (1 - p) / total + z * z / (4 * total * total)) / denominator
    return {"lower": centre - half, "upper": centre + half}


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--outcomes", type=Path, required=True)
    parser.add_argument("--fidelity-summary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    manifest = [json.loads(line) for line in args.manifest.open(encoding="utf-8") if line.strip()]
    outcomes = [json.loads(line) for line in args.outcomes.open(encoding="utf-8") if line.strip()]
    fidelity = json.loads(args.fidelity_summary.read_text(encoding="utf-8"))
    by_refcode = {row["refcode"]: row for row in outcomes}
    if len(by_refcode) != len(outcomes):
        raise SystemExit("outcome refcodes are not unique")

    strata = defaultdict(Counter)
    missing = []
    for row in manifest:
        refcode = row["refcode"]
        outcome = by_refcode.get(refcode)
        if outcome is None:
            missing.append(refcode)
            continue
        counter = strata[row["sample_stratum"]]
        counter["total"] += 1
        if outcome.get("status", "").startswith("EMITTED_"):
            counter["record_emitted"] += 1
        if outcome.get("structural_record"):
            counter["structural_record"] += 1
        counter[outcome.get("status", "UNKNOWN")] += 1
        counter[f"output_scope:{outcome.get('scope', 'unknown')}"] += 1

    stratum_rows = {}
    for name, counter in sorted(strata.items()):
        interval = _wilson(counter["structural_record"], counter["total"])
        stratum_rows[name] = {
            **dict(counter),
            "structural_rate": counter["structural_record"] / counter["total"],
            "wilson_95": interval,
        }
    total = len(manifest)
    structural = sum(counter["structural_record"] for counter in strata.values())
    payload = {
        "manifest_rows": total,
        "outcome_rows": len(outcomes),
        "missing_outcomes": missing,
        "structural_records": structural,
        "structural_rate": structural / total if total else None,
        "overall_wilson_95": _wilson(structural, total),
        "strata": stratum_rows,
        "source_fidelity": fidelity,
        "manifest_sha256": hashlib.sha256(args.manifest.read_bytes()).hexdigest(),
        "outcomes_sha256": hashlib.sha256(args.outcomes.read_bytes()).hexdigest(),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
