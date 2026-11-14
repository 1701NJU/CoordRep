#!/usr/bin/env python3
"""Build the preregistered four-stratum CSD audit pilot manifest."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Any, Dict, List

from ccdc.io import EntryReader

BASE = Path(__file__).resolve().parent.parent
if str(BASE) not in sys.path:
    sys.path.insert(0, str(BASE))

from coordrep.audit.metal_policy import (
    METAL_POLICY_ID,
    METAL_POLICY_SHA256,
    is_target_metal_symbol,
    metal_blocks,
)


DEFAULT_NAMESPACE = "coordrep-all-metal-csd-audit-pilot-20260825-v1"
STRATA = (
    "mono_no_native_pi",
    "mono_native_pi_flag",
    "multi_no_native_pi",
    "multi_native_pi_flag",
)


def _rank(namespace: str, index: int) -> int:
    digest = hashlib.sha256(f"{namespace}|{index}".encode("ascii")).digest()
    return int.from_bytes(digest[:8], "big", signed=False)


def _is_pi(bond: Any) -> bool:
    try:
        text = str(bond.bond_type).lower()
    except Exception:
        text = ""
    return "pi" in text or "deloc" in text


def _classify(entry: Any) -> Dict[str, Any]:
    if not bool(getattr(entry, "has_3d_structure", False)):
        return {"stratum": "non_target"}
    molecule = entry.molecule
    if molecule is None:
        return {"stratum": "non_target"}
    metals = [
        atom
        for atom in molecule.atoms
        if is_target_metal_symbol(getattr(atom, "atomic_symbol", "?"))
    ]
    if not metals:
        return {"stratum": "non_target"}
    has_pi = False
    for metal in metals:
        for bond in metal.bonds:
            if _is_pi(bond):
                has_pi = True
                break
        if has_pi:
            break
    prefix = "mono" if len(metals) == 1 else "multi"
    suffix = "native_pi_flag" if has_pi else "no_native_pi"
    return {
        "stratum": f"{prefix}_{suffix}",
        "n_metals": len(metals),
        "metal_elements": sorted({str(atom.atomic_symbol) for atom in metals}),
        "metal_blocks": metal_blocks(atom.atomic_symbol for atom in metals),
        "has_native_pi_or_delocalized_bond": has_pi,
        "has_disorder_entry_flag": bool(getattr(entry, "has_disorder", False)),
        "is_polymeric": bool(getattr(entry, "is_polymeric", False)),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--database", default="CSD")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--namespace", default=DEFAULT_NAMESPACE)
    parser.add_argument("--per-stratum", type=int, default=400)
    parser.add_argument("--probe-sizes", default="20000,40000,65536,131072")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    reader = EntryReader(args.database)
    total = len(reader)
    probe_sizes = [int(item) for item in args.probe_sizes.split(",")]
    if sorted(probe_sizes) != probe_sizes or probe_sizes[-1] > total:
        raise SystemExit("probe sizes must be increasing and no larger than the database")

    ranked_indices = sorted(range(total), key=lambda index: _rank(args.namespace, index))
    selected: Dict[str, List[Dict[str, Any]]] = {stratum: [] for stratum in STRATA}
    counters: Counter = Counter()
    previous_stop = 0
    terminal_probe_size = None
    for probe_size in probe_sizes:
        for position in range(previous_stop, probe_size):
            index = ranked_indices[position]
            counters["entries_read"] += 1
            try:
                entry = reader[index]
                facts = _classify(entry)
            except Exception:
                counters["entry_or_classification_error"] += 1
                continue
            stratum = facts["stratum"]
            counters[f"observed:{stratum}"] += 1
            if stratum not in selected or len(selected[stratum]) >= args.per_stratum:
                continue
            selected[stratum].append({
                "source_index": index,
                "refcode": str(entry.identifier),
                "sample_role": "probability_stratum",
                "sample_stratum": stratum,
                "sample_rank_u64": _rank(args.namespace, index),
                **{key: value for key, value in facts.items() if key != "stratum"},
            })
        previous_stop = probe_size
        if all(len(selected[stratum]) >= args.per_stratum for stratum in STRATA):
            terminal_probe_size = probe_size
            break
    if terminal_probe_size is None:
        missing = {key: len(rows) for key, rows in selected.items() if len(rows) < args.per_stratum}
        raise SystemExit(f"insufficient preregistered probe pool: {missing}")

    rows = sorted(
        (row for stratum in STRATA for row in selected[stratum]),
        key=lambda row: (row["sample_stratum"], row["sample_rank_u64"], row["refcode"]),
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="\n") as handle:
        for row in rows:
            handle.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")
    summary = {
        "namespace": args.namespace,
        "metal_policy_id": METAL_POLICY_ID,
        "metal_policy_sha256": METAL_POLICY_SHA256,
        "database_argument": args.database,
        "database_file": str(getattr(reader, "file_name", "")),
        "database_entries": total,
        "per_stratum": args.per_stratum,
        "terminal_probe_size": terminal_probe_size,
        "manifest_rows": len(rows),
        "stratum_counts": {key: len(value) for key, value in selected.items()},
        "probe_counters": dict(counters),
        "manifest_sha256": hashlib.sha256(args.output.read_bytes()).hexdigest(),
    }
    summary_path = args.output.with_suffix(args.output.suffix + ".summary.json")
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
