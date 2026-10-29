#!/usr/bin/env python3
"""Audit geometry-assigned local-shell candidates for native orphan sites.

This is a post hoc, provenance-preserving sensitivity analysis.  Native CSD
bonds are never replaced.  A deep copy of the crystal is re-perceived only
for entries containing metal sites with neither donor incidences nor direct
metal--metal relations in the frozen native graph.  Recovered objects are
local-shell candidates, not schema-valid augmented CoordRep records, and must
remain separate from native structural-record coverage.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import sys
from collections import Counter
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

import numpy as np
import ccdc
from ccdc.io import EntryReader


BASE = Path(__file__).resolve().parent.parent
if str(BASE) not in sys.path:
    sys.path.insert(0, str(BASE))

from coordrep.audit.records import ENTRY_RESOLUTIONS
from coordrep.io.tmqm_reader import TRANSITION_METALS


FAILED_STATUSES = {"FAILED_INPUT", "FAILED_PIPELINE"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--database", required=True)
    parser.add_argument("--outcomes", type=Path, required=True)
    parser.add_argument("--output-internal", type=Path, required=True)
    parser.add_argument("--output-summary", type=Path, required=True)
    parser.add_argument("--coordinate-tolerance", type=float, default=1.0e-4)
    parser.add_argument("--expected-database-entries", type=int)
    parser.add_argument("--expected-outcomes", type=int)
    parser.add_argument("--database-sha256")
    parser.add_argument("--outcomes-sha256")
    parser.add_argument("--merged-summary", type=Path)
    return parser.parse_args()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _atomic_json(path: Path, payload: Dict[str, Any]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    os.replace(temporary, path)


def _coordinates(atom: Any) -> Optional[np.ndarray]:
    try:
        value = atom.coordinates
        result = np.asarray([value.x, value.y, value.z], dtype=float)
    except Exception:
        return None
    return result if result.shape == (3,) and np.all(np.isfinite(result)) else None


def _bond_type(bond: Any) -> str:
    try:
        return str(bond.bond_type)
    except Exception:
        return "unknown"


def _native_orphan_metals(entry: Any) -> List[Any]:
    molecule = entry.molecule
    atoms = list(molecule.atoms)
    metals = [atom for atom in atoms if atom.atomic_symbol in TRANSITION_METALS]
    incident = Counter()
    for bond in molecule.bonds:
        endpoints = list(bond.atoms)
        if len(endpoints) != 2:
            continue
        for atom in endpoints:
            if atom.atomic_symbol in TRANSITION_METALS:
                incident[int(atom.index)] += 1
    return [atom for atom in metals if incident[int(atom.index)] == 0]


def _same_environment(candidates: Sequence[Dict[str, Any]]) -> bool:
    if not candidates:
        return False
    signatures = {
        (
            tuple(candidate["neighbor_elements"]),
            tuple(candidate["bond_types"]),
            tuple(candidate["distances"]),
            int(candidate.get("endpoint_mapping_failures", 0)),
        )
        for candidate in candidates
    }
    return len(signatures) == 1


def _other_bond_endpoint(bond: Any, atom: Any) -> Optional[Any]:
    """Return the other endpoint using molecule-local indices.

    CCDC Atom properties create fresh wrapper objects, so Python object
    identity is not a valid test of whether an endpoint denotes ``atom``.
    """
    endpoints = list(bond.atoms)
    if len(endpoints) != 2:
        return None
    try:
        atom_index = int(atom.index)
        matches = [int(endpoint.index) == atom_index for endpoint in endpoints]
    except Exception:
        return None
    if matches == [True, False]:
        return endpoints[1]
    if matches == [False, True]:
        return endpoints[0]
    return None


def _perceived_candidates(
    perceived_atoms: Sequence[Any], source_metal: Any, tolerance: float
) -> List[Dict[str, Any]]:
    source_coord = _coordinates(source_metal)
    if source_coord is None:
        return []
    source_label = str(getattr(source_metal, "label", ""))
    matches = []
    for atom in perceived_atoms:
        if atom.atomic_symbol != source_metal.atomic_symbol:
            continue
        if str(getattr(atom, "label", "")) != source_label:
            continue
        coord = _coordinates(atom)
        if coord is None or float(np.linalg.norm(coord - source_coord)) > tolerance:
            continue
        neighbours = []
        endpoint_mapping_failures = 0
        for bond in atom.bonds:
            other = _other_bond_endpoint(bond, atom)
            if other is None:
                endpoint_mapping_failures += 1
                continue
            other_coord = _coordinates(other)
            distance = None
            if other_coord is not None:
                distance = round(float(np.linalg.norm(other_coord - coord)), 4)
            neighbours.append((
                str(other.atomic_symbol),
                _bond_type(bond),
                distance,
            ))
        neighbours.sort(key=repr)
        matches.append({
            "neighbor_elements": [item[0] for item in neighbours],
            "bond_types": [item[1] for item in neighbours],
            "distances": [item[2] for item in neighbours],
            "endpoint_mapping_failures": endpoint_mapping_failures,
        })
    return matches


def main() -> None:
    args = parse_args()
    if args.coordinate_tolerance <= 0:
        raise SystemExit("coordinate-tolerance must be positive")
    reader = EntryReader(args.database)
    if (
        args.expected_database_entries is not None
        and len(reader) != args.expected_database_entries
    ):
        raise SystemExit(
            f"database contains {len(reader)} entries, expected "
            f"{args.expected_database_entries}"
        )
    database_path = Path(args.database)
    observed_database_sha256 = None
    if args.database_sha256:
        if not database_path.is_file():
            raise SystemExit("database-sha256 requires an explicit database file")
        observed_database_sha256 = _sha256(database_path)
        if observed_database_sha256.lower() != args.database_sha256.lower():
            raise SystemExit(
                f"database SHA-256 mismatch: {observed_database_sha256} != "
                f"{args.database_sha256.lower()}"
            )
    observed_outcomes_sha256 = _sha256(args.outcomes)
    if (
        args.outcomes_sha256
        and observed_outcomes_sha256.lower() != args.outcomes_sha256.lower()
    ):
        raise SystemExit(
            f"outcomes SHA-256 mismatch: {observed_outcomes_sha256} != "
            f"{args.outcomes_sha256.lower()}"
        )
    merged_summary_sha256 = None
    merged_audit_verified = False
    if args.merged_summary is not None:
        merged_summary = json.loads(args.merged_summary.read_text(encoding="utf-8"))
        merged_summary_sha256 = _sha256(args.merged_summary)
        required_values_present = all(
            value is not None
            for value in (
                args.expected_database_entries,
                args.expected_outcomes,
                observed_database_sha256,
                args.outcomes_sha256,
            )
        )
        if not required_values_present:
            raise SystemExit(
                "merged-summary verification also requires both expected counts "
                "and both input SHA-256 arguments"
            )
        if int(merged_summary["database_entries_evaluated"]) != args.expected_database_entries:
            raise SystemExit("merged-summary database denominator mismatch")
        if int(merged_summary["transition_metal_entries"]) != args.expected_outcomes:
            raise SystemExit("merged-summary transition-metal denominator mismatch")
        if str(merged_summary.get("database_sha256", "")).lower() != (
            observed_database_sha256 or ""
        ).lower():
            raise SystemExit("merged-summary database SHA-256 mismatch")
        if str(
            merged_summary.get("artifacts", {}).get("outcomes_internal_sha256", "")
        ).lower() != observed_outcomes_sha256.lower():
            raise SystemExit("merged-summary outcomes SHA-256 mismatch")
        merged_audit_verified = True

    targets: List[Dict[str, Any]] = []
    base = Counter()
    input_refcodes = set()
    nonnegative_source_indices = set()
    saw_negative_source_index = False
    with args.outcomes.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            row = json.loads(line)
            refcode = str(row["refcode"])
            source_index = int(row["source_index"])
            if refcode in input_refcodes:
                raise SystemExit(f"duplicate outcome refcode: {refcode}")
            input_refcodes.add(refcode)
            if source_index >= 0:
                if source_index in nonnegative_source_indices:
                    raise SystemExit(f"duplicate nonnegative source index: {source_index}")
                nonnegative_source_indices.add(source_index)
            else:
                saw_negative_source_index = True

            status = str(row.get("status", ""))
            if status not in ENTRY_RESOLUTIONS | FAILED_STATUSES:
                raise SystemExit(f"unknown terminal status for {refcode}: {status}")
            metal_sites = int(row.get("metal_sites", 0))
            structural_sites = int(row.get("structural_metal_sites", 0))
            audit_sites = int(row.get("audit_only_metal_sites", 0))
            failed_sites = int(row.get("failed_metal_sites", 0))
            if min(metal_sites, structural_sites, audit_sites, failed_sites) < 0:
                raise SystemExit(f"negative site count for {refcode}")
            if structural_sites + audit_sites + failed_sites != metal_sites:
                raise SystemExit(f"metal-site fields do not close for {refcode}")
            if status in ENTRY_RESOLUTIONS:
                if row.get("schema_valid") is not True:
                    raise SystemExit(f"emitted outcome is not schema-valid for {refcode}")
                if failed_sites:
                    raise SystemExit(f"emitted outcome has failed sites for {refcode}")
                if bool(row.get("structural_record")) != (structural_sites > 0):
                    raise SystemExit(f"structural-record flag mismatch for {refcode}")
                if bool(row.get("all_metal_sites_structural")) != (audit_sites == 0):
                    raise SystemExit(f"all-sites flag mismatch for {refcode}")
                if (status == "EMITTED_AUDIT_ONLY") != (structural_sites == 0):
                    raise SystemExit(f"audit-only terminal status mismatch for {refcode}")
            elif structural_sites or audit_sites or failed_sites != metal_sites:
                raise SystemExit(f"failed outcome site fields do not close for {refcode}")

            base["entries"] += 1
            base["metal_sites"] += metal_sites
            base["structural_metal_sites"] += structural_sites
            base["audit_only_metal_sites"] += audit_sites
            base["failed_metal_sites"] += failed_sites
            base["entries_any_structural"] += int(bool(row.get("structural_record")))
            base["entries_all_sites_structural"] += int(
                bool(row.get("all_metal_sites_structural"))
            )
            if status in ENTRY_RESOLUTIONS and audit_sites:
                targets.append({
                    "refcode": refcode,
                    "source_index": source_index,
                    "expected_orphans": audit_sites,
                    "input_status": status,
                    "scope_flags": tuple(row.get("scope_flags", [])),
                    "issue_codes": tuple(row.get("issue_codes", [])),
                    "had_structural_site": structural_sites > 0,
                })

    if args.expected_outcomes is not None and base["entries"] != args.expected_outcomes:
        raise SystemExit(
            f"outcomes contain {base['entries']} rows, expected {args.expected_outcomes}"
        )
    if base["metal_sites"] != (
        base["structural_metal_sites"]
        + base["audit_only_metal_sites"]
        + base["failed_metal_sites"]
    ):
        raise SystemExit("aggregate native metal-site fields do not close")

    args.output_internal.parent.mkdir(parents=True, exist_ok=True)
    aggregate = Counter()
    recovered_by_refcode: Dict[str, int] = {}
    with args.output_internal.open("w", encoding="utf-8", newline="\n") as output:
        for position, target in enumerate(targets, 1):
            refcode = target["refcode"]
            source_index = target["source_index"]
            expected_orphans = target["expected_orphans"]
            if source_index >= 0:
                entry = reader[source_index]
                if str(getattr(entry, "identifier", "")) != refcode:
                    raise SystemExit(
                        f"source-index/refcode mismatch: {source_index} != {refcode}"
                    )
            else:
                entry = reader.entry(refcode)
            if entry is None:
                raise SystemExit(f"entry not found in supplied database: {refcode}")
            orphan_metals = _native_orphan_metals(entry)
            if len(orphan_metals) != expected_orphans:
                raise SystemExit(
                    f"native orphan-site mismatch for {refcode}: reconstructed "
                    f"{len(orphan_metals)}, outcome reports {expected_orphans}"
                )
            aggregate["target_entries"] += 1
            aggregate["expected_orphan_sites"] += expected_orphans
            aggregate["reconstructed_native_orphans"] += len(orphan_metals)
            try:
                perceived = entry.crystal.copy()
                perceived.assign_bonds()
                perceived_atoms = list(perceived.molecule.atoms)
                perception_error = ""
            except Exception as exc:
                perceived_atoms = []
                perception_error = f"{type(exc).__name__}: {str(exc)[:200]}"
                aggregate["perception_exception_entries"] += 1

            recovered = 0
            site_rows = []
            for metal in orphan_metals:
                candidates = _perceived_candidates(
                    perceived_atoms, metal, args.coordinate_tolerance
                )
                accepted = bool(candidates) and _same_environment(candidates)
                representative = candidates[0] if accepted else {
                    "neighbor_elements": [],
                    "bond_types": [],
                    "distances": [],
                    "endpoint_mapping_failures": 0,
                }
                valid_distances = (
                    bool(representative["distances"])
                    and all(
                        distance is not None and distance > args.coordinate_tolerance
                        for distance in representative["distances"]
                    )
                )
                accepted_site = (
                    accepted
                    and bool(representative["neighbor_elements"])
                    and valid_distances
                    and not representative["endpoint_mapping_failures"]
                )
                if accepted_site:
                    recovered += 1
                    aggregate["recovered_orphan_sites"] += 1
                    retains_ambiguity = (
                        target["input_status"] == "EMITTED_AMBIGUOUS"
                        or "disorder" in target["scope_flags"]
                        or len(candidates) > 1
                    )
                    status = (
                        "GEOMETRY_ASSIGNED_CANDIDATE_AMBIGUOUS"
                        if retains_ambiguity
                        else "GEOMETRY_ASSIGNED_CANDIDATE_PARTIAL"
                    )
                    issue = "GEOMETRY_ASSIGNED_FALLBACK"
                else:
                    aggregate["unrecovered_orphan_sites"] += 1
                    status = "GEOMETRY_ASSIGNMENT_UNRECOVERED"
                    issue = "GEOMETRY_ASSIGNMENT_EMPTY_OR_NONUNIQUE"
                source_token = hashlib.sha256(
                    (
                        refcode
                        + "|"
                        + str(getattr(metal, "label", ""))
                        + "|"
                        + str(getattr(metal, "index", ""))
                    ).encode("utf-8")
                ).hexdigest()
                site_rows.append({
                    "source_metal_token_sha256": source_token,
                    "metal_element": str(metal.atomic_symbol),
                    "candidate_matches": len(candidates),
                    "neighbor_elements": representative["neighbor_elements"],
                    "bond_types": representative["bond_types"],
                    "distances": representative["distances"],
                    "endpoint_mapping_failures": representative[
                        "endpoint_mapping_failures"
                    ],
                    "candidate_evidence_level": (
                        "local_shell_only" if accepted_site else "none"
                    ),
                    "resolution": (
                        "ambiguous"
                        if accepted_site and status.endswith("AMBIGUOUS")
                        else "partial" if accepted_site else "unresolved"
                    ),
                    "status": status,
                    "evidence": "CCDC_geometry_assigned" if accepted_site else "none",
                    "issue_code": issue,
                })
            recovered_by_refcode[refcode] = recovered
            output.write(json.dumps({
                "source_index": source_index,
                "refcode": refcode,
                "expected_audit_only_sites": expected_orphans,
                "native_orphan_sites_reconstructed": len(orphan_metals),
                "recovered_sites": recovered,
                "perception_error": perception_error,
                "input_status": target["input_status"],
                "input_scope_flags": target["scope_flags"],
                "input_issue_codes": target["issue_codes"],
                "sites": site_rows,
            }, sort_keys=True, separators=(",", ":")) + "\n")
            if position % 100 == 0:
                print(json.dumps({
                    "entries": position,
                    "recovered_sites": aggregate["recovered_orphan_sites"],
                }), flush=True)

    candidate_sites = base["structural_metal_sites"] + aggregate["recovered_orphan_sites"]
    candidate_any = base["entries_any_structural"]
    candidate_all = base["entries_all_sites_structural"]
    for target in targets:
        refcode = target["refcode"]
        expected_orphans = target["expected_orphans"]
        recovered = recovered_by_refcode.get(refcode, 0)
        if recovered and not target["had_structural_site"]:
            candidate_any += 1
        if recovered == expected_orphans:
            candidate_all += 1

    if not 0 <= candidate_any <= base["entries"]:
        raise SystemExit("candidate entry-any numerator is outside its denominator")
    if not 0 <= candidate_all <= base["entries"]:
        raise SystemExit("candidate entry-all numerator is outside its denominator")
    if not 0 <= candidate_sites <= base["metal_sites"]:
        raise SystemExit("candidate site numerator is outside its denominator")
    summary = {
        "method": "CCDC_geometry_assigned_local_shells_on_deep_crystal_copy_for_native_orphan_sites_only",
        "interpretation": (
            "Sensitivity analysis only: recovered shells are not schema-valid "
            "augmented CoordRep records and are excluded from native coverage."
        ),
        "coordinate_tolerance_angstrom": args.coordinate_tolerance,
        "provenance": {
            "database_argument": args.database,
            "database_entries": len(reader),
            "expected_database_entries": args.expected_database_entries,
            "database_sha256": observed_database_sha256,
            "outcomes_sha256": observed_outcomes_sha256,
            "expected_outcomes": args.expected_outcomes,
            "merged_summary_sha256": merged_summary_sha256,
            "fallback_script_sha256": _sha256(Path(__file__).resolve()),
            "frozen_merged_audit_verified": bool(
                merged_audit_verified
                and not saw_negative_source_index
            ),
            "runtime": {
                "python": sys.version.replace("\n", " "),
                "implementation": platform.python_implementation(),
                "platform": platform.platform(),
                "ccdc": str(getattr(ccdc, "__version__", "unknown")),
                "numpy": str(np.__version__),
            },
        },
        "native": {
            "entries": base["entries"],
            "metal_sites": base["metal_sites"],
            "entries_any_structural": base["entries_any_structural"],
            "entries_all_sites_structural": base["entries_all_sites_structural"],
            "structural_metal_sites": base["structural_metal_sites"],
            "audit_only_metal_sites": base["audit_only_metal_sites"],
            "failed_metal_sites": base["failed_metal_sites"],
        },
        "fallback": dict(aggregate),
        "geometry_assigned_candidate_sensitivity": {
            "entries_any_candidate_or_native_shell": candidate_any,
            "entry_any_candidate_or_native_fraction": (
                candidate_any / base["entries"] if base["entries"] else None
            ),
            "entries_all_sites_candidate_or_native_shell": candidate_all,
            "entry_all_sites_candidate_or_native_fraction": (
                candidate_all / base["entries"] if base["entries"] else None
            ),
            "candidate_or_native_metal_sites": candidate_sites,
            "candidate_or_native_metal_site_fraction": (
                candidate_sites / base["metal_sites"] if base["metal_sites"] else None
            ),
        },
        "internal_output_sha256": hashlib.sha256(
            args.output_internal.read_bytes()
        ).hexdigest(),
    }
    args.output_summary.parent.mkdir(parents=True, exist_ok=True)
    _atomic_json(args.output_summary, summary)
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
