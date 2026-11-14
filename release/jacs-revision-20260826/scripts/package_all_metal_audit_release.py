#!/usr/bin/env python3
"""Build the public aggregate evidence package for the all-metal CSD audit."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any, Dict, Iterable, List


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--census", type=Path, required=True)
    parser.add_argument("--merged-summary", type=Path, required=True)
    parser.add_argument("--source-fidelity", type=Path, required=True)
    parser.add_argument("--claim-ready", type=Path, required=True)
    parser.add_argument("--policy-probe", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read(path: Path) -> Dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_json(path: Path, payload: Dict[str, Any]) -> None:
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def _percent(numerator: int, denominator: int) -> str:
    return f"{100.0 * numerator / denominator:.4f}%" if denominator else ""


def _write_csv(path: Path, fieldnames: Iterable[str], rows: List[Dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    inputs = {
        "metal_domain_census": args.census,
        "record_audit_summary": args.merged_summary,
        "source_fidelity_summary": args.source_fidelity,
        "claim_ready_summary": args.claim_ready,
        "metal_policy_probe": args.policy_probe,
    }
    census = _read(args.census)
    merged = _read(args.merged_summary)
    fidelity = _read(args.source_fidelity)
    claim = _read(args.claim_ready)
    probe = _read(args.policy_probe)

    target = int(census["metal_containing_3d_entries"])
    emitted = int(merged["records_emitted"])
    if not (
        target == int(merged["target_entries"]) == int(claim["target_entries"])
        and emitted == target
        and int(merged["failed_outcomes"]) == 0
        and int(fidelity["records"]) == target
        and int(fidelity["passed"]) == target
        and int(fidelity["failed"]) == 0
    ):
        raise SystemExit("audit, census, and source-fidelity denominators do not close")
    for key in ("database_sha256", "metal_policy_id", "metal_policy_sha256"):
        claim_value = (
            claim.get("provenance", {}).get("database_sha256")
            if key == "database_sha256"
            else claim.get(key)
        )
        values = {census.get(key), merged.get(key), claim_value}
        if key == "database_sha256":
            values.add(fidelity.get(key))
        if len(values) != 1:
            raise SystemExit(f"cross-artifact provenance mismatch for {key}")
    if not all(claim.get("accounting_checks", {}).values()):
        raise SystemExit("one or more claim-ready accounting checks failed")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    protected = [
        args.output_dir / "METAL_DOMAIN_CENSUS_PUBLIC.json",
        args.output_dir / "STRUCTURAL_RECORD_AUDIT_PUBLIC.json",
        args.output_dir / "SOURCE_FIDELITY_PUBLIC.json",
        args.output_dir / "ALL_METAL_CLAIM_READY_SUMMARY.json",
        args.output_dir / "RUN_PROVENANCE.json",
        args.output_dir / "TABLE_S9_SOURCE.csv",
        args.output_dir / "FIGURE5_SOURCE.csv",
        args.output_dir / "README.md",
        args.output_dir / "SHA256SUMS.txt",
    ]
    if any(path.exists() for path in protected):
        raise SystemExit("public release output already exists")

    _write_json(protected[0], census)
    _write_json(protected[1], merged)
    fidelity_public = {
        key: fidelity.get(key)
        for key in (
            "records",
            "passed",
            "failed",
            "pass_rate",
            "shards",
            "validator_sha256",
            "database_sha256",
            "metal_policy_id",
            "metal_policy_sha256",
            "runtime",
        )
    }
    _write_json(protected[2], fidelity_public)
    _write_json(protected[3], claim)

    total_csd = int(census["database_entries"])
    metal_all = int(census["metal_containing_entries"])
    metal_no3d = int(census["metal_containing_without_3d_entries"])
    observed_metal_elements = len(census["metal_element_entry_counts"])
    observed_3d_metal_elements = len(census["metal_3d_element_entry_counts"])
    structural_entries = int(
        claim["coverage"]["entry_any_structural"]["numerator"]
    )
    all_sites_entries = int(
        claim["coverage"]["entry_all_metal_sites_structural"]["numerator"]
    )
    total_sites = int(merged["metal_sites_accounted"])
    structural_sites = int(merged["structural_metal_sites"])
    audit_sites = int(merged["audit_only_metal_sites"])
    status = claim["status_counts"]
    chemistry = claim["chemical_scope"]

    table_rows = [
        {"Analysis": "April 2025 CSD", "Quantity": "All entries", "Count": total_csd, "Fraction": "100.0000%", "Interpretation": "Frozen database denominator"},
        {"Analysis": "Frozen metal-domain census", "Quantity": "Metal-containing entries", "Count": metal_all, "Fraction": _percent(metal_all, total_csd), "Interpretation": "At least one symbol in the frozen 96-symbol CCDC-derived policy"},
        {"Analysis": "Frozen metal-domain census", "Quantity": "Metal-containing entries without 3D", "Count": metal_no3d, "Fraction": _percent(metal_no3d, metal_all), "Interpretation": "Census only; no 3D structural record assigned"},
        {"Analysis": "Structural-audit target", "Quantity": "Three-dimensional metal-containing entries", "Count": target, "Fraction": _percent(target, metal_all), "Interpretation": "Operational structure-record denominator"},
        {"Analysis": "Audit records", "Quantity": "Schema-valid emitted records", "Count": emitted, "Fraction": _percent(emitted, target), "Interpretation": "Closed typed audit outcomes"},
        {"Analysis": "Structural records", "Quantity": "Entries with any structural metal site", "Count": structural_entries, "Fraction": _percent(structural_entries, target), "Interpretation": "At least one native-source donor or metal relation"},
        {"Analysis": "Structural records", "Quantity": "Entries with every metal site structural", "Count": all_sites_entries, "Fraction": _percent(all_sites_entries, target), "Interpretation": "No audit-only site in the entry"},
        {"Analysis": "Metal-site records", "Quantity": "All in-scope metal sites", "Count": total_sites, "Fraction": "100.0000%", "Interpretation": "Site-level denominator"},
        {"Analysis": "Metal-site records", "Quantity": "Structural metal sites", "Count": structural_sites, "Fraction": _percent(structural_sites, total_sites), "Interpretation": "Native-source structural relation retained"},
        {"Analysis": "Metal-site records", "Quantity": "Audit-only metal sites", "Count": audit_sites, "Fraction": _percent(audit_sites, total_sites), "Interpretation": "No qualifying native-source structural relation"},
        {"Analysis": "Independent source reread", "Quantity": "Records passing source-signature validation", "Count": int(fidelity["passed"]), "Fraction": _percent(int(fidelity["passed"]), int(fidelity["records"])), "Interpretation": "Agreement under the same frozen metal policy"},
    ]
    _write_csv(
        protected[5],
        ("Analysis", "Quantity", "Count", "Fraction", "Interpretation"),
        table_rows,
    )

    figure_rows = [
        {"Panel": "B", "Metric": "All April 2025 CSD entries", "Count": total_csd, "Denominator": total_csd, "Percent": "100.0000%", "Note": "Frozen release"},
        {"Panel": "B", "Metric": "Metal-containing census", "Count": metal_all, "Denominator": total_csd, "Percent": _percent(metal_all, total_csd), "Note": "Frozen 96-symbol policy"},
        {"Panel": "B", "Metric": "3D structural-audit target", "Count": target, "Denominator": metal_all, "Percent": _percent(target, metal_all), "Note": "Operational target"},
        {"Panel": "B", "Metric": "Nonpolymeric single-metal entries", "Count": int(chemistry["strict_nonpolymeric_single_metal_center_entries"]), "Denominator": target, "Percent": _percent(int(chemistry["strict_nonpolymeric_single_metal_center_entries"]), target), "Note": "Mutually exclusive target class"},
        {"Panel": "B", "Metric": "Nonpolymeric multimetal entries", "Count": int(chemistry["nonpolymeric_multiple_metal_centers_or_external_relation_entries"]), "Denominator": target, "Percent": _percent(int(chemistry["nonpolymeric_multiple_metal_centers_or_external_relation_entries"]), target), "Note": "Mutually exclusive target class"},
        {"Panel": "B", "Metric": "Polymeric entries", "Count": int(chemistry["polymeric_entries"]), "Denominator": target, "Percent": _percent(int(chemistry["polymeric_entries"]), target), "Note": "Mutually exclusive target class"},
        {"Panel": "C", "Metric": "Entries with any structural site", "Count": structural_entries, "Denominator": target, "Percent": _percent(structural_entries, target), "Note": "Entry-level coverage"},
        {"Panel": "C", "Metric": "Entries with every site structural", "Count": all_sites_entries, "Denominator": target, "Percent": _percent(all_sites_entries, target), "Note": "Entry-level all-site coverage"},
        {"Panel": "C", "Metric": "Structural metal sites", "Count": structural_sites, "Denominator": total_sites, "Percent": _percent(structural_sites, total_sites), "Note": "Site-level coverage"},
        {"Panel": "C", "Metric": "Resolved structural-object outcomes", "Count": int(status["EMITTED_RESOLVED"]), "Denominator": target, "Percent": _percent(int(status["EMITTED_RESOLVED"]), target), "Note": "Does not imply exact ID or CShM"},
        {"Panel": "C", "Metric": "Partial outcomes", "Count": int(status["EMITTED_PARTIAL"]), "Denominator": target, "Percent": _percent(int(status["EMITTED_PARTIAL"]), target), "Note": "Explicitly retained"},
        {"Panel": "C", "Metric": "Disorder-ambiguous outcomes", "Count": int(status["EMITTED_AMBIGUOUS"]), "Denominator": target, "Percent": _percent(int(status["EMITTED_AMBIGUOUS"]), target), "Note": "Explicitly retained"},
        {"Panel": "C", "Metric": "All-sites audit-only outcomes", "Count": int(status["EMITTED_AUDIT_ONLY"]), "Denominator": target, "Percent": _percent(int(status["EMITTED_AUDIT_ONLY"]), target), "Note": "No native structural relation at any metal site"},
        {"Panel": "D", "Metric": "Shared-donor entries", "Count": int(chemistry["bridged_entries"]), "Denominator": target, "Percent": _percent(int(chemistry["bridged_entries"]), target), "Note": "Nonexclusive flag"},
        {"Panel": "D", "Metric": "Confirmed haptic/pi entries", "Count": int(chemistry["haptic_entries"]), "Denominator": target, "Percent": _percent(int(chemistry["haptic_entries"]), target), "Note": "Nonexclusive flag"},
        {"Panel": "D", "Metric": "Entries with a resolved nonzero translation edge", "Count": int(chemistry["nonzero_translation_edge_entries"]), "Denominator": target, "Percent": _percent(int(chemistry["nonzero_translation_edge_entries"]), target), "Note": "Nonexclusive flag"},
        {"Panel": "D", "Metric": "Ambiguous collective-pi candidates", "Count": int(chemistry["ambiguous_pi_candidate_entries"]), "Denominator": target, "Percent": _percent(int(chemistry["ambiguous_pi_candidate_entries"]), target), "Note": "No forced hapticity"},
    ]
    _write_csv(
        protected[6],
        ("Panel", "Metric", "Count", "Denominator", "Percent", "Note"),
        figure_rows,
    )

    provenance = {
        "release_kind": "public aggregate evidence; no licensed CSD rows",
        "source_release": census["source_release"],
        "database_sha256": census["database_sha256"],
        "metal_policy_id": census["metal_policy_id"],
        "metal_policy_sha256": census["metal_policy_sha256"],
        "input_artifact_sha256": {name: _sha256(path) for name, path in inputs.items()},
        "internal_ledger_commitments": merged.get("artifacts", {}),
        "packaging_script_sha256": _sha256(Path(__file__).resolve()),
        "policy_probe": probe,
    }
    _write_json(protected[4], provenance)

    readme = f"""# April 2025 CSD metal-containing structural-record audit

This directory contains aggregate public evidence for the frozen
`{merged['protocol_id']}` audit. The April 2025 CSD contained {total_csd:,}
entries. Under the frozen 96-symbol CCDC-derived metal policy, {metal_all:,}
entries are metal-containing; {target:,} have three-dimensional structure and
form the operational structural-audit target, while {metal_no3d:,} remain in
the census without a 3D structural record.

Of the 96 frozen policy symbols, {observed_metal_elements} occur in this CSD
snapshot and {observed_3d_metal_elements} occur in the 3D target. Policy
membership is therefore not evidence that every listed synthetic element was
empirically exercised by this release.

All {target:,} target entries emitted schema-valid audit records. At least one
CSD-native structural relation was retained in {structural_entries:,} entries,
and all in-scope metal sites were structural in {all_sites_entries:,} entries.
At the site level, {structural_sites:,} of {total_sites:,} metal sites were
structural. The independent native-object reread passed for
{int(fidelity['passed']):,}/{int(fidelity['records']):,} emitted records.

These rates describe source-faithful structural transcription. They do not
imply that every site has a unique canonical CoordRep identity, a supported
CShM vector, or a chemically complete distance-derived first sphere. Donor,
bridge, haptic/pi, and periodic relations are derived from the CSD-native
molecular bond graph; no distance-derived contact is added.

Row-level records and outcomes are licensed CSD derivatives and are not
redistributed. Their SHA-256 commitments are retained in the public summaries.
"""
    protected[7].write_text(readme, encoding="utf-8")

    checksum_targets = [path for path in protected[:-1] if path.exists()]
    protected[8].write_text(
        "".join(f"{_sha256(path)}  {path.name}\n" for path in checksum_targets),
        encoding="utf-8",
    )
    print(json.dumps({
        "output_dir": str(args.output_dir),
        "files": [path.name for path in protected],
        "target_entries": target,
        "source_fidelity_passed": int(fidelity["passed"]),
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
