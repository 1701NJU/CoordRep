#!/usr/bin/env python3
"""Integrity and numerical-lock checks for the current JACS release."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RELEASE = ROOT / "release" / "jacs-revision-20260826"
REQUIRED = [
    "README.md",
    "PROVENANCE.json",
    "RELEASE_MANIFEST.json",
    "PUBLIC_ARTIFACT_MANIFEST.json",
    "RELEASE_SHA256SUMS.txt",
    "coordrep/__init__.py",
    "canonicalization/general_invariance/summary.json",
    "canonicalization/legal_orbit_challenge/summary.json",
    "audits/full_csd/ALL_METAL_CLAIM_READY_SUMMARY.json",
    "audits/full_csd/METAL_DOMAIN_CENSUS_PUBLIC.json",
    "audits/full_csd/STRUCTURAL_RECORD_AUDIT_PUBLIC.json",
    "audits/full_csd/SOURCE_FIDELITY_PUBLIC.json",
    "audits/full_csd/FIGURE5_SOURCE.csv",
    "audits/full_csd/TABLE_S9_SOURCE.csv",
    "figures/Figure2/Figure2_Canonicalization_Ties_Specificity_20260822_source_data.csv",
    "figures/Figure3/Figure3_summary_and_internal_QC.csv",
    "figures/Figure5/Figure5_FullRelease_MetalContaining_Audit_20260826.svg",
    "figures/Figure5/Figure5_source_data_all_metal_20260826.csv",
    "figures/Figure5/Figure5_caption_all_metal_20260826.txt",
    "figures/Figure6/Figure6_source_data_20260823.csv",
    "figures/Supplementary/Supplementary_Figure_S1_source_data.csv",
    "figures/Supplementary/Supplementary_Figure_S2_DeltaS_Threshold_Sensitivity_source_data.csv",
    "supporting_information/Supplementary_Methods_1_and_Tables_S9A_S9B_20260826.docx",
    "protocols/CSD_MOF_PERIODIC_CANONICAL_LOCAL_SITE_PROTOCOL_v3.json",
]
BINARY_SUFFIXES = {
    ".docx", ".jpg", ".jpeg", ".pdf", ".png", ".svg", ".tif", ".tiff", ".xlsx", ".zip"
}


def fail(message: str) -> None:
    raise AssertionError(message)


def require_file(relative: str) -> Path:
    path = RELEASE / relative
    if not path.is_file() or path.stat().st_size == 0:
        fail(f"missing or empty: {relative}")
    return path


def load_json(relative: str) -> dict:
    return json.loads(require_file(relative).read_text(encoding="utf-8"))


def csv_rows(relative: str) -> list[dict[str, str]]:
    with require_file(relative).open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def canonical_bytes(path: Path) -> bytes:
    data = path.read_bytes()
    if path.suffix.lower() not in BINARY_SUFFIXES:
        data = data.replace(b"\r\n", b"\n")
    return data


def sha256(path: Path) -> str:
    return hashlib.sha256(canonical_bytes(path)).hexdigest()


def verify_release_checksums() -> int:
    checksum_path = require_file("RELEASE_SHA256SUMS.txt")
    observed_paths: set[str] = set()
    for line in checksum_path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        expected, relative = line.split("  ", 1)
        if relative in observed_paths:
            fail(f"duplicate checksum target: {relative}")
        observed_paths.add(relative)
        path = RELEASE / relative
        if not path.is_file():
            fail(f"checksum target missing: {relative}")
        if sha256(path) != expected:
            fail(f"checksum mismatch: {relative}")
    expected_paths = {
        path.relative_to(RELEASE).as_posix()
        for path in RELEASE.rglob("*")
        if path.is_file()
        and path.name != "RELEASE_SHA256SUMS.txt"
        and "__pycache__" not in path.parts
        and ".pytest_cache" not in path.parts
        and path.suffix.lower() not in {".pyc", ".pyo"}
    }
    if observed_paths != expected_paths:
        missing = sorted(expected_paths - observed_paths)
        extra = sorted(observed_paths - expected_paths)
        fail(f"checksum coverage mismatch; missing={missing[:5]}, extra={extra[:5]}")
    return len(observed_paths)


def verify_nested_checksums() -> int:
    """Verify focused evidence-tree SHA256SUMS files recursively."""
    checked = 0
    for checksum_path in sorted(RELEASE.rglob("SHA256SUMS.txt")):
        seen: set[str] = set()
        for line in checksum_path.read_text(encoding="utf-8").splitlines():
            if not line.strip():
                continue
            expected, relative = line.split("  ", 1)
            if relative in seen:
                fail(f"duplicate nested checksum target: {checksum_path}: {relative}")
            seen.add(relative)
            target = checksum_path.parent / relative
            if not target.is_file():
                fail(f"nested checksum target missing: {target.relative_to(RELEASE)}")
            if sha256(target) != expected:
                fail(f"nested checksum mismatch: {target.relative_to(RELEASE)}")
            checked += 1
    return checked


def main() -> int:
    if not RELEASE.is_dir():
        fail(f"release directory not found: {RELEASE}")
    for relative in REQUIRED:
        require_file(relative)

    manifest = load_json("RELEASE_MANIFEST.json")
    if manifest["release_id"] != "jacs-revision-20260826" or manifest["status"] != "current":
        fail("unexpected release identity/status")
    if manifest["record_language_version"] != "1.1.2rc3":
        fail("unexpected representation version")
    if manifest["audit_protocol"]["runtime"]["ccdc"] != "3.6.0":
        fail("full-CSD runtime must be CCDC Python API 3.6.0")

    audit = manifest["full_csd_audit"]
    expected_scalars = {
        "all_entries": 1371757,
        "metal_containing_entries": 783263,
        "metal_containing_census_only_without_3d": 36752,
        "target_entries": 746511,
        "schema_valid_records": 746511,
    }
    for key, expected in expected_scalars.items():
        if audit[key] != expected:
            fail(f"manifest lock mismatch: {key}")
    if audit["source_fidelity"] != {"passed": 746511, "failed": 0}:
        fail("source-fidelity lock mismatch")
    if audit["entry_coverage"] != {
        "at_least_one_structural_site": 741402,
        "all_metal_sites_structural": 733004,
    }:
        fail("entry-coverage lock mismatch")
    if audit["site_coverage"] != {"structural": 2574081, "audit_only": 34367, "total": 2608448}:
        fail("site-accounting lock mismatch")
    if sum(audit["mutually_exclusive_target_partition"].values()) != 746511:
        fail("target partition does not close")
    if sum(audit["mutually_exclusive_metal_block_signatures"].values()) != 746511:
        fail("metal-block partition does not close")
    if sum(audit["terminal_entry_outcomes"].values()) != 746511:
        fail("terminal outcomes do not close")
    if sum(audit["geometry_fields"].values()) != 2608448:
        fail("geometry-field accounting does not close")

    fidelity = load_json("audits/full_csd/SOURCE_FIDELITY_PUBLIC.json")
    if (fidelity["records"], fidelity["passed"], fidelity["failed"]) != (746511, 746511, 0):
        fail("public source-fidelity evidence mismatch")
    record_audit = load_json("audits/full_csd/STRUCTURAL_RECORD_AUDIT_PUBLIC.json")
    if record_audit["metal_containing_entries"] != 783263:
        fail("public census evidence mismatch")
    if record_audit["entry_outcomes"] != 746511:
        fail("public target evidence mismatch")
    if record_audit["structural_metal_sites"] != 2574081:
        fail("public site evidence mismatch")

    fig5 = csv_rows("figures/Figure5/Figure5_source_data_all_metal_20260826.csv")
    observed5 = {(r["metric"], r["count"], r["denominator"]) for r in fig5}
    for item in {
        ("Metal-containing census", "783263", "1371757"),
        ("Predefined audit corpus", "746511", "783263"),
        ("Independent source-fidelity matches", "746511", "746511"),
        ("Structural metal sites", "2574081", "2608448"),
        ("Shared-donor entries", "213595", "746511"),
        ("Confirmed haptic/pi entries", "113726", "746511"),
        ("Resolved nonzero-translation-edge entries", "129092", "746511"),
    }:
        if item not in observed5:
            fail(f"Figure 5 lock missing: {item}")

    general = load_json("canonicalization/general_invariance/summary.json")
    if not general["full"]["publication_pass"] or general["full"]["coordrep_mismatch_rows"] != 0:
        fail("general canonical-invariance challenge did not pass")
    orbit = load_json("canonicalization/legal_orbit_challenge/summary.json")
    if orbit["full_exact_coordrep"]["exact_variant_matches"] != 100000:
        fail("legal-orbit exact-match lock mismatch")

    nested_checked = verify_nested_checksums()
    checked = verify_release_checksums()
    print(f"PASS: {RELEASE.relative_to(ROOT)}")
    print("PASS: all-metal CSD numerical lock and Figure 5 source table")
    print(f"PASS: {nested_checked} nested evidence files covered by SHA-256")
    print(f"PASS: {checked} public release files covered by SHA-256")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
