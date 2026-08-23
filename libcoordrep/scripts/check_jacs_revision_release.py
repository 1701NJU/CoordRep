#!/usr/bin/env python3
"""Integrity, scope, and numerical-lock checks for the current JACS release."""

from __future__ import annotations

import csv
import hashlib
import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RELEASE = ROOT / "release" / "jacs-revision-20260823"
REQUIRED = [
    "README.md",
    "PROVENANCE.json",
    "RELEASE_MANIFEST.json",
    "PUBLIC_ARTIFACT_MANIFEST.json",
    "RELEASE_SHA256SUMS.txt",
    "coordrep/__init__.py",
    "coordrep/audit/csd_adapter.py",
    "coordrep/audit/records.py",
    "canonicalization/general_invariance/summary.json",
    "canonicalization/legal_orbit_challenge/summary.json",
    "canonicalization/expanded_record_regressions/summary.json",
    "audits/full_csd/FULL_CSD_CLAIM_READY_SUMMARY.json",
    "audits/full_csd/SOURCE_FIDELITY_PUBLIC_SUMMARY.json",
    "figures/Figure1/Figure1_Cisplatin_3C_20260822.svg",
    "figures/Figure2/Figure2_Canonicalization_Ties_Specificity_20260822_source_data.csv",
    "figures/Figure3/Figure3_summary_and_internal_QC.csv",
    "figures/Figure4/Figure4_Compositional_coordination_states_JACS_20260806_source_data.csv",
    "figures/Figure5/Figure5_ScopeAware_source_data_20260818.csv",
    "figures/Figure5/Figure5_ScopeAware_caption_20260818.txt",
    "figures/Figure5/Figure5_Full_CSD_Audit_20260823.svg",
    "figures/Figure6/Figure6_source_data_20260823.csv",
    "figures/Supplementary/Supplementary_Figure_S1_DeltaS_Threshold_Sensitivity_source_data.csv",
    "figures/Supplementary/Supplementary_Figure_S2_source_data.csv",
    "supporting_information/Supplementary Information - synchronized 20260823.docx",
    "supporting_information/Supplementary Information - synchronized 20260823.pdf",
    "protocols/CSD_MOF_PERIODIC_CANONICAL_LOCAL_SITE_PROTOCOL_v3.json",
]
RAW_SUFFIXES = {".cif", ".mol", ".mol2", ".xyz", ".gjf", ".hkl", ".res", ".fcf"}
BINARY_SUFFIXES = {
    ".docx",
    ".jpg",
    ".jpeg",
    ".pdf",
    ".png",
    ".svg",
    ".tif",
    ".tiff",
    ".xlsx",
    ".zip",
}


def fail(message: str) -> None:
    raise AssertionError(message)


def require_file(relative: str) -> Path:
    path = RELEASE / relative
    if not path.is_file() or path.stat().st_size == 0:
        fail(f"missing or empty: {relative}")
    return path


def csv_rows(relative: str) -> list[dict[str, str]]:
    with require_file(relative).open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def canonical_bytes(path: Path) -> bytes:
    if path.suffix.lower() in BINARY_SUFFIXES:
        relative = path.relative_to(ROOT).as_posix()
        try:
            return subprocess.check_output(
                ["git", "show", f":{relative}"],
                cwd=ROOT,
                stderr=subprocess.DEVNULL,
            )
        except (FileNotFoundError, subprocess.CalledProcessError):
            return path.read_bytes()
    return path.read_bytes().replace(b"\r\n", b"\n")


def sha256(path: Path) -> str:
    return hashlib.sha256(canonical_bytes(path)).hexdigest()


def verify_release_checksums() -> int:
    checksum_path = require_file("RELEASE_SHA256SUMS.txt")
    checked = 0
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
        observed = sha256(path)
        if observed != expected:
            fail(f"checksum mismatch: {relative}: {observed} != {expected}")
        checked += 1
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
    if checked < 100:
        fail(f"unexpectedly small checksum manifest: {checked} files")
    return checked


def main() -> int:
    if not RELEASE.is_dir():
        fail(f"release directory not found: {RELEASE}")
    for relative in REQUIRED:
        require_file(relative)

    manifest = json.loads(require_file("RELEASE_MANIFEST.json").read_text(encoding="utf-8"))
    if manifest["release_id"] != "jacs-revision-20260823" or manifest["status"] != "current":
        fail("unexpected release identity/status")
    if manifest["record_language_version"] != "1.1.2rc3":
        fail("unexpected representation version")
    if manifest["audit_protocol"]["runtime"]["ccdc"] != "3.6.0":
        fail("full-CSD runtime must be CCDC Python API 3.6.0")
    if "Sc–Zn, Y–Cd, La, and Hf–Hg" not in manifest["audit_protocol"]["target_definition"]:
        fail("frozen in-domain element set is not explicit")

    audit = manifest["full_csd_audit"]
    if (audit["target_entries"], audit["schema_valid_records"]) != (602116, 602116):
        fail("manifest full-CSD emission lock mismatch")
    if audit["entry_coverage"]["all_in_domain_sites_structural"] != 601402:
        fail("manifest all-sites structural lock mismatch")
    if audit["site_coverage"] != {"structural": 1984062, "audit_only": 2136, "total": 1986198}:
        fail("manifest site accounting mismatch")
    for key, value in {
        "in_domain_shared_donor_entries": 147348,
        "confirmed_haptic_pi_entries": 102241,
        "resolved_nonzero_translation_edge_entries": 94452,
    }.items():
        if audit["nonexclusive_relations"][key] != value:
            fail(f"manifest relation lock mismatch: {key}")

    claim = json.loads(require_file("audits/full_csd/FULL_CSD_CLAIM_READY_SUMMARY.json").read_text(encoding="utf-8"))
    if claim["emitted_entries"] != 602116 or claim["failed_entries"] != 0:
        fail("claim-ready emission lock mismatch")
    if claim["coverage"]["entry_all_metal_sites_structural"]["numerator"] != 601402:
        fail("claim-ready all-sites structural lock mismatch")
    if claim["coverage"]["metal_site_structural"]["numerator"] != 1984062:
        fail("claim-ready structural-site lock mismatch")
    if claim["chemical_scope"]["bridged_entries"] != 147348:
        fail("claim-ready shared-donor lock mismatch")
    if claim["chemical_scope"]["haptic_entries"] != 102241:
        fail("claim-ready haptic lock mismatch")

    fidelity = json.loads(require_file("audits/full_csd/SOURCE_FIDELITY_PUBLIC_SUMMARY.json").read_text(encoding="utf-8"))
    if (fidelity["passed"], fidelity["failed"]) != (602116, 0):
        fail("source-fidelity lock mismatch")

    fig2 = csv_rows("figures/Figure2/Figure2_Canonicalization_Ties_Specificity_20260822_source_data.csv")
    observed2 = {(r["panel"], r["record"], r["quantity"], r["value"]) for r in fig2}
    for item in {
        ("B", "CoordRep-State · all arms", "exact variant matches", "300000"),
        ("C", "real legal-orbit screen", "challenge-eligible records", "5173"),
        ("C", "invariant-key sort only", "exact variant matches", "64834"),
        ("C", "full exact legal-tie search", "exact variant matches", "100000"),
        ("C", "non-automorphic attachment control", "false merges", "0"),
        ("D", "idealized fac/mer control", "fac relation inventory", "3 × N-Cl trans"),
        ("D", "idealized fac/mer control", "mer relation inventory", "N-N + N-Cl + Cl-Cl trans"),
    }:
        if item not in observed2:
            fail(f"Figure 2 lock missing: {item}")

    general = json.loads(require_file("canonicalization/general_invariance/summary.json").read_text(encoding="utf-8"))
    if not general["full"]["publication_pass"] or general["full"]["coordrep_mismatch_rows"] != 0:
        fail("general canonical-invariance challenge did not pass")
    orbit = json.loads(require_file("canonicalization/legal_orbit_challenge/summary.json").read_text(encoding="utf-8"))
    if orbit["source_population"]["challenge_eligible_records"] != 5173:
        fail("legal-orbit screening lock mismatch")
    if orbit["full_exact_coordrep"]["exact_variant_matches"] != 100000:
        fail("legal-orbit exact-match lock mismatch")

    fig3 = csv_rows("figures/Figure3/Figure3_summary_and_internal_QC.csv")
    observed3 = {r["scope"]: r for r in fig3}
    if observed3["CN4-CN6"]["n"] != "33863" or observed3["CN4-CN6"]["boundary_margin_lt_1_n"] != "1026":
        fail("Figure 3 cohort lock mismatch")

    fig5 = csv_rows("figures/Figure5/Figure5_ScopeAware_source_data_20260818.csv")
    observed5 = {(row["metric"], row["numerator"], row["denominator"]) for row in fig5}
    for item in {
        ("schema-valid emission", "602116", "602116"),
        ("source-fidelity matches", "602116", "602116"),
        ("all sites structural", "601402", "602116"),
        ("native structural metal sites", "1984062", "1986198"),
        ("in-domain shared donor", "147348", "602116"),
        ("explicit haptic/pi site", "102241", "602116"),
        ("resolved nonzero translation edge", "94452", "602116"),
    }:
        if item not in observed5:
            fail(f"Figure 5 lock missing: {item}")
    caption5 = require_file("figures/Figure5/Figure5_ScopeAware_caption_20260818.txt").read_text(encoding="utf-8")
    for phrase in ("CuN₄Cl₂", "τ = (0, 0, 0)", "τ = (0, 0, 1)", "15,905/15,906"):
        if phrase not in caption5:
            fail(f"Figure 5 caption lock missing: {phrase}")

    fig6 = csv_rows("figures/Figure6/Figure6_source_data_20260823.csv")
    observed6 = {(r["series"], r["category"], r["value"]) for r in fig6}
    for item in {
        ("mean pairwise agreement", "L3", "97.370813"),
        ("XEYVEC hierarchy", "L0", "5"),
        ("XEYVEC hierarchy", "L3", "1"),
        ("identity diversity", "metric IDs", "36552"),
        ("identity diversity", "stereo IDs", "60139"),
        ("canonical-invariance audit", "entries", "400"),
        ("canonical-invariance audit", "sites", "5920"),
    }:
        if item not in observed6:
            fail(f"Figure 6 lock missing: {item}")

    s2 = csv_rows("figures/Supplementary/Supplementary_Figure_S2_source_data.csv")
    observed_s2 = {(r["panel"], r["series"], r["category"], r["value"]) for r in s2}
    for item in {
        ("A", "legal-orbit population", "2", "2652"),
        ("A", "legal-orbit population", "120", "98"),
        ("B", "Full exact search", "Median", "0.7409999961964786"),
        ("B", "Full exact search", "Maximum", "63.83759999880567"),
    }:
        if item not in observed_s2:
            fail(f"Supplementary Figure S2 lock missing: {item}")

    fig4_text = require_file(
        "figures/Figure4/Figure4_Compositional_coordination_states_JACS_20260806_source_data.csv"
    ).read_text(encoding="utf-8").lower()
    for model in ("gine_wide", "schnet_3d", "visnet_3d"):
        if model not in fig4_text:
            fail(f"Figure 4 source table does not contain model label {model}")

    forbidden = []
    for path in RELEASE.rglob("*"):
        if path.is_file() and (
            path.suffix.lower() in RAW_SUFFIXES or "internal_licensed" in path.name.lower()
        ):
            forbidden.append(path.relative_to(ROOT).as_posix())
    if forbidden:
        fail("redistribution-prohibited files: " + ", ".join(forbidden[:10]))

    artifact_manifest = json.loads(require_file("PUBLIC_ARTIFACT_MANIFEST.json").read_text(encoding="utf-8"))
    if artifact_manifest["release_id"] != "jacs-revision-20260823":
        fail("artifact-manifest release ID mismatch")
    if artifact_manifest["algorithm"] != "SHA-256 over canonical-LF text bytes and raw binary bytes":
        fail("unexpected artifact-checksum policy")
    checked = verify_release_checksums()

    print(f"PASS: current JACS release integrity ({RELEASE.relative_to(ROOT)})")
    print("PASS: rc3 canonicalization locks and current Figure 1-6/S1-S2 evidence")
    print("PASS: full-CSD 602116/602116 source fidelity; 601402 all-sites structural entries")
    print("PASS: no licensed raw CSD files in the public package")
    print(f"PASS: {checked} release files verified by SHA-256")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (AssertionError, KeyError, json.JSONDecodeError, ValueError) as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        raise SystemExit(1)
