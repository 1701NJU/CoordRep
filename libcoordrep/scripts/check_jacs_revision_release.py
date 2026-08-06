#!/usr/bin/env python3
"""Integrity checks for the current JACS revision release.

This checker intentionally scopes itself to ``release/jacs-revision-20260806``.
The older ``libcoordrep/revision_results`` tree is retained as provenance and
contains historical outputs that are not part of the current evidence package.
"""

from __future__ import annotations

import csv
import json
import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RELEASE = ROOT / "release" / "jacs-revision-20260806"
REQUIRED = [
    "README.md",
    "RELEASE_MANIFEST.json",
    "pyproject.toml",
    "coordrep/__init__.py",
    "figures/Figure4/Figure4_Compositional_coordination_states_JACS_20260806_source_data.csv",
    "figures/Figure5/Figure5_CoordRep_current_CSD_source_data_20260802.csv",
    "figures/Figure6/Figure6_source_data_20260806.csv",
    "protocols/CSD_MOF_PERIODIC_CANONICAL_LOCAL_SITE_PROTOCOL_v3.json",
]
RAW_SUFFIXES = {".cif", ".mol", ".mol2", ".xyz", ".gjf", ".hkl", ".res", ".fcf"}


def fail(message: str) -> None:
    raise AssertionError(message)


def require_file(relative: str) -> Path:
    path = RELEASE / relative
    if not path.is_file() or path.stat().st_size == 0:
        fail(f"missing or empty: {relative}")
    return path


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def main() -> int:
    if not RELEASE.is_dir():
        fail(f"release directory not found: {RELEASE}")
    for relative in REQUIRED:
        require_file(relative)

    manifest = json.loads((RELEASE / "RELEASE_MANIFEST.json").read_text(encoding="utf-8"))
    if manifest["release_id"] != "jacs-revision-20260806":
        fail("unexpected release_id")
    if manifest["record_language_version"] != "1.1.2rc2":
        fail("unexpected record-language version")

    init_text = (RELEASE / "coordrep/__init__.py").read_text(encoding="utf-8")
    if '"1.1.2rc2"' not in init_text:
        fail("CoordRep package version is not 1.1.2rc2")

    fig6 = rows(RELEASE / "figures/Figure6/Figure6_source_data_20260806.csv")
    expected_fig6 = {
        ("processing success", "15905"),
        ("state-bearing entries", "10948"),
        ("states with nonzero translation edge", "74354"),
        ("unique metric IDs", "36552"),
        ("unique stereo IDs", "60139"),
    }
    observed_fig6 = {(row.get("metric", ""), row.get("numerator", "")) for row in fig6}
    if not expected_fig6 <= observed_fig6:
        fail(f"Figure 6 current-v3 lock not found: {sorted(expected_fig6 - observed_fig6)}")

    fig5 = rows(RELEASE / "figures/Figure5/Figure5_CoordRep_current_CSD_source_data_20260802.csv")
    observed_fig5 = {(row.get("key", ""), row.get("value", "")) for row in fig5}
    for key, value in {
        ("all_csd", "1371757"),
        ("tm_3d", "602116"),
        ("strict_core", "101878"),
        ("ligand_topology", "84453"),
        ("linked", "52760"),
        ("haptic", "38203"),
    }:
        if not any(k == key and v.startswith(value) for k, v in observed_fig5):
            fail(f"Figure 5 current CSD value missing: {key}={value}")

    fig4_text = (RELEASE / "figures/Figure4/Figure4_Compositional_coordination_states_JACS_20260806_source_data.csv").read_text(encoding="utf-8")
    for model in ("gine_wide", "schnet_3d", "visnet_3d"):
        if model not in fig4_text.lower():
            fail(f"Figure 4 source table does not contain model label {model}")

    forbidden = []
    for path in RELEASE.rglob("*"):
        if path.is_file() and (path.suffix.lower() in RAW_SUFFIXES or "internal_licensed" in path.name.lower()):
            forbidden.append(path.relative_to(ROOT).as_posix())
    if forbidden:
        fail("redistribution-prohibited files in current release: " + ", ".join(forbidden[:10]))

    print(f"PASS: current JACS release integrity ({RELEASE.relative_to(ROOT)})")
    print("PASS: CoordRep 1.1.2rc2; current Figure 4/5/6 source locks; no raw CSD files")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (AssertionError, KeyError, json.JSONDecodeError) as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        raise SystemExit(1)
