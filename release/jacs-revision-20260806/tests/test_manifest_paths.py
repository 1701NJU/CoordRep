"""Validate that the current release map points to existing artifacts."""

import json
from pathlib import Path


RELEASE = Path(__file__).resolve().parent.parent
ROOT = RELEASE.parent.parent


def test_release_manifest_exists_and_is_current():
    manifest = RELEASE / "RELEASE_MANIFEST.json"
    assert manifest.exists() and manifest.stat().st_size > 100
    payload = json.loads(manifest.read_text(encoding="utf-8"))
    assert payload["release_id"] == "jacs-revision-20260806"
    assert payload["record_language_version"] == "1.1.2rc2"


def test_figure_map_targets_exist():
    figure_map = RELEASE / "figures" / "FIGURE_MAP.md"
    assert figure_map.exists()
    for figure in range(1, 7):
        assert (RELEASE / "figures" / f"Figure{figure}").is_dir()


def test_current_periodic_protocol_exists():
    protocol = RELEASE / "protocols" / "CSD_MOF_PERIODIC_CANONICAL_LOCAL_SITE_PROTOCOL_v3.json"
    payload = json.loads(protocol.read_text(encoding="utf-8"))
    assert payload["protocol_version"] == "CSD_MOF_PERIODIC_CANONICAL_LOCAL_SITE_v3"


def test_no_raw_csd_files_in_current_release():
    raw_exts = {".cif", ".mol", ".mol2", ".xyz", ".gjf", ".hkl", ".res", ".fcf"}
    found = [p for p in RELEASE.rglob("*") if p.is_file() and p.suffix.lower() in raw_exts]
    assert not found, "raw coordinate files found: " + ", ".join(str(p) for p in found[:10])
