"""Test that all paths referenced in DATA_MANIFEST.md exist."""

import re
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent.parent  # coordrep-release/
LIB_ROOT = REPO_ROOT / "libcoordrep"
MANIFEST = REPO_ROOT / "DATA_MANIFEST.md"


def _extract_manifest_paths():
    """Extract file paths from DATA_MANIFEST.md table rows."""
    if not MANIFEST.exists():
        pytest.skip("DATA_MANIFEST.md not found")

    text = MANIFEST.read_text()
    paths = set()
    # Match paths in table cells: `revision_results/...` or `scripts/...`
    for match in re.finditer(
        r"`((?:revision_results|scripts)/[^`]+)`", text
    ):
        p = match.group(1)
        # Skip directory-only references and wildcards
        if p.endswith("/") or "*" in p:
            continue
        paths.add(p)
    # Checkpoint paths are at repo root, not under libcoordrep
    for match in re.finditer(r"`(checkpoints/[^`]+)`", text):
        p = match.group(1)
        if p.endswith("/") or "*" in p:
            continue
        full = REPO_ROOT / p
        if not full.exists():
            # Not a test failure — checkpoints may be large/external
            pass
        paths.add("../" + p)  # relative from libcoordrep
    return sorted(paths)


class TestManifestPaths:
    """Every file path in DATA_MANIFEST.md must exist on disk."""

    @pytest.fixture(scope="class")
    def manifest_paths(self):
        return _extract_manifest_paths()

    def test_manifest_exists(self):
        assert MANIFEST.exists(), "DATA_MANIFEST.md not found at repo root"

    def test_manifest_not_empty(self):
        assert MANIFEST.stat().st_size > 1000, "DATA_MANIFEST.md too small"

    def test_all_paths_exist(self, manifest_paths):
        missing = []
        for p in manifest_paths:
            full = LIB_ROOT / p
            if not full.exists():
                missing.append(p)
        assert not missing, (
            f"{len(missing)} manifest paths missing:\n"
            + "\n".join(f"  - {m}" for m in missing[:20])
        )

    def test_no_raw_csd_files(self):
        raw_exts = {".cif", ".hkl", ".res", ".fcf"}
        found = []
        for ext in raw_exts:
            found.extend(LIB_ROOT.rglob(f"*{ext}"))
        assert not found, (
            f"Raw CSD files found:\n"
            + "\n".join(str(f) for f in found[:10])
        )
