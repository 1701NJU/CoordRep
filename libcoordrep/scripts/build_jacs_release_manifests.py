#!/usr/bin/env python3
"""Build the public artifact manifest and whole-release SHA-256 lock."""

from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RELEASE = ROOT / "release" / "jacs-revision-20260819"
ARTIFACT_MANIFEST = RELEASE / "PUBLIC_ARTIFACT_MANIFEST.json"
CHECKSUMS = RELEASE / "RELEASE_SHA256SUMS.txt"
EXCLUDED_NAMES = {ARTIFACT_MANIFEST.name, CHECKSUMS.name}
EXCLUDED_PARTS = {"__pycache__", ".pytest_cache"}
EXCLUDED_SUFFIXES = {".pyc", ".pyo"}
BINARY_SUFFIXES = {".pdf", ".png", ".svg", ".tif", ".tiff"}


def canonical_bytes(path: Path) -> bytes:
    """Return raw binary bytes or LF-normalized text bytes.

    Git may materialize text files with CRLF in a Windows checkout even though
    GitHub archives contain LF.  The public checksum contract therefore uses
    canonical LF bytes for text and untouched bytes for publication artwork.
    """
    if path.suffix.lower() in BINARY_SUFFIXES:
        # Binary artwork is stored with -text in .gitattributes. Read the
        # staged blob when Git is available so a racy Windows worktree stat or
        # an editor's newline rewrite cannot alter the published lock.
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


def public_files(*, include_artifact_manifest: bool) -> list[Path]:
    files = []
    for path in RELEASE.rglob("*"):
        if not path.is_file():
            continue
        relative = path.relative_to(RELEASE)
        if EXCLUDED_PARTS.intersection(relative.parts) or path.suffix.lower() in EXCLUDED_SUFFIXES:
            continue
        if path.name == CHECKSUMS.name:
            continue
        if not include_artifact_manifest and path.name == ARTIFACT_MANIFEST.name:
            continue
        files.append(path)
    return sorted(files, key=lambda p: p.relative_to(RELEASE).as_posix())


def role(relative: str) -> str:
    if relative.startswith("figures/"):
        return "figure artwork, caption, or public source data"
    if relative.startswith("audits/full_csd/"):
        return "full-CSD aggregate evidence"
    if relative.startswith("coordrep/") or relative.startswith("coordrep_tools/"):
        return "CoordRep implementation"
    if relative.startswith("tests/"):
        return "regression test or small fixture"
    if relative.startswith("protocols/"):
        return "frozen protocol"
    if relative.startswith("models/") or relative.startswith("brain/"):
        return "model implementation or compact public result"
    return "release documentation or supporting source"


def main() -> None:
    if not RELEASE.is_dir():
        raise SystemExit(f"release directory not found: {RELEASE}")

    artifacts = []
    for path in public_files(include_artifact_manifest=False):
        relative = path.relative_to(RELEASE).as_posix()
        artifacts.append(
            {
                "path": relative,
                "bytes": len(canonical_bytes(path)),
                "sha256": sha256(path),
                "role": role(relative),
                "publication_status": "current" if not relative.startswith("archive/retired/") else "historical/superseded",
            }
        )

    payload = {
        "release_id": "jacs-revision-20260819",
        "manifest_date": "2026-08-19",
        "algorithm": "SHA-256 over canonical-LF text bytes and raw binary bytes",
        "scope": "all public release files except this manifest and RELEASE_SHA256SUMS.txt",
        "n_files": len(artifacts),
        "total_bytes": sum(item["bytes"] for item in artifacts),
        "files": artifacts,
    }
    ARTIFACT_MANIFEST.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
        newline="\n",
    )

    lines = []
    for path in public_files(include_artifact_manifest=True):
        relative = path.relative_to(RELEASE).as_posix()
        lines.append(f"{sha256(path)}  {relative}")
    CHECKSUMS.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")

    print(f"Wrote {ARTIFACT_MANIFEST.relative_to(ROOT)} ({len(artifacts)} artifacts)")
    print(f"Wrote {CHECKSUMS.relative_to(ROOT)} ({len(lines)} checksums)")


if __name__ == "__main__":
    main()
