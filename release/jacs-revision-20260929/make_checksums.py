"""Write or verify the addendum's SHA-256 file manifest."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parent
MANIFEST = ROOT / "SHA256SUMS.txt"


def digest(path: Path) -> str:
    hasher = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


def expected_lines() -> list[str]:
    files = sorted(
        path for path in ROOT.rglob("*")
        if path.is_file()
        and path != MANIFEST
        and "__pycache__" not in path.parts
        and ".pytest_cache" not in path.parts
    )
    return [f"{digest(path)}  {path.relative_to(ROOT).as_posix()}" for path in files]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true", help="regenerate SHA256SUMS.txt")
    args = parser.parse_args()
    expected = expected_lines()
    if args.write:
        MANIFEST.write_bytes(("\n".join(expected) + "\n").encode("utf-8"))
        print(f"Wrote {len(expected)} file digests to {MANIFEST}")
    else:
        actual = MANIFEST.read_text(encoding="utf-8").splitlines()
        if actual != expected:
            raise SystemExit("SHA256SUMS.txt is incomplete or contains a mismatch")
        print(f"PASS: {len(expected)} addendum file digests verified")


if __name__ == "__main__":
    main()
