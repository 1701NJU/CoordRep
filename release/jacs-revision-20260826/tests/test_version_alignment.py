"""Release metadata must agree between the package and pyproject."""

from __future__ import annotations

import tomllib
from pathlib import Path

import coordrep


def test_runtime_version_matches_pyproject() -> None:
    root = Path(__file__).resolve().parents[1]
    pyproject = tomllib.loads((root / "pyproject.toml").read_text(encoding="utf-8"))
    assert coordrep.__version__ == pyproject["project"]["version"]
