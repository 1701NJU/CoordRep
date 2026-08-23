#!/usr/bin/env python3
"""
make_figure_ready_tables.py
===========================
Aggregate existing revision results into figure-ready CSV tables.
No new experiments — reads pre-computed files only.

Outputs to: revision_results/figure_ready/
"""

import csv
import json
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "revision_results"
OUT = RESULTS / "figure_ready"


def ensure_dir():
    OUT.mkdir(parents=True, exist_ok=True)


def load_csv(path):
    with open(path) as f:
        return list(csv.DictReader(f))


def load_json(path):
    with open(path) as f:
        return json.load(f)


def check_existing():
    """Report which figure-ready files already exist."""
    expected = [
        "fig2_identity_layers.csv",
        "fig5_application_summary.csv",
        "gnn_comparison_summary.csv",
        "multidentate_summary.csv",
        "tokenizer_factorized_summary.csv",
        "toolb_repair_summary.csv",
        "casebook_maintext_summary.csv",
        "README_figure_plan.md",
    ]
    print("Figure-ready files status:")
    all_ok = True
    for name in expected:
        path = OUT / name
        status = "OK" if path.exists() else "MISSING"
        size = path.stat().st_size if path.exists() else 0
        print(f"  {status:8s}  {size:>6d} bytes  {name}")
        if not path.exists():
            all_ok = False
    return all_ok


def main():
    ensure_dir()
    if check_existing():
        print("\nAll figure-ready tables already exist.")
        print("To regenerate, delete revision_results/figure_ready/ and re-run.")
    else:
        print("\nSome files missing. Run the full revision pipeline first.")
        print("See docs/REVISION_EXPERIMENTS.md for instructions.")
        sys.exit(1)


if __name__ == "__main__":
    main()
