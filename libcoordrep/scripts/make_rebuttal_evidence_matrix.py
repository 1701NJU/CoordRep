#!/usr/bin/env python3
"""
make_rebuttal_evidence_matrix.py
================================
Verify and report on the rebuttal evidence matrix files.
No new experiments — reads pre-computed files only.

Outputs to: revision_results/rebuttal/
"""

import csv
import json
import os
import sys
from pathlib import Path
from collections import Counter

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "revision_results"
OUT = RESULTS / "rebuttal"


def check_existing():
    """Report which rebuttal files exist and their contents."""
    expected = [
        "evidence_matrix.csv",
        "evidence_matrix.md",
        "response_letter_bullets.md",
        "manuscript_revision_map.md",
        "reviewer_priority_summary.json",
    ]
    print("Rebuttal files status:")
    all_ok = True
    for name in expected:
        path = OUT / name
        status = "OK" if path.exists() else "MISSING"
        size = path.stat().st_size if path.exists() else 0
        print(f"  {status:8s}  {size:>6d} bytes  {name}")
        if not path.exists():
            all_ok = False
    return all_ok


def summarize_matrix():
    """Print summary statistics from evidence_matrix.csv."""
    path = OUT / "evidence_matrix.csv"
    if not path.exists():
        return
    with open(path) as f:
        rows = list(csv.DictReader(f))
    print(f"\nEvidence matrix: {len(rows)} comments")
    reviewer_counts = Counter(r['reviewer'] for r in rows)
    status_counts = Counter(r['status'] for r in rows)
    severity_counts = Counter(r['severity'] for r in rows)
    print(f"  By reviewer: {dict(reviewer_counts)}")
    print(f"  By status:   {dict(status_counts)}")
    print(f"  By severity: {dict(severity_counts)}")


def summarize_priority():
    """Print summary from reviewer_priority_summary.json."""
    path = OUT / "reviewer_priority_summary.json"
    if not path.exists():
        return
    with open(path) as f:
        d = json.load(f)
    print(f"\nPriority summary:")
    for k, v in d.items():
        print(f"  {k}: {len(v)} items")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    if check_existing():
        print("\nAll rebuttal files exist.")
        summarize_matrix()
        summarize_priority()
    else:
        print("\nSome files missing. See previous task outputs.")
        sys.exit(1)


if __name__ == "__main__":
    main()
