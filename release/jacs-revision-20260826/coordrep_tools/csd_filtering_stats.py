"""
csd_filtering_stats.py
======================
Waterfall statistics tracker for CSD filtering pipeline.
Accumulates counts at each stage and produces summary tables.
"""

from __future__ import annotations

import csv
import json
import os
from collections import Counter, defaultdict
from typing import Dict, List


class WaterfallTracker:
    """Track pass/fail counts through the filtering waterfall."""

    def __init__(self):
        self.n_scanned = 0
        self.n_has_3d = 0
        self.n_transition_metal = 0
        self.n_mononuclear = 0
        self.n_eta1 = 0          # no hapticity
        self.n_cn_ok = 0
        self.n_no_disorder = 0
        self.n_no_polymer = 0
        self.n_donor_ok = 0
        self.n_valid_smiles = 0
        self.n_valid_coordrep = 0

        self.rejection_reasons = Counter()
        self.retained_by_metal = Counter()
        self.retained_by_cn = Counter()
        self.boundary_count = 0
        self.total_retained = 0

    def record_scan(self):
        self.n_scanned += 1

    def record_filter_result(self, passed: bool, reason: str,
                             metal: str = "", cn: int = 0):
        if not passed:
            self.rejection_reasons[reason] += 1
        else:
            self.retained_by_metal[metal] += 1
            self.retained_by_cn[cn] += 1
            self.total_retained += 1

    def record_waterfall_stage(self, stage: str):
        """Increment the counter for a named stage."""
        attr = f"n_{stage}"
        if hasattr(self, attr):
            setattr(self, attr, getattr(self, attr) + 1)

    def record_boundary(self):
        self.boundary_count += 1

    def summary_dict(self) -> dict:
        top_rej = self.rejection_reasons.most_common(10)
        return {
            "n_csd_scanned": self.n_scanned,
            "n_has_3d": self.n_has_3d,
            "n_transition_metal": self.n_transition_metal,
            "n_mononuclear": self.n_mononuclear,
            "n_eta1": self.n_eta1,
            "n_cn_ok": self.n_cn_ok,
            "n_no_disorder": self.n_no_disorder,
            "n_no_polymer": self.n_no_polymer,
            "n_donor_ok": self.n_donor_ok,
            "n_valid_smiles": self.n_valid_smiles,
            "n_valid_coordrep": self.n_valid_coordrep,
            "retention_total": self.total_retained,
            "boundary_fraction": round(self.boundary_count / max(self.total_retained, 1), 4),
            "top_rejection_reasons": {k: v for k, v in top_rej},
        }

    def write_all(self, out_dir: str):
        """Write all CSV/JSON outputs."""
        os.makedirs(out_dir, exist_ok=True)

        # summary.json
        with open(os.path.join(out_dir, "summary.json"), "w") as f:
            json.dump(self.summary_dict(), f, indent=2)

        # csd_filter_waterfall.csv
        stages = [
            ("csd_scanned", self.n_scanned),
            ("has_3d", self.n_has_3d),
            ("transition_metal", self.n_transition_metal),
            ("mononuclear", self.n_mononuclear),
            ("eta1_only", self.n_eta1),
            ("cn_in_range", self.n_cn_ok),
            ("no_disorder", self.n_no_disorder),
            ("no_polymer", self.n_no_polymer),
            ("donor_ok", self.n_donor_ok),
            ("valid_smiles", self.n_valid_smiles),
            ("valid_coordrep", self.n_valid_coordrep),
        ]
        with open(os.path.join(out_dir, "csd_filter_waterfall.csv"), "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["stage", "count", "retention_pct"])
            for name, cnt in stages:
                pct = round(100 * cnt / max(self.n_scanned, 1), 2)
                w.writerow([name, cnt, pct])

        # csd_rejection_reasons.csv
        with open(os.path.join(out_dir, "csd_rejection_reasons.csv"), "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["reason", "count", "pct_of_scanned"])
            for reason, cnt in self.rejection_reasons.most_common():
                pct = round(100 * cnt / max(self.n_scanned, 1), 4)
                w.writerow([reason, cnt, pct])

        # csd_retained_by_metal.csv
        with open(os.path.join(out_dir, "csd_retained_by_metal.csv"), "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["metal", "count"])
            for m, cnt in self.retained_by_metal.most_common():
                w.writerow([m, cnt])

        # csd_retained_by_cn.csv
        with open(os.path.join(out_dir, "csd_retained_by_cn.csv"), "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["cn", "count"])
            for cn, cnt in sorted(self.retained_by_cn.items()):
                w.writerow([cn, cnt])

        # csd_boundary_fraction.csv
        with open(os.path.join(out_dir, "csd_boundary_fraction.csv"), "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["metric", "value"])
            w.writerow(["total_retained", self.total_retained])
            w.writerow(["boundary_count", self.boundary_count])
            w.writerow(["boundary_fraction",
                         round(self.boundary_count / max(self.total_retained, 1), 4)])
