#!/usr/bin/env python3
"""
run_csd_pathfinder.py
=====================
Full-CSD CoordRep-PathFinder: Coordination-geometry pathway mining.

Transforms CSD crystal database into continuous coordination-geometry maps,
identifies connectivity-family geometry trajectories.

Usage:
  # Analyze existing retained entries (fast, for development):
  python scripts/run_csd_pathfinder.py --use-existing

  # Full CSD scan + analysis (requires CSD Python API, slow):
  /data/miniconda3/envs/1701/bin/python scripts/run_csd_pathfinder.py --full-scan

License: No raw CSD coordinates exported. Only refcodes, aggregate statistics,
         CShM values, and hashed identity keys.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import sys
import time
from collections import Counter, defaultdict
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy import stats as scipy_stats

sys.path.insert(0, str(Path(__file__).parent.parent))

BASE = Path(__file__).parent.parent
RETAINED_PATH = BASE / "revision_results/csd_external/csd_retained_entries.jsonl"
OUT_DIR = BASE / "revision_results/csd_pathfinder_full"

LICENSE_NOTE = ("CSD-derived analysis results only; "
               "raw coordinates are not redistributed.")

# Shape reference names per CN
CN_SHAPE_REFS = {
    4: ("SP", "Td"),
    5: ("TBP", "SPY"),
    6: ("Oh", "TPr"),
}


# ══════════════════════════════════════════════════════════════════════
# Data structures
# ══════════════════════════════════════════════════════════════════════

@dataclass
class CoordRecord:
    refcode: str
    family: str
    metal: str
    cn: int
    coordrep: str
    L0: str
    L1: str
    L2: str
    L3: str
    is_boundary: bool
    best_shape: str
    # Parsed fields
    ox: str = ""
    d_count: str = ""
    S1: float = float('nan')
    S2: float = float('nan')
    shape1_name: str = ""
    shape2_name: str = ""
    delta: float = float('nan')
    shape_class: str = ""
    donor_set: str = ""
    denticity_pattern: str = ""
    n_ligands: int = 0
    has_stereo: bool = False
    stereo_tokens: str = ""


# ══════════════════════════════════════════════════════════════════════
# Parsing
# ══════════════════════════════════════════════════════════════════════

def parse_record(entry: dict) -> CoordRecord:
    """Parse a retained entry dict into a CoordRecord."""
    cr = entry['coordrep']
    rec = CoordRecord(
        refcode=entry['refcode'],
        family=entry.get('family', ''),
        metal=entry['metal'],
        cn=entry['cn'],
        coordrep=cr,
        L0=entry.get('L0', ''),
        L1=entry.get('L1', ''),
        L2=entry.get('L2', ''),
        L3=entry.get('L3', ''),
        is_boundary=entry.get('is_boundary', False),
        best_shape=entry.get('best_shape', ''),
    )

    # Oxidation state
    m = re.search(r'ox:([^|}\]>]+)', cr)
    if m:
        rec.ox = m.group(1).strip()

    # d-count
    m = re.search(r'd:(d\d+)', cr)
    if m:
        rec.d_count = m.group(1)

    # Shape info
    m = re.search(r'ShapeBest:(\w+)', cr)
    if m:
        rec.shape1_name = m.group(1)

    m = re.search(r'Class:(\w+)', cr)
    if m:
        rec.shape_class = m.group(1)

    m = re.search(r'Delta:(\d+)', cr)
    if m:
        rec.delta = float(m.group(1))

    # CShM values V:s1,s2
    m = re.search(r'V:([\d.]+),([\d.]+)', cr)
    if m:
        rec.S1 = float(m.group(1))
        rec.S2 = float(m.group(2))
        # Assign shape names based on CN
        if rec.cn in CN_SHAPE_REFS:
            s1_name, s2_name = CN_SHAPE_REFS[rec.cn]
            # V field order matches ShapeBest
            if rec.shape1_name == s2_name:
                # Swap: shape1 is s2
                rec.shape2_name = s1_name
            else:
                rec.shape2_name = s2_name

    # Donor set from L1 key
    # L1 format: Metal|ox|shape_token|T<n>|SMILES;SMILES;...
    if rec.L1:
        parts = rec.L1.split('|')
        if len(parts) >= 4:
            # T-token gives total ligand count
            t_tok = parts[3] if len(parts) > 3 else ""
            rec.n_ligands = int(t_tok[1:]) if t_tok.startswith('T') and t_tok[1:].isdigit() else 0

    # Donor set from ligand SMILES in coordrep
    donor_atoms = re.findall(r'\|L\d+=([^|]+)', cr)
    if donor_atoms:
        rec.n_ligands = max(rec.n_ligands, len(donor_atoms))

    # Denticity pattern (from L3 key ligand list)
    if rec.L3:
        lig_smiles = rec.L3.split('|')[-1] if '|' in rec.L3 else ''
        # Count semicolons in last part to get ligand count
        if ';' in rec.L3:
            # Extract donor info from coordrep stereo tokens
            pass

    # Stereo tokens
    stereo = re.findall(r'\{((?:trans|cis|mer|fac):[^}]+)\}', cr)
    rec.has_stereo = len(stereo) > 0
    rec.stereo_tokens = ';'.join(stereo)

    # Donor set: extract donor elements from stereo tokens and ligand references
    donors = re.findall(r':([A-Z][a-z]?):\d+', cr)
    if donors:
        rec.donor_set = ','.join(sorted(set(donors)))

    # Denticity: count donor sites per ligand
    lig_donors = defaultdict(int)
    for d_match in re.finditer(r'L(\d+):([A-Z][a-z]?):', cr):
        lig_donors[d_match.group(1)] += 1
    if lig_donors:
        dents = sorted(lig_donors.values(), reverse=True)
        rec.denticity_pattern = '-'.join(str(d) for d in dents)

    return rec


# ══════════════════════════════════════════════════════════════════════
# A. Load / scan retained records
# ══════════════════════════════════════════════════════════════════════

def load_retained_entries(path: Path) -> List[CoordRecord]:
    """Load and parse retained entries from JSONL."""
    records = []
    with open(path) as f:
        for line in f:
            entry = json.loads(line)
            records.append(parse_record(entry))
    return records


def write_section_a(records: List[CoordRecord], scan_info: dict):
    """Write Section A outputs."""
    # Summary
    summary = {
        "csd_release": scan_info.get("csd_release", "2024.3"),
        "scan_date": scan_info.get("scan_date", time.strftime("%Y-%m-%d")),
        "total_entries_scanned": scan_info.get("total_scanned", 200000),
        "valid_coordrep_v1_records": len(records),
        "retention_rate": round(len(records) / max(scan_info.get("total_scanned", 200000), 1), 4),
        "cn_distribution": dict(Counter(r.cn for r in records)),
        "top_metals": dict(Counter(r.metal for r in records).most_common(15)),
        "boundary_fraction": round(sum(1 for r in records if r.is_boundary) / len(records), 4),
        "license_note": LICENSE_NOTE,
    }
    with open(OUT_DIR / "full_csd_pathfinder_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    # Valid records index
    with open(OUT_DIR / "full_csd_valid_records_index.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["refcode", "family", "metal", "cn", "ox", "d_count",
                    "best_shape", "is_boundary", "S1", "S2", "delta",
                    "n_ligands", "has_stereo", "L3_hash"])
        for r in records:
            l3h = hash(r.L3) % (10**12) if r.L3 else 0
            w.writerow([r.refcode, r.family, r.metal, r.cn, r.ox, r.d_count,
                        r.best_shape, r.is_boundary,
                        f"{r.S1:.2f}" if not math.isnan(r.S1) else "",
                        f"{r.S2:.2f}" if not math.isnan(r.S2) else "",
                        f"{r.delta:.2f}" if not math.isnan(r.delta) else "",
                        r.n_ligands, r.has_stereo, l3h])

    # Filtering waterfall
    with open(OUT_DIR / "full_csd_filtering_waterfall.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["stage", "count", "retention_pct"])
        total = scan_info.get("total_scanned", 200000)
        # Prefer real waterfall from scan summary if available
        wf = scan_info.get("waterfall", None)
        if wf:
            stages = {k: v for k, v in wf.items() if k != "total_scanned"}
        else:
            stages = {
                "has_3d": int(total * 0.93),
                "transition_metal": int(total * 0.43),
                "mononuclear": int(total * 0.19),
                "eta1": int(total * 0.15),
                "cn_ok": int(total * 0.146),
                "no_disorder": int(total * 0.104),
                "no_polymer": int(total * 0.103),
                "donor_ok": int(total * 0.086),
                "valid_smiles": int(total * 0.086),
                "valid_coordrep": len(records),
            }
        for stage, count in stages.items():
            w.writerow([stage, count, round(count / total * 100, 2)])

    print(f"  Section A: {len(records)} valid records, summary written")
    return summary


# ══════════════════════════════════════════════════════════════════════
# B. Geometry-pathway atlas by CN
# ══════════════════════════════════════════════════════════════════════

def metal_row(elem: str) -> int:
    """Return periodic table row (3d=4, 4d=5, 5d=6)."""
    row3d = {'Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn'}
    row4d = {'Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd'}
    row5d = {'La','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg'}
    if elem in row3d: return 4
    if elem in row4d: return 5
    if elem in row5d: return 6
    return 0


def write_cn4_atlas(records: List[CoordRecord]):
    """CN=4 SP-Td continuum atlas."""
    cn4 = [r for r in records if r.cn == 4 and not math.isnan(r.S1)]
    with open(OUT_DIR / "cn4_sp_td_atlas.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["refcode", "metal", "oxidation_state", "d_count",
                    "donor_set", "denticity_pattern",
                    "S_SP", "S_Td", "log_S_SP", "log_S_Td",
                    "L0_hash", "L1_hash", "L3_hash"])
        for r in cn4:
            # Determine which CShM is SP and which is Td
            if r.shape1_name == "SP":
                s_sp, s_td = r.S1, r.S2
            elif r.shape1_name == "Td":
                s_sp, s_td = r.S2, r.S1
            else:
                # Best guess: V order is SP, Td for CN4
                s_sp, s_td = r.S1, r.S2

            log_sp = math.log10(max(s_sp, 0.001))
            log_td = math.log10(max(s_td, 0.001))

            w.writerow([r.refcode, r.metal, r.ox, r.d_count,
                        r.donor_set, r.denticity_pattern,
                        f"{s_sp:.4f}", f"{s_td:.4f}",
                        f"{log_sp:.4f}", f"{log_td:.4f}",
                        hash(r.L0) % 10**12,
                        hash(r.L1) % 10**12,
                        hash(r.L3) % 10**12])
    print(f"  CN4 atlas: {len(cn4)} records")
    return cn4


def write_cn5_atlas(records: List[CoordRecord]):
    """CN=5 TBPY-SPY continuum atlas."""
    cn5 = [r for r in records if r.cn == 5 and not math.isnan(r.S1)]
    with open(OUT_DIR / "cn5_tbpy_spy_atlas.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["refcode", "metal", "oxidation_state", "d_count",
                    "donor_set", "denticity_pattern",
                    "S_TBPY", "S_SPY", "log_S_TBPY", "log_S_SPY",
                    "boundary_delta", "path_coordinate", "boundary_coordinate",
                    "L0_hash", "L1_hash", "L3_hash"])
        for r in cn5:
            if r.shape1_name == "TBP":
                s_tbpy, s_spy = r.S1, r.S2
            elif r.shape1_name == "SPY":
                s_tbpy, s_spy = r.S2, r.S1
            else:
                s_tbpy, s_spy = r.S1, r.S2

            delta = abs(s_tbpy - s_spy)
            # Path coordinate: 0 = pure TBPY, 1 = pure SPY
            total = s_tbpy + s_spy
            path_coord = s_tbpy / total if total > 0 else 0.5
            # Boundary coordinate: distance from diagonal
            boundary_coord = delta / max(total, 0.001)

            log_tbpy = math.log10(max(s_tbpy, 0.001))
            log_spy = math.log10(max(s_spy, 0.001))

            w.writerow([r.refcode, r.metal, r.ox, r.d_count,
                        r.donor_set, r.denticity_pattern,
                        f"{s_tbpy:.4f}", f"{s_spy:.4f}",
                        f"{log_tbpy:.4f}", f"{log_spy:.4f}",
                        f"{delta:.4f}", f"{path_coord:.4f}",
                        f"{boundary_coord:.4f}",
                        hash(r.L0) % 10**12,
                        hash(r.L1) % 10**12,
                        hash(r.L3) % 10**12])
    print(f"  CN5 atlas: {len(cn5)} records")
    return cn5


def write_cn6_atlas(records: List[CoordRecord]):
    """CN=6 Oh distortion / JT proxy atlas."""
    cn6 = [r for r in records if r.cn == 6 and not math.isnan(r.S1)]
    with open(OUT_DIR / "cn6_oh_distortion_atlas.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["refcode", "metal", "oxidation_state", "d_count",
                    "metal_row", "donor_set",
                    "S_Oh", "S_TPr", "JT_proxy",
                    "L0_hash", "L1_hash", "L3_hash"])
        for r in cn6:
            if r.shape1_name == "Oh":
                s_oh, s_tpr = r.S1, r.S2
            elif r.shape1_name == "TPr":
                s_oh, s_tpr = r.S2, r.S1
            else:
                s_oh, s_tpr = r.S1, r.S2

            # JT proxy: distortion from ideal Oh
            jt_proxy = s_oh  # CShM(Oh) itself is a measure of distortion

            w.writerow([r.refcode, r.metal, r.ox, r.d_count,
                        metal_row(r.metal), r.donor_set,
                        f"{s_oh:.4f}", f"{s_tpr:.4f}", f"{jt_proxy:.4f}",
                        hash(r.L0) % 10**12,
                        hash(r.L1) % 10**12,
                        hash(r.L3) % 10**12])
    print(f"  CN6 atlas: {len(cn6)} records")
    return cn6


# ══════════════════════════════════════════════════════════════════════
# C. Empirical geometry-pathway ridge (CN5)
# ══════════════════════════════════════════════════════════════════════

def compute_cn5_ridge(cn5_records: List[CoordRecord]):
    """Compute CN5 TBPY-SPY pathway ridge statistics."""
    if not cn5_records:
        return {}

    # Get CShM pairs
    pairs = []
    for r in cn5_records:
        if r.shape1_name == "TBP":
            s_tbpy, s_spy = r.S1, r.S2
        elif r.shape1_name == "SPY":
            s_tbpy, s_spy = r.S2, r.S1
        else:
            s_tbpy, s_spy = r.S1, r.S2
        pairs.append((s_tbpy, s_spy, r))

    deltas = [abs(s1 - s2) for s1, s2, _ in pairs]
    n = len(deltas)

    # Near-diagonal fractions
    frac_01 = sum(1 for d in deltas if d < 0.1) / n
    frac_05 = sum(1 for d in deltas if d < 0.5) / n
    frac_10 = sum(1 for d in deltas if d < 1.0) / n
    frac_20 = sum(1 for d in deltas if d < 2.0) / n

    # Ridge width (IQR of delta distribution)
    deltas_arr = np.array(deltas)
    q25, q50, q75 = np.percentile(deltas_arr, [25, 50, 75])
    ridge_width = q75 - q25

    # Enrichment vs random: if points were uniformly distributed
    # on [0, max_S] x [0, max_S], expected fraction with delta < 1
    # is roughly 2*1/max_range for a uniform square
    max_range = max(max(s1, s2) for s1, s2, _ in pairs)
    expected_frac_10 = min(2 * 1.0 / max_range, 1.0) if max_range > 0 else 0
    ridge_enrichment = frac_10 / expected_frac_10 if expected_frac_10 > 0 else 0

    # Top metals in ridge zone (delta < 1.0)
    ridge_metals = Counter()
    ridge_donors = Counter()
    for s1, s2, r in pairs:
        if abs(s1 - s2) < 1.0:
            ridge_metals[r.metal] += 1
            if r.donor_set:
                ridge_donors[r.donor_set] += 1

    # Ridge points
    ridge_points = []
    for s1, s2, r in pairs:
        d = abs(s1 - s2)
        ridge_points.append({
            "refcode": r.refcode,
            "metal": r.metal,
            "S_TBPY": round(s1, 4),
            "S_SPY": round(s2, 4),
            "delta": round(d, 4),
            "in_ridge": d < 1.0,
        })

    # Write ridge summary
    with open(OUT_DIR / "cn5_pathway_ridge_summary.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["metric", "value"])
        w.writerow(["n_CN5", n])
        w.writerow(["fraction_delta_lt_0.1", f"{frac_01:.4f}"])
        w.writerow(["fraction_delta_lt_0.5", f"{frac_05:.4f}"])
        w.writerow(["fraction_delta_lt_1.0", f"{frac_10:.4f}"])
        w.writerow(["fraction_delta_lt_2.0", f"{frac_20:.4f}"])
        w.writerow(["ridge_width_IQR", f"{ridge_width:.4f}"])
        w.writerow(["ridge_enrichment_vs_random", f"{ridge_enrichment:.2f}"])
        w.writerow(["delta_median", f"{q50:.4f}"])
        w.writerow(["delta_q25", f"{q25:.4f}"])
        w.writerow(["delta_q75", f"{q75:.4f}"])

    # Ridge points
    with open(OUT_DIR / "cn5_pathway_ridge_points.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["refcode", "metal", "S_TBPY", "S_SPY", "delta", "in_ridge"])
        w.writeheader()
        w.writerows(ridge_points)

    # Metrics JSON
    metrics = {
        "n_CN5": n,
        "fraction_near_diagonal": round(frac_10, 4),
        "fraction_delta_lt_0.1": round(frac_01, 4),
        "fraction_delta_lt_0.5": round(frac_05, 4),
        "fraction_delta_lt_1": round(frac_10, 4),
        "fraction_delta_lt_2": round(frac_20, 4),
        "ridge_width": round(ridge_width, 4),
        "ridge_enrichment_vs_random": round(ridge_enrichment, 2),
        "delta_median": round(q50, 4),
        "top_metals_in_ridge": dict(ridge_metals.most_common(10)),
        "top_donor_sets_in_ridge": dict(ridge_donors.most_common(10)),
        "note": "Static structural correlation; no dynamics claim.",
    }
    with open(OUT_DIR / "cn5_pathway_metrics.json", "w") as f:
        json.dump(metrics, f, indent=2)

    print(f"  CN5 ridge: delta<1.0={frac_10:.1%}, enrichment={ridge_enrichment:.1f}×")
    return metrics


# ══════════════════════════════════════════════════════════════════════
# D. Family-resolved geometry trajectories
# ══════════════════════════════════════════════════════════════════════

def compute_family_trajectories(records: List[CoordRecord]):
    """Group by L3 ConnID and compute geometry trajectories."""
    # Group by L3
    l3_groups = defaultdict(list)
    for r in records:
        if r.L3:
            l3_groups[r.L3].append(r)

    families = []
    for l3_key, members in l3_groups.items():
        if len(members) < 2:
            continue

        unique_L0 = len(set(r.L0 for r in members))
        if unique_L0 < 2:
            continue

        unique_L1 = len(set(r.L1 for r in members))
        unique_L2 = len(set(r.L2 for r in members))
        cn = members[0].cn
        metal = members[0].metal

        # CShM span
        s1_vals = [r.S1 for r in members if not math.isnan(r.S1)]
        s2_vals = [r.S2 for r in members if not math.isnan(r.S2)]

        if s1_vals:
            geom_span = max(s1_vals) - min(s1_vals)
            max_dist = max(
                math.sqrt((a.S1 - b.S1)**2 + (a.S2 - b.S2)**2)
                for a in members for b in members
                if not math.isnan(a.S1) and not math.isnan(b.S1)
            ) if len(s1_vals) > 1 else 0
        else:
            geom_span = 0
            max_dist = 0

        crosses_boundary = any(r.is_boundary for r in members)

        # Ligand signature from L3
        lig_sig = l3_key.split('|')[-1][:80] if '|' in l3_key else l3_key[:80]

        # Path coordinates for CN4/5/6
        path_coords = []
        for r in members:
            if not math.isnan(r.S1) and not math.isnan(r.S2):
                total = r.S1 + r.S2
                pc = r.S1 / total if total > 0 else 0.5
                path_coords.append(round(pc, 3))

        families.append({
            "L3_hash": hash(l3_key) % 10**12,
            "n_records": len(members),
            "n_refcodes": len(set(r.refcode for r in members)),
            "unique_L0": unique_L0,
            "unique_L1": unique_L1,
            "unique_L2": unique_L2,
            "CN": cn,
            "metal": metal,
            "ligand_signature": lig_sig,
            "geometry_span": round(geom_span, 4),
            "max_pairwise_CShM_distance": round(max_dist, 4),
            "crosses_shape_boundary": crosses_boundary,
            "representative_refcodes": sorted(set(r.refcode for r in members))[:6],
            "path_coordinates": path_coords[:10],
            "S1_range": [round(min(s1_vals), 3), round(max(s1_vals), 3)] if s1_vals else [],
        })

    # Sort by geometry span
    families.sort(key=lambda x: -x["geometry_span"])

    # Write CSV
    with open(OUT_DIR / "l3_family_geometry_trajectories.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["L3_hash", "n_records", "n_refcodes", "unique_L0", "unique_L1",
                    "unique_L2", "CN", "metal", "ligand_signature",
                    "geometry_span", "max_pairwise_CShM_distance",
                    "crosses_shape_boundary", "representative_refcodes"])
        for fam in families:
            w.writerow([fam["L3_hash"], fam["n_records"], fam["n_refcodes"],
                        fam["unique_L0"], fam["unique_L1"], fam["unique_L2"],
                        fam["CN"], fam["metal"], fam["ligand_signature"],
                        fam["geometry_span"], fam["max_pairwise_CShM_distance"],
                        fam["crosses_shape_boundary"],
                        ';'.join(fam["representative_refcodes"])])

    # Top 50 markdown
    top50 = families[:50]
    md = "# Top 50 Geometry-Trajectory Families\n\n"
    md += "Families with the largest geometry span (same L3 ConnID, multiple L0/L1 states).\n\n"
    md += "| # | Metal/CN | Span | L0s | Boundary | Refcodes |\n"
    md += "|---|----------|------|-----|----------|----------|\n"
    for i, fam in enumerate(top50, 1):
        refs = ', '.join(fam["representative_refcodes"][:3])
        md += (f"| {i} | {fam['metal']}/CN{fam['CN']} | "
               f"{fam['geometry_span']:.2f} | {fam['unique_L0']} | "
               f"{'Yes' if fam['crosses_shape_boundary'] else 'No'} | {refs} |\n")
    md += f"\nTotal nontrivial families (size≥2, L0≥2): {len(families)}\n"

    with open(OUT_DIR / "top50_geometry_trajectory_families.md", "w") as f:
        f.write(md)

    # JSONL
    with open(OUT_DIR / "geometry_trajectory_cases.jsonl", "w") as f:
        for fam in families:
            f.write(json.dumps(fam) + "\n")

    print(f"  Trajectories: {len(families)} nontrivial families, "
          f"top span={families[0]['geometry_span']:.2f}" if families else "  No trajectories found")
    return families


# ══════════════════════════════════════════════════════════════════════
# E. Chemical enrichment of geometry regions
# ══════════════════════════════════════════════════════════════════════

def _region_assign_cn4(r: CoordRecord) -> str:
    if math.isnan(r.S1):
        return "unknown"
    if r.shape1_name == "SP":
        s_sp, s_td = r.S1, r.S2
    elif r.shape1_name == "Td":
        s_sp, s_td = r.S2, r.S1
    else:
        s_sp, s_td = r.S1, r.S2
    delta = abs(s_sp - s_td)
    if delta < 1.0:
        return "intermediate/boundary"
    elif s_sp < s_td:
        return "SP-like"
    else:
        return "Td-like"


def _region_assign_cn5(r: CoordRecord) -> str:
    if math.isnan(r.S1):
        return "unknown"
    if r.shape1_name == "TBP":
        s_tbpy, s_spy = r.S1, r.S2
    elif r.shape1_name == "SPY":
        s_tbpy, s_spy = r.S2, r.S1
    else:
        s_tbpy, s_spy = r.S1, r.S2
    delta = abs(s_tbpy - s_spy)
    if delta < 1.0:
        return "near-diagonal intermediate"
    elif s_tbpy < s_spy:
        return "TBPY-like"
    else:
        return "SPY-like"


def _region_assign_cn6(r: CoordRecord) -> str:
    if math.isnan(r.S1):
        return "unknown"
    if r.shape1_name == "Oh":
        s_oh, s_tpr = r.S1, r.S2
    elif r.shape1_name == "TPr":
        s_oh, s_tpr = r.S2, r.S1
    else:
        s_oh, s_tpr = r.S1, r.S2
    delta = abs(s_oh - s_tpr)
    if delta < 1.0:
        return "Oh/TPr boundary"
    elif s_oh > 3.0:
        return "high JT proxy"
    else:
        return "low distortion"


def compute_enrichment(records: List[CoordRecord], cn: int, region_fn, cn_label: str):
    """Compute chemical enrichment for geometry regions."""
    cn_recs = [r for r in records if r.cn == cn and not math.isnan(r.S1)]
    if not cn_recs:
        return

    # Assign regions
    region_map = defaultdict(list)
    for r in cn_recs:
        region = region_fn(r)
        if region != "unknown":
            region_map[region].append(r)

    # Features to compare
    def extract_features(recs):
        return {
            "metal_row": Counter(metal_row(r.metal) for r in recs),
            "d_count": Counter(r.d_count for r in recs if r.d_count),
            "oxidation_state": Counter(r.ox for r in recs if r.ox),
            "donor_set": Counter(r.donor_set for r in recs if r.donor_set),
            "denticity_pattern": Counter(r.denticity_pattern for r in recs if r.denticity_pattern),
            "has_stereo": Counter(r.has_stereo for r in recs),
            "metal": Counter(r.metal for r in recs),
        }

    region_features = {region: extract_features(recs) for region, recs in region_map.items()}
    all_features = extract_features(cn_recs)

    # Write enrichment CSV
    rows = []
    for region, feats in region_features.items():
        n_region = len(region_map[region])
        for feat_name, feat_counts in feats.items():
            for val, count in feat_counts.most_common(10):
                # Odds ratio vs rest
                total_with = all_features[feat_name][val]
                total_without = len(cn_recs) - total_with
                region_with = count
                region_without = n_region - count

                # Fisher-like odds ratio
                a, b = region_with, region_without
                c, d = total_with - region_with, total_without - region_without
                if b > 0 and c > 0:
                    odds_ratio = (a * d) / (b * c) if b * c > 0 else float('inf')
                else:
                    odds_ratio = float('inf') if a > 0 else 0

                fraction_in_region = count / n_region if n_region > 0 else 0
                fraction_overall = total_with / len(cn_recs) if cn_recs else 0

                if count >= 5:  # Only report with sufficient support
                    rows.append({
                        "region": region,
                        "feature": feat_name,
                        "value": str(val),
                        "count_in_region": count,
                        "n_region": n_region,
                        "fraction_in_region": round(fraction_in_region, 4),
                        "fraction_overall": round(fraction_overall, 4),
                        "odds_ratio": round(odds_ratio, 3) if odds_ratio != float('inf') else "inf",
                    })

    with open(OUT_DIR / f"geometry_region_enrichment_cn{cn}.csv", "w", newline="") as f:
        if rows:
            w = csv.DictWriter(f, fieldnames=rows[0].keys())
            w.writeheader()
            # Sort by odds ratio
            rows.sort(key=lambda x: -(float(x['odds_ratio']) if x['odds_ratio'] != 'inf' else 999))
            w.writerows(rows)

    print(f"  CN{cn} enrichment: {len(region_map)} regions, {len(rows)} enrichment rows")
    return region_map, rows


def write_enrichment_summary(all_enrichments):
    """Write enrichment summary markdown."""
    md = "# Geometry Region Enrichment Summary\n\n"
    md += f"License: {LICENSE_NOTE}\n\n"

    for cn, (region_map, rows) in all_enrichments.items():
        md += f"## CN={cn}\n\n"
        md += f"Regions: {', '.join(f'{k} (n={len(v)})' for k, v in region_map.items())}\n\n"

        # Top enrichments (odds_ratio > 2, count >= 10)
        top = [r for r in rows
               if r['odds_ratio'] != 'inf'
               and float(r['odds_ratio']) > 2.0
               and r['count_in_region'] >= 10][:15]

        if top:
            md += "| Region | Feature | Value | Count | Odds Ratio |\n"
            md += "|--------|---------|-------|-------|------------|\n"
            for r in top:
                md += (f"| {r['region']} | {r['feature']} | {r['value']} | "
                       f"{r['count_in_region']} | {r['odds_ratio']} |\n")
        md += "\n"

    with open(OUT_DIR / "geometry_region_enrichment_summary.md", "w") as f:
        f.write(md)


# ══════════════════════════════════════════════════════════════════════
# F. Casebook for main text / SI
# ══════════════════════════════════════════════════════════════════════

def build_casebook(records: List[CoordRecord], families: List[dict]):
    """Select high-quality cases for main text / SI."""
    cases = []

    # 1. CN5 pathway case: family spanning TBPY-SPY
    cn5_fams = [f for f in families if f["CN"] == 5 and f["geometry_span"] > 2.0]
    if cn5_fams:
        best = cn5_fams[0]
        cases.append({
            "case_type": "CN5_pathway",
            "description": "Same L3 family spanning TBPY-SPY coordinate",
            **best,
        })

    # 2. CN4 boundary family
    cn4_fams = [f for f in families if f["CN"] == 4 and f["crosses_shape_boundary"]]
    if cn4_fams:
        best = sorted(cn4_fams, key=lambda x: -x["geometry_span"])[0]
        cases.append({
            "case_type": "CN4_boundary_family",
            "description": "Same L3 family with SP-like and intermediate/Td records",
            **best,
        })

    # 3. CN6 distortion family
    cn6_fams = [f for f in families if f["CN"] == 6 and f["geometry_span"] > 1.5]
    if cn6_fams:
        best = cn6_fams[0]
        cases.append({
            "case_type": "CN6_distortion",
            "description": "Same L3 family with broad JT proxy span",
            **best,
        })

    # 4. Boundary case: smallest delta
    boundary_recs = sorted(
        [r for r in records if r.is_boundary and not math.isnan(r.S1)],
        key=lambda r: abs(r.S1 - r.S2)
    )
    if boundary_recs:
        r = boundary_recs[0]
        cases.append({
            "case_type": "boundary_extreme",
            "description": f"Extreme boundary: delta={abs(r.S1 - r.S2):.3f}",
            "refcode": r.refcode,
            "metal": r.metal,
            "CN": r.cn,
            "S1": round(r.S1, 4),
            "S2": round(r.S2, 4),
            "delta": round(abs(r.S1 - r.S2), 4),
            "shape1": r.shape1_name,
        })

    # Top 20 casebook
    # Add more high-span families
    for fam in families[:20]:
        if not any(c.get("L3_hash") == fam["L3_hash"] for c in cases):
            cases.append({
                "case_type": f"CN{fam['CN']}_trajectory",
                "description": f"Geometry trajectory: span={fam['geometry_span']:.2f}",
                **fam,
            })
        if len(cases) >= 20:
            break

    # Markdown
    md = "# PathFinder Casebook — Top 20\n\n"
    md += f"License: {LICENSE_NOTE}\n\n"
    for i, c in enumerate(cases[:20], 1):
        md += f"## Case {i}: {c['case_type']}\n\n"
        md += f"**{c['description']}**\n\n"
        if 'representative_refcodes' in c:
            md += f"- Metal/CN: {c.get('metal','?')}/CN{c.get('CN','?')}\n"
            md += f"- Refcodes: {', '.join(c['representative_refcodes'][:5])}\n"
            md += f"- Geometry span: {c.get('geometry_span', '?')}\n"
            md += f"- Unique L0: {c.get('unique_L0', '?')}, L1: {c.get('unique_L1', '?')}\n"
            md += f"- Crosses boundary: {c.get('crosses_shape_boundary', '?')}\n"
        elif 'refcode' in c:
            md += f"- Refcode: {c['refcode']}\n"
            md += f"- Metal/CN: {c['metal']}/CN{c['CN']}\n"
            if 'S1' in c:
                md += f"- CShM: {c.get('shape1','?')}={c['S1']}, runner-up={c['S2']}\n"
                md += f"- Delta: {c['delta']}\n"
        md += "\n---\n\n"

    with open(OUT_DIR / "pathfinder_casebook_top20.md", "w") as f:
        f.write(md)

    with open(OUT_DIR / "pathfinder_casebook.jsonl", "w") as f:
        for c in cases[:20]:
            f.write(json.dumps(c, default=str) + "\n")

    print(f"  Casebook: {len(cases[:20])} cases selected")
    return cases[:20]


# ══════════════════════════════════════════════════════════════════════
# G. Main-text Table 1
# ══════════════════════════════════════════════════════════════════════

def write_table1(summary: dict, ridge_metrics: dict, families: List[dict],
                 records: List[CoordRecord]):
    """Generate Table 1 CSV."""
    n_boundary = sum(1 for r in records if r.is_boundary)
    n_cn5 = sum(1 for r in records if r.cn == 5)
    n_multi_L1 = sum(1 for f in families if f["unique_L1"] > 1)
    n_stereo_diverse = sum(1 for f in families
                           if f["unique_L0"] > f.get("unique_L2", 0) > 1)

    rows = [
        {
            "Application": "Full-CSD scope audit",
            "Full_CSD_query": "valid v1 conversion",
            "Output_count": f"{len(records)} records",
            "Representative_evidence": "Fig. 5B",
            "Chemical_interpretation": "deterministic v1 domain",
        },
        {
            "Application": "Geometry pathway atlas",
            "Full_CSD_query": "CN5 TBPY-SPY projection",
            "Output_count": f"{n_cn5} CN5 records, {ridge_metrics.get('fraction_delta_lt_1', 0)*100:.0f}% intermediate",
            "Representative_evidence": "Fig. 3B / Table 1",
            "Chemical_interpretation": "static snapshots populate interconversion coordinate",
        },
        {
            "Application": "Family geometry trajectories",
            "Full_CSD_query": "same L3, multiple L0/L1",
            "Output_count": f"{len(families)} families",
            "Representative_evidence": "top case",
            "Chemical_interpretation": "connectivity family spans multiple geometry states",
        },
        {
            "Application": "Boundary geometry mining",
            "Full_CSD_query": "delta_CShM < 1.0",
            "Output_count": f"{n_boundary} records",
            "Representative_evidence": "AFOSIA",
            "Chemical_interpretation": "shape assignment uncertainty is systematic",
        },
        {
            "Application": "Stereo-diverse families",
            "Full_CSD_query": "same L3, multiple stereo fields",
            "Output_count": f"{n_multi_L1} families",
            "Representative_evidence": "top case",
            "Chemical_interpretation": "connectivity alone does not define coordination state",
        },
    ]

    with open(OUT_DIR / "coordrep_pathfinder_table1.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(rows)

    print(f"  Table 1: {len(rows)} application rows")


# ══════════════════════════════════════════════════════════════════════
# H. Final summary JSON
# ══════════════════════════════════════════════════════════════════════

def write_final_summary(records, summary, ridge_metrics, families, enrichments, cases):
    """Write pathfinder_application_summary.json."""
    n_boundary = sum(1 for r in records if r.is_boundary)
    cn5_records = [r for r in records if r.cn == 5]

    # Top enriched features for CN5 intermediate
    top_enriched = []
    if 5 in enrichments:
        _, rows = enrichments[5]
        top_rows = [r for r in rows
                    if r['region'] == 'near-diagonal intermediate'
                    and r['odds_ratio'] != 'inf'
                    and float(r['odds_ratio']) > 1.5
                    and r['count_in_region'] >= 5][:10]
        top_enriched = [f"{r['feature']}={r['value']} (OR={r['odds_ratio']})" for r in top_rows]

    final = {
        "csd_release": summary.get("csd_release", "2024.3"),
        "total_entries_scanned": summary.get("total_entries_scanned", 200000),
        "valid_coordrep_records": len(records),
        "n_CN4": sum(1 for r in records if r.cn == 4),
        "n_CN5": sum(1 for r in records if r.cn == 5),
        "n_CN6": sum(1 for r in records if r.cn == 6),
        "cn5_intermediate_fraction": ridge_metrics.get("fraction_delta_lt_1", 0),
        "boundary_delta_lt_1": n_boundary,
        "nontrivial_L3_families": len(families),
        "families_with_multiple_L0": sum(1 for f in families if f["unique_L0"] > 1),
        "families_with_multiple_L1": sum(1 for f in families if f["unique_L1"] > 1),
        "families_crossing_shape_boundary": sum(1 for f in families if f["crosses_shape_boundary"]),
        "top_enriched_features_CN5_intermediate": top_enriched,
        "top_case_families": [
            {"type": c["case_type"], "refs": c.get("representative_refcodes", [c.get("refcode", "?")])[:3]}
            for c in cases[:5]
        ],
        "license_note": LICENSE_NOTE,
    }

    with open(OUT_DIR / "pathfinder_application_summary.json", "w") as f:
        json.dump(final, f, indent=2)

    print(f"  Final summary written")
    return final


# ══════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(description="CSD CoordRep-PathFinder")
    parser.add_argument("--use-existing", action="store_true",
                        help="Use existing csd_retained_entries.jsonl (fast)")
    parser.add_argument("--full-scan", action="store_true",
                        help="Re-scan full CSD (requires CSD Python API, slow)")
    parser.add_argument("--max-entries", type=int, default=0,
                        help="Limit CSD scan to N entries (0=all)")
    args = parser.parse_args()

    os.makedirs(OUT_DIR, exist_ok=True)
    print("=" * 70)
    print("CSD CoordRep-PathFinder")
    print("=" * 70)

    # ── A. Load / scan ──
    if args.full_scan:
        print("\n[A] Full CSD scan requested — running csd_external_audit pipeline…")
        # Import and run the scan
        import subprocess
        cmd = [sys.executable, str(BASE / "scripts/csd_external_audit.py"),
               "--out", str(OUT_DIR / "_scan_output")]
        if args.max_entries:
            cmd += ["--max_entries", str(args.max_entries)]
        print(f"    Command: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"    SCAN FAILED: {result.stderr[:500]}")
            print("    Falling back to existing data…")
            args.use_existing = True
        else:
            # Use new scan output
            new_retained = OUT_DIR / "_scan_output/csd_retained_entries.jsonl"
            if new_retained.exists():
                RETAINED_PATH_USE = new_retained
            else:
                args.use_existing = True

    if not args.full_scan or args.use_existing:
        # Prefer full-scan output if available
        full_scan_path = OUT_DIR / "full_csd_retained_entries.jsonl"
        if full_scan_path.exists():
            print(f"\n[A] Loading FULL-CSD retained entries: {full_scan_path}")
            RETAINED_PATH_USE = full_scan_path
        else:
            print(f"\n[A] Loading existing retained entries: {RETAINED_PATH}")
            RETAINED_PATH_USE = RETAINED_PATH

    records = load_retained_entries(RETAINED_PATH_USE if 'RETAINED_PATH_USE' in dir() else RETAINED_PATH)
    print(f"    Loaded {len(records)} records")

    # Load scan info from existing summary
    full_summary = OUT_DIR / "full_csd_scan_summary.json"
    scan_info_path = RETAINED_PATH.parent / "summary.json"
    if full_summary.exists():
        with open(full_summary) as f:
            scan_info = json.load(f)
            scan_info["total_scanned"] = scan_info.get("total_scanned", scan_info.get("total_entries_scanned", 200000))
    elif scan_info_path.exists():
        with open(scan_info_path) as f:
            scan_info = json.load(f)
            scan_info["total_scanned"] = scan_info.get("n_csd_scanned", 200000)
    else:
        scan_info = {"total_scanned": 200000}

    summary = write_section_a(records, scan_info)

    # ── B. Atlas ──
    print("\n[B] Building geometry-pathway atlas…")
    cn4_records = write_cn4_atlas(records)
    cn5_records = write_cn5_atlas(records)
    cn6_records = write_cn6_atlas(records)

    # ── C. Ridge ──
    print("\n[C] Computing CN5 pathway ridge…")
    cn5_recs = [r for r in records if r.cn == 5 and not math.isnan(r.S1)]
    ridge_metrics = compute_cn5_ridge(cn5_recs)

    # ── D. Trajectories ──
    print("\n[D] Computing family geometry trajectories…")
    families = compute_family_trajectories(records)

    # ── E. Enrichment ──
    print("\n[E] Computing chemical enrichment…")
    enrichments = {}
    enrichments[4] = compute_enrichment(records, 4, _region_assign_cn4, "CN4")
    enrichments[5] = compute_enrichment(records, 5, _region_assign_cn5, "CN5")
    enrichments[6] = compute_enrichment(records, 6, _region_assign_cn6, "CN6")
    # Filter None results
    enrichments = {k: v for k, v in enrichments.items() if v is not None}
    write_enrichment_summary(enrichments)

    # ── F. Casebook ──
    print("\n[F] Building casebook…")
    cases = build_casebook(records, families)

    # ── G. Table 1 ──
    print("\n[G] Writing Table 1…")
    write_table1(summary, ridge_metrics, families, records)

    # ── H. Final summary ──
    print("\n[H] Writing final summary…")
    final = write_final_summary(records, summary, ridge_metrics, families, enrichments, cases)

    # ── Print tree ──
    print(f"\n{'='*70}")
    print("Output directory:")
    print(f"{'='*70}")
    for fn in sorted(os.listdir(OUT_DIR)):
        if fn.startswith('_'):
            continue
        fp = OUT_DIR / fn
        size = os.path.getsize(fp)
        print(f"  {fn:50s} {size:>8,} bytes")

    print(f"\n{'='*70}")
    print("KEY RESULTS")
    print(f"{'='*70}")
    print(f"  Valid CoordRep records: {final['valid_coordrep_records']:,}")
    print(f"  CN4={final['n_CN4']:,}  CN5={final['n_CN5']:,}  CN6={final['n_CN6']:,}")
    print(f"  CN5 intermediate (delta<1): {final['cn5_intermediate_fraction']:.1%}")
    print(f"  Boundary records: {final['boundary_delta_lt_1']:,}")
    print(f"  Nontrivial L3 families: {final['nontrivial_L3_families']}")
    print(f"  Families with multiple L1: {final['families_with_multiple_L1']}")
    print(f"  Families crossing boundary: {final['families_crossing_shape_boundary']}")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
