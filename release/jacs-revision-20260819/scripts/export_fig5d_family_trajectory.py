#!/usr/bin/env python3
"""
export_fig5d_family_trajectory.py
==================================
Export figure-ready data for Fig. 5D:
  - Top CN4 family trajectories overlay
  - Family centroids
  - Family-to-family edges
  - LINMOL XYZ thumbnails
  - Preview plot

Run:  /data/miniconda3/envs/1701/bin/python scripts/export_fig5d_family_trajectory.py
Out:  revision_results/fig5d_family_trajectory/
"""
from __future__ import annotations
import csv, json, math, os, random, sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))
# Import coordrep BEFORE ccdc to avoid numpy/segfault conflict
from coordrep_tools.csd_adapter import _count_metals, _get_donor_neighbors
from ccdc.io import EntryReader

BASE = Path(__file__).parent.parent
SRC  = BASE / "revision_results" / "csd_pathfinder_full"
OUT  = BASE / "revision_results" / "fig5d_family_trajectory"
os.makedirs(OUT, exist_ok=True)

LICENSE = "CSD-derived coordinates for internal figure rendering only; do not redistribute raw XYZ."
TOP_N = 15  # number of top families to export

LINMOL_L3 = "463888590564"

# ══════════════════════════════════════════════════════════════
# XYZ helpers
# ══════════════════════════════════════════════════════════════
def write_xyz(atoms, comment, fp):
    with open(fp, "w") as f:
        f.write(f"{len(atoms)}\n{comment}\n")
        for e,x,y,z in atoms:
            f.write(f"{e:3s} {x:12.6f} {y:12.6f} {z:12.6f}\n")

def mol_no_H(mol):
    return [(a.atomic_symbol,a.coordinates.x,a.coordinates.y,a.coordinates.z)
            for a in mol.atoms if a.coordinates and a.atomic_symbol != "H"]

def first_sphere(mol):
    metals, nm = _count_metals(mol)
    if nm == 0: return []
    m = metals[0]; donors = _get_donor_neighbors(mol, m)
    al = list(mol.atoms); keep = {al.index(m)}
    for d in donors:
        di = al.index(d); keep.add(di)
        for b in d.bonds:
            o = b.atoms[0] if b.atoms[1]==d else b.atoms[1]
            if o in al: keep.add(al.index(o))
    return [(al[i].atomic_symbol,al[i].coordinates.x,al[i].coordinates.y,al[i].coordinates.z)
            for i in sorted(keep) if al[i].atomic_symbol!="H" and al[i].coordinates]

def polyhedron(mol):
    metals, nm = _count_metals(mol)
    if nm == 0: return []
    m = metals[0]; donors = _get_donor_neighbors(mol, m)
    al = list(mol.atoms); idxs = [al.index(m)] + [al.index(d) for d in donors if d in al]
    return [(al[i].atomic_symbol,al[i].coordinates.x,al[i].coordinates.y,al[i].coordinates.z)
            for i in idxs if al[i].coordinates]

# ══════════════════════════════════════════════════════════════
# 1. Load data
# ══════════════════════════════════════════════════════════════
print("="*60)
print("Fig. 5D — Family-to-Family Trajectory Export")
print("="*60)

# --- Load family trajectories ---
fam_all = []
with open(SRC / "l3_family_geometry_trajectories.csv") as f:
    for r in csv.DictReader(f):
        if int(r["CN"])==4 and int(r["n_records"])>=2 and int(r["unique_L0"])>=2:
            fam_all.append(r)
fam_all.sort(key=lambda r: -float(r["geometry_span"]))

# Ensure LINMOL is in top set
top_fams = fam_all[:TOP_N]
linmol_in = any(r["L3_hash"]==LINMOL_L3 for r in top_fams)
if not linmol_in:
    for r in fam_all:
        if r["L3_hash"]==LINMOL_L3:
            top_fams.append(r); break
top_hashes = {r["L3_hash"] for r in top_fams}

print(f"\n[1] Loaded {len(fam_all)} CN4 families, exporting top {len(top_fams)}")

# --- Load CN4 atlas ---
atlas_by_l3 = defaultdict(list)
cn4_all = []
with open(SRC / "cn4_sp_td_atlas.csv") as f:
    for r in csv.DictReader(f):
        cn4_all.append(r)
        if r["L3_hash"] in top_hashes:
            atlas_by_l3[r["L3_hash"]].append(r)

print(f"    CN4 atlas: {len(cn4_all)} total records")

# --- Family labels (first refcode in representative list) ---
fam_label = {}
fam_meta  = {}
for r in top_fams:
    h = r["L3_hash"]
    refs = r["representative_refcodes"].split(";")
    fam_label[h] = refs[0] if refs else h
    fam_meta[h] = r

# ══════════════════════════════════════════════════════════════
# A. Top family trajectories overlay CSV
# ══════════════════════════════════════════════════════════════
print("\n[A] Writing overlay CSV...")

overlay_fields = [
    "family_rank","family_label","L3_key_hash","metal","CN",
    "n_records","n_refcodes","n_unique_L0","n_unique_L1","trajectory_span",
    "refcode","L0_key_hash","L1_shape_id","S_SP","S_Td",
    "log_S_SP","log_S_Td","trajectory_scalar",
    "is_LINMOL","is_selected_for_thumbnail",
]

# For LINMOL: pick 3 members for thumbnails (it has 3)
linmol_members = atlas_by_l3.get(LINMOL_L3, [])
for m in linmol_members:
    m["_s_sp"] = float(m["S_SP"]); m["_s_td"] = float(m["S_Td"])
    m["_scalar"] = m["_s_sp"] - m["_s_td"]
linmol_members.sort(key=lambda m: m["_scalar"])
linmol_thumb_refs = {m["refcode"] for m in linmol_members}  # all 3

overlay_rows = 0
with open(OUT / "fig5D_top_family_trajectories_overlay.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=overlay_fields)
    w.writeheader()
    for rank, fam in enumerate(top_fams, 1):
        h = fam["L3_hash"]
        members = atlas_by_l3.get(h, [])
        for m in members:
            m["_s_sp"] = float(m["S_SP"]); m["_s_td"] = float(m["S_Td"])
            m["_scalar"] = m["_s_sp"] - m["_s_td"]
        members.sort(key=lambda m: m["_scalar"])
        for m in members:
            w.writerow({
                "family_rank": rank,
                "family_label": fam_label[h],
                "L3_key_hash": h,
                "metal": fam["metal"],
                "CN": 4,
                "n_records": fam["n_records"],
                "n_refcodes": fam["n_refcodes"],
                "n_unique_L0": fam["unique_L0"],
                "n_unique_L1": fam["unique_L1"],
                "trajectory_span": fam["geometry_span"],
                "refcode": m["refcode"],
                "L0_key_hash": m["L0_hash"],
                "L1_shape_id": m["L1_hash"],
                "S_SP": f'{m["_s_sp"]:.4f}',
                "S_Td": f'{m["_s_td"]:.4f}',
                "log_S_SP": f'{math.log10(m["_s_sp"]+0.1):.4f}',
                "log_S_Td": f'{math.log10(m["_s_td"]+0.1):.4f}',
                "trajectory_scalar": f'{m["_scalar"]:.4f}',
                "is_LINMOL": h == LINMOL_L3,
                "is_selected_for_thumbnail": m["refcode"] in linmol_thumb_refs,
            })
            overlay_rows += 1

print(f"  → fig5D_top_family_trajectories_overlay.csv ({overlay_rows} rows, {len(top_fams)} families)")

# ══════════════════════════════════════════════════════════════
# B. Family centroids CSV
# ══════════════════════════════════════════════════════════════
print("\n[B] Writing centroids CSV...")

centroid_fields = [
    "family_rank","family_label","L3_key_hash","metal","CN",
    "n_records","n_refcodes","n_unique_L0","n_unique_L1","trajectory_span",
    "centroid_log_S_SP","centroid_log_S_Td","span_log_S_SP","span_log_S_Td",
    "crosses_shape_boundary","dominant_L1_shape","geometry_region",
]

with open(OUT / "fig5D_family_centroids.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=centroid_fields)
    w.writeheader()
    for rank, fam in enumerate(top_fams, 1):
        h = fam["L3_hash"]
        members = atlas_by_l3.get(h, [])
        if not members:
            continue
        log_sps = [math.log10(float(m["S_SP"])+0.1) for m in members]
        log_tds = [math.log10(float(m["S_Td"])+0.1) for m in members]
        centroid_sp = sum(log_sps)/len(log_sps)
        centroid_td = sum(log_tds)/len(log_tds)
        span_sp = max(log_sps) - min(log_sps)
        span_td = max(log_tds) - min(log_tds)
        # Dominant L1
        l1_counts = defaultdict(int)
        for m in members:
            l1_counts[m["L1_hash"]] += 1
        dominant_l1 = max(l1_counts, key=l1_counts.get)
        # Geometry region
        avg_delta = sum(abs(float(m["S_SP"])-float(m["S_Td"])) for m in members)/len(members)
        if avg_delta < 1.0:
            region = "boundary"
        elif centroid_sp < centroid_td:
            region = "SP-dominated"
        elif centroid_sp > centroid_td:
            region = "Td-dominated"
        else:
            region = "mixed"
        w.writerow({
            "family_rank": rank,
            "family_label": fam_label[h],
            "L3_key_hash": h,
            "metal": fam["metal"],
            "CN": 4,
            "n_records": fam["n_records"],
            "n_refcodes": fam["n_refcodes"],
            "n_unique_L0": fam["unique_L0"],
            "n_unique_L1": fam["unique_L1"],
            "trajectory_span": fam["geometry_span"],
            "centroid_log_S_SP": f"{centroid_sp:.4f}",
            "centroid_log_S_Td": f"{centroid_td:.4f}",
            "span_log_S_SP": f"{span_sp:.4f}",
            "span_log_S_Td": f"{span_td:.4f}",
            "crosses_shape_boundary": fam["crosses_shape_boundary"],
            "dominant_L1_shape": dominant_l1,
            "geometry_region": region,
        })

print(f"  → fig5D_family_centroids.csv ({len(top_fams)} rows)")

# ══════════════════════════════════════════════════════════════
# C. Family-to-family edges
# ══════════════════════════════════════════════════════════════
print("\n[C] Computing family-to-family edges...")

def donor_set_distance(ds1, ds2):
    """Symmetric difference size between two comma-separated donor sets."""
    s1 = set(ds1.split(",")) if ds1 else set()
    s2 = set(ds2.split(",")) if ds2 else set()
    return len(s1.symmetric_difference(s2))

def ligand_sig_similarity(sig1, sig2):
    """Jaccard on semicolon-split ligand tokens."""
    t1 = set(sig1.split(";")) if sig1 else set()
    t2 = set(sig2.split(";")) if sig2 else set()
    if not t1 and not t2:
        return 1.0
    if not t1 or not t2:
        return 0.0
    return len(t1 & t2) / len(t1 | t2)

# Gather donor_set per family from atlas members
fam_donor_set = {}
for h in top_hashes:
    members = atlas_by_l3.get(h, [])
    ds_counts = defaultdict(int)
    for m in members:
        ds_counts[m["donor_set"]] += 1
    fam_donor_set[h] = max(ds_counts, key=ds_counts.get) if ds_counts else ""

edge_fields = [
    "source_L3_key_hash","target_L3_key_hash",
    "source_family_label","target_family_label",
    "edge_reason","ligand_similarity","donor_set_distance",
    "same_metal","same_CN",
]

edges = []
top_list = list(top_fams)
for i in range(len(top_list)):
    for j in range(i+1, len(top_list)):
        fi, fj = top_list[i], top_list[j]
        hi, hj = fi["L3_hash"], fj["L3_hash"]
        same_metal = fi["metal"] == fj["metal"]
        same_cn = True  # all CN=4
        lig_sim = ligand_sig_similarity(fi["ligand_signature"], fj["ligand_signature"])
        ds_dist = donor_set_distance(fam_donor_set.get(hi,""), fam_donor_set.get(hj,""))

        reasons = []
        if lig_sim >= 0.75:
            reasons.append("ligand_similarity>=0.75")
        if ds_dist <= 1:
            reasons.append("donor_set_distance<=1")

        if reasons and same_metal:
            edges.append({
                "source_L3_key_hash": hi,
                "target_L3_key_hash": hj,
                "source_family_label": fam_label[hi],
                "target_family_label": fam_label[hj],
                "edge_reason": "; ".join(reasons),
                "ligand_similarity": f"{lig_sim:.3f}",
                "donor_set_distance": ds_dist,
                "same_metal": same_metal,
                "same_CN": same_cn,
            })

with open(OUT / "fig5D_family_edges.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=edge_fields)
    w.writeheader()
    w.writerows(edges)

print(f"  → fig5D_family_edges.csv ({len(edges)} edges)")

# ══════════════════════════════════════════════════════════════
# D. LINMOL XYZ thumbnails
# ══════════════════════════════════════════════════════════════
print("\n[D] Exporting LINMOL XYZ thumbnails...")

csd = EntryReader("CSD")
linmol_members.sort(key=lambda m: m["_scalar"])

for i, m in enumerate(linmol_members):
    rc = m["refcode"]
    sn = f"{i+1:02d}"
    try:
        entry = csd.entry(rc); mol = entry.molecule
        write_xyz(mol_no_H(mol), f"CSD_code={rc} | LINMOL state{sn} | no H",
                  OUT / f"LINMOL_state{sn}_full.xyz")
        write_xyz(first_sphere(mol), f"CSD_code={rc} | first_sphere | no H",
                  OUT / f"LINMOL_state{sn}_first_sphere.xyz")
        write_xyz(polyhedron(mol), f"CSD_code={rc} | polyhedron",
                  OUT / f"LINMOL_state{sn}_polyhedron.xyz")
        nf = len(mol_no_H(mol)); nfs = len(first_sphere(mol)); np_ = len(polyhedron(mol))
        print(f"    state{sn}: {rc:12s}  full={nf} fs={nfs} poly={np_}")
    except Exception as e:
        print(f"    state{sn}: {rc:12s}  ERROR: {e}")

# ══════════════════════════════════════════════════════════════
# E. Family statistics CSV (carried forward)
# ══════════════════════════════════════════════════════════════
with open(OUT / "fig5D_family_statistics.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["metric","value","label_for_figure"])
    w.writerow(["nontrivial_L3_families",6358,"6,358 nontrivial L3 families"])
    w.writerow(["multiple_L1_shape_families",3739,"3,739 families with multiple L1 shapes"])
    w.writerow(["shape_boundary_crossing_families",1810,"1,810 families crossing shape boundary"])

# ══════════════════════════════════════════════════════════════
# F. CN4 background sample (8 K)
# ══════════════════════════════════════════════════════════════
print("\n[E/F] Background sample + statistics...")
random.seed(42)
bg = random.sample(cn4_all, min(8000, len(cn4_all)))
with open(OUT / "fig5D_cn4_background_sample.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["refcode","metal","CN","S_SP","S_Td","log_S_SP","log_S_Td"])
    w.writeheader()
    for r in bg:
        sp = float(r["S_SP"]); td = float(r["S_Td"])
        w.writerow({"refcode":r["refcode"],"metal":r["metal"],"CN":4,
                     "S_SP":f"{sp:.4f}","S_Td":f"{td:.4f}",
                     "log_S_SP":f"{math.log10(sp+0.1):.4f}",
                     "log_S_Td":f"{math.log10(td+0.1):.4f}"})

# ══════════════════════════════════════════════════════════════
# G. Summary JSON
# ══════════════════════════════════════════════════════════════
summary = {
    "top_trajectory_families_exported": len(top_fams),
    "nontrivial_L3_families": 6358,
    "multiple_L1_shape_families": 3739,
    "shape_boundary_crossing_families": 1810,
    "top_family": "LINMOL",
    "top_family_metal": "Pd",
    "top_family_CN": 4,
    "top_family_span": 10.10,
    "families_exported": [
        {"rank": i+1, "label": fam_label[r["L3_hash"]], "metal": r["metal"],
         "span": float(r["geometry_span"]), "n_members": int(r["n_records"]),
         "n_L0": int(r["unique_L0"]), "n_L1": int(r["unique_L1"]),
         "boundary": r["crosses_shape_boundary"]=="True"}
        for i, r in enumerate(top_fams)
    ],
    "n_edges": len(edges),
    "license_note": LICENSE,
}
json.dump(summary, open(OUT / "fig5D_family_to_family_summary.json", "w"), indent=2)
print(f"\n  → fig5D_family_to_family_summary.json")

print("\n[H] Preview plot will be generated by separate matplotlib script.")

# ══════════════════════════════════════════════════════════════
# README
# ══════════════════════════════════════════════════════════════
readme = f"""# Fig. 5D — Family-to-Family Trajectory Visualization Data

## Overview

This directory contains figure-ready data for Fig. 5D, showing how
L3 connectivity families traverse the CN=4 SP–Td continuous shape measure space.

## Primary showcase family

- **LINMOL** (Pd, CN=4): 3 members, span = 10.10, boundary-crossing
- L3 ConnID hash: {LINMOL_L3}

## Files

### Trajectory data
- **`fig5D_top_family_trajectories_overlay.csv`** — {overlay_rows} point-rows across
  {len(top_fams)} top families; each row is one observed geometry state within a family.
- **`fig5D_family_centroids.csv`** — One row per family with centroid, span, region.
- **`fig5D_family_edges.csv`** — {len(edges)} family-to-family edges (same metal,
  ligand similarity ≥ 0.75 or donor-set distance ≤ 1).

### LINMOL XYZ thumbnails
- `LINMOL_state01_{{full,first_sphere,polyhedron}}.xyz` — trajectory start
- `LINMOL_state02_{{full,first_sphere,polyhedron}}.xyz` — intermediate
- `LINMOL_state03_{{full,first_sphere,polyhedron}}.xyz` — trajectory end
- **Recommended for rendering: `first_sphere.xyz`**

### Background & statistics
- **`fig5D_cn4_background_sample.csv`** — 8,000 CN4 atlas points for gray scatter.
- **`fig5D_family_statistics.csv`** — Key family counts for annotation.
- **`fig5D_family_to_family_summary.json`** — Master summary.

### Preview
- **`fig5D_trajectory_preview.png`** — Quick-look overlay of all exported families.

## Data provenance

All numbers from the **full-CSD run** (1,413,222 entries → 124,837 retained).
Trajectory scalar = S_SP − S_Td.

## License

{LICENSE}
No CIF, no reversible coordinate archives.
XYZ files for internal rendering only.
"""

with open(OUT / "README.md", "w") as f:
    f.write(readme)

# ══════════════════════════════════════════════════════════════
# Final verification
# ══════════════════════════════════════════════════════════════
req = [
    "fig5D_top_family_trajectories_overlay.csv",
    "fig5D_family_centroids.csv",
    "fig5D_family_edges.csv",
    "fig5D_family_to_family_summary.json",
    "fig5D_family_statistics.csv",
    "fig5D_cn4_background_sample.csv",
    "README.md",
]
for i in range(len(linmol_members)):
    sn = f"{i+1:02d}"
    req += [f"LINMOL_state{sn}_full.xyz",
            f"LINMOL_state{sn}_first_sphere.xyz",
            f"LINMOL_state{sn}_polyhedron.xyz"]

print(f"\n{'='*60}")
print("Output verification")
print(f"{'='*60}")
ok = 0
for fn in req:
    p = OUT / fn
    if p.exists():
        print(f"  ✓ {fn:55s} {os.path.getsize(p):>10,} B")
        ok += 1
    else:
        print(f"  ✗ {fn:55s} MISSING")

prev = OUT / "fig5D_trajectory_preview.png"
if prev.exists():
    print(f"  ✓ {'fig5D_trajectory_preview.png':55s} {os.path.getsize(prev):>10,} B  (preview)")

print(f"\nRequired: {ok}/{len(req)}")
print("="*60)
