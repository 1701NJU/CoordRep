#!/usr/bin/env python3
"""
screen_visual_families.py
==========================
Screen CN4 L3 families for visual interpretability by computing
real donor-metal-donor angle sets and geometric metrics via CSD API.

Run:  /data/miniconda3/envs/1701/bin/python scripts/screen_visual_families.py
Out:  revision_results/fig5d_family_trajectory/
"""
from __future__ import annotations
import csv, json, math, os, sys, itertools
from collections import defaultdict
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent))
from coordrep.io.tmqm_reader import TRANSITION_METALS
from coordrep.graph.neighbors import DONOR_ELEMENTS
from ccdc.io import EntryReader

BASE = Path(__file__).parent.parent
SRC  = BASE / "revision_results" / "csd_pathfinder_full"
OUT  = BASE / "revision_results" / "fig5d_family_trajectory"
os.makedirs(OUT, exist_ok=True)

LICENSE = "CSD-derived; internal figure rendering only."
MAX_CANDIDATES = 40  # screen top N families by span

# ══════════════════════════════════════════════════════════════
# Geometry computation helpers
# ══════════════════════════════════════════════════════════════

def get_metal_and_donors_csd(mol):
    """Return (metal_atom, [donor_atoms]) from a CSD molecule."""
    metals = [a for a in mol.atoms if a.atomic_symbol in TRANSITION_METALS]
    if len(metals) != 1:
        return None, []
    metal = metals[0]
    donors = []
    for b in metal.bonds:
        other = b.atoms[0] if b.atoms[1] == metal else b.atoms[1]
        if other.atomic_symbol != "H" and other.atomic_symbol in DONOR_ELEMENTS:
            donors.append(other)
    return metal, donors


def atom_coords(atom):
    c = atom.coordinates
    return np.array([c.x, c.y, c.z])


def angle_deg(v1, v2):
    """Angle between two vectors in degrees."""
    cos = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-12)
    cos = np.clip(cos, -1.0, 1.0)
    return math.degrees(math.acos(cos))


def compute_dmd_angles(metal, donors):
    """Compute all donor-metal-donor angles. Returns sorted list."""
    mc = atom_coords(metal)
    dc = [atom_coords(d) for d in donors]
    angles = []
    for i in range(len(dc)):
        for j in range(i+1, len(dc)):
            a = angle_deg(dc[i] - mc, dc[j] - mc)
            angles.append(round(a, 2))
    angles.sort()
    return angles


def compute_md_distances(metal, donors):
    """Compute metal-donor distances."""
    mc = atom_coords(metal)
    return sorted([round(np.linalg.norm(atom_coords(d) - mc), 4) for d in donors])


def angle_rmsd(a1, a2):
    """RMSD between two sorted angle lists."""
    if len(a1) != len(a2):
        return 999.0
    return math.sqrt(sum((x - y)**2 for x, y in zip(a1, a2)) / len(a1))


def dist_rmsd(d1, d2):
    if len(d1) != len(d2):
        return 999.0
    return math.sqrt(sum((x - y)**2 for x, y in zip(d1, d2)) / len(d1))


def trans_angle_count(angles, threshold=160.0):
    """Count angles >= threshold (likely trans)."""
    return sum(1 for a in angles if a >= threshold)


def cis_angle_count(angles, threshold=100.0):
    """Count angles <= threshold (likely cis)."""
    return sum(1 for a in angles if a <= threshold)


def classify_pattern(angles):
    """Classify as SP-like, Td-like, seesaw, etc."""
    n_trans = trans_angle_count(angles)
    n_cis = cis_angle_count(angles)
    if n_trans == 2:
        return "SP"
    elif n_trans == 0 and all(90 < a < 130 for a in angles):
        return "Td"
    elif n_trans == 1:
        return "seesaw"
    else:
        return "other"


# ══════════════════════════════════════════════════════════════
# XYZ helpers
# ══════════════════════════════════════════════════════════════

def write_xyz(atoms, comment, fp):
    with open(fp, "w") as f:
        f.write(f"{len(atoms)}\n{comment}\n")
        for e,x,y,z in atoms:
            f.write(f"{e:3s} {x:12.6f} {y:12.6f} {z:12.6f}\n")


def first_sphere_from_mol(mol, metal, donors):
    al = list(mol.atoms)
    keep = set()
    mi = al.index(metal); keep.add(mi)
    for d in donors:
        di = al.index(d); keep.add(di)
        for b in d.bonds:
            o = b.atoms[0] if b.atoms[1]==d else b.atoms[1]
            if o in al: keep.add(al.index(o))
    return [(al[i].atomic_symbol, al[i].coordinates.x, al[i].coordinates.y, al[i].coordinates.z)
            for i in sorted(keep) if al[i].atomic_symbol!="H" and al[i].coordinates]


def polyhedron_from_mol(mol, metal, donors):
    result = [(metal.atomic_symbol, metal.coordinates.x, metal.coordinates.y, metal.coordinates.z)]
    for d in donors:
        if d.coordinates:
            result.append((d.atomic_symbol, d.coordinates.x, d.coordinates.y, d.coordinates.z))
    return result


# ══════════════════════════════════════════════════════════════
# 1. Load candidate families
# ══════════════════════════════════════════════════════════════

print("=" * 60)
print("Visual Family Screening")
print("=" * 60)

fam_rows = []
with open(SRC / "l3_family_geometry_trajectories.csv") as f:
    for r in csv.DictReader(f):
        if (int(r["CN"]) == 4 and int(r["n_records"]) >= 3 and
            int(r["unique_L1"]) >= 2 and r["crosses_shape_boundary"] == "True"):
            fam_rows.append(r)
fam_rows.sort(key=lambda r: -float(r["geometry_span"]))
fam_rows = fam_rows[:MAX_CANDIDATES]

# Load CN4 atlas
atlas_by_l3 = defaultdict(list)
with open(SRC / "cn4_sp_td_atlas.csv") as f:
    for r in csv.DictReader(f):
        atlas_by_l3[r["L3_hash"]].append(r)

print(f"\nScreening {len(fam_rows)} candidate families...")

# ══════════════════════════════════════════════════════════════
# 2. Compute geometry metrics via CSD API
# ══════════════════════════════════════════════════════════════

csd = EntryReader("CSD")

results = []

for fi, fam in enumerate(fam_rows):
    h = fam["L3_hash"]
    refs_str = fam["representative_refcodes"]
    members = atlas_by_l3.get(h, [])
    if not members:
        continue

    label = refs_str.split(";")[0]
    refcodes = [m["refcode"] for m in members]

    # Compute angles/distances for each member
    member_data = []
    for rc in refcodes:
        try:
            entry = csd.entry(rc)
            mol = entry.molecule
            metal, donors = get_metal_and_donors_csd(mol)
            if metal is None or len(donors) != 4:
                continue
            angles = compute_dmd_angles(metal, donors)
            dists = compute_md_distances(metal, donors)
            pattern = classify_pattern(angles)
            s_sp = float([m for m in members if m["refcode"]==rc][0]["S_SP"])
            s_td = float([m for m in members if m["refcode"]==rc][0]["S_Td"])
            delta = abs(s_sp - s_td)
            member_data.append({
                "refcode": rc,
                "angles": angles,
                "dists": dists,
                "pattern": pattern,
                "S_SP": s_sp, "S_Td": s_td, "delta": delta,
                "scalar": s_sp - s_td,
                "is_boundary": delta < 1.0,
                "n_trans": trans_angle_count(angles),
                "max_angle": max(angles),
                "min_angle": min(angles),
            })
        except Exception as e:
            pass

    if len(member_data) < 3:
        continue

    member_data.sort(key=lambda m: m["scalar"])

    # Pairwise angle RMSD
    pw_angle_rmsds = []
    for i in range(len(member_data)):
        for j in range(i+1, len(member_data)):
            pw_angle_rmsds.append(angle_rmsd(member_data[i]["angles"], member_data[j]["angles"]))
    max_angle_rmsd = max(pw_angle_rmsds) if pw_angle_rmsds else 0

    # Pairwise distance RMSD
    pw_dist_rmsds = []
    for i in range(len(member_data)):
        for j in range(i+1, len(member_data)):
            pw_dist_rmsds.append(dist_rmsd(member_data[i]["dists"], member_data[j]["dists"]))
    max_dist_rmsd = max(pw_dist_rmsds) if pw_dist_rmsds else 0

    # Trans/cis pattern change
    patterns = [m["pattern"] for m in member_data]
    unique_patterns = len(set(patterns))
    max_trans_start = member_data[0]["n_trans"]
    max_trans_end = member_data[-1]["n_trans"]
    trans_change = abs(max_trans_end - max_trans_start)

    # Max angle range: difference between largest max_angle and smallest min_angle across members
    max_angle_range = max(m["max_angle"] for m in member_data) - min(m["min_angle"] for m in member_data)

    # Start vs end angle signature difference
    start_end_rmsd = angle_rmsd(member_data[0]["angles"], member_data[-1]["angles"])

    # Has boundary member?
    has_boundary = any(m["is_boundary"] for m in member_data)

    results.append({
        "L3_hash": h,
        "family_label": label,
        "metal": fam["metal"],
        "n_members": len(member_data),
        "n_unique_L1": int(fam["unique_L1"]),
        "trajectory_span": float(fam["geometry_span"]),
        "has_boundary": has_boundary,
        "max_pairwise_angle_RMSD": round(max_angle_rmsd, 2),
        "mean_pairwise_angle_RMSD": round(sum(pw_angle_rmsds)/len(pw_angle_rmsds), 2) if pw_angle_rmsds else 0,
        "max_pairwise_dist_RMSD": round(max_dist_rmsd, 4),
        "start_end_angle_RMSD": round(start_end_rmsd, 2),
        "unique_angle_patterns": unique_patterns,
        "trans_count_change": trans_change,
        "max_angle_range": round(max_angle_range, 2),
        "start_pattern": patterns[0],
        "end_pattern": patterns[-1],
        "start_angles": member_data[0]["angles"],
        "end_angles": member_data[-1]["angles"],
        "member_data": member_data,
    })

    status = "✓" if max_angle_rmsd > 10 else " "
    print(f"  {status} {fi+1:2d}/{len(fam_rows)} {label:12s} {fam['metal']:3s}  "
          f"n={len(member_data)}  angle_RMSD={max_angle_rmsd:5.1f}  "
          f"s-e_RMSD={start_end_rmsd:5.1f}  patterns={'/'.join(patterns)}  "
          f"trans_Δ={trans_change}  bd={'Y' if has_boundary else 'N'}")

# ══════════════════════════════════════════════════════════════
# 3. Rank by visual interpretability
# ══════════════════════════════════════════════════════════════

# Score: prioritize angle RMSD, pattern diversity, trans change, boundary
for r in results:
    r["visual_score"] = (
        r["max_pairwise_angle_RMSD"] * 2.0 +
        r["start_end_angle_RMSD"] * 1.5 +
        r["unique_angle_patterns"] * 10.0 +
        r["trans_count_change"] * 15.0 +
        (5.0 if r["has_boundary"] else 0) +
        r["n_unique_L1"] * 3.0
    )

results.sort(key=lambda r: -r["visual_score"])

print(f"\n{'='*60}")
print("Top 15 by visual score")
print(f"{'='*60}")
for i, r in enumerate(results[:15]):
    print(f"  {i+1:2d}. score={r['visual_score']:5.1f}  {r['family_label']:12s} {r['metal']:3s}  "
          f"n={r['n_members']}  L1={r['n_unique_L1']}  "
          f"angle_RMSD={r['max_pairwise_angle_RMSD']:5.1f}  "
          f"s-e={r['start_end_angle_RMSD']:5.1f}  "
          f"patterns={r['start_pattern']}→{r['end_pattern']}  "
          f"trans_Δ={r['trans_count_change']}  bd={'Y' if r['has_boundary'] else 'N'}")
    print(f"      start_angles={r['start_angles']}")
    print(f"      end_angles  ={r['end_angles']}")

# ══════════════════════════════════════════════════════════════
# 4. Export top_visual_family_candidates.csv
# ══════════════════════════════════════════════════════════════

csv_fields = [
    "visual_rank", "family_label", "L3_key_hash", "metal", "CN",
    "n_members", "n_unique_L1", "trajectory_span",
    "has_boundary_member",
    "max_pairwise_angle_RMSD", "mean_pairwise_angle_RMSD",
    "start_end_angle_RMSD", "max_pairwise_dist_RMSD",
    "unique_angle_patterns", "trans_count_change", "max_angle_range",
    "start_pattern", "end_pattern",
    "start_angles", "end_angles",
    "visual_score",
]

with open(OUT / "top_visual_family_candidates.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=csv_fields)
    w.writeheader()
    for i, r in enumerate(results[:20]):
        w.writerow({
            "visual_rank": i+1,
            "family_label": r["family_label"],
            "L3_key_hash": r["L3_hash"],
            "metal": r["metal"],
            "CN": 4,
            "n_members": r["n_members"],
            "n_unique_L1": r["n_unique_L1"],
            "trajectory_span": r["trajectory_span"],
            "has_boundary_member": r["has_boundary"],
            "max_pairwise_angle_RMSD": r["max_pairwise_angle_RMSD"],
            "mean_pairwise_angle_RMSD": r["mean_pairwise_angle_RMSD"],
            "start_end_angle_RMSD": r["start_end_angle_RMSD"],
            "max_pairwise_dist_RMSD": r["max_pairwise_dist_RMSD"],
            "unique_angle_patterns": r["unique_angle_patterns"],
            "trans_count_change": r["trans_count_change"],
            "max_angle_range": r["max_angle_range"],
            "start_pattern": r["start_pattern"],
            "end_pattern": r["end_pattern"],
            "start_angles": r["start_angles"],
            "end_angles": r["end_angles"],
            "visual_score": round(r["visual_score"], 1),
        })

print(f"\n  → top_visual_family_candidates.csv (top 20)")

# ══════════════════════════════════════════════════════════════
# 5. Export XYZ for top 5 candidates
# ══════════════════════════════════════════════════════════════

print(f"\n{'='*60}")
print("Exporting XYZ for top 5 visual families")
print(f"{'='*60}")

for rank, r in enumerate(results[:5], 1):
    label = r["family_label"]
    md = r["member_data"]
    # Pick start, intermediate (boundary-like or middle), end
    start = md[0]
    end = md[-1]
    # Find best intermediate: prefer boundary, else middle
    boundary_members = [m for m in md[1:-1] if m["is_boundary"]]
    if boundary_members:
        mid = min(boundary_members, key=lambda m: m["delta"])
    else:
        mid = md[len(md)//2]

    picks = [("start", start), ("mid", mid), ("end", end)]

    print(f"\n  [{rank}] {label} ({r['metal']}, score={r['visual_score']:.1f})")

    for tag, m in picks:
        rc = m["refcode"]
        prefix = f"vis{rank}_{label}_{tag}"
        try:
            entry = csd.entry(rc)
            mol = entry.molecule
            metal, donors = get_metal_and_donors_csd(mol)
            if metal is None:
                print(f"      {tag}: {rc} — no metal found, skipping")
                continue
            fs = first_sphere_from_mol(mol, metal, donors)
            poly = polyhedron_from_mol(mol, metal, donors)
            write_xyz(fs,
                      f"CSD={rc} | {label} {tag} | first_sphere | angles={m['angles']}",
                      OUT / f"{prefix}_first_sphere.xyz")
            write_xyz(poly,
                      f"CSD={rc} | {label} {tag} | polyhedron | pattern={m['pattern']}",
                      OUT / f"{prefix}_polyhedron.xyz")
            print(f"      {tag}: {rc:12s}  pattern={m['pattern']:6s}  "
                  f"angles={m['angles']}  fs={len(fs)} poly={len(poly)}")
        except Exception as e:
            print(f"      {tag}: {rc} ERROR: {e}")

# ══════════════════════════════════════════════════════════════
# 6. Summary JSON
# ══════════════════════════════════════════════════════════════

summary = {
    "description": "Visual family screening: top CN4 families ranked by donor-geometry diversity",
    "n_screened": len(fam_rows),
    "n_with_valid_geometry": len(results),
    "scoring_formula": "angle_RMSD*2 + start_end_RMSD*1.5 + unique_patterns*10 + trans_change*15 + boundary*5 + n_L1*3",
    "top_5": [
        {"rank": i+1, "label": r["family_label"], "metal": r["metal"],
         "visual_score": round(r["visual_score"], 1),
         "angle_RMSD": r["max_pairwise_angle_RMSD"],
         "start_pattern": r["start_pattern"], "end_pattern": r["end_pattern"],
         "trans_change": r["trans_count_change"],
         "n_members": r["n_members"], "has_boundary": r["has_boundary"]}
        for i, r in enumerate(results[:5])
    ],
    "license_note": LICENSE,
}
json.dump(summary, open(OUT / "visual_family_screening_summary.json", "w"), indent=2)

# Verification
print(f"\n{'='*60}")
print("Verification")
print(f"{'='*60}")
for fn in ["top_visual_family_candidates.csv", "visual_family_screening_summary.json"]:
    p = OUT / fn
    if p.exists():
        print(f"  ✓ {fn:50s} {os.path.getsize(p):>8,} B")
xyz_count = 0
for f in OUT.glob("vis*_*.xyz"):
    xyz_count += 1
print(f"  ✓ {xyz_count} XYZ files exported")
print("=" * 60)
