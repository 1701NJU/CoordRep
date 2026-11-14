#!/usr/bin/env python3
"""
export_fig5d_gazkun_trajectory.py
==================================
Production export for Fig. 5D: GAZKUN family (Cu, CN4, N4 donors)
as the visually interpretable representative L3 family trajectory.

GAZKUN family: 8 members, 5 L1 shapes, boundary-crossing,
pattern: distorted → Td → SP (angles dramatically change).

Run:  /data/miniconda3/envs/1701/bin/python scripts/export_fig5d_gazkun_trajectory.py
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

LICENSE = "CSD-derived coordinates for internal figure rendering only; do not redistribute raw XYZ."

GAZKUN_L3 = "895774452902"
FAMILY_LABEL = "GAZKUN"
FAMILY_METAL = "Cu"
FAMILY_CN = 4

# ══════════════════════════════════════════════════════════════
# Geometry helpers
# ══════════════════════════════════════════════════════════════

def get_metal_and_donors(mol):
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

def atom_coords(a):
    return np.array([a.coordinates.x, a.coordinates.y, a.coordinates.z])

def angle_deg(v1, v2):
    cos = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-12)
    return math.degrees(math.acos(np.clip(cos, -1, 1)))

def compute_dmd_angles(metal, donors):
    mc = atom_coords(metal)
    dc = [atom_coords(d) for d in donors]
    angles = []
    for i in range(len(dc)):
        for j in range(i+1, len(dc)):
            angles.append(round(angle_deg(dc[i]-mc, dc[j]-mc), 2))
    angles.sort()
    return angles

def compute_md_distances(metal, donors):
    mc = atom_coords(metal)
    return sorted([round(np.linalg.norm(atom_coords(d)-mc), 4) for d in donors])

def classify_pattern(angles):
    n_trans = sum(1 for a in angles if a >= 160)
    if n_trans == 2:
        return "SP"
    elif n_trans == 0 and all(90 < a < 130 for a in angles):
        return "Td"
    elif n_trans == 1:
        return "seesaw"
    else:
        return "distorted"

# ══════════════════════════════════════════════════════════════
# XYZ helpers
# ══════════════════════════════════════════════════════════════

def write_xyz(atoms, comment, fp):
    with open(fp, "w") as f:
        f.write(f"{len(atoms)}\n{comment}\n")
        for e,x,y,z in atoms:
            f.write(f"{e:3s} {x:12.6f} {y:12.6f} {z:12.6f}\n")

def mol_no_H(mol):
    return [(a.atomic_symbol, a.coordinates.x, a.coordinates.y, a.coordinates.z)
            for a in mol.atoms if a.coordinates and a.atomic_symbol != "H"]

def first_sphere_atoms(mol, metal, donors):
    al = list(mol.atoms); keep = {al.index(metal)}
    for d in donors:
        di = al.index(d); keep.add(di)
        for b in d.bonds:
            o = b.atoms[0] if b.atoms[1]==d else b.atoms[1]
            if o in al: keep.add(al.index(o))
    return [(al[i].atomic_symbol, al[i].coordinates.x, al[i].coordinates.y, al[i].coordinates.z)
            for i in sorted(keep) if al[i].atomic_symbol!="H" and al[i].coordinates]

def polyhedron_atoms(mol, metal, donors):
    r = [(metal.atomic_symbol, metal.coordinates.x, metal.coordinates.y, metal.coordinates.z)]
    for d in donors:
        if d.coordinates:
            r.append((d.atomic_symbol, d.coordinates.x, d.coordinates.y, d.coordinates.z))
    return r

# ══════════════════════════════════════════════════════════════
# 1. Load GAZKUN family from CN4 atlas
# ══════════════════════════════════════════════════════════════

print("=" * 60)
print("Fig. 5D — GAZKUN Family Production Export")
print("=" * 60)

cn4_all = []
family_atlas = []
with open(SRC / "cn4_sp_td_atlas.csv") as f:
    for r in csv.DictReader(f):
        cn4_all.append(r)
        if r["L3_hash"] == GAZKUN_L3:
            family_atlas.append(r)

print(f"\nGAZKUN family: {len(family_atlas)} members")

# ══════════════════════════════════════════════════════════════
# 2. Compute geometry for each member via CSD API
# ══════════════════════════════════════════════════════════════

csd = EntryReader("CSD")
member_data = []

for ar in family_atlas:
    rc = ar["refcode"]
    s_sp = float(ar["S_SP"]); s_td = float(ar["S_Td"])
    try:
        entry = csd.entry(rc); mol = entry.molecule
        metal, donors = get_metal_and_donors(mol)
        if metal is None or len(donors) != 4:
            print(f"  SKIP {rc}: metal/donors issue")
            continue
        angles = compute_dmd_angles(metal, donors)
        dists = compute_md_distances(metal, donors)
        pattern = classify_pattern(angles)
        delta = abs(s_sp - s_td)
        scalar = s_sp - s_td

        member_data.append({
            "refcode": rc,
            "metal_elem": metal.atomic_symbol,
            "donor_elems": [d.atomic_symbol for d in donors],
            "S_SP": s_sp, "S_Td": s_td,
            "log_S_SP": math.log10(s_sp + 0.1),
            "log_S_Td": math.log10(s_td + 0.1),
            "delta": delta,
            "trajectory_scalar": scalar,
            "is_boundary": delta < 1.0,
            "angles": angles,
            "dists": dists,
            "pattern": pattern,
            "n_trans": sum(1 for a in angles if a >= 160),
            "max_angle": max(angles),
            "min_angle": min(angles),
            "L0_hash": ar["L0_hash"],
            "L1_hash": ar["L1_hash"],
            "L3_hash": ar["L3_hash"],
            "donor_set": ar["donor_set"],
            "denticity_pattern": ar["denticity_pattern"],
        })
    except Exception as e:
        print(f"  ERROR {rc}: {e}")

member_data.sort(key=lambda m: m["trajectory_scalar"])
n_members = len(member_data)

print(f"\n{'='*60}")
print(f"GAZKUN family: {n_members} members, sorted by trajectory scalar")
print(f"{'='*60}")
for i, m in enumerate(member_data):
    bd = "BD" if m["is_boundary"] else "  "
    print(f"  {i+1}. {m['refcode']:12s} S_SP={m['S_SP']:6.2f} S_Td={m['S_Td']:6.2f} "
          f"Δ={m['delta']:5.2f} {bd} pattern={m['pattern']:10s} "
          f"angles={m['angles']}")

# ══════════════════════════════════════════════════════════════
# 3. Select 5 representative states
# ══════════════════════════════════════════════════════════════

# Strategy: pick at ~0%, 25%, 50%, 75%, 100% quantiles of trajectory
# Ensure boundary member is included
indices = [0]  # start
n = len(member_data)
q25 = max(1, round(n * 0.25))
q50 = max(1, round(n * 0.50))
q75 = min(n-2, round(n * 0.75))
end_idx = n - 1

# Find boundary member closest to mid-trajectory
boundary_members = [(i, m) for i, m in enumerate(member_data) if m["is_boundary"]]

candidate_indices = sorted(set([0, q25, q50, q75, end_idx]))
# Ensure at least one boundary member
if boundary_members:
    best_bd_idx = boundary_members[0][0]
    if best_bd_idx not in candidate_indices:
        # Replace the closest candidate
        dists_to_bd = [(abs(ci - best_bd_idx), ci) for ci in candidate_indices if ci != 0 and ci != end_idx]
        if dists_to_bd:
            dists_to_bd.sort()
            replace = dists_to_bd[0][1]
            candidate_indices.remove(replace)
            candidate_indices.append(best_bd_idx)
            candidate_indices.sort()

# Trim to exactly 5
if len(candidate_indices) > 5:
    candidate_indices = candidate_indices[:5]
while len(candidate_indices) < 5 and len(candidate_indices) < n:
    for ci in range(n):
        if ci not in candidate_indices:
            candidate_indices.append(ci)
            candidate_indices.sort()
            break
    if len(candidate_indices) >= 5:
        break

selected = [member_data[i] for i in candidate_indices[:5]]

# Assign reasons
reasons = {}
for i, m in enumerate(selected):
    idx = candidate_indices[i]
    if idx == 0:
        reasons[m["refcode"]] = "trajectory start (boundary, most Td-like)"
    elif idx == end_idx:
        reasons[m["refcode"]] = "trajectory end (most SP-like)"
    elif m["is_boundary"]:
        reasons[m["refcode"]] = f"boundary state (delta={m['delta']:.2f})"
    elif m["pattern"] == "Td":
        reasons[m["refcode"]] = "near-ideal tetrahedral"
    elif m["pattern"] == "SP":
        reasons[m["refcode"]] = "near-square-planar intermediate"
    else:
        reasons[m["refcode"]] = "intermediate geometry"

print(f"\n[3] Selected 5 render states:")
for i, m in enumerate(selected):
    print(f"  render {i+1}: {m['refcode']:12s}  pattern={m['pattern']:10s}  "
          f"Δ={m['delta']:.2f}  → {reasons[m['refcode']]}")

# ══════════════════════════════════════════════════════════════
# 4. Write trajectory points CSV (all members)
# ══════════════════════════════════════════════════════════════

selected_refs = {m["refcode"] for m in selected}
traj_fields = [
    "point_rank","refcode","family_label","L3_key_hash","L0_key_hash",
    "L1_shape_id","metal","CN",
    "S_SP","S_Td","log_S_SP","log_S_Td",
    "delta_SP_Td","shape_boundary_flag","angle_pattern","trajectory_scalar",
    "dmd_angles","md_distances","n_trans_angles","max_angle","min_angle",
    "selected_for_render",
]

with open(OUT / "fig5D_GAZKUN_trajectory_points.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=traj_fields)
    w.writeheader()
    for i, m in enumerate(member_data):
        w.writerow({
            "point_rank": i+1,
            "refcode": m["refcode"],
            "family_label": FAMILY_LABEL,
            "L3_key_hash": m["L3_hash"],
            "L0_key_hash": m["L0_hash"],
            "L1_shape_id": m["L1_hash"],
            "metal": FAMILY_METAL,
            "CN": FAMILY_CN,
            "S_SP": f'{m["S_SP"]:.4f}',
            "S_Td": f'{m["S_Td"]:.4f}',
            "log_S_SP": f'{m["log_S_SP"]:.4f}',
            "log_S_Td": f'{m["log_S_Td"]:.4f}',
            "delta_SP_Td": f'{m["delta"]:.4f}',
            "shape_boundary_flag": m["is_boundary"],
            "angle_pattern": m["pattern"],
            "trajectory_scalar": f'{m["trajectory_scalar"]:.4f}',
            "dmd_angles": ";".join(f"{a:.1f}" for a in m["angles"]),
            "md_distances": ";".join(f"{d:.4f}" for d in m["dists"]),
            "n_trans_angles": m["n_trans"],
            "max_angle": m["max_angle"],
            "min_angle": m["min_angle"],
            "selected_for_render": m["refcode"] in selected_refs,
        })

print(f"\n  → fig5D_GAZKUN_trajectory_points.csv ({n_members} rows)")

# ══════════════════════════════════════════════════════════════
# 5. Write selected render states CSV
# ══════════════════════════════════════════════════════════════

render_fields = [
    "render_order","refcode","family_label","L3_key_hash","L0_key_hash",
    "L1_shape_id","metal","CN",
    "S_SP","S_Td","delta_SP_Td","angle_pattern",
    "dmd_angles","trajectory_scalar","selection_reason",
    "xyz_full","xyz_first_sphere","xyz_polyhedron",
]

with open(OUT / "fig5D_GAZKUN_render_states.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=render_fields)
    w.writeheader()
    for i, m in enumerate(selected):
        sn = f"{i+1:02d}"
        w.writerow({
            "render_order": i+1,
            "refcode": m["refcode"],
            "family_label": FAMILY_LABEL,
            "L3_key_hash": m["L3_hash"],
            "L0_key_hash": m["L0_hash"],
            "L1_shape_id": m["L1_hash"],
            "metal": FAMILY_METAL,
            "CN": FAMILY_CN,
            "S_SP": f'{m["S_SP"]:.4f}',
            "S_Td": f'{m["S_Td"]:.4f}',
            "delta_SP_Td": f'{m["delta"]:.4f}',
            "angle_pattern": m["pattern"],
            "dmd_angles": ";".join(f"{a:.1f}" for a in m["angles"]),
            "trajectory_scalar": f'{m["trajectory_scalar"]:.4f}',
            "selection_reason": reasons[m["refcode"]],
            "xyz_full": f"GAZKUN_state{sn}_full.xyz",
            "xyz_first_sphere": f"GAZKUN_state{sn}_first_sphere.xyz",
            "xyz_polyhedron": f"GAZKUN_state{sn}_polyhedron.xyz",
        })

print(f"  → fig5D_GAZKUN_render_states.csv (5 states)")

# ══════════════════════════════════════════════════════════════
# 6. Export XYZ files
# ══════════════════════════════════════════════════════════════

print(f"\n[6] Exporting XYZ files...")
for i, m in enumerate(selected):
    rc = m["refcode"]; sn = f"{i+1:02d}"
    prefix = f"GAZKUN_state{sn}"
    try:
        entry = csd.entry(rc); mol = entry.molecule
        metal, donors = get_metal_and_donors(mol)
        if metal is None:
            print(f"    state{sn} {rc}: no metal, skip"); continue

        full_a = mol_no_H(mol)
        fs_a = first_sphere_atoms(mol, metal, donors)
        poly_a = polyhedron_atoms(mol, metal, donors)

        write_xyz(full_a,
                  f"CSD={rc} | GAZKUN state{sn} | {m['pattern']} | full no H",
                  OUT / f"{prefix}_full.xyz")
        write_xyz(fs_a,
                  f"CSD={rc} | GAZKUN state{sn} | {m['pattern']} | first_sphere | angles={m['angles']}",
                  OUT / f"{prefix}_first_sphere.xyz")
        write_xyz(poly_a,
                  f"CSD={rc} | GAZKUN state{sn} | {m['pattern']} | Cu+4N polyhedron",
                  OUT / f"{prefix}_polyhedron.xyz")

        print(f"    state{sn}: {rc:12s}  pattern={m['pattern']:10s}  "
              f"full={len(full_a)} fs={len(fs_a)} poly={len(poly_a)}")
    except Exception as e:
        print(f"    state{sn}: {rc} ERROR: {e}")

# ══════════════════════════════════════════════════════════════
# 7. Trajectory summary JSON
# ══════════════════════════════════════════════════════════════

unique_L0 = len(set(m["L0_hash"] for m in member_data))
unique_L1 = len(set(m["L1_hash"] for m in member_data))
span = member_data[-1]["trajectory_scalar"] - member_data[0]["trajectory_scalar"]

traj_summary = {
    "family_label": FAMILY_LABEL,
    "metal": FAMILY_METAL,
    "CN": FAMILY_CN,
    "donor_set": "N4",
    "n_family_members": n_members,
    "trajectory_span": round(span, 2),
    "trajectory_axis_definition": "S_SP - S_Td",
    "n_unique_L0": unique_L0,
    "n_unique_L1": unique_L1,
    "crosses_shape_boundary": any(m["is_boundary"] for m in member_data),
    "selected_render_count": 5,
    "L3_key_hash": GAZKUN_L3,
    "pattern_sequence": [m["pattern"] for m in member_data],
    "visual_justification": (
        "GAZKUN family shows dramatic visual geometry change: "
        "GAZKUN (boundary, S_Td=0.44, flattened Td-like) → "
        "LOYZII (near-ideal Td, all angles ~109°) → "
        "TUGKOX (near-ideal SP, two trans angles ~176°). "
        "This makes the SP↔Td interconversion visually obvious "
        "in both polyhedron and first_sphere renderings."
    ),
    "why_not_IKEXEC": (
        "IKEXEC/MIGTEE family, despite high CShM span, "
        "shows minimal visual difference in first_sphere XYZ "
        "because all members have similar Ni-donor angles near "
        "the distorted region. GAZKUN shows actual SP↔Td pattern "
        "change with angle RMSD=42.5° vs IKEXEC's ~5°."
    ),
    "representative_refcodes": [m["refcode"] for m in member_data],
    "recommended_render_file": "first_sphere.xyz",
    "license_note": LICENSE,
}
json.dump(traj_summary, open(OUT / "fig5D_GAZKUN_trajectory_summary.json", "w"), indent=2)
print(f"\n  → fig5D_GAZKUN_trajectory_summary.json")

# ══════════════════════════════════════════════════════════════
# Verification
# ══════════════════════════════════════════════════════════════

req = [
    "fig5D_GAZKUN_trajectory_points.csv",
    "fig5D_GAZKUN_render_states.csv",
    "fig5D_GAZKUN_trajectory_summary.json",
    "top_visual_family_candidates.csv",
]
for i in range(5):
    sn = f"{i+1:02d}"
    req += [f"GAZKUN_state{sn}_full.xyz",
            f"GAZKUN_state{sn}_first_sphere.xyz",
            f"GAZKUN_state{sn}_polyhedron.xyz"]

print(f"\n{'='*60}\nVerification\n{'='*60}")
ok = 0
for fn in req:
    p = OUT / fn
    if p.exists():
        print(f"  ✓ {fn:55s} {os.path.getsize(p):>8,} B")
        ok += 1
    else:
        print(f"  ✗ {fn:55s} MISSING")
print(f"\n{ok}/{len(req)} present.")
print("=" * 60)
