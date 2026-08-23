#!/usr/bin/env python3
"""
export_fig5d_maintext_subset.py
================================
Main-text-ready subset for Fig. 5D:
  5 top-span CN4 families (diverse metals) + IKEXEC 5-state thumbnail index.

Two steps:
  1) CSD API step: export IKEXEC 5×3 XYZ (run with 1701 env)
  2) Pure-Python CSV/JSON (no CSD needed)

Run:
  /data/miniconda3/envs/1701/bin/python scripts/export_fig5d_maintext_subset.py

Output → revision_results/fig5d_family_trajectory/
"""
from __future__ import annotations
import csv, json, math, os, sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))
from coordrep_tools.csd_adapter import _count_metals, _get_donor_neighbors
from ccdc.io import EntryReader

BASE = Path(__file__).parent.parent
SRC  = BASE / "revision_results" / "csd_pathfinder_full"
OUT  = BASE / "revision_results" / "fig5d_family_trajectory"
os.makedirs(OUT, exist_ok=True)

LICENSE = "CSD-derived; internal figure rendering only."

# ── Selected 5 families (diverse metals) ─────────────────────
SELECTED_FAMILIES = {
    "463888590564": {"label": "LINMOL",   "metal": "Pd", "rank": 1},
    "371011783775": {"label": "IKEXEC",   "metal": "Ni", "rank": 2},
    "918517807621": {"label": "TOPHIP",   "metal": "Cu", "rank": 3},
    "99885397104":  {"label": "YOSNAV01", "metal": "Pt", "rank": 4},
    "426961457385": {"label": "COSALP10", "metal": "Co", "rank": 5},
}

LINMOL_L3  = "463888590564"
IKEXEC_L3  = "371011783775"

# IKEXEC 5 selected refcodes (drop MIGTEE, keep boundary MIGTEE01)
IKEXEC_SELECTED_5 = ["IKEXEC", "MIGTEE01", "SONXUP", "MIGTEE02", "MIGTEE03"]
IKEXEC_SELECT_REASONS = {
    "IKEXEC":   "trajectory start (most SP-like)",
    "MIGTEE01": "boundary state (delta = 0.90)",
    "SONXUP":   "early Td-leaning intermediate",
    "MIGTEE02": "late intermediate",
    "MIGTEE03": "trajectory end (most Td-like)",
}

# ── XYZ helpers ──────────────────────────────────────────────
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


print("=" * 60)
print("Fig. 5D — Main-text subset export")
print("=" * 60)

# ═══════════════════════════════════════════════════════════
# 1. Load CN4 atlas for selected families
# ═══════════════════════════════════════════════════════════
atlas_by_l3 = defaultdict(list)
with open(SRC / "cn4_sp_td_atlas.csv") as f:
    for r in csv.DictReader(f):
        if r["L3_hash"] in SELECTED_FAMILIES:
            atlas_by_l3[r["L3_hash"]].append(r)

total_pts = sum(len(v) for v in atlas_by_l3.values())
print(f"\n[1] Loaded {total_pts} atlas points across {len(atlas_by_l3)} families")

# ═══════════════════════════════════════════════════════════
# 2. fig5D_maintext_overlay_subset.csv
# ═══════════════════════════════════════════════════════════
fields = [
    "family_rank", "family_label", "L3_key_hash", "metal", "CN",
    "n_family_members", "trajectory_span",
    "refcode", "L0_key_hash", "L1_shape_id",
    "S_SP", "S_Td", "log_S_SP", "log_S_Td",
    "delta_SP_Td", "shape_boundary_flag", "shape_label",
    "trajectory_scalar",
    "is_LINMOL", "is_IKEXEC_thumbnail",
    "render_order", "selection_reason",
]

ikexec_set = set(IKEXEC_SELECTED_5)
rows_written = 0

with open(OUT / "fig5D_maintext_overlay_subset.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fields)
    w.writeheader()

    for l3_hash in sorted(SELECTED_FAMILIES, key=lambda h: SELECTED_FAMILIES[h]["rank"]):
        info = SELECTED_FAMILIES[l3_hash]
        members = atlas_by_l3.get(l3_hash, [])
        for m in members:
            m["_s_sp"] = float(m["S_SP"])
            m["_s_td"] = float(m["S_Td"])
            m["_scalar"] = m["_s_sp"] - m["_s_td"]
            m["_delta"] = abs(m["_s_sp"] - m["_s_td"])
        members.sort(key=lambda m: m["_scalar"])

        for m in members:
            if m["_delta"] < 1.0:
                slabel = "boundary"
            elif m["_s_sp"] < m["_s_td"]:
                slabel = "SP-like"
            else:
                slabel = "Td-like"

            is_thumb = (l3_hash == IKEXEC_L3 and m["refcode"] in ikexec_set)
            rorder = ""
            sreason = ""
            if is_thumb:
                rorder = IKEXEC_SELECTED_5.index(m["refcode"]) + 1
                sreason = IKEXEC_SELECT_REASONS[m["refcode"]]

            w.writerow({
                "family_rank": info["rank"],
                "family_label": info["label"],
                "L3_key_hash": l3_hash,
                "metal": info["metal"],
                "CN": 4,
                "n_family_members": len(members),
                "trajectory_span": round(max(m2["_scalar"] for m2 in members) -
                                         min(m2["_scalar"] for m2 in members), 2),
                "refcode": m["refcode"],
                "L0_key_hash": m["L0_hash"],
                "L1_shape_id": m["L1_hash"],
                "S_SP": f'{m["_s_sp"]:.4f}',
                "S_Td": f'{m["_s_td"]:.4f}',
                "log_S_SP": f'{math.log10(m["_s_sp"]+0.1):.4f}',
                "log_S_Td": f'{math.log10(m["_s_td"]+0.1):.4f}',
                "delta_SP_Td": f'{m["_delta"]:.4f}',
                "shape_boundary_flag": m["_delta"] < 1.0,
                "shape_label": slabel,
                "trajectory_scalar": f'{m["_scalar"]:.4f}',
                "is_LINMOL": l3_hash == LINMOL_L3,
                "is_IKEXEC_thumbnail": is_thumb,
                "render_order": rorder,
                "selection_reason": sreason,
            })
            rows_written += 1

print(f"\n[2] fig5D_maintext_overlay_subset.csv — {rows_written} rows, {len(SELECTED_FAMILIES)} families")

# ═══════════════════════════════════════════════════════════
# 3. Export IKEXEC 5×3 XYZ via CSD API
# ═══════════════════════════════════════════════════════════
print("\n[3] Exporting IKEXEC 5-state XYZ thumbnails...")

csd = EntryReader("CSD")

render_index = []
for i, rc in enumerate(IKEXEC_SELECTED_5):
    sn = f"{i+1:02d}"
    prefix = f"IKEXEC_state{sn}"
    print(f"    state{sn}: {rc} ...", end=" ")
    try:
        entry = csd.entry(rc)
        mol = entry.molecule

        full_a  = mol_no_H(mol)
        fs_a    = first_sphere(mol)
        poly_a  = polyhedron(mol)

        write_xyz(full_a,
                  f"CSD_code={rc} | IKEXEC family state{sn} | no H",
                  OUT / f"{prefix}_full.xyz")
        write_xyz(fs_a,
                  f"CSD_code={rc} | first_sphere | Ni + donors + backbone | no H",
                  OUT / f"{prefix}_first_sphere.xyz")
        write_xyz(poly_a,
                  f"CSD_code={rc} | polyhedron | Ni + donors",
                  OUT / f"{prefix}_polyhedron.xyz")

        # atlas data for this refcode
        atlas_row = None
        for m in atlas_by_l3.get(IKEXEC_L3, []):
            if m["refcode"] == rc:
                atlas_row = m; break

        render_index.append({
            "render_order": i + 1,
            "refcode": rc,
            "family_label": "IKEXEC",
            "metal": "Ni",
            "CN": 4,
            "S_SP": f'{atlas_row["_s_sp"]:.4f}' if atlas_row else "",
            "S_Td": f'{atlas_row["_s_td"]:.4f}' if atlas_row else "",
            "delta": f'{atlas_row["_delta"]:.4f}' if atlas_row else "",
            "shape_label": ("boundary" if atlas_row and atlas_row["_delta"] < 1.0
                            else "SP-like" if atlas_row and atlas_row["_s_sp"] < atlas_row["_s_td"]
                            else "Td-like"),
            "selection_reason": IKEXEC_SELECT_REASONS[rc],
            "xyz_full": f"{prefix}_full.xyz",
            "xyz_first_sphere": f"{prefix}_first_sphere.xyz",
            "xyz_polyhedron": f"{prefix}_polyhedron.xyz",
            "n_atoms_full": len(full_a),
            "n_atoms_first_sphere": len(fs_a),
            "n_atoms_polyhedron": len(poly_a),
        })
        print(f"full={len(full_a)} fs={len(fs_a)} poly={len(poly_a)}")
    except Exception as e:
        print(f"ERROR: {e}")

# ═══════════════════════════════════════════════════════════
# 4. Render index CSV + JSON
# ═══════════════════════════════════════════════════════════
ri_fields = [
    "render_order","refcode","family_label","metal","CN",
    "S_SP","S_Td","delta","shape_label","selection_reason",
    "xyz_full","xyz_first_sphere","xyz_polyhedron",
    "n_atoms_full","n_atoms_first_sphere","n_atoms_polyhedron",
]
with open(OUT / "fig5D_IKEXEC_render_index.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=ri_fields)
    w.writeheader()
    w.writerows(render_index)

json.dump(render_index, open(OUT / "fig5D_IKEXEC_render_index.json", "w"), indent=2)
print(f"\n[4] fig5D_IKEXEC_render_index.csv/json — {len(render_index)} states")

# ═══════════════════════════════════════════════════════════
# 5. Maintext summary JSON
# ═══════════════════════════════════════════════════════════
msummary = {
    "description": "Main-text Fig. 5D: 5 selected CN4 L3 families + IKEXEC 5-state thumbnails",
    "selected_families": [
        {"rank": info["rank"], "label": info["label"], "metal": info["metal"],
         "L3_hash": h, "n_members": len(atlas_by_l3.get(h, [])),
         "span": round(max(m["_scalar"] for m in atlas_by_l3[h]) -
                        min(m["_scalar"] for m in atlas_by_l3[h]), 2)
                 if atlas_by_l3.get(h) else 0}
        for h, info in sorted(SELECTED_FAMILIES.items(), key=lambda x: x[1]["rank"])
    ],
    "thumbnail_family": "IKEXEC",
    "thumbnail_metal": "Ni",
    "thumbnail_CN": 4,
    "thumbnail_n_states": 5,
    "thumbnail_span": 8.00,
    "nontrivial_L3_families": 6358,
    "multiple_L1_shape_families": 3739,
    "shape_boundary_crossing_families": 1810,
    "top_family": "LINMOL",
    "top_family_span": 10.10,
    "recommended_render_file": "first_sphere.xyz",
    "license_note": LICENSE,
}
json.dump(msummary, open(OUT / "fig5D_maintext_summary.json", "w"), indent=2)
print(f"[5] fig5D_maintext_summary.json")

# ═══════════════════════════════════════════════════════════
# Verification
# ═══════════════════════════════════════════════════════════
req = [
    "fig5D_maintext_overlay_subset.csv",
    "fig5D_IKEXEC_render_index.csv",
    "fig5D_IKEXEC_render_index.json",
    "fig5D_maintext_summary.json",
]
for i in range(5):
    sn = f"{i+1:02d}"
    req += [f"IKEXEC_state{sn}_full.xyz",
            f"IKEXEC_state{sn}_first_sphere.xyz",
            f"IKEXEC_state{sn}_polyhedron.xyz"]

print(f"\n{'='*60}\nNew files verification\n{'='*60}")
ok = 0
for fn in req:
    p = OUT / fn
    if p.exists():
        print(f"  ✓ {fn:50s} {os.path.getsize(p):>8,} B")
        ok += 1
    else:
        print(f"  ✗ {fn:50s} MISSING")
print(f"\n{ok}/{len(req)} files present.")
print("=" * 60)
