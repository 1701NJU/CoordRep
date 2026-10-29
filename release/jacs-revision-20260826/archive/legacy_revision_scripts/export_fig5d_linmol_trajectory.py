#!/usr/bin/env python3
"""
export_fig5d_linmol_trajectory.py
==================================
Export figure-ready data for Fig. 5D: LINMOL family trajectory
in the CN=4 SP–Td CShM space.

Run with the 1701 conda env (has CSD Python API):
  /data/miniconda3/envs/1701/bin/python scripts/export_fig5d_linmol_trajectory.py

Outputs → revision_results/fig5d_family_trajectory/

License: XYZ files are CSD-derived, for internal figure rendering ONLY.
         Do NOT redistribute raw coordinates.
"""

from __future__ import annotations

import csv
import json
import math
import os
import random
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent))

from ccdc.io import EntryReader

from coordrep.io.tmqm_reader import TRANSITION_METALS
from coordrep.graph.neighbors import DONOR_ELEMENTS
from coordrep_tools.csd_adapter import _count_metals, _get_donor_neighbors

BASE = Path(__file__).parent.parent
SRC = BASE / "revision_results" / "csd_pathfinder_full"
OUT = BASE / "revision_results" / "fig5d_family_trajectory"
os.makedirs(OUT, exist_ok=True)

LICENSE_NOTE = (
    "CSD-derived coordinates for internal figure rendering only; "
    "do not redistribute raw XYZ."
)

# ══════════════════════════════════════════════════════════════
# LINMOL family definition
# ══════════════════════════════════════════════════════════════

LINMOL_L3_HASH = 463888590564
FAMILY_LABEL = "LINMOL"
FAMILY_METAL = "Pd"
FAMILY_CN = 4

# ══════════════════════════════════════════════════════════════
# XYZ I/O helpers (reused from export_fig5c_case_xyz.py)
# ══════════════════════════════════════════════════════════════

def write_xyz(atoms, comment, filepath):
    """Write a list of (element, x, y, z) tuples to XYZ format."""
    with open(filepath, "w") as f:
        f.write(f"{len(atoms)}\n")
        f.write(f"{comment}\n")
        for elem, x, y, z in atoms:
            f.write(f"{elem:3s} {x:12.6f} {y:12.6f} {z:12.6f}\n")


def csd_mol_to_atom_list(mol):
    """CSD molecule → list of (element, x, y, z)."""
    atoms = []
    for a in mol.atoms:
        c = a.coordinates
        if c is None:
            continue
        atoms.append((a.atomic_symbol, c.x, c.y, c.z))
    return atoms


def filter_no_H(atoms):
    """Remove hydrogen atoms."""
    return [(e, x, y, z) for e, x, y, z in atoms if e != "H"]


def get_metal_and_donors(mol):
    """Return (metal_atom, donor_atoms, metal_idx, donor_indices)."""
    metals, n_metals = _count_metals(mol)
    if n_metals == 0:
        return None, [], -1, []
    metal_atom = metals[0]
    donors = _get_donor_neighbors(mol, metal_atom)
    atom_list = list(mol.atoms)
    metal_idx = atom_list.index(metal_atom) if metal_atom in atom_list else -1
    donor_indices = []
    for d in donors:
        if d in atom_list:
            donor_indices.append(atom_list.index(d))
    return metal_atom, donors, metal_idx, donor_indices


def first_sphere_atoms(mol):
    """Metal + donor atoms + atoms directly bonded to donors (backbone), no H."""
    metal_atom, donors, metal_idx, donor_indices = get_metal_and_donors(mol)
    if metal_atom is None:
        return []
    keep_set = set()
    atom_list = list(mol.atoms)
    keep_set.add(metal_idx)
    for di in donor_indices:
        keep_set.add(di)
    for d_atom in donors:
        for b in d_atom.bonds:
            other = b.atoms[0] if b.atoms[1] == d_atom else b.atoms[1]
            if other in atom_list:
                idx = atom_list.index(other)
                keep_set.add(idx)
    result = []
    for idx in sorted(keep_set):
        a = atom_list[idx]
        if a.atomic_symbol == "H":
            continue
        c = a.coordinates
        if c is None:
            continue
        result.append((a.atomic_symbol, c.x, c.y, c.z))
    return result


def polyhedron_atoms(mol):
    """Metal + donor atoms only (for polyhedron rendering)."""
    metal_atom, donors, metal_idx, donor_indices = get_metal_and_donors(mol)
    if metal_atom is None:
        return []
    atom_list = list(mol.atoms)
    result = []
    for idx in [metal_idx] + donor_indices:
        a = atom_list[idx]
        c = a.coordinates
        if c is None:
            continue
        result.append((a.atomic_symbol, c.x, c.y, c.z))
    return result


# ══════════════════════════════════════════════════════════════
# A. Load LINMOL family from CN4 atlas
# ══════════════════════════════════════════════════════════════

print("=" * 60)
print("Fig. 5D — LINMOL Family Trajectory Export")
print("=" * 60)

# Read CN4 atlas for LINMOL members
cn4_atlas_path = SRC / "cn4_sp_td_atlas.csv"
linmol_rows = []
cn4_all_rows = []

with open(cn4_atlas_path) as f:
    reader = csv.DictReader(f)
    for row in reader:
        cn4_all_rows.append(row)
        if row["L3_hash"] == str(LINMOL_L3_HASH):
            linmol_rows.append(row)

print(f"\n[A] LINMOL family: {len(linmol_rows)} members found (L3_hash={LINMOL_L3_HASH})")
for r in linmol_rows:
    print(f"    {r['refcode']:12s}  S_SP={float(r['S_SP']):6.2f}  S_Td={float(r['S_Td']):6.2f}  L0={r['L0_hash']}  L1={r['L1_hash']}")

# ══════════════════════════════════════════════════════════════
# B. Compute trajectory scalar and build trajectory points
# ══════════════════════════════════════════════════════════════

# trajectory_scalar = S_SP - S_Td  (measures position along SP→Td axis)
# Positive = more Td-distorted; negative = more SP-distorted
# Alternative: use the principal axis. With only 3 points we use the simpler definition.

for r in linmol_rows:
    s_sp = float(r["S_SP"])
    s_td = float(r["S_Td"])
    r["_s_sp"] = s_sp
    r["_s_td"] = s_td
    r["_log_sp"] = math.log10(s_sp + 0.1)
    r["_log_td"] = math.log10(s_td + 0.1)
    r["_delta"] = abs(s_sp - s_td)
    r["_trajectory_scalar"] = s_sp - s_td  # SP–Td difference

# Sort by trajectory_scalar
linmol_rows.sort(key=lambda r: r["_trajectory_scalar"])

# Determine shape labels
for r in linmol_rows:
    if r["_s_sp"] < r["_s_td"]:
        r["_shape_label"] = "SP-like"
    elif abs(r["_s_sp"] - r["_s_td"]) < 1.0:
        r["_shape_label"] = "boundary"
    else:
        r["_shape_label"] = "Td-like"
    r["_boundary_flag"] = r["_delta"] < 1.0

# Mark selected for render
# With only 3 members, select ALL 3
n_members = len(linmol_rows)
for r in linmol_rows:
    r["_selected"] = True

print(f"\n[B] Trajectory scalar: S_SP - S_Td")
for i, r in enumerate(linmol_rows):
    print(f"    rank {i+1}: {r['refcode']:12s}  scalar={r['_trajectory_scalar']:+6.2f}  "
          f"shape={r['_shape_label']}  delta={r['_delta']:.2f}")

# ── B.1: Trajectory points CSV ──────────────────────────────

traj_fieldnames = [
    "point_rank", "refcode", "family_label", "L3_key_hash", "L0_key_hash",
    "L1_shape_id", "metal", "CN", "S_SP", "S_Td", "log_S_SP", "log_S_Td",
    "shape_boundary_flag", "shape_label", "trajectory_scalar", "selected_for_render",
]

with open(OUT / "fig5D_LINMOL_trajectory_points.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=traj_fieldnames)
    w.writeheader()
    for i, r in enumerate(linmol_rows):
        w.writerow({
            "point_rank": i + 1,
            "refcode": r["refcode"],
            "family_label": FAMILY_LABEL,
            "L3_key_hash": r["L3_hash"],
            "L0_key_hash": r["L0_hash"],
            "L1_shape_id": r["L1_hash"],
            "metal": r["metal"],
            "CN": FAMILY_CN,
            "S_SP": f"{r['_s_sp']:.4f}",
            "S_Td": f"{r['_s_td']:.4f}",
            "log_S_SP": f"{r['_log_sp']:.4f}",
            "log_S_Td": f"{r['_log_td']:.4f}",
            "shape_boundary_flag": r["_boundary_flag"],
            "shape_label": r["_shape_label"],
            "trajectory_scalar": f"{r['_trajectory_scalar']:.4f}",
            "selected_for_render": r["_selected"],
        })

print(f"\n  → fig5D_LINMOL_trajectory_points.csv written ({n_members} rows)")

# ── B.2: Trajectory summary JSON ────────────────────────────

unique_L0 = len(set(r["L0_hash"] for r in linmol_rows))
unique_L1 = len(set(r["L1_hash"] for r in linmol_rows))
crosses_boundary = any(r["_boundary_flag"] for r in linmol_rows)
span = max(r["_trajectory_scalar"] for r in linmol_rows) - min(r["_trajectory_scalar"] for r in linmol_rows)

summary = {
    "family_label": FAMILY_LABEL,
    "metal": FAMILY_METAL,
    "CN": FAMILY_CN,
    "n_family_members": n_members,
    "trajectory_span": round(span, 2),
    "trajectory_axis_definition": "S_SP - S_Td (continuous shape measure difference; SP→Td axis)",
    "n_unique_L0": unique_L0,
    "n_unique_L1": unique_L1,
    "crosses_shape_boundary": crosses_boundary,
    "selected_render_count": sum(1 for r in linmol_rows if r["_selected"]),
    "L3_key_hash": LINMOL_L3_HASH,
    "representative_refcodes": [r["refcode"] for r in linmol_rows],
    "note": (
        "LINMOL family has 3 members (not 5). All 3 are selected for rendering. "
        "Each member has a distinct L0 (geometry-hash) and L1 (shape-ID), "
        "demonstrating clear geometry-state diversity within a single connectivity family."
    ),
    "license_note": LICENSE_NOTE,
}

json.dump(summary, open(OUT / "fig5D_LINMOL_trajectory_summary.json", "w"), indent=2)
print(f"  → fig5D_LINMOL_trajectory_summary.json written")


# ══════════════════════════════════════════════════════════════
# C. Selected render states CSV
# ══════════════════════════════════════════════════════════════

# Selection reasons for 3 members
selection_reasons = {
    0: "trajectory start (lowest S_SP - S_Td)",
    1: "near-boundary midpoint",
    2: "trajectory end (highest S_SP - S_Td)",
}
if n_members == 3:
    # Adjust reasons based on actual geometry
    for i, r in enumerate(linmol_rows):
        if r["_boundary_flag"]:
            selection_reasons[i] = "boundary-like state (delta < 1.0)"

render_fieldnames = [
    "render_order", "refcode", "family_label", "L3_key_hash", "L0_key_hash",
    "L1_shape_id", "shape_label", "shape_boundary_flag", "S_SP", "S_Td",
    "log_S_SP", "log_S_Td", "trajectory_scalar", "selection_reason",
    "xyz_file", "xyz_first_sphere_file", "xyz_polyhedron_file",
]

selected_rows = [r for r in linmol_rows if r["_selected"]]

with open(OUT / "fig5D_selected_render_states.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=render_fieldnames)
    w.writeheader()
    for i, r in enumerate(selected_rows):
        state_num = f"{i+1:02d}"
        w.writerow({
            "render_order": i + 1,
            "refcode": r["refcode"],
            "family_label": FAMILY_LABEL,
            "L3_key_hash": r["L3_hash"],
            "L0_key_hash": r["L0_hash"],
            "L1_shape_id": r["L1_hash"],
            "shape_label": r["_shape_label"],
            "shape_boundary_flag": r["_boundary_flag"],
            "S_SP": f"{r['_s_sp']:.4f}",
            "S_Td": f"{r['_s_td']:.4f}",
            "log_S_SP": f"{r['_log_sp']:.4f}",
            "log_S_Td": f"{r['_log_td']:.4f}",
            "trajectory_scalar": f"{r['_trajectory_scalar']:.4f}",
            "selection_reason": selection_reasons.get(i, "intermediate"),
            "xyz_file": f"LINMOL_state{state_num}_full.xyz",
            "xyz_first_sphere_file": f"LINMOL_state{state_num}_first_sphere.xyz",
            "xyz_polyhedron_file": f"LINMOL_state{state_num}_polyhedron.xyz",
        })

print(f"\n[C] fig5D_selected_render_states.csv written ({len(selected_rows)} states)")


# ══════════════════════════════════════════════════════════════
# D. Export XYZ files via CSD API
# ══════════════════════════════════════════════════════════════

print("\n[D] Exporting XYZ files via CSD API...")

reader = EntryReader("CSD")

for i, r in enumerate(selected_rows):
    refcode = r["refcode"]
    state_num = f"{i+1:02d}"
    prefix = f"LINMOL_state{state_num}"

    print(f"    State {i+1}: {refcode} ...")

    try:
        entry = reader.entry(refcode)
        mol = entry.molecule

        # Full structure (no H for cleaner render)
        all_atoms = csd_mol_to_atom_list(mol)
        clean_atoms = filter_no_H(all_atoms)
        write_xyz(clean_atoms,
                  f"CSD_code={refcode} | family={FAMILY_LABEL} | state{state_num} | no H",
                  OUT / f"{prefix}_full.xyz")

        # First coordination sphere
        fs_atoms = first_sphere_atoms(mol)
        write_xyz(fs_atoms,
                  f"CSD_code={refcode} | first_sphere | Pd + 4 donors + backbone | no H",
                  OUT / f"{prefix}_first_sphere.xyz")

        # Polyhedron (metal + donors only)
        poly_atoms = polyhedron_atoms(mol)
        write_xyz(poly_atoms,
                  f"CSD_code={refcode} | polyhedron | Pd + donors only",
                  OUT / f"{prefix}_polyhedron.xyz")

        print(f"      full={len(clean_atoms)} atoms, "
              f"first_sphere={len(fs_atoms)} atoms, "
              f"polyhedron={len(poly_atoms)} atoms")

    except Exception as e:
        print(f"      ERROR: {e}")


# ══════════════════════════════════════════════════════════════
# E. CN4 background sample
# ══════════════════════════════════════════════════════════════

print("\n[E] Sampling CN4 background points...")

random.seed(42)
n_bg = min(8000, len(cn4_all_rows))
bg_sample = random.sample(cn4_all_rows, n_bg)

bg_fieldnames = ["refcode", "metal", "CN", "S_SP", "S_Td", "log_S_SP", "log_S_Td"]

with open(OUT / "fig5D_cn4_background_sample.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=bg_fieldnames)
    w.writeheader()
    for r in bg_sample:
        s_sp = float(r["S_SP"])
        s_td = float(r["S_Td"])
        w.writerow({
            "refcode": r["refcode"],
            "metal": r["metal"],
            "CN": 4,
            "S_SP": f"{s_sp:.4f}",
            "S_Td": f"{s_td:.4f}",
            "log_S_SP": f"{math.log10(s_sp + 0.1):.4f}",
            "log_S_Td": f"{math.log10(s_td + 0.1):.4f}",
        })

print(f"  → fig5D_cn4_background_sample.csv written ({n_bg} rows)")


# ══════════════════════════════════════════════════════════════
# F. Family statistics CSV
# ══════════════════════════════════════════════════════════════

stats = [
    ("nontrivial_L3_families", 6358, "6,358 nontrivial L3 families"),
    ("multiple_L1_shape_families", 3739, "3,739 families with multiple L1 shapes"),
    ("shape_boundary_crossing_families", 1810, "1,810 families crossing shape boundary"),
]

with open(OUT / "fig5D_family_statistics.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["metric", "value", "label_for_figure"])
    for m, v, l in stats:
        w.writerow([m, v, l])

print(f"\n[F] fig5D_family_statistics.csv written")


# ══════════════════════════════════════════════════════════════
# G. README.md
# ══════════════════════════════════════════════════════════════

readme_text = f"""# Fig. 5D — LINMOL Family Trajectory Data

## Family Definition

- **Family label:** LINMOL
- **L3 ConnID hash:** {LINMOL_L3_HASH}
- **Metal:** Pd
- **CN:** 4
- **Members:** LINMOL, LINMOL01, LINMOL02 (3 entries)
- **Trajectory span:** {round(span, 2)} (largest among all full-CSD L3 families)

The LINMOL family is identified by its L3 ConnID (connectivity-level identity key).
All three refcodes share the same L3 hash but have distinct L0 (geometry-hash)
and L1 (ShapeID), demonstrating clear geometry-state diversity.

## Trajectory scalar

**Definition:** `trajectory_scalar = S_SP − S_Td`

This measures position along the SP→Td interconversion axis in CN=4 CShM space.
- Negative values → closer to ideal square-planar (SP)
- Positive values → closer to ideal tetrahedral (Td)
- Near zero → boundary / intermediate geometry

## Selected render states

Since the LINMOL family has 3 members (not 5), **all 3** are selected for rendering:

| State | Refcode  | S_SP  | S_Td  | Scalar | Shape    |
|-------|----------|-------|-------|--------|----------|
"""

for i, r in enumerate(selected_rows):
    readme_text += (
        f"| {i+1}     | {r['refcode']:8s} | {r['_s_sp']:5.2f} | {r['_s_td']:5.2f} "
        f"| {r['_trajectory_scalar']:+6.2f} | {r['_shape_label']:8s} |\n"
    )

readme_text += f"""
## XYZ file versions

For each state, three XYZ files are provided:

1. **`_full.xyz`** — Complete molecular structure (H atoms removed for clarity)
2. **`_first_sphere.xyz`** — Metal + donor atoms + donor-bonded backbone (no H).
   **Recommended for main-text rendering.**
3. **`_polyhedron.xyz`** — Metal + donor atoms only (5 atoms for CN=4).
   For minimal geometric illustration.

## Rendering recommendation

Use `first_sphere.xyz` for main-text Fig. 5D structure thumbnails.
The polyhedron version is a backup for extremely simplified renderings.

## Background sample

`fig5D_cn4_background_sample.csv` contains {n_bg:,} randomly sampled points
from the full-CSD CN=4 atlas (47,187 total) for plotting as gray background
scatter / density in the SP–Td plane.

## License

- XYZ files are CSD-derived and intended for **internal figure rendering only**.
- **Do NOT redistribute** raw coordinates, CIF, or XYZ files.
- No CIF or reversible CSD coordinate archives are exported.
- Only refcodes, hashed keys, aggregate statistics, CShM values, and field summaries
  appear in CSV/JSON outputs.

## Files

```
fig5d_family_trajectory/
├── README.md
├── fig5D_LINMOL_trajectory_points.csv
├── fig5D_LINMOL_trajectory_summary.json
├── fig5D_selected_render_states.csv
├── fig5D_cn4_background_sample.csv
├── fig5D_family_statistics.csv
├── fig5D_trajectory_preview.png
├── LINMOL_state01_full.xyz
├── LINMOL_state01_first_sphere.xyz
├── LINMOL_state01_polyhedron.xyz
├── LINMOL_state02_full.xyz
├── LINMOL_state02_first_sphere.xyz
├── LINMOL_state02_polyhedron.xyz
├── LINMOL_state03_full.xyz
├── LINMOL_state03_first_sphere.xyz
└── LINMOL_state03_polyhedron.xyz
```
"""

with open(OUT / "README.md", "w") as f:
    f.write(readme_text)

print(f"[G] README.md written")


# ══════════════════════════════════════════════════════════════
# H. Quick-look preview plot
# ══════════════════════════════════════════════════════════════

print("\n[H] Generating preview plot...")

ORANGE = "#E87D2F"
GRAY_DARK = "#333333"
GRAY_MED = "#888888"
GRAY_LIGHT = "#CCCCCC"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 10,
    "axes.edgecolor": GRAY_DARK,
    "axes.linewidth": 0.8,
})

fig, ax = plt.subplots(figsize=(7, 6))

# Background: CN4 atlas sample
bg_x = [math.log10(float(r["S_SP"]) + 0.1) for r in bg_sample]
bg_y = [math.log10(float(r["S_Td"]) + 0.1) for r in bg_sample]
ax.scatter(bg_x, bg_y, s=2, c=GRAY_LIGHT, alpha=0.3, rasterized=True, label="CN4 background")

# Diagonal line
mn = min(min(bg_x), min(bg_y)) - 0.1
mx = max(max(bg_x), max(bg_y)) + 0.1
ax.plot([mn, mx], [mn, mx], "--", color=GRAY_MED, lw=1, alpha=0.5)

# LINMOL trajectory
traj_x = [r["_log_sp"] for r in linmol_rows]
traj_y = [r["_log_td"] for r in linmol_rows]

# Sort by trajectory scalar for line plot
sorted_pts = sorted(zip(traj_x, traj_y, linmol_rows), key=lambda t: t[2]["_trajectory_scalar"])
sx = [p[0] for p in sorted_pts]
sy = [p[1] for p in sorted_pts]

# Draw trajectory line
ax.plot(sx, sy, "-", color=ORANGE, lw=2.5, zorder=5, label="LINMOL trajectory")

# Draw points with labels
for i, (x, y, r) in enumerate(sorted_pts):
    state_num = i + 1
    color = "red" if r["_boundary_flag"] else ORANGE
    ax.scatter([x], [y], s=120, c=color, edgecolors=GRAY_DARK, linewidths=1.2, zorder=6)
    ax.annotate(
        f"  {r['refcode']}\n  S_SP={r['_s_sp']:.1f}\n  S_Td={r['_s_td']:.1f}\n  Δ={r['_delta']:.1f}",
        (x, y),
        fontsize=7, color=GRAY_DARK,
        xytext=(12, -5 + i * 18), textcoords="offset points",
        bbox=dict(facecolor="white", edgecolor=GRAY_LIGHT, alpha=0.9, pad=2),
        arrowprops=dict(arrowstyle="-", color=GRAY_MED, lw=0.5),
        zorder=7,
    )

# Annotation box
ax.text(0.03, 0.97,
        f"LINMOL family (Pd, CN=4)\n"
        f"L3 hash: {LINMOL_L3_HASH}\n"
        f"Members: {n_members}\n"
        f"Span: {round(span, 2)}\n"
        f"Unique L0: {unique_L0}  L1: {unique_L1}\n"
        f"Boundary crossing: {'Yes' if crosses_boundary else 'No'}",
        transform=ax.transAxes, va="top", fontsize=8,
        bbox=dict(facecolor="white", edgecolor=GRAY_MED, alpha=0.9))

ax.set_xlabel("log₁₀(S_SP + 0.1)")
ax.set_ylabel("log₁₀(S_Td + 0.1)")
ax.set_title("Fig. 5D Preview — LINMOL Family Trajectory (CN=4 SP–Td)",
             fontsize=11, fontweight="bold", color=GRAY_DARK)
ax.legend(loc="lower right", fontsize=8)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

fig.tight_layout()
fig.savefig(OUT / "fig5D_trajectory_preview.png", dpi=200, facecolor="white")
plt.close(fig)
print(f"  → fig5D_trajectory_preview.png written")


# ══════════════════════════════════════════════════════════════
# Final verification
# ══════════════════════════════════════════════════════════════

print(f"\n{'='*60}")
print("Output verification")
print(f"{'='*60}")

required_files = [
    "fig5D_LINMOL_trajectory_points.csv",
    "fig5D_LINMOL_trajectory_summary.json",
    "fig5D_selected_render_states.csv",
    "fig5D_cn4_background_sample.csv",
    "fig5D_family_statistics.csv",
    "README.md",
]

# XYZ files
for i in range(len(selected_rows)):
    sn = f"{i+1:02d}"
    required_files.extend([
        f"LINMOL_state{sn}_full.xyz",
        f"LINMOL_state{sn}_first_sphere.xyz",
        f"LINMOL_state{sn}_polyhedron.xyz",
    ])

ok = 0
for fn in required_files:
    p = OUT / fn
    if p.exists():
        sz = os.path.getsize(p)
        print(f"  ✓ {fn:50s} {sz:>8,} B")
        ok += 1
    else:
        print(f"  ✗ {fn:50s} MISSING")

# Preview
preview = OUT / "fig5D_trajectory_preview.png"
if preview.exists():
    print(f"  ✓ {'fig5D_trajectory_preview.png':50s} {os.path.getsize(preview):>8,} B  (preview)")

print(f"\nRequired: {ok}/{len(required_files)} files present.")
print(f"{'='*60}")
