# Fig. 5D — Family-to-Family Trajectory Visualization Data

## Overview

This directory contains figure-ready data for Fig. 5D, showing how
L3 connectivity families traverse the CN=4 SP–Td continuous shape measure space.

## Primary showcase family

- **LINMOL** (Pd, CN=4): 3 members, span = 10.10, boundary-crossing
- L3 ConnID hash: 463888590564

## Files

### Trajectory data
- **`fig5D_top_family_trajectories_overlay.csv`** — 41 point-rows across
  15 top families; each row is one observed geometry state within a family.
- **`fig5D_family_centroids.csv`** — One row per family with centroid, span, region.
- **`fig5D_family_edges.csv`** — 5 family-to-family edges (same metal,
  ligand similarity ≥ 0.75 or donor-set distance ≤ 1).

### LINMOL XYZ thumbnails
- `LINMOL_state01_{full,first_sphere,polyhedron}.xyz` — trajectory start
- `LINMOL_state02_{full,first_sphere,polyhedron}.xyz` — intermediate
- `LINMOL_state03_{full,first_sphere,polyhedron}.xyz` — trajectory end
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

CSD-derived coordinates for internal figure rendering only; do not redistribute raw XYZ.
No CIF, no reversible coordinate archives.
XYZ files for internal rendering only.
