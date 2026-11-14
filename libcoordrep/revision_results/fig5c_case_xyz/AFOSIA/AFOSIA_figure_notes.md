# AFOSIA — Boundary Geometry (Fig. 5C, Case B)

## What this case demonstrates

AFOSIA is a Cr(III) CN=6 complex where CShM(TPr)=4.49 and CShM(Oh)=4.48.
The gap Δ=0.01 is far below the boundary threshold (1.0), meaning any discrete
shape label assignment ("octahedral" vs "trigonal prismatic") is essentially
arbitrary — the geometry sits exactly on the boundary.

CoordRep's continuous shape annotation preserves **both** CShM values and
flags the boundary explicitly via `TPr/Oh_boundary` in the L1 identity token,
rather than committing to a single brittle label.

## Rendering guidance

- Use `AFOSIA_polyhedron.xyz` (metal + 6 donors) for a clean CN6 polyhedron.
- Use `AFOSIA_first_sphere.xyz` for a slightly richer ball-and-stick view.
- Highlight Cr in a distinct color; label donor atoms (O, Cl, I).
- The full molecule (`_clean_for_render.xyz`) has 3 THF rings + Cl + I — may be too
  crowded for a small panel.

## Key numbers for annotation

| Metric | Value |
|--------|-------|
| CShM(TPr) | 4.49 |
| CShM(Oh)  | 4.48 |
| Δ         | 0.01 |
| boundary_thresh | 1.0 |
| L1 token  | `TPr/Oh_boundary.dist.D0` |
