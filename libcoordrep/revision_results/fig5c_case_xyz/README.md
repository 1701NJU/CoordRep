# Fig. 5C — CSD Casebook XYZ Files

## Purpose
Local rendering assets for the three CSD casebook panels in Fig. 5C.

## ⚠️ License
All XYZ files in this directory are **CSD-derived** and are for
**internal figure rendering only**. Do NOT redistribute raw coordinates
or upload to GitHub / SI data packages.

What CAN be published:
- Rendered images (PNG/SVG)
- Aggregated statistics (CShM values, validator status, ranker scores)
- Case summaries and CoordRep field excerpts
- Refcode identifiers

## Cases

| Case | Refcode | Metal/CN | Purpose | Recommended XYZ |
|------|---------|----------|---------|-----------------|
| Boundary | AFOSIA | Cr/CN6 | Boundary-aware shape identity avoids brittle discrete shape … | `AFOSIA_first_sphere.xyz or AFOSIA_polyhedron.xyz` |
| Repair | AGOTIA | Re/CN7 | Explicit grammar enables deterministic syntax repair and fie… | `AGOTIA_first_sphere.xyz` |
| Stereo | CIJWUO | Cu/CN5 | CoordRep fields enable semantic consistency checking beyond … | `CIJWUO_first_sphere.xyz or CIJWUO_stereo_view.xyz` |

## Directory structure

```
fig5c_case_xyz/
├── README.md
├── case_xyz_index.csv
├── AFOSIA/
│   ├── AFOSIA_raw.xyz
│   ├── AFOSIA_clean_for_render.xyz
│   ├── AFOSIA_first_sphere.xyz
│   ├── AFOSIA_polyhedron.xyz
│   ├── AFOSIA_render_metadata.json
│   └── AFOSIA_figure_notes.md
├── AGOTIA/
│   ├── AGOTIA_raw.xyz
│   ├── AGOTIA_clean_for_render.xyz
│   ├── AGOTIA_first_sphere.xyz
│   ├── AGOTIA_repair_excerpt.txt
│   ├── AGOTIA_render_metadata.json
│   └── AGOTIA_figure_notes.md
└── CIJWUO/
    ├── CIJWUO_raw.xyz
    ├── CIJWUO_clean_for_render.xyz
    ├── CIJWUO_first_sphere.xyz
    ├── CIJWUO_stereo_view.xyz
    ├── CIJWUO_stereo_tokens.txt
    ├── CIJWUO_render_metadata.json
    └── CIJWUO_figure_notes.md
```

## Rendering tips

1. **Software**: ChimeraX, PyMOL, or Olex2
2. **Background**: white, no axes, no labels (except optional metal/CN)
3. **Hydrogens**: hidden (use `_clean_for_render.xyz` or `_first_sphere.xyz`)
4. **Metal**: highlighted in distinct color
5. **Resolution**: ≥2000×2000 px for publication
6. **Scale**: consistent across all three panels
7. For crowded structures, prefer `_first_sphere.xyz`
