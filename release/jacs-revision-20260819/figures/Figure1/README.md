# Figure 1C CoSyMLib polyhedra

## Recommended files

- `Fig1C_Cosymlib_ATIZEK_triptych_20260805.svg`: compact Figure 1C draft containing the nearest ideal Oh reference, the PBE ATIZEK RhN3O3 first sphere, and the nearest ideal TPr reference.
- `Fig1C_ATIZEK_nearest_OC-6_icon_20260805.svg`: transparent Oh icon aligned and scaled to ATIZEK by CoSyMLib.
- `Fig1C_ATIZEK_nearest_TPR-6_icon_20260805.svg`: transparent TPr icon aligned and scaled to ATIZEK by CoSyMLib.
- `Fig1C_Cosymlib_ideal_reference_pair_20260805.svg`: the two fitted ideal references without the central ATIZEK structure.
- `Fig1C_OC-6_cosymlib_canonical_icon_20260805.svg` and `Fig1C_TPR-6_cosymlib_canonical_icon_20260805.svg`: canonical, unaligned CoSyMLib references for general reuse.

All SVG files are vector-only and contain no embedded raster image. Matching PDF and 600 dpi PNG files are provided for checking and Word placement.

## Scientific provenance

The official CoSyMLib CN = 6 labels were used:

- `OC-6`: octahedron, Oh symmetry;
- `TPR-6`: trigonal prism, D3h symmetry.

The source geometry is `xyz/ATIZEK.xyz` in the locally frozen official tmQMg archive. Its SHA-256 is:

`e3d68ed616618aa43c5f248fc4c969f7d9d79f6633b39ad60b79265731a07825`

The Rh center and six audited donors were selected in the exact CoordRep order `[48, 44, 42, 46, 39, 40, 41]` using zero-based indices. CoSyMLib reproduced:

- `S(OC-6) = 1.204080932467511`
- `S(TPR-6) = 12.311628666783502`

The ideal polyhedra used in the ATIZEK artwork were generated with `Shape.structure(label, central_atom=1)`. They are therefore exact ideal OC-6/TPR-6 objects after the optimum CoSyMLib translation, rotation, scale, and donor permutation for the real ATIZEK first coordination sphere. They were not drawn by eye.

Coordinates and complete metadata are recorded in:

- `Fig1C_Cosymlib_coordinates_20260805.csv`
- `Fig1C_Cosymlib_manifest_20260805.json`
- the five `.xyz` files in this directory.

## Reproduction

The tested Windows environment is Python 3.10 with CoSyMLib 0.11.1, NumPy 1.26.4, SciPy 1.10.1, wfnsympy 0.3.5, and symgroupy 0.5.9. CoSyMLib 0.12.0/0.12.1 was not used because the Windows wheels downloaded during this run contained a CPython 3.14 extension despite their CPython 3.12 filename tags.

From the JACS project root:

```powershell
& "tmp\cosymlib_py310_env\python.exe" `
  "revision\figure1c_cosymlib_polyhedra_20260805\generate_fig1c_cosymlib_polyhedra.py"
```

The environment can be recreated with:

```powershell
conda create -y -p "tmp\cosymlib_py310_env" python=3.10 pip
& "tmp\cosymlib_py310_env\python.exe" -m pip install `
  "numpy==1.26.4" "scipy==1.10.1" "matplotlib==3.8.4" `
  "PyYAML==6.0.3" "wfnsympy==0.3.5" "symgroupy==0.5.9" `
  "huckelpy==0.3" "pointgroup==0.4.4" "cosymlib==0.11.1"
```

Software: https://github.com/GrupEstructuraElectronicaSimetria/cosymlib
CoSyMLib archive DOI: https://doi.org/10.5281/zenodo.4925766

## Suggested caption text

> **(C) Continuous geometry.** The PBE-optimized RhN3O3 first coordination sphere of the ATIZEK snapshot shown in panels A and B is evaluated against the nearest ideal octahedral (OC-6, Oh) and trigonal-prismatic (TPR-6, D3h) reference structures generated with CoSyMLib. The resulting shape vector, [S(Oh), S(TPr)] = [1.20, 12.31], identifies an Oh-like coordination environment while retaining its measured deviation from both reference geometries. Ligand backbones are omitted for clarity.
