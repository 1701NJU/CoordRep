# CoordRep — JACS revision reproducibility package

This `jacs-revision` branch corresponds to the current CoordRep manuscript
revision. The authoritative publication package is
[`release/jacs-revision-20260806`](release/jacs-revision-20260806/), which
contains the `CoordRep 1.1.2rc2` source, current Figure 1–6 artwork and public
source tables, the modern E(3)-equivariant controls, frozen protocols, and
release-level integrity metadata.

## What is current

- **3C record language:** canonicalization, continuous CShM coordinates, and
  compositional metal/ligand/donor-relation fields, including finite
  multinuclear and haptic/π coordination-site records.
- **Molecular CSD survey:** 1,371,757 total entries; 602,116 3D transition-
  metal entries; 101,878 conservative mononuclear records; 84,453
  ligand-topology records; 52,760 linked-multimetal candidates; 38,203
  haptic/π candidates.
- **Property benchmark:** frozen 48,057-record graph-grouped cohort with GINE,
  CoordRep relation/shape ablations, SchNet, and ViSNet controls.
- **Periodic extension:** all 15,906 entries of the public CSD MOF Collection;
  15,905 processed entries, 172,332 local states from 10,948 entries, 43.15%
  of local states carrying a nonzero translation-labelled edge, 36,552 metric
  IDs, and 60,139 stereo IDs.

The periodic result is a canonical local first-sphere/state audit; it is not a
claim of complete MOF-framework identity or universal predictive superiority
over graph or Cartesian geometric neural networks.

## Start here

```bash
git clone https://github.com/1701NJU/CoordRep.git
cd CoordRep
git checkout jacs-revision
python -m pip install -e "release/jacs-revision-20260806[dev]"
python libcoordrep/scripts/check_jacs_revision_release.py
```

For CSD-dependent reproduction, install the CCDC Python API and use the local
licensed database. Raw CSD coordinates, CIF/MOL/MOL2 files, and commercial
database exports are intentionally excluded; see
[`CSD_REDISTRIBUTION_NOTICE.md`](CSD_REDISTRIBUTION_NOTICE.md).

## Repository map

| Item | Location |
|---|---|
| Current release README and manifest | `release/jacs-revision-20260806/` |
| Figure-to-manuscript map | `release/jacs-revision-20260806/figures/FIGURE_MAP.md` |
| CoordRep 1.1.2rc2 source | `release/jacs-revision-20260806/{coordrep,coordrep_tools,brain}` |
| E(3)-equivariant baselines | `release/jacs-revision-20260806/models/E3_baselines/` |
| Periodic-v3 and ML protocols | `release/jacs-revision-20260806/protocols/` |
| Historical/retired result tree | `libcoordrep/revision_results/` (not current evidence) |

## License and citation

Code is released under the MIT License. CSD-derived outputs are subject to
the restrictions in `CSD_REDISTRIBUTION_NOTICE.md`.

```bibtex
@software{coordrep2026,
  title  = {CoordRep: A Canonical, Continuous, Compositional, and Multi-Resolution Representation Layer for Machine Learning in Coordination Chemistry},
  author = {Li, Cheng-Hui and Luo, Wen-Lin},
  year   = {2026},
  url    = {https://github.com/1701NJU/CoordRep}
}
```
