# CoordRep — JACS revision release (2026-08-06)

This directory is the publication-facing release for the `jacs-revision`
branch. It couples the audited `CoordRep 1.1.2rc2` source with the public
source tables, vector figures, model baselines, and protocols used for the
current manuscript and Supporting Information.

The release is deliberately separated from the historical material retained
under `libcoordrep/revision_results/`. Files in the historical tree may refer
to earlier serializers, the retired Figure 5G/J package, or the superseded
15,901-entry periodic pipeline; they are not evidence for the current claims.

## Current evidence package

| Scope | Current release evidence |
|---|---|
| CoordRep implementation | `coordrep/`, `coordrep_tools/`, `brain/`, `scripts/`, `tests/` |
| Figure 1 | `figures/Figure1/` — 3C record anatomy, CoSyMLib shape references, and molecular/periodic scope assets |
| Figure 2 | `figures/Figure2/` — canonicalization stress test and L0–L3 identity source files |
| Figure 3 | `figures/Figure3/` — public CShM atlas tables and manuscript preview; the final composite is supplied in the submission artwork |
| Figure 4 | `figures/Figure4/` — relation/shape field-learning comparison, modern E(3)-equivariant controls, vector artwork, and source tables |
| Figure 5 | `figures/Figure5/` — current CSD molecular-scope survey and public source table |
| Figure 6 | `figures/Figure6/` — matched periodic-v3 vector artwork, caption, and public collection census |
| Baselines | `models/E3_baselines/` — SchNet, ViSNet, and PaiNN-compatible runners/model code plus compact OOF audit outputs |
| Frozen protocols | `protocols/` — periodic canonical-v3, corrected CShM hybrid, and ViSNet OOF protocols |

## Numbers locked for the current manuscript

The current periodic-v3 audit covers all 15,906 records in the public CSD MOF
Collection: 15,905 processed entries, 10,948 state-bearing entries, 172,332
local metal-site states, 74,354 states with at least one nonzero translation
edge (43.15%), 36,552 metric local-state IDs, and 60,139 stereo-distinct IDs.
The family-disjoint audit contains 400 entries and 5,920 local sites. These
values are reproduced in `figures/Figure6/Figure6_source_data_20260806.csv`.

The molecular CSD survey contains 1,371,757 entries, 602,116 3D transition-
metal entries, 101,878 conservative mononuclear records, 84,453 ligand-
topology records, 52,760 linked-multimetal candidates, and 38,203 haptic/π
candidates. The public source table is
`figures/Figure5/Figure5_CoordRep_current_CSD_source_data_20260802.csv`.

The graph-degenerate property benchmark uses the frozen 48,057-record cohort,
five graph-grouped folds, and seeds 11/22/33. The current Figure 4 source table
reports the relation and shape ablations together with SchNet and ViSNet
controls; no claim of universal superiority over geometric neural networks is
made.

## Installation and tests

```bash
cd release/jacs-revision-20260806
python -m pip install -e ".[dev]"
pytest -q
```

The optional `ml` and `csd` extras are intentionally not installed by the
default test command. CSD scripts require a licensed CCDC/CSD Python API and
local access to the user's database; raw CSD coordinates and structure files
are not redistributed here.

Run the release integrity check from the repository root:

```bash
python libcoordrep/scripts/check_jacs_revision_release.py
```

## Data and licensing boundary

The repository contains code, aggregate statistics, figure-source tables,
refcodes/hashes, and audit metadata only. It does not contain CSD CIF/MOL/MOL2
files, coordinate dumps, or commercial database exports. See the root
`CSD_REDISTRIBUTION_NOTICE.md` before using any CSD-dependent workflow.

## Citation

```bibtex
@software{coordrep2026,
  title  = {CoordRep: A Canonical, Continuous, Compositional, and Multi-Resolution Representation Layer for Machine Learning in Coordination Chemistry},
  author = {Li, Cheng-Hui and Luo, Wen-Lin},
  year   = {2026},
  url    = {https://github.com/1701NJU/CoordRep}
}
```
