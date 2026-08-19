# CoordRep — JACS revision reproducibility package

This `jacs-revision` branch contains the public code and evidence for the
current CoordRep manuscript revision. The authoritative package is
[`release/jacs-revision-20260819`](release/jacs-revision-20260819/). It
supersedes the 2026-08-06 package and synchronizes the corrected molecular
figures with the release-wide April 2025 CSD audit.

## Current evidence lock

- **Full-CSD audit:** 1,371,757 entries were examined. The operational target
  comprises 602,116 entries with a three-dimensional structure and at least
  one in-domain d-block center.
- **Record emission and source fidelity:** 602,116/602,116 targets yielded
  schema-valid audit records, and 602,116/602,116 passed an independent reread
  of the native CSD object.
- **Structural coverage:** every in-domain site is structural in
  601,402/602,116 entries (99.8814%); 1,984,062/1,986,198 individual sites are
  structurally recorded (99.8925%).
- **Relations retained:** 147,348 entries contain an in-domain shared donor
  group, 102,241 contain a confirmed haptic/π site, and 94,452 contain a
  resolved nonzero lattice-translation edge. These are nonexclusive counts.
- **Periodic extension:** 15,905/15,906 public CSD MOF Collection entries were
  processed. The 10,948 state-bearing entries yielded 172,332 local states;
  74,354 states (43.15%) contain a nonzero translation-labelled edge.

The 101,878-record strict mononuclear η¹ core and 84,453-record
connectivity-supported tier are downstream analysis subsets, not estimates of
overall representational coverage. Likewise, full-CSD emission and source
fidelity do not imply exact chemical-species identity: 305,845 audit records
are exact label-free, while 296,271 conservatively retain a source-order tie.

## Start here

```bash
git clone https://github.com/1701NJU/CoordRep.git
cd CoordRep
git checkout jacs-revision
python -m pip install -e "release/jacs-revision-20260819[dev]"
python libcoordrep/scripts/check_jacs_revision_release.py
```

The full CSD audit requires a licensed CCDC installation and local CSD access.
Raw CSD coordinates, structure files, and licensed row-level exports are not
redistributed; see
[`CSD_REDISTRIBUTION_NOTICE.md`](CSD_REDISTRIBUTION_NOTICE.md).

## Repository map

| Item | Location |
|---|---|
| Current release and integrity lock | `release/jacs-revision-20260819/` |
| Full-CSD aggregate evidence | `release/jacs-revision-20260819/audits/full_csd/` |
| Current Figure 1–6 artifacts | `release/jacs-revision-20260819/figures/` |
| CoordRep 1.1.2rc2 source | `release/jacs-revision-20260819/{coordrep,coordrep_tools,brain}` |
| Full-CSD audit implementation | `release/jacs-revision-20260819/{coordrep/audit,scripts}/` |
| E(3)-equivariant baselines | `release/jacs-revision-20260819/models/E3_baselines/` |
| Periodic-v3 and ML protocols | `release/jacs-revision-20260819/protocols/` |
| Superseded 2026-08-06 entry point | `release/jacs-revision-20260806/` |
| Historical result tree | `libcoordrep/revision_results/` |

## License and citation

Code is released under the MIT License. CSD-derived outputs remain subject to
the restrictions in `CSD_REDISTRIBUTION_NOTICE.md`.

```bibtex
@software{coordrep2026,
  title  = {CoordRep: A Canonical, Continuous, Compositional, and Multi-Resolution Representation Layer for Machine Learning in Coordination Chemistry},
  author = {Li, Cheng-Hui and Luo, Wen-Lin},
  year   = {2026},
  url    = {https://github.com/1701NJU/CoordRep}
}
```
