# CoordRep — JACS revision reproducibility package

This `jacs-revision` branch contains the public code and evidence for the
current manuscript, **CoordRep: A Canonical, Continuous, and Compositional
Representation for Coordination Chemistry**. The authoritative package is
[`release/jacs-revision-20260826`](release/jacs-revision-20260826/).

## Current evidence lock

- **April 2025 CSD census:** 1,371,757 entries were examined with CCDC Python
  API 3.6.0 under a frozen 96-symbol, CCDC-derived metal policy. In total,
  783,263 entries are metal-containing; 746,511 deposited structures meet the
  three-dimensional structural-audit criterion, while 36,752 remain in the
  census without a 3D structural record.
- **Record emission and source fidelity:** all 746,511 targets yielded a
  schema-valid audit record, and 746,511/746,511 passed an independent reread
  of the native CSD object.
- **Structural coverage:** at least one structural metal-site record was
  retained in 741,402/746,511 entries; every in-scope metal site is structural
  in 733,004/746,511 entries. At site level, 2,574,081/2,608,448 metal sites
  are structural and 34,367 are explicitly audit-only.
- **Target classes:** 298,932 nonpolymeric single-metal, 295,164
  nonpolymeric multimetal, and 152,415 polymeric entries form a mutually
  exclusive partition of the 746,511-entry structural-audit target.
- **Typed relations:** 213,595 entries contain a shared donor group, 113,726
  contain a confirmed haptic/π site, and 129,092 contain a resolved nonzero
  lattice-translation edge. These flags are nonexclusive; 63 entries retain
  an ambiguous collective-π candidate without forcing hapticity.
- **Periodic collection:** 15,905/15,906 public CSD MOF Collection entries were
  processed. The 10,948 state-bearing entries yielded 172,332 local states;
  74,354 (43.15%) contain a nonzero translation-labelled edge.

These are source-faithful structural-transcription results, not a claim that
every site has a unique exact CoordRep identity, a supported CShM vector, or a
distance-completed first coordination sphere. Donor, bridge, haptic/π, and
periodic relations in the release-wide audit are derived from the CSD-native
molecular bond graph; no distance-derived contact is added. Raw licensed CSD
coordinates and row-level ledgers are not redistributed.

## Canonicalization version boundary

CoordRep **1.1.2rc3** adds ligand-local attachment-set keys, donor-attachment
orbits, exact residual-orbit enumeration, and fail-closed handling when the
metadata required to establish exchangeability is unavailable. Figure 2 and
Supplementary Figure S2 are the locked rc3 canonicalization evidence.

This change can alter molecular strings, token sequences, and L0–L3 hashes
relative to rc2. Figure 4 and the molecular-family statistics shown in Figure
6 remain explicitly labelled frozen rc2 evidence; they are not claimed as rc3
reruns. Release-wide CSD counts are typed-record census quantities and do not
depend on final rc3 string ordering. Periodic MID/SID results use the separate
frozen periodic-v3 identity construction.

## Start here

```bash
git clone https://github.com/1701NJU/CoordRep.git
cd CoordRep
git checkout jacs-revision
python -m pip install -e "release/jacs-revision-20260826[dev]"
python libcoordrep/scripts/check_jacs_revision_release.py
```

The full database rerun requires a licensed CCDC installation and local April
2025 CSD access. See
[`CSD_REDISTRIBUTION_NOTICE.md`](CSD_REDISTRIBUTION_NOTICE.md).

## Repository map

| Item | Location |
|---|---|
| Current frozen release | `release/jacs-revision-20260826/` |
| All-metal CSD aggregate evidence | `release/jacs-revision-20260826/audits/full_csd/` |
| Current Figure 1–6 artwork and source tables | `release/jacs-revision-20260826/figures/` |
| Supplementary Figures S1–S2 | `release/jacs-revision-20260826/figures/Supplementary/` |
| rc3 canonicalization evidence | `release/jacs-revision-20260826/canonicalization/` |
| CoordRep 1.1.2rc3 source | `release/jacs-revision-20260826/{coordrep,coordrep_tools,brain}/` |
| Current audit-methods/Table S9 excerpt | `release/jacs-revision-20260826/supporting_information/` |
| Superseded frozen releases | `release/jacs-revision-20260823/`, `release/jacs-revision-20260819/`, `release/jacs-revision-20260806/` |

## License and citation

Code is released under the MIT License. CSD-derived outputs remain subject to
the restrictions in `CSD_REDISTRIBUTION_NOTICE.md`.

```bibtex
@software{coordrep2026,
  title  = {CoordRep: A Canonical, Continuous, and Compositional Representation for Coordination Chemistry},
  author = {Luo, Wen-Lin and Li, Cheng-Hui},
  year   = {2026},
  url    = {https://github.com/1701NJU/CoordRep}
}
```
