# CoordRep — JACS revision reproducibility package

This `jacs-revision` branch contains the public code and evidence for the
current manuscript, **CoordRep: A Canonical, Continuous, and Compositional
Representation for Coordination Chemistry**. The authoritative frozen package
is [`release/jacs-revision-20260823`](release/jacs-revision-20260823/).

## Current evidence lock

- **April 2025 CSD census:** 1,371,757 entries were examined with CCDC Python
  API 3.6.0. The predefined structural-audit target is 602,116 entries with a
  three-dimensional structure and at least one atom in the frozen in-domain
  set Sc–Zn, Y–Cd, La, and Hf–Hg.
- **Record emission and source fidelity:** 602,116/602,116 targets yielded a
  schema-valid audit outcome, and 602,116/602,116 passed an independent reread
  of the native CSD object.
- **Structural coverage:** every in-domain site is structural in
  601,402/602,116 entries (99.8814%); 1,984,062/1,986,198 in-domain sites are
  structurally recorded (99.8925%).
- **Typed relations:** 147,348 entries contain an in-domain shared donor group,
  102,241 contain a confirmed haptic/π site, and 94,452 contain a resolved
  nonzero lattice-translation edge. These flags are nonexclusive.
- **Periodic collection:** 15,905/15,906 public CSD MOF Collection entries were
  processed. The 10,948 state-bearing entries yielded 172,332 local states;
  74,354 (43.15%) contain a nonzero translation-labelled edge.

The value 602,116 is not an all-metal CSD count. Entries containing only
Ce–Lu, actinides, or main-group metals are outside this denominator; native
relations to such centers are retained when they occur in mixed-metal target
entries. Raw licensed CSD coordinates and row-level structure ledgers are not
redistributed.

## Canonicalization version boundary

CoordRep **1.1.2rc3** adds ligand-local attachment-set keys, donor-attachment
orbits, exact residual-orbit enumeration, and fail-closed handling when the
metadata required to establish exchangeability is unavailable. Figure 2 and
Supplementary Figure S2 are the locked rc3 canonicalization evidence.

This change can alter molecular strings, token sequences, and L0–L3 hashes
relative to rc2. The Figure 4 property benchmarks and the molecular-family
statistics shown in Figure 6 remain the explicitly labelled frozen rc2 evidence
used in the manuscript; this release does not claim that those downstream
cohorts were rerun under rc3. The release-wide CSD audit counts are
source-derived typed-record census quantities and do not depend on the final
ordering of an rc3 molecular string. Periodic MID/SID results use their
separate frozen periodic-v3 identity construction.

## Start here

```bash
git clone https://github.com/1701NJU/CoordRep.git
cd CoordRep
git checkout jacs-revision
python -m pip install -e "release/jacs-revision-20260823[dev]"
python libcoordrep/scripts/check_jacs_revision_release.py
```

The full database rerun requires a licensed CCDC installation and local April
2025 CSD access. See
[`CSD_REDISTRIBUTION_NOTICE.md`](CSD_REDISTRIBUTION_NOTICE.md).

## Repository map

| Item | Location |
|---|---|
| Current frozen release | `release/jacs-revision-20260823/` |
| Full-CSD aggregate evidence | `release/jacs-revision-20260823/audits/full_csd/` |
| Current Figure 1–6 artwork and source tables | `release/jacs-revision-20260823/figures/` |
| Supplementary Figures S1–S2 | `release/jacs-revision-20260823/figures/Supplementary/` |
| rc3 canonicalization evidence | `release/jacs-revision-20260823/canonicalization/` |
| CoordRep 1.1.2rc3 source | `release/jacs-revision-20260823/{coordrep,coordrep_tools,brain}/` |
| Supporting Information | `release/jacs-revision-20260823/supporting_information/` |
| Superseded frozen releases | `release/jacs-revision-20260819/`, `release/jacs-revision-20260806/` |

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
