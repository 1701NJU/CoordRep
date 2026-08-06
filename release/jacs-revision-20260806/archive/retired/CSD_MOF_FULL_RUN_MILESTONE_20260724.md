# CoordRep periodic MOF application: frozen full-run milestone

Date: 2026-07-24
Protocol: `CSD_MOF_PERIODIC_LOCAL_SITE_v2`
Protocol SHA256: `A522BC5F9DC7D45D804D01BC733D9873F1E97BEE23BCC2777B53CDF19728C1E3`
Run signature SHA256: `978A7622BC61707F2627DCB570847D14925A928F39384AD57261D1BAF9A7EBA3`

## Milestone status

The frozen pipeline completed an unlimited run over all 15,906 entries in the public CSD MOF Collection. The terminal publication-integrity gate passed all 35 checks. The run produced 15,906 unique entry records and 224,342 unique in-domain-metal site records in 2,947.2 s. All 12 output hashes in `SHA256SUMS.csv` were independently rechecked after completion and matched.

This application should be presented as a periodic-aware first-sphere adapter into a CoordRep-compatible molecular local-state layer. It is not a claim that the local-state token encodes complete MOF framework topology.

## Full-collection deployment results

| Stage | Result |
|---|---:|
| CIF parsing | 15,906/15,906 (100.000%) |
| Complete pipeline processing | 15,901/15,906 (99.969%) |
| Processing failures | 5/15,906 (0.031%; all ambiguous torus clustering) |
| Processed entries containing a predefined target metal | 12,802/15,901 (80.511%) |
| Target-metal entries with at least one CN 2–6 site | 11,641/12,802 (90.931%) |
| Target-metal entries with at least one emitted local state | 10,957/12,802 (85.588%) |
| Target-metal sites | 224,342 |
| Chemical-domain and vector-witness-passing emitted sites | 172,783 |
| Primary internal-consistency subset | 170,021/172,783 emitted sites (98.401%) |
| Local-state token parse/serialize round-trip | 172,783/172,783 (100.000%) |

The emitted-site atlas is dominated by CN 6 (66,649), CN 4 (57,319), and CN 5 (38,110). The most frequent idealized shapes are Oh (64,093), Td (49,586), and SPY (25,048). The most represented metals are Zn (50,606), Cu (42,294), and Co (20,225).

## Evidence that the periodic adapter performs substantive work

In the supplied P1 cell representation, 74,483/172,783 emitted sites (43.108%) contained at least one explicit coordination edge with a nonzero lattice translation. The equal-family mean of the within-entry fraction was 53.280% (bootstrap stability interval 52.670–53.907%). This metric is cell/origin dependent and must not be described as an intrinsic topological invariant, but it demonstrates that periodic translation recovery is used extensively rather than being a cosmetic code path.

The narrower multi-image-union QC flag occurred for 2,760/172,783 emitted sites (1.597%). It should not replace the nonzero-translation result as the main evidence for periodic processing.

## Family-disjoint confirmatory results

The confirmatory cohort contains 12,607 entries whose refcode families do not overlap the 3,299-entry/2,396-family development set or any family encountered in prior sampled audits.

| Confirmatory criterion | Observed | Frozen threshold | Outcome |
|---|---:|---:|---|
| CIF parse success | 100.000% | ≥99% | Met |
| Pipeline processing success | 99.968% | ≥99% | Met |
| Chemical-domain vector witness, equal-family mean | 100.000% | ≥99% | Met |
| Local-state round-trip, equal-family mean | 100.000% | ≥99% | Met |
| Complete-explicit-image availability, equal-family mean | 97.9907% | ≥98% | **Narrowly missed** |
| Best-shape agreement when a complete image is available | 100.000% | ≥99% | Met |
| Absolute CShM difference, p95 | 6.53 × 10⁻¹⁴ | ≤0.001 | Met |

The complete-image availability criterion missed its development-informed threshold by 0.0093 percentage points. It must not be silently rounded to a pass or used to change the frozen threshold. Report instead:

> Complete explicit-image comparison was available for 125,091/127,285 emitted confirmatory sites (98.276% by raw site count; 97.991% equal-family mean, 95% bootstrap stability interval 97.777–98.190%). Among comparable sites, best-shape agreement was 100%, and the p95 absolute CShM difference was 6.53 × 10⁻¹⁴.

This narrow miss does not invalidate the main atlas because the atlas was prospectively defined over all chemical-domain/vector-witness-passing emitted sites. It only prevents a claim of at least 98% family-equal complete-image availability.

## Release/version-shift audit

Of the 12,607 confirmatory entries, 975 have no exact refcode match in the local CSD v6.00 installation. The stricter comparison retains 944 entries whose entire refcode family is also absent from the exact-refcode-present collection cohort.

For this stricter family-unmatched cohort:

- Pipeline processing differed from the present cohort by −0.085 percentage points (family-bootstrap stability interval −0.329 to +0.055), providing no evidence of a material processing collapse.
- The probability of containing a predefined target metal was 3.951 percentage points higher (interval +1.238 to +6.514), demonstrating a composition shift.
- Conditional on successful processing and an in-domain metal, representation-emission coverage was 6.324 percentage points lower (interval −9.408 to −3.359).
- The analogous primary internal-consistency coverage difference was −6.434 percentage points (interval −9.560 to −3.460).

These are coverage/domain-shift results, not accuracy estimates and not a fully prospective time split. The correct label is “exact-refcode-unmatched release/version-shift cohort.”

## Independent distance sensitivity

The fixed validation sample contains 400 unique, development-family-disjoint families: 200 exact-refcode-present and 200 exact-refcode-absent. Every selected entry parsed and processed, every emitted validation site completed the bounded lattice search, no search was nonconvergent, and the maximum translation radius used was 2.

The distance route agreed exactly with explicit connectivity for 88.642% of present-cohort emitted sites and 91.896% of absent-cohort emitted sites. Among distance-shape-comparable sites, best-shape agreement was 94.150% and 95.401%, respectively. These imperfect but high agreements support using distance adjacency as sensitivity/QC, not as chemical ground truth or the primary inclusion rule.

## Manuscript-level interpretation

This full application is substantially stronger than the 4.88% overall HOMO–LUMO-gap MAE improvement because it establishes an orthogonal contribution:

1. deployment from molecular coordination complexes to 15,906 periodic MOFs;
2. explicit recovery of periodic translation and translated multiedges;
3. a 172,783-site, coordinate-free local-state atlas;
4. family-disjoint validation and a strict release/version-shift coverage audit; and
5. transparent failure and domain-boundary reporting.

The GNN comparison should remain a bounded complementarity result. It should say that CoordRep supplies interpretable, serializable local-state information that can improve specific paired material-gap comparisons, while GNNs/SchNet remain stronger on some tasks. It should not claim universal predictive superiority.

## Recommended next milestone

1. Replace the old Rosetta-centered Figure 5 with a vector four-panel MOF application figure:
   - full-collection coverage waterfall;
   - metal × geometry atlas;
   - periodic-translation and internal-consistency evidence;
   - present versus strict family-unmatched coverage-difference forest plot.
2. Write the exact replacement paragraphs and caption against `新建文件夹/manuscript - 20260701.docx`, without restructuring the rest of the paper.
3. Add a response-to-reviewer section that separates the bounded GNN comparison from the new periodic-MOF deployment.
4. Preserve the complete-image threshold miss explicitly in the Supporting Information and rebuttal.

## Evidence files

- `results/csd_mof_periodic_application_v2/SUMMARY_PUBLIC_AGGREGATE.json`
- `results/csd_mof_periodic_application_v2/RUN_MANIFEST_PUBLIC.json`
- `results/csd_mof_periodic_application_v2/SHA256SUMS.csv`
- `CSD_MOF_PERIODIC_LOCAL_SITE_PROTOCOL_v2.json`
