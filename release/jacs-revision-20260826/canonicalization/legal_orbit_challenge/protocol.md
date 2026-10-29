# Pre-registered protocol: real legal-orbit canonicalization challenge

Date locked: 2026-08-22

Status: completed against the frozen six-file attachment-aware source now
integrated as CoordRep 1.1.2rc3. This protocol was written before the all-6,427
scan and K = 100 challenge results were inspected. It is not a CSD-prevalence
experiment.

## Source population

The source population is exactly the 6,427 graph-disjoint eligible molecular
graphs established by `canonical_stress_20260822`: one transition-metal center,
CN 4--6, complete successful rc2 pipeline status, no rc2 validation issues,
SMILES ligand payloads only, and official tmQMg/PBE XYZ payload hashes verified
against the frozen ZIP. If more than one eligible row carries the same graph
hash, one representative is chosen by the smallest SHA-256 score of
`20260822|unique_graph_representative|graph_sha256|mol_id`.

## Source-index-free challenge screen

Each representative source snapshot is encoded once by the frozen candidate.
The screen uses only ligand-payload fields, ligand-local coloured-graph
attachment keys, donor-attachment orbit keys, and quantized local donor
geometry already present in the encoded record. Global/source atom indices are
not used.

A record is challenge-eligible if the exact canonicalizer has more than one
legal candidate after initial invariant refinement, because at least one of
the following is present:

1. two or more ligand instances have the same payload, joint attachment-set
   key, and donor-attachment-key multiset; or
2. a same-element donor-rank orbit has size greater than one after topology
   and quantized local-geometry refinement.

The exact candidate count is the product of the factorials of all legal
ligand-instance group sizes and donor-rank orbit sizes. Repeated payloads that
are not attachment-equivalent are reported but do not by themselves make a
record challenge-eligible. Missing attachment identity in a repeated-payload
group is recorded as a fail-closed case.

## Selection if the challenge pool exceeds 1,000

The executed cohort is capped at N=1,000. CN-specific quotas are assigned in
proportion to the complete challenge-pool counts by Hamilton largest remainder
allocation. Within each CN, records are ranked by the smallest SHA-256 score of
`20260822|challenge_selection|CN|graph_sha256|mol_id`. No outcome from the
K=100 trial is used for selection. Membership in the previously locked N=1,000
general stress cohort is retained only as a descriptive manifest field.

## K=100 atom-row reindexing challenge

For every selected record, 100 deterministic uniform permutations of all atom
rows are generated from NumPy `SeedSequence([151120260, selection_order,
variant_index])`. Any supplied connection table is permuted consistently; the
official XYZ inputs in this population do not carry a supplied bond-order
matrix, so their molecular graph is rebuilt from the unchanged fixed
coordinates after each row permutation. Coordinates are not perturbed, rotated,
or reflected in this experiment.

Two record constructions are compared from the same encoded object:

- **signature-only ablation**: validate attachment metadata, sort ligands by
  the candidate invariant sort key, assign L1...Ln, rewrite and sort donor
  relations, round the shape field, and serialize; deliberately skip exact
  equivalent-ligand/donor-orbit enumeration and whole-record minimization;
- **full exact CoordRep**: run the complete frozen candidate canonicalizer,
  including the Cartesian product of all legal ligand-instance and donor-rank
  permutations and selection of the lexicographically minimum complete record.

Primary outcomes are exact variant matches to the unpermuted reference and the
number of structures for which all 100 variants collapse to that reference.
Secondary outcomes are grammar parse/re-serialization round trips, exceptions,
candidate-count stability, and signature-only/full-exact wall times. Timing is
descriptive and is summarized by median, P95, and maximum; it is not a
cross-machine performance benchmark.

## Specificity and fail-closed controls

An explicit non-automorphic negative control compares adjacent (1,2) and
opposite (1,4) two-site attachment sets on the same cyclohexane parent graph.
The attachment-set keys and complete canonical strings must remain distinct
under exhaustive legal relabeling. A repeated-payload record with absent
attachment-set identity must raise the documented fail-closed error. These
controls define the reported false-merge count; cross-record collisions in the
tmQMg cohort are not labelled false merges because distinct source graphs can
legitimately differ only outside CoordRep's recorded first-sphere state.

## Record round trip

The harness independently parses the emitted mononuclear grammar into metal,
shape/status, ordered relation, and ligand-dictionary blocks, validates field
syntax and sequential ligand identifiers, and reconstructs the string from the
parsed fields. A round trip passes only on byte-identical reconstruction.

## Locked pass conditions

The full-exact arm passes only if all selected structures and all K=100
variants match their unpermuted reference, all emitted strings pass the
independent grammar round trip, all per-variant legal candidate counts equal
their reference counts, no variant errors occur, the non-automorphic control
remains distinct, and the missing-metadata control fails closed. The
signature-only result is an ablation outcome and is not a pass gate.
