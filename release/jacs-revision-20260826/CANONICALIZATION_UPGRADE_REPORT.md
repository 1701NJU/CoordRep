# CoordRep canonicalization upgrade: implementation and claim boundary

Date: 2026-08-03
Working copy: `revision_experiments/coordrep_v2beta_full_canonical_upgrade_20260803`
Source archive: the read-only local source archive used during the revision audit (not redistributed)

## Bottom line

The label/permutation defect can be repaired without changing the three-part
Results architecture.  The upgraded code now provides exact canonical
labelling of the information contained in a **geometry-resolved coordination
record graph** for:

- the mononuclear validated core at CN <= 6;
- finite multinuclear metal--site--ligand record graphs with at most 12 metals;
- haptic/pi coordination-site objects treated as unordered donor atom sets.

This is not yet a full chemical-species identifier.  The current CSD v2-beta
adapter still emits element-only or `[*]` ligand placeholders in some paths;
no canonicalization algorithm can recover ligand connectivity, bond order,
stereochemistry, or charge that the adapter did not encode.

## Defects reproduced before the patch

### Symmetric multinuclear graph

The former fixed-round WL sort leaves automorphic metal vertices tied and then
inherits Python input order.  A homogeneous Cu4 cycle gives:

- 24 exhaustive input metal permutations;
- 3 distinct serialized metal graphs;
- 8 permutations in each output class.

### Multiple duplicate-ligand groups in the mononuclear core

`coordrep/canonical/canonicalize.py::_apply_equivalent_relabeling` formerly
optimized one duplicate-SMILES group at a time.  It did not evaluate the
Cartesian product of two or more group permutations.  A CN4 record containing
two equivalent N ligands, two equivalent chloride ligands, and one trans
relation produced two strings under the four within-group traversal orders.
The same routine also silently truncated groups larger than six ligands.

### Source atom labels in haptic records

`CoordinationSite.site_signature()` and the haptic serializer included raw CSD
atom labels.  Renaming a Cp-ring atom therefore changed the site ordering and
L0 string even though the abstract eta5 site was unchanged.

### Validation false positive

The old `metal_order_invariant` gate created `shuffled_order` but never applied
it to `rec3.metals`.  It tested a relabeling in the original list order, not an
actual metal-order perturbation.

## Implemented changes

### Expanded molecular records

File: `coordrep/v2beta/canonicalize.py`

- Added exact individualization--refinement over metal vertices.
- Refines colors with typed metal edges and complete metal/site/ligand
  incidence context.
- Scores every discrete leaf with a source-label-free code and selects the
  lexicographic minimum.
- Uses exact transposition-automorphism checks only for safe symmetry pruning.
- Canonicalizes ligand and site objects without redistributing sites between
  duplicate ligands.
- Replaces source donor labels by site-local canonical tokens while preserving
  original labels in `site.meta['source_donor_labels']`.
- Normalizes haptic centroid labels and the single-metal target label.

File: `coordrep/v2beta/core.py`

- `canonical_atom_key()` no longer uses `donor_atoms` source labels.
- It uses optional label-free `donor_canonical_keys`; otherwise it uses the
  unordered donor-element multiset, which is the strongest invariant present
  in the v2-beta record.

File: `coordrep/v2beta/serialize.py`

- Added a complete `[Sites:]` block for multinuclear records.  Earlier strings
  contained terminal-site labels but omitted their full ligand/donor/target
  incidence, so distinct in-memory record graphs could share an L0 string.

File: `coordrep/v2beta/csd_v2beta_adapter.py`

- Added label-free atom-graph refinement signatures based on element, charge,
  degree, bond type, and recursively refined neighboring colors.
- Source CSD labels are lookup/provenance fields only, not identity fields.
- This does not replace placeholder ligands with full canonical ligand graphs.

File: `coordrep/v2beta/validate.py`

- The metal list is now actually shuffled and arbitrarily relabeled.
- Donor atoms are shuffled and renamed.
- Eta sites are explicitly rotated and reversed before relabeling.
- Referential integrity of metal, site, ligand, edge, and bridge references is
  checked.

File: `coordrep/v2beta/cases_multi.py`

- Corrected `multi_028` (Re2Cl8(2-)): four sites referenced an absent `L2`.
  The corrected record contains eight chloride ligand objects, one per site.

### Mononuclear core

File: `coordrep/canonical/canonicalize.py`

- Duplicate-ligand groups are defined by the full ligand record rather than
  SMILES alone.
- The exact Cartesian product of all within-group permutations is evaluated.
- The former `len(group) > 6` silent truncation was removed.
- Searches above 100,000 candidates fail explicitly as outside the validated
  CN <= 6 core instead of returning an approximate, input-dependent string.

## Verification results

### Expanded records audit

Command:

```powershell
python .\scripts\audit_v2beta_exact_canonicalization.py
```

Results:

| Audit | Result |
|---|---:|
| Cu4 exhaustive permutations | 24/24 collapsed to one string |
| Legacy Cu4 distinct outputs | 3 |
| Upgraded Cu4 distinct outputs | 1 |
| Curated multinuclear records | 55 |
| Random nuisance variants per multinuclear record | 100 |
| Multinuclear nuisance trials | 5,500/5,500 passed |
| Curated haptic records | 33 |
| Rotation/reversal/rename variants per haptic record | 100 |
| Haptic nuisance trials | 3,300/3,300 passed |
| Non-isomorphic Cu4 cycle/path | different L0 keys |
| Multinuclear validation gates | 55/55 passed |
| Haptic validation gates | 33/33 passed |

The machine-readable result is
`revision_results/canonical_invariance_upgrade/canonical_invariance_audit_summary.json`.

### Core and regression tests

Command:

```powershell
python -m pytest -q tests\test_core_exact_equivalent_relabeling.py tests\test_v2beta_exact_canonicalization.py tests\test_v2beta_records.py tests\test_basic_encoding.py tests\test_candidate_repairs.py tests\test_identity_keys.py tests\test_csd_hapticity_spellings.py tests\test_version_alignment.py
```

Result: **68/68 passed**.

The core red-team specifically includes:

- CN4 with two simultaneous duplicate-ligand groups: 4/4 traversals collapse;
- CN6 with three duplicate-ligand groups: 8/8 within-group traversals collapse;
- 100 random global ligand-list orders and arbitrary source-ID renamings:
  100/100 collapse;
- CN9 factorial stress: rejected explicitly, never approximated.

An additional negative-control test flips the two internal N donor ranks of
each ligand in a symmetric Co(en)3 record.  The eight flip patterns still give
two strings.  This is expected from the present data model: `LigandModule`
does not store atom-graph automorphism orbits, so canonicalization cannot tell
whether two same-element donors are symmetry-equivalent or constitutionally
distinct.  Canonical identity therefore requires the adapter's fixed donor-rank
convention; arbitrary internal donor-rank renaming is not yet supported.

The earlier Windows fixture errors were unrelated to the algorithm:
`tests/test_v2beta_records.py` used locale-default `open()`, causing GBK to
decode UTF-8 JSONL.  Both fixture readers now specify `encoding='utf-8'`.

A full repository test attempt reached the unrelated coordinate-CShM tests and
then the local MKL build aborted inside `numpy.linalg.svd`.  The canonical and
selected regression suites above complete normally.

## Bounded real-CSD ligand-payload feasibility audit

Two additional scripts tested whether the 55 multinuclear and 33 haptic
templates could immediately be upgraded from element/placeholder ligands to
full atom-mapped ligand objects:

```powershell
`python scripts/audit_csd_ligand_components_phase1.py` (run inside a licensed CCDC environment)
python .\scripts\audit_rdkit_ligand_components_phase2.py
```

Phase 1 used the local April 2025 CSD installation / API 3.4.1
(1,371,757 entries), removed all transition
metals, and treated each donor-bearing non-metal connected component as one
ligand object.  The result is a decisive provenance limitation:

| Template scope | Requested | Refcode found | Expected metal set matched |
|---|---:|---:|---:|
| Multinuclear | 55 | 3 | 0 |
| Haptic | 33 | 3 | 2 |

Thus, the 55/33 collection is a programmatic template set, not 88 verified
current-CSD records.  Most identifiers are absent; several identifiers that do
exist refer to structures with a different metal count/composition.

The two metal-scope-matched entries were FEROCE01 and CPMNCO.  They yielded six
donor-bearing components.  RDKit successfully produced canonical payloads and
source-label-free donor symmetry-rank multisets for the three carbonyl
components, with 60/60 atom-order permutation trials passing.  All three eta5
Cp components failed Mol2/SMILES sanitization because metal removal leaves a
neutral five-membered aromatic fragment (`c1cccc1`) without the charge/bonding
assignment needed for a chemically valid standalone ligand.  Consequently:

- matched components: 6;
- full payload + atom-map invariant: 3;
- eta5 sanitization failures: 3;
- completely successful matched entries: 0/2;
- multinuclear full-payload success rate: not estimable (0 matching entries).

Machine-readable outputs are in
`revision_results/csd_ligand_component_feasibility/phase1_summary.json` and
`phase2_summary.json`.

This is a **NO-GO** for claiming full ligand chemical identity from the current
55/33 evidence.  Repairing Cp fragments by assigning charge or altered bond
orders would be a chemical inference and must be specified, validated, and
audited rather than silently forced by the serializer.

## Manuscript-ready scope statement

Recommended Results/SI text:

> For the expanded molecular scope, each finite complex was represented as a
> colored metal--coordination-site--ligand incidence graph. Exact
> individualization--refinement was used to resolve residual automorphism
> classes, after which the lexicographically minimal serialization defined the
> geometry-resolved record state. Source metal, site, ligand, and donor-atom
> labels were treated as provenance rather than identity fields. In an
> exhaustive symmetric Cu4 test, all 24 metal input permutations collapsed to
> one record, whereas the preliminary WL-only beta implementation produced
> three. Across 55 curated multinuclear and 33 haptic record templates, 8,800
> nuisance variants combining container shuffling, arbitrary source-label
> renaming, and eta-site rotation or reversal yielded bit-identical records and
> L0 keys.

Recommended reviewer response:

> We thank the reviewer for prompting a stricter treatment of graph
> canonicalization. During revision we identified two unresolved automorphism
> cases: the preliminary multinuclear implementation inherited input order
> when WL refinement left symmetric metal vertices tied, and the mononuclear
> implementation optimized multiple duplicate-ligand groups independently.
> We replaced these steps by exact global canonical searches within the stated
> record domain and expanded the invariance tests to include symmetric metal
> graphs, simultaneous duplicate-ligand groups, arbitrary identifier renaming,
> and cyclic eta-site rotations/reversals. The revised claim is canonical
> identity of a well-posed geometry-resolved coordination record, not identity
> of a chemical species in the InChI sense.

## Claims that remain unsafe

Do not write any of the following without an additional adapter/data audit:

- "canonical chemical-species identifier";
- "complete ligand identity" for CSD adapter outputs containing `[ * ]` or
  element-only placeholders;
- "all multinuclear and haptic CSD chemistry is canonical" based only on the
  55/33 curated templates;
- invariance to alternate tautomer, protonation, disorder, occupancy, bonding,
  or oxidation-state assignments;
- internal donor-rank invariance for an asymmetric multidentate ligand unless
  donor ranks are derived from a canonical atom-mapped ligand graph;
- arbitrary internal donor-rank invariance even for a symmetric ligand until
  atom-graph automorphism orbits are stored by the adapter;
- CN > 6 mononuclear-core coverage under the exact duplicate-ligand search.

## Work still required for a full ligand-chemical-identity claim

1. Replace `_safe_smiles()` placeholders with a canonical, atom-mapped ligand
   graph containing bond orders, formal charges, isotopes, and stereochemistry.
2. Derive donor ranks from canonical ligand-graph atom orbits rather than CSD
   labels or traversal order.
3. Re-run the actual CSD multinuclear and haptic cohorts under random atom,
   metal, ligand, and container permutations; the existing 218/218 and 398/398
   topology counts are not substitutes for this invariance audit.
4. Report conversion failures and unsupported ambiguity classes separately.

Until those steps are complete, the strongest defensible and still substantial
claim is: **CoordRep provides exact, auditable canonical serialization of the
geometry-resolved coordination record graph across its validated core and
expanded molecular record domains.**
