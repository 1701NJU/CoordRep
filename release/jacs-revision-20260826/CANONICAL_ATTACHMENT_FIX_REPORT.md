# Mononuclear atom-renumbering canonicalization repair

Date: 2026-08-22

Status: integrated as CoordRep 1.1.2rc3 in the frozen 2026-08-23 public
release. This report preserves the pre-integration repair rationale and test
results.

## Failure being repaired

The previous exact duplicate-ligand search used
`tuple(sorted(zip(lig.donor_elements, lig.attach_atoms)))` in its equivalence
key. `attach_atoms` contains source/global atom indices. Two graph-equivalent
ligand instances therefore ceased to be exchangeable after a harmless atom
renumbering. Simply deleting `attach_atoms` would be unsafe: two copies of the
same unmarked ligand SMILES can coordinate through non-automorphic atom sets.

There was a second, independent identity leak. The serializer emitted only
the unmarked parent ligand SMILES. Even if the canonicalizer correctly refused
to exchange two non-equivalent attachment sets, the final strings could still
be identical when no relation token exposed the difference.

A subsequent smoke test exposed a third ordering leak: molecular weight was
placed before canonical SMILES in the ligand sort key. TOHTIV contains two
nonidentical N,O-bidentate ligands with the same formula. Atom traversal changed
the last floating-point bit of the RDKit mass sum (`131.175` versus
`131.17499999999998`), reversing their L1/L2 order. Molecular weight is now
removed from the identity sort; the discrete payload and attachment fields
provide the deterministic total order.

## Repair design

Each graph-supported `LigandModule` now carries three source-index-free fields.

1. `attachment_set_key`: canonical ligand graph with every coordinated atom
   coloured. This is a joint set key, so it distinguishes, for example,
   adjacent and opposite two-site attachment on the same ring even when all
   individual atoms belong to one orbit in the unmarked graph.
2. `donor_attachment_keys`: one canonical coloured-graph key per donor. All
   attachment atoms remain coloured and the queried donor receives a second
   colour. Equal keys therefore define the donor orbits under the stabilizer
   of the complete attachment set.
3. `donor_rank_orbits`: same-element donor ranks that remain tied after both
   the donor-orbit key and the quantized local geometric key are applied.

The ligand-instance equivalence key is now:

```text
(payload fields, attachment_set_key, sorted donor_attachment_keys)
```

It contains no global atom index. Exact canonicalization enumerates the
Cartesian product of:

- legal permutations of equivalent ligand instances; and
- legal rank permutations within declared `donor_rank_orbits`.

The lexicographically smallest **whole serialized record** is retained. Donor
ranks are never freely permuted across different elements, different
attachment-set stabilizer orbits, or resolved geometric keys.

The serializer now uses the donor-marked canonical ligand graph as the SMILES
payload whenever `attachment_set_key` is available. Thus attachment identity
is present in the final record, not only in an internal comparison key.

For repeated ligand payloads without an attachment-set key, canonicalization
fails closed with an explicit `ValueError`. This avoids silently merging a
positional coordination isomer or claiming an invariance that the adapter has
not established. Unique formula-fallback ligands retain the existing explicit
fallback path.

## Minimal source diff

Only the following mononuclear-core files were changed:

- `coordrep/core.py`: adds the three ligand-local attachment/orbit fields and
  places `attachment_set_key` in the deterministic ligand sort key; removes
  floating molecular weight from that key.
- `coordrep/graph/ligand_module.py`: generates joint attachment-set and
  per-donor orbit keys from canonical coloured RDKit ligand graphs.
- `coordrep/graph/donor_sites.py`: orders donors by ligand-local topology then
  the existing quantized local geometric information, and records only legal
  unresolved tie orbits.
- `coordrep/geometry/rel_config.py`: carries the attachment descriptors through
  the alternate donor-site construction route.
- `coordrep/canonical/canonicalize.py`: removes global indices from ligand
  equivalence, validates missing metadata, and performs the combined exact
  ligand/donor-orbit Cartesian search.
- `coordrep/serialize/to_string.py`: serializes donor-marked attachment
  identity so the negative control cannot pass by output omission.

No multimetal, haptic/pi, periodic, or identity-level implementation file was
changed. The six mononuclear-core files listed above are the complete rc3
source diff relative to the frozen rc2 release.

## Regression tests

New tests are in `tests/test_attachment_orbit_canonicalization.py`; the former
internal-donor negative control in
`tests/test_core_exact_equivalent_relabeling.py` is upgraded to a positive
orbit-aware test.

Required cases and observed results:

| Case | Result |
|---|---:|
| ABAZAH raw molecule, 16 seeded full-atom permutations | 16/16 bit-identical canonical strings |
| TOHTIV raw molecule, 20 seeded full-atom permutations | 20/20 bit-identical canonical strings |
| Equivalent ligands with different global attachment indices | same equivalence key and same canonical string |
| Same parent payload, adjacent vs opposite non-orbit attachment sets | different attachment keys and different final strings |
| Three symmetric bidentate ligands, all 2^3 internal donor-rank flips | one canonical string |
| Symmetric same-orbit/same-geometry donor pair in the adapter | one legal two-rank orbit declared |
| Repeated payload with missing attachment identity | explicit fail-closed error |
| Prior simultaneous duplicate-group regressions | all pass |

Targeted suite: `12 passed`.

Copied publication suite excluding its repository-layout-only manifest test:
`110 passed, 2 skipped`.

The omitted `test_manifest_paths.py` expects `DATA_MANIFEST.md` two directories
above the copied package and is not a source-code failure. Its data-backed
v2-beta tests were included after copying the unchanged `revision_results`
fixtures into the candidate package used before public integration.

## ABAZAH fixture provenance

- Local test fixture (not redistributed): `tests/fixtures/xyz/ABAZAH.xyz`
- Fixture SHA-256:
  `dfe24592992d4dc83ff9e1fa2fbe07dbe84e135fcf87bf11c0fd398f383db0ad`
- Extracted unchanged from:
  `revision_experiments/source_cache/official_tmqmg_github_xyz/tmQMg_xyz.zip`
- Source ZIP SHA-256:
  `e0d15a70bcba294717cd9f9792e7fac99ef0c5c61c3a6e08dcc8a8643f53660a`
- ZIP member: `xyz/ABAZAH.xyz`

The TOHTIV regression fixture was extracted unchanged from the same source
ZIP (`xyz/TOHTIV.xyz`) and has SHA-256
`6d0189688fb5bbf09125c42673b80ed34b686e63ee00ff4878a7c5f4ed22f918`.

## Integration and evidence boundary

The reviewed six-file repair is integrated in CoordRep 1.1.2rc3. The planned
N = 1,000, K = 100 general-invariance experiment and the separate 100,000-
variant legal-orbit challenge were subsequently completed and are archived in
`canonicalization/general_invariance/` and
`canonicalization/legal_orbit_challenge/`. The repair can change molecular
strings and L0–L3 hashes relative to rc2; the frozen Figure 4 property
benchmarks and molecular-family statistics were not rerun under rc3 and remain
explicitly labelled rc2 downstream evidence.
