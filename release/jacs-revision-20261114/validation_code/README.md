# Executable chemical-validation input routes

This directory supplies the final versioned XYZ and conditional-CIF algorithms,
the native CSD adapter, configurable execution wrappers, complete-set comparison,
and fresh-process replay. It contains no deposited CIF, CSD coordinate file,
third-party paper, or exported crystal atom mapping. Generated local outputs can
contain licensed coordinates and must not be uploaded without authorization.

The 84-case register and reported results are in the sibling chemical-validation
data package. Executing an encoder is not itself a chemical-validation success.
The public molecular atom-index references permit direct scoring of the 21
literature-qualified Fe structures. The remaining external molecular structures
stay in the register; encoding them does not promote them to complete-reference
successes. Crystal reference comparisons require the user's authorized source
structures and independently documented reference assignments.

## Environment and quick software controls

Tested with Python 3.9.23 and the versions in `requirements-tested.txt`. CCDC is
not required for the constructed controls, external XYZ route, or source-CIF
coordinate route. A locally authorized CCDC Python API and CSD installation are
required for the `native` route. The crystal panel used April 2025 CSD v6.00;
changing the database version can change deposited chemistry or labels.

```powershell
python -m pip install -r requirements-tested.txt
python -B synthetic_controls.py --output local-software-controls.json
```

The controls generate invented structures in a temporary directory. They test
rotation/translation/atom-order invariance, Si donor admission, simultaneous
heavy-atom screening, exclusion of H as a donor/blocker, periodic donor images,
alternative disorder conditions, MID/SID parity, and the fixed 0.005 Å set
comparison. They do not represent experiments on deposited compounds.

## Obtain and validate external molecular inputs

`public_xyz_sources_44.csv` records all 44 input basenames, checksums, and source
URLs at Molassembler commit
`839c5024787724df5b13850f41fe631b80017074`. The helper fetches those input files
directly from that immutable upstream source into a new user-owned directory:

```powershell
python -B fetch_public_xyz.py --directory local-xyz
python -B run_panel.py --manifest ../chemical_validation/case_register_84.csv --input-directory local-xyz --output-directory local-panel --database CSD --source-release CSD-April-2025-v6.00
python -B validate_molecular_references.py --references ../chemical_validation/molecular_public_reference_indices.json --input-directory local-xyz --output-directory local-molecular-validation
```

If a local Python installation cannot load its system certificate store, pass
`--ca-file trusted-ca-bundle.pem` to the fetch command. TLS certificate and
hostname verification remain enabled; do not bypass them to obtain the inputs.

Adjust these relative paths to the installed chemical-validation package if the
directories are not siblings. `run_panel.py` attempts every register row,
including the 40 crystal rows; without an authorized CSD it retains their input
failures rather than deleting rows. To execute only the 21 directly scored
external molecular cases, use `validate_molecular_references.py` alone. It
requires the exact input checksum before applying any reference indices. The
encoder receives only the XYZ, not the expected coordination number or donor
set. Atom indices are zero-based in the public reference JSON.

The source papers are Proppe and Reiher, *J. Chem. Theory Comput.* **2017**, DOI
10.1021/acs.jctc.7b00235, and Sobez and Reiher, *J. Chem. Inf. Model.* **2020**,
DOI 10.1021/acs.jcim.0c00503. Original supporting-information metadata identifies
CC BY-NC 4.0 terms. The upstream repository has its own BSD-3-Clause code
license. Neither third-party coordinates nor their rights are relicensed under
this directory's MIT license. Generated upstream MOL graphs are not treated as
independent chemical ground truth.

## Execute one input and check a fresh-process replay

```powershell
python -B run_input.py xyz --input local-xyz/01.xyz --metal-element Fe --case-id 01 --output local-01.json
python -B replay_check.py --directory local-replay-01 -- xyz --input local-xyz/01.xyz --metal-element Fe --case-id 01
python -B run_input.py native --refcode TAXDOL --database CSD --source-release CSD-April-2025-v6.00 --case-id own-TAXDOL --output local-TAXDOL.json
python -B run_input.py cif --input user-owned-source.cif --block example_block --case-id own-case --output local-source-conditions.json
```

Use the exact CIF block name, without the leading `data_` delimiter, or omit
`--block` only for a single-block CIF. Outputs never overwrite existing files.
The replay runs the same command in two separate Python processes and compares
the complete generated payloads, retaining failures and logs.

The XYZ route supplies no bond orders, donor hints, expected CN, or oxidation
and spin assignments. It builds local metal-centred records; it does not claim
one whole-assembly canonical identity for a multimetal molecule. `--metal-element`
selects centres, not a donor set. The native route instead serializes deposited
CCDC graph relations without distance completion; its output cannot establish
independent chemical correctness without a separate reference comparison.

## Conditional source-CIF reconstruction

The source-CIF route reads cell, symmetry operations, atom elements/labels,
fractional coordinates, occupancy, and disorder assembly/group selectors. It
does not read deposited bond tables, angle tables, formula annotations, or
paper text as coordination answers. Integer translations reconstruct physical
donor images. Alternative positive groups in the same disorder assembly are
processed separately, with common atoms retained. Per-atom occupancy is retained;
it is not converted to a global conformer probability. Unassigned or unsupported
local disorder context remains explicit.

Use `conditional_source_inputs.csv` from the data package with user-acquired,
checksum-matching deposited CIFs:

```powershell
python -B run_panel.py --manifest ../chemical_validation/conditional_source_inputs.csv --input-directory local-deposited-cifs --output-directory local-conditional-run
```

For the fixed D017 and D019 examples this executes the source-derived conditions,
not just the active CSD disorder group. The reported chemical comparison covers
the 16 reference-qualified configurations and their 72 donor images. D019 also
has two emitted configurations associated with the third BF4 orientation that
do not have equally complete independent chemical references. These are retained
as emitted, unscored configurations, not silently discarded or counted as
agreement. Licensed reference atom mappings are not bundled in this release.

## Compare privately held complete crystal references

`compare_local_sets.py` is a separate, generic comparison tool, not a new source
of gold labels. A reference JSON has a `references` array. Each row contains a
`reference_id`, `target_label`, and complete `donor_images` array; each donor has
`element`, metal-centred `vector_A`, and optionally an exact
`relative_translation`. Source-CIF rows additionally give the exact `condition`
mapping. Native rows can set `include_metal_neighbours: true` if the predeclared
comparison includes metal-metal contacts. XYZ references use `metal_index0` and
`donor_indices0` (or `expected_donor_indices0`).

```powershell
python -B compare_local_sets.py --output-record local-source-conditions.json --references user-owned-complete-references.json --output local-comparison.json
```

Matching is one-to-one at the fixed 0.005 Å positional tolerance. Both missing and
additional donors are failures of complete-set equality. A unique target and
condition must be supplied in advance; the tool never chooses the most favorable
output. It reports CIF typed-incidence counts but does not replace the detailed
published atom-mapping adjudication with a count-only chemical test. Public
aggregate tables alone cannot reproduce restricted crystal comparisons.

## Algorithm provenance and wrapper scope

`ALGORITHM_PROVENANCE.json` records per-file SHA-256 values for the 44 copied
algorithm files. They are byte-identical to the authors' final xyz3/CIF-condition
validation snapshot, including CR-PLS/2. No perception, canonicalization, or
comparison threshold was changed for this publication wrapper.

The preserved XYZ package identifies itself as `1.1.2rc3+xyz3.20261113`. Its
unchanged neighbour submodule retains `xyz2.20261113`; the conditional-CIF module
retains `cifcfg1.20261113`. These are source version labels, not publication-date
rewrites. Admission uses 1.3 times the sum of RDKit covalent radii followed by
simultaneous screening against all original non-H atoms. H is excluded as donor
and blocker. The radius-boundary and Gabriel tolerances remain 1e-12 Å and
1e-12 Å². CR-PLS/2 retains its quantized Gram and signed triple-product records
(six decimals) and exhaustive CN 2-6 donor permutations.

New top-level scripts are `publication-wrapper/1`: they replace private directory
defaults with explicit user inputs, provide process receipts and non-overwriting
output, and keep input/reference comparison separate. Wrapper result JSON fields
are not asserted to be byte-identical to earlier private orchestration logs.
The native CSD adapter source is unchanged from the deposited-graph audit route.
No source structure, saved source result, model weight, or original validation
status is altered by this package.

`WRAPPER_QA_REPORT.json` records publication-wrapper acceptance: all 51 Python
files compiled; 18 constructed controls passed; all 44 external XYZ inputs were
freshly downloaded and checksum-verified; all 21 complete public molecular donor
references matched using those downloads; and a two-process replay agreed.
The licensed native route was also executed as a one-case smoke test, not scored
as an independent chemical validation. Raw QA outputs are not distributed here.
