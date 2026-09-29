# Periodic-v3 MID/SID public-ZIP runner and worked example

This **public code-only addendum** offers both a one-deposition `ABADUG`
worked example and a batch runner
for the CCDC CSD MOF Collection ZIP through a locally authorized CCDC Python
API. It writes no CIF, coordinates, metal–donor site rows, or refcode ledger.
The batch runner prints only aggregate JSON and safe progress/error counters;
MID/SID hashes stay in memory until aggregation.

The collection archive is available for non-commercial use under
[CC BY-NC-SA 4.0 through CCDC](https://www.ccdc.cam.ac.uk/support-and-resources/Downloads/).
The exact frozen archive has SHA-256
`ABF9F2A7288FD83E2A59F19B717C202414E744A1B2C8C1F9EEC1F909554346FA`.
If CCDC has replaced that download, obtain the frozen version from an
authorized source before comparing against these commitments; do not bypass
the checksum.

## Scope and provenance

`periodic_v3/canonical_periodic_site.py` and
`periodic_v3/canonical_csd_quotient.py` are byte-for-byte copies of the frozen
research-run sources, with SHA-256 values
`292AFD4DB89D1DD552F9A03433B207FD648BB30E9DED717BFB8DBCE17D4694AC`
and `8C9DE81E1EBABDECE083BD6393786357BDC35CE3D6AD3142AAB205E2ECB19EEA`.
The runner checks both before processing. Eight vendored canonicalizer tests exercise
donor permutations, origin shifts, image gauges, quotient-orbit renaming,
unimodular basis changes, proper rotations, reflections, orbit partitions,
genuine image changes, and entry-multiset multiplicity. Two further pure
tests check the public-summary comparator.

`run_abadug.py` is a deliberately narrow **single-entry extraction adapter**.
It reads the
public archive's `Framework details.csv` to locate `abadug_P1.cif`, calls
`ccdc.crystal.Crystal.from_string`, clusters symmetry-expanded images at the
frozen 0.002 Å tolerance, uses the CCDC **explicit molecular bond graph** for
metal–donor incidences, reconstructs translation-labelled quotient edges,
checks the 0.001 Å vector-witness and 0.60–3.50 Å distance bounds, and emits
CR-PLS/2 IDs. It does **not** add distance-inferred metal–donor contacts.
The adapter is limited to ABADUG, for which all six sites are emitted under
the frozen chemical gates. The separate `run_batch.py` plus
`periodic_v3/adapter.py` port the frozen primary site-emission gates across
the whole pinned MOF archive: metal/donor domain, CN 2–6, haptic/π and
connected/repeated carbon exclusions, metadata's unreliable-chemistry flag,
explicit distance bounds, bond-role conflicts, and a Cartesian vector
witness. It omits corrected CShM and optional distance-QC reporting because
those are not primary MID/SID emission gates. The private v2 ledger is not
read for the public calculation. The 400-entry re-expression challenge and
secondary by-shape analysis remain outside this bundle.

MID is the SHA-256 of a canonical JSON payload containing metal/CN,
donor-element/bond-role/site labels, the donor-orbit equality partition, and
the six-decimal quantized labelled Gram matrix
`G_ij = v_i · v_j` (Å²). SID includes those fields **plus the complete
six-decimal quantized signed triple products**
`T_ijk = (v_i × v_j) · v_k` (Å³). Both payloads are minimized over all legal
donor permutations (at most 720 for CN ≤ 6). MID is parity-free under `O(3)`;
SID distinguishes local orientations under proper rotations `SO(3)` when the
geometry and labels permit. The exact algorithm uses signed magnitudes, **not
only a tuple of signs**. These are local first-sphere keys, not framework-net
or macroscopic chirality assignments.

## Requirements and commands

- Python ≥ 3.9 and NumPy ≥ 1.24. Pure unit tests ran with Python 3.12.7 and
  NumPy 1.26.4.
- A locally authorized CCDC Python API installation to parse the CIF and its
  explicit bond graph. The successful integration run used Python 3.9.23,
  CCDC API 3.4.0, and NumPy 2.0.2 in the original `csd_env`. The public bundle
  does not include or install CCDC software or a CSD database.
- The user's own pinned `CSD_MOF_Collection.zip`; it is **not** bundled here.

In PowerShell, from this directory:

```powershell
python -B -m unittest discover -s tests -v
$pythonCsd = 'C:\path\to\authorized_ccdc_env\python.exe'
$mofZip = 'C:\path\to\CSD_MOF_Collection.zip'
& $pythonCsd -B run_abadug.py --mof-zip $mofZip
& $pythonCsd -B run_batch.py --mof-zip $mofZip --limit 10
& $pythonCsd -B run_batch.py --mof-zip $mofZip --start-index 1000 --limit 100 --workers 2
$env:COORDREP_MOF_ZIP = $mofZip
& $pythonCsd -B -m unittest discover -s tests -v
```

Set the two paths to your own authorized interpreter and download. The first
command runs the ten pure-Python/NumPy tests and skips the
two optional CIF integration tests. The last command runs all twelve tests when
`COORDREP_MOF_ZIP` is set in the CCDC environment.

For an explicitly requested full census, use `--all` instead of `--limit`;
the runner will verify the pinned archive and canonicalizer hashes first.
Use `--workers 1` by default, or 2–4 only when the local CCDC licence and RAM
permit independent parser processes. No private ledger is required:

```powershell
& $pythonCsd -B run_batch.py --mof-zip $mofZip --all --workers 2 |
    Set-Content -Encoding utf8 batch_aggregate.json
& $pythonCsd -B compare_public_summary.py --batch-json batch_aggregate.json `
    --public-summary 'C:\path\to\MOF_stereochemical_pairing_summary.json'
```

The public summary is in the CoordRep `jacs-revision` release's
`figures/Figure6/` directory. The comparison covers the primary atlas and
pairing fields, **not** `shape_conflicting_metric_ids` because this minimal
runner does not compute corrected CShM. The batch output itself does not
write a file; the example PowerShell pipeline saves only aggregate JSON.

The optional `--compare-frozen-entries PATH` is for the original lab's local
QA only. It compares per-entry status, orbit/metal/emission/periodic counts,
and MID/SID multiset hashes, but prints only exact-match totals. Do not
publish that private ledger or any row-level output. The runner uses the
metadata's explicit `Yes` as the Sohncke flag; every other value is labelled
**not flagged Sohncke**, not asserted to be crystallographically non-Sohncke.
The reported gate counts are nonexclusive.

On the pinned input, the adapter exits zero and reports:

```json
{
  "status": "PASS_ABADUG_public_ZIP_one_CIF_matches_frozen_v3_commitments",
  "refcode": "ABADUG",
  "result": {
    "periodic_atom_orbits": 282,
    "coordination_numbers": {"5": 6},
    "emitted_sites": 6,
    "sites_with_nonzero_translation": 3,
    "entry_mid_multiset_sha256": "fd1969008e6147f3077418b22588d2084bfb4767452f2afce82ef71e7fadbe61",
    "entry_sid_multiset_sha256": "b8d4823c52c4435e0be8052177f78b50e7187e48fc528ebf31b4552a9bb67f8d"
  }
}
```

The observed ABADUG run matched all five frozen commitments. In separate
batch smoke checks against the original lab's private frozen entry ledger,
the first 100, metadata indices 1000–1099, and indices 10000–10099 each
matched **100/100 entries exactly** across the seven per-entry checks.

## Full-archive CIF rerun

`FULL_ZIP_AGGREGATE_20260929.json` is the aggregate-only output of a
completed `run_batch.py --all --workers 2` execution on the pinned public
15,906-CIF archive, using the same frozen canonicalizer and emission gates
as the worked example. The run used CCDC Python API 3.4.0 and NumPy 2.0.2;
it did **not** read the private output ledger (`private_frozen_QA` is
`not_requested`). It processed 15,905 entries and recorded one
`site_extract:ValueError`, consistent with the frozen processed-entry
denominator. The CIF rerun emitted 172,332 local sites, with 74,354
nonzero-translation states, 36,552 unique MIDs, 60,139 unique SIDs, and
23,587 MIDs having two observed SIDs. It also recovered 7,721 paired
entries among 9,189 state-bearing entries not flagged Sohncke and none of
1,759 flagged Sohncke entries.

The aggregate file retains the runner's conservative
`full_zip_run_requires_claim_review` status. The separate public-summary
comparison is the claim check: all **11** primary Figure 6 top-level fields,
including nested pairing and family counts, match the frozen public source
with no mismatch. Run it from the repository root:

```powershell
python -B release/jacs-revision-20260929/periodic/compare_public_summary.py `
  --batch-json release/jacs-revision-20260929/periodic/FULL_ZIP_AGGREGATE_20260929.json `
  --public-summary release/jacs-revision-20260826/figures/Figure6/MOF_stereochemical_pairing_summary.json
```

This is a full-archive **raw-CIF rerun of the primary MID/SID atlas**, not
merely reaggregation of pre-existing output. Its scope still excludes the
secondary corrected-CShM `shape_conflicting_metric_ids` calculation and
the separate 400-entry representation-change challenge. The public file
contains only collection-level aggregates and no CIF, coordinates,
per-site rows, refcode ledger, or private path.
