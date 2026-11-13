# Literature-referenced chemical validation

This package contains the authors' derived comparison results for the fixed
84-structure panel: all 44 published Fe test structures, 20 periodic structures,
and 20 disorder-labelled structures. Every case remains in the register.
Processing, numerical reproducibility, and agreement with a complete chemical
reference are reported as separate endpoints.

| Endpoint | Result | Meaning |
|---|---:|---|
| External molecular inputs | 44 structures; 47 local Fe records | All inputs processed and reproduced in a fresh process. |
| Published single-Fe coordination numbers | 40/41 agree | The three two-Fe inputs are outside this denominator. Case 28 compares different optimized and experimental states. |
| Complete same-optimized-state donor references | 21/21 local sets; 103/103 assigned donors | No missing or additional donor in these qualified sets (Table S16). |
| Complete crystallographic local references | 17/19 exact sets from 12 entries | All 107 assigned nonmetal contacts retained; four additional longer Cu-O contacts occur across two Cu sites (Table S17). |
| Metal-metal relations | 6/6 directed Au-Tl relations | Counted separately from the 107 nonmetal contacts. |
| Complete conditional disorder references | 16/16 local alternatives; 72/72 physical donor images | Twelve ZnN2I2 and four NaF6 alternatives; zero missing/additional donor images in the qualified comparisons (Table S18). |
| Additional BF4 orientation | 2 local candidates retained, not scored | The third BF4 orientation lacks an equally complete independent local-donor reference. |
| Molecular re-expression controls | 470/470 centre-level trials; another 2,000/2,000 for case 13 | Stability tests, not additional chemical-reference determinations. |
| Conditional-CIF re-expression controls | 24 adapter executions; 216/216 local-record comparisons | Distinct execution and comparison units; not 216 structures. |

The 84 input cases are not 84 complete reference determinations. References
were assembled by the authors from the cited papers and deposited structures;
this is not a new human-blind adjudication or an unseen holdout benchmark. The
complete-reference endpoints above use different qualified subsets and are not
pooled into one accuracy estimate.

## Files

- `case_register_84.csv`: every fixed case, public identifiers, input selectors,
  reference qualification, and outcome. XYZ basenames and SHA-256 identify the
  external inputs without redistributing their coordinates.
- `molecular44_execution_and_CN.csv`: all 44 molecular inputs, 47 local records,
  and the 41 single-centre coordination-number comparisons.
- `Table_S16_molecular_donor_references.csv`: the 21 qualified optimized-state
  primary-donor comparisons.
- `molecular_public_reference_indices.json`: zero-based Fe and expected-donor
  indices for those 21 independently assigned references, tied by SHA-256 to
  the external published XYZ files. These are not native CSD atom mappings.
- `Table_S17_crystallographic_references.csv`: the 12 entries providing 19
  complete local comparisons, including both Cu differences.
- `Table_S18_disorder_conditional_references.csv`: the two source-CIF cases,
  their 16 qualified local alternatives and 72 donor images.
- `conditional_source_inputs.csv`: source-CIF basename, block selector and
  SHA-256 for separately acquired inputs; no CIF content is included.
- `references_S1_S16.csv` and `REFERENCES.txt`: short references and DOI links.
  These S-numbers follow the SI order; they are not replacements for the
  manuscript's existing numbered Molassembler and CSD references.
- `SCOPE_AND_PROTOCOL.md`: comparison rules, endpoint definitions and the
  interpretation of case 28, Cu contacts and disorder alternatives.
- `validation_summary.json`: machine-readable numerical results.
- `source_commitments.json`: SHA-256 commitments to archived author-side
  comparison reports, not substitutes for the licensed source files.
- `verify_summary_tables.py`: standard-library consistency and checksum test.
- `SHA256SUMS.txt`: file hashes for this derived-results package.

Run `python verify_summary_tables.py` to verify the public tables and file
hashes. This checks the reported arithmetic and internal consistency; it does
not reconstruct atoms or independently establish chemical assignments. The
companion executable implementation is in `../validation_code/`.
The 21 external molecular donor sets can be recomputed against the public
index references after acquiring the matching published XYZ inputs. Raw-input
reruns require the separately acquired XYZ/CIF inputs or an
appropriately licensed CSD installation, as specified by that implementation.

## Redistribution

Only public identifiers, bibliographic metadata, author-derived chemical
descriptions, counts, scope labels and comparison results, and the 21 external
public-XYZ index references are included. The
tables are the authors' comparison results, not an official CCDC or publisher
reference dataset. No raw CIF, third-party or CSD coordinates, publisher PDF,
native-CSD atom-label/index mapping, private local path, credential or licensed
site-level structural ledger is distributed. The repository-wide
`CSD_REDISTRIBUTION_NOTICE.md` also applies. Obtain source structures through
their original providers and comply with their terms.
