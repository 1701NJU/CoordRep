# AGOTIA — Grammar Repair (Fig. 5C, Case C)

## What this case demonstrates

AGOTIA is a Re(VII) CN=7 complex. A single bracket is removed from the
CoordRep string (position 129 inside L1 SMILES), making the
string **grammatically invalid** (bracket imbalance).

CoordRep's explicit bracket/delimiter/field grammar enables a **deterministic
rule-based repair** that:
1. Detects the bracket imbalance
2. Re-inserts the missing bracket
3. Restores syntactic validity
4. Preserves 5/6 semantic fields (metal, CN, ox, shape, stereo)

The only field not perfectly recovered is the ligand L1 SMILES (bracket
position may shift), which is expected for a character-level deletion.

## Rendering guidance

- Use `AGOTIA_first_sphere.xyz` for a small Re/CN7 thumbnail in the card corner.
- The main evidence is the **repair excerpt** (`AGOTIA_repair_excerpt.txt`),
  not the 3D structure.
- Highlight Re in a distinct color.

## Key numbers

| Metric | Value |
|--------|-------|
| Corruption | missing `)` at position 129 |
| Validator before | INVALID |
| Validator after  | VALID |
| Fields preserved | 5/6 |
