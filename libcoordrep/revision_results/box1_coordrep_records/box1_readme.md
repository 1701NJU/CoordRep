# Box 1: Representative Complete CoordRep Records

## Purpose

These records demonstrate the full CoordRep-v1 representation for three
classes of coordination compounds, directly addressing the reviewer concern
that "full syntax is insufficiently documented."

## Data provenance

- **Source**: Cambridge Structural Database (CSD)
- **Serializer**: `coordrep_v1` (`coordrep.encode.encode_molecule` →
  `canonicalize_complex` → `to_string()`)
- **Validation**: `coordrep_tools.validate.is_valid_coordrep(strict=True)`
- **Identity keys**: `coordrep.identity.extract_identity_keys()`

## Selection criteria

| Example | Criteria |
|---------|----------|
| A | CN=4, Pt, shape_best=SP, all monodentate, has trans tokens, not boundary |
| B | CN=6, Co, shape_best=Oh, has bidentate ligand(s), has trans tokens |
| C | CN=6, has {fm:fac} or {fm:mer} token, Oh preferred |

## Commands to reproduce

```bash
/data/miniconda3/envs/1701/bin/python scripts/export_box1_coordrep_records.py
```

## Guarantees

- All records are **direct serializer output**; no hand-written strings.
- Every record passes `parse_valid=True` and `roundtrip_valid=True`.
- No placeholders ("...", "TBD", "unknown") appear in any output.
- Raw CSD coordinates are **not** redistributed.

## Files

| File | Description |
|------|-------------|
| `box1_record_index.csv` | Summary index for all examples |
| `box1_display_records.md` | Main-text Box 1 formatted for typesetting |
| `box1_full_records.jsonl` | Machine-readable full records |
| `box1_validation_report.csv` | Per-record validation gate results |
| `box1_si_full_records.txt` | Supporting Information full records |
| `box1_readme.md` | This file |
