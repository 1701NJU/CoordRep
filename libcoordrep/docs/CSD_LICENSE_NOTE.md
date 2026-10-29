# CSD License Note

## Data Redistribution Policy

Raw CSD (Cambridge Structural Database) data is **not redistributed** in this repository due to CCDC licensing restrictions.

## What Is Included

- **Aggregate statistics**: filter waterfall, rejection reasons, identity layer match rates
- **Anonymized summaries**: metal/CN/shape distributions without identifiable coordinates
- **Scripts**: all CSD processing scripts are included and documented

## What Is NOT Included

- `csd_retained_entries.jsonl` (full CoordRep strings derived from CSD)
- Raw CSD coordinates (`.cif`, `.mol2`, `.xyz`)
- Bulk refcode-to-structure mappings
- Any data from which individual CSD structures could be reconstructed

## Requirements for CSD Scripts

To run CSD-related scripts, you need:

1. A valid CCDC license
2. Local CSD Python API installation (`ccdc` package)
3. Access to the CSD database files

```bash
# Example: run CSD external audit (requires licensed CSD installation)
python scripts/csd_external_audit.py --max_entries 1000
```

If you do not have a CSD license, all CSD-derived **aggregate results** are provided in:
```
revision_results/csd_external_summary_only/
```

## Contact

For questions about CSD data access, contact the CCDC: https://www.ccdc.cam.ac.uk/
