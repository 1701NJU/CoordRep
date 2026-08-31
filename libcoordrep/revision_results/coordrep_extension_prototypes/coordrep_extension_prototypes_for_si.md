# CoordRep Extension Prototypes — Supplementary Feasibility Tests

## Scope

CoordRep v1, as described in the main text, is restricted to **mononuclear,
atom-resolved η1 coordination snapshots**. The prototypes below demonstrate
that the CoordRep grammar is **extensible** to multinuclear and haptic/π-bonded
systems, but do not claim full production support.

These tests respond to reviewer concerns about multinuclear complexes and
π/haptic ligands.

**Important:**
- This is a feasibility demonstration only.
- We do not claim full support for MOFs, clusters, metallocenes, or all
  organometallics in CoordRep v1.
- No raw CSD coordinates are exported.
- Grammar extensions are prototypes (v0) and subject to future refinement.

## Part A: Multinuclear η1 / Bridging-Ligand Prototype (CoordRep-Multi-v0)

**15 cases** covering dinuclear (Cu, Fe, Rh, Co, Mn, Pt, Zn,
Mo, Ni) and tri-/tetranuclear (Cr, Fe, Cu, Mn) complexes.

- Bridge types: μ2-Cl, μ2-O, μ2-OH, μ2-OAc, μ2-OPh, μ2-dppm, μ3-O, μ3-OH, μ3-OMe
- Metal counts: 2–4
- All records parse-valid: **True**
- All records roundtrip-valid: **True**

See `multinuclear_examples_for_si.md` for full records.

## Part B: π / Haptic Coordination-Site Prototype (CoordRep-Haptic-v0)

**8 cases** covering:

- η5-Cp (ferrocene, Cp2TiCl2, CpMn(CO)3)
- η2-alkene (Zeise's salt, Ni(C2H4)(PPh3)2)
- η3-allyl (Pd(allyl)Cl2)
- η6-arene (Cr(benzene)(CO)3)
- Mixed hapticity (Ru(Cp*)(cod)Cl)

- All records parse-valid: **True**
- All records roundtrip-valid: **True**

See `haptic_examples_for_si.md` for full records.

## Validation Summary

| Property | Multinuclear | Haptic |
|---|---|---|
| Cases | 15 | 8 |
| All parse-valid | True | True |
| All roundtrip-valid | True | True |

## Files

- `multinuclear_case_index.csv` — case index with validation flags
- `multinuclear_records.jsonl` — full prototype records (line-broken + unwrapped)
- `multinuclear_validation_report.csv` — per-case validation checks
- `multinuclear_examples_for_si.md` — formatted examples for SI text
- `haptic_case_index.csv` — case index with validation flags
- `haptic_records.jsonl` — full prototype records
- `haptic_validation_report.csv` — per-case validation checks
- `haptic_examples_for_si.md` — formatted examples for SI text
- `coordrep_extension_prototype_summary.json` — machine-readable summary
- `coordrep_extension_prototypes_for_si.md` — this file
- `README.md` — project README
