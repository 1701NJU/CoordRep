# CoordRep-v2 Beta — Methods

## Record Construction

All v2-beta records are generated programmatically by the serializer
(`coordrep.v2beta.serialize`).  No records are hand-written.

## Canonicalization

Metal-node ordering uses Weisfeiler-Lehman-like iterative refinement
on the metal graph with chemical signatures (element, oxidation, CN,
local site types, bridge degree).  Site ordering uses canonical
site_signature tuples.  Ligand ordering uses SMILES + site-list hash.

## Validation Gates (14 checks per record)

1. parse_valid
2. roundtrip_valid
3. no_placeholder_tokens
4. metal_order_invariant
5. atom_order_invariant
6. ligand_order_invariant
7. site_atom_order_invariant
8. bridge_consistency_valid
9. eta_mu_consistency_valid
10. local_sphere_consistency_valid
11. no_duplicate_ligand_for_same_bridge
12. no_raw_coordinates_exported
13. identity_keys_present
14. human_readable_si_example

## Tier 2 Random Audit

Random CSD audit (n=500 multinuclear, n=300 haptic) is projected from
curated validation rates and known category heterogeneity.  Actual
execution requires CSD Python API access (not available in this environment).

## Coverage Estimate

If v2-beta achieves 72% on multinuclear and 68% on haptic random subsets,
the estimated additional coverage is ~284,725
entries (+228.1% over the current 124,837 v1 records).
