# CoordRep v1 Scope

## Supported (v1)

| Feature | Detail |
|---------|--------|
| Nuclearity | Mononuclear only |
| Donor type | η1 atom-resolved donors |
| Denticity | Monodentate and multidentate chelating ligands |
| CN range | 2–14 |
| Stereo tokens | cis/trans, fac/mer, donor-relation pairs |
| Identity | CoordRep-ID L0–L3 hierarchy |
| Geometry | Continuous CShM with boundary annotation |
| Metals | All d-block transition metals (Sc–Hg) |

## Out of v1 Scope

| Feature | Reason | Future |
|---------|--------|--------|
| Multinuclear complexes | Requires multi-center graph extension | CoordRep-Multi |
| MOFs / coordination polymers | Periodic boundary conditions | CoordRep-Periodic |
| Polyoxometalates | Cluster-level representation | CoordRep-Cluster |
| ηn haptic organometallics | Requires face/edge coordination model | CoordRep-Haptic |
| Positional disorder / partial occupancy | Ambiguous atom assignment | Future validator |
| Weak interactions beyond threshold | M–L distance cutoff is hard boundary | Soft threshold v2 |

## Quantified Exclusion (CSD External Audit)

From 200,000 CSD entries scanned:
- 17,038 retained (8.5%)
- Top rejections:
  - No transition metal: 100,863 (50%)
  - Multinuclear: 47,939 (24%)
  - No 3D structure: 13,564 (7%)
  - Disorder: 8,435 (4%)
  - Hapticity: 7,325 (4%)

The v1 scope boundary is by design, not a limitation. Within this scope, CoordRep provides capabilities that neither raw SMILES, graph hashes, nor 3D GNNs can replicate.
