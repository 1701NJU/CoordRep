# Box 1 — CoordRep-v2 beta Extension Examples: Selection Notes

These two examples are clearly labeled as **CoordRep-v2 beta extension
examples** and are distinct from production v1 records (Examples 1–3).


## Example 4: Molecular Multinuclear Record Graph

**Example ID:** example_4
**Refcode:** ETOREN
**Record type:** molecular_multinuclear_graph

### Why selected

ETOREN is Cu₂(μ-Cl)₂Cl₂(PR₃)₂ — a textbook di-μ-chloro bridged
dicopper(I) complex. It was selected because:

- **Finite molecular complex** (not polymeric, not a MOF/coordination polymer).
- **Two metal centers** — simplest possible multinuclear graph.
- **Two clearly assigned μ₂-Cl bridges** — the most recognizable
  bridging mode in introductory inorganic chemistry.
- **Mixed ligand environment** (bridging Cl, terminal Cl, terminal P)
  shows that the record graph correctly handles both bridge and
  terminal coordination in the same complex.
- **Symmetric Cu₂** pair makes canonical metal ordering straightforward.

### What v2-beta feature it demonstrates

CoordRep-v2 beta represents multinuclear structures as **metal-centered
record graphs** rather than forcing all donors into one mononuclear
coordination sphere. Key features shown:

1. **Independent local spheres per metal**: M1={S1,S2,S3,S5}, M2={S1,S3,S4,S6}.
   Each metal's CN is derived from its own local site count.
2. **Bridging ligand appears once**: S1 and S3 are single site objects
   with mu=2 and target_metals=[M1,M2]. They are NOT duplicated as
   two separate terminal ligands.
3. **MetalGraph edge**: M1--M2 relation is "bridged" with MM_bond="no",
   correctly noting that the Cu–Cu connection is via bridging Cl, not a
   direct metal–metal bond.
4. **Graph-level identity keys** (L0–L3) + per-metal LocalIDs enable
   hierarchical identity matching.

### Chemical validation

- 2 Cu centers, each with CN_site=4 (two bridging Cl + one terminal Cl
  + one phosphine)
- Bridge sites S1 and S3: mu=2, target_metals=[M1,M2] — consistent with
  mu₂ bridging
- No metal–metal bond claimed (correct for d¹⁰ Cu(I) at typical
  Cu–Cu distances in di-μ-chloro complexes)

### Record validation

All 14 validation gates passed:
- parse_valid, roundtrip_valid, no_placeholder_tokens
- metal_order_invariant, atom_order_invariant, ligand_order_invariant
- site_atom_order_invariant, bridge_consistency_valid
- eta_mu_consistency_valid, local_sphere_consistency_valid
- no_duplicate_ligand_for_same_bridge
- no_raw_coordinates_exported, identity_keys_present
- human_readable_si_example

### Known limitations

- Oxidation state is "unk" — the v2-beta pipeline does not yet assign
  oxidation states for multinuclear complexes.
- local_shape_best is empty — CShM is not yet implemented for
  multinuclear metals in v2-beta.
- Ligand SMILES are atom-level fragments (e.g. [Cl], [P]) rather than
  full molecular SMILES of the entire phosphine ligand.


---


## Example 5: Haptic/π Coordination-Site Object

**Example ID:** example_5
**Refcode:** FEROCE01
**Record type:** haptic_pi_site_object

### Why selected

FEROCE01 is ferrocene, Fe(η⁵-Cp)₂ — the most iconic organometallic
sandwich complex and arguably the most widely recognized haptic system
in chemistry. It was selected because:

- **Mononuclear** — simplest possible haptic record.
- **Two η⁵-Cp sites** — the canonical haptic ligand example.
- **No ambiguity** in hapticity assignment (all five carbons of each Cp
  ring are π-bonded to Fe).
- **D₅d/D₅h symmetry** makes canonical ordering unambiguous.
- Universally familiar to any inorganic/organometallic chemist.

### What v2-beta feature it demonstrates

CoordRep-v2 beta represents haptic/π ligands as **coordination-site
objects** rather than decomposing them into individual η¹ donor atoms.
Key features shown:

1. **Site objects**: S1 and S2 are type="pi_fragment" with eta=5,
   containing five donor atoms each. The Cp ring is ONE site, not
   five independent C donors.
2. **CN_site vs eta_sum distinction**: CN_site=2 (two coordination
   sites), eta_sum=10 (total atom contacts 5+5). This avoids the
   ambiguity of "CN=10" in traditional atom-resolved notation.
3. **Contiguous atom sets**: C1-C5 and C1B-C5B are contiguous
   five-membered rings — the BFS fragment grouping correctly identifies
   connected π systems.
4. **Haptic identity keys** (L0_HapticState through L3_Connectivity)
   provide hierarchical identity matching for organometallic complexes.

### Chemical validation

- Fe center with CN_site=2 (two η⁵-Cp sites)
- eta_sum=10 matches: 5 atoms (S1) + 5 atoms (S2) = 10
- Both sites are pi_fragment type with eta=5, mode=eta5
- Atom sets {C1,C2,C3,C4,C5} and {C1B,C2B,C3B,C4B,C5B} are
  contiguous five-membered carbon rings

### Record validation

All 14 validation gates passed:
- parse_valid, roundtrip_valid, no_placeholder_tokens
- metal_order_invariant (trivially true: single metal)
- atom_order_invariant, ligand_order_invariant, site_atom_order_invariant
- bridge_consistency_valid (N/A: no bridges)
- eta_mu_consistency_valid (eta_sum=10 consistent)
- local_sphere_consistency_valid (cn_site=2 = len(sites)=2)
- no_duplicate_ligand_for_same_bridge
- no_raw_coordinates_exported, identity_keys_present
- human_readable_si_example

### Known limitations

- Oxidation state is "unk" — the v2-beta pipeline does not yet assign
  oxidation states for haptic complexes (Fe²⁺ d⁶ would be expected).
- local_shape_best is empty — no site-level geometry label is assigned.
- Ligand SMILES are disconnected atom fragments [C].[C].[C].[C].[C]
  rather than a cyclopentadienyl ring SMILES (c1cccc1 or similar).
- centroid_label is empty — centroid coordinates are not stored (no raw
  coordinates exported by design).


---


## Data Provenance

- Example 4 (ETOREN) is from the **multinuclear_molecular_subset_audit**
  in `revision_results/full_csd_v2beta_postfix_audit/`, which audited
  4,579 molecular multinuclear CSD entries with 100% conversion,
  roundtrip, and invariance pass rates.

- Example 5 (FEROCE01) is from the **haptic curated case set** in
  `revision_results/coordrep_v2_beta_extension/`, case haptic_001,
  which passed all validation gates.

Both records were generated by the production pipeline
(`coordrep.v2beta.csd_v2beta_adapter`) and validated by the production
validator (`coordrep.v2beta.validate`). No manual editing was applied.
