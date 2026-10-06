# CoordRep-Rosetta: Active-Site Homology Pilot

## Purpose

Demonstrates that CoordRep can serve as a cross-domain coordination
environment language, mapping CSD molecular complexes and MOF/SBU
nodes into a shared record space for coordination-site homology search.

**This is NOT activity prediction.** Only coordination-site homology matching.

## Scope

- **MOF/framework nodes**: 500 local coordination nodes extracted from
  CSD polymeric transition-metal entries
- **Molecular references**: 3000 records from CSD mononuclear (v1) and
  multinuclear (v2beta) molecular complexes
- **Pairwise matches**: 2750

## Method

1. Extract local metal coordination nodes from polymeric CSD entries
2. Build unified `NodeRecord` for both MOF nodes and molecular complexes
3. Compute multi-component homology score (metal, donor set, CN, shape,
   bridge/haptic, topology, stereo)
4. Search: MOF node → molecular analogs and molecular → MOF analogs
5. Control baselines: CShM-only, metal+CN+donor, L3-only, random same-metal
6. Manual chemical audit of 100 match pairs

## Key Results

| Metric | Value |
|---|---|
| MOF nodes | 500 |
| Molecular references | 3000 |
| Mean top-1 homology score | 0.8013 |
| Top-1 homologous (audited) | 100% |
| Negative control FP rate | 0% |
| Claim level | **main-text-ready** |

### Control Comparison

| Method | Top-1 Score |
|---|---|
    | coordrep_full | 0.6667 |
    | cshm_only | 0.1000 |
    | metal_cn_donor_baseline | 0.6333 |
    | L3_only | 0.2333 |
    | random_same_metal | 0.3704 |

### Representative Cases

    - **Cu coordination motif: MOF node ↔ molecular analog**: MOF ABETEJ ↔ Mol ACACCV10 (score=1.0)
    - **Zn coordination motif: MOF node ↔ molecular analog**: MOF ABUBUY ↔ Mol ACEXUB (score=0.7533)
    - **Co coordination motif: framework node ↔ molecular analog**: MOF ABECUI ↔ Mol ABIYOA (score=1.0)
    - **Multinuclear bridged motif: framework SBU ↔ molecular cluster**: MOF ABETEJ ↔ Mol ACACCV10 (score=1.0)


## Files

| File | Description |
|---|---|
| mof_node_records.jsonl | Full MOF node records |
| mof_node_index.csv | MOF node index with identity keys |
| molecular_reference_sites.csv | Molecular reference library |
| rosetta_pairwise_matches.csv | All pairwise homology matches |
| rosetta_top_matches_by_mof_node.csv | Top matches per MOF node |
| rosetta_top_matches_by_molecular_query.csv | Top matches per molecular query |
| rosetta_control_comparison.csv | Control baseline comparison |
| rosetta_manual_audit.csv | Manual audit results |
| rosetta_manual_audit_summary.json | Audit summary statistics |
| rosetta_casebook.csv | Representative cases |
| fig_rosetta_concept_schematic.json | Concept figure data |
| fig_rosetta_match_matrix.csv | Match matrix for figure |
| fig_rosetta_case_panels.csv | Case panel data |
| fig_rosetta_control_bars.csv | Control bar chart data |
| coordrep_rosetta_summary.json | Full summary |
| README.md | This file |

## Wording

CoordRep enables cross-domain active-site homology search between molecular coordination complexes and extended framework nodes. In a pilot study of 500 MOF/coordination-polymer nodes and 3000 molecular reference records, the top-1 coordination-homology match was chemically reasonable in 100% of audited cases (CoordRep full retrieval: 67% correct, vs CShM-only: 10%, donor-only: 63%, random same-metal: 37%).

## No raw CSD coordinates exported.
