# Protocol v1 retired before outcome analysis

Protocol v1 and its exact code/test sources are retained only as an audit
record. The first full-collection process was stopped after 1,000 of 15,906
CIFs, before aggregation or inspection of validation outcomes, because an
independent methods audit identified:

1. validation refcodes were not fully disjoint from development CSD families;
2. distance-based QC was incorrectly allowed to determine inclusion for sites
   lacking a complete explicit image;
3. site-level rates did not yet use refcode-clustered reporting;
4. the public-output gate did not prove a complete, internally consistent run.

The only inspected full-input quantity was exact-refcode presence in the local
CSD snapshot (14,719 present; 1,187 absent), which had already been established
before protocol v1. No v1 application outcome is used in the manuscript.

Exact retired sources:

- `run_csd_mof_periodic_application_v1_retired.py`
  SHA-256 `82DEE66E95A163BF4A5C2CF2E5ACB3B3F61352D8D83473F1BAA369989A4022B6`
- `test_csd_mof_periodic_application_v1_retired.py`
  SHA-256 `E464E4C2782B28C660F60CF3C15D99A6814FAD5414C25F4F8D0AC5558C412400`

Protocol v2 supersedes v1 and is frozen before its full run.
