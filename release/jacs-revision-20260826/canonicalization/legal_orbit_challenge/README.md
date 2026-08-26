# Legal residual-orbit challenge

This directory contains the locked 5,173/6,427 legal-orbit screen and the
1,000-record, 100-reindexing challenge used in Figure 2 and Supplementary
Figure S2. Signature-only sorting reproduced 64,834/100,000 variants and fully
collapsed 380/1,000 records. Full legal-orbit enumeration followed by
minimum-whole-record selection reproduced 100,000/100,000 variants and fully
collapsed 1,000/1,000 records, with stable candidate counts, successful grammar
round trips, and no fatal errors.

The records are fixed tmQMg/PBE CN = 4–6 snapshots and are not an April 2025
CSD prevalence sample. `run_orbit_challenge.py` is the exact archived harness;
its internal task-local status labels record the pre-integration execution
context. The tested implementation is now frozen as CoordRep 1.1.2rc3.
