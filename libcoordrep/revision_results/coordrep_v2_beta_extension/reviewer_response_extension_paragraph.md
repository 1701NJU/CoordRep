## Response to Reviewer — Multinuclear and Haptic Extensions

We thank the reviewer for raising the important question of coverage
beyond mononuclear η1 complexes.  We have addressed this concern as
follows:

1. **Denominator-corrected scope waterfall.**  The headline figure of
   8.83% refers to all 1,413,222 CSD entries; among the 615,498
   transition-metal candidates, CoordRep v1 covers 20.3%, and among
   the 126,197 entries within its intended scope (mononuclear, η1,
   CN 2–6, no disorder), 98.9% are successfully converted.

2. **Multinuclear systems** are the largest outside-scope class
   (346,468 entries, 56.3% of TM candidates).  We implemented
   CoordRep-Multi-v2beta as a metal-centered record graph that
   encodes local spheres, bridging ligands, and metal–metal
   relations in a single canonical record.

3. **Haptic/π systems** (51,867 entries, 8.4% of TM candidates)
   are addressed by CoordRep-Haptic-v2beta, which introduces
   coordination-site objects (atom_set / pi_fragment / centroid)
   with canonical atom ordering and site-level geometry.

4. **Validation.**  On 55 curated multinuclear cases and 33 curated
   haptic cases, all records pass parse, roundtrip, and invariance
   checks (14 validation gates each).  Projected random-CSD-subset
   pass rates are ~72% (multinuclear) and ~68% (haptic).

5. **Scope.**  We do not claim full production support for these
   extensions in the present v1 benchmark.  The main text states:
   "Prototype SI tests further show that the CoordRep grammar can
   be extended beyond v1 through metal-centered record graphs for
   multinuclear systems and coordination-site objects for haptic
   ligands.  These tests are not included in the present v1
   benchmark."
