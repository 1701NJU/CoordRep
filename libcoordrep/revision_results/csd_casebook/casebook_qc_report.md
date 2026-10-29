# CSD Casebook QC Report

## QC1: Case A — ACUWOK (Refcode-Family Identity Ladder)

**Verdict: ✓ PASS — suitable for main text**

| Check | Result |
|-------|--------|
| Family size | 12 entries (ACUWOK … ACUWOK11) |
| All same family | ✓ all share family prefix ACUWOK |
| L0 unique count | **12/12** — every entry has a distinct L0 |
| L0 pairwise matches | **0/66** — no two L0s are identical |
| L1 unique count | 2 (Td/SP.good.D3 and SP/Td.ideal.D2) |
| L2 unique count | 2 |
| L3 unique count | **1** — all 12 share the same L3 |
| Metal | Pt |
| CN | 4 |
| Shapes | {Td, SP} — polymorphic square-planar / tetrahedral |
| is_boundary | all False |

**L3 common key:**
`Pt|CN4|Cl;Cl;[H]C1([H])SC([H])([H])C([H])([H])SC([H])([H])C([H])([H])SC1([H])[H]`

**Why L0 differs:** CShM values differ across 12 conformers (geometric detail: V values vary).
**Why L3 links:** Metal + CN + canonical ligand set is invariant → same chemical identity.

---

## QC2: Case B — AFOSIA (Boundary Geometry)

**Verdict: ✓ PASS — suitable for main text**

| Check | Result |
|-------|--------|
| Refcode | AFOSIA |
| Metal | Cr, CN = 6 |
| is_boundary | **True** |
| ShapeBest | TPr |
| CShM(TPr) = V1 | **4.49** |
| CShM(Oh) = V2 | **4.48** |
| Δ = |V1 − V2| | **0.01** ← extreme boundary |
| Class | dist |
| L1 token | `Cr|?|TPr/Oh_boundary.dist.D0|T3|…` |

**Boundary rule verified:**
`is_boundary = (delta < boundary_thresh)` where `boundary_thresh = 1.0` (from `shape_binning.py` line 57).
Here delta = 0.01 ≪ 1.0, so `is_boundary = True` ✓

**Why this matters:** CShM(TPr) and CShM(Oh) are nearly identical. Any discrete shape label assignment is arbitrary — CoordRep's continuous boundary annotation captures this ambiguity.

---

## QC3: Case C — ADIZUI (Tool B Repair)

**Verdict: ✗ FAIL as originally framed — ADIZUI must NOT go in main text**

### Field-level comparison (ADIZUI, truncation corruption)

| Field | Clean | Corrupted | Rule | MLM |
|-------|-------|-----------|------|-----|
| metal | Co | Co | Co ✓ | ? ✗ |
| cn | 6 | 6 | 6 ✓ | ? ✗ |
| ox | +3 | +3 | +3 ✓ | ? ✗ |
| shape | Oh | Oh | Oh ✓ | ? ✗ |
| ligands | 4 ligs | 4 ligs | 4 ligs ✗ | 4 ligs ✗ |
| stereo | 2 rels | 2 rels | 2 rels ✓ | 2 rels ✓ |
| **Total** | **(ref)** | **5/6** | **5/6** | **1/6** |

MLM reformats the header (`[Metal:Co|ox:+3|...]` → `[Co];ox=+3;...`), destroying field parsability. Over 15 candidates tested, **MLM never outperformed rule repair at field level.**

### Systematic survey (500 entries × 2 corruption types)

| Method | Typical field score |
|--------|-------------------|
| No repair (corrupted) | 5/6 (header preserved) |
| Rule repair | 5–6/6 |
| Char bigram | 5–6/6 (identical to rule) |
| CoordRep-MLM | 1–2/6 (header destroyed) |

### Reframing recommendation

Case C should be reframed from "MLM beats rules" to:

> **CoordRep's explicit grammar enables highly effective rule-based repair.**
> The structured field–delimiter–bracket grammar allows deterministic rules to
> recover 100% syntactic validity on CSD-derived corruption, while preserving
> 5/6 semantic fields. This repair is a direct consequence of CoordRep's
> design — raw SMILES or unstructured text would not support such grammar-
> aware recovery.

**Proposed main-text candidate:** BEMVIY (Co/CN6, missing_bracket)
- Bracket removed from ligand L1 SMILES position 111: `C([H])` → `C([H]`
- Corrupted string: **invalid** (bracket imbalance)
- Rule repair: restores validity using CoordRep bracket grammar
- Preserves 5/6 fields (metal, CN, ox, shape, stereo all correct; ligand L1 has bracket at wrong position)

---

## QC4: Case D — CIJWUO (Stereo Semantic Consistency)

**Verdict: ✓ PASS — suitable for main text**

### Stereo decoy construction

| Item | Value |
|------|-------|
| Refcode | CIJWUO |
| Metal | Cu, CN = 5 |
| Real stereo | `{trans:L1:N:1--L2:N:1}` |
| Decoy stereo | `{cis:L1:N:1--L2:N:1}` |
| Decoy grammar-valid | **True** (passes `is_valid_coordrep`) |
| Character differences | 197 (many due to different tokenization) |

### Ranker scoring — all 5 candidates

| Refcode | Metal/CN | Flip | Score(real) | Score(decoy) | Margin | Correct |
|---------|----------|------|-------------|-------------|--------|---------|
| CIJWUO | Cu/5 | trans→cis | −2.479 | −3.164 | **+0.685** | ✓ |
| BATTIZ | Pd/4 | trans→cis | −4.653 | −4.979 | **+0.326** | ✓ |
| CIRRUS | Pt/4 | trans→cis | +3.932 | +3.592 | **+0.340** | ✓ |
| CEHZIZ | Mn/6 | trans→cis | −0.717 | −0.813 | **+0.097** | ✓ |
| BECDUK | Cu/2 | trans→cis | −10.531 | −11.170 | **+0.639** | ✓ |

**5/5 correct** — CoordRep-Ranker consistently scores real stereo higher than flipped decoy.

**Key point:** The stereo decoy is grammar-valid (passes all structural checks). Only a learned model that understands field-level semantic compatibility can distinguish real from decoy. This is a semantic-level test, not a syntax test.
