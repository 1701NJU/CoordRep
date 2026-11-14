# Supplementary Information

**CoordRep: A Canonical Text Representation for Coordination Complexes Enabling Masked-Field Learning of Chemical Grammar**

---

## Supplementary Methods 1. Data sources, scope definition, and filtering rules

### SM1.1 Data sources and version identifiers

CoordRep draws structures from three databases:

| Source | Version | Records available | Coverage purpose |
|--------|---------|-------------------|------------------|
| Cambridge Structural Database (CSD) | 2024.3 (1,413,222 entries) | Single-crystal XRD | Primary high-quality source |
| Crystallography Open Database (COD) | 2024-01 snapshot | Single-crystal XRD | Open-access complement |
| tmQM | v1 (86,665 complexes) | DFT-optimized | Computed geometries for training diversity |

### SM1.2 Scope definition

CoordRep-v1 targets mononuclear, η¹-only transition-metal coordination complexes:

**Supported features:**
- Mononuclear transition-metal complexes (3d/4d/5d)
- Atom-resolved η¹ donors (N, O, S, P, C, Cl, Br, I, F, Se, Te, As)
- Monodentate, bidentate, and higher-denticity η¹ chelating ligands
- Coordination number 2–6 (primary), 7–14 (shape vectors only)
- cis/trans, fac/mer, and donor-relation stereochemical tokens
- CoordRep-ID hierarchy L0–L3

**Explicitly excluded:**
- Multinuclear complexes (cluster, polynuclear)
- MOFs and coordination polymers
- Polyoxometalates
- η^n organometallics (ferrocene, arene complexes)
- Positional disorder / partial occupancy entries
- Ambiguous coordination spheres (uncertain CN)

### SM1.3 Waterfall filtering procedure

Starting from 1,413,222 CSD entries, we apply sequential hard filters. Each stage is independently verifiable:

**Supplementary Table S1 | Full-CSD waterfall filter stages.**

| Stage | Criterion | Retained | Removed |
|-------|-----------|----------|---------|
| 0 | Total CSD entries scanned | 1,413,222 | — |
| 1 | Has 3D structure | 1,335,141 | 78,081 |
| 2 | Contains transition metal | 615,498 | 719,643 |
| 3 | Mononuclear | 269,030 | 346,468 |
| 4 | η¹ only (no hapticity) | 217,163 | 51,867 |
| 5 | CN ∈ {2,3,4,5,6} | 209,826 | 7,337 |
| 6 | No crystallographic disorder | 145,728 | 64,098 |
| 7 | Non-polymeric | 144,271 | 1,457 |
| 8 | Valid RawMolecule construction | 126,197 | 18,074 |
| 9 | Valid SMILES for all ligands | 126,197 | 0 |
| 10 | Valid CoordRep serialization | 124,837 | 1,360 |

Final retention rate: 8.83% (124,837 / 1,413,222).

Of the 124,837 retained entries, 17,791 (14.3%) are flagged as geometry-boundary records (ΔCShM < 1.0 between top-2 reference shapes).

### SM1.4 Top rejection reasons

| Reason | Count |
|--------|-------|
| No transition metal | 719,643 |
| Multinuclear | 346,468 |
| No 3D coordinates | 78,081 |
| Crystallographic disorder | 64,098 |
| Hapticity (η^n) | 51,867 |
| RawMolecule construction failure | 18,074 |
| CN outside 2–6 | 7,337 |
| Polymeric | 1,457 |
| Invalid CoordRep | 1,360 |

### SM1.5 Stratified data split

The pretraining and evaluation corpus is split 80/10/10 (train/val/test) with stratification by metal element, coordination number, and data source. Split statistics are verified to preserve coverage:

**Supplementary Table S2** | Stratified split statistics by metal series (3d/4d/5d). [Keep existing content.]

**Supplementary Table S3** | Stratified split statistics by coordination number. [Keep existing content.]

---

## Supplementary Methods 2. CoordRep grammar, tokenization, and complete serializer outputs

### SM2.1 Token grammar specification

A complete CoordRep-v1 string follows the grammar:

```
<record>   ::= <metal> <shape> <stereo>* <ligands> <trailing>?
<metal>    ::= '[Metal:' ELEMENT '|ox:' OX '|d:' DCOUNT '|CN:' INT ']'
<shape>    ::= '<ShapeBest:' SHAPE_LABEL '|Class:' CLASS '|Delta:' INT '|V:' FLOAT ',' FLOAT '>'
<stereo>   ::= '{trans:' DONOR_REF '--' DONOR_REF '}' | '{cis:' ... '}' | '{fm:' FAC_MER '}'
<ligands>  ::= ('|L' INT '=' SMILES)+ '|'
<donor_ref>::= 'L' INT ':' ELEMENT ':' INT
```

Where:
- `ELEMENT` ∈ {H, He, ..., Og} (periodic table symbols)
- `OX` ∈ {+1, +2, ..., +8, 0, -1, ...} (formal oxidation state)
- `DCOUNT` ∈ {d0, d1, ..., d10} (d-electron count)
- `SHAPE_LABEL` ∈ {SP, Td, Oh, TPr, SPY, TBPY, L, TY, TP, ...}
- `CLASS` ∈ {ideal, dist, boundary} (geometry classification)
- `FAC_MER` ∈ {fac, mer}

### SM2.2 Token types

**Supplementary Table S4 | Token types in the CoordRep vocabulary.**

| Token class | Example | Count in vocab | Role |
|-------------|---------|----------------|------|
| Metal header | `[Metal:Pt\|ox:+2\|d:d8\|CN:4]` | ~480 composite | Encode metal state |
| Shape token | `<ShapeBest:SP\|Class:ideal\|Delta:2\|V:3.03,0.18>` | continuous | Geometry encoding |
| Stereo constraint | `{trans:L1:N:1--L2:N:1}` | variable | Stereochemical relations |
| fac/mer | `{fm:fac}`, `{fm:mer}` | 2 | Isomer distinction |
| Ligand SMILES | `\|L1=[H]N([H])[H]\|` | open vocabulary | Ligand structure |
| Donor reference | `L1:N:1` | derived | Donor site identifier |

### SM2.3 Composite vs factorized tokenization

In the default CoordRep-v1 tokenizer, the metal header is a **composite token**: `[Metal:Fe|ox:+2|d:d6|CN:6]` is a single vocabulary entry. This design choice enables:
- Efficient single-token prediction of the full metal state
- Natural co-occurrence statistics during training

The **factorized** variant splits this into independent tokens: `[Metal:Fe]`, `[ox:+2]`, `[d:d6]`, `[CN:6]`. This enables independent recovery assessment (see Supplementary Methods 6).

### SM2.4 Complete serializer outputs (Box 1 records)

The following records are direct outputs of the CoordRep serializer applied to CSD crystal structures. They are reproduced here in full, unwrapped form for machine readability. Line-broken display versions appear in the main text Box 1.

#### Example A: cis-[Pt(MeNH₂)₂I₂] (CSD: IYIPOV)

```
[Metal:Pt|ox:+2|d:d8|CN:4]<ShapeBest:SP|Class:ideal|Delta:2|V:3.25,0.47>{trans:L1:N:1--L3:I:1}{trans:L2:N:1--L4:I:1}|L1=[H]N([H])C([H])([H])[H]|L2=[H]N([H])C([H])([H])[H]|L3=I|L4=I|
```

Identity keys:
- L0 (StateKey): full string above (hash: `ed29845b`)
- L1 (ShapeID): `Pt|+2|SP/Td.ideal.D2|T2|I;I;[H]N([H])C([H])([H])[H];[H]N([H])...`
- L2 (TopoID): `Pt|CN4|SP|I;I;[H]N([H])C([H])([H])[H];[H]N([H])C([H])([H])[H]`
- L3 (ConnID): `Pt|CN4|I;I;[H]N([H])C([H])([H])[H];[H]N([H])C([H])([H])[H]`

Validation: parse_valid=True, roundtrip_valid=True.

#### Example B: [Co(en)₃]³⁺ (CSD: JOXSOI)

```
[Metal:Co|ox:+3|d:d6|CN:6]<ShapeBest:Oh|Class:dist|Delta:2|V:3.64,7.04>{trans:L1:N:1--L2:N:2}{trans:L1:N:2--L3:N:1}{trans:L2:N:1--L3:N:2}|L1=C2D8N2|L2=C2D8N2|L3=C2D8N2|
```

Identity keys:
- L0 hash: `71024466`
- L1: `Co|+3|Oh/TPr.dist.D2|T3|C2D8N2;C2D8N2;C2D8N2`
- L2: `Co|CN6|Oh|C2D8N2;C2D8N2;C2D8N2`
- L3: `Co|CN6|C2D8N2;C2D8N2;C2D8N2`

Key observation: Each ethylenediamine (en) ligand contributes **two** donor sites (L1:N:1, L1:N:2). The three trans pairs cross different ligands, encoding the Δ-/Λ-independent octahedral topology.

Validation: parse_valid=True, roundtrip_valid=True.

#### Example C (fac): fac-[Co(dien)(CN)₃] (CSD: IDARAG)

```
[Metal:Co|ox:+3|d:d6|CN:6]<ShapeBest:Oh|Class:dist|Delta:1|V:5.76,6.98>{trans:L1:N:1--L2:C:1}{trans:L1:N:2--L3:C:1}{trans:L1:N:3--L4:C:1}{fm:fac}|L1=[H]N([H])C([H])([H])C([H])([H])N([H])C([H])([H])C([H])([H])C([H])([H])N([H])[H]|L2=C#N|L3=C#N|L4=C#N|
```

Identity keys:
- L0 hash: `83759f7d`
- L1: `Co|+3|Oh/TPr.dist.D1|T3_fac|C#N;C#N;C#N;[H]N([H])C([H])([H])...`
- L2: `Co|CN6|Oh|C#N;C#N;C#N;[H]N([H])C([H])([H])C([H])([H])N([H])C...`
- L3: `Co|CN6|C#N;C#N;C#N;[H]N([H])C([H])([H])C([H])([H])N([H])C([H])([H])C([H])([H])C([H])([H])N([H])[H]`

The `{fm:fac}` token explicitly encodes that the three CN⁻ ligands occupy a facial arrangement. All three N-donors from dien are each trans to one CN.

#### Example D (mer): mer-[Co(dien)(CN)₃] (CSD: IDAREK)

```
[Metal:Co|ox:+3|d:d6|CN:6]<ShapeBest:TPr|Class:dist|Delta:0|V:5.67,5.34>{trans:L1:N:1--L1:N:2}{trans:L1:N:3--L2:C:1}{trans:L3:C:1--L4:C:1}{fm:mer}|L1=[H]N([H])C([H])([H])C([H])([H])N([H])C([H])([H])C([H])([H])C([H])([H])N([H])[H]|L2=C#N|L3=C#N|L4=C#N|
```

Identity keys:
- L0 hash: `6a1a2c27`
- L1: `Co|+3|TPr/Oh_boundary.dist.D0|T3_mer|...`
- L3: `Co|CN6|C#N;C#N;C#N;[H]N([H])C([H])([H])C([H])([H])N([H])C([H])([H])C([H])([H])C([H])([H])N([H])[H]`

Key contrast with Example C: **L3 is identical** (same connectivity), but L0, L1, and stereo tokens differ. The `{fm:mer}` token and the trans pattern `{trans:L1:N:1--L1:N:2}` (two N-donors from the same ligand trans to each other) distinguish the mer isomer. Additionally, the ShapeBest is TPr with Delta:0 (boundary), reflecting the geometric distortion that accompanies meridional coordination.

Validation: parse_valid=True, roundtrip_valid=True.

---

## Supplementary Methods 3. Canonicalization and CoordRep-ID construction

### SM3.1 Definition of determinism

A CoordRep serialization is **deterministic** if the same coordination complex produces a bit-identical output string regardless of:
- Input atom ordering
- Coordinate frame (rotation, translation)
- Traversal order of the molecular graph
- Minor numerical noise (within rounding tolerance)

### SM3.2 Deterministic pipeline overview

The canonicalization pipeline proceeds as:

1. **Metal identification** → unique metal center (mononuclear constraint)
2. **Coordination sphere detection** → distance + bond-order criterion
3. **Graph segmentation** → ligand extraction via metal-deletion in molecular graph
4. **Canonical ligand ordering** → deterministic sort by:
   - (a) RDKit canonical SMILES
   - (b) Donor element alphabetical
   - (c) Donor index within ligand (graph distance from attachment)
5. **Shape vector** → permutation-invariant CShM (Hungarian assignment)
6. **Constraint canonicalization** → trans-pairs and fac/mer sorted by ligand/donor indices
7. **String assembly** → concatenate in fixed field order

### SM3.3 Ligand segmentation

Ligand extraction is performed by:
1. Remove the metal atom from the molecular graph
2. Identify connected components → each component is one ligand
3. For each ligand, identify donor atoms = atoms that were bonded to the metal
4. Generate canonical SMILES for each ligand fragment (RDKit Chem.MolToSmiles with canonical=True)

### SM3.4 Canonical ligand ordering with deterministic tie-breakers

Ligands are sorted by a lexicographic key:
```
sort_key(L) = (canonical_SMILES, donor_element_tuple, n_donors, denticity)
```

When two ligands have identical SMILES and donor elements (e.g., two NH₃ in cisplatin), disambiguation uses the stereochemical constraint tokens: the ligand participating in the first (lowest-index) trans-pair receives the lower index.

### SM3.5 CoordRep-ID hierarchy construction

The four-level identity hierarchy provides multi-resolution complex identification:

| Level | Name | Content | Distinguishes |
|-------|------|---------|---------------|
| L0 | StateKey | Full canonical CoordRep string | Everything (unique state) |
| L1 | ShapeID | Metal + ox + shape_bin + topology + stereo_sig + truncated ligands | Stereoisomers within same topology |
| L2 | TopoID | Metal + CN + shape_best + sorted full ligand SMILES | Same-connectivity different geometry |
| L3 | ConnID | Metal + CN + sorted ligand molecular formulas | Connectivity family (ignores geometry and SMILES detail) |

**Hash construction:**
- L0: the full string itself (or SHA-256 hash for storage)
- L1: constructed by concatenating metal, oxidation, shape classification (e.g., `Oh/TPr.dist.D2`), topology signature (e.g., `T3_fac`), and truncated sorted ligand SMILES
- L2: `metal|CN{n}|shape_best|sorted_ligand_smiles_joined_by_semicolon`
- L3: `metal|CN{n}|sorted_ligand_formula_joined_by_semicolon`

### SM3.6 Empirical determinism verification

We verify determinism by:
1. **Permutation test**: For 1,000 complexes, randomly permute input atom order 10× each → check L0 identity. Result: 100.0% identical across all permutations.
2. **Rotation test**: Apply random 3D rotations → check L0 identity. Result: 100.0% identical.
3. **Noise test**: Add Gaussian noise (σ = 0.001 Å) → check CShM rounding stability. Result: 100.0% identical for CShM rounded to 2 decimal places.

---

## Supplementary Methods 4. Coordination Identity Challenge

### SM4.1 Motivation

We define the **Coordination Identity Challenge** to quantify which operations a coordination representation can and cannot support. Existing representations (SMILES, InChI, raw coordinates, CShM vectors) each have blind spots for coordination chemistry. The challenge provides a structured comparison.

### SM4.2 Operations tested

Six fundamental operations required for coordination chemistry informatics:

| Operation | Definition | Why it matters |
|-----------|-----------|----------------|
| Collapse invariance | Same complex → same identifier regardless of input format | Database deduplication |
| Stereo separation | Different stereoisomers → different identifiers | cis/trans, fac/mer discrimination |
| Family linking | Same connectivity, different geometry → linked as family | Structural surveys, polymorph tracking |
| Boundary detection | Ambiguous geometry (ΔCShM < 1) → flagged, not force-labelled | Honest uncertainty handling |
| Grammar validation | Syntactically invalid record → detected | Quality control |
| Multi-resolution ID | Query at different granularity (exact vs family) | Flexible retrieval |

### SM4.3 Representations evaluated

| Representation | Dimensions | Metal-aware? | Stereo tokens? |
|----------------|-----------|--------------|----------------|
| Canonical SMILES (multiset) | 1D string | Partial | No |
| InChI/InChIKey | 1D string | No (organic) | No |
| Raw 3D coordinates | 3N floats | Implicit | Implicit |
| CShM vector | 2–3 floats | No | No |
| CoordRep-v1 (full) | 1D string | Yes | Yes |
| CoordRep-ID (L0–L3) | 4-level hash | Yes | Yes |

### SM4.4 Challenge pair construction

For each operation, we construct **20 challenge pairs** from the retained CSD corpus:

- **Invariance**: 20 polymorph pairs sharing L3 + constraint signature (same complex, different crystal packing). Pool: 4,950 pairs.
- **Stereo separation**: 20 pairs with identical L3 ConnID but different cis/trans or fac/mer tokens. Pool: 745 pairs.
- **Family linking**: 20 pairs with same L3 but different L0/L1 (geometry state). Pool: 3,739 pairs.
- **Boundary detection**: 20 entries with ΔCShM < 1.0. Pool: 17,791 entries.

### SM4.5 Results matrix

**Supplementary Table S5 | Representation capability matrix.**

| Representation | Collapse | Stereo | Family | Boundary | Grammar | Multi-res |
|----------------|----------|--------|--------|----------|---------|-----------|
| Canonical SMILES | Partial | ✗ | ✗ | ✗ | ✗ | ✗ |
| InChI/InChIKey | Partial | ✗ | ✗ | ✗ | Partial | ✗ |
| Raw 3D coords | ✗ | Partial | ✗ | ✗ | ✗ | ✗ |
| CShM vector | ✗ | ✗ | ✗ | Partial | ✗ | ✗ |
| CoordRep full | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| CoordRep-ID L0–L3 | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |

---

## Supplementary Methods 5. Continuous geometry and CShM analyses

### SM5.1 Mathematical definition

The Continuous Shape Measure (CShM) quantifies deviation of an observed coordination geometry from an ideal reference polyhedron P₀:

$$S(P_0) = \min_{R, \sigma} \frac{\sum_{i=1}^{N} \|q_i - R \cdot p_{\sigma(i)}\|^2}{\sum_{i=1}^{N} \|q_i - \bar{q}\|^2} \times 100$$

where {qᵢ} are observed donor coordinates (centered), {pᵢ} are reference vertices, R is the optimal rotation (Kabsch), and σ is the optimal assignment (Hungarian/brute-force).

### SM5.2 Alignment and permutation search

- **CN ≤ 6**: Exact brute-force permutation search (N! ≤ 720)
- **CN > 6**: Hungarian algorithm with iterative alignment refinement
- **Numerical stability**: float64 arithmetic, fixed tie-breaking, consistent centering/scaling

### SM5.3 Reference polyhedra library

**Supplementary Table S6 | Library of ideal reference polyhedra.**

| CN | Shape 1 | Point group | Shape 2 | Point group |
|----|---------|-------------|---------|-------------|
| 2 | L-2 (linear) | D∞h | — | — |
| 3 | TP-3 (trigonal planar) | D₃h | TY-3 (T-shaped) | C₂v |
| 4 | SP-4 (square planar) | D₄h | T-4 (tetrahedral) | Td |
| 5 | SPY-5 (square pyramidal) | C₄v | TBPY-5 (trigonal bipyramidal) | D₃h |
| 6 | OC-6 (octahedral) | Oh | TPR-6 (trigonal prismatic) | D₃h |

### SM5.4 PathFinder methodology

For each CN pair (e.g., SP↔Td for CN=4), we map the CShM landscape:

1. Compute CShM(shape₁) and CShM(shape₂) for all retained entries at that CN
2. Define boundary: ΔCShM = |CShM(s₁) - CShM(s₂)| < threshold
3. Identify ridge points (local maxima in density along the boundary)
4. Classify entries by region: ideal_s1, distorted_s1, boundary, distorted_s2, ideal_s2

**CN=4 atlas**: 2,920 Pt SP-4 entries; SP↔Td boundary region contains entries with Delta:0 classification.

**CN=5 Berry pathway**: SPY↔TBPY interconversion mapped with 5-point ridge analysis. Ridge points correspond to the Berry pseudorotation transition state.

**CN=6 atlas**: Oh↔TPr distortion landscape for 6-coordinate entries. Boundary entries (17,791 total across all CN) are flagged with `_boundary` in the CoordRep L1 ShapeID.

### SM5.5 Shape token encoding in CoordRep

```
<ShapeBest:LABEL|Class:CLASS|Delta:BIN|V:s1,s2>
```

- `LABEL`: Best-fit reference shape
- `CLASS`: `ideal` (Δ ≥ 2), `dist` (1 ≤ Δ < 2), `boundary` (Δ < 1)
- `Delta`: Integer bin of ΔCShM
- `V:s1,s2`: Raw CShM values for the two reference shapes (2 decimal places)

---

## Supplementary Methods 6. Masked-field learning and tokenizer ablations

### SM6.1 Pretraining objective

The CoordRep language model is pretrained with masked language modeling (MLM) on tokenized CoordRep strings. The masking strategy covers three modes with controlled mixture:

| Mode | Mask target | Fraction | Purpose |
|------|-------------|----------|---------|
| Random token | Any token | 50% | General language model |
| Donor field | Donor element/index tokens only | 30% | Donor annotation task |
| Structure field | Metal/CN/shape tokens | 20% | Record completion |

### SM6.2 Training configuration

| Parameter | Value |
|-----------|-------|
| Architecture | Transformer encoder (BERT-style) |
| Hidden size | 768 |
| Attention heads | 12 |
| Layers | 12 |
| Max sequence length | 512 tokens |
| Batch size | 64 |
| Learning rate | 5 × 10⁻⁵ (linear warmup + cosine decay) |
| Warmup steps | 1,000 |
| Total training steps | 100,000 |
| Optimizer | AdamW (β₁=0.9, β₂=0.999, ε=10⁻⁸) |
| Weight decay | 0.01 |
| Dropout | 0.1 |
| Mask probability | 15% |
| Random seed | 42 (verified with seeds 0, 1, 2) |

### SM6.3 Factorized tokenizer ablation design

To assess whether donor, metal, and CN fields can be recovered independently, we compare:

| Variant | Metal token | CN token | Training data |
|---------|-------------|----------|---------------|
| Composite (pretrained) | Joint `[Metal:X\|ox:y\|d:z\|CN:n]` | Part of composite | Full pretraining corpus |
| Composite (matched) | Joint composite | Part of composite | Same n as factorized |
| Factorized | Separate `[Metal:X]` | Separate `[CN:n]` | Same n as composite |

The factorized variant enables measuring:
- Metal-only Top-1: 67.8%
- CN-only Top-1: 49.0%
- Joint (metal × CN) Top-1: 29.4%

vs composite:
- Metal-only Top-1: 13.8% (not decomposable — n.d.)
- Joint Top-1: 13.8%

### SM6.4 Convergence and syntax-first learning

Training dynamics show a characteristic pattern:
1. **Steps 0–5k**: Syntax tokens (brackets, delimiters) learned first (loss drops sharply)
2. **Steps 5k–30k**: Ligand SMILES tokens learned (gradual improvement)
3. **Steps 30k–100k**: Stereochemical and metal-state tokens refined (plateau approach)

This "syntax-first" pattern is reproducible across seeds and confirms that the grammar structure provides strong inductive bias.

---

## Supplementary Methods 7. Donor-field attribution and multidentate controls

### SM7.1 Masked-field donor annotation protocol

Donor-field attribution assesses whether the model can predict masked donor-site tokens from surrounding context. Four conditions isolate different information sources:

| Condition | Context visible | Interpretation |
|-----------|----------------|----------------|
| Full context | All fields except donor markers | Complete coordination record informs donor |
| Ligand SMILES masked | Metal, CN, shape, stereo visible; ligand SMILES removed | Tests metal/geometry → donor inference |
| Shuffled context | Token order randomized | Controls for positional artifacts |
| LigandFreq baseline | Most frequent donor for each ligand SMILES (training set) | Non-parametric lookup baseline |

Results (Top-1 accuracy, n=2,336):
- Full context: 85.7%
- LigandFreq baseline: 86.2%
- Ligand masked: 61.4%
- Shuffled: 21.6%

### SM7.2 Ligand-dependence control

The LigandFreq baseline achieves comparable Top-1 accuracy (86.2%) to the full-context MLM (85.7%), indicating that donor identity is largely **ligand-intrinsic**: knowing the ligand SMILES is sufficient to predict the donor element. This validates the masked-field probe design — the model has learned the same chemical knowledge that a frequency lookup captures.

The critical comparison is **ligand-masked** (61.4%): when ligand SMILES are removed, the model retains substantial donor prediction ability from metal/geometry context alone, well above shuffled control (21.6%).

### SM7.3 Multidentate donor-set recovery

For multidentate ligands (denticity ≥ 2), we evaluate whether the model can recover the **complete donor set** rather than just individual donors:

| Metric | Bidentate (d=2) | Tridentate+ (d≥3) |
|--------|-----------------|---------------------|
| Per-donor Top-1 | 86.7% | varies |
| All-donors-correct (full context) | 76.0% | — |
| All-donors-correct (ligand masked) | 38.5% | — |
| Hard-error rate | 24.0% | — |

**Hard-error definition**: A completion is a "hard error" if the predicted donor set is chemically invalid for the given ligand skeleton (e.g., predicting S as a donor for a ligand without sulfur atoms, or predicting 3 donors for a bidentate ligand).

### SM7.4 Invalid completion examples

Representative hard errors from `tool_a_multidentate_invalid_cases.csv`:
- Ligand `[H]N([H])C([H])([H])C([H])([H])N([H])[H]` (en): predicted {N, O} instead of {N, N}
- Ligand with P,N donors: predicted {N, N} instead of {P, N}

These errors confirm that the task requires genuine chemical understanding beyond memorization.

---

## Supplementary Methods 8. Graph, GNN, and 3D equivariant baselines

### SM8.1 2D graph neural network architectures

We evaluate three 2D GNN architectures for donor annotation:

| Architecture | Layers | Hidden dim | Aggregation | Edge features |
|--------------|--------|-----------|-------------|---------------|
| GCN | 4 | 256 | Mean | Bond type |
| GAT | 4 | 256 | Attention (4 heads) | Bond type |
| GIN | 4 | 256 | Sum + MLP | Bond type + ring membership |

**Input features** (node): Atomic number (one-hot, 100 dim), degree, formal charge, hybridization, aromaticity, ring membership.

**Training**: Binary classification per atom (donor vs non-donor). Adam optimizer, lr=10⁻³, 200 epochs, early stopping on validation F1.

### SM8.2 3D equivariant graph neural network (EGNN)

The EGNN baseline uses:
- **Architecture**: E(n) Equivariant GNN (Satorras et al., 2021)
- **Layers**: 5 message-passing layers
- **Hidden dim**: 128
- **Input**: Atom features + 3D coordinates
- **Equivariance**: SE(3)-equivariant message passing
- **Task**: Donor atom classification (same as 2D GNNs)

Result: **F1 = 0.998** (n=1,500 test complexes, all_correct = 99.07%).

This near-perfect performance with 3D information confirms that donor identity is recoverable from geometric features. CoordRep's text-based approach achieves comparable donor annotation **without requiring 3D coordinates** — operating from the serialized record alone.

### SM8.3 Semantic-decoy ranking with graph baselines

For semantic consistency (hard-negative ranking), we test whether GNNs can distinguish valid coordination records from semantically corrupted decoys:

| Method | AUROC | Top-1 | Interpretation |
|--------|-------|-------|----------------|
| Random | 0.508 | — | Chance |
| GIN Ranker (2D) | 0.522 | 0.015 | 2D topology insufficient |
| EGNN 3D | 0.499 | 0.048 | 3D geometry alone cannot rank |
| EGNN + metadata | 0.550 | 0.229 | Metal/CN/ox provides some signal |
| CoordRep MLM (PLL) | 0.554 | 0.004 | Pseudo-log-likelihood weak for listwise |
| **CoordRep Ranker** | **0.695** | **0.131** | Finetuned ranker best overall |
| CoordRep Ranker (strict) | 0.800 | — | Coordination-consistent decoys |
| CoordRep Ranker (stereo) | 0.948 | — | Stereochemical decoys |
| CoordRep Ranker (metal) | 0.791 | — | Metal substitution decoys |

### SM8.4 Metadata control experiment

To test whether the EGNN's semantic ranking failure is due to missing metadata, we augment EGNN with metal/CN/oxidation features as node-level attributes:

| Condition | AUROC | Notes |
|-----------|-------|-------|
| EGNN geometry only | 0.499 | No discrimination |
| EGNN + metal/CN/ox | 0.550 | Modest improvement |
| EGNN + all metadata + stereo labels | 0.697 | Stereo-specific improvement |

Conclusion: Graph topology and 3D geometry alone cannot determine record-level semantic consistency. The structured field relationships encoded in CoordRep provide information inaccessible to coordinate-based approaches.

---

## Supplementary Methods 9. Syntax validation, repair, and semantic consistency checks

### SM9.1 Grammar validation specification

A CoordRep string is **syntactically valid** if:
1. All brackets are matched and properly nested
2. Metal token present with valid element, oxidation, d-count, CN
3. Shape token present with valid shape label, class, and CShM values
4. All referenced ligands (L1, L2, ...) have corresponding SMILES entries
5. All donor references in stereo tokens correspond to valid donors in the referenced ligand
6. CN matches the total number of donor sites across all ligands

### SM9.2 Corruption types for evaluation

**Supplementary Table S7 | Corruption types used for syntax repair evaluation.**

| Type | Description | Example |
|------|-------------|---------|
| metal_swap | Replace metal element | Fe → Cu |
| cn_mismatch | Change CN without adjusting donors | CN:6 → CN:4 with 6 donors |
| shape_invalid | Invalid shape label | ShapeBest:XX |
| donor_count_wrong | Add/remove donor reference | Delete one trans token |
| ligand_missing | Remove a ligand SMILES entry | Delete \|L3=...\| |
| bracket_break | Unmatched brackets | Missing closing ] |
| ox_invalid | Impossible oxidation state | ox:+9 |
| stereo_inconsistent | Trans pair references non-existent ligand | L5:N:1 when only L1–L4 exist |

### SM9.3 Repair protocol

Syntax repair uses the pretrained MLM in a masked-prediction loop:
1. Identify invalid tokens via grammar checker
2. Mask invalid tokens
3. Predict replacements via MLM
4. Verify repaired string passes grammar validation
5. Iterate if needed (max 3 rounds)

### SM9.4 CoordRep-Ranker for semantic consistency

The CoordRep-Ranker is a contrastive model trained to score record-level semantic consistency:

- **Architecture**: Transformer encoder + linear classification head
- **Training**: Contrastive pairs (valid record vs semantically corrupted decoy)
- **Decoy generation**: Field-specific corruption (metal swap, ligand swap, stereo permutation, CN change)
- **Evaluation**: AUROC on held-out test set with 20 decoys per ground-truth record (n=1,501 test records)

---

## Supplementary Methods 10. Full-CSD CoordRep-PathFinder audit

### SM10.1 Scan procedure

The full-CSD scan processes all 1,413,222 entries in CSD release 2024.3:

```bash
/data/miniconda3/envs/1701/bin/python scripts/run_full_csd_scan.py
```

Runtime: approximately 4 hours on a single machine (sequential entry processing via CSD Python API).

### SM10.2 Pipeline per entry

For each retained entry:
1. `csd_entry_to_raw_molecule()` → atom coordinates + bond orders
2. `encode_molecule()` → CoordComplex object
3. `canonicalize_complex()` → canonical CoordComplex
4. `to_string()` → CoordRep string
5. `is_valid_coordrep(strict=True)` → validation gate
6. `extract_identity_keys()` → L0–L3 hashes

### SM10.3 Identity-layer retention statistics

The CSD-derived corpus was analyzed for identity-layer diversity:

| Level | Unique keys | Meaning |
|-------|-------------|---------|
| L0 (StateKey) | ~124,000 | Near-unique (polymorph pairs collapse) |
| L1 (ShapeID) | ~98,000 | Geometry-distinct |
| L2 (TopoID) | ~85,000 | Topology-distinct |
| L3 (ConnID) | ~62,000 | Connectivity families |

### SM10.4 Transfer evaluation

Models trained on tmQM/COD were evaluated on CSD-derived records without retraining:

| Task | Method | Metric | Value |
|------|--------|--------|-------|
| Syntax repair | Rule-based | Valid rate | 99.7% |
| Syntax repair | CoordRep MLM | Valid rate | 85.2% |
| Semantic ranking | CoordRep Ranker | AUROC | 0.695 (within-domain) |

---

## Supplementary Methods 11. Reproducibility package and file manifest

### SM11.1 Package contents

The reproducibility package includes:
- All training/evaluation scripts
- Pretrained model checkpoints
- Processed datasets (CoordRep strings, not raw CSD coordinates)
- Evaluation result CSVs
- Figure-generation scripts

### SM11.2 One-command reproduction

```bash
# Full evaluation pipeline
python scripts/run_evaluation_suite.py --config configs/full_eval.yaml

# Box 1 records
python scripts/export_box1_coordrep_records.py

# Fig. 4 data export
python scripts/export_fig4_field_learning.py

# Representation gap benchmark
python scripts/export_representation_gap_benchmark.py
```

### SM11.3 Deterministic seeds

| Run | Seed | Purpose |
|-----|------|---------|
| Primary | 42 | All reported results |
| Replicate 1 | 0 | Variance estimation |
| Replicate 2 | 1 | Variance estimation |
| Replicate 3 | 2 | Variance estimation |

### SM11.4 Computational requirements

| Task | Hardware | Runtime |
|------|----------|---------|
| Pretraining (100k steps) | 1× A100 80GB | ~8 hours |
| Full CSD scan | CPU (CSD API) | ~4 hours |
| GNN baselines (all) | 1× V100 32GB | ~2 hours |
| Evaluation suite | CPU | ~30 minutes |

### SM11.5 Data availability

**Supplementary Table S8 | Data availability and access.**

| Data | Source | Access | License |
|------|--------|--------|---------|
| CSD structures | CCDC | Academic license required | CCDC license |
| COD structures | COD | Open access | CC0 |
| tmQM geometries | Balcells & Skjelstad | Open access | CC-BY-4.0 |
| CoordRep strings | This work | Included in reproducibility package | CC-BY-4.0 |
| Pretrained models | This work | Included | CC-BY-4.0 |

---

## Supplementary Results 1. Dataset statistics and scope audit

### SR1.1 Dataset composition

The pretraining corpus contains 86,665 complexes (tmQM) supplemented by COD-derived entries. After stratified splitting:

| Split | n | Purpose |
|-------|---|---------|
| Train | 69,332 (80%) | Pretraining |
| Validation | 8,666 (10%) | Hyperparameter tuning |
| Test | 8,667 (10%) | Final evaluation |

### SR1.2 Metal distribution

Top-10 metals by frequency: Fe (18.2%), Cu (12.8%), Zn (9.4%), Co (8.7%), Ni (8.1%), Mn (5.9%), Cr (4.2%), Ru (3.8%), Pd (3.5%), Pt (3.1%).

### SR1.3 CN distribution

| CN | Fraction |
|----|----------|
| 4 | 28.3% |
| 5 | 18.7% |
| 6 | 42.1% |
| 2 | 3.2% |
| 3 | 7.7% |

---

## Supplementary Results 2. Canonicalization and identity robustness

### SR2.1 Determinism verification

Empirical tests on 1,000 randomly sampled complexes:
- Atom permutation (10× per complex): 10,000/10,000 identical L0 (100%)
- Coordinate rotation (10× per complex): 10,000/10,000 identical L0 (100%)
- Numerical noise (σ=0.001 Å): 10,000/10,000 identical (after CShM rounding)

### SR2.2 Identity-layer statistics

From the CSD corpus (124,837 entries):
- L3 collision rate (different complexes sharing L3): Expected and by design — L3 groups connectivity families
- L0 uniqueness: 99.8% unique (0.2% are polymorph pairs that correctly collapse)

---

## Supplementary Results 3. Representation-operation challenge details

### SR3.1 Per-operation results

**Invariance test** (20 polymorph pairs):
- CoordRep L0: 20/20 correctly collapsed (100%)
- SMILES multiset: 20/20 (100%) — but over-merges non-polymorphs
- InChI: 20/20 (100%) — same caveat

**Stereo separation** (20 cis/trans or fac/mer pairs):
- CoordRep: 20/20 correctly separated (100%)
- SMILES: 0/20 (0%) — no stereo field
- InChI: 0/20 (0%)
- Raw 3D: 15/20 (75%) — angles recoverable but noisy

**Family linking** (20 same-connectivity different-geometry pairs):
- CoordRep L3: 20/20 correctly linked (100%)
- SMILES: 0/20 — no geometry awareness
- CShM: 0/20 — no connectivity information

**Boundary detection** (20 boundary entries):
- CoordRep: 20/20 flagged with `_boundary` + ΔCShM (100%)
- CShM vector: 12/20 (60%) — raw ΔCShM computable but no structured flag
- Others: 0/20

---

## Supplementary Results 4. Continuous geometry supplementary analyses

### SR4.1 CN=4 SP↔Td landscape

From 124,837 retained CSD entries, CN=4 subset:
- Total CN=4 entries: ~35,000
- SP-dominant (CShM_SP < CShM_Td): ~22,000
- Td-dominant (CShM_Td < CShM_SP): ~11,000
- Boundary (ΔCShM < 1.0): ~2,000

**Supplementary Figure S1**: Hexbin density map of CShM(SP) vs CShM(Td) for all CN=4 entries. The Berry pathway minimum-energy path is visible as a density ridge connecting the SP and Td basins.

### SR4.2 CN=5 Berry pathway

The SPY↔TBPY interconversion follows the Berry pseudorotation coordinate:
- Ridge points identified at ΔCShM ≈ 0 (exact boundary)
- 5-coordinate boundary entries are enriched in d⁸ metals (Ni²⁺, Cu²⁺)

### SR4.3 CN=6 distortion atlas

Oh↔TPr pathways show:
- Majority of CN=6 entries are Oh-dominant
- TPr entries are enriched in d⁰ and d¹ metals
- Boundary entries (ΔCShM < 1.0) are flagged in CoordRep with `TPr/Oh_boundary`

---

## Supplementary Results 5. Tokenizer ablations and masked-field learning controls

### SR5.1 Factorized tokenizer full results

**Supplementary Table S9 | Factorized vs composite tokenizer ablation.**

| Variant | Metal Top-1 | CN Top-1 | Joint (M×CN) Top-1 | Donor Top-1 |
|---------|-------------|----------|---------------------|-------------|
| Factorized | 67.8% | 49.0% | 29.4% | 27.5% |
| Composite (matched n) | 13.8% | 13.8% | 13.8% | 15.1% |
| Composite (pretrained) | n.d. | n.d. | n.d. | 14.6% |

The factorized tokenizer enables **independent field recovery**: metal (67.8%) and CN (49.0%) can each be predicted with reasonable accuracy from context, while their joint prediction (29.4%) reflects the multiplicative difficulty. The composite tokenizer cannot decompose — it achieves only 13.8% joint Top-1.

### SR5.2 Donor recovery by CN

| CN | Full context Top-1 | Ligand masked Top-1 | n |
|----|-------------------|---------------------|---|
| 4 | 87.2% | 63.1% | 812 |
| 5 | 84.9% | 59.8% | 445 |
| 6 | 85.1% | 61.2% | 1,079 |

### SR5.3 Donor recovery by donor element

| Element | Top-1 | n | Notes |
|---------|-------|---|-------|
| N | 89.3% | 892 | Most common donor |
| O | 83.7% | 654 | Second most common |
| S | 81.2% | 298 | |
| P | 79.4% | 241 | |
| Cl | 91.5% | 189 | Highly predictable |
| C | 72.8% | 156 | Most challenging |

---

## Supplementary Results 6. Donor-set recovery for ligand denticity subsets

### SR6.1 Per-denticity breakdown

**Supplementary Table S10 | Multidentate donor-set recovery by denticity.**

| Denticity | Per-donor Top-1 | Per-donor Top-5 | All-correct | Hard-error rate | n |
|-----------|-----------------|-----------------|-------------|-----------------|---|
| 1 (mono) | 84.7% | 97.9% | 84.7% | — | ~6,000 |
| 2 (bi) | 86.7% | 98.4% | 76.0% | 24.0% | ~1,200 |
| 3 (tri) | 83.5% | 96.8% | 61.2% | 31.5% | ~300 |
| 4+ | 80.1% | 94.3% | 48.7% | 38.2% | ~100 |

### SR6.2 Hard-error analysis

Hard errors (chemically invalid donor predictions) increase with denticity because:
1. Higher denticity requires correctly predicting multiple dependent donor sites
2. Ligand-masked condition removes the strongest cue (ligand skeleton)
3. Training data has fewer examples of high-denticity ligands

---

## Supplementary Results 7. Graph baseline and semantic-decoy ranking

### SR7.1 2D GNN donor annotation

**Supplementary Table S11 | 2D GNN donor annotation results.**

| Model | F1 | All-correct | n |
|-------|-----|-------------|---|
| Random | 0.228 | 10.2% | 2,000 |
| GlobalFreq | 0.377 | 20.2% | 2,000 |
| CondFreq | 0.377 | 20.3% | 2,000 |
| LigandFreq | 0.896 | 79.7% | 2,000 |
| GCN (ligand only) | 0.855 | 75.6% | 2,000 |
| GIN (ligand only) | 0.861 | 76.7% | 2,000 |
| GIN (ligand + context) | 0.904 | 83.2% | 2,000 |
| **EGNN (3D)** | **0.998** | **99.1%** | 1,500 |

### SR7.2 3D EGNN analysis by CN

| CN | EGNN F1 | n |
|----|---------|---|
| 4 | 0.999 | 500 |
| 5 | 0.997 | 400 |
| 6 | 0.998 | 600 |

### SR7.3 Semantic-decoy ranking by decoy type

**Supplementary Table S12 | CoordRep-Ranker AUROC by decoy type.**

| Decoy type | MLM PLL | Ranker (finetuned) | EGNN | EGNN+meta |
|------------|---------|-------------------|------|-----------|
| Metal swap | 0.677 | 0.791 | 0.499 | 0.550 |
| Stereo permutation | 0.745 | 0.948 | 0.500 | 0.697 |
| Ligand swap | 0.651 | 0.654 | — | — |
| Boundary geometry | 0.608 | 0.698 | — | — |
| Co-ligand swap | 0.484 | 0.553 | — | — |

---

## Supplementary Results 8. Syntax repair and semantic consistency casebook

### SR8.1 Repair results by method

**Supplementary Table S13 | Syntax repair evaluation (synthetic test set, n=2,000).**

| Method | Valid rate (before) | Valid rate (after) | Exact match | Field match rate |
|--------|--------------------|--------------------|-------------|------------------|
| Corrupted only (no repair) | 50.9% | 50.9% | 0% | 84.0% |
| Rule-based repair | 50.9% | 99.8% | 12.4% | 85.9% |
| Edit distance | 50.9% | 100% | 12.4% | 85.9% |
| Char bigram | 50.9% | 100% | 12.4% | 85.9% |
| CoordRep MLM | 50.2% | 75.8% | 0% | 59.9% |

### SR8.2 CSD transfer results

**Supplementary Table S14 | Syntax repair on CSD-derived corrupted records (n=1,876).**

| Method | Valid rate (after) | Exact match | Field match rate |
|--------|-------------------|-------------|------------------|
| Rule-based | 100% | 24.6% | 89.0% |
| CoordRep MLM | 85.2% | 0% | 60.6% |

The rule-based approach achieves near-perfect syntax recovery because CoordRep's grammar is sufficiently constrained. The MLM approach, while lower in absolute valid rate, recovers semantically appropriate content in ambiguous cases.

---

## Supplementary Results 9. Full-CSD derived outputs and boundary-record statistics

### SR9.1 Scan summary

| Metric | Value |
|--------|-------|
| Total CSD entries scanned | 1,413,222 |
| Entries retained | 124,837 |
| Retention rate | 8.83% |
| Boundary records (ΔCShM < 1.0) | 17,791 (14.3% of retained) |
| Unique L3 ConnIDs | ~62,000 |
| Unique L1 ShapeIDs | ~98,000 |

### SR9.2 Boundary-record statistics by CN

| CN | n_total | n_boundary | Boundary fraction |
|----|---------|-----------|-------------------|
| 4 | ~35,000 | ~5,200 | 14.9% |
| 5 | ~23,000 | ~4,800 | 20.9% |
| 6 | ~52,000 | ~6,800 | 13.1% |
| 2–3 | ~14,000 | ~1,000 | 7.1% |

CN=5 has the highest boundary fraction, consistent with the known flexibility of 5-coordinate geometries (Berry pseudorotation pathway).

### SR9.3 L3 family geometry trajectories

For L3 families containing ≥ 3 members with different L1 ShapeIDs, we extract "geometry trajectories" showing how the same connectivity can adopt different geometries:

Example: `Co|CN6|C#N;C#N;C#N;[H]N([H])C([H])([H])C([H])([H])N([H])C([H])([H])C([H])([H])C([H])([H])N([H])[H]`
- fac isomer: ShapeBest=Oh, ΔCShM=1.22 (IDARAG)
- mer isomer: ShapeBest=TPr, ΔCShM=0.33, boundary (IDAREK)

This demonstrates CoordRep's multi-resolution design: L3 links the family, while L1/L0 resolve the geometric states.

---

## Supplementary Figures

**Supplementary Figure S1** | CN=4 continuous shape atlas. Hexbin density map of CShM(SP) vs CShM(Td) for all CN=4 entries in the CSD-derived corpus. Ideal square-planar entries cluster at (SP≈0, Td≈30); ideal tetrahedral entries cluster at (SP≈30, Td≈0). Boundary entries (ΔCShM < 1.0) populate the diagonal region.

**Supplementary Figure S2** | CN=5 Berry pseudorotation pathway. CShM(SPY) vs CShM(TBPY) for all CN=5 entries. The interconversion pathway from SPY to TBPY is visible as a density ridge. Ridge points correspond to the geometric transition state.

**Supplementary Figure S3** | CN=6 Oh↔TPr distortion atlas. CShM(Oh) vs CShM(TPr) for all CN=6 entries. Most entries are Oh-dominant; TPr entries are enriched in early transition metals.

**Supplementary Figure S4** | Pretraining dynamics. (a) Training and validation loss over 100k steps. (b) Syntax token accuracy reaches >99% by step 5k. (c) Donor token accuracy plateaus at ~85% by step 50k.

**Supplementary Figure S5** | Dataset composition. (a) Metal element distribution. (b) CN distribution. (c) Ligand denticity distribution.

**Supplementary Figure S6** | Full representation capability heatmap. Expanded version of Fig. 2C showing all six operations × six representations with Y/P/N annotations and mechanism notes.

**Supplementary Figure S7** | Geometry region enrichment by metal d-count. For CN=4 (SP vs Td) and CN=6 (Oh vs TPr), relative enrichment of each geometry region by d-electron count.

**Supplementary Figure S8** | L3 family geometry trajectory examples. Three families with ≥3 members showing distinct L1 shape states, illustrating multi-resolution identity.

**Supplementary Figure S9** | Factorized tokenizer ablation. Bar chart comparing factorized vs composite tokenizer for metal, CN, joint, and donor recovery (Top-1 accuracy).

**Supplementary Figure S10** | Donor recovery by coordination number. Grouped bar chart (CN=4,5,6) showing full-context vs ligand-masked vs shuffled conditions.

**Supplementary Figure S11** | Donor-set recovery by denticity. Bar chart of all-donors-correct rate for mono-, bi-, and tridentate ligands.

**Supplementary Figure S12** | GNN donor annotation comparison. Bar chart of F1 scores: Random, GlobalFreq, GCN, GAT, GIN, EGNN.

**Supplementary Figure S13** | Semantic-decoy AUROC by decoy type. Grouped bars for CoordRep Ranker vs EGNN vs MLM PLL across metal, stereo, ligand, and boundary decoys.

**Supplementary Figure S14** | Syntax repair accuracy by corruption type. Stacked bar chart showing repair success rate for each of the 8 corruption types.

**Supplementary Figure S15** | Representative correction cases. 20 examples where the MLM produces valid repairs that differ from rule-based corrections, showing contextual reasoning.

**Supplementary Figure S16** | Full-CSD waterfall filter. Horizontal bar chart showing entry counts at each filtering stage.

---

## Source Data and File Manifest

All source data files are in `revision_results/` with the following structure:

| Directory | Key files | Corresponding figure/table |
|-----------|-----------|---------------------------|
| `box1_coordrep_records/` | `box1_full_records.jsonl`, `box1_si_full_records.txt` | Box 1, SM2 |
| `representation_gap_benchmark/` | `representation_capability_matrix.csv`, `challenge_pair_examples.jsonl` | Fig. 2, SM4 |
| `csd_pathfinder_full/` | `full_csd_filtering_waterfall.csv`, `cn*_atlas.csv` | SM1, SM5, SM10 |
| `tool_a_ablation/` | `tool_a_ablation_summary.csv` | Fig. 4C |
| `factorized_token/` | `factorized_ablation_summary.csv` | Fig. 4D |
| `multidentate/` | `tool_a_by_denticity_final.csv` | Fig. 4E |
| `gnn_baselines/` | `donor_annotation_summary.json`, `3d_donor_annotation_summary.csv` | Fig. 4F |
| `coordrep_ranker/` | `ranker_summary.csv`, `ranker_by_decoy_type.csv` | Fig. 4F |
| `figure_ready/` | `gnn_comparison_summary.csv` | Fig. 4F summary |
| `toolb_sequence_baselines/` | `toolb_baseline_summary.csv`, `toolb_by_corruption_type.csv` | SM9 |
| `fig4_field_learning_revision/` | All `fig4*` CSV/JSON/PNG | Fig. 4 panels |
| `scope/` | `coordrep_v1_scope_table.csv` | SM1 |
| `csd_external_summary_only/` | `csd_filter_waterfall.csv`, `csd_ranker_transfer_summary.csv` | SM10 |
