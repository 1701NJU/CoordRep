# Response Letter Bullets

## Reviewer 1

### R1.1 — Scope: multinuclear / MOF / ηn exclusion

- We agree CoordRep v1 has an explicit scope boundary.
- We added a CSD external audit: 200,000 entries scanned, 17,038 retained (8.5%).
- Top rejection reasons quantified: no transition metal 50%, multinuclear 24%, no 3D 7%, disorder 4%.
- v1 scope is mononuclear η1 with CN 2–6. This is stated in the revised Abstract, Methods, and Discussion.
- Future extensions (CoordRep-Multi, CoordRep-Periodic, CoordRep-Haptic) are noted in Discussion.

### R1.2 — Composite metal/CN token

- We discovered a tokenizer implementation mismatch: `[Metal:Fe|ox:+2|d:d6|CN:6]` was treated as a single composite token despite factorized tokens being defined in code.
- We adopted the factorized metal-field tokenizer as the revised default.
- Matched training shows metal Top-1 improves from 13.8% to 67.8%, CN from 13.8% to 49.0%.
- Caveat: this is matched short training, not directly comparable to the fully pretrained model.

### R1.3 — Tool A physical grounding

- We reframe Tool A from "inverse ligand design" to donor annotation and field attribution.
- Ablation: full_context MLM = 85.7%, LigandFreq = 86.2%, no_ligand_SMILES ≈ 61%, no_geometry_stereo ≈ full_context, shuffled control = 21.6%.
- Conclusion: donor identity is ligand-intrinsic; Tool A is a donor annotator, not a stability predictor.
- EGNN baseline (F1 = 0.998) confirms donor assignment is geometry-accessible.
- We no longer claim thermodynamic or synthetic validation.

### R1.4 — Bidentate / multidentate

- 76.1% of training complexes contain at least one multidentate ligand.
- Bidentate per-donor Top-1 = 86.7%, comparable to monodentate 84.7%.
- Validator hard-error rate = 0.46%.
- Three annotated examples added: Pt/CN4 monodentate, Pd bidentate, Ir octahedral multi-stereo.
- CoordRep v1 is mononuclear η1, not monodentate-only.

---

## Reviewer 2

### R2.1 + R2.2 — GNN baselines

- We added 2D GNN (GCN, GIN) and 3D EGNN baselines for three tasks.
- **Donor annotation**: EGNN F1 = 0.998, GIN+context F1 = 0.904, LigandFreq F1 = 0.896. EGNN dominates.
- **Hard-negative semantic ranking**: EGNN AUROC = 0.499 (random-level), GIN AUROC = 0.522. CoordRep-Ranker strict semantic AUROC = 0.800, stereo_hard = 0.948.
- **CSD identity**: WL hash recall = 0.983, precision = 0.800; CoordRep L3 recall = 0.830, precision = 0.816.
- Conclusion: 3D GNNs solve geometry-local tasks; CoordRep handles field-level semantic records. Complementary.

### R2.3 — Borderline coordination

- CoordRep-ID includes boundary labels based on CShM competition (threshold Δ < 1.0).
- AFOSIA case: CShM(TPr) = 4.49, CShM(Oh) = 4.48, Δ = 0.01 — extreme boundary.
- 13.4% of CSD retained entries are boundary cases.
- ηn haptic motifs are explicitly excluded in v1.

### R2.4 — Fraction excluded

- CSD audit quantifies: 200k → 17k, 8.5% retention.
- Rejection breakdown provided in both main text and SI.

### R2.5 — Subtle geometry distinction

- CoordRep-ID defines L0 (full string), L1 (shape + constraints), L2 (metal + CN + shape + ligands), L3 (metal + CN + ligands).
- Perturbation: at σ = 0.01 Å, L0 retention = 33%, L3 = 94%.
- CSD families: L0 match = 21.5%, L3 match = 88.9%.
- ACUWOK: 12 entries, 12 unique L0, single L3.

### R2.6 + R2.7 — MLM vs GNN

- MLM does not outperform GNN for donor annotation. We acknowledge this.
- CoordRep advantage is in the representation layer: field-level semantic consistency, auditable records, multi-resolution identity.
- CoordRep-Ranker strict semantic AUROC = 0.800 vs EGNN = 0.499 demonstrates this.

### R2.8 — Tokenizer sensitivity

- Composite→factorized ablation demonstrates sensitivity.
- Factorized tokenizer adopted as revised default.

### R2.9 — Statistical correlations vs chemistry

- Field ablation: removing stereo tokens → stereo_hard AUROC drops from 0.948 to 0.500.
- Removing metal tokens → metal_hard AUROC drops from 0.792 to 0.500.
- Each CoordRep field contributes to its corresponding decoy type discrimination.

### R2.10 — No chemical insight

- Four QC-verified CSD casebook cases demonstrate chemical insight:
  - ACUWOK: identity linking across conformers.
  - AFOSIA: boundary-aware geometry annotation.
  - AGOTIA: deterministic grammar repair.
  - CIJWUO: stereo semantic consistency checking.

### R2.11 — Standard molecular graph comparison

- SMILES multiset, WL hash, ECFP, Metal+SMILES all included.
- CoordRep provides auditable multi-resolution identity, not just a single hash.

### R2.12 — Discovery / generation lack validation

- De novo generation claim removed from Abstract and throughout.
- Replaced with CSD application and curation workflow.
- Tool A reframed as donor annotation.

### R2.13 — Electronic structure effects

- Explicitly stated: CoordRep does not model electronic structure or thermodynamic stability.
- Future integration with DFT labels / electronic descriptors noted.

---

## Reviewer 3

### R3.1 — Missing literature

- References added: Pidko (coordination catalyst databases), Corminboeuf (ligand descriptors), Kulik (ML for transition metals), hemilability and isomerism datasets.
- CoordRep positioned relative to existing coordination chemistry data infrastructure.

### R3.2 — Geometry indices

- References added: Fey (ligand knowledge bases), Cundari (computational TM chemistry), CShM origin.
- CoordRep builds on, not replaces, coordination geometry indices.

### R3.3 — WL hash

- WL hash included in CSD identity baselines: recall = 0.983, precision = 0.800.
- CoordRep L3: recall = 0.830, precision = 0.816. Multi-resolution hierarchy is the differentiator.

### R3.4 — Spin-state caveat

- Caveat added: 3d metal geometry distributions may reflect spin-state effects not explicitly separated in CoordRep v1.
- Spin-state labels not systematically available in tmQM; future stratification possible.

### R3.5 — GitHub repo

- Repository URL corrected; code and scripts uploaded.
- CSD raw coordinates not redistributable (CCDC license); derived statistics and anonymized excerpts provided.

---

## Reviewer 4

### R4.1 — Canonicalization claim too strong

- Revised from "InChI-like identifier" to "multi-resolution identity hierarchy".
- L0 is geometry-state sensitive (21.5% CSD family match); L3 is robust chemical identity (88.9%).
- ACUWOK: 12 entries share single L3 despite 12 unique L0.
- Perturbation robustness quantified.

### R4.2 — CoordRep syntax documentation

- Three fully annotated examples added: monodentate Pt/CN4, bidentate Pd, octahedral Ir with multi-stereo/multidentate.
- Field-by-field annotation in SI.

### R4.3 — Tool A ligand-SMILES ablation

- Ablation: full = 85.7%, LigandFreq = 86.2%, no_lig_SMILES = 61.4%, no_geometry_stereo ≈ full, shuffled = 21.6%.
- Donor identity is ligand-intrinsic. Geometry/stereo fields minimally affect donor prediction.
- Tool A reframed as donor annotation not inverse design.

### R4.4 — Multidentate tracking

- 76.1% of complexes contain multidentate ligands.
- Bidentate Top-1 = 86.7%; validator error = 0.46%.
- CoordRep v1 = mononuclear η1, not monodentate-only.

### R4.5 — Tool A overpromotion

- Inverse ligand design claim removed.
- Tool A reframed as donor annotation and field attribution.
- GNN baselines and LigandFreq included as comparators.
- Ranker framed as semantic consistency probe, not stability predictor.
