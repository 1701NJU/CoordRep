# Evidence Matrix — Reviewer Comment ↔ Revision Evidence Mapping

## Reviewer 1

### R1.1 — Scope: multinuclear / MOF / ηn exclusion [critical → resolved]

| Item | Detail |
|------|--------|
| Concern | Pipeline rejects multi-metal, hapticity >1, MOFs, POMs, bimetallic catalysts, metallocenes |
| Evidence | CSD audit: 200k → 17k (8.5%); top rejections: no_TM 50%, multinuclear 24%, hapticity quantified |
| Key result | v1 = mononuclear η1 CN 2–6; scope is explicit, not universal |
| Manuscript | Discussion/Limitations: v1/v2 scope table |
| SI | Scope table with rejection counts |
| Files | `fig5_application_summary.csv`, `csd_external/summary.json` |

### R1.2 — Composite metal/CN token [critical → resolved]

| Item | Detail |
|------|--------|
| Concern | Fusing metal+CN prevents independent exploration |
| Evidence | Tokenizer audit found composite mismatch; factorized default adopted |
| Key result | metal Top-1: 13.8%→67.8%; CN: 13.8%→49.0%; joint: 13.8%→29.4% |
| Manuscript | Methods: tokenizer section revised |
| SI | Tokenizer audit detail + by-task breakdown |
| Files | `tokenizer_factorized_summary.csv`, `factorized_token/tokenizer_audit_summary.json` |

### R1.3 — Tool A physical grounding / DFT / statistical bias [critical → reframed]

| Item | Detail |
|------|--------|
| Concern | Tool A may only learn frequency; lacks physical validation |
| Evidence | Tool A ablation: full=85.7%, LigandFreq=86.2%, no_lig=61.4%, no_geom≈full |
| Key result | Donor identity is ligand-intrinsic; EGNN F1=0.998 confirms geometry-accessible |
| Manuscript | Results: reframe as donor annotation; add ablation table |
| SI | Full ablation by CN/denticity/donor element |
| Files | `tool_a_ablation/tool_a_ablation_summary.csv`, `gnn_baselines/donor_annotation_summary.csv` |

### R1.4 — Bidentate / multidentate support [major → resolved]

| Item | Detail |
|------|--------|
| Concern | Need scaling beyond monodentate |
| Evidence | 76.1% complexes multidentate; bidentate Top-1=86.7%; validator error=0.46% |
| Manuscript | Results: multidentate paragraph + 3 annotated examples |
| SI | Coverage by CN/metal/source; validator failures |
| Files | `multidentate_summary.csv`, `multidentate/tool_a_multidentate_summary.json` |

---

## Reviewer 2

### R2.1 — Why not GNN / equivariant GNN? [critical → resolved]

| Item | Detail |
|------|--------|
| Concern | Modern GNNs operate directly on 3D structures |
| Evidence | EGNN donor F1=0.998; EGNN hard-neg AUROC=0.499; CoordRep-Ranker strict=0.800 |
| Manuscript | Results: GNN comparison subsection |
| SI | Full GNN tables |
| Files | `gnn_comparison_summary.csv`, `gnn_baselines/3d_hard_negative_summary.csv` |

### R2.2 — No GNN benchmark [critical → resolved]

| Item | Detail |
|------|--------|
| Evidence | Three tasks: donor (EGNN 0.998), hard-neg (EGNN 0.499 vs CoordRep 0.800), identity (WL/ECFP/L3) |
| Files | `gnn_comparison_summary.csv` |

### R2.3 — Borderline coordination / ηn [major → resolved for boundary]

| Item | Detail |
|------|--------|
| Evidence | AFOSIA Δ=0.01; 13.4% CSD boundary; ηn future scope |
| Files | `casebook_maintext_summary.csv`, `csd_external/summary.json` |

### R2.4 — Fraction excluded [major → resolved]

| Item | Detail |
|------|--------|
| Evidence | 200k→17k; rejection breakdown quantified |
| Files | `fig5_application_summary.csv` |

### R2.5 — Subtle geometry distinction [major → resolved]

| Item | Detail |
|------|--------|
| Evidence | L0 21.5%, L3 88.9%; perturbation σ=0.01: L0=33% L3=94%; ACUWOK 12/1 |
| Files | `fig2_identity_layers.csv` |

### R2.6 — Reason for MLM [moderate → reframed]

| Item | Detail |
|------|--------|
| Evidence | MLM operates on structured fields; GNNs included; syntax by rules, semantics by models |

### R2.7 — MLM advantage over GNN? [major → reframed]

| Item | Detail |
|------|--------|
| Evidence | Not for donor; yes for strict semantic (0.80 vs 0.50 AUROC) |
| Files | `gnn_comparison_summary.csv` |

### R2.8 — Tokenizer sensitivity [moderate → resolved]

| Item | Detail |
|------|--------|
| Evidence | Composite→factorized: metal 67.8%, CN 49.0% |
| Files | `tokenizer_factorized_summary.csv` |

### R2.9 — Statistical correlations vs chemistry [major → resolved]

| Item | Detail |
|------|--------|
| Evidence | Field ablation: no_stereo→stereo AUROC 0.948→0.500; no_metal→0.792→0.500 |
| Files | `gnn_baselines/taskB_coordrep_field_ablation.csv` |

### R2.10 — No chemical insight [major → resolved]

| Item | Detail |
|------|--------|
| Evidence | 4 QC-verified CSD cases: ACUWOK, AFOSIA, AGOTIA, CIJWUO |
| Files | `casebook_maintext_summary.csv`, `csd_casebook/casebook_qc_report.md` |

### R2.11 — Standard molecular graph comparison [moderate → resolved]

| Item | Detail |
|------|--------|
| Evidence | SMILES/WL/ECFP/CoordRep L3 compared |
| Files | `gnn_comparison_summary.csv`, `gnn_baselines/csd_identity_baselines.csv` |

### R2.12 — Discovery / generation lacks validation [critical → reframed]

| Item | Detail |
|------|--------|
| Evidence | De novo claim removed; CSD curation workflow replaces |

### R2.13 — Electronic structure effects [moderate → future_scope]

| Item | Detail |
|------|--------|
| Evidence | Explicit limitation; future DFT integration |

---

## Reviewer 3

### R3.1 — Missing literature [moderate → partially_resolved]

| Item | Detail |
|------|--------|
| Evidence | References to add: Pidko, Corminboeuf, Kulik, hemilability datasets |
| Status | Requires manuscript citation update |

### R3.2 — Geometry indices / Fey / Cundari [minor → partially_resolved]

| Item | Detail |
|------|--------|
| Evidence | Add CShM/geometry descriptor references |

### R3.3 — WL hash discussion [moderate → resolved]

| Item | Detail |
|------|--------|
| Evidence | WL recall=0.983 prec=0.800; CoordRep L3 recall=0.830 prec=0.816 |
| Files | `gnn_comparison_summary.csv` |

### R3.4 — Spin-state caveat [minor → partially_resolved]

| Item | Detail |
|------|--------|
| Evidence | Add caveat about 3d metal geometry and spin-state convolution |

### R3.5 — GitHub repo missing [moderate → partially_resolved]

| Item | Detail |
|------|--------|
| Evidence | Fix URL; code uploaded; CSD data license caveat |

---

## Reviewer 4

### R4.1 — Canonicalization claim too strong [critical → resolved]

| Item | Detail |
|------|--------|
| Evidence | CoordRep-ID L0–L3; L0=21.5%, L3=88.9%; ACUWOK 12/1 |
| Manuscript | Revise to "multi-resolution identity hierarchy" |
| Files | `fig2_identity_layers.csv` |

### R4.2 — Full syntax insufficiently documented [moderate → resolved]

| Item | Detail |
|------|--------|
| Evidence | 3 annotated examples: Pt monodentate, Pd bidentate, Ir octahedral multi-stereo |

### R4.3 — Tool A ligand-SMILES ablation [major → resolved]

| Item | Detail |
|------|--------|
| Evidence | full=85.7%, LigandFreq=86.2%, no_lig=61.4%, shuffled=21.6% |
| Files | `tool_a_ablation/tool_a_ablation_summary.csv` |

### R4.4 — Multidentate ligand tracking [major → resolved]

| Item | Detail |
|------|--------|
| Evidence | 76.1% multidentate; bidentate Top-1=86.7% |
| Files | `multidentate_summary.csv` |

### R4.5 — Tool A overpromotion [critical → resolved]

| Item | Detail |
|------|--------|
| Evidence | Inverse design claim removed; reframed as donor annotation + GNN comparison |
