# Supplementary Information (SI) Materials

## File Structure

```
SI/
├── Supplementary_Information.md    # Main document (S1-S10, ~1000 lines)
├── README.md                       # This file
│
├── figures/                        # Figures (7 panels)
│   ├── Figure_S1_training_curves.pdf       # Training curves
│   ├── Figure_S2_data_distribution.pdf     # Data distribution
│   ├── Figure_S3_cshm_distribution.pdf     # CShM distribution
│   ├── Figure_S4_delta_heatmap.pdf         # ΔTop-1 heatmap
│   ├── Figure_S5_correction_cases.pdf      # Correction case studies
│   ├── Figure_S6_tokenizer_analysis.pdf    # Tokenizer analysis
│   ├── Figure_S7_corruption_results.pdf    # Corruption repair results
│   └── README.md                           # Figure descriptions
│
├── tables/                         # Tabular data (5 files)
│   ├── Table_S1_reference_polyhedra.csv    # CShM reference polyhedron coordinates
│   ├── Table_S2_donor_priority.csv         # Donor priority list
│   ├── Table_S3_data_filtering.csv         # Data filtering statistics
│   ├── Table_S4_hyperparameters.csv        # Training hyperparameters
│   └── Table_S5_corruption_suite.csv       # Synthetic corruption statistics
│
└── data/                           # Reproducible data (25+ files)
    ├── training_curves.csv         # Training curves
    ├── figure_s1_*.csv             # Figure S1 data
    ├── figure_s2*.csv              # Figure S2 data
    ├── figure_s3*.csv              # Figure S3 data
    ├── figure_s4*.csv              # Figure S4 data
    ├── figure_s5*.csv              # Figure S5 data
    ├── figure_s6*.csv              # Figure S6 data
    ├── figure_s7*.csv              # Figure S7 data
    └── panel_*.csv                 # Fig5 panel data
```

## SI Sections Overview

| Section | Title | Content |
|---------|-------|---------|
| S1 | Data Sources & Statistics | tmQM/COD versions, filtering rules, data flow |
| S2 | Conversion Pipeline | Structure → CoordRep full pipeline |
| S3 | CShM Definitions | Mathematical definition, reference polyhedron coordinates |
| S4 | Canonicalization | Ligand ordering rules, stereochemistry handling |
| S5 | Tokenizer Specification | Token types, grammar, known limitations |
| S6 | Masking Protocol | Tool A/B masking strategies |
| S7 | Synthetic Corruption | Corruption type definitions, dataset statistics |
| S8 | Baselines | Random/Freq/CondFreq definitions |
| S9 | Training Details | Hardware, hyperparameters, training curves |
| S10 | Reproducibility Package | Code structure, reproduction commands, data availability |

## Correspondence with Main Text

| Main Text Section | SI Reference |
|-------------------|-------------|
| 2.1 Data & Representation | S1, S2, S3 |
| 2.2 Canonicalization | S4 |
| 2.3 Language Model | S5, S6, S9 |
| 2.4 Tool A (Donor) | S6, S8, Fig5 B/C |
| 2.5 Tool B (Repair) | S7, Fig5 D/E |
| Methods | S9, S10 |

---

*Last updated: 2026-01-26*
