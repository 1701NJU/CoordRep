# Supplementary Figures

## Figure List

| Figure | Filename | Content |
|--------|----------|---------|
| **S1** | `Figure_S1_training_curves.pdf` | Training curves (Loss, Perplexity, Validity, LR) |
| **S2** | `Figure_S2_data_distribution.pdf` | Data distribution (Metal, CN, Donor, Split) |
| **S3** | `Figure_S3_cshm_distribution.pdf` | CShM distribution (2D density plots for CN=4/5/6) |
| **S4** | `Figure_S4_delta_heatmap.pdf` | ΔTop-1 heatmap (30 metals × 9 CNs) |
| **S5** | `Figure_S5_correction_cases.pdf` | Correction case studies (20 representative cases) |
| **S6** | `Figure_S6_tokenizer_analysis.pdf` | Tokenizer analysis (Token types + sequence lengths) |
| **S7** | `Figure_S7_corruption_results.pdf` | Synthetic corruption repair results |

## Figure Descriptions

### Figure S1: Training Curves
- **(a)** Train/Val Loss vs Epoch
- **(b)** Validation Perplexity (log scale)
- **(c)** Syntax Validity (%) — annotated at ≥99.9% achievement point
- **(d)** Learning Rate Schedule (warmup + cosine decay)

### Figure S2: Data Distribution
- **(a)** Metal Element Distribution (grouped by 3d/4d/5d, log scale)
- **(b)** Coordination Number Distribution
- **(c)** Donor Atom Distribution
- **(d)** Train/Val/Test Split (80/10/10)

### Figure S3: CShM Distribution
- **(a)** CN=4: S(Td) vs S(SP) density plot
- **(b)** CN=5: S(TBP) vs S(SPY) density plot (Berry pathway annotated)
- **(c)** CN=6: S(Oh) vs S(TPr) density plot

### Figure S4: Extended ΔTop-1 Analysis
- 30 metals × 9 CNs heatmap
- Color scale: red (−40pp) → white (0) → blue (+40pp)
- Dividing lines for 3d/4d/5d

### Figure S5: Correction Case Studies
- 20 cases where Model is correct and CondFreq is incorrect
- Columns: Metal, CN, Target, CondFreq (wrong), Model (correct), Confidence

### Figure S6: Tokenizer Analysis
- **(a)** Token type distribution (7 categories, log scale)
- **(b)** Sequence length distribution (histogram)

### Figure S7: Synthetic Corruption Results
- **(a)** Before/After grouped bar chart
- **(b)** Uplift waterfall chart

## Reproduction

```bash
cd libcoordrep
python scripts/generate_si_figures.py
```

All figures are saved to this directory; data files are saved to `../data/`.

## Corresponding Data Files

| Figure | Data File(s) |
|--------|-------------|
| S1 | `figure_s1_training_curves.csv` |
| S2 | `figure_s2a_metal_distribution.csv`, `figure_s2b_cn_distribution.csv`, `figure_s2c_donor_distribution.csv`, `figure_s2d_split.csv` |
| S3 | `figure_s3a_cshm_cn4.csv`, `figure_s3b_cshm_cn5.csv`, `figure_s3c_cshm_cn6.csv` |
| S4 | `figure_s4_delta_full.csv`, `figure_s4_delta_matrix.csv` |
| S5 | `figure_s5_correction_cases.csv` |
| S6 | `figure_s6a_token_types.csv`, `figure_s6b_sequence_lengths.csv` |
| S7 | `figure_s7_corruption_results.csv` |
