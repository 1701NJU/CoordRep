# Fig. 5 Reproducible Data

This directory contains the complete reproducible data for all Fig. 5 panels.

## Quick Reproduction

```bash
python plot_fig5_from_csv.py
```

Generated plots are saved to the `plots/` subdirectory.

---

## Data Files

### Panel B: ΔTop-1 Distribution (Beyond Conditional Lookup)

| File | Description |
|------|-------------|
| `panel_b_delta_histogram.csv` | Histogram bin data (17 bins) |
| `panel_b_delta_by_metal_cn.csv` | Full data for 179 (metal, CN) groups |
| `panel_b_top_gains.csv` | Top 20 groups with largest improvements |
| `panel_b_casecards.csv` | 30 cases where the model corrects CondFreq |
| `panel_b_summary.csv` | Statistical summary |

**Key columns (`panel_b_delta_by_metal_cn.csv`):**
- `metal`: Metal element
- `cn`: Coordination number
- `n_samples`: Number of samples
- `model_top1_acc`: Model Top-1 accuracy
- `condfreq_top1_acc`: Conditional frequency baseline Top-1 accuracy
- `delta_top1_pp`: ΔTop-1 (percentage points)

---

### Panel C: Confidence–Coverage Curve

| File | Description |
|------|-------------|
| `panel_c_confidence_coverage.csv` | Coverage/accuracy at 23 threshold points |
| `panel_c_donor_type_delta.csv` | Δ data for 7 donor types |

**Key columns (`panel_c_confidence_coverage.csv`):**
- `threshold`: Confidence threshold
- `coverage_pct`: Coverage (%)
- `accuracy_pct`: Accuracy (%)
- `n_samples`: Number of samples

**Key columns (`panel_c_donor_type_delta.csv`):**
- `donor_atom`: Donor atom type (C/N/O/P/F/S/I)
- `model_top1_acc_pct`: Model accuracy (%)
- `delta_top1_pp`: Δ vs CondFreq (percentage points)
- `x_random`: Fold-improvement vs Random baseline

---

### Panel D: Validity Uplift

| File | Description |
|------|-------------|
| `panel_d_validity_uplift.csv` | Before/after validity rates for 6 corruption types |

**Key columns:**
- `corruption_type`: Corruption type
- `valid_after_corruption`: Validity rate after corruption
- `valid_after_repair`: Validity rate after repair
- `uplift`: Uplift magnitude

---

### Panel E: CN-Stratified Structure Recovery

| File | Description |
|------|-------------|
| `panel_e_cn_stratified.csv` | Recovery rates for 13 CN values |
| `panel_e_summary.csv` | Overall statistics |

**Key columns:**
- `cn`: Coordination number
- `n_samples`: Number of samples
- `top1_acc_pct`: Top-1 recovery rate (%)
- `top5_acc_pct`: Top-5 recovery rate (%)

---

## Key Results Summary

### Panel B
- **Total groups**: 179 (metal, CN) combinations
- **Model wins**: 120 groups (67%)
- **CondFreq wins**: 30 groups (17%)
- **Mean ΔTop-1**: +4.26pp

### Panel C
- **τ = 0.85**: 95.6% accuracy @ 74.3% coverage
- **τ = 0.95**: 97.5% accuracy @ 61.9% coverage

### Panel D
- **Overall uplift**: +26.8pp
- **Missing bracket**: 0% → 58.8% (+58.8pp)
- **Truncation**: 36.8% → 100% (+63.3pp)

### Panel E
- **Overall**: 97.79% Top-1, 98.56% Top-5
- **CN = 6**: 97.98% (n = 2877)
- **CN = 8**: 99.27% (n = 1782)

---

## Dependencies

```
pandas
numpy
matplotlib
```
