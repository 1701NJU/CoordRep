# Evaluator QC Report

## Root Cause of Accuracy Discrepancy

The factorized ablation script (`factorized_token_ablation.py`)
and the standard Tool A evaluator use **fundamentally different
masking protocols**:

| Aspect | Standard Tool A | Factorized Ablation |
|--------|----------------|---------------------|
| Test set | fig5_tasks.jsonl (19,992 samples) | Random 2,000 from pipeline |
| Masked token | Single SMILES character (donor atom) | Donor marker (`:N:1` etc.) in constraint block |
| y_true type | Element char: C, N, O (from ~7 classes) | Full token (from 657+ vocab) |
| Candidates | ~7 donor elements | Entire vocabulary |
| Masks per sample | 1 | All donors simultaneously |
| Tokenization | Pre-tokenized with composite metal block | Re-tokenized at runtime |

## Impact on Reported Accuracy

- **Standard Tool A Top-1 = 85.7%**: Predicts which of ~7 elements
  (C, N, O, S, P, F, I) fills the masked SMILES position.
- **Factorized ablation Top-1 = 14.6–27.5%**: Predicts which of 657+
  vocabulary tokens fills the masked donor-marker position.

These numbers are **not comparable**. The factorized ablation metric
is valid for relative comparison (composite vs factorized), but the
absolute values cannot be compared with Tool A.

## Token Type at Masked Positions (fig5_tasks)

- **smiles_element**: 19992 (100.0%)

## y_true Distribution (fig5_tasks)

- **C**: 16031 (80.2%)
- **N**: 1586 (7.9%)
- **O**: 1579 (7.9%)
- **P**: 290 (1.5%)
- **F**: 256 (1.3%)
- **S**: 223 (1.1%)
- **I**: 27 (0.1%)

## Resolution

The standard Tool A evaluation must be re-run using:
1. The same fig5_tasks.jsonl test set
2. The same single-donor masking per sample
3. For factorized: re-tokenize the string, remap mask position
4. Restrict Top-k prediction to donor element tokens only
5. Stratify by denticity, CN, and ablation mode
