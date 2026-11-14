# Tokenizer Fix: Composite → Factorized

## Problem Identified

The CoordRep serializer emits metal blocks using **pipe** separators:
```
[Metal:Fe|ox:+2|d:d6|CN:6]
```

The tokenizer was designed to split on **semicolons**:
```
[Fe;ox=+2;d=6;CN=6] → [Fe] ;ox=+2 ;d=6 ;CN=6
```

Because of this mismatch, the entire pipe-separated block was treated as a **single composite token** (191 fused metal+CN tokens in the vocabulary). This prevented the model from independently learning metal, oxidation state, and coordination number.

## Fix Applied

`brain/tokenizer.py` → `_tokenize_metal()` now detects both formats and **always emits factorized tokens**:

```python
# Input (pipe format):  [Metal:Fe|ox:+2|d:d6|CN:6]
# Input (semicolon):    [Fe;ox=+2;d=6;CN=6]
# Output (always):      ["[Fe]", ";ox=+2", ";d=6", ";CN=6"]
```

The legacy composite tokenizer is preserved in `coordrep_tools/tokenizer_variants.py` for ablation comparison only.

## Ablation Results (Matched Short Training)

| Task | Composite | Factorized | Improvement |
|------|-----------|------------|-------------|
| Metal Top-1 | 13.8% | 67.8% | +54.0 pts |
| CN Top-1 | 13.8% | 49.0% | +35.3 pts |
| Joint Metal+CN Top-1 | 13.8% | 29.4% | +15.6 pts |
| Donor Top-1 | 15.1% | 27.5% | +12.4 pts |

## Caveat

These numbers are from matched short training (same epochs, same data). They are **not directly comparable** to the fully pretrained Tool A model. The main conclusion is that factorized tokenization enables independent field recovery and improves tokenizer interpretability.

## Revised Default

The factorized tokenizer is now the **default** for all new training and inference. The composite variant is retained only as an ablation baseline.
