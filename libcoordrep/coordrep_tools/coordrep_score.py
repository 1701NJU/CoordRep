"""
coordrep_score.py
=================
Pseudo-log-likelihood (PLL) scoring of CoordRep strings using the
pretrained MLM.  Two strategies:

* **mask-one-token** (token-PLL): mask each non-special token once,
  sum log P(t_i | t_{\\i}).  Slow but gold-standard.
* **mask-field**: mask all tokens in a structural field at once and
  sum log probs.  Fast, also gives *field-wise* decomposition.

Public API
----------
compute_pseudo_log_likelihood(model, tokenizer, tokens, device, mode)
compute_fieldwise_scores(model, tokenizer, tokens, device)
batch_pll(model, tokenizer, token_lists, device, batch_size, mode)
"""

from __future__ import annotations

import math
import re
from typing import Dict, List, Optional, Tuple

import torch
import numpy as np

# ── field detection (reuse ablation_masking helpers) ─────────────
_METAL_RE = re.compile(r'^\[Metal:')
_LIG_ID_RE = re.compile(r'^L\d+$')
_SPECIAL = {"[CLS]", "[SEP]", "[PAD]", "[MASK]", "[UNK]"}

FIELD_NAMES = ["metal", "shape", "stereo", "ligand", "donor", "full"]


def _classify_token_fields(tokens: List[str]) -> Dict[str, List[int]]:
    """Return {field_name: [token indices]} for the major CoordRep fields."""
    n = len(tokens)
    fields: Dict[str, List[int]] = {
        "metal": [], "shape": [], "stereo": [],
        "ligand_smiles": [], "donor_marker": [],
        "delimiter": [], "other": [],
    }

    i = 0
    while i < n:
        tok = tokens[i]

        # Metal token
        if _METAL_RE.match(tok):
            fields["metal"].append(i)
            i += 1
            continue

        # Shape block:  <Shape: ... >
        if tok == "<Shape:":
            fields["shape"].append(i)
            i += 1
            while i < n and tokens[i] != ">":
                fields["shape"].append(i)
                i += 1
            if i < n:
                fields["shape"].append(i)  # the ">"
            i += 1
            continue

        # Stereo block:  {trans: ... } or {cis: ... }
        if tok in ("{trans:", "{cis:"):
            fields["stereo"].append(i)
            i += 1
            while i < n and tokens[i] != "}":
                fields["stereo"].append(i)
                i += 1
            if i < n:
                fields["stereo"].append(i)  # the "}"
            i += 1
            continue

        # Ligand block:  | Lid = SMILES ... |
        if tok == "|":
            fields["delimiter"].append(i)
            # Check if next is a ligand ID
            nxt = i + 1
            is_lig_start = False
            if nxt < n:
                if _LIG_ID_RE.match(tokens[nxt]):
                    is_lig_start = True
                elif tokens[nxt] == "L" and nxt + 1 < n and tokens[nxt + 1].isdigit():
                    is_lig_start = True
            if is_lig_start:
                i += 1
                # Lid token(s)
                if _LIG_ID_RE.match(tokens[i]):
                    fields["other"].append(i)  # Lid itself
                    i += 1
                else:
                    # split: "L" then digit
                    fields["other"].append(i)
                    i += 1
                    if i < n and tokens[i].isdigit():
                        fields["other"].append(i)
                        i += 1
                # = sign
                if i < n and tokens[i] == "=":
                    fields["delimiter"].append(i)
                    i += 1
                # SMILES tokens until next |
                while i < n and tokens[i] != "|":
                    # donor markers :X:n
                    if re.match(r'^:[A-Z][a-z]?:\d+$', tokens[i]):
                        fields["donor_marker"].append(i)
                    else:
                        fields["ligand_smiles"].append(i)
                    i += 1
                # Don't consume the closing |; loop will pick it up
                continue
            else:
                i += 1
                continue

        # Special / other
        if tok in _SPECIAL:
            fields["other"].append(i)
        elif tok in ("=", "--"):
            fields["delimiter"].append(i)
        else:
            fields["other"].append(i)
        i += 1

    return fields


# ── token → id conversion ───────────────────────────────────────

def _tokens_to_ids(tokens, tokenizer):
    mask_id = tokenizer.token2id.get("[MASK]", 3)
    unk_id = tokenizer.token2id.get("[UNK]", 1)
    ids = []
    for t in tokens:
        if t == "[MASK]":
            ids.append(mask_id)
        elif t in tokenizer.token2id:
            ids.append(tokenizer.token2id[t])
        else:
            ids.append(unk_id)
    return ids


# ── core PLL ─────────────────────────────────────────────────────

@torch.no_grad()
def compute_pseudo_log_likelihood(
    model,
    tokenizer,
    tokens: List[str],
    device: str = "cuda",
    mode: str = "mask_one_token",
    positions: Optional[List[int]] = None,
) -> float:
    """
    Compute pseudo-log-likelihood of *tokens*.

    Parameters
    ----------
    mode : str
        ``"mask_one_token"`` – mask each scorable position individually.
        ``"mask_field"`` – mask all *positions* at once and score.
    positions : list[int] or None
        If ``mode="mask_field"``, the token indices to mask.
        If ``None`` with ``mask_one_token``, score all non-special tokens.

    Returns
    -------
    float  – sum of log-probabilities (higher = more compatible).
    """
    mask_id = tokenizer.token2id.get("[MASK]", 3)
    unk_id = tokenizer.token2id.get("[UNK]", 1)
    base_ids = _tokens_to_ids(tokens, tokenizer)
    max_len = 512

    if mode == "mask_one_token":
        # Determine which positions to score
        if positions is None:
            positions = [i for i, t in enumerate(tokens)
                         if t not in _SPECIAL and i < max_len]
        if not positions:
            return 0.0

        # Batch: create one copy per position with that position masked
        batch_ids = []
        for pos in positions:
            ids_copy = list(base_ids[:max_len])
            if pos < len(ids_copy):
                ids_copy[pos] = mask_id
            batch_ids.append(ids_copy)

        # Pad and run in sub-batches
        pll = 0.0
        sub_bs = 128
        for start in range(0, len(batch_ids), sub_bs):
            end = min(start + sub_bs, len(batch_ids))
            chunk = batch_ids[start:end]
            chunk_pos = positions[start:end]
            ml = max(len(ids) for ids in chunk)
            padded = [ids + [0] * (ml - len(ids)) for ids in chunk]
            masks = [[1] * len(ids) + [0] * (ml - len(ids)) for ids in chunk]

            inp = torch.tensor(padded, device=device)
            att = torch.tensor(masks, device=device)
            logits, _ = model(inp, attention_mask=att)

            for j, pos in enumerate(chunk_pos):
                if pos >= logits.shape[1]:
                    continue
                probs = torch.softmax(logits[j, pos, :], dim=-1)
                true_id = base_ids[pos] if pos < len(base_ids) else unk_id
                p = probs[true_id].item()
                pll += math.log(max(p, 1e-20))

        return pll

    elif mode == "mask_chunked":
        # Partition positions into chunks of ~15% of tokens.
        # Each chunk is scored with the rest visible. ~7 forward passes.
        if positions is None:
            positions = [i for i, t in enumerate(tokens)
                         if t not in _SPECIAL and i < max_len]
        if not positions:
            return 0.0
        n_chunks = 7
        chunk_size = max(1, len(positions) // n_chunks)
        # Shuffle positions for even distribution, then partition
        pos_list = list(positions)
        chunks = [pos_list[i:i + chunk_size]
                  for i in range(0, len(pos_list), chunk_size)]

        pll = 0.0
        for chunk in chunks:
            ids_copy = list(base_ids[:max_len])
            valid_chunk = [p for p in chunk if p < len(ids_copy)]
            for p in valid_chunk:
                ids_copy[p] = mask_id
            if not valid_chunk:
                continue
            inp = torch.tensor([ids_copy], device=device)
            att = torch.ones_like(inp)
            logits, _ = model(inp, attention_mask=att)
            for p in valid_chunk:
                if p >= logits.shape[1]:
                    continue
                probs = torch.softmax(logits[0, p, :], dim=-1)
                true_id = base_ids[p] if p < len(base_ids) else unk_id
                prob = probs[true_id].item()
                pll += math.log(max(prob, 1e-20))
        return pll

    elif mode == "mask_all":
        # Approximate PLL: mask ALL scorable tokens at once.
        # Single forward pass — fast but over-estimates uncertainty.
        if positions is None:
            positions = [i for i, t in enumerate(tokens)
                         if t not in _SPECIAL and i < max_len]
        if not positions:
            return 0.0
        ids_copy = list(base_ids[:max_len])
        valid_pos = [p for p in positions if p < len(ids_copy)]
        for p in valid_pos:
            ids_copy[p] = mask_id
        inp = torch.tensor([ids_copy], device=device)
        att = torch.ones_like(inp)
        logits, _ = model(inp, attention_mask=att)
        pll = 0.0
        for p in valid_pos:
            if p >= logits.shape[1]:
                continue
            probs = torch.softmax(logits[0, p, :], dim=-1)
            true_id = base_ids[p] if p < len(base_ids) else unk_id
            prob = probs[true_id].item()
            pll += math.log(max(prob, 1e-20))
        return pll

    elif mode == "mask_field":
        if positions is None:
            return 0.0
        ids_copy = list(base_ids[:max_len])
        valid_pos = [p for p in positions if p < len(ids_copy)]
        for p in valid_pos:
            ids_copy[p] = mask_id
        if not valid_pos:
            return 0.0

        inp = torch.tensor([ids_copy], device=device)
        att = torch.ones_like(inp)
        logits, _ = model(inp, attention_mask=att)

        pll = 0.0
        for p in valid_pos:
            if p >= logits.shape[1]:
                continue
            probs = torch.softmax(logits[0, p, :], dim=-1)
            true_id = base_ids[p] if p < len(base_ids) else unk_id
            prob = probs[true_id].item()
            pll += math.log(max(prob, 1e-20))
        return pll

    else:
        raise ValueError(f"Unknown mode: {mode}")


def compute_fieldwise_scores(
    model,
    tokenizer,
    tokens: List[str],
    device: str = "cuda",
    mode: str = "mask_field",
) -> Dict[str, float]:
    """
    Decompose PLL into per-field contributions.

    Returns dict with keys: metal, shape, stereo, ligand_smiles,
    donor_marker, full, n_tokens.
    """
    fields = _classify_token_fields(tokens)

    scores = {}
    total = 0.0
    for fname in ("metal", "shape", "stereo", "ligand_smiles", "donor_marker"):
        positions = fields.get(fname, [])
        if positions:
            s = compute_pseudo_log_likelihood(
                model, tokenizer, tokens, device, mode=mode, positions=positions,
            )
        else:
            s = 0.0
        scores[fname] = s
        total += s

    # "full" = score everything scorable
    all_pos = []
    for fname in ("metal", "shape", "stereo", "ligand_smiles",
                   "donor_marker", "delimiter", "other"):
        all_pos.extend(fields.get(fname, []))
    all_pos = [p for p in all_pos if tokens[p] not in _SPECIAL]

    if mode == "mask_field":
        # For full, use mask-one-token to avoid masking everything
        scores["full"] = compute_pseudo_log_likelihood(
            model, tokenizer, tokens, device,
            mode="mask_one_token", positions=sorted(all_pos),
        )
    else:
        scores["full"] = compute_pseudo_log_likelihood(
            model, tokenizer, tokens, device,
            mode=mode, positions=sorted(all_pos),
        )

    scores["n_tokens"] = len(tokens)
    scores["n_scored"] = len(all_pos)
    return scores


# ── batched convenience ──────────────────────────────────────────

def compute_normalized_pll(
    model,
    tokenizer,
    tokens: List[str],
    device: str = "cuda",
    mode: str = "mask_one_token",
) -> Tuple[float, float, int]:
    """
    Compute both raw and length-normalized PLL.

    Returns (raw_pll, per_token_pll, n_scored).
    """
    positions = [i for i, t in enumerate(tokens)
                 if t not in _SPECIAL and i < 512]
    n_scored = len(positions)
    if n_scored == 0:
        return 0.0, 0.0, 0
    raw = compute_pseudo_log_likelihood(
        model, tokenizer, tokens, device, mode=mode, positions=positions,
    )
    return raw, raw / n_scored, n_scored


def batch_pll(
    model,
    tokenizer,
    token_lists: List[List[str]],
    device: str = "cuda",
    batch_size: int = 32,
    mode: str = "mask_one_token",
    normalize: bool = True,
) -> List[float]:
    """Compute PLL for many sequences. Returns list of scores (normalized if requested)."""
    scores = []
    for tl in token_lists:
        raw, normed, n = compute_normalized_pll(model, tokenizer, tl, device, mode=mode)
        scores.append(normed if normalize else raw)
    return scores


# ── conditional field scoring ────────────────────────────────────

@torch.no_grad()
def score_field(
    model,
    tokenizer,
    tokens: List[str],
    field: str,
    device: str = "cuda",
) -> Dict[str, float]:
    """
    Mask a target field and compute log P(true_field | remaining context).

    Parameters
    ----------
    field : str
        One of "metal", "shape", "stereo", "ligand_smiles", "donor_marker".

    Returns
    -------
    dict with keys:
        raw   – sum of log P(token_i | context) for masked positions
        norm  – raw / n_positions
        n     – number of masked positions
    """
    fields = _classify_token_fields(tokens)
    positions = fields.get(field, [])
    if not positions:
        return {"raw": 0.0, "norm": 0.0, "n": 0}
    positions = [p for p in positions if p < 512]
    if not positions:
        return {"raw": 0.0, "norm": 0.0, "n": 0}

    raw = compute_pseudo_log_likelihood(
        model, tokenizer, tokens, device,
        mode="mask_field", positions=positions,
    )
    n = len(positions)
    return {"raw": raw, "norm": raw / n if n > 0 else 0.0, "n": n}


@torch.no_grad()
def score_field_marginal(
    model,
    tokenizer,
    tokens: List[str],
    field: str,
    device: str = "cuda",
) -> float:
    """
    Compute marginal log P(field) by masking ALL other tokens (context-free).

    This gives the prior probability of the field tokens, used for PMI.
    """
    fields = _classify_token_fields(tokens)
    target_pos = set(fields.get(field, []))
    if not target_pos:
        return 0.0

    mask_id = tokenizer.token2id.get("[MASK]", 3)
    unk_id = tokenizer.token2id.get("[UNK]", 1)
    base_ids = _tokens_to_ids(tokens, tokenizer)
    max_len = 512

    # Mask everything EXCEPT the target field and special tokens
    ids_copy = list(base_ids[:max_len])
    for i in range(len(ids_copy)):
        t = tokens[i] if i < len(tokens) else "[PAD]"
        if i not in target_pos and t not in _SPECIAL:
            ids_copy[i] = mask_id

    # Now mask target field too and score
    ids_masked = list(ids_copy)
    valid_pos = [p for p in target_pos if p < len(ids_masked)]
    for p in valid_pos:
        ids_masked[p] = mask_id

    if not valid_pos:
        return 0.0

    inp = torch.tensor([ids_masked], device=device)
    att = torch.ones_like(inp)
    logits, _ = model(inp, attention_mask=att)

    pll = 0.0
    for p in valid_pos:
        if p >= logits.shape[1]:
            continue
        probs = torch.softmax(logits[0, p, :], dim=-1)
        true_id = base_ids[p] if p < len(base_ids) else unk_id
        pll += math.log(max(probs[true_id].item(), 1e-20))
    return pll


def compute_pmi_score(
    model,
    tokenizer,
    tokens: List[str],
    field: str,
    device: str = "cuda",
) -> Dict[str, float]:
    """
    PMI-style score: log P(field | context) - log P(field).

    Positive PMI means the field is more likely given context than alone,
    indicating learned compatibility.

    Returns dict with raw, norm, pmi, n.
    """
    cond = score_field(model, tokenizer, tokens, field, device)
    marginal = score_field_marginal(model, tokenizer, tokens, field, device)
    n = cond["n"]
    pmi = cond["raw"] - marginal if n > 0 else 0.0
    pmi_norm = pmi / n if n > 0 else 0.0
    return {
        "raw": cond["raw"],
        "norm": cond["norm"],
        "marginal": marginal,
        "pmi": pmi,
        "pmi_norm": pmi_norm,
        "n": n,
    }
