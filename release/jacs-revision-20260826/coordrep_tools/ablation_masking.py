#!/usr/bin/env python3
"""
ablation_masking.py
===================

Token-level masking functions for the legacy donor-field recovery diagnostic.

Each function operates on a **token list** (as produced by CoordRepTokenizer
or stored in fig5_tasks.jsonl ``tokens``/``masked_tokens`` fields) and
returns a new token list with the appropriate spans replaced by ``[MASK]``
or special placeholder tokens.

Ablation modes
--------------
A.  ``full_context``              – mask donor marker only
B.  ``no_ligand_smiles_keep_length`` – mask donor + SMILES chars, keep length
C.  ``no_ligand_smiles_collapsed``   – mask donor + collapse SMILES to [LIGMASK]
D.  ``no_geometry_stereo``        – mask shape + constraint blocks, keep SMILES
E.  ``metal_cn_only``             – keep metal/CN only; mask SMILES + geom + stereo
F.  ``shuffled_ligand_smiles_control`` – randomly swap ligand SMILES spans
"""

from __future__ import annotations

import random
import re
from copy import deepcopy
from typing import List, Optional, Tuple


# ──────────────────────────────────────────────────────────────────
# Token classification helpers
# ──────────────────────────────────────────────────────────────────

_METAL_RE = re.compile(r'^\[Metal:')
_SHAPE_START_TOKENS = {"<Shape:"}
_SHAPE_END_TOKENS = {">"}
_CONSTRAINT_START_TOKENS = {"{trans:", "{cis:"}
_CONSTRAINT_END_TOKENS = {"}"}
_LIG_ID_RE = re.compile(r'^L\d+$')
_DONOR_MARKER_RE = re.compile(r'^:[A-Z][a-z]?:\d+$')  # :N:1, :O:2, :Cl:1
_SMILES_ATOM_CHARS = set("BCNOPSFIHnoscpbfli0123456789()[]=#@+\\/-.")
_DELIMITER_TOKENS = {"|", "=", "--"}

_SPECIAL_TOKENS = {"[CLS]", "[SEP]", "[PAD]", "[MASK]", "[UNK]"}


def _is_metal_token(tok: str) -> bool:
    return bool(_METAL_RE.match(tok))


def _is_shape_token(tok: str) -> bool:
    """True for <Shape:, shape labels, V_xxx, =, ,  and >."""
    if tok in _SHAPE_START_TOKENS or tok in _SHAPE_END_TOKENS:
        return True
    if tok.startswith("V_"):
        return True
    if tok in {"Td", "SP", "Oh", "TP", "TBP", "SPY", "TPr", "L", ",", "="}:
        return True  # will be refined by context in block-aware functions
    return False


def _is_constraint_token(tok: str) -> bool:
    if tok in _CONSTRAINT_START_TOKENS or tok in _CONSTRAINT_END_TOKENS:
        return True
    if tok == "--":
        return True
    if _DONOR_MARKER_RE.match(tok):
        return True
    if _LIG_ID_RE.match(tok):
        return True  # L1:N:1 is split as "L1", ":", "N", ":", "1" or single token
    return False


def _is_donor_marker(tok: str) -> bool:
    return bool(_DONOR_MARKER_RE.match(tok))


def _is_lig_id(tok: str) -> bool:
    return bool(_LIG_ID_RE.match(tok))


# ──────────────────────────────────────────────────────────────────
# Block detection on token lists
# ──────────────────────────────────────────────────────────────────

def _detect_blocks(tokens: List[str]) -> dict:
    """
    Identify index ranges for metal, shape, constraint, and ligand-dict
    blocks in a token list.

    Returns dict with keys:
        metal_range:       (start, end) or None
        shape_range:       (start, end) or None
        constraint_ranges: list of (start, end)
        ligand_blocks:     list of {start, end, lid_pos, eq_pos, smiles_start, smiles_end}
    """
    result = {
        "metal_range": None,
        "shape_range": None,
        "constraint_ranges": [],
        "ligand_blocks": [],
    }

    n = len(tokens)

    # Metal: single token matching [Metal:...]
    for i, tok in enumerate(tokens):
        if _is_metal_token(tok):
            result["metal_range"] = (i, i + 1)
            break

    # Shape block: <Shape: ... >
    shape_start = None
    for i, tok in enumerate(tokens):
        if tok in _SHAPE_START_TOKENS:
            shape_start = i
        elif tok in _SHAPE_END_TOKENS and shape_start is not None:
            result["shape_range"] = (shape_start, i + 1)
            shape_start = None
            break

    # Constraint blocks: {trans: ... } or {cis: ... }
    con_start = None
    for i, tok in enumerate(tokens):
        if tok in _CONSTRAINT_START_TOKENS:
            con_start = i
        elif tok in _CONSTRAINT_END_TOKENS and con_start is not None:
            result["constraint_ranges"].append((con_start, i + 1))
            con_start = None

    # Ligand dict blocks: | Lid = <smiles chars> |
    # Handles both single-token ligand IDs (e.g. "L1") and split
    # ligand IDs (e.g. "L", "2") which occur in some tokenizations.
    # NOTE: the closing "|" of one block may be the opening "|" of
    # the next, so after consuming a block we re-check the closing
    # "|" position rather than skipping past it.
    i = 0
    while i < n:
        if tokens[i] == "|" and i + 1 < n:
            # Case 1: single-token Lid  →  | L1 = ...
            # Case 2: split Lid         →  | L 2 = ...
            lid_pos = None
            lid_end = None  # index after the Lid token(s)
            if _is_lig_id(tokens[i + 1]):
                lid_pos = i + 1
                lid_end = i + 2
            elif (tokens[i + 1] == "L" and i + 2 < n
                  and tokens[i + 2].isdigit()):
                lid_pos = i + 1
                lid_end = i + 3  # skip "L" and digit

            if lid_pos is not None:
                block_start = i
                eq_pos = None
                smiles_start = None
                j = lid_end
                # find '='
                if j < n and tokens[j] == "=":
                    eq_pos = j
                    smiles_start = j + 1
                # find closing '|'
                j = smiles_start if smiles_start else j
                while j < n and tokens[j] != "|":
                    j += 1
                smiles_end = j  # exclusive

                result["ligand_blocks"].append({
                    "start": block_start,
                    "end": j + 1 if j < n else j,
                    "lid_pos": lid_pos,
                    "eq_pos": eq_pos,
                    "smiles_start": smiles_start,
                    "smiles_end": smiles_end,
                })
                # Re-check the closing "|" — it may be the opening
                # "|" of the next ligand block.
                i = j
                continue
        i += 1

    return result


# ──────────────────────────────────────────────────────────────────
# Core masking functions
# ──────────────────────────────────────────────────────────────────

def mask_donor_markers(tokens: List[str], mask_token: str = "[MASK]") -> Tuple[List[str], List[int]]:
    """
    Mask donor-marker tokens (e.g. ``:N:1``, ``:O:2``).

    Returns (masked_tokens, mask_positions).
    """
    out = tokens[:]
    positions = []
    for i, tok in enumerate(out):
        if _is_donor_marker(tok):
            out[i] = mask_token
            positions.append(i)
    return out, positions


def mask_ligand_smiles(
    tokens: List[str],
    mode: str = "keep_length",
    mask_token: str = "[MASK]",
    collapsed_token: str = "[LIGMASK]",
) -> List[str]:
    """
    Mask ligand SMILES characters within ligand-dict blocks.

    mode='keep_length':  replace each SMILES char with [MASK], preserving length
    mode='collapsed':    replace entire SMILES span with single [LIGMASK]
    """
    blocks = _detect_blocks(tokens)
    out = tokens[:]

    if mode == "keep_length":
        for lb in blocks["ligand_blocks"]:
            for j in range(lb["smiles_start"], lb["smiles_end"]):
                out[j] = mask_token
        return out

    elif mode == "collapsed":
        # Build output with collapsed spans; indices shift
        result = []
        skip_until = -1
        for i, tok in enumerate(tokens):
            if i < skip_until:
                continue
            # Check if this is the smiles_start of a ligand block
            in_lig = False
            for lb in blocks["ligand_blocks"]:
                if i == lb["smiles_start"]:
                    result.append(collapsed_token)
                    skip_until = lb["smiles_end"]
                    in_lig = True
                    break
            if not in_lig:
                result.append(tok)
        return result

    else:
        raise ValueError(f"Unknown mode: {mode}")


def mask_geometry_stereo(tokens: List[str], mask_token: str = "[MASK]") -> List[str]:
    """
    Mask shape block and constraint blocks, preserving ligand SMILES.
    """
    blocks = _detect_blocks(tokens)
    out = tokens[:]

    if blocks["shape_range"]:
        s, e = blocks["shape_range"]
        for j in range(s, e):
            out[j] = mask_token

    for cs, ce in blocks["constraint_ranges"]:
        for j in range(cs, ce):
            out[j] = mask_token

    return out


def shuffle_ligand_smiles(
    tokens: List[str],
    seed: int = None,
) -> List[str]:
    """
    Randomly permute ligand SMILES spans while keeping metal/geometry/stereo
    and ligand slot structure intact.
    """
    blocks = _detect_blocks(tokens)
    if len(blocks["ligand_blocks"]) < 2:
        return tokens[:]

    # Extract SMILES sub-lists
    smiles_spans = []
    for lb in blocks["ligand_blocks"]:
        smiles_spans.append(tokens[lb["smiles_start"]:lb["smiles_end"]])

    # Shuffle
    rng = random.Random(seed)
    indices = list(range(len(smiles_spans)))
    rng.shuffle(indices)
    shuffled_spans = [smiles_spans[i] for i in indices]

    # Reconstruct
    out = tokens[:]
    for lb, new_smiles in zip(blocks["ligand_blocks"], shuffled_spans):
        s, e = lb["smiles_start"], lb["smiles_end"]
        # Replace in-place (lengths may differ — rebuild)
        pass

    # Safer: full rebuild
    result = []
    lb_idx = 0
    i = 0
    while i < len(tokens):
        if lb_idx < len(blocks["ligand_blocks"]):
            lb = blocks["ligand_blocks"][lb_idx]
            if i == lb["smiles_start"]:
                result.extend(shuffled_spans[lb_idx])
                i = lb["smiles_end"]
                lb_idx += 1
                continue
        result.append(tokens[i])
        i += 1
    return result


# ──────────────────────────────────────────────────────────────────
# Composite ablation builder
# ──────────────────────────────────────────────────────────────────

ABLATION_MODES = [
    "full_context",
    "no_ligand_smiles_keep_length",
    "no_ligand_smiles_collapsed",
    "no_geometry_stereo",
    "metal_cn_only",
    "shuffled_ligand_smiles_control",
]


def make_ablation_input(
    tokens: List[str],
    donor_mask_positions: List[int],
    ablation_mode: str,
    mask_token: str = "[MASK]",
    collapsed_token: str = "[LIGMASK]",
    shuffle_seed: int = None,
) -> Tuple[List[str], List[int]]:
    """
    Build an ablated token list for a given mode.

    Parameters
    ----------
    tokens : list[str]
        Original (unmasked) token list.
    donor_mask_positions : list[int]
        Positions of donor-marker tokens to predict.
    ablation_mode : str
        One of ABLATION_MODES.
    mask_token, collapsed_token : str
        Replacement tokens.
    shuffle_seed : int, optional
        Random seed for shuffled control.

    Returns
    -------
    (ablated_tokens, new_mask_positions)
        new_mask_positions accounts for length changes in collapsed mode.
    """
    if ablation_mode not in ABLATION_MODES:
        raise ValueError(f"Unknown ablation_mode={ablation_mode!r}. "
                         f"Choose from {ABLATION_MODES}")

    out = tokens[:]

    if ablation_mode == "full_context":
        # Mode A: mask donor markers only (current Tool A)
        for p in donor_mask_positions:
            if p < len(out):
                out[p] = mask_token
        return out, donor_mask_positions

    elif ablation_mode == "no_ligand_smiles_keep_length":
        # Mode B: mask donor markers + mask SMILES chars (keep length)
        for p in donor_mask_positions:
            if p < len(out):
                out[p] = mask_token
        out = mask_ligand_smiles(out, mode="keep_length", mask_token=mask_token)
        return out, donor_mask_positions

    elif ablation_mode == "no_ligand_smiles_collapsed":
        # Mode C: mask donor markers + collapse SMILES
        # First mark donors, then collapse
        for p in donor_mask_positions:
            if p < len(out):
                out[p] = mask_token

        blocks = _detect_blocks(out)
        result = []
        new_positions = []
        skip_until = -1
        position_map = {}  # old_idx -> new_idx
        for i, tok in enumerate(out):
            if i < skip_until:
                continue
            in_lig = False
            for lb in blocks["ligand_blocks"]:
                if i == lb["smiles_start"]:
                    # Check if any donor mask positions fall in this span
                    for p in donor_mask_positions:
                        if lb["smiles_start"] <= p < lb["smiles_end"]:
                            # Donor was inside SMILES — put a MASK then LIGMASK
                            new_positions.append(len(result))
                            result.append(mask_token)
                    result.append(collapsed_token)
                    skip_until = lb["smiles_end"]
                    in_lig = True
                    break
            if not in_lig:
                if tok == mask_token and i in donor_mask_positions:
                    new_positions.append(len(result))
                result.append(tok)

        # Deduplicate
        if not new_positions:
            # Fallback: find all [MASK] in result
            new_positions = [j for j, t in enumerate(result) if t == mask_token]

        return result, new_positions

    elif ablation_mode == "no_geometry_stereo":
        # Mode D: mask shape + constraints, keep SMILES
        out = mask_geometry_stereo(out, mask_token=mask_token)
        for p in donor_mask_positions:
            if p < len(out):
                out[p] = mask_token
        return out, donor_mask_positions

    elif ablation_mode == "metal_cn_only":
        # Mode E: keep metal/CN + ligand slot structure only
        # Mask: SMILES chars + geometry + stereo
        out = mask_geometry_stereo(out, mask_token=mask_token)
        out = mask_ligand_smiles(out, mode="keep_length", mask_token=mask_token)
        for p in donor_mask_positions:
            if p < len(out):
                out[p] = mask_token
        return out, donor_mask_positions

    elif ablation_mode == "shuffled_ligand_smiles_control":
        # Mode F: shuffle SMILES spans, then mask donor markers
        out = shuffle_ligand_smiles(out, seed=shuffle_seed)
        for p in donor_mask_positions:
            if p < len(out):
                out[p] = mask_token
        return out, donor_mask_positions

    else:
        raise ValueError(f"Unknown ablation_mode: {ablation_mode}")


# ──────────────────────────────────────────────────────────────────
# Convenience: apply ablation to a fig5_tasks sample dict
# ──────────────────────────────────────────────────────────────────

def ablate_sample(
    sample: dict,
    ablation_mode: str,
    shuffle_seed: int = None,
) -> dict:
    """
    Take a ``fig5_tasks.jsonl`` sample (with ``tokens``, ``mask_pos``,
    ``masked_tokens``) and return a new sample with ablated tokens.
    """
    tokens = sample["tokens"]
    mask_pos = sample.get("mask_pos")

    if mask_pos is None or (isinstance(mask_pos, list) and len(mask_pos) == 0):
        # Infer from masked_tokens
        mtoks = sample.get("masked_tokens", [])
        mask_pos_list = [i for i in range(min(len(tokens), len(mtoks)))
                         if tokens[i] != mtoks[i]]
    elif isinstance(mask_pos, int):
        mask_pos_list = [mask_pos]
    else:
        mask_pos_list = list(mask_pos)

    ablated, new_pos = make_ablation_input(
        tokens, mask_pos_list, ablation_mode,
        shuffle_seed=shuffle_seed,
    )

    new_sample = dict(sample)
    new_sample["ablated_tokens"] = ablated
    new_sample["ablated_mask_positions"] = new_pos
    new_sample["ablation_mode"] = ablation_mode
    return new_sample
