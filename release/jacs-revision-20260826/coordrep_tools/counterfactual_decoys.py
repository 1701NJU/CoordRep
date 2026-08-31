"""
counterfactual_decoys.py
========================
Generate syntactically valid but chemically counterfactual CoordRep
decoys by swapping individual fields.

Decoy types
-----------
1. **metal_swap** – replace metal element (and optionally CN)
2. **shape_swap** – replace shape label with another valid one
3. **stereo_swap** – flip trans↔cis or add/remove constraints
4. **ligand_swap** – replace one or more ligand SMILES blocks
5. **donor_pattern_swap** – permute donor markers within a ligand
6. **boundary** – near-valid perturbations (CN ±1, single-atom donor edit)

Each generator takes a token list and a pool of donor materials
(other complexes' tokens) and returns a list of decoy token lists.
Every decoy is checked with `is_valid_coordrep` before inclusion.
"""

from __future__ import annotations

import random
import re
from copy import deepcopy
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

from .validate import is_valid_coordrep
from .ablation_masking import _detect_blocks

# ── helpers ──────────────────────────────────────────────────────

_METAL_RE = re.compile(r'^\[Metal:(\w+)\|.*CN:(\d+)\]$')
_LIG_ID_RE = re.compile(r'^L\d+$')
_DONOR_RE = re.compile(r'^:[A-Z][a-z]?:\d+$')


def _parse_metal_token(tok: str) -> Optional[Tuple[str, int]]:
    m = _METAL_RE.match(tok)
    if m:
        return m.group(1), int(m.group(2))
    return None


def _rebuild_metal_token(tok: str, new_metal: str = None,
                         new_cn: int = None) -> str:
    """Replace metal element and/or CN in a metal token string."""
    s = tok
    if new_metal:
        s = re.sub(r'(\[Metal:)\w+', rf'\g<1>{new_metal}', s)
    if new_cn is not None:
        s = re.sub(r'CN:\d+', f'CN:{new_cn}', s)
    return s


def _tokens_to_string(tokens: List[str]) -> str:
    return "".join(tokens)


def _valid(tokens: List[str]) -> bool:
    s = _tokens_to_string(tokens)
    return is_valid_coordrep(s)


def _extract_ligand_block_tokens(tokens, block_info):
    """Extract the SMILES tokens from a ligand block."""
    return tokens[block_info["smiles_start"]:block_info["smiles_end"]]


# ── pool builder ─────────────────────────────────────────────────

class DecoyPool:
    """Pre-compute swap materials from a corpus of real complexes."""

    def __init__(self, all_token_lists: List[List[str]], seed: int = 42):
        self.rng = random.Random(seed)
        self.metals: List[str] = []
        self.metal_tokens: List[str] = []  # full [Metal:...|CN:...] tokens
        self.ligand_blocks: List[List[str]] = []  # SMILES token lists
        self.stereo_blocks: List[List[str]] = []  # full stereo token spans
        self.cn_values: List[int] = []

        metal_set = set()
        for toks in all_token_lists:
            blocks = _detect_blocks(toks)

            # Metals
            if blocks["metal_range"] is not None:
                mt = toks[blocks["metal_range"][0]]
                parsed = _parse_metal_token(mt)
                if parsed:
                    metal_set.add(parsed[0])
                    self.metal_tokens.append(mt)
                    self.cn_values.append(parsed[1])

            # Ligand SMILES
            for lb in blocks["ligand_blocks"]:
                if lb.get("smiles_start") is None or lb.get("smiles_end") is None:
                    continue
                smi_toks = toks[lb["smiles_start"]:lb["smiles_end"]]
                if smi_toks:
                    self.ligand_blocks.append(smi_toks)

            # Stereo
            for sr in blocks.get("constraint_ranges", []):
                st_toks = toks[sr[0]:sr[1]]
                if st_toks:
                    self.stereo_blocks.append(st_toks)

        self.metals = sorted(metal_set)

    def random_metal(self, exclude: str = None) -> str:
        choices = [m for m in self.metals if m != exclude]
        return self.rng.choice(choices) if choices else exclude

    def random_metal_token(self, exclude_metal: str = None) -> str:
        choices = [t for t in self.metal_tokens
                   if _parse_metal_token(t)[0] != exclude_metal]
        return self.rng.choice(choices) if choices else self.metal_tokens[0]

    def random_ligand_smiles(self, exclude: List[str] = None) -> List[str]:
        if exclude:
            ex_str = "".join(exclude)
            choices = [lb for lb in self.ligand_blocks
                       if "".join(lb) != ex_str]
        else:
            choices = self.ligand_blocks
        return list(self.rng.choice(choices)) if choices else list(exclude or [])

    def random_cn(self, exclude: int = None) -> int:
        choices = [c for c in self.cn_values if c != exclude]
        return self.rng.choice(choices) if choices else (exclude or 6)


# ── decoy generators ────────────────────────────────────────────

def generate_metal_swap_decoys(
    tokens: List[str],
    pool: DecoyPool,
    n: int = 10,
) -> List[Tuple[List[str], str]]:
    """Swap the metal element, keeping CN. Returns [(decoy_tokens, reason)]."""
    blocks = _detect_blocks(tokens)
    if blocks["metal_range"] is None:
        return []

    mp = blocks["metal_range"][0]
    orig_tok = tokens[mp]
    parsed = _parse_metal_token(orig_tok)
    if not parsed:
        return []
    orig_metal, orig_cn = parsed

    decoys = []
    tried = set()
    for _ in range(n * 3):
        new_metal = pool.random_metal(exclude=orig_metal)
        if new_metal in tried or new_metal == orig_metal:
            continue
        tried.add(new_metal)

        new_tok = _rebuild_metal_token(orig_tok, new_metal=new_metal)
        d = list(tokens)
        d[mp] = new_tok
        if _valid(d):
            decoys.append((d, f"metal_swap:{orig_metal}->{new_metal}"))
        if len(decoys) >= n:
            break

    return decoys[:n]


def generate_shape_swap_decoys(
    tokens: List[str],
    pool: DecoyPool,
    n: int = 5,
) -> List[Tuple[List[str], str]]:
    """Swap shape region content. Only for complexes with non-empty shapes."""
    blocks = _detect_blocks(tokens)
    sr = blocks.get("shape_range")
    if not sr:
        return []

    s_start, s_end = sr
    orig_shape = tokens[s_start:s_end]

    # Skip empty-shape complexes: inserting shape is not a counterfactual
    if len(orig_shape) <= 2:
        return []

    decoys = []
    # Decoy: clear the shape
    d = list(tokens)
    d[s_start:s_end] = ["<Shape:", ">"]
    if _valid(d):
        decoys.append((d, "shape_swap:clear"))

    return decoys[:n]


def generate_stereo_swap_decoys(
    tokens: List[str],
    pool: DecoyPool,
    n: int = 5,
) -> List[Tuple[List[str], str]]:
    """Flip trans↔cis, remove, or add stereo constraints."""
    blocks = _detect_blocks(tokens)
    constraint_ranges = blocks.get("constraint_ranges", [])

    decoys = []

    if constraint_ranges:
        # Decoy 1: remove all stereo constraints
        d = list(tokens)
        # Remove from end to start to preserve indices
        for sr_start, sr_end in reversed(constraint_ranges):
            d[sr_start:sr_end] = []
        if _valid(d):
            decoys.append((d, "stereo_swap:remove_all"))

        # Decoy 2: flip trans↔cis
        d = list(tokens)
        for sr_start, sr_end in constraint_ranges:
            if d[sr_start] == "{trans:":
                d[sr_start] = "{cis:"
            elif d[sr_start] == "{cis:":
                d[sr_start] = "{trans:"
        if _valid(d):
            decoys.append((d, "stereo_swap:flip"))

        # Decoy 3: remove first constraint only
        if len(constraint_ranges) > 1:
            d = list(tokens)
            sr_start, sr_end = constraint_ranges[0]
            d[sr_start:sr_end] = []
            if _valid(d):
                decoys.append((d, "stereo_swap:remove_first"))

    else:
        # No stereo → add random stereo from pool
        if pool.stereo_blocks:
            for _ in range(min(n, len(pool.stereo_blocks))):
                stereo = pool.rng.choice(pool.stereo_blocks)
                d = list(tokens) + stereo
                if _valid(d):
                    decoys.append((d, "stereo_swap:add_random"))
                if len(decoys) >= n:
                    break

    return decoys[:n]


def generate_ligand_swap_decoys(
    tokens: List[str],
    pool: DecoyPool,
    n: int = 10,
) -> List[Tuple[List[str], str]]:
    """Replace one ligand's SMILES with a same-length ligand to keep sequence length constant."""
    blocks = _detect_blocks(tokens)
    lig_blocks = blocks["ligand_blocks"]
    if not lig_blocks:
        return []

    # Build a length-indexed pool for same-length matching
    len_pool = defaultdict(list)
    for lb_smi in pool.ligand_blocks:
        len_pool[len(lb_smi)].append(lb_smi)

    # Filter to blocks with valid SMILES ranges
    lig_blocks = [lb for lb in lig_blocks
                  if lb.get("smiles_start") is not None and lb.get("smiles_end") is not None]
    if not lig_blocks:
        return []

    decoys = []
    for _ in range(n * 5):
        lb = pool.rng.choice(lig_blocks)
        orig_smi = tokens[lb["smiles_start"]:lb["smiles_end"]]
        orig_len = len(orig_smi)
        orig_str = "".join(orig_smi)

        # Find a ligand with same token length (no padding needed)
        candidates = [s for s in len_pool.get(orig_len, [])
                      if "".join(s) != orig_str]
        if not candidates:
            # Allow ±2 tokens
            for delta in [1, -1, 2, -2]:
                candidates = [s for s in len_pool.get(orig_len + delta, [])
                              if "".join(s) != orig_str]
                if candidates:
                    break
        if not candidates:
            continue

        new_smi = list(pool.rng.choice(candidates))

        d = list(tokens)
        d[lb["smiles_start"]:lb["smiles_end"]] = new_smi
        if _valid(d):
            lid = tokens[lb["lid_pos"]]
            if lid == "L":
                lid = f"L{tokens[lb['lid_pos']+1]}"
            decoys.append((d, f"ligand_swap:{lid}"))
        if len(decoys) >= n:
            break

    return decoys[:n]


def generate_donor_pattern_swap_decoys(
    tokens: List[str],
    pool: DecoyPool,
    n: int = 5,
) -> List[Tuple[List[str], str]]:
    """Permute donor markers within a ligand (e.g., swap :N:1 and :O:2)."""
    blocks = _detect_blocks(tokens)
    lig_blocks = blocks["ligand_blocks"]
    if not lig_blocks:
        return []

    decoys = []
    for lb in lig_blocks:
        if lb.get("smiles_start") is None or lb.get("smiles_end") is None:
            continue
        # Find donor markers in this block
        donor_positions = []
        for i in range(lb["smiles_start"], lb["smiles_end"]):
            if i < len(tokens) and _DONOR_RE.match(tokens[i]):
                donor_positions.append(i)

        if len(donor_positions) < 2:
            continue

        # Generate permutations
        donor_tokens = [tokens[p] for p in donor_positions]
        for _ in range(n * 2):
            shuffled = list(donor_tokens)
            pool.rng.shuffle(shuffled)
            if shuffled == donor_tokens:
                continue
            d = list(tokens)
            for pos, new_tok in zip(donor_positions, shuffled):
                d[pos] = new_tok
            if _valid(d):
                lid = tokens[lb["lid_pos"]]
                if lid == "L":
                    lid = f"L{tokens[lb['lid_pos']+1]}"
                decoys.append((d, f"donor_pattern_swap:{lid}"))
            if len(decoys) >= n:
                break
        if len(decoys) >= n:
            break

    return decoys[:n]


def generate_boundary_decoys(
    tokens: List[str],
    pool: DecoyPool,
    n: int = 5,
) -> List[Tuple[List[str], str]]:
    """
    Near-valid perturbations:
    - CN ± 1
    - Single SMILES atom substitution (C↔N, O↔S, etc.)
    """
    blocks = _detect_blocks(tokens)
    decoys = []

    # CN ± 1
    if blocks["metal_range"] is not None:
        mp = blocks["metal_range"][0]
        parsed = _parse_metal_token(tokens[mp])
        if parsed:
            orig_metal, orig_cn = parsed
            for delta in [-1, +1, -2, +2]:
                new_cn = orig_cn + delta
                if new_cn < 1 or new_cn > 16:
                    continue
                new_tok = _rebuild_metal_token(tokens[mp], new_cn=new_cn)
                d = list(tokens)
                d[mp] = new_tok
                if _valid(d):
                    decoys.append((d, f"boundary:CN_{orig_cn}->{new_cn}"))
                if len(decoys) >= n:
                    break

    # Single-atom substitution in SMILES
    ATOM_SWAPS = {"C": "N", "N": "C", "O": "S", "S": "O", "P": "N"}
    lig_blocks = blocks["ligand_blocks"]
    for lb in lig_blocks:
        if lb.get("smiles_start") is None or lb.get("smiles_end") is None:
            continue
        for i in range(lb["smiles_start"], min(lb["smiles_end"], len(tokens))):
            tok = tokens[i]
            if tok in ATOM_SWAPS:
                d = list(tokens)
                d[i] = ATOM_SWAPS[tok]
                if _valid(d):
                    decoys.append((d, f"boundary:atom_{tok}->{ATOM_SWAPS[tok]}_@{i}"))
                if len(decoys) >= n:
                    return decoys[:n]
        if len(decoys) >= n:
            break

    return decoys[:n]


# ── unified generator ────────────────────────────────────────────

DECOY_TYPES = [
    "metal_swap", "shape_swap", "stereo_swap",
    "ligand_swap", "donor_pattern_swap", "boundary",
]

_GENERATORS = {
    "metal_swap": generate_metal_swap_decoys,
    "shape_swap": generate_shape_swap_decoys,
    "stereo_swap": generate_stereo_swap_decoys,
    "ligand_swap": generate_ligand_swap_decoys,
    "donor_pattern_swap": generate_donor_pattern_swap_decoys,
    "boundary": generate_boundary_decoys,
}


def generate_all_decoys(
    tokens: List[str],
    pool: DecoyPool,
    n_per_type: int = 10,
    types: Optional[List[str]] = None,
) -> List[Tuple[List[str], str]]:
    """Generate decoys of all (or specified) types."""
    if types is None:
        types = DECOY_TYPES
    all_decoys = []
    for dtype in types:
        gen = _GENERATORS.get(dtype)
        if gen:
            decoys = gen(tokens, pool, n=n_per_type)
            all_decoys.extend(decoys)
    return all_decoys
