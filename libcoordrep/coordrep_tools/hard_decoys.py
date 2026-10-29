"""
hard_decoys.py
==============
Hard-negative counterfactual decoy generators.

Every decoy MUST pass the full CoordRep validator.  The decoys are
designed to be **chemically plausible** so that a simple rule-based
validator cannot distinguish them from the real complex.

Hard decoy types
----------------
1. **ligand_replacement_hard** – swap a ligand SMILES with another
   real ligand of the same token-length (±5%), keeping all other
   fields unchanged.
2. **metal_replacement_hard** – swap metal to another metal that
   appears with the same CN in the training corpus.
3. **stereo_rearrangement_hard** – flip trans↔cis / fac↔mer only
   for complexes that have stereo constraints.
4. **co_ligand_set_hard** – replace the full ligand set with
   another real complex's ligand set that has the same CN and
   number of ligands.
5. **boundary_hard** – swap metal to an adjacent-group or same-
   group metal that shares the same CN in training data.
"""

from __future__ import annotations

import math
import random
import re
from collections import Counter, defaultdict
from copy import deepcopy
from typing import Dict, List, Optional, Set, Tuple

from .ablation_masking import _detect_blocks
from .validate import is_valid_coordrep

# ── helpers ──────────────────────────────────────────────────────

_METAL_RE = re.compile(r'^\[Metal:(\w+)\|.*CN:(\d+)\]$')
_DONOR_IN_CONSTRAINT = re.compile(r'^(L\d+):([A-Z][a-z]?):(\d+)$')


def _parse_metal_token(tok: str) -> Optional[Tuple[str, int]]:
    m = _METAL_RE.match(tok)
    if m:
        return m.group(1), int(m.group(2))
    return None


def _rebuild_metal_token(tok: str, new_metal: str = None,
                         new_cn: int = None) -> str:
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
    return is_valid_coordrep(s, strict=True)


def _extract_ligand_info(toks, blocks):
    """Return list of dicts with SMILES str, token-length, etc."""
    infos = []
    for lb in blocks["ligand_blocks"]:
        s, e = lb.get("smiles_start"), lb.get("smiles_end")
        if s is None or e is None:
            continue
        smi_toks = toks[s:e]
        smi_str = "".join(smi_toks)
        infos.append({
            "smi_str": smi_str,
            "smi_toks": list(smi_toks),
            "n_toks": len(smi_toks),
            "start": s,
            "end": e,
            "block": lb,
        })
    return infos


# ── Pool of hard-swap materials ──────────────────────────────────

class HardDecoyPool:
    """Pre-compute swap materials indexed for hard-negative generation."""

    def __init__(self, all_token_lists: List[List[str]], seed: int = 42):
        self.rng = random.Random(seed)

        # Metal → set of CNs seen in training
        self.metal_cns: Dict[str, Set[int]] = defaultdict(set)
        # CN → set of metals
        self.cn_metals: Dict[int, Set[str]] = defaultdict(set)
        # (metal, CN) count
        self.metal_cn_count: Counter = Counter()

        # Ligand SMILES token-lists indexed by token-length
        self.lig_by_len: Dict[int, List[List[str]]] = defaultdict(list)
        # Unique SMILES strings for dedup
        self.lig_str_set: Set[str] = set()

        # Full ligand-set signatures: (metal, CN, n_ligs) → list of full lig_block_toks
        self.full_lig_sets: Dict[Tuple, List[List[List[str]]]] = defaultdict(list)

        # Metal row mapping
        self._3d = {"Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn"}
        self._4d = {"Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd"}
        self._5d = {"La","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg"}

        # Same-group metals (periodic table neighbors)
        self.same_group = defaultdict(set)
        groups = [
            ["Sc","Y","La"], ["Ti","Zr","Hf"], ["V","Nb","Ta"],
            ["Cr","Mo","W"], ["Mn","Tc","Re"], ["Fe","Ru","Os"],
            ["Co","Rh","Ir"], ["Ni","Pd","Pt"], ["Cu","Ag","Au"],
            ["Zn","Cd","Hg"],
        ]
        for g in groups:
            for m in g:
                self.same_group[m] = set(g) - {m}

        for toks in all_token_lists:
            blocks = _detect_blocks(toks)

            # Metal/CN
            if blocks["metal_range"] is not None:
                mt = toks[blocks["metal_range"][0]]
                parsed = _parse_metal_token(mt)
                if parsed:
                    metal, cn = parsed
                    self.metal_cns[metal].add(cn)
                    self.cn_metals[cn].add(metal)
                    self.metal_cn_count[(metal, cn)] += 1

            # Ligand SMILES by length
            lig_infos = _extract_ligand_info(toks, blocks)
            lig_set = []
            for info in lig_infos:
                smi_str = info["smi_str"]
                smi_toks = info["smi_toks"]
                if smi_str not in self.lig_str_set:
                    self.lig_str_set.add(smi_str)
                    self.lig_by_len[len(smi_toks)].append(smi_toks)
                lig_set.append(smi_toks)

            # Full ligand sets
            if blocks["metal_range"] is not None and lig_infos:
                parsed = _parse_metal_token(toks[blocks["metal_range"][0]])
                if parsed:
                    key = (parsed[0], parsed[1], len(lig_infos))
                    self.full_lig_sets[key].append(lig_set)

    def metals_for_cn(self, cn: int, exclude: str = None) -> List[str]:
        choices = [m for m in self.cn_metals.get(cn, set()) if m != exclude]
        return choices

    def same_group_metals_for_cn(self, metal: str, cn: int) -> List[str]:
        """Metals from the same periodic-table group with same CN."""
        group = self.same_group.get(metal, set())
        return [m for m in group if cn in self.metal_cns.get(m, set())]

    def adjacent_metals_for_cn(self, metal: str, cn: int) -> List[str]:
        """Same-group + neighboring metals with same CN."""
        candidates = self.same_group_metals_for_cn(metal, cn)
        # Also add metals from same row with same CN
        row = None
        for r_set, r_name in [(self._3d, "3d"), (self._4d, "4d"), (self._5d, "5d")]:
            if metal in r_set:
                row = r_set
                break
        if row:
            for m in row:
                if m != metal and cn in self.metal_cns.get(m, set()):
                    if m not in candidates:
                        candidates.append(m)
        return candidates

    def same_length_ligands(self, orig_smi_toks: List[str],
                            tolerance: float = 0.05) -> List[List[str]]:
        """Find ligands with similar token count (within tolerance)."""
        orig_len = len(orig_smi_toks)
        orig_str = "".join(orig_smi_toks)
        lo = max(1, int(orig_len * (1 - tolerance)))
        hi = int(orig_len * (1 + tolerance)) + 1
        candidates = []
        for length in range(lo, hi + 1):
            for smi in self.lig_by_len.get(length, []):
                if "".join(smi) != orig_str:
                    candidates.append(smi)
        return candidates

    def co_ligand_sets(self, metal: str, cn: int, n_ligs: int,
                       exclude_set: Optional[str] = None) -> List[List[List[str]]]:
        """Full ligand sets from other complexes with same (metal, CN, n_ligs) or (*, CN, n_ligs)."""
        key_exact = (metal, cn, n_ligs)
        candidates = []
        for ls in self.full_lig_sets.get(key_exact, []):
            sig = "|".join("".join(l) for l in ls)
            if sig != exclude_set:
                candidates.append(ls)
        # Also allow different metals with same CN + n_ligs
        if len(candidates) < 5:
            for m2 in self.cn_metals.get(cn, set()):
                if m2 == metal:
                    continue
                for ls in self.full_lig_sets.get((m2, cn, n_ligs), []):
                    sig = "|".join("".join(l) for l in ls)
                    if sig != exclude_set:
                        candidates.append(ls)
                    if len(candidates) >= 20:
                        break
        return candidates


# ── Hard decoy generators ────────────────────────────────────────

def generate_metal_replacement_hard(
    tokens: List[str],
    pool: HardDecoyPool,
    n: int = 8,
) -> List[Tuple[List[str], str]]:
    """Replace metal with another metal that appears with same CN in corpus."""
    blocks = _detect_blocks(tokens)
    if blocks["metal_range"] is None:
        return []

    mp = blocks["metal_range"][0]
    parsed = _parse_metal_token(tokens[mp])
    if not parsed:
        return []
    orig_metal, orig_cn = parsed

    # Priority: same-group metals > same-row > any same-CN
    candidates = pool.same_group_metals_for_cn(orig_metal, orig_cn)
    if len(candidates) < n:
        candidates = pool.adjacent_metals_for_cn(orig_metal, orig_cn)
    if len(candidates) < n:
        extra = pool.metals_for_cn(orig_cn, exclude=orig_metal)
        for m in extra:
            if m not in candidates:
                candidates.append(m)

    decoys = []
    pool.rng.shuffle(candidates)
    for new_metal in candidates[:n * 2]:
        if new_metal == orig_metal:
            continue
        new_tok = _rebuild_metal_token(tokens[mp], new_metal=new_metal)
        d = list(tokens)
        d[mp] = new_tok
        if _valid(d):
            decoys.append((d, f"metal_hard:{orig_metal}->{new_metal}"))
        if len(decoys) >= n:
            break
    return decoys[:n]


def generate_ligand_replacement_hard(
    tokens: List[str],
    pool: HardDecoyPool,
    n: int = 8,
) -> List[Tuple[List[str], str]]:
    """Replace one ligand SMILES with a same-length ligand from the corpus."""
    blocks = _detect_blocks(tokens)
    lig_infos = _extract_ligand_info(tokens, blocks)
    if not lig_infos:
        return []

    decoys = []
    for _ in range(n * 5):
        info = pool.rng.choice(lig_infos)
        candidates = pool.same_length_ligands(info["smi_toks"], tolerance=0.0)
        if not candidates:
            candidates = pool.same_length_ligands(info["smi_toks"], tolerance=0.05)
        if not candidates:
            continue

        new_smi = pool.rng.choice(candidates)
        d = list(tokens)
        d[info["start"]:info["end"]] = new_smi
        if _valid(d):
            # Get lid label
            lb = info["block"]
            lid = tokens[lb["lid_pos"]]
            if lid == "L":
                lid = f"L{tokens[lb['lid_pos']+1]}"
            decoys.append((d, f"ligand_hard:{lid}"))
        if len(decoys) >= n:
            break
    return decoys[:n]


def generate_stereo_rearrangement_hard(
    tokens: List[str],
    pool: HardDecoyPool,
    n: int = 5,
) -> List[Tuple[List[str], str]]:
    """Flip trans↔cis / fac↔mer in existing constraints."""
    blocks = _detect_blocks(tokens)
    constraint_ranges = blocks.get("constraint_ranges", [])
    if not constraint_ranges:
        return []

    decoys = []

    # Decoy 1: flip all trans↔cis
    d = list(tokens)
    flipped = False
    for cr_start, cr_end in constraint_ranges:
        if d[cr_start] == "{trans:":
            d[cr_start] = "{cis:"
            flipped = True
        elif d[cr_start] == "{cis:":
            d[cr_start] = "{trans:"
            flipped = True
        elif d[cr_start] == "{fac:":
            d[cr_start] = "{mer:"
            flipped = True
        elif d[cr_start] == "{mer:":
            d[cr_start] = "{fac:"
            flipped = True
    if flipped and _valid(d):
        decoys.append((d, "stereo_hard:flip_all"))

    # Decoy 2: flip only first constraint
    if len(constraint_ranges) >= 2:
        d = list(tokens)
        cr_start = constraint_ranges[0][0]
        if d[cr_start] == "{trans:":
            d[cr_start] = "{cis:"
        elif d[cr_start] == "{cis:":
            d[cr_start] = "{trans:"
        if _valid(d):
            decoys.append((d, "stereo_hard:flip_first"))

    # Decoy 3: swap donor references between two constraints
    if len(constraint_ranges) >= 2:
        d = list(tokens)
        cr1 = constraint_ranges[0]
        cr2 = constraint_ranges[1]
        # Swap the content (excluding the type token and closing bracket)
        body1 = d[cr1[0]+1:cr1[1]-1]
        body2 = d[cr2[0]+1:cr2[1]-1]
        if len(body1) == len(body2):  # same length ensures valid splice
            d2 = list(tokens)
            d2[cr1[0]+1:cr1[1]-1] = body2
            d2[cr2[0]+1:cr2[1]-1] = body1
            if _valid(d2):
                decoys.append((d2, "stereo_hard:swap_refs"))

    return decoys[:n]


def generate_co_ligand_set_hard(
    tokens: List[str],
    pool: HardDecoyPool,
    n: int = 5,
) -> List[Tuple[List[str], str]]:
    """Replace full ligand set with another complex's ligand set (same CN, n_ligs)."""
    blocks = _detect_blocks(tokens)
    if blocks["metal_range"] is None:
        return []
    parsed = _parse_metal_token(tokens[blocks["metal_range"][0]])
    if not parsed:
        return []
    metal, cn = parsed
    lig_infos = _extract_ligand_info(tokens, blocks)
    if not lig_infos:
        return []

    n_ligs = len(lig_infos)
    exclude_sig = "|".join(info["smi_str"] for info in lig_infos)
    candidates = pool.co_ligand_sets(metal, cn, n_ligs, exclude_set=exclude_sig)
    if not candidates:
        return []

    decoys = []
    pool.rng.shuffle(candidates)
    for new_set in candidates[:n * 3]:
        if len(new_set) != n_ligs:
            continue
        d = list(tokens)
        # Replace ligand SMILES from end to start to preserve indices
        valid_replace = True
        for i in range(n_ligs - 1, -1, -1):
            info = lig_infos[i]
            new_smi = new_set[i]
            d[info["start"]:info["end"]] = new_smi
        if _valid(d):
            decoys.append((d, f"co_ligand_hard:n={n_ligs}"))
        if len(decoys) >= n:
            break
    return decoys[:n]


def generate_boundary_hard(
    tokens: List[str],
    pool: HardDecoyPool,
    n: int = 5,
) -> List[Tuple[List[str], str]]:
    """Same-group metal swap (hardest: very similar chemistry)."""
    blocks = _detect_blocks(tokens)
    if blocks["metal_range"] is None:
        return []
    mp = blocks["metal_range"][0]
    parsed = _parse_metal_token(tokens[mp])
    if not parsed:
        return []
    orig_metal, orig_cn = parsed

    # Only same-group metals with same CN
    candidates = pool.same_group_metals_for_cn(orig_metal, orig_cn)
    decoys = []
    for new_metal in candidates:
        if new_metal == orig_metal:
            continue
        new_tok = _rebuild_metal_token(tokens[mp], new_metal=new_metal)
        d = list(tokens)
        d[mp] = new_tok
        if _valid(d):
            decoys.append((d, f"boundary_hard:{orig_metal}->{new_metal}"))
        if len(decoys) >= n:
            break
    return decoys[:n]


# ── unified generator ────────────────────────────────────────────

HARD_DECOY_TYPES = [
    "metal_hard", "ligand_hard", "stereo_hard",
    "co_ligand_hard", "boundary_hard",
]

_HARD_GENERATORS = {
    "metal_hard": generate_metal_replacement_hard,
    "ligand_hard": generate_ligand_replacement_hard,
    "stereo_hard": generate_stereo_rearrangement_hard,
    "co_ligand_hard": generate_co_ligand_set_hard,
    "boundary_hard": generate_boundary_hard,
}


def generate_all_hard_decoys(
    tokens: List[str],
    pool: HardDecoyPool,
    n_per_type: int = 8,
    types: Optional[List[str]] = None,
) -> List[Tuple[List[str], str]]:
    """Generate hard decoys of all (or specified) types."""
    if types is None:
        types = HARD_DECOY_TYPES
    all_decoys = []
    for dtype in types:
        gen = _HARD_GENERATORS.get(dtype)
        if gen:
            decoys = gen(tokens, pool, n=n_per_type)
            all_decoys.extend(decoys)
    return all_decoys
