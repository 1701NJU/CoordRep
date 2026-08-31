"""
identity_keys.py
================

Parse a full CoordRep string (or a CoordComplex object) and return four
hierarchical identity keys:

    L0  StateKey  – the full canonical string itself (exact snapshot)
    L1  ShapeID   – metal + oxidation + binned CShM + constraint types + ligands
    L2  TopoID    – metal + CN + best-shape label + sorted ligand SMILES
    L3  ConnID    – metal + CN + sorted ligand SMILES (no geometry at all)

Each successive level is strictly coarser than the previous one:
    L0 == L0  ⟹  L1 == L1  ⟹  L2 == L2  ⟹  L3 == L3
"""

from __future__ import annotations

import hashlib
import re
from dataclasses import dataclass
from typing import List, Optional, Tuple

from .shape_binning import bin_shape_from_string, BinnedShape, DEFAULT_SHAPE_BIN_CONFIG


# ──────────────────────────────────────────────────────────────────
# Result container
# ──────────────────────────────────────────────────────────────────

@dataclass
class IdentityKeys:
    """Container holding all four identity layers for one complex."""
    L0_StateKey: str
    L1_ShapeID: str
    L2_TopoID: str
    L3_ConnID: str

    # Parsed sub-fields (for downstream analysis)
    metal: str = ""
    oxidation: str = ""
    cn: int = 0
    best_shape: str = ""
    binned_shape: Optional[BinnedShape] = None
    ligand_smiles: Tuple[str, ...] = ()
    constraint_signature: str = ""

    def as_dict(self) -> dict:
        return {
            "L0_StateKey": self.L0_StateKey,
            "L1_ShapeID": self.L1_ShapeID,
            "L2_TopoID": self.L2_TopoID,
            "L3_ConnID": self.L3_ConnID,
            "metal": self.metal,
            "oxidation": self.oxidation,
            "cn": self.cn,
            "best_shape": self.best_shape,
            "ligand_smiles": list(self.ligand_smiles),
            "constraint_signature": self.constraint_signature,
        }


# ──────────────────────────────────────────────────────────────────
# String-level parser helpers
# ──────────────────────────────────────────────────────────────────

_RE_METAL = re.compile(
    r'\[Metal:(?P<elem>\w+)'
    r'(?:\|ox:(?P<ox>[+\-]?\d+))?'
    r'(?:\|d:d(?P<dcount>\d+))?'
    r'\|CN:(?P<cn>\d+)\]'
)

_RE_SHAPE = re.compile(
    r'<ShapeBest:(?P<best>\w[\w-]*)'
    r'\|Class:(?P<cls>\w+)'
    r'\|Delta:(?P<dbin>\d+)'
    r'\|V:(?P<vals>[^>]+)>'
)

_RE_LIGAND = re.compile(r'\|(?P<lid>L\d+)=(?P<smi>[^|]+)')

_RE_TRANS = re.compile(r'\{trans:([^}]+)\}')
_RE_CIS = re.compile(r'\{cis:([^}]+)\}')
_RE_FM = re.compile(r'\{fm:(\w+)\}')


def _parse_metal(s: str) -> Tuple[str, str, int]:
    """Return (element, oxidation_str, CN)."""
    m = _RE_METAL.search(s)
    if not m:
        return ("?", "?", 0)
    elem = m.group("elem")
    ox = m.group("ox") or "?"
    cn = int(m.group("cn"))
    return elem, ox, cn


def _parse_ligands(s: str) -> Tuple[str, ...]:
    """Return sorted tuple of ligand SMILES."""
    ligs = _RE_LIGAND.findall(s)
    smiles_list = [smi for _, smi in ligs]
    return tuple(sorted(smiles_list))


def _parse_constraint_signature(s: str) -> str:
    """
    Build a constraint *type* signature that ignores exact site indices.

    Examples:
        "{trans:L1:N:1--L2:Cl:1}{fm:fac}" → "T1_fac"
        "{trans:L1:N:1--L1:N:2}{trans:L2:O:1--L3:O:1}" → "T2"
    """
    n_trans = len(_RE_TRANS.findall(s))
    n_cis = len(_RE_CIS.findall(s))
    fm = _RE_FM.search(s)
    parts = []
    if n_trans:
        parts.append(f"T{n_trans}")
    if n_cis:
        parts.append(f"C{n_cis}")
    if fm:
        parts.append(fm.group(1))
    return "_".join(parts) if parts else "none"


# ──────────────────────────────────────────────────────────────────
# Public API: from string
# ──────────────────────────────────────────────────────────────────

def extract_identity_keys(
    coordrep_str: str,
    shape_bin_config: dict | None = None,
) -> IdentityKeys:
    """
    Extract all four identity layers from a canonical CoordRep string.

    Parameters
    ----------
    coordrep_str : str
        Full canonical CoordRep string.
    shape_bin_config : dict, optional
        Override default shape-binning thresholds (see shape_binning.py).

    Returns
    -------
    IdentityKeys
    """
    cfg = shape_bin_config or DEFAULT_SHAPE_BIN_CONFIG

    # ── Parse sub-fields ──────────────────────────────────────────
    metal, ox, cn = _parse_metal(coordrep_str)
    ligands = _parse_ligands(coordrep_str)
    constraint_sig = _parse_constraint_signature(coordrep_str)
    binned = bin_shape_from_string(coordrep_str, cfg)

    # ── L0: full string ───────────────────────────────────────────
    l0 = coordrep_str

    # ── L1: ShapeID ───────────────────────────────────────────────
    #   metal + ox + binned_shape_token + constraint_type_sig + ligands
    lig_block = ";".join(ligands)
    l1 = f"{metal}|{ox}|{binned.token}|{constraint_sig}|{lig_block}"

    # ── L2: TopoID ────────────────────────────────────────────────
    #   metal + CN + best_shape_label + sorted ligands
    best = binned.best_shape if binned else "?"
    l2 = f"{metal}|CN{cn}|{best}|{lig_block}"

    # ── L3: ConnID ────────────────────────────────────────────────
    #   metal + CN + sorted ligands  (no geometry)
    l3 = f"{metal}|CN{cn}|{lig_block}"

    return IdentityKeys(
        L0_StateKey=l0,
        L1_ShapeID=l1,
        L2_TopoID=l2,
        L3_ConnID=l3,
        metal=metal,
        oxidation=ox,
        cn=cn,
        best_shape=best,
        binned_shape=binned,
        ligand_smiles=ligands,
        constraint_signature=constraint_sig,
    )


# ──────────────────────────────────────────────────────────────────
# Public API: from CoordComplex object
# ──────────────────────────────────────────────────────────────────

def extract_identity_keys_from_cc(
    cc,
    shape_bin_config: dict | None = None,
) -> IdentityKeys:
    """
    Extract identity keys directly from a ``CoordComplex`` object.

    The complex is canonicalized and serialized first.
    """
    cc_canon = cc.canonicalize()
    coordrep_str = cc_canon.to_string()
    return extract_identity_keys(coordrep_str, shape_bin_config)
