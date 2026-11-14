"""
CoordRep v2-beta core data structures.

Two principal abstractions:
  A. MultiMetalRecord  – metal-centered record graph for multinuclear systems
  B. HapticRecord      – coordination-site object for π/haptic systems
"""

from __future__ import annotations

import hashlib
import re
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple


# ════════════════════════════════════════════════════════════════════
# Shared primitives
# ════════════════════════════════════════════════════════════════════

@dataclass
class MetalCenter:
    label: str                          # M1, M2, …
    element: str
    oxidation: Optional[int] = None     # None → unk
    dcount: Optional[int] = None
    cn_site: int = 0                    # coordination number counted by sites
    cn_atom: int = 0                    # coordination number counted by atoms
    eta_sum: int = 0                    # total η across all haptic sites
    local_shape_best: str = ""
    local_shape_class: str = ""
    local_shape_delta: int = 0
    local_shape_vals: str = ""
    local_sites: List[str] = field(default_factory=list)   # site labels
    meta: Dict[str, Any] = field(default_factory=dict)

    @property
    def ox_str(self) -> str:
        if self.oxidation is None:
            return "unk"
        return f"+{self.oxidation}" if self.oxidation >= 0 else str(self.oxidation)

    def chemical_signature(self) -> tuple:
        return (
            self.element,
            self.ox_str,
            self.dcount or -1,
            self.cn_site,
            self.cn_atom,
            self.eta_sum,
            self.local_shape_best,
            tuple(sorted(self.local_sites)),
        )

    def to_dict(self) -> dict:
        return {
            "label": self.label,
            "element": self.element,
            "oxidation": self.oxidation,
            "dcount": self.dcount,
            "cn_site": self.cn_site,
            "cn_atom": self.cn_atom,
            "eta_sum": self.eta_sum,
            "local_shape_best": self.local_shape_best,
        }


@dataclass
class CoordinationSite:
    """A coordination site — atom, atom_set, pi_fragment, or centroid."""
    label: str                          # S1, S2, …
    site_type: str                      # atom | atom_set | pi_fragment | centroid
    ligand_label: str                   # L1, L2, …
    donor_atoms: List[str]              # atom labels within ligand
    donor_elements: List[str]           # element symbols
    eta: int = 1
    mu: int = 1
    target_metals: List[str] = field(default_factory=list)  # M1, M2, …
    mode: str = ""                      # e.g. Cp, allyl, alkene, mu2-O, …
    centroid_label: str = ""            # for pi_fragment / centroid types
    meta: Dict[str, Any] = field(default_factory=dict)

    def canonical_atom_key(self) -> tuple:
        return tuple(sorted(zip(self.donor_elements, self.donor_atoms)))

    def site_signature(self) -> tuple:
        return (
            self.site_type,
            self.canonical_atom_key(),
            self.eta,
            self.mu,
            tuple(sorted(self.target_metals)),
            self.mode,
        )

    def to_dict(self) -> dict:
        return {
            "label": self.label,
            "site_type": self.site_type,
            "ligand_label": self.ligand_label,
            "donor_atoms": self.donor_atoms,
            "donor_elements": self.donor_elements,
            "eta": self.eta,
            "mu": self.mu,
            "target_metals": self.target_metals,
            "mode": self.mode,
            "centroid_label": self.centroid_label,
        }


@dataclass
class Ligand:
    label: str                          # L1, L2, …
    smiles: str
    dent: int = 1
    charge: Optional[int] = None
    sites: List[str] = field(default_factory=list)  # site labels
    meta: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict:
        return {
            "label": self.label,
            "smiles": self.smiles,
            "dent": self.dent,
            "charge": self.charge,
            "sites": self.sites,
        }


@dataclass
class MetalEdge:
    """An edge in the metal graph."""
    m1: str
    m2: str
    relation: str       # bridged_only | direct_MM_bond | contact | unsupported
    d_mm: Optional[float] = None
    mm_bond: str = "unknown"    # yes | no | contact | unknown
    bridges: List[str] = field(default_factory=list)  # bridge site labels
    confidence: str = "medium"

    def canonical_key(self) -> tuple:
        a, b = sorted([self.m1, self.m2])
        return (a, b, self.relation)

    def to_dict(self) -> dict:
        return {
            "m1": self.m1, "m2": self.m2,
            "relation": self.relation,
            "d_mm": self.d_mm,
            "mm_bond": self.mm_bond,
            "bridges": self.bridges,
            "confidence": self.confidence,
        }


# ════════════════════════════════════════════════════════════════════
# A. MultiMetalRecord
# ════════════════════════════════════════════════════════════════════

@dataclass
class MultiIdentityKeys:
    L0_GlobalState: str = ""
    L1_GlobalShape: str = ""
    L2_MetalGraphTopo: str = ""
    L3_Connectivity: str = ""
    local_ids: Dict[str, str] = field(default_factory=dict)  # M_label → local ID

    def to_dict(self) -> dict:
        return {
            "L0_GlobalState": self.L0_GlobalState,
            "L1_GlobalShape": self.L1_GlobalShape,
            "L2_MetalGraphTopo": self.L2_MetalGraphTopo,
            "L3_Connectivity": self.L3_Connectivity,
            "local_ids": self.local_ids,
        }


@dataclass
class MultiMetalRecord:
    """CoordRep-Multi-v2beta record."""
    case_id: str
    refcode: str
    metals: List[MetalCenter]
    metal_edges: List[MetalEdge]
    sites: List[CoordinationSite]
    ligands: List[Ligand]
    identity: MultiIdentityKeys = field(default_factory=MultiIdentityKeys)
    source_database: str = "CSD_curated"
    description: str = ""
    notes: str = ""

    def to_dict(self) -> dict:
        return {
            "case_id": self.case_id,
            "refcode": self.refcode,
            "description": self.description,
            "metals": [m.to_dict() for m in self.metals],
            "metal_edges": [e.to_dict() for e in self.metal_edges],
            "sites": [s.to_dict() for s in self.sites],
            "ligands": [lg.to_dict() for lg in self.ligands],
            "identity": self.identity.to_dict(),
        }


# ════════════════════════════════════════════════════════════════════
# B. HapticRecord
# ════════════════════════════════════════════════════════════════════

@dataclass
class HapticIdentityKeys:
    L0_HapticState: str = ""
    L1_HapticShape: str = ""
    L2_SiteTopo: str = ""
    L3_Connectivity: str = ""

    def to_dict(self) -> dict:
        return {
            "L0_HapticState": self.L0_HapticState,
            "L1_HapticShape": self.L1_HapticShape,
            "L2_SiteTopo": self.L2_SiteTopo,
            "L3_Connectivity": self.L3_Connectivity,
        }


@dataclass
class HapticRecord:
    """CoordRep-Haptic-v2beta record."""
    case_id: str
    refcode: str
    metal: MetalCenter
    sites: List[CoordinationSite]
    ligands: List[Ligand]
    identity: HapticIdentityKeys = field(default_factory=HapticIdentityKeys)
    source_database: str = "CSD_curated"
    description: str = ""
    notes: str = ""

    def to_dict(self) -> dict:
        return {
            "case_id": self.case_id,
            "refcode": self.refcode,
            "description": self.description,
            "metal": self.metal.to_dict(),
            "sites": [s.to_dict() for s in self.sites],
            "ligands": [lg.to_dict() for lg in self.ligands],
            "identity": self.identity.to_dict(),
        }
