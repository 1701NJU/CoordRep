"""
minphys_features.py
===================
Minimal-physics diagnostic features for CoordRep-Ranker ablation.

Three features per ligand:
  A. Heavy atom count bin (small / medium / large)
  B. Formal charge bin (anionic / neutral / cationic / mixed)
  C. Donor hard-soft class (hard / borderline / soft / mixed)

Plus a global summary across all ligands in the complex.
"""

from __future__ import annotations

import re
from collections import Counter
from typing import Dict, List, Optional, Tuple

import numpy as np

try:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Descriptors
    RDLogger.logger().setLevel(RDLogger.ERROR)
    _HAS_RDKIT = True
except ImportError:
    _HAS_RDKIT = False

from .ablation_masking import _detect_blocks

# ── donor hard-soft classification ────────────────────────────────

_DONOR_CLASS = {
    "O": "hard", "F": "hard",
    "N": "borderline", "Cl": "borderline", "Br": "borderline",
    "P": "soft", "S": "soft", "C": "soft", "I": "soft",
    "Se": "soft", "As": "soft", "Te": "soft",
}

_DONOR_RE = re.compile(r'^L\d+:([A-Z][a-z]?):\d+$')


def _classify_donor(atom: str) -> str:
    return _DONOR_CLASS.get(atom, "borderline")


# ── SMILES → heavy atom count ────────────────────────────────────

def _smiles_heavy_atoms(smiles: str) -> Optional[int]:
    """Return heavy atom count, or None if SMILES is unparseable."""
    if not _HAS_RDKIT or not smiles:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        # try sanitize=False
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is None:
            return None
        try:
            Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_FINDRADICALS |
                             Chem.SanitizeFlags.SANITIZE_SETAROMATICITY |
                             Chem.SanitizeFlags.SANITIZE_SETCONJUGATION |
                             Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION)
        except Exception:
            pass
    return mol.GetNumHeavyAtoms()


def _smiles_formal_charge(smiles: str) -> Optional[int]:
    """Return formal charge, or None if unparseable."""
    if not _HAS_RDKIT or not smiles:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is None:
        return None
    return Chem.GetFormalCharge(mol)


# ── extract ligand info from token list ───────────────────────────

def extract_ligand_features(tokens: List[str]) -> List[dict]:
    """
    Extract per-ligand features from a CoordRep token list.

    Returns list of dicts, one per ligand block:
      {lid, smiles, heavy_atoms, formal_charge, donor_atoms, donor_classes}
    """
    blocks = _detect_blocks(tokens)
    ligands = []

    # collect donor atoms per ligand
    donor_map: Dict[str, List[str]] = {}
    for tok in tokens:
        m = _DONOR_RE.match(tok)
        if m:
            lid = tok.split(":")[0]  # e.g. "L1"
            atom = m.group(1)
            donor_map.setdefault(lid, []).append(atom)

    for lb in blocks.get("ligand_blocks", []):
        lid_tok = tokens[lb["lid_pos"]]
        # normalise lid
        if lid_tok == "L" and lb["lid_pos"] + 1 < len(tokens):
            lid = "L" + tokens[lb["lid_pos"] + 1]
        else:
            lid = lid_tok

        smiles = ""
        if lb.get("smiles_start") is not None and lb.get("smiles_end") is not None:
            smiles = "".join(tokens[lb["smiles_start"]:lb["smiles_end"]])

        ha = _smiles_heavy_atoms(smiles)
        fc = _smiles_formal_charge(smiles)
        donors = donor_map.get(lid, [])

        ligands.append({
            "lid": lid,
            "smiles": smiles,
            "heavy_atoms": ha,
            "formal_charge": fc,
            "donor_atoms": donors,
            "donor_classes": list(set(_classify_donor(a) for a in donors)) if donors else [],
        })

    return ligands


# ── binning ───────────────────────────────────────────────────────

class MinPhysBinner:
    """Fit HA quantile bins on train data, then bin any complex."""

    def __init__(self):
        self.ha_q33 = 10
        self.ha_q66 = 20
        self._fitted = False

    def fit(self, train_token_lists: List[List[str]]):
        """Compute 33/66 quantiles from train ligand heavy atom counts."""
        all_ha = []
        for toks in train_token_lists:
            for lig in extract_ligand_features(toks):
                if lig["heavy_atoms"] is not None:
                    all_ha.append(lig["heavy_atoms"])
        if all_ha:
            self.ha_q33 = int(np.percentile(all_ha, 33))
            self.ha_q66 = int(np.percentile(all_ha, 66))
        self._fitted = True

    def bin_ha(self, ha: Optional[int]) -> str:
        if ha is None:
            return "unk"
        if ha <= self.ha_q33:
            return "small"
        if ha <= self.ha_q66:
            return "medium"
        return "large"

    @staticmethod
    def bin_charge(fc: Optional[int]) -> str:
        if fc is None:
            return "unk"
        if fc < 0:
            return "anionic"
        if fc > 0:
            return "cationic"
        return "neutral"

    @staticmethod
    def bin_donor_hardsoft(donor_classes: List[str]) -> str:
        if not donor_classes:
            return "unk"
        s = set(donor_classes)
        if len(s) == 1:
            return s.pop()
        return "mixed"

    def featurise_complex(self, tokens: List[str]) -> dict:
        """
        Return per-ligand and global MinPhys features.

        Returns:
          {
            "per_ligand": [{lid, ha_bin, charge_bin, donor_hs}, ...],
            "global": {size_load, charge_sum, hardsoft}
          }
        """
        ligs = extract_ligand_features(tokens)
        per_lig = []
        all_ha_bins = []
        all_charges = []
        all_donor_classes = set()

        for lig in ligs:
            ha_bin = self.bin_ha(lig["heavy_atoms"])
            charge_bin = self.bin_charge(lig["formal_charge"])
            donor_hs = self.bin_donor_hardsoft(lig["donor_classes"])

            per_lig.append({
                "lid": lig["lid"],
                "ha_bin": ha_bin,
                "charge_bin": charge_bin,
                "donor_hs": donor_hs,
            })
            all_ha_bins.append(ha_bin)
            if lig["formal_charge"] is not None:
                all_charges.append(lig["formal_charge"])
            all_donor_classes.update(lig["donor_classes"])

        # global bins
        ha_counts = Counter(all_ha_bins)
        most_common_ha = ha_counts.most_common(1)[0][0] if ha_counts else "unk"

        total_charge = sum(all_charges) if all_charges else None
        global_charge_bin = self.bin_charge(total_charge)

        global_hs = self.bin_donor_hardsoft(list(all_donor_classes))

        return {
            "per_ligand": per_lig,
            "global": {
                "size_load": most_common_ha,
                "charge_sum": global_charge_bin,
                "hardsoft": global_hs,
            },
        }
