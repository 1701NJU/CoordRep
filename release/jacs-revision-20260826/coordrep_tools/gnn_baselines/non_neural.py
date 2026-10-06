"""
non_neural.py
=============
Non-neural baselines for donor annotation and identity matching.

Donor annotation baselines:
- Random: random donor assignment
- GlobalFreq: global donor element frequency
- CondFreq(metal, CN): conditional frequency given metal + CN
- LigandFreq(SMILES): majority donor set from training SMILES

Identity baselines:
- WL hash: Weisfeiler-Leman graph hash
- Canonical SMILES multiset key
- ECFP Tanimoto similarity
"""

from __future__ import annotations

import hashlib
import re
from collections import Counter, defaultdict
from typing import Dict, List, Optional, Set, Tuple

import numpy as np

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit import DataStructs
except ImportError:
    Chem = None


# ── Donor annotation baselines ────────────────────────────

ALL_DONOR_ELEMENTS = ['C', 'N', 'O', 'S', 'P', 'F', 'Cl', 'Br', 'I', 'Se']


class RandomDonorBaseline:
    """Random k-of-n donor assignment."""

    def __init__(self, seed: int = 42):
        self.rng = np.random.RandomState(seed)

    def predict(self, mol_atoms: List[str], denticity: int) -> List[int]:
        """Return indices predicted as donor atoms."""
        heavy = [i for i, a in enumerate(mol_atoms) if a != 'H']
        if not heavy:
            return []
        k = min(denticity, len(heavy))
        return sorted(self.rng.choice(heavy, size=k, replace=False).tolist())


class GlobalFreqDonorBaseline:
    """Predict most-frequent donor elements globally."""

    def __init__(self):
        self.element_freq = Counter()

    def fit(self, train_samples: List[dict]):
        for s in train_samples:
            for d in s.get('donor_elements', []):
                self.element_freq[d] += 1
        return self

    def predict_elements(self, denticity: int) -> List[str]:
        """Return top-denticity donor elements."""
        ranked = [e for e, _ in self.element_freq.most_common()]
        return ranked[:denticity]


class CondFreqDonorBaseline:
    """Predict donor elements conditioned on (metal, CN)."""

    def __init__(self):
        self.freq = defaultdict(Counter)

    def fit(self, train_samples: List[dict]):
        for s in train_samples:
            key = (s.get('metal', '?'), s.get('cn', 0))
            for d in s.get('donor_elements', []):
                self.freq[key][d] += 1
        return self

    def predict_elements(self, metal: str, cn: int, denticity: int) -> List[str]:
        key = (metal, cn)
        if key in self.freq:
            ranked = [e for e, _ in self.freq[key].most_common()]
        else:
            # Fallback: marginal over metal
            merged = Counter()
            for (m, c), counts in self.freq.items():
                if m == metal:
                    merged += counts
            ranked = [e for e, _ in merged.most_common()]
        return ranked[:denticity]


class LigandFreqDonorBaseline:
    """Predict donor set from most-frequent training assignment for same SMILES."""

    def __init__(self):
        self.smiles_donors = defaultdict(Counter)
        self.fallback = GlobalFreqDonorBaseline()

    def fit(self, train_samples: List[dict]):
        for s in train_samples:
            smi = s.get('smiles', '')
            donor_key = tuple(sorted(s.get('donor_elements', [])))
            if smi and donor_key:
                self.smiles_donors[smi][donor_key] += 1
        self.fallback.fit(train_samples)
        return self

    def predict_elements(self, smiles: str, denticity: int) -> List[str]:
        if smiles in self.smiles_donors:
            best_set = self.smiles_donors[smiles].most_common(1)[0][0]
            return list(best_set)[:denticity]
        return self.fallback.predict_elements(denticity)

    def has_smiles(self, smiles: str) -> bool:
        return smiles in self.smiles_donors


# ── Identity baselines ────────────────────────────────────

def canonical_smiles_multiset_key(ligand_smiles: List[str]) -> str:
    """Canonical identity key from sorted ligand SMILES."""
    if Chem is None:
        canon = sorted(ligand_smiles)
    else:
        canon = []
        for s in ligand_smiles:
            mol = Chem.MolFromSmiles(s)
            if mol:
                canon.append(Chem.MolToSmiles(mol))
            else:
                canon.append(s)
        canon.sort()
    return '|'.join(canon)


def wl_hash(smiles: str, iterations: int = 3) -> str:
    """Weisfeiler-Leman hash of a molecular graph from SMILES."""
    if Chem is None:
        return hashlib.md5(smiles.encode()).hexdigest()
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is None:
            return hashlib.md5(smiles.encode()).hexdigest()

        # Initialize labels as element symbols
        labels = {i: atom.GetSymbol() for i, atom in enumerate(mol.GetAtoms())}

        for _ in range(iterations):
            new_labels = {}
            for i, atom in enumerate(mol.GetAtoms()):
                neighbors = sorted(labels[n.GetIdx()] for n in atom.GetNeighbors())
                new_labels[i] = hashlib.md5(
                    f"{labels[i]}|{'|'.join(neighbors)}".encode()
                ).hexdigest()[:8]
            labels = new_labels

        sorted_labels = sorted(labels.values())
        return hashlib.md5('|'.join(sorted_labels).encode()).hexdigest()
    except Exception:
        return hashlib.md5(smiles.encode()).hexdigest()


def ligand_set_wl_key(ligand_smiles: List[str], metal: str = '', cn: int = 0) -> str:
    """WL hash of the multiset of ligand graphs + metal + CN."""
    hashes = sorted(wl_hash(s) for s in ligand_smiles)
    combined = f"{metal}|{cn}|{'|'.join(hashes)}"
    return hashlib.md5(combined.encode()).hexdigest()


def ecfp_fingerprint(smiles: str, radius: int = 2, n_bits: int = 2048):
    """Morgan fingerprint (ECFP) from SMILES."""
    if Chem is None:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)


def tanimoto_similarity(fp1, fp2) -> float:
    """Tanimoto similarity between two fingerprints."""
    if fp1 is None or fp2 is None:
        return 0.0
    return DataStructs.TanimotoSimilarity(fp1, fp2)


def ligand_set_tanimoto(smiles_a: List[str], smiles_b: List[str]) -> float:
    """
    Max-weighted Tanimoto similarity between two ligand sets.
    Average of best-match Tanimoto for each ligand in set A.
    """
    fps_a = [ecfp_fingerprint(s) for s in smiles_a]
    fps_b = [ecfp_fingerprint(s) for s in smiles_b]

    if not fps_a or not fps_b:
        return 0.0

    scores = []
    for fa in fps_a:
        best = max(tanimoto_similarity(fa, fb) for fb in fps_b)
        scores.append(best)
    return float(np.mean(scores))
