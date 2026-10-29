"""
stability_consensus.py
======================

For each structure, run N perturbation trials at small σ values and compute
consensus statistics on shape assignment.  Structures where consensus drops
below a configurable threshold are flagged as **geometry-boundary** rather
than being assigned a single deterministic shape label.

Output per molecule
-------------------
- best_shape_consensus   : fraction of trials agreeing on best-shape label
- top2_rank_consensus    : fraction of trials with identical (top1, top2) ordering
- shape_id_consensus     : fraction of trials producing the same L1 ShapeID token
- is_geometry_boundary   : True if any consensus < threshold
- dominant_best_shape    : most frequent best-shape label across trials
- dominant_shape_token   : most frequent L1 ShapeID token across trials
"""

from __future__ import annotations

import warnings
from collections import Counter
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np

from .shape_binning import bin_shape_from_string, DEFAULT_SHAPE_BIN_CONFIG

warnings.filterwarnings("ignore")


@dataclass
class ConsensusResult:
    """Consensus statistics for one molecule."""
    mol_id: str
    sigma: float
    n_trials: int
    n_success: int

    best_shape_consensus: float
    top2_rank_consensus: float
    shape_id_consensus: float

    dominant_best_shape: str
    dominant_shape_token: str
    is_geometry_boundary: bool

    best_shape_counts: Dict[str, int] = field(default_factory=dict)
    shape_token_counts: Dict[str, int] = field(default_factory=dict)

    def as_dict(self) -> dict:
        return {
            "mol_id": self.mol_id,
            "sigma": self.sigma,
            "n_trials": self.n_trials,
            "n_success": self.n_success,
            "best_shape_consensus": round(self.best_shape_consensus, 4),
            "top2_rank_consensus": round(self.top2_rank_consensus, 4),
            "shape_id_consensus": round(self.shape_id_consensus, 4),
            "dominant_best_shape": self.dominant_best_shape,
            "dominant_shape_token": self.dominant_shape_token,
            "is_geometry_boundary": self.is_geometry_boundary,
        }


# ──────────────────────────────────────────────────────────────────
# Noise utility (local copy to avoid circular import)
# ──────────────────────────────────────────────────────────────────

def _add_noise(atoms, sigma: float, seed: int):
    """Return a copy of *atoms* with Gaussian coordinate noise."""
    from coordrep.io.tmqm_reader import Atom
    rng = np.random.RandomState(seed)
    new_atoms = []
    for a in atoms:
        dx, dy, dz = rng.normal(0, sigma, 3)
        new_atoms.append(Atom(
            index=a.index, element=a.element,
            x=a.x + dx, y=a.y + dy, z=a.z + dz,
        ))
    return new_atoms


# ──────────────────────────────────────────────────────────────────
# Core consensus function
# ──────────────────────────────────────────────────────────────────

def compute_consensus(
    mol,
    sigma: float = 0.005,
    n_trials: int = 50,
    consensus_threshold: float = 0.90,
    shape_bin_config: dict | None = None,
    config=None,
) -> ConsensusResult:
    """
    Run *n_trials* Gaussian perturbations and measure shape-assignment
    consensus.

    Parameters
    ----------
    mol : RawMolecule
        Input molecule.
    sigma : float
        Standard deviation of coordinate noise (Å).
    n_trials : int
        Number of perturbation trials.
    consensus_threshold : float
        Below this, the structure is flagged as geometry-boundary.
    shape_bin_config : dict, optional
        Shape-binning thresholds.
    config : CoordRepConfig, optional
        Encoding config.  Defaults to ``CoordRepConfig.default()``.

    Returns
    -------
    ConsensusResult
    """
    from coordrep import encode_molecule, CoordRepConfig
    from coordrep.io.tmqm_reader import RawMolecule

    cfg = config or CoordRepConfig.default()
    sbcfg = shape_bin_config or DEFAULT_SHAPE_BIN_CONFIG

    best_shapes: List[str] = []
    top2_ranks: List[Tuple[str, str]] = []
    shape_tokens: List[str] = []

    for trial in range(n_trials):
        noisy_atoms = _add_noise(mol.atoms, sigma, seed=trial)
        noisy_mol = RawMolecule(
            mol_id=mol.mol_id,
            atoms=noisy_atoms,
            bond_orders=mol.bond_orders,
            properties=getattr(mol, "properties", None),
        )
        try:
            cc = encode_molecule(noisy_mol, cfg)
            s = cc.canonicalize().to_string()
        except Exception:
            continue

        if not s:
            continue

        binned = bin_shape_from_string(s, sbcfg)
        best_shapes.append(binned.best_shape)
        top2_ranks.append((binned.best_shape, binned.second_shape))
        shape_tokens.append(binned.token)

    n_success = len(best_shapes)
    if n_success == 0:
        return ConsensusResult(
            mol_id=mol.mol_id, sigma=sigma, n_trials=n_trials,
            n_success=0,
            best_shape_consensus=0.0,
            top2_rank_consensus=0.0,
            shape_id_consensus=0.0,
            dominant_best_shape="?",
            dominant_shape_token="?",
            is_geometry_boundary=True,
        )

    bs_counts = Counter(best_shapes)
    t2_counts = Counter(top2_ranks)
    st_counts = Counter(shape_tokens)

    bs_consensus = bs_counts.most_common(1)[0][1] / n_success
    t2_consensus = t2_counts.most_common(1)[0][1] / n_success
    st_consensus = st_counts.most_common(1)[0][1] / n_success

    is_boundary = min(bs_consensus, t2_consensus, st_consensus) < consensus_threshold

    return ConsensusResult(
        mol_id=mol.mol_id,
        sigma=sigma,
        n_trials=n_trials,
        n_success=n_success,
        best_shape_consensus=bs_consensus,
        top2_rank_consensus=t2_consensus,
        shape_id_consensus=st_consensus,
        dominant_best_shape=bs_counts.most_common(1)[0][0],
        dominant_shape_token=st_counts.most_common(1)[0][0],
        is_geometry_boundary=is_boundary,
        best_shape_counts=dict(bs_counts),
        shape_token_counts=dict(st_counts),
    )


# ──────────────────────────────────────────────────────────────────
# Batch API
# ──────────────────────────────────────────────────────────────────

def compute_consensus_batch(
    molecules,
    sigmas: List[float] = (0.005, 0.01),
    n_trials: int = 50,
    consensus_threshold: float = 0.90,
    shape_bin_config: dict | None = None,
    config=None,
    verbose: bool = True,
) -> List[ConsensusResult]:
    """
    Run consensus analysis on a list of molecules at multiple σ values.
    """
    results = []
    total = len(molecules) * len(sigmas)
    done = 0
    for mol in molecules:
        for sigma in sigmas:
            cr = compute_consensus(
                mol, sigma=sigma, n_trials=n_trials,
                consensus_threshold=consensus_threshold,
                shape_bin_config=shape_bin_config,
                config=config,
            )
            results.append(cr)
            done += 1
            if verbose and done % 20 == 0:
                print(f"  consensus: {done}/{total}")
    return results
