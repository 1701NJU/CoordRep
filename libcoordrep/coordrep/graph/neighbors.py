"""
Coordination Neighbors Determination

1. BO
2.
3. (H donor)
4.
"""

import numpy as np
from typing import List, Tuple, Dict, Set
from dataclasses import dataclass

from ..io.tmqm_reader import RawMolecule, Atom, TRANSITION_METALS


DONOR_ELEMENTS = {
 'N', 'O', 'S', 'P', 'Cl', 'Br', 'I', 'F', # donor
 'C', # η-
 'Se', 'Te', 'As', #
}

METAL_DONOR_CUTOFFS = {
    # (metal, donor) -> max_distance
    'default': 2.8,
}

ATOMIC_NUMBERS = {
    'H': 1, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'P': 15, 'S': 16, 'Cl': 17,
    'Br': 35, 'I': 53, 'Se': 34, 'Te': 52, 'As': 33,
}


@dataclass
class NeighborInfo:
    """Neighbor information."""
    atom_idx: int
    element: str
    distance: float
    bond_order: float
    score: float
    
    def get_sort_key(self) -> tuple:
        """
 key
        
 1. score desc ()
        2. distance asc
        3. atom_Z desc
        4. original_index asc (tie-break)
        """
        atom_z = ATOMIC_NUMBERS.get(self.element, 0)
        return (
            -self.score,      # desc
            self.distance,    # asc
            -atom_z,          # desc
            self.atom_idx,    # asc (tie-break)
        )


def get_coordination_neighbors(mol: RawMolecule,
                               metal_idx: int,
                               bo_threshold: float = 0.3,
                               distance_cutoff: float = 2.8) -> List[NeighborInfo]:
    """
    
 1. BO
 2.
 3. donor H donor
 4.
    
    Args:
 mol:
 metal_idx:
 bo_threshold: BO
 distance_cutoff:
    
    Returns:
    """
    metal_coord = mol.atoms[metal_idx].coords
    
    candidates = {}  # atom_idx -> NeighborInfo
    
    if mol.bond_orders is not None:
        for i, bo in enumerate(mol.bond_orders[metal_idx]):
            if i == metal_idx:
                continue
            
            atom = mol.atoms[i]
            
            if atom.element not in DONOR_ELEMENTS:
                continue
            
            if bo > bo_threshold:
                dist = np.linalg.norm(atom.coords - metal_coord)
                candidates[i] = NeighborInfo(
                    atom_idx=i,
                    element=atom.element,
                    distance=dist,
                    bond_order=bo,
                    score=bo
                )
    
    for i, atom in enumerate(mol.atoms):
        if i == metal_idx:
            continue
        
        if atom.element not in DONOR_ELEMENTS:
            continue
        
        dist = np.linalg.norm(atom.coords - metal_coord)
        
        if dist < distance_cutoff and i not in candidates:
            bo = 0.0
            if mol.bond_orders is not None:
                bo = mol.bond_orders[metal_idx, i]
            
            score = 1.0 / max(dist, 0.5)
            
            candidates[i] = NeighborInfo(
                atom_idx=i,
                element=atom.element,
                distance=dist,
                bond_order=bo,
                score=score
            )
    
    for info in candidates.values():
        bo_score = info.bond_order * 2.0
        dist_score = 1.0 / max(info.distance, 0.5)
        info.score = bo_score + dist_score
    
    neighbors = sorted(candidates.values(), key=lambda x: x.get_sort_key())
    
    return neighbors


def stabilize_neighbors(neighbors: List[NeighborInfo],
                       expected_cn: int = None) -> List[NeighborInfo]:
    """
    
 expected_cn top-k
    """
    if expected_cn is None:
        return neighbors
    
    return neighbors[:expected_cn]


def get_donor_indices(neighbors: List[NeighborInfo]) -> List[int]:
    """Get donor atom indices."""
    return [n.atom_idx for n in neighbors]


def validate_donors(mol: RawMolecule, donor_indices: List[int]) -> List[str]:
    """Validate donors."""
    issues = []
    
    for idx in donor_indices:
        elem = mol.atoms[idx].element
        if elem == 'H':
            issues.append(f"H atom at index {idx} should not be a donor!")
        if elem not in DONOR_ELEMENTS:
            issues.append(f"Unusual donor element: {elem} at index {idx}")
    
    return issues





