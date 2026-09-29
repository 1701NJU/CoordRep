"""
Graph Builder

 BO/
"""

import numpy as np
import networkx as nx
from typing import List, Tuple, Optional, Dict
from dataclasses import dataclass

from ..core import (
    MetalState, LigandModule, AssemblyGraph, 
    CoordRepIssue, IssueCodes, CoordRepConfig
)
from ..io.tmqm_reader import RawMolecule, Atom, TRANSITION_METALS


METAL_INFO = {
    # 3d (row 4)
    'Sc': {'row': 4, 'group': 3, 'en': 1.36},
    'Ti': {'row': 4, 'group': 4, 'en': 1.54},
    'V':  {'row': 4, 'group': 5, 'en': 1.63},
    'Cr': {'row': 4, 'group': 6, 'en': 1.66},
    'Mn': {'row': 4, 'group': 7, 'en': 1.55},
    'Fe': {'row': 4, 'group': 8, 'en': 1.83},
    'Co': {'row': 4, 'group': 9, 'en': 1.88},
    'Ni': {'row': 4, 'group': 10, 'en': 1.91},
    'Cu': {'row': 4, 'group': 11, 'en': 1.90},
    'Zn': {'row': 4, 'group': 12, 'en': 1.65},
    # 4d (row 5)
    'Y':  {'row': 5, 'group': 3, 'en': 1.22},
    'Zr': {'row': 5, 'group': 4, 'en': 1.33},
    'Nb': {'row': 5, 'group': 5, 'en': 1.6},
    'Mo': {'row': 5, 'group': 6, 'en': 2.16},
    'Tc': {'row': 5, 'group': 7, 'en': 1.9},
    'Ru': {'row': 5, 'group': 8, 'en': 2.2},
    'Rh': {'row': 5, 'group': 9, 'en': 2.28},
    'Pd': {'row': 5, 'group': 10, 'en': 2.20},
    'Ag': {'row': 5, 'group': 11, 'en': 1.93},
    'Cd': {'row': 5, 'group': 12, 'en': 1.69},
    # 5d (row 6)
    'La': {'row': 6, 'group': 3, 'en': 1.1},
    'Hf': {'row': 6, 'group': 4, 'en': 1.3},
    'Ta': {'row': 6, 'group': 5, 'en': 1.5},
    'W':  {'row': 6, 'group': 6, 'en': 2.36},
    'Re': {'row': 6, 'group': 7, 'en': 1.9},
    'Os': {'row': 6, 'group': 8, 'en': 2.2},
    'Ir': {'row': 6, 'group': 9, 'en': 2.2},
    'Pt': {'row': 6, 'group': 10, 'en': 2.28},
    'Au': {'row': 6, 'group': 11, 'en': 2.54},
    'Hg': {'row': 6, 'group': 12, 'en': 2.0},
}


@dataclass
class GraphBuildResult:
    """Graph build result."""
    metal_idx: int
    metal_state: MetalState
    donor_indices: List[int]
    donor_bond_orders: Dict[int, float]
    molecular_graph: nx.Graph
    issues: List[CoordRepIssue]


class GraphBuilder:
    """Graph builder."""
    
    def __init__(self, config: CoordRepConfig = None):
        self.config = config or CoordRepConfig.default()
    
    def _find_metal(self, mol: RawMolecule) -> Tuple[Optional[int], List[CoordRepIssue]]:
        """Find metal center."""
        issues = []
        metal_indices = [i for i, a in enumerate(mol.atoms) if a.element in TRANSITION_METALS]
        
        if len(metal_indices) == 0:
            issues.append(CoordRepIssue(
                code=IssueCodes.NO_METAL,
                message="No transition metal found",
                severity="error"
            ))
            return None, issues
        
        if len(metal_indices) > 1 and not self.config.allow_multimetal:
            issues.append(CoordRepIssue(
                code=IssueCodes.MULTI_METAL,
                message=f"Multiple metals found: {metal_indices}. Using first.",
                severity="warning"
            ))
        
        return metal_indices[0], issues
    
    def _get_bo_threshold(self, mol: RawMolecule, metal_idx: int) -> float:
        """Get BO threshold."""
        if self.config.bo_threshold_mode == "fixed":
            return self.config.bo_threshold_fixed
        
        if mol.bond_orders is None:
            return 0.3
        
        metal_bos = mol.bond_orders[metal_idx]
        nonzero_bos = metal_bos[metal_bos > 0.01]
        
        if len(nonzero_bos) == 0:
            return 0.3
        
        threshold = np.percentile(nonzero_bos, 25)
        return max(0.1, min(0.5, threshold))
    
    def _find_donors(self, mol: RawMolecule, metal_idx: int) -> Tuple[List[int], Dict[int, float], List[CoordRepIssue]]:
        """Find donor atoms."""
        issues = []
        donors = []
        bond_orders = {}
        
        threshold = self._get_bo_threshold(mol, metal_idx)
        
        if mol.bond_orders is not None:
            for i, bo in enumerate(mol.bond_orders[metal_idx]):
                if i != metal_idx and bo > threshold:
                    donors.append(i)
                    bond_orders[i] = bo
                    
                    if bo < 0.5:
                        issues.append(CoordRepIssue(
                            code=IssueCodes.LOW_BOND_ORDER,
                            message=f"Low BO ({bo:.2f}) for donor {i}",
                            severity="info",
                            context={'atom': i, 'bo': bo}
                        ))
        else:
            metal_coord = mol.atoms[metal_idx].coords
            for i, atom in enumerate(mol.atoms):
                if i == metal_idx:
                    continue
                dist = np.linalg.norm(atom.coords - metal_coord)
                if dist < 2.8:
                    donors.append(i)
                    bond_orders[i] = 1.0 / dist
        
        if len(donors) == 0:
            issues.append(CoordRepIssue(
                code=IssueCodes.NO_DONORS,
                message="No donor atoms found",
                severity="error"
            ))
        
        return donors, bond_orders, issues
    
    def _estimate_oxidation_state(self, mol: RawMolecule, metal_idx: int, 
                                  donors: List[int]) -> Optional[int]:
        """Estimate oxidation state."""
        metal = mol.atoms[metal_idx].element
        cn = len(donors)
        
        common_ox = {
            'Fe': {4: 2, 5: 2, 6: 2},
            'Co': {4: 2, 5: 2, 6: 3},
            'Ni': {4: 2, 5: 2, 6: 2},
            'Cu': {4: 2, 5: 2, 6: 2},
            'Zn': {4: 2, 5: 2, 6: 2},
            'Pd': {4: 2, 5: 2, 6: 4},
            'Pt': {4: 2, 5: 2, 6: 4},
            'Ru': {5: 2, 6: 2},
            'Rh': {4: 1, 5: 3, 6: 3},
            'Ir': {4: 1, 5: 3, 6: 3},
        }
        
        if metal in common_ox:
            return common_ox[metal].get(cn, 2)
        
        return None
    
    def _get_d_count(self, metal: str, oxidation: Optional[int]) -> Optional[int]:
        """Get d-electron count."""
        if oxidation is None:
            return None
        
        info = METAL_INFO.get(metal)
        if not info:
            return None
        
        group = info['group']
        d_count = group - oxidation
        return max(0, min(10, d_count))
    
    def _create_metal_state(self, mol: RawMolecule, metal_idx: int, 
                           donors: List[int]) -> MetalState:
        """Create metal state."""
        metal_atom = mol.atoms[metal_idx]
        element = metal_atom.element
        
        oxidation = self._estimate_oxidation_state(mol, metal_idx, donors)
        d_count = self._get_d_count(element, oxidation)
        
        info = METAL_INFO.get(element, {})
        phys_vec = np.array([
            info.get('row', 0),
            info.get('group', 0),
            info.get('en', 0),
        ])
        
        return MetalState(
            element=element,
            oxidation=oxidation,
            spin=None,
            dcount=d_count,
            phys_vec=phys_vec
        )
    
    def _build_molecular_graph(self, mol: RawMolecule) -> nx.Graph:
        """Build molecular graph."""
        G = nx.Graph()
        
        for atom in mol.atoms:
            G.add_node(atom.index, element=atom.element, coords=atom.coords)
        
        if mol.bond_orders is not None:
            n = len(mol.atoms)
            for i in range(n):
                for j in range(i + 1, n):
                    bo = mol.bond_orders[i, j]
                    if bo > 0.1:
                        G.add_edge(i, j, bond_order=bo)
        else:
            coords = mol.get_coords()
            n = len(mol.atoms)
            for i in range(n):
                for j in range(i + 1, n):
                    dist = np.linalg.norm(coords[i] - coords[j])
                    if dist < 1.8:
                        G.add_edge(i, j, bond_order=1.0)
        
        return G
    
    def build(self, mol: RawMolecule) -> GraphBuildResult:
        """Build graph for molecule."""
        all_issues = []
        
        metal_idx, issues = self._find_metal(mol)
        all_issues.extend(issues)
        
        if metal_idx is None:
            raise ValueError(f"No metal found in {mol.mol_id}")
        
        donors, bond_orders, issues = self._find_donors(mol, metal_idx)
        all_issues.extend(issues)
        
        metal_state = self._create_metal_state(mol, metal_idx, donors)
        
        mol_graph = self._build_molecular_graph(mol)
        
        return GraphBuildResult(
            metal_idx=metal_idx,
            metal_state=metal_state,
            donor_indices=donors,
            donor_bond_orders=bond_orders,
            molecular_graph=mol_graph,
            issues=all_issues
        )





