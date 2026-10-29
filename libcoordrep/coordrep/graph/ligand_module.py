"""
Ligand Module

SMILES + attachment points + metadata
"""

import networkx as nx
from typing import List, Dict, Optional, Tuple
from collections import defaultdict

from ..core import LigandModule, CoordRepIssue, IssueCodes
from ..io.tmqm_reader import RawMolecule, Atom


try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False


class LigandExtractor:
    """"""
    
    def __init__(self, strict_sanitize: bool = False):
        self.strict_sanitize = strict_sanitize
    
    def extract_ligands(self, 
                       mol: RawMolecule,
                       mol_graph: nx.Graph,
                       metal_idx: int,
                       donor_indices: List[int]) -> Tuple[List[LigandModule], List[CoordRepIssue]]:
        """
        
        1.
        2. →
 3. canonical SMILES
 4. attachment points (donors)
        """
        issues = []
        
        ligand_graph = mol_graph.copy()
        if metal_idx in ligand_graph:
            ligand_graph.remove_node(metal_idx)
        
        components = list(nx.connected_components(ligand_graph))
        
        donor_to_component = {}
        for donor_idx in donor_indices:
            for i, comp in enumerate(components):
                if donor_idx in comp:
                    donor_to_component[donor_idx] = i
                    break
        
        component_donors = defaultdict(list)
        for donor_idx, comp_idx in donor_to_component.items():
            component_donors[comp_idx].append(donor_idx)
        
        ligands = []
        lig_counter = 1
        
        for comp_idx, donors in sorted(component_donors.items()):
            comp = components[comp_idx]
            
            lig_atoms = [mol.atoms[i] for i in sorted(comp)]
            
            smiles, mol_weight, rdkit_issue = self._generate_smiles(lig_atoms, mol_graph, comp)
            if rdkit_issue:
                issues.append(rdkit_issue)
            
            donor_elements = [mol.atoms[d].element for d in donors]
            
            eta = self._detect_hapticity(donors, lig_atoms, mol_graph, comp)
            
            ligand = LigandModule(
                lig_id=f"L{lig_counter}",
                smiles=smiles,
                attach_atoms=donors, #
                donor_elements=donor_elements,
                dent=len(donors),
                eta=eta,
                charge=None, # v1
                meta={'mw': mol_weight, 'atoms': list(comp)}
            )
            
            ligands.append(ligand)
            lig_counter += 1
        
        return ligands, issues
    
    def _generate_smiles(self, atoms: List[Atom], 
                        full_graph: nx.Graph,
                        atom_indices: set) -> Tuple[str, float, Optional[CoordRepIssue]]:
        """Generate canonical SMILES."""
        issue = None
        mol_weight = 0.0
        
        if not RDKIT_AVAILABLE:
            element_counts = defaultdict(int)
            for atom in atoms:
                element_counts[atom.element] += 1
            
            parts = []
            for elem in sorted(element_counts.keys()):
                count = element_counts[elem]
                if count == 1:
                    parts.append(elem)
                else:
                    parts.append(f"{elem}{count}")
            
            smiles = ''.join(parts)
            return smiles, mol_weight, issue
        
        try:
            rw_mol = Chem.RWMol()
            
            old_to_new = {}
            
            for i, atom in enumerate(atoms):
                rdkit_atom = Chem.Atom(atom.element)
                new_idx = rw_mol.AddAtom(rdkit_atom)
                old_to_new[atom.index] = new_idx
            
            added_bonds = set()
            for old_idx in atom_indices:
                if old_idx not in old_to_new:
                    continue
                new_i = old_to_new[old_idx]
                
                for neighbor in full_graph.neighbors(old_idx):
                    if neighbor not in old_to_new:
                        continue
                    new_j = old_to_new[neighbor]
                    
                    bond_key = (min(new_i, new_j), max(new_i, new_j))
                    if bond_key in added_bonds:
                        continue
                    
                    bo = full_graph[old_idx][neighbor].get('bond_order', 1.0)
                    
                    if bo > 2.5:
                        bond_type = Chem.BondType.TRIPLE
                    elif bo > 1.5:
                        bond_type = Chem.BondType.DOUBLE
                    else:
                        bond_type = Chem.BondType.SINGLE
                    
                    rw_mol.AddBond(new_i, new_j, bond_type)
                    added_bonds.add(bond_key)
            
            mol = rw_mol.GetMol()
            
            # Sanitize
            try:
                if self.strict_sanitize:
                    Chem.SanitizeMol(mol)
                else:
                    Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ 
                                    Chem.SanitizeFlags.SANITIZE_KEKULIZE)
            except Exception as e:
                issue = CoordRepIssue(
                    code=IssueCodes.RDKIT_SANITIZE_FAILED,
                    message=f"RDKit sanitize failed: {str(e)}",
                    severity="warning"
                )
            
            smiles = Chem.MolToSmiles(mol, canonical=True)
            mol_weight = Descriptors.MolWt(mol) if mol else 0.0
            
            return smiles, mol_weight, issue
            
        except Exception as e:
            # Fallback
            issue = CoordRepIssue(
                code=IssueCodes.RDKIT_SANITIZE_FAILED,
                message=f"SMILES generation failed: {str(e)}",
                severity="warning"
            )
            
            element_counts = defaultdict(int)
            for atom in atoms:
                element_counts[atom.element] += 1
            
            parts = []
            for elem in sorted(element_counts.keys()):
                count = element_counts[elem]
                if count == 1:
                    parts.append(elem)
                else:
                    parts.append(f"{elem}{count}")
            
            return ''.join(parts), 0.0, issue
    
    def _detect_hapticity(self, donors: List[int], 
                         atoms: List[Atom],
                         full_graph: nx.Graph,
                         atom_indices: set) -> Optional[int]:
        """Detect hapticity."""
        if len(donors) < 2:
            return None
        
        donor_set = set(donors)
        
        subgraph = full_graph.subgraph(atom_indices)
        
        all_adjacent = True
        for i, d1 in enumerate(donors):
            for d2 in donors[i+1:]:
                if d2 not in subgraph.neighbors(d1):
                    all_adjacent = False
                    break
            if not all_adjacent:
                break
        
        if all_adjacent and len(donors) > 1:
            donor_elements = [full_graph.nodes[d].get('element', '') for d in donors]
            if all(e == 'C' for e in donor_elements):
                return len(donors)
        
        return None





