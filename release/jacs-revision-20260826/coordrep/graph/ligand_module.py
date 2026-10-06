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
            
            (
                payload,
                mol_weight,
                payload_provenance,
                connectivity_status,
                rdkit_issue,
            ) = self._generate_smiles(
                lig_atoms,
                mol_graph,
                comp,
                has_bond_orders=mol.bond_orders is not None,
            )
            if rdkit_issue:
                issues.append(rdkit_issue)

            attachment_set_key = None
            donor_attachment_keys = []
            if payload_provenance == "SMILES":
                (
                    attachment_set_key,
                    donor_key_by_global_index,
                ) = self._generate_attachment_descriptors(
                    lig_atoms,
                    mol_graph,
                    comp,
                    donors,
                )
                donor_attachment_keys = [
                    donor_key_by_global_index.get(donor_idx, "")
                    for donor_idx in donors
                ]
            
            donor_elements = [mol.atoms[d].element for d in donors]
            
            eta = self._detect_hapticity(donors, lig_atoms, mol_graph, comp)
            
            ligand = LigandModule(
                lig_id=f"L{lig_counter}",
                smiles=payload,
                attach_atoms=donors, #
                donor_elements=donor_elements,
                dent=len(donors),
                attachment_set_key=attachment_set_key,
                donor_attachment_keys=donor_attachment_keys,
                eta=eta,
                charge=None, # v1
                meta={'mw': mol_weight, 'atoms': list(comp)},
                payload_provenance=payload_provenance,
                connectivity_status=connectivity_status,
            )
            
            ligands.append(ligand)
            lig_counter += 1
        
        return ligands, issues
    
    def _generate_smiles(self, atoms: List[Atom], 
                        full_graph: nx.Graph,
                        atom_indices: set,
                        has_bond_orders: bool) -> Tuple[
                            str, float, str, str, Optional[CoordRepIssue]
                        ]:
        """Generate a canonical ligand payload with explicit provenance.

        Coordinate-only graphs can provide a useful nominal SMILES payload,
        but their bond orders were inferred rather than supplied.  Such
        payloads are marked low confidence.  If SMILES generation fails, the
        fallback is explicitly a molecular formula and is never presented as
        connectivity.
        """
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
            
            formula = ''.join(parts)
            issue = CoordRepIssue(
                code=IssueCodes.FORMULA_FALLBACK,
                message="RDKit unavailable; emitted formula-only ligand payload",
                severity="warning",
                context={
                    "payload_provenance": "FORMULA",
                    "connectivity_status": "unsupported_formula_fallback",
                },
            )
            return (
                formula,
                mol_weight,
                "FORMULA",
                "unsupported_formula_fallback",
                issue,
            )
        
        try:
            mol, _ = self._build_rdkit_component(atoms, full_graph, atom_indices)
            
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
                    severity="warning",
                    context={
                        "payload_provenance": "SMILES",
                        "connectivity_status": "low_confidence_sanitize_warning",
                    },
                )
            
            smiles = Chem.MolToSmiles(mol, canonical=True)
            mol_weight = Descriptors.MolWt(mol) if mol else 0.0
            
            if issue is not None:
                connectivity_status = "low_confidence_sanitize_warning"
            elif has_bond_orders:
                connectivity_status = "supported_bond_order_graph"
            else:
                connectivity_status = "low_confidence_distance_inferred"

            return (
                smiles,
                mol_weight,
                "SMILES",
                connectivity_status,
                issue,
            )
            
        except Exception as e:
            # Fallback
            issue = CoordRepIssue(
                code=IssueCodes.FORMULA_FALLBACK,
                message=f"SMILES generation failed; emitted formula only: {str(e)}",
                severity="warning",
                context={
                    "payload_provenance": "FORMULA",
                    "connectivity_status": "unsupported_formula_fallback",
                },
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
            
            return (
                ''.join(parts),
                0.0,
                "FORMULA",
                "unsupported_formula_fallback",
                issue,
            )

    def _build_rdkit_component(self,
                               atoms: List[Atom],
                               full_graph: nx.Graph,
                               atom_indices: set):
        """Build an RDKit ligand graph and its global-to-local atom map."""
        rw_mol = Chem.RWMol()
        old_to_new = {}

        for atom in atoms:
            new_idx = rw_mol.AddAtom(Chem.Atom(atom.element))
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

        return rw_mol.GetMol(), old_to_new

    def _generate_attachment_descriptors(self,
                                         atoms: List[Atom],
                                         full_graph: nx.Graph,
                                         atom_indices: set,
                                         donor_indices: List[int]):
        """Return exact ligand-local attachment-set and donor-orbit keys.

        RDKit canonical SMILES is evaluated on two coloured versions of the
        same ligand graph.  Atom-map 1 marks every coordinated atom.  For a
        per-donor key, atom-map 2 singles out the queried donor while all
        other coordinated atoms retain map 1.  Consequently, two donors get
        the same key exactly when an automorphism preserving the complete
        attachment set can exchange them.  Global/source atom numbers are
        never serialized into either key.
        """
        if not RDKIT_AVAILABLE:
            return None, {}

        try:
            mol, old_to_new = self._build_rdkit_component(
                atoms, full_graph, atom_indices
            )
            try:
                if self.strict_sanitize:
                    Chem.SanitizeMol(mol)
                else:
                    Chem.SanitizeMol(
                        mol,
                        sanitizeOps=(
                            Chem.SanitizeFlags.SANITIZE_ALL
                            ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE
                        ),
                    )
            except Exception:
                # MolToSmiles can still canonicalize some graphs that fail a
                # strict valence check; if it cannot, the outer guard returns
                # unresolved metadata rather than a misleading graph key.
                pass

            local_donors = [old_to_new[idx] for idx in donor_indices]

            def coloured_smiles(target_local_idx=None):
                coloured = Chem.Mol(mol)
                for local_idx in local_donors:
                    coloured.GetAtomWithIdx(local_idx).SetAtomMapNum(1)
                if target_local_idx is not None:
                    coloured.GetAtomWithIdx(target_local_idx).SetAtomMapNum(2)
                return Chem.MolToSmiles(
                    coloured,
                    canonical=True,
                    isomericSmiles=True,
                )

            attachment_set_key = coloured_smiles()
            donor_keys = {
                global_idx: coloured_smiles(old_to_new[global_idx])
                for global_idx in donor_indices
            }
            return attachment_set_key, donor_keys
        except Exception:
            return None, {}
    
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



