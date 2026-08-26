"""
Donor Site Identification - Sprint 2 Version

1. donor_id canonical SMILES +
2. RDKit atom index
3. RDKit canonical ranking

: {LigID}:{Element}:{donor_order}:{canonical_rank}
: L1:N:1:42, L2:Cl:1:0
"""

import numpy as np
from typing import List, Dict, Tuple
from dataclasses import dataclass
from collections import defaultdict

from ..core import DonorSite, LigandModule
from ..io.tmqm_reader import RawMolecule

try:
    from rdkit import Chem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False


ATOMIC_NUMBERS = {
    'H': 1, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'P': 15, 'S': 16, 'Cl': 17,
    'Br': 35, 'I': 53, 'Se': 34, 'Te': 52, 'As': 33,
}


def compute_donor_canonical_key(mol: RawMolecule,
                                donor_idx: int,
                                metal_idx: int,
                                ligand_smiles: str = None) -> tuple:
    """
    donor canonical rank
    
    1.
 2. multiset
 3.
    """
    atom = mol.atoms[donor_idx]
    element = atom.element
    
    metal_coord = mol.atoms[metal_idx].coords
    distance = np.linalg.norm(atom.coords - metal_coord)
    dist_key = int(round(distance * 10000))
    
    atom_z = ATOMIC_NUMBERS.get(element, 0)
    
    bo = 0.0
    if mol.bond_orders is not None:
        bo = mol.bond_orders[metal_idx, donor_idx]
    bo_key = int(round(bo * 1000))
    
    neighbor_elements = []
    if mol.bond_orders is not None:
        for i, b in enumerate(mol.bond_orders[donor_idx]):
            if i != metal_idx and i != donor_idx and b > 0.1:
                neighbor_elements.append(ATOMIC_NUMBERS.get(mol.atoms[i].element, 0))

    return dist_key, atom_z, bo_key, tuple(sorted(neighbor_elements))


def compute_donor_canonical_rank(mol: RawMolecule,
                                 donor_idx: int,
                                 metal_idx: int,
                                 ligand_smiles: str = None) -> int:
    """Backward-compatible integer view of the local geometric key."""
    dist_key, atom_z, bo_key, neighbor_elements = compute_donor_canonical_key(
        mol, donor_idx, metal_idx, ligand_smiles
    )
    env_hash = 0
    for z in neighbor_elements:
        env_hash = env_hash * 100 + z
    return (
        dist_key * (10**12)
        + atom_z * (10**9)
        + bo_key * (10**6)
        + (env_hash % (10**6))
    )


def assign_donor_ranks(mol: RawMolecule,
                      ligand: LigandModule,
                      metal_idx: int) -> Dict[int, DonorSite]:
    """
    donor canonical rank
    
    donor_idx -> DonorSite
    """
    donor_indices = ligand.attach_atoms
    
    if not donor_indices:
        return {}
    
    donor_geometry_keys = {
        donor_idx: compute_donor_canonical_key(
            mol, donor_idx, metal_idx, ligand.smiles
        )
        for donor_idx in donor_indices
    }
    donor_topology_keys = {
        donor_idx: ligand.donor_attachment_keys[position]
        for position, donor_idx in enumerate(donor_indices)
        if position < len(ligand.donor_attachment_keys)
        and ligand.donor_attachment_keys[position]
    }
    
    element_groups = defaultdict(list)
    for donor_idx in donor_indices:
        element = mol.atoms[donor_idx].element
        element_groups[element].append((
            donor_idx,
            donor_topology_keys.get(donor_idx, ""),
            donor_geometry_keys[donor_idx],
        ))
    
    result = {}
    
    ligand.donor_rank_orbits = []
    for element, donor_rows in element_groups.items():
        # The ligand-local topological key is primary.  Geometry only orders
        # atoms that remain in the same attachment-set stabilizer orbit.
        sorted_donors = sorted(donor_rows, key=lambda x: (x[1], x[2]))
        tie_groups = defaultdict(list)

        for donor_order, (donor_idx, topology_key, geometry_key) in enumerate(
            sorted_donors, 1
        ):
            result[donor_idx] = DonorSite(
                lig_id=ligand.lig_id,
                donor_element=element,
                donor_rank=donor_order,
            )
            if topology_key:
                tie_groups[(topology_key, geometry_key)].append(donor_order)

        for tied_ranks in tie_groups.values():
            if len(tied_ranks) > 1:
                ligand.donor_rank_orbits.append([
                    (element, rank) for rank in tied_ranks
                ])
    
    return result


def build_all_donor_sites(mol: RawMolecule,
                         ligands: List[LigandModule],
                         metal_idx: int,
                         mol_graph=None) -> Dict[int, DonorSite]:
    """
    donor site
    """
    all_sites = {}
    
    for ligand in ligands:
        sites = assign_donor_ranks(mol, ligand, metal_idx)
        all_sites.update(sites)
    
    return all_sites


def validate_donor_sites(donor_sites: Dict[int, DonorSite],
                        mol: RawMolecule) -> List[str]:
    """Validate donor sites."""
    issues = []
    
    for donor_idx, site in donor_sites.items():
        actual_element = mol.atoms[donor_idx].element
        
        if site.donor_element != actual_element:
            issues.append(f"Element mismatch: {site} has element {actual_element}")
        
        if actual_element == 'H':
            issues.append(f"H cannot be a donor: {site}")
    
    return issues
