"""
datasets.py
===========
Unified graph dataset construction for GNN baselines.

Builds three types of graph data objects:
1. LigandGraph — 2D molecular graph of a ligand (atoms + bonds)
2. CoordinationGraph — metal-centred graph (metal node + donor nodes + ligand subgraphs)
3. ComplexGraph3D — full 3D structure from tmQM XYZ

All datasets share the same train/val/test split as the CoordRep MLM,
derived from the fig5_tasks IDs (test split) and pipeline data.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

try:
    import torch
    from torch_geometric.data import Data, InMemoryDataset
except ImportError:
    raise ImportError("torch and torch_geometric required. Install via pip.")

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors
except ImportError:
    Chem = None


# ── Atom featurization ────────────────────────────────────

ATOM_FEATURES = ['C', 'N', 'O', 'S', 'P', 'F', 'Cl', 'Br', 'I', 'Se',
                 'Si', 'B', 'H', 'other']
HYBRIDIZATION = [
    Chem.rdchem.HybridizationType.SP,
    Chem.rdchem.HybridizationType.SP2,
    Chem.rdchem.HybridizationType.SP3,
    Chem.rdchem.HybridizationType.SP3D,
    Chem.rdchem.HybridizationType.SP3D2,
] if Chem else []


def one_hot(val, choices):
    encoding = [0] * (len(choices) + 1)
    try:
        idx = choices.index(val)
        encoding[idx] = 1
    except ValueError:
        encoding[-1] = 1
    return encoding


def atom_features(atom) -> List[float]:
    """Compute atom feature vector from RDKit atom."""
    symbol = atom.GetSymbol()
    feats = one_hot(symbol, ATOM_FEATURES)
    feats += one_hot(atom.GetHybridization(), HYBRIDIZATION)
    feats += [
        atom.GetDegree() / 6.0,
        atom.GetFormalCharge() / 3.0,
        atom.GetNumRadicalElectrons() / 2.0,
        int(atom.GetIsAromatic()),
        int(atom.IsInRing()),
        atom.GetTotalNumHs() / 4.0,
    ]
    return feats


def mol_to_graph(mol, donor_indices: Optional[List[int]] = None) -> Optional[Data]:
    """Convert RDKit Mol to PyG Data object with donor labels."""
    if mol is None:
        return None
    n_atoms = mol.GetNumAtoms()
    if n_atoms == 0:
        return None

    # Node features
    x = []
    for atom in mol.GetAtoms():
        x.append(atom_features(atom))
    x = torch.tensor(x, dtype=torch.float)

    # Edge index (undirected)
    edge_index = []
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        edge_index.append([i, j])
        edge_index.append([j, i])

    if edge_index:
        edge_index = torch.tensor(edge_index, dtype=torch.long).t().contiguous()
    else:
        edge_index = torch.zeros((2, 0), dtype=torch.long)

    # Donor labels (binary per atom)
    y = torch.zeros(n_atoms, dtype=torch.float)
    if donor_indices:
        for idx in donor_indices:
            if 0 <= idx < n_atoms:
                y[idx] = 1.0

    data = Data(x=x, edge_index=edge_index, y=y)
    data.num_nodes = n_atoms
    return data


# ── SMILES → graph with donor identification ─────────────

def smiles_to_graph_with_donors(
    smiles: str,
    donor_elements: List[str],
    denticity: int,
) -> Optional[Data]:
    """
    Build a ligand graph from SMILES or molecular formula.
    Always uses formula_to_graph for consistent feature dim across dataset.
    """
    return formula_to_graph(smiles, donor_elements, denticity)


def formula_to_graph(
    formula: str,
    donor_elements: List[str],
    denticity: int,
) -> Optional[Data]:
    """
    Build a pseudo-graph from a molecular formula (e.g. C5H4NO).
    Each element becomes one node; edges connect all heavy atoms.
    Donor labels are assigned to elements matching donor_elements.
    """
    import re as _re
    # Parse formula: C5H4NO -> [(C,5),(H,4),(N,1),(O,1)]
    parts = _re.findall(r'([A-Z][a-z]?)(\d*)', formula)
    atoms = []
    for elem, cnt in parts:
        if not elem:
            continue
        n = int(cnt) if cnt else 1
        atoms.extend([elem] * n)

    if not atoms:
        return None

    n_atoms = len(atoms)
    ELEM_LIST = ['C', 'N', 'O', 'S', 'P', 'F', 'Cl', 'Br', 'I', 'Se', 'H', 'other']

    x = []
    for a in atoms:
        feat = one_hot(a, ELEM_LIST)
        feat.append(0.0)  # placeholder for degree
        x.append(feat)
    x = torch.tensor(x, dtype=torch.float)

    # Edges: connect all heavy atoms sequentially
    heavy = [i for i, a in enumerate(atoms) if a != 'H']
    edge_index = []
    for i in range(len(heavy) - 1):
        edge_index.append([heavy[i], heavy[i + 1]])
        edge_index.append([heavy[i + 1], heavy[i]])
    if edge_index:
        edge_index = torch.tensor(edge_index, dtype=torch.long).t().contiguous()
    else:
        edge_index = torch.zeros((2, 0), dtype=torch.long)

    # Donor labels
    y = torch.zeros(n_atoms, dtype=torch.float)
    donor_count = {}
    for d in donor_elements:
        donor_count[d] = donor_count.get(d, 0) + 1
    assigned = {}
    for i, a in enumerate(atoms):
        if a in donor_count and donor_count[a] > 0 and assigned.get(a, 0) < donor_count[a]:
            y[i] = 1.0
            assigned[a] = assigned.get(a, 0) + 1

    data = Data(x=x, edge_index=edge_index, y=y)
    data.num_nodes = n_atoms
    data.smiles = formula
    data.denticity = denticity
    return data


# ── Pipeline data loading ─────────────────────────────────

PIPELINE_PATH = "inputs/tmqm/results.jsonl"
FIG5_PATH = "inputs/property_benchmarks/fig5_tasks.jsonl"
TMQM_XYZ_DIR = "inputs/tmQM"
CSD_PATH = "inputs/csd/csd_retained_entries.jsonl"


def load_pipeline_data(max_n: Optional[int] = None) -> List[dict]:
    """Load processed CoordRep pipeline entries."""
    data = []
    with open(PIPELINE_PATH) as f:
        for line in f:
            d = json.loads(line)
            data.append(d)
            if max_n and len(data) >= max_n:
                break
    return data


def get_test_ids() -> set:
    """Get test split mol_ids from fig5_tasks."""
    ids = set()
    with open(FIG5_PATH) as f:
        for line in f:
            d = json.loads(line)
            mol_id = d.get('id', '')
            # Extract mol_id: tmqm_0_donor_143 -> tmqm_0
            match = re.match(r'^(tmqm_\d+)', mol_id)
            if match:
                ids.add(match.group(1))
    return ids


def deterministic_split(mol_id: str, train_frac=0.8, val_frac=0.1) -> str:
    """Deterministic split based on hash of mol_id."""
    h = int(hashlib.md5(mol_id.encode()).hexdigest(), 16) % 1000
    if h < train_frac * 1000:
        return 'train'
    elif h < (train_frac + val_frac) * 1000:
        return 'val'
    else:
        return 'test'


def build_unified_split(entries: List[dict]) -> Dict[str, List[dict]]:
    """
    Split pipeline data into train/val/test.
    Test set = fig5 mol_ids; rest split 80/10 by hash.
    """
    test_ids = get_test_ids()
    splits = {'train': [], 'val': [], 'test': []}

    for entry in entries:
        mol_id = entry.get('mol_id', '')
        # fig5 test set uses tmqm_{idx}
        mol_key = f"tmqm_{mol_id}" if not mol_id.startswith('tmqm') else mol_id
        if mol_key in test_ids or mol_id in test_ids:
            splits['test'].append(entry)
        else:
            sp = deterministic_split(mol_id)
            splits[sp].append(entry)

    return splits


# ── Ligand Graph Dataset ──────────────────────────────────

class LigandGraphDataset:
    """
    Constructs ligand-level graph data with donor annotation labels.

    Each sample is one ligand in one complex:
    - Graph: 2D molecular graph from SMILES
    - Labels: binary donor annotation per atom
    - Meta: complex_id, metal, CN, denticity, source
    """

    def __init__(self, entries: List[dict], split: str = 'all'):
        self.samples = []
        self._build(entries, split)

    def _build(self, entries, split):
        for entry in entries:
            mol_id = entry.get('mol_id', '')
            metal = entry.get('metal', '?')
            cn = entry.get('cn', 0)
            lig_smiles_list = entry.get('ligand_smiles', [])
            lig_dents = entry.get('ligand_dents', [])
            donor_elements = entry.get('donor_elements', [])

            if not lig_smiles_list:
                continue

            # Map donors to ligands
            donor_offset = 0
            for lig_idx, (smi, dent) in enumerate(zip(lig_smiles_list, lig_dents)):
                if dent <= 0:
                    continue
                lig_donors = donor_elements[donor_offset:donor_offset + dent]
                donor_offset += dent

                self.samples.append({
                    'complex_id': mol_id,
                    'ligand_id': f'L{lig_idx+1}',
                    'smiles': smi,
                    'denticity': dent,
                    'donor_elements': lig_donors,
                    'metal': metal,
                    'cn': cn,
                    'source': 'tmQM',
                })

    def __len__(self):
        return len(self.samples)

    def __getitem__(self, idx):
        return self.samples[idx]

    def to_pyg_list(self) -> List[Data]:
        """Convert to list of PyG Data objects."""
        pyg_data = []
        for sample in self.samples:
            data = smiles_to_graph_with_donors(
                sample['smiles'],
                sample['donor_elements'],
                sample['denticity'],
            )
            if data is not None:
                data.complex_id = sample['complex_id']
                data.metal = sample['metal']
                data.cn = sample['cn']
                data.denticity_val = sample['denticity']
                data.ligand_id = sample['ligand_id']
                pyg_data.append(data)
        return pyg_data


# ── 3D Complex Dataset (from tmQM XYZ) ───────────────────

def parse_tmqm_xyz(xyz_path: str, target_mol_ids: Optional[set] = None) -> Dict[str, dict]:
    """
    Parse tmQM XYZ file, returning {mol_id: {atoms: [...], coords: [...]}}.
    """
    molecules = {}
    with open(xyz_path) as f:
        while True:
            line = f.readline()
            if not line:
                break
            try:
                n_atoms = int(line.strip())
            except ValueError:
                continue
            comment = f.readline().strip()
            # Extract CSD code
            match = re.search(r'CSD_code\s*=\s*(\S+)', comment)
            mol_id = match.group(1) if match else f'unk_{len(molecules)}'

            if target_mol_ids and mol_id not in target_mol_ids:
                # Skip
                for _ in range(n_atoms):
                    f.readline()
                continue

            atoms = []
            coords = []
            for _ in range(n_atoms):
                parts = f.readline().split()
                if len(parts) >= 4:
                    atoms.append(parts[0])
                    coords.append([float(parts[1]), float(parts[2]), float(parts[3])])

            molecules[mol_id] = {'atoms': atoms, 'coords': coords, 'n_atoms': n_atoms}

    return molecules


def build_3d_graph(mol_data: dict, metal_element: str,
                   donor_elements: Optional[List[str]] = None) -> Optional[Data]:
    """Build 3D graph from XYZ coordinates with metal-centric edges."""
    atoms = mol_data['atoms']
    coords = np.array(mol_data['coords'])
    n_atoms = len(atoms)

    if n_atoms == 0:
        return None

    # Find metal atom
    metal_idx = None
    for i, a in enumerate(atoms):
        if a == metal_element:
            metal_idx = i
            break
    if metal_idx is None:
        return None

    # Simple distance-based edges (within 3.0 Å)
    cutoff = 3.0
    edge_index = []
    for i in range(n_atoms):
        for j in range(i + 1, n_atoms):
            dist = np.linalg.norm(coords[i] - coords[j])
            if dist < cutoff:
                edge_index.append([i, j])
                edge_index.append([j, i])

    # Node features: one-hot element + is_metal flag
    ELEMENTS_3D = ['H', 'C', 'N', 'O', 'S', 'P', 'F', 'Cl', 'Br', 'Se'] + [metal_element]
    x = []
    for i, a in enumerate(atoms):
        feat = one_hot(a, ELEMENTS_3D)
        feat.append(1.0 if i == metal_idx else 0.0)  # is_metal
        x.append(feat)

    x = torch.tensor(x, dtype=torch.float)
    pos = torch.tensor(coords, dtype=torch.float)

    if edge_index:
        edge_index_t = torch.tensor(edge_index, dtype=torch.long).t().contiguous()
    else:
        edge_index_t = torch.zeros((2, 0), dtype=torch.long)

    # Donor labels: atoms bonded to metal with matching element
    y = torch.zeros(n_atoms, dtype=torch.float)
    metal_coord = coords[metal_idx]
    for i in range(n_atoms):
        if i == metal_idx:
            continue
        dist = np.linalg.norm(coords[i] - metal_coord)
        if dist < 2.8 and atoms[i] != 'H':  # within coordination shell
            if donor_elements is None or atoms[i] in donor_elements:
                y[i] = 1.0

    data = Data(x=x, edge_index=edge_index_t, pos=pos, y=y)
    data.num_nodes = n_atoms
    data.metal_idx = metal_idx
    return data


# ── Leakage diagnostics ──────────────────────────────────

def leakage_diagnostic(splits: Dict[str, List[dict]]) -> dict:
    """Check for mol_id leakage between splits."""
    train_ids = {e['mol_id'] for e in splits['train']}
    val_ids = {e['mol_id'] for e in splits['val']}
    test_ids = {e['mol_id'] for e in splits['test']}

    tv_leak = train_ids & val_ids
    tt_leak = train_ids & test_ids
    vt_leak = val_ids & test_ids

    # Also check SMILES leakage
    train_smi = set()
    for e in splits['train']:
        for s in e.get('ligand_smiles', []):
            train_smi.add(s)
    test_smi = set()
    for e in splits['test']:
        for s in e.get('ligand_smiles', []):
            test_smi.add(s)
    smi_overlap = train_smi & test_smi

    return {
        'train_val_mol_id_leak': len(tv_leak),
        'train_test_mol_id_leak': len(tt_leak),
        'val_test_mol_id_leak': len(vt_leak),
        'train_size': len(splits['train']),
        'val_size': len(splits['val']),
        'test_size': len(splits['test']),
        'ligand_smiles_train_test_overlap': len(smi_overlap),
        'ligand_smiles_train_total': len(train_smi),
        'ligand_smiles_test_total': len(test_smi),
    }
