"""
tmQM Dataset Reader

Parses tmQM XYZ, BO (sparse format), and CSV property files into
RawMolecule objects suitable for CoordRep encoding.
"""

import re
import gzip
from pathlib import Path
from dataclasses import dataclass
from typing import Dict, List, Optional, Iterator, Tuple
import numpy as np


@dataclass
class Atom:
    """Atom data."""
    index: int
    element: str
    x: float
    y: float
    z: float

    @property
    def coords(self) -> np.ndarray:
        return np.array([self.x, self.y, self.z])


@dataclass
class RawMolecule:
    """Raw molecule data (unprocessed)."""
    mol_id: str
    atoms: List[Atom]
    bond_orders: Optional[np.ndarray] = None
    properties: Dict = None

    @property
    def n_atoms(self) -> int:
        return len(self.atoms)

    def get_coords(self) -> np.ndarray:
        """Get all atomic coordinates as an (N, 3) array."""
        return np.array([a.coords for a in self.atoms])

    def get_elements(self) -> List[str]:
        """Get all element symbols."""
        return [a.element for a in self.atoms]


TRANSITION_METALS = {
    # 3d
    'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
    # 4d
    'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
    # 5d
    'La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
}


class TMQMReader:
    """Reader for the tmQM dataset."""

    def __init__(self, data_dir: str):
        self.data_dir = Path(data_dir)
        self._xyz_data = {}
        self._bo_data = {}
        self._prop_data = {}
        self._loaded = False

    def _load_xyz_file(self, filepath: Path) -> Dict[str, List[Atom]]:
        """Load an XYZ file containing multiple molecules."""
        molecules = {}

        if filepath.suffix == '.gz':
            open_func = lambda p: gzip.open(p, 'rt')
        else:
            open_func = lambda p: open(p, 'r')

        with open_func(filepath) as f:
            lines = f.readlines()

        i = 0
        while i < len(lines):
            try:
                n_atoms = int(lines[i].strip())
            except (ValueError, IndexError):
                i += 1
                continue

            comment = lines[i + 1].strip()
            mol_id_match = re.search(r'CSD_code\s*=\s*(\w+)', comment)
            if not mol_id_match:
                mol_id_match = re.search(r'CSD-(\w+)', comment)

            if mol_id_match:
                mol_id = mol_id_match.group(1)
            else:
                mol_id = f"MOL_{i}"

            atoms = []
            for j in range(n_atoms):
                line = lines[i + 2 + j].strip().split()
                if len(line) >= 4:
                    atoms.append(Atom(
                        index=j,
                        element=line[0],
                        x=float(line[1]),
                        y=float(line[2]),
                        z=float(line[3])
                    ))

            if atoms:
                molecules[mol_id] = atoms

            i += n_atoms + 2

        return molecules

    def _load_bo_file(self, filepath: Path) -> Dict[str, np.ndarray]:
        """
        Load a Bond Order file in tmQM sparse format.

        Format:
            CSD_code = WELROW | 2020-2024 CSD
                 1  La  3.031        Se   4 0.429    Se   5 0.411    ...
                 2  Se  1.967        P    8 1.313    La   1 0.285    ...
        """
        molecules = {}

        if filepath.suffix == '.gz':
            open_func = lambda p: gzip.open(p, 'rt')
        else:
            open_func = lambda p: open(p, 'r')

        with open_func(filepath) as f:
            lines = f.readlines()

        current_mol_id = None
        current_bo_entries = []

        for line in lines:
            line = line.strip()
            if not line:
                continue

            if line.startswith('CSD_code'):
                if current_mol_id and current_bo_entries:
                    bo_matrix = self._build_bo_matrix(current_bo_entries)
                    if bo_matrix is not None:
                        molecules[current_mol_id] = bo_matrix

                match = re.search(r'CSD_code\s*=\s*(\w+)', line)
                if match:
                    current_mol_id = match.group(1)
                else:
                    current_mol_id = None
                current_bo_entries = []
            else:
                parts = line.split()
                if len(parts) >= 3:
                    try:
                        atom_idx = int(parts[0]) - 1  # convert to 0-based
                        atom_elem = parts[1]
                        total_bo = float(parts[2])

                        neighbors = []
                        i = 3
                        while i + 2 < len(parts):
                            neighbor_elem = parts[i]
                            try:
                                neighbor_idx = int(parts[i + 1]) - 1
                                bo_value = float(parts[i + 2])
                                neighbors.append((neighbor_idx, bo_value))
                            except (ValueError, IndexError):
                                pass
                            i += 3

                        current_bo_entries.append((atom_idx, neighbors))
                    except (ValueError, IndexError):
                        continue

        if current_mol_id and current_bo_entries:
            bo_matrix = self._build_bo_matrix(current_bo_entries)
            if bo_matrix is not None:
                molecules[current_mol_id] = bo_matrix

        return molecules

    def _build_bo_matrix(self, entries: List[Tuple[int, List[Tuple[int, float]]]]) -> Optional[np.ndarray]:
        """Build a dense BO matrix from sparse entries."""
        if not entries:
            return None

        max_idx = 0
        for atom_idx, neighbors in entries:
            max_idx = max(max_idx, atom_idx)
            for neighbor_idx, _ in neighbors:
                max_idx = max(max_idx, neighbor_idx)

        n = max_idx + 1
        bo_matrix = np.zeros((n, n), dtype=np.float32)

        for atom_idx, neighbors in entries:
            for neighbor_idx, bo_value in neighbors:
                if 0 <= atom_idx < n and 0 <= neighbor_idx < n:
                    bo_matrix[atom_idx, neighbor_idx] = bo_value
                    bo_matrix[neighbor_idx, atom_idx] = bo_value

        return bo_matrix

    def _load_properties(self, filepath: Path) -> Dict[str, Dict]:
        """Load property CSV."""
        import csv

        properties = {}
        with open(filepath, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                mol_id = row.get('CSD_code', row.get('csd_code', ''))
                if mol_id:
                    properties[mol_id] = {
                        'hl_gap': float(row.get('Electronic_E_gap_eV', row.get('hl_gap', 0)) or 0),
                        'dipole': float(row.get('Dipole_M_Debye', row.get('dipole_moment', 0)) or 0),
                        'homo': float(row.get('HOMO_eV', row.get('homo_energy', 0)) or 0),
                        'lumo': float(row.get('LUMO_eV', row.get('lumo_energy', 0)) or 0),
                    }

        return properties

    def load(self):
        """Load all XYZ, BO, and property data."""
        if self._loaded:
            return

        for pattern in ['tmQM_X*.xyz', '*.xyz']:
            for fp in sorted(self.data_dir.glob(pattern)):
                if '.gz' not in fp.name:
                    self._xyz_data.update(self._load_xyz_file(fp))

        for pattern in ['tmQM_X*.BO', '*.BO']:
            for fp in sorted(self.data_dir.glob(pattern)):
                if '.gz' not in fp.name:
                    self._bo_data.update(self._load_bo_file(fp))

        for pattern in ['tmQM_y.csv', '*.csv']:
            for fp in self.data_dir.glob(pattern):
                self._prop_data.update(self._load_properties(fp))

        self._loaded = True
        print(f"Loaded: {len(self._xyz_data)} XYZ, {len(self._bo_data)} BO, {len(self._prop_data)} properties")

    def get_molecule(self, mol_id: str) -> Optional[RawMolecule]:
        """Get a single molecule by ID."""
        self.load()

        if mol_id not in self._xyz_data:
            return None

        atoms = self._xyz_data[mol_id]
        bo = self._bo_data.get(mol_id)
        props = self._prop_data.get(mol_id, {})

        return RawMolecule(
            mol_id=mol_id,
            atoms=atoms,
            bond_orders=bo,
            properties=props
        )

    def iter_molecules(self, limit: int = None) -> Iterator[RawMolecule]:
        """Iterate over all molecules."""
        self.load()

        count = 0
        for mol_id in self._xyz_data:
            if limit and count >= limit:
                break

            mol = self.get_molecule(mol_id)
            if mol:
                yield mol
                count += 1

    def get_mol_ids(self) -> List[str]:
        """Get all molecule IDs."""
        self.load()
        return list(self._xyz_data.keys())

    def __len__(self) -> int:
        self.load()
        return len(self._xyz_data)


def find_metal_index(atoms: List[Atom]) -> Optional[int]:
    """Find the index of the transition-metal atom."""
    for i, atom in enumerate(atoms):
        if atom.element in TRANSITION_METALS:
            return i
    return None


def find_all_metal_indices(atoms: List[Atom]) -> List[int]:
    """Find indices of all transition-metal atoms."""
    return [i for i, atom in enumerate(atoms) if atom.element in TRANSITION_METALS]
