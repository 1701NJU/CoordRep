"""
CIF Reader for COD (Crystallography Open Database)

Extracts atomic coordinates from CIF files and converts fractional
coordinates to Cartesian coordinates for CoordRep encoding.
"""

import re
import math
import numpy as np
from pathlib import Path
from dataclasses import dataclass
from typing import List, Dict, Optional, Iterator, Tuple

from .tmqm_reader import Atom, RawMolecule, TRANSITION_METALS


@dataclass
class CIFData:
    """Parsed CIF file data."""
    cod_id: str
    formula: str
    cell_params: Dict[str, float]  # a, b, c, alpha, beta, gamma
    atoms: List[Atom]
    space_group: str = ""


class CIFReader:
    """Reader for CIF files."""

    def __init__(self, cif_dir: str):
        self.cif_dir = Path(cif_dir)
        self.files = []
        self._loaded = False

    def load(self, limit: int = None):
        """Load the list of CIF files."""
        if self._loaded:
            return

        self.files = sorted(self.cif_dir.glob("*.cif"))
        if limit:
            self.files = self.files[:limit]

        self._loaded = True
        print(f"Found {len(self.files)} CIF files")

    def parse_cif(self, filepath: Path) -> Optional[CIFData]:
        """Parse a single CIF file."""
        try:
            with open(filepath, 'r', errors='ignore') as f:
                content = f.read()

            cod_id = filepath.stem

            formula_match = re.search(r"_chemical_formula_sum\s+['\"]?([^'\"\n]+)['\"]?", content)
            formula = formula_match.group(1).strip() if formula_match else ""

            if not any(m in formula for m in TRANSITION_METALS):
                return None

            cell_params = {}
            for param in ['a', 'b', 'c', 'alpha', 'beta', 'gamma']:
                pattern = rf"_cell_(?:length_{param}|angle_{param})\s+([\d.]+)"
                match = re.search(pattern, content)
                if match:
                    cell_params[param] = float(match.group(1))

            required = ['a', 'b', 'c', 'alpha', 'beta', 'gamma']
            if not all(p in cell_params for p in required):
                return None

            atoms = self._parse_atoms(content, cell_params)
            if not atoms:
                return None

            sg_match = re.search(r"_symmetry_space_group_name_H-M\s+['\"]?([^'\"\n]+)['\"]?", content)
            space_group = sg_match.group(1).strip() if sg_match else ""

            return CIFData(
                cod_id=cod_id,
                formula=formula,
                cell_params=cell_params,
                atoms=atoms,
                space_group=space_group
            )

        except Exception as e:
            return None

    def _parse_atoms(self, content: str, cell_params: Dict) -> List[Atom]:
        """Parse atom positions and convert to Cartesian coordinates."""
        atoms = []

        lines = content.split('\n')

        in_loop = False
        header_cols = []
        data_start = -1

        for i, line in enumerate(lines):
            line = line.strip()

            if line.startswith('loop_'):
                in_loop = True
                header_cols = []
                continue

            if in_loop:
                if line.startswith('_atom_site_'):
                    header_cols.append(line)
                elif line and not line.startswith('_') and header_cols:
                    if '_atom_site_fract_x' in ' '.join(header_cols):
                        data_start = i
                        break
                    else:
                        in_loop = False
                        header_cols = []

        if data_start < 0 or not header_cols:
            return atoms

        col_map = {}
        for j, col in enumerate(header_cols):
            if 'label' in col:
                col_map['label'] = j
            elif 'type_symbol' in col:
                col_map['symbol'] = j
            elif 'fract_x' in col:
                col_map['x'] = j
            elif 'fract_y' in col:
                col_map['y'] = j
            elif 'fract_z' in col:
                col_map['z'] = j

        if not all(k in col_map for k in ['x', 'y', 'z']):
            return atoms

        frac_to_cart = self._get_frac_to_cart_matrix(cell_params)

        for i in range(data_start, len(lines)):
            line = lines[i].strip()
            if not line or line.startswith('_') or line.startswith('loop_') or line.startswith('#'):
                break

            parts = line.split()
            if len(parts) < max(col_map.values()) + 1:
                continue

            try:
                if 'symbol' in col_map:
                    element = parts[col_map['symbol']]
                elif 'label' in col_map:
                    label = parts[col_map['label']]
                    element = ''.join(c for c in label if c.isalpha())[:2]
                else:
                    continue

                element = element.capitalize()
                if len(element) > 1:
                    element = element[0] + element[1].lower()

                frac_x = self._parse_coord(parts[col_map['x']])
                frac_y = self._parse_coord(parts[col_map['y']])
                frac_z = self._parse_coord(parts[col_map['z']])

                if frac_x is None or frac_y is None or frac_z is None:
                    continue

                frac = np.array([frac_x, frac_y, frac_z])
                cart = frac_to_cart @ frac

                atoms.append(Atom(
                    index=len(atoms),
                    element=element,
                    x=cart[0],
                    y=cart[1],
                    z=cart[2]
                ))

            except (ValueError, IndexError):
                continue

        return atoms

    def _parse_coord(self, s: str) -> Optional[float]:
        """Parse a coordinate value (may include parenthesized uncertainty)."""
        s = re.sub(r'\([^)]*\)', '', s)
        try:
            return float(s)
        except ValueError:
            return None

    def _get_frac_to_cart_matrix(self, cell: Dict) -> np.ndarray:
        """Build the fractional-to-Cartesian transformation matrix."""
        a = cell['a']
        b = cell['b']
        c = cell['c']
        alpha = math.radians(cell['alpha'])
        beta = math.radians(cell['beta'])
        gamma = math.radians(cell['gamma'])

        cos_alpha = math.cos(alpha)
        cos_beta = math.cos(beta)
        cos_gamma = math.cos(gamma)
        sin_gamma = math.sin(gamma)

        volume_factor = math.sqrt(
            1 - cos_alpha**2 - cos_beta**2 - cos_gamma**2 +
            2 * cos_alpha * cos_beta * cos_gamma
        )

        matrix = np.array([
            [a, b * cos_gamma, c * cos_beta],
            [0, b * sin_gamma, c * (cos_alpha - cos_beta * cos_gamma) / sin_gamma],
            [0, 0, c * volume_factor / sin_gamma]
        ])

        return matrix

    def to_raw_molecule(self, cif_data: CIFData) -> RawMolecule:
        """Convert CIFData to a RawMolecule."""
        return RawMolecule(
            mol_id=cif_data.cod_id,
            atoms=cif_data.atoms,
            bond_orders=None,
            properties={'formula': cif_data.formula, 'space_group': cif_data.space_group}
        )

    def iter_molecules(self, limit: int = None) -> Iterator[RawMolecule]:
        """Iterate over molecules parsed from CIF files."""
        self.load(limit=None)

        count = 0
        for filepath in self.files:
            if limit and count >= limit:
                break

            cif_data = self.parse_cif(filepath)
            if cif_data:
                yield self.to_raw_molecule(cif_data)
                count += 1

    def get_statistics(self, limit: int = 1000) -> Dict:
        """Get dataset statistics."""
        self.load()

        stats = {
            'total_files': len(self.files),
            'with_metals': 0,
            'parsed_ok': 0,
            'metal_counts': {},
            'atom_counts': [],
            'formulas': [],
        }

        for i, fp in enumerate(self.files[:limit]):
            cif_data = self.parse_cif(fp)
            if cif_data:
                stats['with_metals'] += 1
                stats['parsed_ok'] += 1
                stats['atom_counts'].append(len(cif_data.atoms))
                stats['formulas'].append(cif_data.formula)

                for atom in cif_data.atoms:
                    if atom.element in TRANSITION_METALS:
                        stats['metal_counts'][atom.element] = stats['metal_counts'].get(atom.element, 0) + 1

        return stats


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--cif-dir", required=True)
    parser.add_argument("--limit", type=int, default=100)
    args = parser.parse_args()

    reader = CIFReader(args.cif_dir)
    stats = reader.get_statistics(args.limit)

    print(f"Total files: {stats['total_files']}")
    print(f"With metals: {stats['with_metals']}")
    print(f"Parsed OK: {stats['parsed_ok']}")
    print(f"Metal distribution: {stats['metal_counts']}")
    print(f"Avg atoms: {np.mean(stats['atom_counts']) if stats['atom_counts'] else 0:.1f}")
