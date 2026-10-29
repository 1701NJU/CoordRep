"""Order-independent quotient-site grouping for CSD P1 periodic molecules.

The frozen MOF census grouped periodic images sequentially and retained the
first encountered wrapped fractional coordinate as the group representative.
That is adequate for reconstructing an edge witness, but it is not canonical
when symmetry-expanded coordinates differ slightly because of CIF rounding.

Here, same-element images are clustered by an order-independent union-find on
their minimum-image Cartesian separation.  Each quotient site's coordinate is
then the arithmetic mean in a unique local lift of the torus.  For a compact
cluster this mean is equivariant to an origin shift and to an exact unimodular
change of lattice basis, while eliminating atom-row encounter order.
"""

from __future__ import annotations

import itertools
import math
from collections import defaultdict
from typing import Any, Dict, Iterable, List, Mapping, MutableMapping, Sequence, Tuple

import numpy as np


DEFAULT_IMAGE_CLUSTER_TOLERANCE_ANGSTROM = 0.002


class _UnionFind:
    def __init__(self, size: int):
        self.parent = list(range(size))
        self.rank = [0] * size

    def find(self, value: int) -> int:
        parent = self.parent[value]
        if parent != value:
            self.parent[value] = self.find(parent)
        return self.parent[value]

    def union(self, left: int, right: int) -> None:
        left_root = self.find(left)
        right_root = self.find(right)
        if left_root == right_root:
            return
        if self.rank[left_root] < self.rank[right_root]:
            left_root, right_root = right_root, left_root
        self.parent[right_root] = left_root
        if self.rank[left_root] == self.rank[right_root]:
            self.rank[left_root] += 1


def _fractional(atom: Any) -> np.ndarray:
    coordinates = atom.fractional_coordinates
    if coordinates is None:
        raise ValueError(f"Atom {atom.index} lacks fractional coordinates")
    value = np.array([coordinates.x, coordinates.y, coordinates.z], dtype=float)
    if not np.all(np.isfinite(value)):
        raise ValueError(f"Atom {atom.index} has nonfinite fractional coordinates")
    return value


def _cartesian(atom: Any) -> np.ndarray:
    coordinates = atom.coordinates
    if coordinates is None:
        raise ValueError(f"Atom {atom.index} lacks Cartesian coordinates")
    value = np.array([coordinates.x, coordinates.y, coordinates.z], dtype=float)
    if not np.all(np.isfinite(value)):
        raise ValueError(f"Atom {atom.index} has nonfinite Cartesian coordinates")
    return value


_NEIGHBOR_OFFSETS = tuple(itertools.product((-1, 0, 1), repeat=3))


def _minimum_image(
    fractional_difference: Sequence[float], cell_matrix: np.ndarray
) -> Tuple[np.ndarray, np.ndarray, float]:
    """Return residual, subtracted lattice integer, and Cartesian norm."""

    difference = np.asarray(fractional_difference, dtype=float)
    centre = np.rint(difference).astype(int)
    best_residual = None
    best_integer = None
    best_norm = float("inf")
    for offset in _NEIGHBOR_OFFSETS:
        integer = centre + np.asarray(offset, dtype=int)
        residual = difference - integer
        norm = float(np.linalg.norm(residual @ cell_matrix))
        candidate_key = (norm, tuple(int(value) for value in integer))
        best_key = (
            best_norm,
            tuple(int(value) for value in best_integer)
            if best_integer is not None
            else (0, 0, 0),
        )
        if best_residual is None or candidate_key < best_key:
            best_residual = residual
            best_integer = integer
            best_norm = norm
    assert best_residual is not None and best_integer is not None
    return best_residual, best_integer, best_norm


def canonical_periodic_groups(
    molecule: Any,
    cell_matrix: Sequence[Sequence[float]],
    *,
    tolerance_angstrom: float = DEFAULT_IMAGE_CLUSTER_TOLERANCE_ANGSTROM,
) -> Tuple[
    Dict[int, Dict[str, Any]],
    Dict[int, int],
    Dict[int, Dict[str, Any]],
]:
    """Group periodic atom images without atom encounter-order dependence."""

    matrix = np.asarray(cell_matrix, dtype=float)
    if matrix.shape != (3, 3) or not np.all(np.isfinite(matrix)):
        raise ValueError("cell_matrix must be finite and 3 by 3")
    singular_values = np.linalg.svd(matrix, compute_uv=False)
    sigma_min = float(np.min(singular_values))
    if sigma_min <= 0.0 or not math.isfinite(sigma_min):
        raise ValueError("cell_matrix is singular")
    tolerance = float(tolerance_angstrom)
    if tolerance <= 0.0 or not math.isfinite(tolerance):
        raise ValueError("tolerance_angstrom must be positive and finite")

    records = []
    for atom in molecule.atoms:
        fractional = _fractional(atom)
        records.append(
            {
                "atom": atom,
                "atom_index": int(atom.index),
                "element": str(atom.atomic_symbol),
                "fractional": fractional,
                "wrapped": np.mod(fractional, 1.0),
                "cartesian": _cartesian(atom),
                "label": str(atom.label),
            }
        )
    union = _UnionFind(len(records))

    # A Cartesian tolerance implies this conservative component-wise
    # fractional bound.  Binning only prunes comparisons; the final decision
    # always uses the minimum-image Cartesian norm.
    fractional_bound = min(0.49, tolerance / sigma_min)
    scale = max(1, int(math.floor(1.0 / max(fractional_bound, 1.0e-12))))
    bins: MutableMapping[Tuple[Any, ...], List[int]] = defaultdict(list)
    ordered_indices = sorted(
        range(len(records)),
        key=lambda index: (
            records[index]["element"],
            tuple(float(value) for value in records[index]["wrapped"]),
            records[index]["label"],
        ),
    )
    for index in ordered_indices:
        record = records[index]
        wrapped = record["wrapped"]
        bin_index = tuple(
            int(math.floor(float(value) * scale)) % scale for value in wrapped
        )
        candidates = set()
        for offset in _NEIGHBOR_OFFSETS:
            key = (
                record["element"],
                (bin_index[0] + offset[0]) % scale,
                (bin_index[1] + offset[1]) % scale,
                (bin_index[2] + offset[2]) % scale,
            )
            candidates.update(bins.get(key, ()))
        for candidate in candidates:
            if records[candidate]["element"] != record["element"]:
                continue
            _, _, distance = _minimum_image(
                record["fractional"] - records[candidate]["fractional"],
                matrix,
            )
            if distance <= tolerance:
                union.union(index, candidate)
        bins[
            (
                record["element"],
                bin_index[0],
                bin_index[1],
                bin_index[2],
            )
        ].append(index)

    components: MutableMapping[int, List[int]] = defaultdict(list)
    for index in range(len(records)):
        components[union.find(index)].append(index)

    component_payloads = []
    for member_indices in components.values():
        ordered_members = sorted(
            member_indices,
            key=lambda index: (
                tuple(float(value) for value in records[index]["wrapped"]),
                records[index]["label"],
            ),
        )
        anchor = np.asarray(records[ordered_members[0]]["fractional"], dtype=float)
        lifted = []
        for index in ordered_members:
            residual, _, _ = _minimum_image(
                records[index]["fractional"] - anchor, matrix
            )
            lifted.append(anchor + residual)
        mean_lift = np.mean(np.asarray(lifted, dtype=float), axis=0)
        wrapped_mean = np.mod(mean_lift, 1.0)
        maximum_residual = 0.0
        translations: Dict[int, Tuple[int, int, int]] = {}
        for index in ordered_members:
            residual, integer, distance = _minimum_image(
                records[index]["fractional"] - wrapped_mean, matrix
            )
            if distance > tolerance + 1.0e-12:
                raise ValueError(
                    "A quotient component exceeds the canonical lift tolerance"
                )
            maximum_residual = max(maximum_residual, distance)
            translations[index] = tuple(int(value) for value in integer)
        component_payloads.append(
            {
                "element": records[ordered_members[0]]["element"],
                "wrapped_fractional": wrapped_mean,
                "member_indices": ordered_members,
                "translations": translations,
                "maximum_image_residual_angstrom": maximum_residual,
            }
        )

    # IDs are deterministic for a fixed origin/basis, but downstream identity
    # uses only the induced equality partition, never these integers.
    component_payloads.sort(
        key=lambda component: (
            component["element"],
            tuple(float(value) for value in component["wrapped_fractional"]),
        )
    )
    atom_info: Dict[int, Dict[str, Any]] = {}
    atom_to_group: Dict[int, int] = {}
    groups: Dict[int, Dict[str, Any]] = {}
    for group_id, component in enumerate(component_payloads):
        groups[group_id] = {
            "group_id": group_id,
            "element": component["element"],
            "wrapped_fractional": component["wrapped_fractional"],
            "atom_indices": [],
            "labels": [],
            "maximum_image_residual_angstrom": component[
                "maximum_image_residual_angstrom"
            ],
        }
        for record_index in component["member_indices"]:
            record = records[record_index]
            atom_index = record["atom_index"]
            translation = component["translations"][record_index]
            atom_info[atom_index] = {
                "atom_index": atom_index,
                "element": record["element"],
                "fractional": record["fractional"],
                "wrapped_fractional": component["wrapped_fractional"],
                "translation": translation,
                "cartesian": record["cartesian"],
                "label": record["label"],
            }
            atom_to_group[atom_index] = group_id
            groups[group_id]["atom_indices"].append(atom_index)
            groups[group_id]["labels"].append(record["label"])
    return atom_info, atom_to_group, groups
