"""
Shape Vector (CShM) - Sprint 2 Assignment-Invariant Version

1.
2. Hungarian assignment tie-break
3. donor

CShM
CShM(P, Q) = 100 * min_permutation(Σ|D_P[i,j] - D_Q[σ(i),σ(j)]|²) / Σ|D_P[i,j]|²
"""

import numpy as np
from typing import List, Tuple, Optional, Dict
from dataclasses import dataclass, field
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist


def _normalize_template(coords: np.ndarray) -> np.ndarray:
 """Normalize + 1"""
    coords = coords - coords.mean(axis=0)
    avg_dist = np.linalg.norm(coords, axis=1).mean()
    if avg_dist > 1e-6:
        coords = coords / avg_dist
    return coords


# CN=2
_LINEAR_2 = _normalize_template(np.array([
    [0, 0, 1],
    [0, 0, -1]
]))

# CN=3
_TRIGONAL_PLANAR_3 = _normalize_template(np.array([
    [1, 0, 0],
    [-0.5, np.sqrt(3)/2, 0],
    [-0.5, -np.sqrt(3)/2, 0]
]))

# CN=4
_TETRAHEDRAL_4 = _normalize_template(np.array([
    [1, 1, 1],
    [1, -1, -1],
    [-1, 1, -1],
    [-1, -1, 1]
]))

_SQUARE_PLANAR_4 = _normalize_template(np.array([
    [1, 0, 0],
    [-1, 0, 0],
    [0, 1, 0],
    [0, -1, 0]
]))

# CN=5
_TRIGONAL_BIPYRAMIDAL_5 = _normalize_template(np.array([
    [0, 0, 1],
    [0, 0, -1],
    [1, 0, 0],
    [-0.5, np.sqrt(3)/2, 0],
    [-0.5, -np.sqrt(3)/2, 0]
]))

_SQUARE_PYRAMIDAL_5 = _normalize_template(np.array([
    [0, 0, 1],
    [1, 0, 0],
    [-1, 0, 0],
    [0, 1, 0],
    [0, -1, 0]
]))

# CN=6
_OCTAHEDRAL_6 = _normalize_template(np.array([
    [1, 0, 0],
    [-1, 0, 0],
    [0, 1, 0],
    [0, -1, 0],
    [0, 0, 1],
    [0, 0, -1]
]))

_TRIGONAL_PRISMATIC_6 = _normalize_template(np.array([
    [1, 0, 0.5],
    [-0.5, np.sqrt(3)/2, 0.5],
    [-0.5, -np.sqrt(3)/2, 0.5],
    [1, 0, -0.5],
    [-0.5, np.sqrt(3)/2, -0.5],
    [-0.5, -np.sqrt(3)/2, -0.5]
]))


IDEAL_GEOMETRIES = {
    2: {'L': _LINEAR_2},
    3: {'TP': _TRIGONAL_PLANAR_3},
    4: {'Td': _TETRAHEDRAL_4, 'SP': _SQUARE_PLANAR_4},
    5: {'TBP': _TRIGONAL_BIPYRAMIDAL_5, 'SPY': _SQUARE_PYRAMIDAL_5},
    6: {'Oh': _OCTAHEDRAL_6, 'TPr': _TRIGONAL_PRISMATIC_6},
}


def compute_atom_signature(sorted_distances: np.ndarray) -> int:
    """
 tie-break
    
    """
    discretized = np.round(sorted_distances * 10000).astype(np.int64)
    
    signature = 0
 for i, d in enumerate(discretized[:6]): # 6
 signature = signature * 100003 + int(d) #
    
    return signature % (10**15)


def compute_cshm_with_tiebreak(P: np.ndarray, Q_template: np.ndarray, 
                               atom_signatures: np.ndarray) -> Tuple[float, np.ndarray]:
    """
 CShM - tie-break
    
 Sprint 2
 1.
 2.
 3. assignment
    
 assignment
    """
    N = len(P)
    if N != len(Q_template):
        return float('inf'), np.arange(N)
    
    D_P = cdist(P, P, metric='euclidean')
    D_Q = cdist(Q_template, Q_template, metric='euclidean')
    
    
    
    P_features = np.array([np.sort(D_P[i]) for i in range(N)])
    Q_features = np.array([np.sort(D_Q[j]) for j in range(N)])
    
    sig_order = np.argsort(atom_signatures)
    
    assigned = set()
    assignment = np.zeros(N, dtype=int)
    
    for i in sig_order:
        best_j = -1
        best_cost = float('inf')
        
        for j in range(N):
            if j in assigned:
                continue
            cost = np.sum((P_features[i] - Q_features[j]) ** 2)
 cost += 1e-10 * j # j
            if cost < best_cost:
                best_cost = cost
                best_j = j
        
        assignment[i] = best_j
        assigned.add(best_j)
    
    D_Q_permuted = D_Q[assignment][:, assignment]
    numerator = np.sum((D_P - D_Q_permuted) ** 2)
    denominator = np.sum(D_P ** 2)
    
    if denominator < 1e-10:
        return float('inf'), assignment
    
    cshm = 100 * numerator / denominator
    
    return cshm, assignment


@dataclass
class ShapeResult:
 """Compute results"""
    cn: int
    ref_shapes: List[str]
 values: np.ndarray # CShM
 values_rounded: np.ndarray #
 assignments: Dict[str, np.ndarray] #
    
    best_shape: str = ""
    best_cshm: float = 0.0
    delta_to_second: float = 0.0
    shape_class: str = ""  # "ideal" / "distorted" / "ambiguous"
    
    def to_canonical_string(self) -> str:
 """"""
        if not self.ref_shapes:
            return ""
        
        if self.delta_to_second < 0.5:
            delta_bin = 0
        elif self.delta_to_second < 2:
            delta_bin = 1
        elif self.delta_to_second < 5:
            delta_bin = 2
        else:
            delta_bin = 3
        
        return f"<ShapeBest:{self.best_shape}|Delta:{delta_bin}>"
    
    def to_dict(self) -> dict:
        return {
            'cn': self.cn,
            'ref_shapes': self.ref_shapes,
            'values_rounded': self.values_rounded.tolist(),
            'best_shape': self.best_shape,
            'best_cshm': round(self.best_cshm, 2),
            'delta_to_second': round(self.delta_to_second, 2),
        }


class ShapeCalculator:
 """ (Assignment-Invariant, Rotation-Invariant)"""
    
    def __init__(self, rounding_decimals: int = 2):
        self.rounding_decimals = rounding_decimals
    
    def _compute_atom_signatures(self, donor_coords: np.ndarray, metal_coord: np.ndarray) -> np.ndarray:
        """
 donor
        
 + donors
        """
        n = len(donor_coords)
        signatures = np.zeros(n, dtype=np.int64)
        
        distances_to_metal = np.linalg.norm(donor_coords - metal_coord, axis=1)
        
        D = cdist(donor_coords, donor_coords, metric='euclidean')
        
        for i in range(n):
            feature = np.concatenate([
                [distances_to_metal[i]],
                np.sort(D[i])
            ])
            signatures[i] = compute_atom_signature(feature)
        
        return signatures
    
    def _canonicalize_donor_order(self, donor_coords: np.ndarray, 
                                  metal_coord: np.ndarray,
                                  signatures: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """
 donor
        
        """
        n = len(donor_coords)
        if n <= 1:
            return donor_coords, signatures
        
        sorted_indices = np.argsort(signatures)
        
        return donor_coords[sorted_indices], signatures[sorted_indices]
    
    def compute(self, donor_coords: np.ndarray, metal_coord: np.ndarray) -> ShapeResult:
        """
        """
        cn = len(donor_coords)
        
        if cn not in IDEAL_GEOMETRIES:
            return ShapeResult(
                cn=cn,
                ref_shapes=[],
                values=np.array([]),
                values_rounded=np.array([]),
                assignments={},
                best_shape=f"CN{cn}",
                best_cshm=0.0,
                delta_to_second=0.0,
                shape_class="unknown"
            )
        
        templates = IDEAL_GEOMETRIES[cn]
        ref_shapes = list(templates.keys())
        
        signatures = self._compute_atom_signatures(donor_coords, metal_coord)
        
        donor_coords_sorted, signatures_sorted = self._canonicalize_donor_order(
            donor_coords, metal_coord, signatures
        )
        
        P = donor_coords_sorted - metal_coord
        
        avg_dist = np.linalg.norm(P, axis=1).mean()
        if avg_dist > 1e-6:
            P = P / avg_dist
        
        values = []
        assignments = {}
        
        for shape_name in ref_shapes:
            Q = templates[shape_name]
            cshm, assignment = compute_cshm_with_tiebreak(P, Q, signatures_sorted)
            values.append(cshm)
            assignments[shape_name] = assignment
        
        values = np.array(values)
        values_rounded = np.round(values, self.rounding_decimals)
        
        sorted_indices = np.argsort(values)
        best_idx = sorted_indices[0]
        best_shape = ref_shapes[best_idx]
        best_cshm = values[best_idx]
        
        if len(values) > 1:
            second_idx = sorted_indices[1]
            delta_to_second = values[second_idx] - values[best_idx]
        else:
            delta_to_second = float('inf')
        
        if best_cshm < 1.0:
            shape_class = "ideal"
        elif best_cshm < 5.0:
            shape_class = "distorted"
        else:
            shape_class = "irregular"
        
        if delta_to_second < 1.0:
            shape_class = "ambiguous"
        
        return ShapeResult(
            cn=cn,
            ref_shapes=ref_shapes,
            values=values,
            values_rounded=values_rounded,
            assignments=assignments,
            best_shape=best_shape,
            best_cshm=best_cshm,
            delta_to_second=delta_to_second,
            shape_class=shape_class
        )
