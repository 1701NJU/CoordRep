"""
Continuous Shape Measures (CShM)


v1
- CN=4: Td, SP
- CN=6: Oh, TP (trigonal prism)
"""

import numpy as np
from typing import List, Tuple, Optional
from scipy.spatial.distance import cdist
from scipy.optimize import linear_sum_assignment
from scipy.spatial.transform import Rotation

from ..core import ShapeVector, CoordRepConfig


IDEAL_GEOMETRIES = {
    # CN = 2
    'L-2': np.array([[0, 0, 1], [0, 0, -1]]),  # Linear
    
    # CN = 3
    'TP-3': np.array([  # Trigonal planar
        [1, 0, 0],
        [-0.5, np.sqrt(3)/2, 0],
        [-0.5, -np.sqrt(3)/2, 0]
    ]),
    
    # CN = 4
    'Td': np.array([  # Tetrahedral
        [1, 1, 1],
        [1, -1, -1],
        [-1, 1, -1],
        [-1, -1, 1]
    ]) / np.sqrt(3),
    
    'SP': np.array([  # Square planar
        [1, 0, 0],
        [-1, 0, 0],
        [0, 1, 0],
        [0, -1, 0]
    ]),
    
    # CN = 5
    'TBP': np.array([  # Trigonal bipyramidal
        [0, 0, 1],
        [0, 0, -1],
        [1, 0, 0],
        [-0.5, np.sqrt(3)/2, 0],
        [-0.5, -np.sqrt(3)/2, 0]
    ]),
    
    'SPY': np.array([  # Square pyramidal
        [0, 0, 1],
        [1, 0, 0],
        [-1, 0, 0],
        [0, 1, 0],
        [0, -1, 0]
    ]),
    
    # CN = 6
    'Oh': np.array([  # Octahedral
        [1, 0, 0],
        [-1, 0, 0],
        [0, 1, 0],
        [0, -1, 0],
        [0, 0, 1],
        [0, 0, -1]
    ]),
    
    'TP': np.array([  # Trigonal prismatic
        [1, 0, 0.5],
        [-0.5, np.sqrt(3)/2, 0.5],
        [-0.5, -np.sqrt(3)/2, 0.5],
        [1, 0, -0.5],
        [-0.5, np.sqrt(3)/2, -0.5],
        [-0.5, -np.sqrt(3)/2, -0.5]
    ]),
}

for name in IDEAL_GEOMETRIES:
    template = IDEAL_GEOMETRIES[name]
    template = template - template.mean(axis=0)
    dists = np.linalg.norm(template, axis=1)
    template = template / dists.mean()
    IDEAL_GEOMETRIES[name] = template


def kabsch_rmsd(P: np.ndarray, Q: np.ndarray) -> Tuple[float, np.ndarray]:
    """
 Kabsch RMSD
    
    Args:
 P: (N, 3)
 Q: (N, 3)
    
    Returns:
 rmsd: RMSD
 R:
    """
    P_centered = P - P.mean(axis=0)
    Q_centered = Q - Q.mean(axis=0)
    
    H = P_centered.T @ Q_centered
    
    # SVD
    U, S, Vt = np.linalg.svd(H)
    
    R = Vt.T @ U.T
    
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T
    
    P_rotated = P_centered @ R
    rmsd = np.sqrt(np.mean(np.sum((P_rotated - Q_centered) ** 2, axis=1)))
    
    return rmsd, R


def compute_shape_measure(coords: np.ndarray, 
                         template_name: str) -> float:
    """
 CShM
    
    CShM = 100 * min_permutation(RMSD^2) / <r^2>
    
    Args:
 coords: (N, 3)
 template_name:
    
    Returns:
 CShM
    """
    template = IDEAL_GEOMETRIES.get(template_name)
    if template is None:
        return float('inf')
    
    if len(coords) != len(template):
        return float('inf')
    
    coords_centered = coords - coords.mean(axis=0)
    r_mean = np.linalg.norm(coords_centered, axis=1).mean()
    if r_mean < 1e-6:
        return float('inf')
    coords_norm = coords_centered / r_mean
    
    cost_matrix = cdist(coords_norm, template)
    row_ind, col_ind = linear_sum_assignment(cost_matrix)
    
    template_permuted = template[col_ind]
    
    # Kabsch RMSD
    rmsd, _ = kabsch_rmsd(coords_norm, template_permuted)
    
    # CShM = 100 * RMSD^2
    cshm = 100 * rmsd ** 2
    
    return cshm


class ShapeCalculator:
 """"""
    
    def __init__(self, config: CoordRepConfig = None):
        self.config = config or CoordRepConfig.default()
    
    def get_reference_shapes(self, cn: int) -> List[str]:
 """ CN """
        if cn == 2:
            return ['L-2']
        elif cn == 3:
            return ['TP-3']
        elif cn == 4:
            return self.config.cn4_ref_shapes
        elif cn == 5:
            return ['TBP', 'SPY']
        elif cn == 6:
            return self.config.cn6_ref_shapes
        else:
            return []
    
    def compute(self, coords: np.ndarray, cn: int) -> ShapeVector:
        """
        
        Args:
 coords:
 cn:
        
        Returns:
 ShapeVector
        """
        ref_shapes = self.get_reference_shapes(cn)
        
        if not ref_shapes:
            return ShapeVector(
                cn=cn,
                ref_shapes=[],
                values=np.array([]),
                values_rounded=np.array([])
            )
        
        values = []
        for shape_name in ref_shapes:
            cshm = compute_shape_measure(coords, shape_name)
            values.append(cshm)
        
        values = np.array(values)
        values_rounded = np.round(values, self.config.rounding_decimals)
        
        return ShapeVector(
            cn=cn,
            ref_shapes=ref_shapes,
            values=values,
            values_rounded=values_rounded
        )
    
    def get_dominant_geometry(self, shape_vector: ShapeVector) -> Optional[str]:
 """CShM """
        if len(shape_vector.ref_shapes) == 0:
            return None
        
        min_idx = np.argmin(shape_vector.values)
        return shape_vector.ref_shapes[min_idx]





