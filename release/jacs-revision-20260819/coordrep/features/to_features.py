"""
Feature Extraction

ML
"""

import numpy as np
from typing import Dict, List, Any

from ..core import CoordComplex


def extract_features(cc: CoordComplex) -> Dict[str, Any]:
    """
    ML

    Returns:
        {
        "tokens": [...], # token
 "geom": np.array(...), #
 "metal": np.array(...), #
 "ligands": [...], #
        }
    """
    features = {}

    features['tokens'] = tokenize_complex(cc)

    features['geom'] = extract_geom_features(cc)

    features['metal'] = extract_metal_features(cc)

    features['ligands'] = extract_ligand_features(cc)

    return features


def tokenize_complex(cc: CoordComplex) -> List[str]:
    """ CoordComplex token """
    tokens = []

    # Metal tokens
    tokens.append(f"[{cc.metal.element}]")
    if cc.metal.oxidation is not None:
        tokens.append(f"ox={cc.metal.oxidation}")
    if cc.metal.dcount is not None:
        tokens.append(f"d{cc.metal.dcount}")

    # CN token
    cn = len(cc.graph.donor_indices) if cc.graph else 0
    tokens.append(f"CN={cn}")

    # Shape tokens
    if cc.shape is not None and len(cc.shape.ref_shapes) > 0:
        min_idx = np.argmin(cc.shape.values)
        dominant = cc.shape.ref_shapes[min_idx]
        tokens.append(f"GEO={dominant}")

        for name, val in zip(cc.shape.ref_shapes, cc.shape.values):
            if val < 1.0:
                tokens.append(f"{name}_low")
            elif val < 5.0:
                tokens.append(f"{name}_mid")
            else:
                tokens.append(f"{name}_high")

    # Constraint tokens
    for s1, s2 in cc.constraints.trans_pairs:
        tokens.append(f"trans:{s1.lig_id}-{s2.lig_id}")

    if 'fac_mer' in cc.constraints.notes:
        tokens.append(f"fm={cc.constraints.notes['fac_mer']}")

    # Ligand tokens
    for lig in cc.ligands:
        tokens.append(f"${lig.lig_id}")
        tokens.append(f"dent={lig.dent}")
        if lig.eta:
            tokens.append(f"eta={lig.eta}")

    return tokens


def extract_geom_features(cc: CoordComplex) -> np.ndarray:
    """"""
    features = []

    # CN (one-hot, 1-9)
    cn = len(cc.graph.donor_indices) if cc.graph else 0
    cn_onehot = [0] * 9
    if 1 <= cn <= 9:
        cn_onehot[cn - 1] = 1
    features.extend(cn_onehot)

    # Shape values (padded to 4)
    shape_vals = [0.0] * 4
    if cc.shape is not None:
        for i, val in enumerate(cc.shape.values[:4]):
            shape_vals[i] = val
    features.extend(shape_vals)

    # Trans pair count
    features.append(len(cc.constraints.trans_pairs))

    # fac/mer (one-hot: none, fac, mer)
    fm = cc.constraints.notes.get('fac_mer', 'none')
    features.extend([
        1 if fm == 'none' else 0,
        1 if fm == 'fac' else 0,
        1 if fm == 'mer' else 0,
    ])

    return np.array(features, dtype=np.float32)


def extract_metal_features(cc: CoordComplex) -> np.ndarray:
    """"""
    features = []

    m = cc.metal

    # Row (one-hot: 4, 5, 6)
    row = m.phys_vec[0] if m.phys_vec is not None else 0
    features.extend([
        1 if row == 4 else 0,
        1 if row == 5 else 0,
        1 if row == 6 else 0,
    ])

    # Group (normalized, 3-12)
    group = m.phys_vec[1] if m.phys_vec is not None else 0
    features.append((group - 3) / 9)

    # Electronegativity
    en = m.phys_vec[2] if m.phys_vec is not None else 0
    features.append(en / 3)  # normalized

    # Oxidation state (normalized, -2 to +6)
    ox = m.oxidation if m.oxidation is not None else 2
    features.append((ox + 2) / 8)

    # d-count (normalized, 0-10)
    dc = m.dcount if m.dcount is not None else 6
    features.append(dc / 10)

    return np.array(features, dtype=np.float32)


def extract_ligand_features(cc: CoordComplex) -> List[Dict]:
    """"""
    lig_features = []

    for lig in cc.ligands:
        feat = {
            'lig_id': lig.lig_id,
            'smiles': lig.smiles,
            'dent': lig.dent,
            'eta': lig.eta,
            'donor_elements': lig.donor_elements,
            'mw': lig.meta.get('mw', 0),
        }
        lig_features.append(feat)

    return lig_features
