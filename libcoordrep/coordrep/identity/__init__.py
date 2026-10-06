"""
CoordRep-ID: Hierarchical Identity Key System
==============================================

Projects full CoordRep strings onto multiple chemical-equivalence levels:

    L0  StateKey   – full canonical CoordRep string (geometry-snapshot level)
    L1  ShapeID    – binned CShM (best-shape + S1_bin + delta_bin), no exact decimals
    L2  TopoID     – metal + best-shape label + sorted ligand SMILES
    L3  ConnID     – metal + CN + sorted ligand SMILES (geometry-free)
"""

from .identity_keys import (
    IdentityKeys,
    extract_identity_keys,
    extract_identity_keys_from_cc,
)
from .shape_binning import (
    ShapeBin,
    BinnedShape,
    bin_shape,
    bin_shape_from_result,
    DEFAULT_SHAPE_BIN_CONFIG,
)
