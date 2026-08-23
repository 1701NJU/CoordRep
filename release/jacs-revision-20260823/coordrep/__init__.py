"""
coordrep -- Canonical Representation for Coordination Complexes

A canonical, rotation-invariant, and permutation-invariant string
representation for transition-metal coordination complexes, with
continuous geometry encoding and modular ligand decomposition.

Usage:
    from coordrep import encode

    cc = encode(xyz_path="complex.xyz", bo_path="complex.BO")
    s = cc.canonicalize().to_string()
    feat = cc.to_features()
"""

from .core import (
    CoordComplex,
    MetalState,
    ShapeVector,
    LigandModule,
    ConstraintSet,
    DonorSite,
    AssemblyGraph,
    CoordRepConfig,
    CoordRepIssue,
    IssueCodes,
)

from .encode import encode, encode_molecule, batch_encode

__version__ = "1.1.2rc3"
__all__ = [
    'encode',
    'encode_molecule',
    'batch_encode',
    'CoordComplex',
    'MetalState',
    'ShapeVector',
    'LigandModule',
    'ConstraintSet',
    'DonorSite',
    'CoordRepConfig',
]
