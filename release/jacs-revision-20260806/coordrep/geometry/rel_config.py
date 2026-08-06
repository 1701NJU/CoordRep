"""
Relative Configuration - Sprint 2 Version

1. trans/cis : TRANS / CIS / OTHER
2. hysteresis
3.
"""

import numpy as np
from typing import List, Tuple, Dict, Optional
from collections import defaultdict
from enum import Enum

from ..core import (
    ConstraintSet, DonorSite, LigandModule,
    CoordRepConfig, CoordRepIssue
)
from ..io.tmqm_reader import RawMolecule


class AngleClass(Enum):
    TRANS = "trans"
    CIS = "cis"
    OTHER = "other"


def classify_angle(angle: float,
                  trans_range: Tuple[float, float] = (165.0, 195.0),
                  cis_range: Tuple[float, float] = (75.0, 105.0)) -> AngleClass:
    """
    -

    TRANS: > 165°
    CIS: 75° - 105°
    OTHER:

     155-165° OTHER
    """
    if angle >= trans_range[0]:
        return AngleClass.TRANS
    elif cis_range[0] <= angle <= cis_range[1]:
        return AngleClass.CIS
    else:
        return AngleClass.OTHER


class ConstraintDetector:
    """ - Sprint 2 """

    def __init__(self, config: CoordRepConfig = None):
        self.config = config or CoordRepConfig.default()

        self.trans_min = 165.0 #
        self.cis_min = 75.0
        self.cis_max = 105.0

    def detect_constraints(self,
                          mol: RawMolecule,
                          metal_idx: int,
                          donor_indices: List[int],
                          ligands: List[LigandModule]) -> ConstraintSet:
        """Detect constraints."""
        donor_sites = self._build_simple_donor_sites(donor_indices, ligands, mol, metal_idx)
        return self.detect_constraints_with_sites(
            mol, metal_idx, donor_indices, ligands, donor_sites
        )

    def detect_constraints_with_sites(self,
                                     mol: RawMolecule,
                                     metal_idx: int,
                                     donor_indices: List[int],
                                     ligands: List[LigandModule],
                                     donor_sites: Dict[int, DonorSite]) -> ConstraintSet:
        """
        donor sites

        TRANS CIS
        """
        trans_pairs = []
        cis_pairs = []

        metal_coord = mol.atoms[metal_idx].coords

        n_donors = len(donor_indices)

        for i in range(n_donors):
            for j in range(i + 1, n_donors):
                d1_idx = donor_indices[i]
                d2_idx = donor_indices[j]

                angle = self._compute_angle(
                    mol.atoms[d1_idx].coords,
                    metal_coord,
                    mol.atoms[d2_idx].coords
                )

                site1 = donor_sites.get(d1_idx)
                site2 = donor_sites.get(d2_idx)

                if site1 is None or site2 is None:
                    continue

                angle_class = classify_angle(
                    angle,
                    trans_range=(self.trans_min, 195.0),
                    cis_range=(self.cis_min, self.cis_max)
                )

                if angle_class == AngleClass.TRANS:
                    trans_pairs.append((site1, site2))
                elif angle_class == AngleClass.CIS:
                    cis_pairs.append((site1, site2))

        constraints = ConstraintSet(
            trans_pairs=trans_pairs,
            cis_pairs=cis_pairs,
            notes={}
        )

        fac_mer = self._detect_fac_mer(trans_pairs, cis_pairs, ligands)
        if fac_mer:
            constraints.notes['fac_mer'] = fac_mer

        return constraints.canonicalize()

    def _build_simple_donor_sites(self,
                                  donor_indices: List[int],
                                  ligands: List[LigandModule],
                                  mol: RawMolecule,
                                  metal_idx: int) -> Dict[int, DonorSite]:
        """Build simple donor sites."""
        from ..graph.donor_sites import assign_donor_ranks

        result = {}
        for lig in ligands:
            valid_attach = [d for d in lig.attach_atoms if d in donor_indices]
            if valid_attach:
                temp_lig = LigandModule(
                    lig_id=lig.lig_id,
                    smiles=lig.smiles,
                    attach_atoms=valid_attach,
                    donor_elements=[mol.atoms[d].element for d in valid_attach],
                    dent=len(valid_attach),
                    eta=lig.eta,
                    charge=lig.charge,
                    meta=lig.meta
                )
                sites = assign_donor_ranks(mol, temp_lig, metal_idx)
                result.update(sites)

        return result

    def _compute_angle(self, p1: np.ndarray, vertex: np.ndarray, p2: np.ndarray) -> float:
        """ p1-vertex-p2 """
        v1 = p1 - vertex
        v2 = p2 - vertex

        v1_norm = np.linalg.norm(v1)
        v2_norm = np.linalg.norm(v2)

        if v1_norm < 1e-6 or v2_norm < 1e-6:
            return 0.0

        cos_angle = np.dot(v1, v2) / (v1_norm * v2_norm)
        cos_angle = np.clip(cos_angle, -1.0, 1.0)

        return np.degrees(np.arccos(cos_angle))

    def _detect_fac_mer(self,
                       trans_pairs: List[Tuple[DonorSite, DonorSite]],
                       cis_pairs: List[Tuple[DonorSite, DonorSite]],
                       ligands: List[LigandModule]) -> Optional[str]:
        """Detect fac/mer isomerism."""
        lig_trans = defaultdict(int)
        lig_cis = defaultdict(int)

        for s1, s2 in trans_pairs:
            if s1.lig_id == s2.lig_id:
                lig_trans[s1.lig_id] += 1

        for s1, s2 in cis_pairs:
            if s1.lig_id == s2.lig_id:
                lig_cis[s1.lig_id] += 1

        for lig in ligands:
            if lig.dent == 3:
                if lig_trans[lig.lig_id] > 0:
                    return 'mer'
                if lig_cis[lig.lig_id] == 3:
                    return 'fac'

        return None
