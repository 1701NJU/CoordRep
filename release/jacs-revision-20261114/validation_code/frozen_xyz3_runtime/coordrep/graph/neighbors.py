"""XYZ-v2 candidate perception, not a claim of general chemical bond truth.

Valid supplied BO: threshold only; no coordinate-based gap filling.
Coordinates only: pairwise covalent-radius cap and simultaneous heavy-atom
Gabriel empty-ball screening. Every excluded candidate has an explicit reason.
"""
from dataclasses import dataclass
from typing import List, Optional
import numpy as np
from rdkit import Chem
from ..io.tmqm_reader import RawMolecule

XYZ_ADAPTER_VERSION = "xyz2.20261113"
COVALENT_RADIUS_FACTOR = 1.3
RADIUS_BOUNDARY_TOLERANCE_A = 1e-12
GABRIEL_MARGIN_TOLERANCE_A2 = 1e-12
BO_SYMMETRY_TOLERANCE = 1e-8
DONOR_ELEMENTS = {"N", "O", "S", "P", "Cl", "Br", "I", "F",
                  "C", "Se", "Te", "As", "Si"}
_PT = Chem.GetPeriodicTable()
ATOMIC_NUMBERS = {s: _PT.GetAtomicNumber(s) for s in DONOR_ELEMENTS | {"H"}}


@dataclass
class NeighborInfo:
    atom_idx: int
    element: str
    distance: float
    bond_order: float
    score: float

    def get_sort_key(self):
        return (-self.score, self.distance, ATOMIC_NUMBERS.get(self.element, 0) * -1, self.atom_idx)


def covalent_radius(element):
    radius = float(_PT.GetRcovalent(element))
    if not np.isfinite(radius) or radius <= 0:
        raise ValueError("No positive covalent radius for " + element)
    return radius


def validated_bond_orders(mol):
    if mol.bond_orders is None:
        return None
    bo = np.asarray(mol.bond_orders, dtype=float)
    n = len(mol.atoms)
    if bo.shape != (n, n):
        raise ValueError("BO matrix must have one row and column per atom")
    if not np.all(np.isfinite(bo)):
        raise ValueError("BO matrix contains nonfinite values")
    if not np.allclose(bo, bo.T, rtol=0, atol=BO_SYMMETRY_TOLERANCE):
        raise ValueError("BO matrix is not symmetric within fixed tolerance")
    if np.any(bo < -BO_SYMMETRY_TOLERANCE):
        raise ValueError("BO matrix contains negative values")
    return bo


def perceive_neighbors(mol: RawMolecule, metal_idx: int, bo_threshold=0.3,
                       distance_cutoff: Optional[float] = None):
    n = len(mol.atoms)
    if not 0 <= metal_idx < n:
        raise IndexError("metal_idx outside atom list")
    for atom in mol.atoms:
        try:
            number = _PT.GetAtomicNumber(atom.element)
        except Exception as exc:
            raise ValueError("Unknown XYZ element: " + str(atom.element)) from exc
        if number <= 0:
            raise ValueError("XYZ requires real chemical elements: " + str(atom.element))
    coords = np.asarray([a.coords for a in mol.atoms], dtype=float)
    if coords.shape != (n, 3) or not np.all(np.isfinite(coords)):
        raise ValueError("All atoms require finite Cartesian coordinates")
    bo = validated_bond_orders(mol)
    center = coords[metal_idx]
    neighbors, decisions = [], []
    heavy = [j for j, a in enumerate(mol.atoms) if j != metal_idx and a.element != "H"]
    metal_radius = covalent_radius(mol.atoms[metal_idx].element)
    for i, atom in enumerate(mol.atoms):
        if i == metal_idx:
            continue
        distance = float(np.linalg.norm(coords[i] - center))
        decision = {"atom_index0": i, "element": atom.element,
                    "distance_A": distance, "accepted": False,
                    "blockers": [], "method": "SUPPLIED_BO_ONLY" if bo is not None else "RADII_GABRIEL_GEOMETRIC_CANDIDATE"}
        if atom.element not in DONOR_ELEMENTS:
            decision["reason"] = "H_EXCLUDED" if atom.element == "H" else "ELEMENT_OUTSIDE_SUPPORTED_DONOR_DOMAIN"
        elif bo is not None:
            order = float(bo[metal_idx, i])
            decision.update(bond_order=order, accepted=order > bo_threshold,
                            reason="BO_ABOVE_THRESHOLD" if order > bo_threshold else "BO_NOT_ABOVE_THRESHOLD_NO_DISTANCE_OVERRIDE")
        else:
            cap = COVALENT_RADIUS_FACTOR * (metal_radius + covalent_radius(atom.element))
            if distance_cutoff is not None:
                cap = min(cap, float(distance_cutoff))
            decision["radius_cap_A"] = cap
            decision["radius_margin_A"] = cap - distance
            if distance <= 1e-12:
                decision["reason"] = "COINCIDENT_COORDINATES_UNRESOLVED"
            elif abs(distance - cap) <= RADIUS_BOUNDARY_TOLERANCE_A:
                decision["reason"] = "NUMERICAL_RADIUS_BOUNDARY_UNRESOLVED"
            elif distance > cap:
                decision["reason"] = "OUTSIDE_PAIR_RADIUS_CAP"
            else:
                # Test the complete original heavy-atom point set, not a
                # sequentially filtered point set that would depend on order.
                for j in heavy:
                    if j == i:
                        continue
                    margin = float(np.dot(coords[j] - center, coords[j] - center)
                                   + np.dot(coords[j] - coords[i], coords[j] - coords[i])
                                   - distance * distance)
                    if margin < -GABRIEL_MARGIN_TOLERANCE_A2:
                        decision["blockers"].append({"atom_index0": j,
                                                     "element": mol.atoms[j].element,
                                                     "empty_ball_margin_A2": margin})
                decision["accepted"] = not decision["blockers"]
                decision["reason"] = "EMPTY_BALL_VISIBLE_CANDIDATE" if decision["accepted"] else "GEOMETRICALLY_OCCLUDED_CANDIDATE_NOT_PROVEN_NONBOND"
        decisions.append(decision)
        if decision["accepted"]:
            order = float(bo[metal_idx, i]) if bo is not None else 0.0
            neighbors.append(NeighborInfo(i, atom.element, distance, order,
                                          2 * order + 1 / max(distance, 0.5)))
    neighbors.sort(key=lambda x: x.get_sort_key())
    return neighbors, {"adapter_version": XYZ_ADAPTER_VERSION,
                       "metal_index0": metal_idx,
                       "bo_input_used": bo is not None,
                       "covalent_radius_factor": COVALENT_RADIUS_FACTOR,
                       "radius_boundary_tolerance_A": RADIUS_BOUNDARY_TOLERANCE_A,
                       "gabriel_margin_tolerance_A2": GABRIEL_MARGIN_TOLERANCE_A2,
                       "hydrogen_donors_and_blockers_excluded": True,
                       "blocker_domain": "ALL_ORIGINAL_NON_H_ATOMS_EXCEPT_CENTER_AND_CANDIDATE_INCLUDING_OTHER_METALS",
                       "decisions": decisions,
                       "warning": "Geometric candidates, not validated chemical bonds. Asymmetric pi contacts and small-bite/weak coordination can be excluded; inspect decisions."}


def get_coordination_neighbors(mol, metal_idx, bo_threshold=0.3, distance_cutoff=None):
    return perceive_neighbors(mol, metal_idx, bo_threshold, distance_cutoff)[0]


def stabilize_neighbors(neighbors: List[NeighborInfo], expected_cn=None):
    # Compatibility only; the XYZ-v2 benchmark never supplies a CN hint.
    return neighbors if expected_cn is None else neighbors[:expected_cn]


def get_donor_indices(neighbors):
    return [n.atom_idx for n in neighbors]


def validate_donors(mol, donor_indices):
    return ["Unusual or excluded donor element: " + mol.atoms[i].element
            for i in donor_indices if mol.atoms[i].element not in DONOR_ELEMENTS]
