"""Primary periodic-v3 site gate and MID/SID extraction from one CCDC crystal.

This is a clean, public-ZIP-facing port of the frozen v1/v2 primary site
eligibility logic. It deliberately omits secondary CShM and distance-QC
reporting, which were not primary emission gates. No source rows are written.
"""

from __future__ import annotations

import math
from collections import Counter, defaultdict
from typing import Any, Dict, List, Mapping, Tuple

import numpy as np

from .canonical_csd_quotient import canonical_periodic_groups
from .canonical_periodic_site import PeriodicDonorImage, canonicalize_fractional_site


TRANSITION_METALS = {
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "La", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
}
DONOR_ELEMENTS = {"N", "O", "S", "P", "Cl", "Br", "I", "F", "C", "Se", "Te", "As"}
DISTANCE_MIN_ANGSTROM = 0.60
DISTANCE_MAX_ANGSTROM = 3.50
VECTOR_WITNESS_TOLERANCE_ANGSTROM = 0.001


class _UnionFind:
    def __init__(self, values):
        self.parent = {value: value for value in values}

    def find(self, value):
        parent = self.parent[value]
        if parent != value:
            self.parent[value] = self.find(parent)
        return self.parent[value]

    def union(self, left, right):
        left_root = self.find(left)
        right_root = self.find(right)
        if left_root != right_root:
            self.parent[right_root] = left_root


def _independent_cell_matrix(crystal: Any) -> np.ndarray:
    lengths = crystal.cell_lengths
    angles = crystal.cell_angles
    a, b, c = float(lengths.a), float(lengths.b), float(lengths.c)
    alpha = math.radians(float(angles.alpha))
    beta = math.radians(float(angles.beta))
    gamma = math.radians(float(angles.gamma))
    sin_gamma = math.sin(gamma)
    if abs(sin_gamma) < 1e-10:
        raise ValueError("Degenerate unit-cell gamma")
    a_vector = np.array([a, 0.0, 0.0], dtype=float)
    b_vector = np.array([b * math.cos(gamma), b * sin_gamma, 0.0], dtype=float)
    c_x = c * math.cos(beta)
    c_y = c * (math.cos(alpha) - math.cos(beta) * math.cos(gamma)) / sin_gamma
    c_z_sq = c * c - c_x * c_x - c_y * c_y
    if c_z_sq < -1e-7:
        raise ValueError("Invalid unit-cell metric")
    c_vector = np.array([c_x, c_y, math.sqrt(max(c_z_sq, 0.0))])
    return np.vstack([a_vector, b_vector, c_vector])


def _edge_key(donor_group: int, donor_translation, metal_translation) -> Tuple[int, int, int, int]:
    return (int(donor_group),) + tuple(
        int(donor_translation[index]) - int(metal_translation[index])
        for index in range(3)
    )


def extract_entry(crystal: Any, *, suspect_chemistry: bool) -> Dict[str, Any]:
    """Return safe in-memory site keys and gate counts for one parsed P1 CIF.

    `states` contains only hashes, CN and a periodicity flag. Callers must not
    publish a per-refcode list without separately considering data permissions.
    """
    molecule = crystal.molecule
    if molecule is None:
        raise ValueError("Parsed crystal lacks molecule")
    cell = np.asarray(crystal.fractional_to_orthogonal.rotation, dtype=float)
    if not np.allclose(_independent_cell_matrix(crystal), cell, rtol=0.0, atol=1e-6):
        raise ValueError("Independent unit-cell matrix differs from CCDC transform")
    atom_info, atom_to_group, groups = canonical_periodic_groups(molecule, cell)
    metal_group_ids = sorted(
        group_id for group_id, group in groups.items()
        if group["element"] in TRANSITION_METALS
    )
    carbon_group_ids = {
        group_id for group_id, group in groups.items() if group["element"] == "C"
    }
    carbon_union = _UnionFind(carbon_group_ids)
    site_edges: Dict[int, Dict[Tuple[int, int, int, int], Dict[str, Any]]] = defaultdict(dict)
    pi_or_delocalized_by_metal_group = defaultdict(bool)

    for bond in molecule.bonds:
        left, right = bond.atoms
        left_group = atom_to_group[left.index]
        right_group = atom_to_group[right.index]
        if (
            left_group in carbon_group_ids
            and right_group in carbon_group_ids
            and left_group != right_group
        ):
            carbon_union.union(left_group, right_group)
        orientations = []
        if left.atomic_symbol in TRANSITION_METALS and right.atomic_symbol in DONOR_ELEMENTS:
            orientations.append((left, right))
        if right.atomic_symbol in TRANSITION_METALS and left.atomic_symbol in DONOR_ELEMENTS:
            orientations.append((right, left))
        if not orientations:
            continue
        bond_type = str(bond.bond_type) if bond.bond_type is not None else ""
        haptic_bond_type = "pi" in bond_type.lower() or "deloc" in bond_type.lower()
        for metal_atom, donor_atom in orientations:
            metal_info = atom_info[metal_atom.index]
            donor_info = atom_info[donor_atom.index]
            metal_group = atom_to_group[metal_atom.index]
            donor_group = atom_to_group[donor_atom.index]
            key = _edge_key(
                donor_group, donor_info["translation"], metal_info["translation"]
            )
            vector = np.asarray(donor_info["cartesian"]) - np.asarray(metal_info["cartesian"])
            edge = site_edges[metal_group].setdefault(
                key,
                {
                    "donor_element": donor_atom.atomic_symbol,
                    "vectors": [],
                    "bond_roles": set(),
                },
            )
            edge["vectors"].append(vector)
            edge["bond_roles"].add("haptic" if haptic_bond_type else "sigma")
            if haptic_bond_type:
                pi_or_delocalized_by_metal_group[metal_group] = True

    states: List[Dict[str, Any]] = []
    gate_counts: Counter[str] = Counter()
    nonzero_translation_sites_all = 0
    for metal_group in metal_group_ids:
        metal_record = groups[metal_group]
        metal_wrapped = np.asarray(metal_record["wrapped_fractional"], dtype=float)
        edges = site_edges.get(metal_group, {})
        sorted_keys = sorted(edges)
        cn = len(sorted_keys)
        has_nonzero_translation = any(
            any(value != 0 for value in key[1:]) for key in sorted_keys
        )
        nonzero_translation_sites_all += int(has_nonzero_translation)
        distances = []
        witness_residuals = []
        donor_elements = []
        for key in sorted_keys:
            donor_wrapped = np.asarray(groups[key[0]]["wrapped_fractional"], dtype=float)
            vector = (donor_wrapped + np.asarray(key[1:], dtype=float) - metal_wrapped) @ cell
            donor_elements.append(str(edges[key]["donor_element"]))
            distances.append(float(np.linalg.norm(vector)))
            observations = np.asarray(edges[key]["vectors"], dtype=float)
            witness_residuals.append(float(np.max(np.linalg.norm(observations - vector[None, :], axis=1))))

        carbon_edges = [key for key in sorted_keys if edges[key]["donor_element"] == "C"]
        carbon_components = Counter(carbon_union.find(key[0]) for key in carbon_edges)
        connected_multi_carbon = any(count >= 2 for count in carbon_components.values())
        repeated_periodic_carbon = len({key[0] for key in carbon_edges}) < len(carbon_edges)
        haptic_suspect = bool(
            pi_or_delocalized_by_metal_group.get(metal_group, False)
            or connected_multi_carbon
            or repeated_periodic_carbon
        )
        extreme_explicit_distance = any(
            distance < DISTANCE_MIN_ANGSTROM or distance > DISTANCE_MAX_ANGSTROM
            for distance in distances
        )
        edge_role_conflict = any(len(edges[key]["bond_roles"]) > 1 for key in sorted_keys)
        vector_witness_pass = max(witness_residuals, default=0.0) <= VECTOR_WITNESS_TOLERANCE_ANGSTROM
        chemical_domain_eligible = bool(
            2 <= cn <= 6
            and not haptic_suspect
            and not suspect_chemistry
            and not extreme_explicit_distance
            and not edge_role_conflict
        )
        emitted = chemical_domain_eligible and vector_witness_pass
        if not emitted:
            gate_counts["nonemitted_sites"] += 1
            if not 2 <= cn <= 6:
                gate_counts["cn_outside_2_6"] += 1
            if haptic_suspect:
                gate_counts["haptic_suspect"] += 1
            if suspect_chemistry:
                gate_counts["suspect_chemistry"] += 1
            if extreme_explicit_distance:
                gate_counts["distance_outside_0p60_3p50"] += 1
            if edge_role_conflict:
                gate_counts["bond_role_conflict"] += 1
            if not vector_witness_pass:
                gate_counts["vector_witness_failed"] += 1
            continue
        donors = [
            PeriodicDonorImage(
                element=donor_elements[index],
                fractional_base=tuple(float(value) for value in groups[key[0]]["wrapped_fractional"]),
                translation=tuple(int(value) for value in key[1:]),
                orbit=int(key[0]),
                bond_role="sigma",
            )
            for index, key in enumerate(sorted_keys)
        ]
        state = canonicalize_fractional_site(
            metal=str(metal_record["element"]),
            metal_fractional_base=tuple(float(value) for value in metal_wrapped),
            metal_translation=(0, 0, 0),
            donors=donors,
            cell_matrix=cell,
            decimals=6,
        )
        states.append(
            {
                "metric_id": state.metric_id,
                "stereo_id": state.stereo_id,
                "cn": cn,
                "metal": str(metal_record["element"]),
                "has_nonzero_translation": has_nonzero_translation,
            }
        )
        gate_counts["emitted_sites"] += 1
    return {
        "periodic_atom_orbits": len(groups),
        "metal_sites": len(metal_group_ids),
        "nonzero_translation_sites_all": nonzero_translation_sites_all,
        "states": states,
        "gate_counts": dict(gate_counts),
    }
