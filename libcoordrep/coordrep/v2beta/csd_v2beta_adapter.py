"""
CSD → CoordRep-v2beta automatic converter.

Takes a CSD entry and produces a MultiMetalRecord or HapticRecord.
No raw CSD coordinates are exported.
"""
from __future__ import annotations

import traceback
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Tuple

from coordrep.io.tmqm_reader import TRANSITION_METALS
from coordrep.graph.neighbors import DONOR_ELEMENTS

from .core import (
    CoordinationSite, HapticRecord, Ligand, MetalCenter, MetalEdge,
    MultiMetalRecord,
)
from .canonicalize import canonicalize_multi, canonicalize_haptic
from .validate import validate_multi, validate_haptic, ValidationResult


_HAPTIC_BOND_TYPES = {"Pi", "Delocalized", "pi", "delocalized"}


@dataclass
class V2BetaConversionResult:
    refcode: str
    record_type: str = ""          # "multi" | "haptic" | ""
    success: bool = False
    record: object = None          # MultiMetalRecord or HapticRecord
    validation: Optional[ValidationResult] = None
    failure_stage: str = ""
    failure_reason: str = ""
    metal_count: int = 0
    metals_str: str = ""
    bridge_count: int = 0
    haptic_class: str = ""
    eta_values: str = ""
    n_sites: int = 0


# ════════════════════════════════════════════════════════════════════
# Helpers
# ════════════════════════════════════════════════════════════════════

def _is_pi_bond(bond) -> bool:
    bt = str(bond.bond_type) if bond.bond_type else ""
    return bt in _HAPTIC_BOND_TYPES or "pi" in bt.lower() or "deloc" in bt.lower()


def _get_metal_atoms(mol) -> list:
    return [a for a in mol.atoms if a.atomic_symbol in TRANSITION_METALS]


def _atom_id(a) -> int:
    return id(a)


def _safe_smiles(mol, atoms_subset=None) -> str:
    """Return a simple element-based placeholder SMILES (no CSD coords)."""
    if atoms_subset and len(atoms_subset) <= 3:
        return ".".join(f"[{a.atomic_symbol}]" for a in atoms_subset)
    return "[*]"


# ════════════════════════════════════════════════════════════════════
# Multinuclear converter
# ════════════════════════════════════════════════════════════════════

def convert_multinuclear(entry) -> V2BetaConversionResult:
    """Convert a multinuclear CSD entry to MultiMetalRecord."""
    refcode = entry.identifier
    res = V2BetaConversionResult(refcode=refcode, record_type="multi")

    try:
        mol = entry.molecule
        if mol is None:
            res.failure_stage = "molecule"
            res.failure_reason = "no_molecule"
            return res

        metal_atoms = _get_metal_atoms(mol)
        n_metals = len(metal_atoms)
        res.metal_count = n_metals
        res.metals_str = "+".join(sorted(a.atomic_symbol for a in metal_atoms))

        if n_metals < 2:
            res.failure_stage = "metal_count"
            res.failure_reason = "not_multinuclear"
            return res

        if n_metals > 12:
            res.failure_stage = "metal_count"
            res.failure_reason = f"too_many_metals_{n_metals}"
            return res

        metal_id_set = {_atom_id(m) for m in metal_atoms}
        metal_by_id = {_atom_id(m): m for m in metal_atoms}

        # Build per-metal neighbor lists
        metal_labels = {}
        for i, m in enumerate(metal_atoms):
            metal_labels[_atom_id(m)] = f"M{i+1}"

        # Find all non-metal neighbors of each metal, and bridging atoms
        metal_neighbors: Dict[int, List] = defaultdict(list)   # metal_id -> [atom, ...]
        atom_metal_map: Dict[int, List[int]] = defaultdict(list)  # atom_id -> [metal_ids]
        metal_metal_bonds: List[Tuple[int, int]] = []

        for m in metal_atoms:
            for b in m.bonds:
                other = b.atoms[0] if b.atoms[1] == m else b.atoms[1]
                if _atom_id(other) in metal_id_set:
                    pair = tuple(sorted([_atom_id(m), _atom_id(other)]))
                    if pair not in metal_metal_bonds:
                        metal_metal_bonds.append(pair)
                else:
                    if not _is_pi_bond(b):
                        metal_neighbors[_atom_id(m)].append(other)
                        atom_metal_map[_atom_id(other)].append(_atom_id(m))

        # Build metals
        metals = []
        for i, m in enumerate(metal_atoms):
            neigh = metal_neighbors[_atom_id(m)]
            metals.append(MetalCenter(
                label=f"M{i+1}",
                element=m.atomic_symbol,
                oxidation=None,
                dcount=None,
                cn_site=len(neigh),
                cn_atom=len(neigh),
                eta_sum=len(neigh),
                local_shape_best="",
            ))

        # Build edges
        edges = []
        # Metal-metal direct bonds
        mm_bond_set = set()
        for (mid1, mid2) in metal_metal_bonds:
            l1 = metal_labels[mid1]
            l2 = metal_labels[mid2]
            mm_bond_set.add((min(l1, l2), max(l1, l2)))
            edges.append(MetalEdge(
                m1=min(l1, l2), m2=max(l1, l2),
                relation="direct_MM_bond", mm_bond="yes",
            ))

        # Bridging connections (atoms bonded to >=2 metals)
        bridge_pairs = set()
        for aid, mids in atom_metal_map.items():
            if len(mids) >= 2:
                for i_m in range(len(mids)):
                    for j_m in range(i_m + 1, len(mids)):
                        l1 = metal_labels[mids[i_m]]
                        l2 = metal_labels[mids[j_m]]
                        pair = (min(l1, l2), max(l1, l2))
                        if pair not in mm_bond_set and pair not in bridge_pairs:
                            bridge_pairs.add(pair)
                            edges.append(MetalEdge(
                                m1=pair[0], m2=pair[1],
                                relation="bridged", mm_bond="no",
                            ))

        # Build sites and ligands
        sites = []
        ligands = []
        lig_counter = 0
        site_counter = 0
        seen_atoms = set()

        for m in metal_atoms:
            mid = _atom_id(m)
            mlabel = metal_labels[mid]
            for neigh_atom in metal_neighbors[mid]:
                aid = _atom_id(neigh_atom)
                if aid in seen_atoms:
                    # Already created a site for this atom
                    # Find it and add this metal to targets
                    for s in sites:
                        if _atom_id(neigh_atom) == s.meta.get("_atom_id"):
                            if mlabel not in s.target_metals:
                                s.target_metals.append(mlabel)
                            break
                    continue

                seen_atoms.add(aid)
                target_metals_for_atom = [
                    metal_labels[m_id] for m_id in atom_metal_map[aid]
                ]
                mu = len(target_metals_for_atom)

                lig_counter += 1
                site_counter += 1
                mode = ""
                if mu > 1:
                    mode = f"mu{mu}-{neigh_atom.atomic_symbol}"

                ligands.append(Ligand(
                    label=f"L{lig_counter}",
                    smiles=f"[{neigh_atom.atomic_symbol}]",
                    dent=1,
                    charge=None,
                ))
                sites.append(CoordinationSite(
                    label=f"S{site_counter}",
                    site_type="atom",
                    ligand_label=f"L{lig_counter}",
                    donor_atoms=[neigh_atom.label],
                    donor_elements=[neigh_atom.atomic_symbol],
                    eta=1, mu=mu,
                    target_metals=sorted(target_metals_for_atom),
                    mode=mode,
                    meta={"_atom_id": aid},
                ))

        res.bridge_count = sum(1 for s in sites if s.mu > 1)
        res.n_sites = len(sites)

        # Fix cn_site on metals after site construction
        for m_obj in metals:
            local = [s for s in sites if m_obj.label in s.target_metals]
            m_obj.cn_site = len(local)
            m_obj.cn_atom = len(local)
            m_obj.eta_sum = len(local)
            m_obj.local_sites = [s.label for s in local]

        if len(sites) == 0:
            res.failure_stage = "site_detection"
            res.failure_reason = "no_sites_found"
            return res

        rec = MultiMetalRecord(
            case_id=f"csd_multi_{refcode}",
            refcode=refcode,
            description=f"CSD {refcode}: {n_metals}-nuclear {res.metals_str}",
            metals=metals, metal_edges=edges, sites=sites, ligands=ligands,
        )

        rec = canonicalize_multi(rec)
        vr = validate_multi(rec)

        res.record = rec
        res.validation = vr
        res.success = vr.all_passed

        if not vr.all_passed:
            failed = [k for k, v in vr.checks.items() if not v]
            res.failure_stage = "validation"
            res.failure_reason = ";".join(failed)

    except Exception as exc:
        res.failure_stage = "exception"
        res.failure_reason = str(exc)[:200]

    return res


# ════════════════════════════════════════════════════════════════════
# Haptic converter
# ════════════════════════════════════════════════════════════════════

def convert_haptic(entry) -> V2BetaConversionResult:
    """Convert a haptic CSD entry to HapticRecord."""
    refcode = entry.identifier
    res = V2BetaConversionResult(refcode=refcode, record_type="haptic")

    try:
        mol = entry.molecule
        if mol is None:
            res.failure_stage = "molecule"
            res.failure_reason = "no_molecule"
            return res

        metal_atoms = _get_metal_atoms(mol)
        if len(metal_atoms) != 1:
            res.failure_stage = "metal_count"
            res.failure_reason = f"expected_1_metal_got_{len(metal_atoms)}"
            return res

        m_atom = metal_atoms[0]
        res.metal_count = 1
        res.metals_str = m_atom.atomic_symbol

        # Separate pi and sigma bonds from metal
        pi_atoms = []
        sigma_atoms = []
        for b in m_atom.bonds:
            other = b.atoms[0] if b.atoms[1] == m_atom else b.atoms[1]
            if other.atomic_symbol in TRANSITION_METALS:
                continue
            if _is_pi_bond(b):
                pi_atoms.append(other)
            else:
                sigma_atoms.append(other)

        if not pi_atoms:
            res.failure_stage = "haptic_detection"
            res.failure_reason = "no_pi_bonds_found"
            return res

        # Group pi atoms into connected fragments
        pi_id_set = {_atom_id(a) for a in pi_atoms}
        visited = set()
        fragments = []

        def _bfs_fragment(start_atom):
            frag = []
            queue = [start_atom]
            while queue:
                cur = queue.pop(0)
                cid = _atom_id(cur)
                if cid in visited:
                    continue
                visited.add(cid)
                frag.append(cur)
                for b in cur.bonds:
                    other = b.atoms[0] if b.atoms[1] == cur else b.atoms[1]
                    oid = _atom_id(other)
                    if oid in pi_id_set and oid not in visited:
                        queue.append(other)
            return frag

        for pa in pi_atoms:
            if _atom_id(pa) not in visited:
                frag = _bfs_fragment(pa)
                if frag:
                    fragments.append(frag)

        # Build sites
        sites = []
        ligands = []
        site_counter = 0
        lig_counter = 0
        eta_values = []

        for frag in fragments:
            eta = len(frag)
            eta_values.append(eta)
            lig_counter += 1
            site_counter += 1

            atoms_labels = [a.label for a in frag]
            atoms_elems = [a.atomic_symbol for a in frag]
            mode = f"eta{eta}" if eta > 1 else ""
            stype = "pi_fragment" if eta > 1 else "atom"

            smiles = ".".join(f"[{a.atomic_symbol}]" for a in frag[:6])

            ligands.append(Ligand(
                label=f"L{lig_counter}",
                smiles=smiles,
                dent=eta,
                charge=None,
            ))
            sites.append(CoordinationSite(
                label=f"S{site_counter}",
                site_type=stype,
                ligand_label=f"L{lig_counter}",
                donor_atoms=atoms_labels,
                donor_elements=atoms_elems,
                eta=eta, mu=1,
                target_metals=["M1"],
                mode=mode,
            ))

        # Add sigma donors
        for sa in sigma_atoms:
            lig_counter += 1
            site_counter += 1
            ligands.append(Ligand(
                label=f"L{lig_counter}",
                smiles=f"[{sa.atomic_symbol}]",
                dent=1, charge=None,
            ))
            sites.append(CoordinationSite(
                label=f"S{site_counter}",
                site_type="atom",
                ligand_label=f"L{lig_counter}",
                donor_atoms=[sa.label],
                donor_elements=[sa.atomic_symbol],
                eta=1, mu=1,
                target_metals=["M1"],
                mode="",
            ))

        eta_sum = sum(s.eta for s in sites)
        cn_site = len(sites)

        metal = MetalCenter(
            label="M1",
            element=m_atom.atomic_symbol,
            oxidation=None,
            dcount=None,
            cn_site=cn_site,
            cn_atom=eta_sum,
            eta_sum=eta_sum,
            local_shape_best="",
        )
        metal.local_sites = [s.label for s in sites]

        res.n_sites = cn_site
        res.eta_values = ",".join(str(e) for e in sorted(eta_values))
        res.haptic_class = f"eta{max(eta_values)}" if eta_values else ""

        rec = HapticRecord(
            case_id=f"csd_haptic_{refcode}",
            refcode=refcode,
            description=f"CSD {refcode}: {m_atom.atomic_symbol} haptic {res.haptic_class}",
            metal=metal, sites=sites, ligands=ligands,
        )

        rec = canonicalize_haptic(rec)
        vr = validate_haptic(rec)

        res.record = rec
        res.validation = vr
        res.success = vr.all_passed

        if not vr.all_passed:
            failed = [k for k, v in vr.checks.items() if not v]
            res.failure_stage = "validation"
            res.failure_reason = ";".join(failed)

    except Exception as exc:
        res.failure_stage = "exception"
        res.failure_reason = str(exc)[:200]

    return res


# ════════════════════════════════════════════════════════════════════
# Classifier
# ════════════════════════════════════════════════════════════════════

def classify_csd_entry(entry) -> str:
    """Classify a CSD entry: 'multi', 'haptic', 'mono_eta1', 'no_tm', 'no_3d', 'other'."""
    if not entry.has_3d_structure:
        return "no_3d"
    mol = entry.molecule
    if mol is None:
        return "no_molecule"

    metals = _get_metal_atoms(mol)
    if len(metals) == 0:
        return "no_tm"
    if len(metals) >= 2:
        return "multi"

    # Check haptic
    m = metals[0]
    for b in m.bonds:
        if _is_pi_bond(b):
            return "haptic"

    return "mono_eta1"
