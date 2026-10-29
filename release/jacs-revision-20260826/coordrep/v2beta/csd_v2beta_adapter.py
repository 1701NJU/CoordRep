"""
CSD → CoordRep-v2beta automatic converter.

Takes a CSD entry and produces a MultiMetalRecord or HapticRecord.
No raw CSD coordinates are exported.
"""
from __future__ import annotations

import hashlib
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


def _atom_id(a) -> str:
    """Use atom label as stable identity (CSD may return different objects for same atom)."""
    return a.label


def _safe_smiles(mol, atoms_subset=None) -> str:
    """Return a simple element-based placeholder SMILES (no CSD coords)."""
    if atoms_subset and len(atoms_subset) <= 3:
        return ".".join(f"[{a.atomic_symbol}]" for a in atoms_subset)
    return "[*]"


def _atom_topology_signatures(mol) -> Dict[str, str]:
    """Build label-free atom-graph refinement signatures.

    CSD atom labels are used only to look signatures up while adapting the
    source object; the signature itself contains element/charge/degree, bond
    types, and recursively refined neighbour colours.  Symmetry-equivalent
    atoms intentionally retain the same signature.  If labels are not unique,
    the function returns an empty mapping and canonicalization falls back to
    the donor-element multiset rather than treating an arbitrary label as
    chemical identity.
    """
    atoms = list(mol.atoms)
    labels = [str(atom.label) for atom in atoms]
    if len(labels) != len(set(labels)):
        return {}

    atom_by_label = {str(atom.label): atom for atom in atoms}

    def formal_charge(atom) -> str:
        for attr in ("formal_charge", "charge"):
            try:
                value = getattr(atom, attr)
            except Exception:
                continue
            if value is not None:
                return str(value)
        return "unk"

    colours = {
        label: repr((
            atom.atomic_symbol,
            formal_charge(atom),
            len(list(atom.bonds)),
        ))
        for label, atom in atom_by_label.items()
    }

    for _ in range(max(2, len(atoms))):
        new_colours = {}
        for label, atom in atom_by_label.items():
            neighbours = []
            for bond in atom.bonds:
                other = bond.atoms[0] if bond.atoms[1] == atom else bond.atoms[1]
                other_label = str(other.label)
                if other_label not in colours:
                    continue
                bond_type = str(bond.bond_type) if bond.bond_type is not None else "unknown"
                neighbours.append((bond_type, colours[other_label]))
            payload = repr((colours[label], tuple(sorted(neighbours))))
            new_colours[label] = hashlib.sha256(payload.encode("utf-8")).hexdigest()
        if new_colours == colours:
            break
        colours = new_colours
    return colours


# ════════════════════════════════════════════════════════════════════
# Multinuclear converter
# ════════════════════════════════════════════════════════════════════

def convert_multinuclear(entry) -> V2BetaConversionResult:
    """Convert a multinuclear CSD entry to MultiMetalRecord.

    Identity strategy:
    - Metals: use Python id() — we iterate metal_atoms directly, guaranteed
      unique objects per metal center.
    - Donor atoms: use atom.label — accessed through bonds, CSD may return
      different Python objects for the same atom, but labels are stable.
      When CSD has duplicate labels (symmetry-generated), donor labels
      are made unique by appending the metal index context.
    """
    refcode = entry.identifier
    res = V2BetaConversionResult(refcode=refcode, record_type="multi")

    try:
        mol = entry.molecule
        if mol is None:
            res.failure_stage = "molecule"
            res.failure_reason = "no_molecule"
            return res

        atom_topology_keys = _atom_topology_signatures(mol)

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

        # Use Python id() for metals (direct iteration → unique objects)
        metal_pyid_set = {id(m) for m in metal_atoms}
        metal_label_map = {}  # id(metal_obj) → "M{i}"
        for i, m in enumerate(metal_atoms):
            metal_label_map[id(m)] = f"M{i+1}"

        # Build metal label set for neighbor-is-metal check via labels
        # (bond traversal may return different objects)
        metal_label_set = {m.label for m in metal_atoms}

        # Find all non-metal neighbors of each metal, and bridging atoms
        # Use label-based keys for donor atoms to handle CSD object aliasing
        metal_neighbors: Dict[int, List] = defaultdict(list)  # id(metal) → [atom, ...]
        atom_metal_map: Dict[str, List[int]] = defaultdict(list)  # donor_label → [id(metal), ...]
        metal_metal_bonds: List[Tuple[int, int]] = []

        for m in metal_atoms:
            mid = id(m)
            for b in m.bonds:
                other = b.atoms[0] if b.atoms[1] == m else b.atoms[1]
                if other.atomic_symbol in TRANSITION_METALS and other.label in metal_label_set:
                    # Metal-metal bond: find the target metal object by label
                    # Match to the correct metal_atom by label; if dup labels,
                    # pick the one that isn't this metal
                    for m2 in metal_atoms:
                        if m2.label == other.label and id(m2) != mid:
                            pair = tuple(sorted([mid, id(m2)]))
                            if pair not in metal_metal_bonds:
                                metal_metal_bonds.append(pair)
                            break
                else:
                    if not _is_pi_bond(b):
                        metal_neighbors[mid].append(other)
                        if mid not in atom_metal_map[other.label]:
                            atom_metal_map[other.label].append(mid)

        # Build metals
        metals = []
        for i, m in enumerate(metal_atoms):
            neigh = metal_neighbors[id(m)]
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
        mm_bond_set = set()
        for (mid1, mid2) in metal_metal_bonds:
            l1 = metal_label_map[mid1]
            l2 = metal_label_map[mid2]
            mm_bond_set.add((min(l1, l2), max(l1, l2)))
            edges.append(MetalEdge(
                m1=min(l1, l2), m2=max(l1, l2),
                relation="direct_MM_bond", mm_bond="yes",
            ))

        # Bridging connections (atoms bonded to >=2 distinct metals)
        bridge_pairs = set()
        for donor_label, mids in atom_metal_map.items():
            unique_mids = list(dict.fromkeys(mids))  # preserve order, dedup
            if len(unique_mids) >= 2:
                for i_m in range(len(unique_mids)):
                    for j_m in range(i_m + 1, len(unique_mids)):
                        l1 = metal_label_map[unique_mids[i_m]]
                        l2 = metal_label_map[unique_mids[j_m]]
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
        seen_donor_labels = set()

        for m in metal_atoms:
            mid = id(m)
            mlabel = metal_label_map[mid]
            for neigh_atom in metal_neighbors[mid]:
                dlabel = neigh_atom.label
                if dlabel in seen_donor_labels:
                    # Already created a site for this donor
                    for s in sites:
                        if s.meta.get("_donor_label") == dlabel:
                            if mlabel not in s.target_metals:
                                s.target_metals.append(mlabel)
                            break
                    continue

                seen_donor_labels.add(dlabel)
                target_metals_for_atom = [
                    metal_label_map[m_id] for m_id in dict.fromkeys(atom_metal_map[dlabel])
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
                    meta={
                        "_donor_label": dlabel,
                        "source_donor_labels": [neigh_atom.label],
                        "donor_canonical_keys": [
                            atom_topology_keys.get(str(neigh_atom.label), "")
                        ],
                    },
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

        atom_topology_keys = _atom_topology_signatures(mol)

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

        # Group pi atoms into connected fragments using label-based BFS
        pi_label_set = {a.label for a in pi_atoms}
        pi_by_label = {a.label: a for a in pi_atoms}
        visited_labels = set()
        fragments = []

        def _bfs_fragment(start_atom):
            frag = []
            queue = [start_atom]
            while queue:
                cur = queue.pop(0)
                clabel = cur.label
                if clabel in visited_labels:
                    continue
                visited_labels.add(clabel)
                frag.append(cur)
                for b in cur.bonds:
                    other = b.atoms[0] if b.atoms[1] == cur else b.atoms[1]
                    olabel = other.label
                    if olabel in pi_label_set and olabel not in visited_labels:
                        # Use canonical atom from pi_by_label
                        queue.append(pi_by_label[olabel])
            return frag

        for pa in pi_atoms:
            if pa.label not in visited_labels:
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
                meta={
                    "source_donor_labels": list(atoms_labels),
                    "donor_canonical_keys": [
                        atom_topology_keys.get(str(label), "")
                        for label in atoms_labels
                    ],
                },
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
                meta={
                    "source_donor_labels": [sa.label],
                    "donor_canonical_keys": [
                        atom_topology_keys.get(str(sa.label), "")
                    ],
                },
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
