"""
csd_adapter.py
==============
Convert CSD entries to RawMolecule objects for the CoordRep pipeline.

Handles:
  - Transition-metal detection
  - Mononuclear check
  - Hapticity / organometallic exclusion (η>1)
  - CN range check
  - Disorder / occupancy check
  - Conversion to RawMolecule → encode_molecule → CoordComplex → CoordRep string

License note: this module never exports raw CSD coordinates or CIF data.
Only aggregated statistics and refcode-level summaries are written.
"""

from __future__ import annotations

import re
import numpy as np
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

from coordrep.io.tmqm_reader import Atom, RawMolecule, TRANSITION_METALS
from coordrep.graph.neighbors import DONOR_ELEMENTS

# Hapticity indicators in CSD bond types
_HAPTIC_BOND_TYPES = {"Pi", "Delocalized", "pi", "delocalized"}


@dataclass
class CSDFilterResult:
    """Result of filtering one CSD entry."""
    refcode: str
    passed: bool
    rejection_reason: str = ""
    metal: str = ""
    cn: int = 0
    n_metals: int = 0
    has_disorder: bool = False
    is_polymeric: bool = False


def _count_metals(molecule) -> Tuple[List, int]:
    """Return (list_of_metal_atoms, count)."""
    metals = []
    for a in molecule.atoms:
        if a.atomic_symbol in TRANSITION_METALS:
            metals.append(a)
    return metals, len(metals)


def _has_hapticity(molecule, metal_atom) -> bool:
    """Check if metal has any η>1 (pi/delocalized) bonds."""
    for b in metal_atom.bonds:
        bt = str(b.bond_type) if b.bond_type else ""
        if bt in _HAPTIC_BOND_TYPES or "pi" in bt.lower():
            return True
    return False


def _get_donor_neighbors(molecule, metal_atom) -> List:
    """Get σ-bonded donor atoms around the metal."""
    donors = []
    for b in metal_atom.bonds:
        other = b.atoms[0] if b.atoms[1] == metal_atom else b.atoms[1]
        if other.atomic_symbol in DONOR_ELEMENTS or other.atomic_symbol == "H":
            # skip H donors for CN count (consistent with CoordRep pipeline)
            if other.atomic_symbol != "H":
                donors.append(other)
    return donors


def filter_csd_entry(entry) -> CSDFilterResult:
    """
    Apply the CoordRep filtering waterfall to one CSD entry.

    Steps:
      1. has 3D structure
      2. detect transition metal
      3. mononuclear check
      4. hapticity exclusion
      5. CN range [2..6]
      6. disorder check
      7. polymeric check

    Returns CSDFilterResult with pass/fail and reason.
    """
    refcode = entry.identifier

    # 1. 3D structure
    if not entry.has_3d_structure:
        return CSDFilterResult(refcode, False, "no_3d_structure")

    mol = entry.molecule
    if mol is None:
        return CSDFilterResult(refcode, False, "no_molecule")

    # 2. Transition metal
    metals, n_metals = _count_metals(mol)
    if n_metals == 0:
        return CSDFilterResult(refcode, False, "no_transition_metal")

    # 3. Mononuclear
    if n_metals > 1:
        return CSDFilterResult(refcode, False, "multinuclear",
                               metal=metals[0].atomic_symbol, n_metals=n_metals)

    metal_atom = metals[0]
    metal_elem = metal_atom.atomic_symbol

    # 4. Hapticity
    if _has_hapticity(mol, metal_atom):
        return CSDFilterResult(refcode, False, "hapticity",
                               metal=metal_elem, n_metals=1)

    # 5. CN range
    donors = _get_donor_neighbors(mol, metal_atom)
    cn = len(donors)
    if cn < 2 or cn > 6:
        return CSDFilterResult(refcode, False, f"cn_out_of_range_{cn}",
                               metal=metal_elem, cn=cn, n_metals=1)

    # 6. Disorder
    has_disorder = entry.has_disorder
    if has_disorder:
        return CSDFilterResult(refcode, False, "disorder",
                               metal=metal_elem, cn=cn, n_metals=1,
                               has_disorder=True)

    # 7. Polymeric
    is_poly = entry.is_polymeric
    if is_poly:
        return CSDFilterResult(refcode, False, "polymeric",
                               metal=metal_elem, cn=cn, n_metals=1,
                               is_polymeric=True)

    return CSDFilterResult(refcode, True, "",
                           metal=metal_elem, cn=cn, n_metals=1)


def csd_entry_to_raw_molecule(entry) -> Optional[RawMolecule]:
    """
    Convert a CSD entry to a RawMolecule.

    Returns None if coordinates are missing.
    Does NOT export raw coordinates — only creates an in-memory object
    for the encoding pipeline.
    """
    mol = entry.molecule
    if mol is None:
        return None

    atoms = []
    for i, a in enumerate(mol.atoms):
        c = a.coordinates
        if c is None:
            return None
        atoms.append(Atom(
            index=i,
            element=a.atomic_symbol,
            x=c.x,
            y=c.y,
            z=c.z,
        ))

    if not atoms:
        return None

    # Build a simple bond-order matrix (all 1.0 for σ-bonds)
    n = len(atoms)
    bo = np.zeros((n, n), dtype=np.float32)
    for b in mol.bonds:
        i0 = b.atoms[0].index
        i1 = b.atoms[1].index
        if i0 < n and i1 < n:
            order = 1.0
            bt = str(b.bond_type) if b.bond_type else ""
            if "Double" in bt or "double" in bt:
                order = 2.0
            elif "Triple" in bt or "triple" in bt:
                order = 3.0
            bo[i0, i1] = order
            bo[i1, i0] = order

    return RawMolecule(
        mol_id=entry.identifier,
        atoms=atoms,
        bond_orders=bo,
    )


def refcode_family(refcode: str) -> str:
    """
    Extract the base refcode family.

    CSD convention: ABCDEF, ABCDEF01, ABCDEF02, ...
    The family is the alphabetic prefix.
    """
    # Strip trailing digits that represent polymorph/redetermination number
    m = re.match(r'^([A-Z]+)(\d*)$', refcode)
    if m:
        return m.group(1)
    return refcode
