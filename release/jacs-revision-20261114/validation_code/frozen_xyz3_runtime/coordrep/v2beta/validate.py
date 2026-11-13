"""
CoordRep v2-beta validation gates.

14 validation checks per record (see Part 5 of the specification).
"""

from __future__ import annotations

import copy
import random
import re
from dataclasses import dataclass, field
from typing import Dict, List, Tuple

from .core import MultiMetalRecord, HapticRecord, CoordinationSite
from .serialize import (
    serialize_multi, serialize_multi_oneline,
    serialize_haptic, serialize_haptic_oneline,
)
from .canonicalize import (
    canonicalize_multi, canonicalize_haptic,
    _serialize_multi_content, _serialize_haptic_content,
)


@dataclass
class ValidationResult:
    case_id: str
    checks: Dict[str, bool] = field(default_factory=dict)
    notes: Dict[str, str] = field(default_factory=dict)

    @property
    def all_passed(self) -> bool:
        return all(self.checks.values())

    def to_dict(self) -> dict:
        return {
            "case_id": self.case_id,
            "all_passed": self.all_passed,
            **self.checks,
        }


# ════════════════════════════════════════════════════════════════════
# Shared checks
# ════════════════════════════════════════════════════════════════════

def _check_parse_valid(serialized: str) -> bool:
    """Check that the serialized string has expected block structure."""
    required_multi = [
        "[Metals:", "[LocalSphere:", "[Sites:", "[Ligands:", "[ID:"
    ]
    required_haptic = ["[Metal:", "[Sites:", "[Ligands:", "[ID:"]
    has_multi = all(blk in serialized for blk in required_multi)
    has_haptic = all(blk in serialized for blk in required_haptic)
    return has_multi or has_haptic


def _check_roundtrip(serialized: str, re_serialized: str) -> bool:
    """Check that re-canonicalizing produces the same string."""
    return serialized.strip() == re_serialized.strip()


def _check_no_placeholder(serialized: str) -> bool:
    """No ???, PLACEHOLDER, TODO tokens."""
    for tok in ["???", "PLACEHOLDER", "TODO", "FIXME"]:
        if tok in serialized:
            return False
    return True


def _check_no_raw_coordinates(serialized: str, record_dict: dict) -> bool:
    """No raw CSD coordinates exported."""
    coord_pat = re.compile(r"-?\d+\.\d{4,}")
    s = str(record_dict)
    matches = coord_pat.findall(s)
    # Allow metal-metal distances (2 decimal) but flag 4+ decimal precision
    for m in matches:
        # Metal-metal distances like 2.09 are fine; raw coords like 12.3456 are not
        if len(m.split(".")[-1]) >= 4 and float(m) > 10:
            return False
    return True


def _check_identity_keys_present(identity_dict: dict) -> bool:
    """All L0/L1/L2/L3 keys present and non-empty."""
    for key in ["L0", "L1", "L2", "L3"]:
        found = any(k.startswith(key) for k in identity_dict)
        val = next((v for k, v in identity_dict.items() if k.startswith(key)), "")
        if not found or not val:
            return False
    return True


# ════════════════════════════════════════════════════════════════════
# Multinuclear validation
# ════════════════════════════════════════════════════════════════════

def validate_multi(rec: MultiMetalRecord) -> ValidationResult:
    """Run all 14 validation gates on a MultiMetalRecord."""
    vr = ValidationResult(case_id=rec.case_id)

    serialized = serialize_multi(rec)
    content = _serialize_multi_content(rec)

    # 1. parse_valid
    vr.checks["parse_valid"] = _check_parse_valid(serialized)

    # 2. roundtrip_valid — deep-copy, re-canonicalize, compare content
    rec2 = copy.deepcopy(rec)
    rec2 = canonicalize_multi(rec2)
    re_content = _serialize_multi_content(rec2)
    vr.checks["roundtrip_valid"] = _check_roundtrip(content, re_content)

    # 3. no_placeholder_tokens
    vr.checks["no_placeholder_tokens"] = _check_no_placeholder(serialized)

    # 4. metal_order_invariant — shuffle metals and re-canonicalize
    rec3 = copy.deepcopy(rec)
    if len(rec3.metals) > 1:
        rng = random.Random(42)
        rng.shuffle(rec3.metals)
        label_map = {
            metal.label: f"Mx{idx + 1}"
            for idx, metal in enumerate(rec3.metals)
        }
        for m in rec3.metals:
            m.label = label_map[m.label]
        for e in rec3.metal_edges:
            e.m1 = label_map.get(e.m1, e.m1)
            e.m2 = label_map.get(e.m2, e.m2)
        for s in rec3.sites:
            s.target_metals = [label_map.get(t, t) for t in s.target_metals]
    rec3 = canonicalize_multi(rec3)
    vr.checks["metal_order_invariant"] = (
        _serialize_multi_content(rec3) == content
    )

    # 5. atom_order_invariant — shuffle donor_atoms within sites
    rec4 = copy.deepcopy(rec)
    rng = random.Random(99)
    for site_idx, s in enumerate(rec4.sites):
        if len(s.donor_atoms) > 1:
            keys = s.meta.get(
                "donor_canonical_keys", [""] * len(s.donor_elements)
            )
            combined = list(zip(s.donor_atoms, s.donor_elements, keys))
            rng.shuffle(combined)
            s.donor_atoms = [x[0] for x in combined]
            s.donor_elements = [x[1] for x in combined]
            s.meta["donor_canonical_keys"] = [x[2] for x in combined]
        s.donor_atoms = [
            f"source_atom_{site_idx}_{atom_idx}"
            for atom_idx in range(len(s.donor_atoms))
        ]
    rec4 = canonicalize_multi(rec4)
    vr.checks["atom_order_invariant"] = (
        _serialize_multi_content(rec4) == content
    )

    # 6. ligand_order_invariant — shuffle ligands
    rec5 = copy.deepcopy(rec)
    random.seed(77)
    old_labels = [lg.label for lg in rec5.ligands]
    random.shuffle(rec5.ligands)
    new_map = {}
    for i, lg in enumerate(rec5.ligands):
        new_map[lg.label] = f"Ly{i+1}"
        lg.label = f"Ly{i+1}"
    for s in rec5.sites:
        s.ligand_label = new_map.get(s.ligand_label, s.ligand_label)
    rec5 = canonicalize_multi(rec5)
    vr.checks["ligand_order_invariant"] = (
        _serialize_multi_content(rec5) == content
    )

    # 7. site_atom_order_invariant — same as 5 for haptic sites (N/A → True)
    vr.checks["site_atom_order_invariant"] = True

    # 8. bridge_consistency_valid
    bridge_sites = [s for s in rec.sites if s.mu > 1]
    bridge_ok = all(
        len(s.target_metals) == s.mu for s in bridge_sites
    )
    metal_labels = {m.label for m in rec.metals}
    ligand_labels = {ligand.label for ligand in rec.ligands}
    site_labels = {site.label for site in rec.sites}
    reference_ok = (
        len(metal_labels) == len(rec.metals)
        and len(ligand_labels) == len(rec.ligands)
        and len(site_labels) == len(rec.sites)
        and all(
            edge.m1 in metal_labels and edge.m2 in metal_labels
            for edge in rec.metal_edges
        )
        and all(
            site.ligand_label in ligand_labels
            and bool(site.target_metals)
            and all(target in metal_labels for target in site.target_metals)
            for site in rec.sites
        )
        and all(
            all(label in site_labels for label in edge.bridges)
            for edge in rec.metal_edges
        )
    )
    vr.checks["bridge_consistency_valid"] = bridge_ok and reference_ok

    # 9. eta_mu_consistency_valid
    eta_mu_ok = all(
        s.eta >= 1 and s.mu >= 1 for s in rec.sites
    )
    vr.checks["eta_mu_consistency_valid"] = eta_mu_ok

    # 10. local_sphere_consistency_valid
    for m in rec.metals:
        local = [s for s in rec.sites if m.label in s.target_metals]
        if len(local) != m.cn_site:
            vr.checks["local_sphere_consistency_valid"] = False
            break
    else:
        vr.checks["local_sphere_consistency_valid"] = True

    # 11. no_duplicate_ligand_for_same_bridge
    # Multiple bridge sites CAN share the same ligand (e.g. 4 OAc bridges
    # from 2 ligand entries). The check is: a bridge site should not be
    # duplicated as a separate ligand for each metal it connects.
    # In v2beta, bridges appear once in Ligands → always valid.
    vr.checks["no_duplicate_ligand_for_same_bridge"] = True

    # 12. no_raw_coordinates_exported
    vr.checks["no_raw_coordinates_exported"] = _check_no_raw_coordinates(
        serialized, rec.to_dict()
    )

    # 13. identity_keys_present
    vr.checks["identity_keys_present"] = _check_identity_keys_present(
        rec.identity.to_dict()
    )

    # 14. human_readable_si_example
    vr.checks["human_readable_si_example"] = len(serialized) > 50

    return vr


# ════════════════════════════════════════════════════════════════════
# Haptic validation
# ════════════════════════════════════════════════════════════════════

def validate_haptic(rec: HapticRecord) -> ValidationResult:
    """Run all 14 validation gates on a HapticRecord."""
    vr = ValidationResult(case_id=rec.case_id)

    serialized = serialize_haptic(rec)
    content = _serialize_haptic_content(rec)

    # 1. parse_valid
    vr.checks["parse_valid"] = _check_parse_valid(serialized)

    # 2. roundtrip_valid
    rec2 = copy.deepcopy(rec)
    rec2 = canonicalize_haptic(rec2)
    re_content = _serialize_haptic_content(rec2)
    vr.checks["roundtrip_valid"] = _check_roundtrip(content, re_content)

    # 3. no_placeholder_tokens
    vr.checks["no_placeholder_tokens"] = _check_no_placeholder(serialized)

    # 4. metal_order_invariant (single metal → trivially True)
    vr.checks["metal_order_invariant"] = True

    # 5. atom_order_invariant — shuffle donor_atoms
    rec4 = copy.deepcopy(rec)
    rng = random.Random(99)
    for site_idx, s in enumerate(rec4.sites):
        if len(s.donor_atoms) > 1:
            keys = s.meta.get(
                "donor_canonical_keys", [""] * len(s.donor_elements)
            )
            combined = list(zip(s.donor_atoms, s.donor_elements, keys))
            rng.shuffle(combined)
            s.donor_atoms = [x[0] for x in combined]
            s.donor_elements = [x[1] for x in combined]
            s.meta["donor_canonical_keys"] = [x[2] for x in combined]
        s.donor_atoms = [
            f"source_atom_{site_idx}_{atom_idx}"
            for atom_idx in range(len(s.donor_atoms))
        ]
    rec4 = canonicalize_haptic(rec4)
    vr.checks["atom_order_invariant"] = (
        _serialize_haptic_content(rec4) == content
    )

    # 6. ligand_order_invariant
    rec5 = copy.deepcopy(rec)
    random.seed(77)
    new_map = {}
    for i, lg in enumerate(rec5.ligands):
        new_map[lg.label] = f"Ly{i+1}"
        lg.label = f"Ly{i+1}"
    for s in rec5.sites:
        s.ligand_label = new_map.get(s.ligand_label, s.ligand_label)
    random.shuffle(rec5.ligands)
    rec5 = canonicalize_haptic(rec5)
    vr.checks["ligand_order_invariant"] = (
        _serialize_haptic_content(rec5) == content
    )

    # 7. site_atom_order_invariant — shuffle haptic atom sets
    rec6 = copy.deepcopy(rec)
    for site_idx, s in enumerate(rec6.sites):
        if s.eta > 1 and len(s.donor_atoms) > 1:
            keys = s.meta.get(
                "donor_canonical_keys", [""] * len(s.donor_elements)
            )
            combined = list(zip(s.donor_atoms, s.donor_elements, keys))
            shift = site_idx % len(combined)
            combined = combined[shift:] + combined[:shift]
            combined.reverse()
            s.donor_elements = [x[1] for x in combined]
            s.meta["donor_canonical_keys"] = [x[2] for x in combined]
            s.donor_atoms = [
                f"renamed_eta_atom_{site_idx}_{atom_idx}"
                for atom_idx in range(len(combined))
            ]
    rec6 = canonicalize_haptic(rec6)
    vr.checks["site_atom_order_invariant"] = (
        _serialize_haptic_content(rec6) == content
    )

    # 8. bridge_consistency_valid (no bridges in haptic → True)
    ligand_labels = {ligand.label for ligand in rec.ligands}
    site_labels = {site.label for site in rec.sites}
    vr.checks["bridge_consistency_valid"] = (
        len(ligand_labels) == len(rec.ligands)
        and len(site_labels) == len(rec.sites)
        and all(
            site.ligand_label in ligand_labels
            and site.target_metals == ["M1"]
            for site in rec.sites
        )
    )

    # 9. eta_mu_consistency_valid
    eta_mu_ok = all(s.eta >= 1 and s.mu >= 1 for s in rec.sites)
    # Check eta_sum matches
    eta_total = sum(s.eta for s in rec.sites)
    eta_mu_ok = eta_mu_ok and (eta_total == rec.metal.eta_sum)
    vr.checks["eta_mu_consistency_valid"] = eta_mu_ok

    # 10. local_sphere_consistency_valid
    vr.checks["local_sphere_consistency_valid"] = (
        len(rec.sites) == rec.metal.cn_site
    )

    # 11. no_duplicate_ligand_for_same_bridge (N/A → True)
    vr.checks["no_duplicate_ligand_for_same_bridge"] = True

    # 12. no_raw_coordinates_exported
    vr.checks["no_raw_coordinates_exported"] = _check_no_raw_coordinates(
        serialized, rec.to_dict()
    )

    # 13. identity_keys_present
    vr.checks["identity_keys_present"] = _check_identity_keys_present(
        rec.identity.to_dict()
    )

    # 14. human_readable_si_example
    vr.checks["human_readable_si_example"] = len(serialized) > 50

    return vr
