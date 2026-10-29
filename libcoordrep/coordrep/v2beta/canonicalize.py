"""
CoordRep v2-beta canonicalization.

Metal-node ordering uses Weisfeiler-Lehman-like iterative refinement.
Site ordering uses site_signature.
Ligand ordering uses graph-hash + site-list signature.
"""

from __future__ import annotations

import hashlib
from typing import Dict, List

from .core import (
    CoordinationSite, HapticRecord, Ligand, MetalCenter, MetalEdge,
    MultiMetalRecord, MultiIdentityKeys, HapticIdentityKeys,
)


def _short_hash(s: str, n: int = 12) -> str:
    return hashlib.sha256(s.encode()).hexdigest()[:n]


# ════════════════════════════════════════════════════════════════════
# Metal node canonical ordering (WL-like)
# ════════════════════════════════════════════════════════════════════

def _wl_metal_signatures(
    metals: List[MetalCenter],
    edges: List[MetalEdge],
    sites: List[CoordinationSite],
    n_iter: int = 3,
) -> Dict[str, str]:
    """Compute stable canonical signatures for metal nodes."""
    # Build adjacency
    adj: Dict[str, List[str]] = {m.label: [] for m in metals}
    for e in edges:
        adj[e.m1].append(e.m2)
        adj[e.m2].append(e.m1)

    # Initial labels
    labels = {}
    for m in metals:
        bridge_sites = [s for s in sites
                        if m.label in s.target_metals and s.mu > 1]
        bridge_sig = tuple(sorted(s.mode for s in bridge_sites))
        local_sites_sorted = tuple(sorted(
            s.site_signature() for s in sites if m.label in s.target_metals
        ))
        labels[m.label] = str((m.chemical_signature(), bridge_sig, local_sites_sorted))

    # WL iterations
    for _ in range(n_iter):
        new_labels = {}
        for m in metals:
            neighbor_labels = tuple(sorted(labels[n] for n in adj[m.label]))
            new_labels[m.label] = str((labels[m.label], neighbor_labels))
        labels = new_labels

    return labels


def _canonical_metal_order(
    metals: List[MetalCenter],
    edges: List[MetalEdge],
    sites: List[CoordinationSite],
) -> List[MetalCenter]:
    """Return metals sorted by WL canonical signature."""
    sigs = _wl_metal_signatures(metals, edges, sites)
    return sorted(metals, key=lambda m: sigs[m.label])


def _relabel_metals(
    metals: List[MetalCenter],
    edges: List[MetalEdge],
    sites: List[CoordinationSite],
    ligands: List[Ligand],
) -> tuple:
    """Relabel metals M1..Mn in canonical order; update all references."""
    ordered = _canonical_metal_order(metals, edges, sites)
    mapping = {m.label: f"M{i+1}" for i, m in enumerate(ordered)}

    for m in ordered:
        m.label = mapping[m.label]
    for e in edges:
        e.m1 = mapping[e.m1]
        e.m2 = mapping[e.m2]
        if e.m1 > e.m2:
            e.m1, e.m2 = e.m2, e.m1
    for s in sites:
        s.target_metals = sorted(mapping[t] for t in s.target_metals)
    return ordered, edges, sites, ligands


# ════════════════════════════════════════════════════════════════════
# Site canonical ordering
# ════════════════════════════════════════════════════════════════════

def _canonical_site_order(sites: List[CoordinationSite]) -> List[CoordinationSite]:
    # First, canonicalize internal atom order within each site
    for s in sites:
        if len(s.donor_atoms) > 1:
            combined = sorted(zip(s.donor_elements, s.donor_atoms))
            s.donor_elements = [x[0] for x in combined]
            s.donor_atoms = [x[1] for x in combined]
    ordered = sorted(sites, key=lambda s: s.site_signature())
    for i, s in enumerate(ordered):
        s.label = f"S{i+1}"
    return ordered


# ════════════════════════════════════════════════════════════════════
# Ligand canonical ordering
# ════════════════════════════════════════════════════════════════════

def _ligand_sort_key(lg: Ligand, sites: List[CoordinationSite]) -> tuple:
    lig_sites = sorted(
        s.site_signature()
        for s in sites if s.ligand_label == lg.label
    )
    return (
        lg.smiles,
        lg.dent,
        lg.charge or 0,
        tuple(lig_sites),
    )


def _canonical_ligand_order(
    ligands: List[Ligand],
    sites: List[CoordinationSite],
) -> List[Ligand]:
    # For groups of duplicate ligands (same smiles/dent/charge),
    # canonically reassign their sites: collect ALL sites for the
    # group, sort them, then distribute to ligands in order.
    from itertools import groupby
    chem_key = lambda lg: (lg.smiles, lg.dent, lg.charge or 0)
    sorted_by_chem = sorted(ligands, key=chem_key)
    for _, grp in groupby(sorted_by_chem, key=chem_key):
        grp_list = list(grp)
        if len(grp_list) <= 1:
            continue
        grp_labels = {lg.label for lg in grp_list}
        grp_sites = sorted(
            [s for s in sites if s.ligand_label in grp_labels],
            key=lambda s: s.site_signature()
        )
        # Count sites per ligand (they should all have similar counts)
        n_per = len(grp_sites) // len(grp_list) if grp_list else 1
        # Assign in order: first n_per sites → first ligand, etc.
        for gi, lg in enumerate(grp_list):
            start = gi * n_per
            end = start + n_per if gi < len(grp_list) - 1 else len(grp_sites)
            for s in grp_sites[start:end]:
                s.ligand_label = lg.label

    ordered = sorted(ligands, key=lambda lg: _ligand_sort_key(lg, sites))
    mapping = {}
    for i, lg in enumerate(ordered):
        old = lg.label
        new = f"L{i+1}"
        mapping[old] = new
        lg.label = new
    for s in sites:
        if s.ligand_label in mapping:
            s.ligand_label = mapping[s.ligand_label]
    return ordered


# ════════════════════════════════════════════════════════════════════
# Identity key generation
# ════════════════════════════════════════════════════════════════════

def _serialize_multi_content(rec: MultiMetalRecord) -> str:
    """Serialize without ID block for stable hashing."""
    from .serialize import serialize_multi
    full = serialize_multi(rec)
    # Remove [ID: ...] block
    idx = full.find("[ID:")
    return full[:idx].strip() if idx >= 0 else full


def _serialize_haptic_content(rec) -> str:
    """Serialize without ID block for stable hashing."""
    from .serialize import serialize_haptic
    full = serialize_haptic(rec)
    idx = full.find("[ID:")
    return full[:idx].strip() if idx >= 0 else full


def _generate_multi_identity(
    rec: MultiMetalRecord,
    serialized: str,
) -> MultiIdentityKeys:
    """Generate L0–L3 identity keys for a multinuclear record."""
    # Use content-only string (no ID block) for L0 to avoid circular dependency
    content = _serialize_multi_content(rec)
    l0 = _short_hash(content, 16)

    # L1: shape-level (includes shape bins)
    shape_parts = []
    for m in rec.metals:
        shape_parts.append(f"{m.element}|{m.local_shape_best}|CN{m.cn_site}")
    l1_raw = ";".join(sorted(shape_parts))
    l1 = _short_hash(l1_raw, 16)

    # L2: metal graph topology
    edge_parts = []
    for e in rec.metal_edges:
        a, b = sorted([e.m1, e.m2])
        edge_parts.append(f"{a}-{b}|{e.relation}")
    metal_parts = [f"{m.element}" for m in rec.metals]
    l2_raw = ";".join(sorted(metal_parts)) + "|" + ";".join(sorted(edge_parts))
    l2 = _short_hash(l2_raw, 16)

    # L3: connectivity only
    lig_parts = sorted(lg.smiles for lg in rec.ligands)
    l3_raw = ";".join(sorted(metal_parts)) + "|" + ";".join(lig_parts)
    l3 = _short_hash(l3_raw, 16)

    # Local IDs
    local_ids = {}
    for m in rec.metals:
        local_sites = sorted(
            s.site_signature()
            for s in rec.sites if m.label in s.target_metals
        )
        lid_raw = f"{m.element}|{m.ox_str}|{m.local_shape_best}|{local_sites}"
        local_ids[m.label] = _short_hash(lid_raw, 12)

    return MultiIdentityKeys(
        L0_GlobalState=l0,
        L1_GlobalShape=l1,
        L2_MetalGraphTopo=l2,
        L3_Connectivity=l3,
        local_ids=local_ids,
    )


def _generate_haptic_identity(
    rec: HapticRecord,
    serialized: str,
) -> HapticIdentityKeys:
    content = _serialize_haptic_content(rec)
    l0 = _short_hash(content, 16)

    m = rec.metal
    shape_raw = f"{m.element}|{m.ox_str}|{m.local_shape_best}|CNsite:{m.cn_site}"
    l1 = _short_hash(shape_raw, 16)

    site_parts = sorted(s.site_signature() for s in rec.sites)
    l2_raw = f"{m.element}|{site_parts}"
    l2 = _short_hash(l2_raw, 16)

    lig_parts = sorted(lg.smiles for lg in rec.ligands)
    l3_raw = f"{m.element}|{';'.join(lig_parts)}"
    l3 = _short_hash(l3_raw, 16)

    return HapticIdentityKeys(
        L0_HapticState=l0,
        L1_HapticShape=l1,
        L2_SiteTopo=l2,
        L3_Connectivity=l3,
    )


# ════════════════════════════════════════════════════════════════════
# Top-level canonicalize
# ════════════════════════════════════════════════════════════════════

def _canonicalize_multi_once(rec: MultiMetalRecord) -> str:
    """One pass of canonical ordering.  Returns content string."""
    metals, edges, sites, ligands = _relabel_metals(
        rec.metals, rec.metal_edges, rec.sites, rec.ligands
    )
    rec.metals = list(metals)
    rec.metal_edges = sorted(edges, key=lambda e: e.canonical_key())
    rec.sites = _canonical_site_order(sites)
    rec.ligands = _canonical_ligand_order(ligands, rec.sites)
    for m in rec.metals:
        m.local_sites = [s.label for s in rec.sites if m.label in s.target_metals]
    return _serialize_multi_content(rec)


def canonicalize_multi(rec: MultiMetalRecord) -> MultiMetalRecord:
    """Canonicalize a MultiMetalRecord to a fixed point."""
    prev = ""
    for _ in range(5):
        cur = _canonicalize_multi_once(rec)
        if cur == prev:
            break
        prev = cur

    from .serialize import serialize_multi
    serialized = serialize_multi(rec)
    rec.identity = _generate_multi_identity(rec, serialized)
    return rec


def _canonicalize_haptic_once(rec: HapticRecord) -> str:
    """One pass."""
    rec.sites = _canonical_site_order(rec.sites)
    rec.ligands = _canonical_ligand_order(rec.ligands, rec.sites)
    rec.metal.local_sites = [s.label for s in rec.sites]
    return _serialize_haptic_content(rec)


def canonicalize_haptic(rec: HapticRecord) -> HapticRecord:
    """Canonicalize a HapticRecord to a fixed point."""
    prev = ""
    for _ in range(5):
        cur = _canonicalize_haptic_once(rec)
        if cur == prev:
            break
        prev = cur

    from .serialize import serialize_haptic
    serialized = serialize_haptic(rec)
    rec.identity = _generate_haptic_identity(rec, serialized)
    return rec
