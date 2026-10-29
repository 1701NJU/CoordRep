"""
CoordRep v2-beta serialization.

Converts MultiMetalRecord / HapticRecord → canonical string.
All records are produced by this serializer — no hand-written records.
"""

from __future__ import annotations

from .core import (
    CoordinationSite, HapticRecord, Ligand, MetalCenter, MetalEdge,
    MultiMetalRecord,
)


# ════════════════════════════════════════════════════════════════════
# Multinuclear serializer
# ════════════════════════════════════════════════════════════════════

def serialize_multi(rec: MultiMetalRecord) -> str:
    """Serialize a MultiMetalRecord to a CoordRep-Multi-v2beta string."""
    parts = []

    # [Metals: …]
    metal_strs = []
    for m in rec.metals:
        d_str = f"|d:d{m.dcount}" if m.dcount is not None else ""
        ox = m.ox_str
        metal_strs.append(
            f"{m.label}=[{m.element}|ox:{ox}{d_str}"
            f"|CNsite:{m.cn_site}|CNatom:{m.cn_atom}]"
        )
    parts.append("[Metals:\n  " + ";\n  ".join(metal_strs) + "]")

    # [MetalGraph: …]
    if rec.metal_edges:
        edge_strs = []
        for e in rec.metal_edges:
            d_str = f"|dMM:{e.d_mm:.2f}" if e.d_mm is not None else ""
            edge_strs.append(
                f"{{{e.m1}--{e.m2}|relation:{e.relation}"
                f"{d_str}|MM_bond:{e.mm_bond}}}"
            )
        parts.append("[MetalGraph:\n  " + ";\n  ".join(edge_strs) + "]")

    # [LocalSphere: …]
    sphere_strs = []
    for m in rec.metals:
        shape_str = ""
        if m.local_shape_best:
            shape_str = (
                f"<ShapeBest:{m.local_shape_best}"
                f"|Class:{m.local_shape_class}"
                f"|Delta:{m.local_shape_delta}"
                f"|V:{m.local_shape_vals}>"
            )
        # Collect sites for this metal
        local_sites = [s for s in rec.sites if m.label in s.target_metals]
        site_labels = ",".join(s.label for s in local_sites)
        sphere_strs.append(
            f"{m.label}={shape_str}{{{site_labels}}}"
        )
    parts.append("[LocalSphere:\n  " + ";\n  ".join(sphere_strs) + "]")

    # [Bridges: …]
    bridge_sites = [s for s in rec.sites if s.mu > 1]
    if bridge_sites:
        br_strs = []
        for s in bridge_sites:
            targets = ",".join(s.target_metals)
            donor_str = s.donor_atoms[0] if s.donor_atoms else "?"
            elem_str = s.donor_elements[0] if s.donor_elements else "?"
            br_strs.append(
                f"{{site:{s.label}|ligand:{s.ligand_label}"
                f"|donor:{elem_str}{donor_str}"
                f"|targets:{targets}|mu:{s.mu}|mode:{s.mode}}}"
            )
        parts.append("[Bridges:\n  " + ";\n  ".join(br_strs) + "]")

    # [Ligands: …]
    lig_strs = [f"{lg.label}={lg.smiles}" for lg in rec.ligands]
    parts.append("[Ligands:\n  " + ";\n  ".join(lig_strs) + "]")

    # [ID: …]
    id_parts = [
        f"L0_GlobalState={rec.identity.L0_GlobalState}",
        f"L1_GlobalShape={rec.identity.L1_GlobalShape}",
        f"L2_MetalGraphTopo={rec.identity.L2_MetalGraphTopo}",
        f"L3_Connectivity={rec.identity.L3_Connectivity}",
    ]
    for m_label, lid in sorted(rec.identity.local_ids.items()):
        id_parts.append(f"LocalID[{m_label}]={lid}")
    parts.append("[ID:\n  " + ";\n  ".join(id_parts) + "]")

    return "\n".join(parts)


def serialize_multi_oneline(rec: MultiMetalRecord) -> str:
    """Single-line version of the multinuclear record."""
    import re
    s = serialize_multi(rec)
    return re.sub(r"\s*\n\s*", " ", s).strip()


# ════════════════════════════════════════════════════════════════════
# Haptic serializer
# ════════════════════════════════════════════════════════════════════

def serialize_haptic(rec: HapticRecord) -> str:
    """Serialize a HapticRecord to a CoordRep-Haptic-v2beta string."""
    parts = []
    m = rec.metal

    # [Metal: …]
    d_str = f"|d:d{m.dcount}" if m.dcount is not None else ""
    ox = m.ox_str
    parts.append(
        f"[Metal:\n  M1=[{m.element}|ox:{ox}{d_str}"
        f"|CNsite:{m.cn_site}|eta_sum:{m.eta_sum}]]"
    )

    # [Sites: …]
    site_strs = []
    for s in rec.sites:
        atoms_str = ",".join(
            f"{e}{a}" for e, a in zip(s.donor_elements, s.donor_atoms)
        )
        centroid_str = f"|centroid:{s.centroid_label}" if s.centroid_label else ""
        mode_str = f"|mode:{s.mode}" if s.mode else ""
        site_strs.append(
            f"{s.label}={{type:{s.site_type}|ligand:{s.ligand_label}"
            f"|atoms:{atoms_str}|eta:{s.eta}{centroid_str}{mode_str}}}"
        )
    parts.append("[Sites:\n  " + ";\n  ".join(site_strs) + "]")

    # [Hapticity: …]
    haptic_sites = [s for s in rec.sites if s.eta > 1]
    if haptic_sites:
        h_strs = []
        for s in haptic_sites:
            targets = ",".join(s.target_metals)
            h_strs.append(
                f"{{{s.label}->{targets}|eta:{s.eta}|mode:eta{s.eta}-{s.mode}}}"
            )
        parts.append("[Hapticity:\n  " + ";\n  ".join(h_strs) + "]")

    # Atom donors (η1 sites)
    atom_sites = [s for s in rec.sites if s.eta == 1]
    if atom_sites:
        a_strs = []
        for s in atom_sites:
            targets = ",".join(s.target_metals)
            elem = s.donor_elements[0] if s.donor_elements else "?"
            atom = s.donor_atoms[0] if s.donor_atoms else "?"
            a_strs.append(
                f"{{{s.label}:{s.ligand_label}:{elem}{atom}->{targets}}}"
            )
        parts.append("[AtomDonors:\n  " + ";\n  ".join(a_strs) + "]")

    # [Ligands: …]
    lig_strs = [f"{lg.label}={lg.smiles}" for lg in rec.ligands]
    parts.append("[Ligands:\n  " + ";\n  ".join(lig_strs) + "]")

    # [Geometry: …]
    geom_parts = [
        f"CNsite={m.cn_site}",
        f"eta_sum={m.eta_sum}",
    ]
    if m.local_shape_best:
        geom_parts.append(f"site_level_shape={m.local_shape_best}")
    parts.append("[Geometry:\n  " + ";\n  ".join(geom_parts) + "]")

    # [ID: …]
    parts.append(
        f"[ID:\n  L0_HapticState={rec.identity.L0_HapticState};\n"
        f"  L1_HapticShape={rec.identity.L1_HapticShape};\n"
        f"  L2_SiteTopo={rec.identity.L2_SiteTopo};\n"
        f"  L3_Connectivity={rec.identity.L3_Connectivity}]"
    )

    return "\n".join(parts)


def serialize_haptic_oneline(rec: HapticRecord) -> str:
    import re
    s = serialize_haptic(rec)
    return re.sub(r"\s*\n\s*", " ", s).strip()
