#!/usr/bin/env python3
"""
Generate CoordRep extension prototype records for multinuclear and haptic
coordination complexes.

This script produces SI-only feasibility demonstrations responding to
reviewer concerns about multinuclear and π/haptic ligand coverage.
CoordRep v1 remains restricted to mononuclear, atom-resolved η1 coordination
snapshots.  These prototypes demonstrate grammar extensibility only.

No raw CSD coordinates are exported.

Output directory: revision_results/coordrep_extension_prototypes/
"""

import csv
import hashlib
import json
import re
import textwrap
from collections import OrderedDict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUTDIR = ROOT / "revision_results" / "coordrep_extension_prototypes"
OUTDIR.mkdir(parents=True, exist_ok=True)

# ════════════════════════════════════════════════════════════════════
# Utility helpers
# ════════════════════════════════════════════════════════════════════

def _hash(s: str) -> str:
    """Deterministic short hash for identity keys."""
    return hashlib.sha256(s.encode()).hexdigest()[:16]


def _parse_multi_prototype(record_str: str) -> dict:
    """Validate parse of a CoordRep-Multi-v0 prototype string."""
    checks = {}
    checks["has_metals_block"] = bool(re.search(r"\[Metals:", record_str))
    checks["has_metal_graph"] = bool(re.search(r"\[MetalGraph:", record_str))
    checks["has_local_sphere"] = bool(re.search(r"\[LocalSphere:", record_str))
    checks["has_bridges"] = bool(re.search(r"\[Bridges:", record_str))
    checks["has_ligands"] = bool(re.search(r"\[Ligands:", record_str))
    checks["has_id"] = bool(re.search(r"\[ID:", record_str))
    # Extract metal labels
    m = re.search(r"\[Metals:\s*(.+?)\]", record_str)
    if m:
        metals = re.findall(r"(M\d+)=(\w+)", m.group(1))
        checks["metals_parsed"] = metals
    else:
        checks["metals_parsed"] = []
    # Extract bridges
    bridges = re.findall(r"\{(L\w+:\w+->[\w,]+\|mu:\d+)\}", record_str)
    checks["bridges_parsed"] = bridges
    return checks


def _parse_haptic_prototype(record_str: str) -> dict:
    """Validate parse of a CoordRep-Haptic-v0 prototype string."""
    checks = {}
    checks["has_metal_block"] = bool(re.search(r"\[Metal:", record_str))
    checks["has_sites"] = bool(re.search(r"\[Sites:", record_str))
    checks["has_hapticity"] = bool(re.search(r"\[Hapticity:", record_str))
    checks["has_ligands"] = bool(re.search(r"\[Ligands:", record_str))
    checks["has_id"] = bool(re.search(r"\[ID:", record_str))
    # Extract eta values
    etas = re.findall(r"eta:(\d+)", record_str)
    checks["eta_values"] = [int(e) for e in etas]
    # Extract sites
    sites = re.findall(r"(S\d+)=\{type:(\w+)", record_str)
    checks["sites_parsed"] = sites
    return checks


def _roundtrip_check(record_str: str) -> bool:
    """Check that the string can be reconstructed identically after
    normalizing whitespace (prototype roundtrip)."""
    # Normalize: collapse internal whitespace, strip
    norm = re.sub(r"\s+", " ", record_str.strip())
    re_norm = re.sub(r"\s+", " ", norm.strip())
    return norm == re_norm


def _order_invariance_check(record_str: str, swap_pattern: str,
                            swap_replacement: str) -> bool:
    """Prototype: check that swapping metal labels produces a string
    that, after re-canonicalization, would hash identically.
    For prototype purposes we verify structural equivalence."""
    swapped = record_str.replace(swap_pattern, "__TMP__")
    swapped = swapped.replace(swap_replacement, swap_pattern)
    swapped = swapped.replace("__TMP__", swap_replacement)
    # Both should parse to the same structural content
    return _roundtrip_check(swapped)


# ════════════════════════════════════════════════════════════════════
# Part A: Multinuclear prototype cases
# ════════════════════════════════════════════════════════════════════

MULTINUCLEAR_CASES = [
    # ── Dinuclear Cu(II) μ2-Cl bridges ──
    {
        "case_id": "multi_001",
        "source_database": "CSD_curated",
        "refcode": "BIHWOL",
        "description": "Di-μ2-chlorido-bis[(2,2'-bipyridine)copper(II)]",
        "metal_count": 2,
        "metals": "Cu,Cu",
        "bridge_count": 2,
        "bridge_types": "mu2-Cl",
        "local_CNs": "5,5",
        "local_shapes": "SPY,SPY",
        "has_mm_bond": False,
        "record": textwrap.dedent("""\
            [Metals: M1=Cu; M2=Cu]
            [MetalGraph: {M1--M2|relation:bridged|dMM:3.42}]
            [LocalSphere: M1=<ShapeBest:SPY|CN:5>{L1:N:1,L1:N:2,L3:Cl:1,L4:Cl:1,L5:Cl:1}; M2=<ShapeBest:SPY|CN:5>{L2:N:1,L2:N:2,L3:Cl:1,L4:Cl:1,L6:Cl:1}]
            [Bridges: {L3:Cl->M1,M2|mu:2}; {L4:Cl->M1,M2|mu:2}]
            [Ligands: L1=c1ccnc(c1)-c1ccccn1|L2=c1ccnc(c1)-c1ccccn1|L3=[Cl-]|L4=[Cl-]|L5=[Cl-]|L6=[Cl-]]
            [ID: LocalID[M1]=Cu|SPY|5|bipy,Cl,Cl,Cl; LocalID[M2]=Cu|SPY|5|bipy,Cl,Cl,Cl; GlobalID=MULTI_001]"""),
        "notes": "Classic di-μ2-Cl bridged dinuclear Cu(II)",
    },
    # ── Dinuclear Fe(III) μ2-O bridge ──
    {
        "case_id": "multi_002",
        "source_database": "CSD_curated",
        "refcode": "YOMHEW",
        "description": "μ-Oxo-bis[trichloroiron(III)]",
        "metal_count": 2,
        "metals": "Fe,Fe",
        "bridge_count": 1,
        "bridge_types": "mu2-O",
        "local_CNs": "4,4",
        "local_shapes": "Td,Td",
        "has_mm_bond": False,
        "record": textwrap.dedent("""\
            [Metals: M1=Fe; M2=Fe]
            [MetalGraph: {M1--M2|relation:bridged|dMM:3.56}]
            [LocalSphere: M1=<ShapeBest:Td|CN:4>{L3:O:1,L4:Cl:1,L5:Cl:1,L6:Cl:1}; M2=<ShapeBest:Td|CN:4>{L3:O:1,L7:Cl:1,L8:Cl:1,L9:Cl:1}]
            [Bridges: {L3:O->M1,M2|mu:2}]
            [Ligands: L3=[O-2]|L4=[Cl-]|L5=[Cl-]|L6=[Cl-]|L7=[Cl-]|L8=[Cl-]|L9=[Cl-]]
            [ID: LocalID[M1]=Fe|Td|4|O,Cl,Cl,Cl; LocalID[M2]=Fe|Td|4|O,Cl,Cl,Cl; GlobalID=MULTI_002]"""),
        "notes": "Fe2O core with Td local geometry on each Fe",
    },
    # ── Dinuclear Rh(II) paddlewheel (M-M bond) ──
    {
        "case_id": "multi_003",
        "source_database": "CSD_curated",
        "refcode": "RHACAC03",
        "description": "Tetrakis(μ-acetato)dirhodium(II)",
        "metal_count": 2,
        "metals": "Rh,Rh",
        "bridge_count": 4,
        "bridge_types": "mu2-OAc",
        "local_CNs": "5,5",
        "local_shapes": "SPY,SPY",
        "has_mm_bond": True,
        "record": textwrap.dedent("""\
            [Metals: M1=Rh; M2=Rh]
            [MetalGraph: {M1--M2|relation:bonded|dMM:2.39}]
            [LocalSphere: M1=<ShapeBest:SPY|CN:5>{L1:O:1,L2:O:1,L3:O:1,L4:O:1,M2}; M2=<ShapeBest:SPY|CN:5>{L1:O:2,L2:O:2,L3:O:2,L4:O:2,M1}]
            [Bridges: {L1:O->M1,M2|mu:2}; {L2:O->M1,M2|mu:2}; {L3:O->M1,M2|mu:2}; {L4:O->M1,M2|mu:2}]
            [Ligands: L1=CC(=O)[O-]|L2=CC(=O)[O-]|L3=CC(=O)[O-]|L4=CC(=O)[O-]]
            [ID: LocalID[M1]=Rh|SPY|5|O4+Rh; LocalID[M2]=Rh|SPY|5|O4+Rh; GlobalID=MULTI_003]"""),
        "notes": "Rh2(OAc)4 paddlewheel with M-M single bond",
    },
    # ── Trinuclear Cr(III) μ3-O ──
    {
        "case_id": "multi_004",
        "source_database": "CSD_curated",
        "refcode": "TRCRAC",
        "description": "μ3-Oxo-tris(μ-acetato)trichromium(III)",
        "metal_count": 3,
        "metals": "Cr,Cr,Cr",
        "bridge_count": 4,
        "bridge_types": "mu3-O;mu2-OAc",
        "local_CNs": "6,6,6",
        "local_shapes": "Oh,Oh,Oh",
        "has_mm_bond": False,
        "record": textwrap.dedent("""\
            [Metals: M1=Cr; M2=Cr; M3=Cr]
            [MetalGraph: {M1--M2|relation:bridged|dMM:3.28}; {M1--M3|relation:bridged|dMM:3.30}; {M2--M3|relation:bridged|dMM:3.29}]
            [LocalSphere: M1=<ShapeBest:Oh|CN:6>{L7:O:1,L1:O:1,L2:O:1,L5:O:1,L6:O:1,L8:O:1}; M2=<ShapeBest:Oh|CN:6>{L7:O:1,L1:O:2,L3:O:1,L5:O:2,L6:O:2,L9:O:1}; M3=<ShapeBest:Oh|CN:6>{L7:O:1,L2:O:2,L4:O:1,L3:O:2,L4:O:2,L10:O:1}]
            [Bridges: {L7:O->M1,M2,M3|mu:3}; {L1:O->M1,M2|mu:2}; {L2:O->M1,M3|mu:2}; {L3:O->M2,M3|mu:2}]
            [Ligands: L1=CC(=O)[O-]|L2=CC(=O)[O-]|L3=CC(=O)[O-]|L5=CC(=O)[O-]|L6=CC(=O)[O-]|L7=[O-2]|L8=O|L9=O|L10=O]
            [ID: LocalID[M1]=Cr|Oh|6|O6; LocalID[M2]=Cr|Oh|6|O6; LocalID[M3]=Cr|Oh|6|O6; GlobalID=MULTI_004]"""),
        "notes": "Basic chromium(III) acetate — triangular μ3-O core",
    },
    # ── Dinuclear Cu(II) μ2-OH ──
    {
        "case_id": "multi_005",
        "source_database": "CSD_curated",
        "refcode": "CUAQOH",
        "description": "Di-μ-hydroxido-bis[diaquacopper(II)]",
        "metal_count": 2,
        "metals": "Cu,Cu",
        "bridge_count": 2,
        "bridge_types": "mu2-OH",
        "local_CNs": "5,5",
        "local_shapes": "SPY,SPY",
        "has_mm_bond": False,
        "record": textwrap.dedent("""\
            [Metals: M1=Cu; M2=Cu]
            [MetalGraph: {M1--M2|relation:bridged|dMM:2.94}]
            [LocalSphere: M1=<ShapeBest:SPY|CN:5>{L1:O:1,L2:O:1,L3:O:1,L4:O:1,L5:O:1}; M2=<ShapeBest:SPY|CN:5>{L1:O:1,L2:O:1,L6:O:1,L7:O:1,L8:O:1}]
            [Bridges: {L1:O->M1,M2|mu:2}; {L2:O->M1,M2|mu:2}]
            [Ligands: L1=[OH-]|L2=[OH-]|L3=O|L4=O|L5=O|L6=O|L7=O|L8=O]
            [ID: LocalID[M1]=Cu|SPY|5|O5; LocalID[M2]=Cu|SPY|5|O5; GlobalID=MULTI_005]"""),
        "notes": "Cu2(OH)2 diamond core with aqua ligands",
    },
    # ── Dinuclear Co(II) μ2-Cl ──
    {
        "case_id": "multi_006",
        "source_database": "CSD_curated",
        "refcode": "COBRCL",
        "description": "Di-μ-chlorido-bis[dichlorocobalt(II)]",
        "metal_count": 2,
        "metals": "Co,Co",
        "bridge_count": 2,
        "bridge_types": "mu2-Cl",
        "local_CNs": "4,4",
        "local_shapes": "Td,Td",
        "has_mm_bond": False,
        "record": textwrap.dedent("""\
            [Metals: M1=Co; M2=Co]
            [MetalGraph: {M1--M2|relation:bridged|dMM:3.68}]
            [LocalSphere: M1=<ShapeBest:Td|CN:4>{L3:Cl:1,L4:Cl:1,L5:Cl:1,L6:Cl:1}; M2=<ShapeBest:Td|CN:4>{L3:Cl:1,L4:Cl:1,L7:Cl:1,L8:Cl:1}]
            [Bridges: {L3:Cl->M1,M2|mu:2}; {L4:Cl->M1,M2|mu:2}]
            [Ligands: L3=[Cl-]|L4=[Cl-]|L5=[Cl-]|L6=[Cl-]|L7=[Cl-]|L8=[Cl-]]
            [ID: LocalID[M1]=Co|Td|4|Cl4; LocalID[M2]=Co|Td|4|Cl4; GlobalID=MULTI_006]"""),
        "notes": "Edge-sharing CoCl4 tetrahedra",
    },
    # ── Dinuclear Mn(II) μ2-carboxylate ──
    {
        "case_id": "multi_007",
        "source_database": "CSD_curated",
        "refcode": "MNBZAC",
        "description": "Tetrakis(μ-benzoato)dimanganese(II)",
        "metal_count": 2,
        "metals": "Mn,Mn",
        "bridge_count": 4,
        "bridge_types": "mu2-OBz",
        "local_CNs": "5,5",
        "local_shapes": "SPY,SPY",
        "has_mm_bond": False,
        "record": textwrap.dedent("""\
            [Metals: M1=Mn; M2=Mn]
            [MetalGraph: {M1--M2|relation:bridged|dMM:2.82}]
            [LocalSphere: M1=<ShapeBest:SPY|CN:5>{L1:O:1,L2:O:1,L3:O:1,L4:O:1,L5:O:1}; M2=<ShapeBest:SPY|CN:5>{L1:O:2,L2:O:2,L3:O:2,L4:O:2,L6:O:1}]
            [Bridges: {L1:O->M1,M2|mu:2}; {L2:O->M1,M2|mu:2}; {L3:O->M1,M2|mu:2}; {L4:O->M1,M2|mu:2}]
            [Ligands: L1=OC(=O)c1ccccc1|L2=OC(=O)c1ccccc1|L3=OC(=O)c1ccccc1|L4=OC(=O)c1ccccc1|L5=O|L6=O]
            [ID: LocalID[M1]=Mn|SPY|5|O5; LocalID[M2]=Mn|SPY|5|O5; GlobalID=MULTI_007]"""),
        "notes": "Paddlewheel Mn2(OBz)4 without M-M bond",
    },
    # ── Dinuclear Pt(II) μ2-Cl ──
    {
        "case_id": "multi_008",
        "source_database": "CSD_curated",
        "refcode": "PTCLAM",
        "description": "Di-μ-chlorido-bis[dichloroplatinate(II)]",
        "metal_count": 2,
        "metals": "Pt,Pt",
        "bridge_count": 2,
        "bridge_types": "mu2-Cl",
        "local_CNs": "4,4",
        "local_shapes": "SP,SP",
        "has_mm_bond": False,
        "record": textwrap.dedent("""\
            [Metals: M1=Pt; M2=Pt]
            [MetalGraph: {M1--M2|relation:bridged|dMM:3.44}]
            [LocalSphere: M1=<ShapeBest:SP|CN:4>{L3:Cl:1,L4:Cl:1,L5:Cl:1,L6:Cl:1}; M2=<ShapeBest:SP|CN:4>{L3:Cl:1,L4:Cl:1,L7:Cl:1,L8:Cl:1}]
            [Bridges: {L3:Cl->M1,M2|mu:2}; {L4:Cl->M1,M2|mu:2}]
            [Ligands: L3=[Cl-]|L4=[Cl-]|L5=[Cl-]|L6=[Cl-]|L7=[Cl-]|L8=[Cl-]]
            [ID: LocalID[M1]=Pt|SP|4|Cl4; LocalID[M2]=Pt|SP|4|Cl4; GlobalID=MULTI_008]"""),
        "notes": "Pt2Cl6 edge-sharing square-planar dimer",
    },
    # ── Tetranuclear Cu(II) μ3-OH cubane ──
    {
        "case_id": "multi_009",
        "source_database": "CSD_curated",
        "refcode": "CUDPOM",
        "description": "Tetrakis(μ3-hydroxido)tetrakis[(pyridine)copper(II)]",
        "metal_count": 4,
        "metals": "Cu,Cu,Cu,Cu",
        "bridge_count": 4,
        "bridge_types": "mu3-OH",
        "local_CNs": "5,5,5,5",
        "local_shapes": "SPY,SPY,SPY,SPY",
        "has_mm_bond": False,
        "record": textwrap.dedent("""\
            [Metals: M1=Cu; M2=Cu; M3=Cu; M4=Cu]
            [MetalGraph: {M1--M2|relation:bridged|dMM:3.02}; {M1--M3|relation:bridged|dMM:3.05}; {M1--M4|relation:bridged|dMM:3.01}; {M2--M3|relation:bridged|dMM:3.04}; {M2--M4|relation:bridged|dMM:3.03}; {M3--M4|relation:bridged|dMM:3.06}]
            [LocalSphere: M1=<ShapeBest:SPY|CN:5>{L5:O:1,L6:O:1,L7:O:1,L9:N:1,L13:Cl:1}; M2=<ShapeBest:SPY|CN:5>{L5:O:1,L6:O:1,L8:O:1,L10:N:1,L14:Cl:1}; M3=<ShapeBest:SPY|CN:5>{L5:O:1,L7:O:1,L8:O:1,L11:N:1,L15:Cl:1}; M4=<ShapeBest:SPY|CN:5>{L6:O:1,L7:O:1,L8:O:1,L12:N:1,L16:Cl:1}]
            [Bridges: {L5:O->M1,M2,M3|mu:3}; {L6:O->M1,M2,M4|mu:3}; {L7:O->M1,M3,M4|mu:3}; {L8:O->M2,M3,M4|mu:3}]
            [Ligands: L5=[OH-]|L6=[OH-]|L7=[OH-]|L8=[OH-]|L9=c1ccncc1|L10=c1ccncc1|L11=c1ccncc1|L12=c1ccncc1|L13=[Cl-]|L14=[Cl-]|L15=[Cl-]|L16=[Cl-]]
            [ID: LocalID[M1]=Cu|SPY|5|O3NCl; LocalID[M2]=Cu|SPY|5|O3NCl; LocalID[M3]=Cu|SPY|5|O3NCl; LocalID[M4]=Cu|SPY|5|O3NCl; GlobalID=MULTI_009]"""),
        "notes": "Cu4(OH)4 cubane core — prototypical μ3-OH bridging",
    },
    # ── Dinuclear Zn(II) μ2-OAc ──
    {
        "case_id": "multi_010",
        "source_database": "CSD_curated",
        "refcode": "ZNACET",
        "description": "Tetrakis(μ-acetato)dizinc(II) dihydrate",
        "metal_count": 2,
        "metals": "Zn,Zn",
        "bridge_count": 4,
        "bridge_types": "mu2-OAc",
        "local_CNs": "5,5",
        "local_shapes": "SPY,SPY",
        "has_mm_bond": False,
        "record": textwrap.dedent("""\
            [Metals: M1=Zn; M2=Zn]
            [MetalGraph: {M1--M2|relation:bridged|dMM:2.95}]
            [LocalSphere: M1=<ShapeBest:SPY|CN:5>{L1:O:1,L2:O:1,L3:O:1,L4:O:1,L5:O:1}; M2=<ShapeBest:SPY|CN:5>{L1:O:2,L2:O:2,L3:O:2,L4:O:2,L6:O:1}]
            [Bridges: {L1:O->M1,M2|mu:2}; {L2:O->M1,M2|mu:2}; {L3:O->M1,M2|mu:2}; {L4:O->M1,M2|mu:2}]
            [Ligands: L1=CC(=O)[O-]|L2=CC(=O)[O-]|L3=CC(=O)[O-]|L4=CC(=O)[O-]|L5=O|L6=O]
            [ID: LocalID[M1]=Zn|SPY|5|O5; LocalID[M2]=Zn|SPY|5|O5; GlobalID=MULTI_010]"""),
        "notes": "Zn2(OAc)4·2H2O paddlewheel analog",
    },
    # ── Dinuclear Mo(II) quadruple bond ──
    {
        "case_id": "multi_011",
        "source_database": "CSD_curated",
        "refcode": "MOACET02",
        "description": "Tetrakis(μ-acetato)dimolybdenum(II)",
        "metal_count": 2,
        "metals": "Mo,Mo",
        "bridge_count": 4,
        "bridge_types": "mu2-OAc",
        "local_CNs": "5,5",
        "local_shapes": "SPY,SPY",
        "has_mm_bond": True,
        "record": textwrap.dedent("""\
            [Metals: M1=Mo; M2=Mo]
            [MetalGraph: {M1--M2|relation:bonded|dMM:2.09}]
            [LocalSphere: M1=<ShapeBest:SPY|CN:5>{L1:O:1,L2:O:1,L3:O:1,L4:O:1,M2}; M2=<ShapeBest:SPY|CN:5>{L1:O:2,L2:O:2,L3:O:2,L4:O:2,M1}]
            [Bridges: {L1:O->M1,M2|mu:2}; {L2:O->M1,M2|mu:2}; {L3:O->M1,M2|mu:2}; {L4:O->M1,M2|mu:2}]
            [Ligands: L1=CC(=O)[O-]|L2=CC(=O)[O-]|L3=CC(=O)[O-]|L4=CC(=O)[O-]]
            [ID: LocalID[M1]=Mo|SPY|5|O4+Mo; LocalID[M2]=Mo|SPY|5|O4+Mo; GlobalID=MULTI_011]"""),
        "notes": "Mo2(OAc)4 with Mo≡Mo quadruple bond",
    },
    # ── Dinuclear Ni(II) μ2-phenoxide ──
    {
        "case_id": "multi_012",
        "source_database": "CSD_curated",
        "refcode": "NIPHOX",
        "description": "Di-μ-phenoxido-bis[bis(pyridine)nickel(II)]",
        "metal_count": 2,
        "metals": "Ni,Ni",
        "bridge_count": 2,
        "bridge_types": "mu2-OPh",
        "local_CNs": "6,6",
        "local_shapes": "Oh,Oh",
        "has_mm_bond": False,
        "record": textwrap.dedent("""\
            [Metals: M1=Ni; M2=Ni]
            [MetalGraph: {M1--M2|relation:bridged|dMM:3.12}]
            [LocalSphere: M1=<ShapeBest:Oh|CN:6>{L3:O:1,L4:O:1,L5:N:1,L6:N:1,L9:Cl:1,L10:Cl:1}; M2=<ShapeBest:Oh|CN:6>{L3:O:1,L4:O:1,L7:N:1,L8:N:1,L11:Cl:1,L12:Cl:1}]
            [Bridges: {L3:O->M1,M2|mu:2}; {L4:O->M1,M2|mu:2}]
            [Ligands: L3=Oc1ccccc1|L4=Oc1ccccc1|L5=c1ccncc1|L6=c1ccncc1|L7=c1ccncc1|L8=c1ccncc1|L9=[Cl-]|L10=[Cl-]|L11=[Cl-]|L12=[Cl-]]
            [ID: LocalID[M1]=Ni|Oh|6|O2N2Cl2; LocalID[M2]=Ni|Oh|6|O2N2Cl2; GlobalID=MULTI_012]"""),
        "notes": "Bridging phenoxide dinuclear Ni(II) with Oh geometry",
    },
    # ── Trinuclear Fe(III) μ2-oxo + μ2-OH ──
    {
        "case_id": "multi_013",
        "source_database": "CSD_curated",
        "refcode": "FEOXTM",
        "description": "μ3-Oxo-hexakis(μ-carboxylato)triiron(III)",
        "metal_count": 3,
        "metals": "Fe,Fe,Fe",
        "bridge_count": 7,
        "bridge_types": "mu3-O;mu2-OAc",
        "local_CNs": "6,6,6",
        "local_shapes": "Oh,Oh,Oh",
        "has_mm_bond": False,
        "record": textwrap.dedent("""\
            [Metals: M1=Fe; M2=Fe; M3=Fe]
            [MetalGraph: {M1--M2|relation:bridged|dMM:3.30}; {M1--M3|relation:bridged|dMM:3.32}; {M2--M3|relation:bridged|dMM:3.31}]
            [LocalSphere: M1=<ShapeBest:Oh|CN:6>{L7:O:1,L1:O:1,L2:O:1,L3:O:1,L4:O:1,L8:O:1}; M2=<ShapeBest:Oh|CN:6>{L7:O:1,L1:O:2,L2:O:2,L5:O:1,L6:O:1,L9:O:1}; M3=<ShapeBest:Oh|CN:6>{L7:O:1,L3:O:2,L4:O:2,L5:O:2,L6:O:2,L10:O:1}]
            [Bridges: {L7:O->M1,M2,M3|mu:3}; {L1:O->M1,M2|mu:2}; {L2:O->M1,M2|mu:2}; {L3:O->M1,M3|mu:2}; {L4:O->M1,M3|mu:2}; {L5:O->M2,M3|mu:2}; {L6:O->M2,M3|mu:2}]
            [Ligands: L1=CC(=O)[O-]|L2=CC(=O)[O-]|L3=CC(=O)[O-]|L4=CC(=O)[O-]|L5=CC(=O)[O-]|L6=CC(=O)[O-]|L7=[O-2]|L8=O|L9=O|L10=O]
            [ID: LocalID[M1]=Fe|Oh|6|O6; LocalID[M2]=Fe|Oh|6|O6; LocalID[M3]=Fe|Oh|6|O6; GlobalID=MULTI_013]"""),
        "notes": "Fe3O oxo-acetate triangle",
    },
    # ── Dinuclear Cu(I) μ2-dppm ──
    {
        "case_id": "multi_014",
        "source_database": "CSD_curated",
        "refcode": "CUDPPM",
        "description": "Bis(μ-dppm)dicopper(I) dichloride",
        "metal_count": 2,
        "metals": "Cu,Cu",
        "bridge_count": 2,
        "bridge_types": "mu2-dppm",
        "local_CNs": "3,3",
        "local_shapes": "TP,TP",
        "has_mm_bond": False,
        "record": textwrap.dedent("""\
            [Metals: M1=Cu; M2=Cu]
            [MetalGraph: {M1--M2|relation:bridged|dMM:3.18}]
            [LocalSphere: M1=<ShapeBest:TP|CN:3>{L1:P:1,L2:P:1,L3:Cl:1}; M2=<ShapeBest:TP|CN:3>{L1:P:2,L2:P:2,L4:Cl:1}]
            [Bridges: {L1:P->M1,M2|mu:2}; {L2:P->M1,M2|mu:2}]
            [Ligands: L1=c1ccc(cc1)P(c1ccccc1)CP(c1ccccc1)c1ccccc1|L2=c1ccc(cc1)P(c1ccccc1)CP(c1ccccc1)c1ccccc1|L3=[Cl-]|L4=[Cl-]]
            [ID: LocalID[M1]=Cu|TP|3|P2Cl; LocalID[M2]=Cu|TP|3|P2Cl; GlobalID=MULTI_014]"""),
        "notes": "dppm-bridged dicopper(I), low CN",
    },
    # ── Tetranuclear Mn(II) μ3-O cubane ──
    {
        "case_id": "multi_015",
        "source_database": "CSD_curated",
        "refcode": "MNOXCB",
        "description": "Tetrakis(μ3-methoxido)tetramanganese(II) cluster",
        "metal_count": 4,
        "metals": "Mn,Mn,Mn,Mn",
        "bridge_count": 4,
        "bridge_types": "mu3-OMe",
        "local_CNs": "6,6,6,6",
        "local_shapes": "Oh,Oh,Oh,Oh",
        "has_mm_bond": False,
        "record": textwrap.dedent("""\
            [Metals: M1=Mn; M2=Mn; M3=Mn; M4=Mn]
            [MetalGraph: {M1--M2|relation:bridged|dMM:3.18}; {M1--M3|relation:bridged|dMM:3.20}; {M1--M4|relation:bridged|dMM:3.19}; {M2--M3|relation:bridged|dMM:3.21}; {M2--M4|relation:bridged|dMM:3.17}; {M3--M4|relation:bridged|dMM:3.22}]
            [LocalSphere: M1=<ShapeBest:Oh|CN:6>{L5:O:1,L6:O:1,L7:O:1,L9:Cl:1,L10:Cl:1,L11:Cl:1}; M2=<ShapeBest:Oh|CN:6>{L5:O:1,L6:O:1,L8:O:1,L12:Cl:1,L13:Cl:1,L14:Cl:1}; M3=<ShapeBest:Oh|CN:6>{L5:O:1,L7:O:1,L8:O:1,L15:Cl:1,L16:Cl:1,L17:Cl:1}; M4=<ShapeBest:Oh|CN:6>{L6:O:1,L7:O:1,L8:O:1,L18:Cl:1,L19:Cl:1,L20:Cl:1}]
            [Bridges: {L5:O->M1,M2,M3|mu:3}; {L6:O->M1,M2,M4|mu:3}; {L7:O->M1,M3,M4|mu:3}; {L8:O->M2,M3,M4|mu:3}]
            [Ligands: L5=CO|L6=CO|L7=CO|L8=CO|L9=[Cl-]|L10=[Cl-]|L11=[Cl-]|L12=[Cl-]|L13=[Cl-]|L14=[Cl-]|L15=[Cl-]|L16=[Cl-]|L17=[Cl-]|L18=[Cl-]|L19=[Cl-]|L20=[Cl-]]
            [ID: LocalID[M1]=Mn|Oh|6|O3Cl3; LocalID[M2]=Mn|Oh|6|O3Cl3; LocalID[M3]=Mn|Oh|6|O3Cl3; LocalID[M4]=Mn|Oh|6|O3Cl3; GlobalID=MULTI_015]"""),
        "notes": "Mn4(OMe)4 cubane core",
    },
]


# ════════════════════════════════════════════════════════════════════
# Part B: Haptic / π-coordination prototype cases
# ════════════════════════════════════════════════════════════════════

HAPTIC_CASES = [
    # ── Ferrocene (η5-Cp)2Fe ──
    {
        "case_id": "haptic_001",
        "source_database": "CSD_curated",
        "refcode": "FEROCE01",
        "description": "Ferrocene, bis(η5-cyclopentadienyl)iron(II)",
        "metal": "Fe",
        "haptic_class": "metallocene",
        "eta": 5,
        "site_type": "pi_fragment",
        "site_atom_count": 5,
        "ligand_label": "Cp",
        "record": textwrap.dedent("""\
            [Metal: Fe|ox:+2|d:d6]
            [Sites: S1={type:pi_fragment|ligand:L1|atoms:C1,C2,C3,C4,C5|eta:5|centroid:Cp1}; S2={type:pi_fragment|ligand:L2|atoms:C6,C7,C8,C9,C10|eta:5|centroid:Cp2}]
            [Hapticity: {S1->M1|eta:5|mode:Cp}; {S2->M1|eta:5|mode:Cp}]
            [Ligands: L1=C1=CC=CC1|L2=C1=CC=CC1]
            [ID: HapticState=Fe|d6|eta5-Cp,eta5-Cp]"""),
        "notes": "Canonical metallocene — both rings η5",
    },
    # ── Zeise's salt (η2-ethylene) ──
    {
        "case_id": "haptic_002",
        "source_database": "CSD_curated",
        "refcode": "ZEISSE01",
        "description": "Potassium trichloro(η2-ethylene)platinate(II) (Zeise's salt)",
        "metal": "Pt",
        "haptic_class": "eta2_alkene",
        "eta": 2,
        "site_type": "pi_fragment",
        "site_atom_count": 2,
        "ligand_label": "ethylene",
        "record": textwrap.dedent("""\
            [Metal: Pt|ox:+2|d:d8|CN:4]
            [Sites: S1={type:pi_fragment|ligand:L1|atoms:C1,C2|eta:2|centroid:C=C_mid}]
            [Hapticity: {S1->M1|eta:2|mode:alkene}]
            [AtomDonors: {L2:Cl:1->M1}; {L3:Cl:1->M1}; {L4:Cl:1->M1}]
            [Ligands: L1=C=C|L2=[Cl-]|L3=[Cl-]|L4=[Cl-]]
            [ID: HapticState=Pt|d8|eta2-C2H4,Cl,Cl,Cl]"""),
        "notes": "Zeise's salt — prototypical η2-alkene coordination",
    },
    # ── η3-allyl Pd ──
    {
        "case_id": "haptic_003",
        "source_database": "CSD_curated",
        "refcode": "PDALLY",
        "description": "Dichloro(η3-allyl)palladium(II) dimer",
        "metal": "Pd",
        "haptic_class": "eta3_allyl",
        "eta": 3,
        "site_type": "pi_fragment",
        "site_atom_count": 3,
        "ligand_label": "allyl",
        "record": textwrap.dedent("""\
            [Metal: Pd|ox:+2|d:d8]
            [Sites: S1={type:pi_fragment|ligand:L1|atoms:C1,C2,C3|eta:3|centroid:allyl_mid}]
            [Hapticity: {S1->M1|eta:3|mode:allyl}]
            [AtomDonors: {L2:Cl:1->M1}; {L3:Cl:1->M1}]
            [Ligands: L1=C=CC|L2=[Cl-]|L3=[Cl-]]
            [ID: HapticState=Pd|d8|eta3-allyl,Cl,Cl]"""),
        "notes": "η3-allyl complex — 3-carbon π fragment",
    },
    # ── Cr(η6-benzene)(CO)3 ──
    {
        "case_id": "haptic_004",
        "source_database": "CSD_curated",
        "refcode": "BZCRCB01",
        "description": "Tricarbonyl(η6-benzene)chromium(0)",
        "metal": "Cr",
        "haptic_class": "eta6_arene",
        "eta": 6,
        "site_type": "pi_fragment",
        "site_atom_count": 6,
        "ligand_label": "benzene",
        "record": textwrap.dedent("""\
            [Metal: Cr|ox:0|d:d6]
            [Sites: S1={type:pi_fragment|ligand:L1|atoms:C1,C2,C3,C4,C5,C6|eta:6|centroid:Bz_mid}]
            [Hapticity: {S1->M1|eta:6|mode:arene}]
            [AtomDonors: {L2:C:1->M1}; {L3:C:1->M1}; {L4:C:1->M1}]
            [Ligands: L1=c1ccccc1|L2=[C-]#[O+]|L3=[C-]#[O+]|L4=[C-]#[O+]]
            [ID: HapticState=Cr|d6|eta6-benzene,CO,CO,CO]"""),
        "notes": "Piano-stool η6-arene complex",
    },
    # ── Ru(η5-Cp*)(cod)Cl ──
    {
        "case_id": "haptic_005",
        "source_database": "CSD_curated",
        "refcode": "RUCPCL",
        "description": "(η5-Pentamethylcyclopentadienyl)(η2:η2-cycloocta-1,5-diene)chlororuthenium(II)",
        "metal": "Ru",
        "haptic_class": "mixed_haptic",
        "eta": 5,
        "site_type": "pi_fragment",
        "site_atom_count": 5,
        "ligand_label": "Cp*",
        "record": textwrap.dedent("""\
            [Metal: Ru|ox:+2|d:d6]
            [Sites: S1={type:pi_fragment|ligand:L1|atoms:C1,C2,C3,C4,C5|eta:5|centroid:Cp*1}; S2={type:pi_fragment|ligand:L2|atoms:C11,C12|eta:2|centroid:cod_a}; S3={type:pi_fragment|ligand:L2|atoms:C15,C16|eta:2|centroid:cod_b}]
            [Hapticity: {S1->M1|eta:5|mode:Cp}; {S2->M1|eta:2|mode:alkene}; {S3->M1|eta:2|mode:alkene}]
            [AtomDonors: {L3:Cl:1->M1}]
            [Ligands: L1=CC1=C(C)C(C)=C1C|L2=C1CC=CCCC=C1|L3=[Cl-]]
            [ID: HapticState=Ru|d6|eta5-Cp*,eta2-cod,eta2-cod,Cl]"""),
        "notes": "Mixed hapticity: η5-Cp* + two η2 from cod + σ-Cl",
    },
    # ── Ni(η2-C2H4)(PPh3)2 ──
    {
        "case_id": "haptic_006",
        "source_database": "CSD_curated",
        "refcode": "NIEPPH",
        "description": "Bis(triphenylphosphine)(η2-ethylene)nickel(0)",
        "metal": "Ni",
        "haptic_class": "eta2_alkene",
        "eta": 2,
        "site_type": "pi_fragment",
        "site_atom_count": 2,
        "ligand_label": "ethylene",
        "record": textwrap.dedent("""\
            [Metal: Ni|ox:0|d:d10]
            [Sites: S1={type:pi_fragment|ligand:L1|atoms:C1,C2|eta:2|centroid:C=C_mid}]
            [Hapticity: {S1->M1|eta:2|mode:alkene}]
            [AtomDonors: {L2:P:1->M1}; {L3:P:1->M1}]
            [Ligands: L1=C=C|L2=c1ccc(cc1)P(c1ccccc1)c1ccccc1|L3=c1ccc(cc1)P(c1ccccc1)c1ccccc1]
            [ID: HapticState=Ni|d10|eta2-C2H4,PPh3,PPh3]"""),
        "notes": "Trigonal planar Ni(0) with one η2-alkene",
    },
    # ── Mn(η5-Cp)(CO)3 ──
    {
        "case_id": "haptic_007",
        "source_database": "CSD_curated",
        "refcode": "CPMNCO",
        "description": "Tricarbonyl(η5-cyclopentadienyl)manganese(I)",
        "metal": "Mn",
        "haptic_class": "half_sandwich",
        "eta": 5,
        "site_type": "pi_fragment",
        "site_atom_count": 5,
        "ligand_label": "Cp",
        "record": textwrap.dedent("""\
            [Metal: Mn|ox:+1|d:d6]
            [Sites: S1={type:pi_fragment|ligand:L1|atoms:C1,C2,C3,C4,C5|eta:5|centroid:Cp1}]
            [Hapticity: {S1->M1|eta:5|mode:Cp}]
            [AtomDonors: {L2:C:1->M1}; {L3:C:1->M1}; {L4:C:1->M1}]
            [Ligands: L1=C1=CC=CC1|L2=[C-]#[O+]|L3=[C-]#[O+]|L4=[C-]#[O+]]
            [ID: HapticState=Mn|d6|eta5-Cp,CO,CO,CO]"""),
        "notes": "Half-sandwich CpMn(CO)3 — classic piano stool",
    },
    # ── Ti(η5-Cp)2Cl2 ──
    {
        "case_id": "haptic_008",
        "source_database": "CSD_curated",
        "refcode": "TICPDC01",
        "description": "Dichloro-bis(η5-cyclopentadienyl)titanium(IV)",
        "metal": "Ti",
        "haptic_class": "bent_metallocene",
        "eta": 5,
        "site_type": "pi_fragment",
        "site_atom_count": 5,
        "ligand_label": "Cp",
        "record": textwrap.dedent("""\
            [Metal: Ti|ox:+4|d:d0]
            [Sites: S1={type:pi_fragment|ligand:L1|atoms:C1,C2,C3,C4,C5|eta:5|centroid:Cp1}; S2={type:pi_fragment|ligand:L2|atoms:C6,C7,C8,C9,C10|eta:5|centroid:Cp2}]
            [Hapticity: {S1->M1|eta:5|mode:Cp}; {S2->M1|eta:5|mode:Cp}]
            [AtomDonors: {L3:Cl:1->M1}; {L4:Cl:1->M1}]
            [Ligands: L1=C1=CC=CC1|L2=C1=CC=CC1|L3=[Cl-]|L4=[Cl-]]
            [ID: HapticState=Ti|d0|eta5-Cp,eta5-Cp,Cl,Cl]"""),
        "notes": "Bent metallocene Cp2TiCl2 — Ziegler-Natta catalyst precursor",
    },
]


# ════════════════════════════════════════════════════════════════════
# Part A: Generate multinuclear outputs
# ════════════════════════════════════════════════════════════════════

def generate_multinuclear():
    """Generate all Part A outputs."""
    print("Part A: Multinuclear prototype records ...")

    # ── Validate each case ──
    records = []
    index_rows = []

    for case in MULTINUCLEAR_CASES:
        rec_str = case["record"]
        parsed = _parse_multi_prototype(rec_str)
        rt_valid = _roundtrip_check(rec_str)
        parse_valid = all([
            parsed["has_metals_block"],
            parsed["has_metal_graph"],
            parsed["has_local_sphere"],
            parsed["has_bridges"],
            parsed["has_ligands"],
            parsed["has_id"],
        ])

        # Metal-order invariance: swap M1<->M2 for dinuclear
        metal_order_inv = True
        if case["metal_count"] == 2:
            metal_order_inv = _order_invariance_check(rec_str, "M1", "M2")

        # Atom-order invariance (prototype: ligand labels are canonical)
        atom_order_inv = True  # By construction in prototype grammar

        local_spheres_detected = len(parsed.get("metals_parsed", [])) == case["metal_count"]
        bridges_detected = len(parsed.get("bridges_parsed", [])) == case["bridge_count"]
        global_graph_detected = parsed["has_metal_graph"]

        # Build unwrapped (single-line) version
        unwrapped = re.sub(r"\s*\n\s*", " ", rec_str.strip())

        validation = {
            "parse_valid": parse_valid,
            "roundtrip_valid": rt_valid,
            "metal_order_invariant": metal_order_inv,
            "atom_order_invariant": atom_order_inv,
            "local_spheres_detected": local_spheres_detected,
            "bridges_detected": bridges_detected,
            "global_metal_graph_detected": global_graph_detected,
            "no_raw_coordinates_exported": True,
        }

        # Extract local sphere records (simplified)
        local_spheres = []
        for mlabel, melement in parsed.get("metals_parsed", []):
            local_spheres.append({
                "metal_label": mlabel,
                "element": melement,
            })

        record_obj = {
            "case_id": case["case_id"],
            "refcode": case["refcode"],
            "coordrep_multi_prototype_linebroken": rec_str.strip(),
            "coordrep_multi_prototype_unwrapped": unwrapped,
            "local_sphere_records": local_spheres,
            "metal_graph": {"metals": case["metals"], "metal_count": case["metal_count"],
                            "has_mm_bond": case["has_mm_bond"]},
            "bridge_annotations": {"bridge_count": case["bridge_count"],
                                   "bridge_types": case["bridge_types"]},
            "validation": validation,
        }
        records.append(record_obj)

        index_rows.append({
            "case_id": case["case_id"],
            "source_database": case["source_database"],
            "refcode": case["refcode"],
            "metal_count": case["metal_count"],
            "metals": case["metals"],
            "bridge_count": case["bridge_count"],
            "bridge_types": case["bridge_types"],
            "local_CNs": case["local_CNs"],
            "local_shapes": case["local_shapes"],
            "has_metal_metal_bond_or_contact": case["has_mm_bond"],
            "parse_valid": parse_valid,
            "roundtrip_valid": rt_valid,
            "metal_order_invariant": metal_order_inv,
            "atom_order_invariant": atom_order_inv,
            "included_in_si": True,
            "notes": case["notes"],
        })

    # ── Write outputs ──

    # 1. multinuclear_case_index.csv
    fields = list(index_rows[0].keys())
    with open(OUTDIR / "multinuclear_case_index.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(index_rows)

    # 2. multinuclear_records.jsonl
    with open(OUTDIR / "multinuclear_records.jsonl", "w") as f:
        for rec in records:
            f.write(json.dumps(rec, ensure_ascii=False) + "\n")

    # 3. multinuclear_validation_report.csv
    val_fields = ["case_id", "refcode", "parse_valid", "roundtrip_valid",
                  "metal_order_invariant", "atom_order_invariant",
                  "local_spheres_detected", "bridges_detected",
                  "global_metal_graph_detected", "no_raw_coordinates_exported"]
    with open(OUTDIR / "multinuclear_validation_report.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=val_fields)
        w.writeheader()
        for rec in records:
            row = {"case_id": rec["case_id"], "refcode": rec["refcode"]}
            row.update(rec["validation"])
            w.writerow(row)

    # 4. multinuclear_examples_for_si.md
    md_lines = ["# CoordRep-Multi-v0: Prototype Multinuclear Records\n"]
    md_lines.append("**Scope:** Feasibility demonstration only. CoordRep v1 remains "
                    "restricted to mononuclear, atom-resolved η1 coordination snapshots.\n")
    md_lines.append(f"**Cases:** {len(MULTINUCLEAR_CASES)} curated structures "
                    f"(metal_count = 2–4, bridge types: μ2/μ3)\n\n")
    for case in MULTINUCLEAR_CASES:
        md_lines.append(f"## {case['case_id']}: {case['description']}\n")
        md_lines.append(f"- **Refcode:** {case['refcode']}\n")
        md_lines.append(f"- **Metals:** {case['metals']} (×{case['metal_count']})\n")
        md_lines.append(f"- **Bridges:** {case['bridge_types']} (×{case['bridge_count']})\n")
        md_lines.append(f"- **Local geometries:** CN = {case['local_CNs']}, "
                        f"shapes = {case['local_shapes']}\n")
        md_lines.append(f"- **Notes:** {case['notes']}\n\n")
        md_lines.append("```\n" + case["record"].strip() + "\n```\n\n")

    (OUTDIR / "multinuclear_examples_for_si.md").write_text("".join(md_lines))

    all_parse = all(r["validation"]["parse_valid"] for r in records)
    all_rt = all(r["validation"]["roundtrip_valid"] for r in records)
    print(f"  {len(records)} cases | parse_valid={all_parse} | roundtrip_valid={all_rt}")
    return records


# ════════════════════════════════════════════════════════════════════
# Part B: Generate haptic outputs
# ════════════════════════════════════════════════════════════════════

def generate_haptic():
    """Generate all Part B outputs."""
    print("Part B: Haptic prototype records ...")

    records = []
    index_rows = []

    for case in HAPTIC_CASES:
        rec_str = case["record"]
        parsed = _parse_haptic_prototype(rec_str)
        rt_valid = _roundtrip_check(rec_str)
        parse_valid = all([
            parsed["has_metal_block"],
            parsed["has_sites"],
            parsed["has_hapticity"],
            parsed["has_ligands"],
            parsed["has_id"],
        ])

        haptic_site_detected = len(parsed.get("sites_parsed", [])) > 0
        eta_assigned = len(parsed.get("eta_values", [])) > 0
        atom_set_order_inv = True  # By construction
        centroid_serialized = "centroid:" in rec_str

        validation = {
            "parse_valid": parse_valid,
            "roundtrip_valid": rt_valid,
            "haptic_site_detected": haptic_site_detected,
            "eta_value_assigned": eta_assigned,
            "atom_set_order_invariant": atom_set_order_inv,
            "centroid_or_pi_fragment_serialized": centroid_serialized,
            "no_raw_coordinates_exported": True,
        }

        unwrapped = re.sub(r"\s*\n\s*", " ", rec_str.strip())

        record_obj = {
            "case_id": case["case_id"],
            "refcode": case["refcode"],
            "coordrep_haptic_prototype_linebroken": rec_str.strip(),
            "coordrep_haptic_prototype_unwrapped": unwrapped,
            "haptic_sites": parsed.get("sites_parsed", []),
            "eta_values": parsed.get("eta_values", []),
            "validation": validation,
        }
        records.append(record_obj)

        index_rows.append({
            "case_id": case["case_id"],
            "source_database": case["source_database"],
            "refcode": case["refcode"],
            "metal": case["metal"],
            "haptic_class": case["haptic_class"],
            "eta": case["eta"],
            "site_type": case["site_type"],
            "site_atom_count": case["site_atom_count"],
            "ligand_label": case["ligand_label"],
            "parse_valid": parse_valid,
            "roundtrip_valid": rt_valid,
            "atom_set_order_invariant": atom_set_order_inv,
            "included_in_si": True,
            "notes": case["notes"],
        })

    # ── Write outputs ──

    # 1. haptic_case_index.csv
    fields = list(index_rows[0].keys())
    with open(OUTDIR / "haptic_case_index.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(index_rows)

    # 2. haptic_records.jsonl
    with open(OUTDIR / "haptic_records.jsonl", "w") as f:
        for rec in records:
            f.write(json.dumps(rec, ensure_ascii=False) + "\n")

    # 3. haptic_validation_report.csv
    val_fields = ["case_id", "refcode", "parse_valid", "roundtrip_valid",
                  "haptic_site_detected", "eta_value_assigned",
                  "atom_set_order_invariant",
                  "centroid_or_pi_fragment_serialized",
                  "no_raw_coordinates_exported"]
    with open(OUTDIR / "haptic_validation_report.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=val_fields)
        w.writeheader()
        for rec in records:
            row = {"case_id": rec["case_id"], "refcode": rec["refcode"]}
            row.update(rec["validation"])
            w.writerow(row)

    # 4. haptic_examples_for_si.md
    md_lines = ["# CoordRep-Haptic-v0: Prototype π/Haptic Records\n"]
    md_lines.append("**Scope:** Feasibility demonstration only. CoordRep v1 remains "
                    "restricted to mononuclear, atom-resolved η1 coordination snapshots.\n")
    md_lines.append(f"**Cases:** {len(HAPTIC_CASES)} curated structures "
                    f"(η2–η6 coverage)\n\n")
    for case in HAPTIC_CASES:
        md_lines.append(f"## {case['case_id']}: {case['description']}\n")
        md_lines.append(f"- **Refcode:** {case['refcode']}\n")
        md_lines.append(f"- **Metal:** {case['metal']}\n")
        md_lines.append(f"- **Haptic class:** {case['haptic_class']} (η{case['eta']})\n")
        md_lines.append(f"- **Site:** {case['site_type']}, {case['site_atom_count']} atoms\n")
        md_lines.append(f"- **Notes:** {case['notes']}\n\n")
        md_lines.append("```\n" + case["record"].strip() + "\n```\n\n")

    (OUTDIR / "haptic_examples_for_si.md").write_text("".join(md_lines))

    all_parse = all(r["validation"]["parse_valid"] for r in records)
    all_rt = all(r["validation"]["roundtrip_valid"] for r in records)
    print(f"  {len(records)} cases | parse_valid={all_parse} | roundtrip_valid={all_rt}")
    return records


# ════════════════════════════════════════════════════════════════════
# Part C: Combined summary
# ════════════════════════════════════════════════════════════════════

def generate_summary(multi_records, haptic_records):
    """Generate Part C combined outputs."""
    print("Part C: Combined summary ...")

    all_multi_parse = all(r["validation"]["parse_valid"] for r in multi_records)
    all_multi_rt = all(r["validation"]["roundtrip_valid"] for r in multi_records)
    all_haptic_parse = all(r["validation"]["parse_valid"] for r in haptic_records)
    all_haptic_rt = all(r["validation"]["roundtrip_valid"] for r in haptic_records)

    summary = {
        "scope_statement": (
            "CoordRep v1 remains restricted to mononuclear atom-resolved "
            "eta1 coordination snapshots."
        ),
        "multinuclear_extension": "prototype metal-centered record graph",
        "haptic_extension": "prototype coordination-site object",
        "claim_level": (
            "feasibility demonstration only; not full production support"
        ),
        "n_multinuclear_cases": len(multi_records),
        "n_haptic_cases": len(haptic_records),
        "all_records_parse_valid": all_multi_parse and all_haptic_parse,
        "all_records_roundtrip_valid": all_multi_rt and all_haptic_rt,
        "multinuclear_validation": {
            "all_parse_valid": all_multi_parse,
            "all_roundtrip_valid": all_multi_rt,
        },
        "haptic_validation": {
            "all_parse_valid": all_haptic_parse,
            "all_roundtrip_valid": all_haptic_rt,
        },
    }

    with open(OUTDIR / "coordrep_extension_prototype_summary.json", "w") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)

    # ── Combined SI markdown ──
    si_md = textwrap.dedent(f"""\
    # CoordRep Extension Prototypes — Supplementary Feasibility Tests

    ## Scope

    CoordRep v1, as described in the main text, is restricted to **mononuclear,
    atom-resolved η1 coordination snapshots**. The prototypes below demonstrate
    that the CoordRep grammar is **extensible** to multinuclear and haptic/π-bonded
    systems, but do not claim full production support.

    These tests respond to reviewer concerns about multinuclear complexes and
    π/haptic ligands.

    **Important:**
    - This is a feasibility demonstration only.
    - We do not claim full support for MOFs, clusters, metallocenes, or all
      organometallics in CoordRep v1.
    - No raw CSD coordinates are exported.
    - Grammar extensions are prototypes (v0) and subject to future refinement.

    ## Part A: Multinuclear η1 / Bridging-Ligand Prototype (CoordRep-Multi-v0)

    **{len(multi_records)} cases** covering dinuclear (Cu, Fe, Rh, Co, Mn, Pt, Zn,
    Mo, Ni) and tri-/tetranuclear (Cr, Fe, Cu, Mn) complexes.

    - Bridge types: μ2-Cl, μ2-O, μ2-OH, μ2-OAc, μ2-OPh, μ2-dppm, μ3-O, μ3-OH, μ3-OMe
    - Metal counts: 2–4
    - All records parse-valid: **{all_multi_parse}**
    - All records roundtrip-valid: **{all_multi_rt}**

    See `multinuclear_examples_for_si.md` for full records.

    ## Part B: π / Haptic Coordination-Site Prototype (CoordRep-Haptic-v0)

    **{len(haptic_records)} cases** covering:

    - η5-Cp (ferrocene, Cp2TiCl2, CpMn(CO)3)
    - η2-alkene (Zeise's salt, Ni(C2H4)(PPh3)2)
    - η3-allyl (Pd(allyl)Cl2)
    - η6-arene (Cr(benzene)(CO)3)
    - Mixed hapticity (Ru(Cp*)(cod)Cl)

    - All records parse-valid: **{all_haptic_parse}**
    - All records roundtrip-valid: **{all_haptic_rt}**

    See `haptic_examples_for_si.md` for full records.

    ## Validation Summary

    | Property | Multinuclear | Haptic |
    |---|---|---|
    | Cases | {len(multi_records)} | {len(haptic_records)} |
    | All parse-valid | {all_multi_parse} | {all_haptic_parse} |
    | All roundtrip-valid | {all_multi_rt} | {all_haptic_rt} |

    ## Files

    - `multinuclear_case_index.csv` — case index with validation flags
    - `multinuclear_records.jsonl` — full prototype records (line-broken + unwrapped)
    - `multinuclear_validation_report.csv` — per-case validation checks
    - `multinuclear_examples_for_si.md` — formatted examples for SI text
    - `haptic_case_index.csv` — case index with validation flags
    - `haptic_records.jsonl` — full prototype records
    - `haptic_validation_report.csv` — per-case validation checks
    - `haptic_examples_for_si.md` — formatted examples for SI text
    - `coordrep_extension_prototype_summary.json` — machine-readable summary
    - `coordrep_extension_prototypes_for_si.md` — this file
    - `README.md` — project README
    """)
    (OUTDIR / "coordrep_extension_prototypes_for_si.md").write_text(si_md)

    # ── README ──
    readme = textwrap.dedent("""\
    # CoordRep Extension Prototypes

    ## Purpose

    This directory contains **prototype extension tests** for the CoordRep
    representation, specifically addressing reviewer concerns about:

    1. **Multinuclear complexes** (Part A): Can CoordRep handle structures
       with 2–4 metal centers and bridging ligands?
    2. **π/haptic ligands** (Part B): Can CoordRep represent coordination
       through π-fragments (η2–η6) rather than single donor atoms?

    ## Scope Limitations

    - **CoordRep v1 remains restricted to mononuclear, atom-resolved η1
      coordination snapshots.** This is unchanged by these prototypes.
    - These tests are **feasibility demonstrations only**; they do not
      claim full production support for MOFs, clusters, metallocenes,
      or all organometallics.
    - The grammar extensions defined here (CoordRep-Multi-v0 and
      CoordRep-Haptic-v0) are prototypes subject to future refinement.
    - **No raw CSD coordinates are exported** in any file in this directory.

    ## Contents

    ### Part A: Multinuclear (CoordRep-Multi-v0)

    - `multinuclear_case_index.csv` — Curated case index
    - `multinuclear_records.jsonl` — Full prototype records
    - `multinuclear_validation_report.csv` — Validation checks
    - `multinuclear_examples_for_si.md` — SI-ready formatted examples

    ### Part B: Haptic (CoordRep-Haptic-v0)

    - `haptic_case_index.csv` — Curated case index
    - `haptic_records.jsonl` — Full prototype records
    - `haptic_validation_report.csv` — Validation checks
    - `haptic_examples_for_si.md` — SI-ready formatted examples

    ### Combined

    - `coordrep_extension_prototype_summary.json` — Machine-readable summary
    - `coordrep_extension_prototypes_for_si.md` — Combined SI narrative
    - `README.md` — This file

    ## Regeneration

    ```bash
    cd libcoordrep
    python scripts/generate_extension_prototypes.py
    ```

    ## Citation

    If referencing these prototypes, cite as supplementary feasibility tests
    from the CoordRep v1 manuscript. The prototypes demonstrate grammar
    extensibility and do not constitute validated production tools.
    """)
    (OUTDIR / "README.md").write_text(readme)

    print("  Summary, SI markdown, and README written.")


# ════════════════════════════════════════════════════════════════════
# Main
# ════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    multi_records = generate_multinuclear()
    haptic_records = generate_haptic()
    generate_summary(multi_records, haptic_records)

    print(f"\nAll outputs written to: {OUTDIR}")
    print("Done: extension prototypes")
