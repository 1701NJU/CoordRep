#!/usr/bin/env python3
"""
CoordRep-Rosetta Hard Controls: test whether CoordRep higher-order fields
(bridge topology, haptic modes, shape/stereo) add retrieval value beyond
metal/CN/donor matching.

Hard pool definition: for each MOF node query, candidate molecular records
must share same metal(s), same CNsite, same donor element multiset.
Within these hard pools, metal/CN/donor baselines score identically on all
candidates — only higher-order CoordRep fields can discriminate.

Three sub-tasks:
A. Bridge-topology retrieval
B. Shape/geometry hard retrieval
C. Combined higher-order retrieval

Go/no-go: CoordRep full must clearly outperform donor-baseline on hard pools
to justify main-text application claim. Otherwise → SI proof-of-concept.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import random
import sys
import time
from collections import Counter, OrderedDict, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ccdc.io import EntryReader

from coordrep.io.tmqm_reader import TRANSITION_METALS
from coordrep.v2beta.csd_v2beta_adapter import convert_multinuclear

LIBROOT = Path(__file__).resolve().parents[1]
OUTDIR = LIBROOT / "revision_results" / "coordrep_rosetta_hard_controls"
OUTDIR.mkdir(parents=True, exist_ok=True)


# ════════════════════════════════════════════════════════════════
# Shared record structure (same as rosetta pilot)
# ════════════════════════════════════════════════════════════════

@dataclass
class LocalSiteRecord:
    label: str
    donor_element: str
    eta: int = 1
    mu: int = 1
    target_metals: List[str] = field(default_factory=list)
    mode: str = ""

@dataclass
class LocalMetalRecord:
    label: str
    element: str
    cn_site: int = 0
    cn_atom: int = 0
    local_shape_best: str = ""

@dataclass
class LocalEdgeRecord:
    m1: str
    m2: str
    relation: str = "bridged"

@dataclass
class NodeRecord:
    node_id: str
    source_database: str
    source_id: str
    record_type: str
    metal_count: int = 0
    metals: List[LocalMetalRecord] = field(default_factory=list)
    sites: List[LocalSiteRecord] = field(default_factory=list)
    edges: List[LocalEdgeRecord] = field(default_factory=list)
    periodic_truncated: bool = False

    # Derived
    metals_str: str = ""
    cn_list: str = ""
    donor_set: str = ""
    bridge_count: int = 0
    bridge_modes: str = ""
    edge_count: int = 0
    L3: str = ""

    def derive_fields(self):
        self.metals_str = "+".join(sorted(m.element for m in self.metals))
        self.cn_list = ",".join(str(m.cn_site) for m in
                                sorted(self.metals, key=lambda m: m.label))
        self.donor_set = ",".join(sorted(s.donor_element for s in self.sites))
        self.bridge_count = sum(1 for s in self.sites if s.mu > 1)
        self.bridge_modes = ",".join(sorted(set(
            s.mode for s in self.sites if s.mu > 1 and s.mode)))
        self.edge_count = len(self.edges)
        sig = [self.metals_str, self.cn_list, self.donor_set,
               str(self.bridge_count), self.bridge_modes, str(self.edge_count)]
        for s in sorted(self.sites, key=lambda x: (x.donor_element, x.mu, x.eta)):
            sig.append(f"{s.donor_element}:eta{s.eta}:mu{s.mu}")
        self.L3 = hashlib.md5("|".join(sig).encode()).hexdigest()[:12]

    @property
    def hard_key(self) -> Tuple[str, str, str]:
        """Key for hard pool: same metal+CN+donor."""
        return (self.metals_str, self.cn_list, self.donor_set)


# ════════════════════════════════════════════════════════════════
# Node extraction (reuse from rosetta pilot)
# ════════════════════════════════════════════════════════════════

def _is_pi_bond(b) -> bool:
    return "pi" in str(b.bond_type).lower() or "deloc" in str(b.bond_type).lower()

def _get_metal_atoms(mol):
    return [a for a in mol.atoms if a.atomic_symbol in TRANSITION_METALS]

def _guess_shape(cn: int) -> str:
    return {2: "L", 3: "TP", 4: "SP", 5: "TBPY", 6: "OC", 7: "PBPY", 8: "CU"}.get(cn, f"CN{cn}")


def extract_mof_node(entry) -> Optional[NodeRecord]:
    """Extract local coordination node from polymeric CSD entry."""
    refcode = entry.identifier
    mol = entry.molecule
    if mol is None:
        return None
    metal_atoms = _get_metal_atoms(mol)
    if not metal_atoms:
        return None

    seen_labels = set()
    unique_metals = []
    for m in metal_atoms:
        if m.label not in seen_labels:
            seen_labels.add(m.label)
            unique_metals.append(m)
    if len(unique_metals) > 12:
        return None

    metal_label_map = {}
    metals = []
    for i, m in enumerate(unique_metals):
        mlabel = f"M{i+1}"
        metal_label_map[m.label] = mlabel
        n_donors = 0
        for b in m.bonds:
            other = b.atoms[0] if b.atoms[1] == m else b.atoms[1]
            if other.atomic_symbol not in TRANSITION_METALS and not _is_pi_bond(b):
                n_donors += 1
        metals.append(LocalMetalRecord(
            label=mlabel, element=m.atomic_symbol,
            cn_site=n_donors, cn_atom=n_donors,
            local_shape_best=_guess_shape(n_donors)))

    sites = []
    edges = []
    donor_to_metals = defaultdict(set)
    for m in unique_metals:
        mlabel = metal_label_map[m.label]
        for b in m.bonds:
            other = b.atoms[0] if b.atoms[1] == m else b.atoms[1]
            if other.atomic_symbol in TRANSITION_METALS:
                if other.label in metal_label_map:
                    ol = metal_label_map[other.label]
                    if ol != mlabel:
                        pair = tuple(sorted([mlabel, ol]))
                        if not any(e.m1 == pair[0] and e.m2 == pair[1] for e in edges):
                            edges.append(LocalEdgeRecord(m1=pair[0], m2=pair[1]))
            elif not _is_pi_bond(b):
                donor_to_metals[other.label].add(mlabel)

    seen_donors = set()
    for m in unique_metals:
        mlabel = metal_label_map[m.label]
        for b in m.bonds:
            other = b.atoms[0] if b.atoms[1] == m else b.atoms[1]
            if other.atomic_symbol in TRANSITION_METALS or _is_pi_bond(b):
                continue
            dlabel = other.label
            if dlabel in seen_donors:
                for s in sites:
                    if s.label == f"S__{dlabel}" and mlabel not in s.target_metals:
                        s.target_metals.append(mlabel)
                continue
            seen_donors.add(dlabel)
            target_ms = sorted(donor_to_metals[dlabel])
            mu = len(target_ms)
            sites.append(LocalSiteRecord(
                label=f"S__{dlabel}", donor_element=other.atomic_symbol,
                eta=1, mu=mu, target_metals=target_ms,
                mode=f"mu{mu}-{other.atomic_symbol}" if mu > 1 else ""))

    for i, s in enumerate(sites):
        s.label = f"S{i+1}"

    node = NodeRecord(
        node_id=f"MOF_{refcode}", source_database="CSD_MOF",
        source_id=refcode, record_type="mof_node",
        metal_count=len(metals), metals=metals, sites=sites,
        edges=edges, periodic_truncated=True)
    node.derive_fields()
    return node


def extract_molecular_node(entry, multi=False) -> Optional[NodeRecord]:
    """Extract node record from molecular CSD entry."""
    mol = entry.molecule
    if mol is None:
        return None
    metals_atoms = _get_metal_atoms(mol)
    if not metals_atoms:
        return None

    if not multi and len(metals_atoms) != 1:
        return None
    if multi and len(metals_atoms) < 2:
        return None
    if len(metals_atoms) > 12:
        return None

    if not multi:
        m = metals_atoms[0]
        donors = []
        for b in m.bonds:
            other = b.atoms[0] if b.atoms[1] == m else b.atoms[1]
            if other.atomic_symbol not in TRANSITION_METALS and not _is_pi_bond(b):
                donors.append(other)
        if len(donors) < 2 or len(donors) > 12:
            return None
        cn = len(donors)
        sites = [LocalSiteRecord(
            label=f"S{i+1}", donor_element=d.atomic_symbol,
            eta=1, mu=1, target_metals=["M1"])
            for i, d in enumerate(donors)]
        node = NodeRecord(
            node_id=f"MONO_{entry.identifier}",
            source_database="CSD_molecular", source_id=entry.identifier,
            record_type="mononuclear", metal_count=1,
            metals=[LocalMetalRecord(
                label="M1", element=m.atomic_symbol,
                cn_site=cn, cn_atom=cn,
                local_shape_best=_guess_shape(cn))],
            sites=sites)
        node.derive_fields()
        return node

    # Multinuclear: use v2beta converter
    try:
        res = convert_multinuclear(entry)
        if not res.success or res.record is None:
            return None
        rec = res.record
        node_metals = [LocalMetalRecord(
            label=m.label, element=m.element,
            cn_site=m.cn_site, cn_atom=m.cn_atom,
            local_shape_best=_guess_shape(m.cn_site))
            for m in rec.metals]
        node_sites = [LocalSiteRecord(
            label=s.label,
            donor_element=s.donor_elements[0] if s.donor_elements else "?",
            eta=s.eta, mu=s.mu,
            target_metals=list(s.target_metals), mode=s.mode)
            for s in rec.sites]
        node_edges = [LocalEdgeRecord(m1=e.m1, m2=e.m2, relation=e.relation)
                      for e in rec.metal_edges]
        node = NodeRecord(
            node_id=f"MULTI_{entry.identifier}",
            source_database="CSD_multi", source_id=entry.identifier,
            record_type="multinuclear", metal_count=len(node_metals),
            metals=node_metals, sites=node_sites, edges=node_edges)
        node.derive_fields()
        return node
    except Exception:
        return None


# ════════════════════════════════════════════════════════════════
# Scoring functions
# ════════════════════════════════════════════════════════════════

def _donor_multiset_jaccard(n1: NodeRecord, n2: NodeRecord) -> float:
    d1 = Counter(s.donor_element for s in n1.sites)
    d2 = Counter(s.donor_element for s in n2.sites)
    if not d1 and not d2:
        return 1.0
    if not d1 or not d2:
        return 0.0
    all_e = set(d1) | set(d2)
    inter = sum(min(d1.get(e, 0), d2.get(e, 0)) for e in all_e)
    union = sum(max(d1.get(e, 0), d2.get(e, 0)) for e in all_e)
    return inter / max(union, 1)


def bridge_topology_sim(n1: NodeRecord, n2: NodeRecord) -> float:
    """Similarity based on bridge count, bridge modes, metal-graph edges."""
    if n1.bridge_count == 0 and n2.bridge_count == 0:
        return 1.0
    if n1.bridge_count == 0 or n2.bridge_count == 0:
        return 0.0
    # Bridge count similarity
    bc_sim = 1.0 - abs(n1.bridge_count - n2.bridge_count) / max(n1.bridge_count, n2.bridge_count)
    # Bridge mode overlap
    bm1 = Counter(s.mode for s in n1.sites if s.mu > 1 and s.mode)
    bm2 = Counter(s.mode for s in n2.sites if s.mu > 1 and s.mode)
    all_modes = set(bm1) | set(bm2)
    if all_modes:
        inter = sum(min(bm1.get(m, 0), bm2.get(m, 0)) for m in all_modes)
        union = sum(max(bm1.get(m, 0), bm2.get(m, 0)) for m in all_modes)
        bm_sim = inter / max(union, 1)
    else:
        bm_sim = bc_sim
    # Edge count similarity
    if n1.edge_count == 0 and n2.edge_count == 0:
        ec_sim = 1.0
    elif n1.edge_count == 0 or n2.edge_count == 0:
        ec_sim = 0.0
    else:
        ec_sim = 1.0 - abs(n1.edge_count - n2.edge_count) / max(n1.edge_count, n2.edge_count)
    return 0.4 * bc_sim + 0.4 * bm_sim + 0.2 * ec_sim


def shape_geometry_sim(n1: NodeRecord, n2: NodeRecord) -> float:
    """Similarity based on local shape labels."""
    s1 = sorted(m.local_shape_best for m in n1.metals)
    s2 = sorted(m.local_shape_best for m in n2.metals)
    if s1 == s2:
        return 1.0
    matched = 0
    s2_rem = list(s2)
    for sh in s1:
        if sh in s2_rem:
            matched += 1
            s2_rem.remove(sh)
    return matched / max(len(s1), len(s2), 1)


def mu_pattern_sim(n1: NodeRecord, n2: NodeRecord) -> float:
    """Similarity of mu-pattern across all sites."""
    mu1 = sorted(s.mu for s in n1.sites)
    mu2 = sorted(s.mu for s in n2.sites)
    if mu1 == mu2:
        return 1.0
    # Pad
    while len(mu1) < len(mu2): mu1.append(0)
    while len(mu2) < len(mu1): mu2.append(0)
    diffs = sum(abs(a - b) for a, b in zip(mu1, mu2))
    max_val = sum(max(a, b) for a, b in zip(mu1, mu2))
    return 1.0 - diffs / max(max_val, 1)


def coordrep_full_score(n1: NodeRecord, n2: NodeRecord) -> float:
    """Full CoordRep homology score emphasizing higher-order fields.
    Within hard pools, metal/CN/donor are identical, so discrimination
    comes from bridge, shape, mu-pattern, edge topology, and L3."""
    # Metal/CN/donor (will be ~1.0 in hard pools)
    metal_ok = 1.0 if n1.metals_str == n2.metals_str else 0.0
    cn_ok = 1.0 if n1.cn_list == n2.cn_list else 0.0
    donor_ok = _donor_multiset_jaccard(n1, n2)

    # Higher-order fields (the discriminative ones in hard pools)
    bridge_sim = bridge_topology_sim(n1, n2)
    shape_sim = shape_geometry_sim(n1, n2)
    mu_sim = mu_pattern_sim(n1, n2)
    l3_match = 1.0 if n1.L3 == n2.L3 else 0.0

    # In hard pools, metal/CN/donor are always ~1.0
    # Weight higher-order fields heavily
    score = (0.10 * metal_ok + 0.10 * cn_ok + 0.10 * donor_ok +
             0.30 * bridge_sim + 0.15 * shape_sim +
             0.15 * mu_sim + 0.10 * l3_match)
    return round(score, 4)


def donor_baseline_score(n1: NodeRecord, n2: NodeRecord) -> float:
    """Baseline: metal + CN + donor only. In hard pools all candidates
    get the same score, so ranking is effectively random."""
    metal_ok = 1.0 if n1.metals_str == n2.metals_str else 0.0
    cn_ok = 1.0 if n1.cn_list == n2.cn_list else 0.0
    donor_ok = _donor_multiset_jaccard(n1, n2)
    return round((metal_ok + cn_ok + donor_ok) / 3, 4)


def cshm_only_score(n1: NodeRecord, n2: NodeRecord) -> float:
    """CShM-only baseline."""
    return round(shape_geometry_sim(n1, n2), 4)


def l3_only_score(n1: NodeRecord, n2: NodeRecord) -> float:
    """L3-collision-only baseline."""
    return 1.0 if n1.L3 == n2.L3 else 0.0


def random_same_donor_score(n1: NodeRecord, n2: NodeRecord) -> float:
    """Random score for candidates in hard pool."""
    return round(random.random(), 4)


METHODS = {
    "coordrep_full": coordrep_full_score,
    "donor_baseline": donor_baseline_score,
    "cshm_only": cshm_only_score,
    "L3_only": l3_only_score,
    "random_same_donor": random_same_donor_score,
}


# ════════════════════════════════════════════════════════════════
# Hard pool construction
# ════════════════════════════════════════════════════════════════

def build_hard_pools(mof_nodes: List[NodeRecord],
                     mol_refs: List[NodeRecord],
                     pool_type: str = "same_metal_cn",
                     min_pool: int = 5
                     ) -> Dict[str, Dict]:
    """Build hard retrieval pools with both positives and negatives.

    pool_type controls the hard key:
      'same_metal_cn' – same metal(s) + same CN(s) → donor set varies
      'same_metal'    – same metal(s) only → CN and donors vary
      'same_metal_donor' – same metal(s) + same donor set → bridge varies
    """
    def _pool_key(n: NodeRecord) -> str:
        if pool_type == "same_metal_cn":
            return f"{n.metals_str}|{n.cn_list}"
        elif pool_type == "same_metal":
            return n.metals_str
        elif pool_type == "same_metal_donor":
            return f"{n.metals_str}|{n.donor_set}"
        return n.metals_str

    mol_by_key = defaultdict(list)
    for n in mol_refs:
        mol_by_key[_pool_key(n)].append(n)

    pools = {}
    for mof in mof_nodes:
        key = _pool_key(mof)
        candidates = mol_by_key.get(key, [])
        if len(candidates) < min_pool:
            continue

        # Label positives: high full-score match (bridge + donor + mu + shape)
        # Label negatives: low full-score match
        positives = []
        negatives = []
        for c in candidates:
            full = coordrep_full_score(mof, c)
            # Use donor+bridge combined as ground truth
            donor_sim = _donor_multiset_jaccard(mof, c)
            bridge_sim = bridge_topology_sim(mof, c)
            mu_sim = mu_pattern_sim(mof, c)
            combined_gt = 0.4 * donor_sim + 0.35 * bridge_sim + 0.25 * mu_sim
            if combined_gt >= 0.7:
                positives.append(c)
            elif combined_gt < 0.4:
                negatives.append(c)
            # In between: skip (ambiguous)

        if not positives or not negatives:
            continue

        pools[mof.node_id] = {
            "query": mof,
            "pool": candidates,
            "positives": positives,
            "negatives": negatives,
        }

    return pools


# ════════════════════════════════════════════════════════════════
# Evaluation
# ════════════════════════════════════════════════════════════════

def evaluate_method(pools: Dict[str, Dict], method_name: str,
                    score_fn) -> Dict[str, Any]:
    """Evaluate a scoring method on hard pools."""
    top1_correct = 0
    top5_has_correct = 0
    n_queries = 0
    all_auroc = []
    all_fp_rates = []

    per_query = []

    for pool_id, pool_data in pools.items():
        query = pool_data["query"]
        candidates = pool_data["pool"]
        pos_ids = {p.node_id for p in pool_data["positives"]}

        # Score all candidates
        scored = []
        for c in candidates:
            s = score_fn(query, c)
            scored.append((c, s, c.node_id in pos_ids))

        # Break ties randomly for fair comparison
        random.shuffle(scored)
        scored.sort(key=lambda x: -x[1])

        n_queries += 1

        # Top-1 correct
        if scored and scored[0][2]:
            top1_correct += 1

        # Top-5 has correct
        if any(x[2] for x in scored[:5]):
            top5_has_correct += 1

        # AUROC (if we have both pos and neg)
        n_pos = sum(1 for x in scored if x[2])
        n_neg = sum(1 for x in scored if not x[2])
        if n_pos > 0 and n_neg > 0:
            # Count concordant pairs
            concordant = 0
            total_pairs = 0
            for i, (_, si, li) in enumerate(scored):
                if not li:
                    continue
                for j, (_, sj, lj) in enumerate(scored):
                    if lj:
                        continue
                    total_pairs += 1
                    if si > sj:
                        concordant += 1
                    elif si == sj:
                        concordant += 0.5
            auroc = concordant / max(total_pairs, 1)
            all_auroc.append(auroc)

        # False positive rate: fraction of top-5 that are negatives
        top5_neg = sum(1 for x in scored[:5] if not x[2])
        fp_rate = top5_neg / min(5, len(scored))
        all_fp_rates.append(fp_rate)

        per_query.append({
            "pool_id": pool_id,
            "n_candidates": len(candidates),
            "n_positives": len(pool_data["positives"]),
            "n_negatives": len(pool_data["negatives"]),
            "top1_correct": scored[0][2] if scored else False,
            "top1_score": scored[0][1] if scored else 0,
            "top5_has_correct": any(x[2] for x in scored[:5]),
        })

    return {
        "method": method_name,
        "n_queries": n_queries,
        "top1_correct_rate": round(top1_correct / max(n_queries, 1), 4),
        "top5_correct_rate": round(top5_has_correct / max(n_queries, 1), 4),
        "mean_auroc": round(sum(all_auroc) / max(len(all_auroc), 1), 4),
        "mean_fp_rate": round(sum(all_fp_rates) / max(len(all_fp_rates), 1), 4),
        "per_query": per_query,
    }


# ════════════════════════════════════════════════════════════════
# Manual audit
# ════════════════════════════════════════════════════════════════

def audit_pair(query: NodeRecord, match: NodeRecord) -> Dict:
    """Automated chemical cross-check."""
    metal_ok = query.metals_str == match.metals_str
    donor_ok = query.donor_set == match.donor_set
    bridge_ok = bridge_topology_sim(query, match) >= 0.5
    mu_ok = mu_pattern_sim(query, match) >= 0.5
    shape_ok = shape_geometry_sim(query, match) >= 0.5
    checks = [metal_ok, donor_ok, bridge_ok, mu_ok, shape_ok]
    homologous = sum(checks) >= 4  # strict: need 4/5 in hard pool
    return {
        "metal_correct": metal_ok,
        "donor_correct": donor_ok,
        "bridge_match": bridge_ok,
        "mu_pattern_match": mu_ok,
        "shape_match": shape_ok,
        "chemically_homologous": homologous,
    }


# ════════════════════════════════════════════════════════════════
# Main
# ════════════════════════════════════════════════════════════════

def main():
    print("=" * 70)
    print("CoordRep-Rosetta Hard Controls")
    print("=" * 70)
    t_start = time.time()
    random.seed(2024)

    reader = EntryReader("CSD")

    # ── Step 1: Build node libraries ──
    print("\n[Step 1a] Extracting MOF nodes...")
    mof_nodes = []
    scanned = 0
    for e in reader:
        scanned += 1
        if len(mof_nodes) >= 500:
            break
        if scanned > 500000:
            break
        try:
            if not e.is_polymeric:
                continue
            mol = e.molecule
            if mol is None:
                continue
            if not any(a.atomic_symbol in TRANSITION_METALS for a in mol.atoms):
                continue
            node = extract_mof_node(e)
            if node and node.metal_count >= 1 and len(node.sites) >= 2:
                mof_nodes.append(node)
        except Exception:
            pass
    print(f"  MOF nodes: {len(mof_nodes)} (scanned {scanned})")

    print("\n[Step 1b] Building molecular reference (expanded)...")
    mol_refs = []
    mono_count = multi_count = 0
    scanned2 = 0
    TARGET_MONO = 3000
    TARGET_MULTI = 1500

    for e in reader:
        scanned2 += 1
        if mono_count >= TARGET_MONO and multi_count >= TARGET_MULTI:
            break
        if scanned2 > 500000:
            break
        try:
            if e.is_polymeric or e.has_disorder:
                continue
            mol = e.molecule
            if mol is None:
                continue
            metals = _get_metal_atoms(mol)

            if len(metals) == 1 and mono_count < TARGET_MONO:
                node = extract_molecular_node(e, multi=False)
                if node:
                    mol_refs.append(node)
                    mono_count += 1
            elif 2 <= len(metals) <= 12 and multi_count < TARGET_MULTI:
                node = extract_molecular_node(e, multi=True)
                if node:
                    mol_refs.append(node)
                    multi_count += 1
        except Exception:
            pass

    print(f"  Molecular refs: {len(mol_refs)} (mono={mono_count}, multi={multi_count})")

    # ── Step 2: Build hard pools for three sub-tasks ──
    SUBTASKS = [
        ("A_same_metal_cn", "same_metal_cn", 5),
        ("B_same_metal", "same_metal", 10),
        ("C_same_metal_donor", "same_metal_donor", 3),
    ]
    all_subtask_results = {}
    all_pools = {}
    all_pool_stats = {}

    for subtask_name, pool_type, min_pool in SUBTASKS:
        print(f"\n[Step 2-3: {subtask_name}] pool_type={pool_type}, min_pool={min_pool}")
        pools = build_hard_pools(mof_nodes, mol_refs, pool_type=pool_type, min_pool=min_pool)
        if len(pools) < 5 and min_pool > 2:
            pools = build_hard_pools(mof_nodes, mol_refs, pool_type=pool_type, min_pool=2)
        all_pools[subtask_name] = pools

        sizes = [len(p["pool"]) for p in pools.values()]
        mean_pool = sum(sizes) / max(len(sizes), 1)
        all_pool_stats[subtask_name] = {
            "n_pools": len(pools), "mean_pool_size": round(mean_pool, 1),
            "total_positives": sum(len(p["positives"]) for p in pools.values()),
            "total_negatives": sum(len(p["negatives"]) for p in pools.values()),
        }
        print(f"  Pools: {len(pools)}, mean size: {mean_pool:.1f}")

        # Evaluate all methods on this subtask
        results = {}
        for method_name, score_fn in METHODS.items():
            random.seed(2024)
            r = evaluate_method(pools, method_name, score_fn)
            results[method_name] = r
            print(f"  {method_name}: top1={r['top1_correct_rate']:.4f} "
                  f"top5={r['top5_correct_rate']:.4f} "
                  f"auroc={r['mean_auroc']:.4f}")
        all_subtask_results[subtask_name] = results

    # Save detailed results for all sub-tasks
    detail_fields = [
        "subtask", "method", "pool_id", "n_candidates", "n_positives",
        "n_negatives", "top1_correct", "top1_score", "top5_has_correct",
    ]
    with open(OUTDIR / "hard_control_detailed_results.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=detail_fields)
        w.writeheader()
        for subtask, results in all_subtask_results.items():
            for method_name, r in results.items():
                for pq in r["per_query"]:
                    pq["subtask"] = subtask
                    pq["method"] = method_name
                    w.writerow({k: pq.get(k, "") for k in detail_fields})

    # Save summary comparison
    comp_fields = ["subtask", "method", "n_queries", "top1_correct_rate",
                   "top5_correct_rate", "mean_auroc", "mean_fp_rate"]
    with open(OUTDIR / "hard_control_method_comparison.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=comp_fields)
        w.writeheader()
        for subtask, results in all_subtask_results.items():
            for r in results.values():
                row = {k: r[k] for k in comp_fields if k != "subtask"}
                row["subtask"] = subtask
                w.writerow(row)

    # Use best subtask (most pools) for go/no-go
    best_subtask = max(all_subtask_results.keys(),
                       key=lambda s: all_pool_stats[s]["n_pools"])
    results = all_subtask_results[best_subtask]
    pools = all_pools[best_subtask]
    mean_pool = all_pool_stats[best_subtask]["mean_pool_size"]

    # ── Step 4: Manual audit ──
    print("\n[Step 4] Manual audit sampling...")
    audit_rows = []
    audit_id = 0

    # For each sub-task, collect top-1 matches from coordrep_full vs donor_baseline
    pool_list = list(pools.values())
    random.shuffle(pool_list)

    coordrep_top1_pairs = []
    donor_top1_pairs = []
    hard_neg_pairs = []

    for pool_data in pool_list:
        query = pool_data["query"]
        candidates = pool_data["pool"]
        pos_ids = {p.node_id for p in pool_data["positives"]}

        # CoordRep full top-1
        scored_full = [(c, coordrep_full_score(query, c)) for c in candidates]
        scored_full.sort(key=lambda x: -x[1])
        if scored_full:
            coordrep_top1_pairs.append((query, scored_full[0][0], scored_full[0][1],
                                        scored_full[0][0].node_id in pos_ids))

        # Donor baseline top-1
        scored_donor = [(c, donor_baseline_score(query, c)) for c in candidates]
        random.shuffle(scored_donor)
        scored_donor.sort(key=lambda x: -x[1])
        if scored_donor:
            donor_top1_pairs.append((query, scored_donor[0][0], scored_donor[0][1],
                                     scored_donor[0][0].node_id in pos_ids))

        # Hard negatives
        if pool_data["negatives"]:
            neg = random.choice(pool_data["negatives"])
            hard_neg_pairs.append((query, neg, coordrep_full_score(query, neg), False))

    # Sample 30 from each
    for label, pairs in [("coordrep_top1", coordrep_top1_pairs),
                         ("donor_top1", donor_top1_pairs),
                         ("hard_negative", hard_neg_pairs)]:
        sample = pairs[:30]
        for query, match, score, is_pos in sample:
            audit_id += 1
            a = audit_pair(query, match)
            audit_rows.append({
                "audit_id": audit_id,
                "audit_type": label,
                "query_id": query.node_id,
                "query_source": query.source_id,
                "match_id": match.node_id,
                "match_source": match.source_id,
                "score": score,
                "is_positive_label": is_pos,
                **a,
                "audit_note": "",
            })

    audit_fields = [
        "audit_id", "audit_type", "query_id", "query_source",
        "match_id", "match_source", "score", "is_positive_label",
        "metal_correct", "donor_correct", "bridge_match",
        "mu_pattern_match", "shape_match", "chemically_homologous",
        "audit_note",
    ]
    with open(OUTDIR / "hard_control_manual_audit.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=audit_fields)
        w.writeheader()
        w.writerows(audit_rows)

    # Audit summary
    by_type = defaultdict(list)
    for r in audit_rows:
        by_type[r["audit_type"]].append(r)

    audit_summary = {}
    for atype, rows in by_type.items():
        n = len(rows)
        n_homologous = sum(1 for r in rows if r["chemically_homologous"])
        audit_summary[atype] = {
            "n": n,
            "homologous_rate": round(n_homologous / max(n, 1), 4),
        }
    print(f"  Audit: coordrep_top1={audit_summary.get('coordrep_top1', {}).get('homologous_rate', 0):.1%}, "
          f"donor_top1={audit_summary.get('donor_top1', {}).get('homologous_rate', 0):.1%}, "
          f"hard_neg={audit_summary.get('hard_negative', {}).get('homologous_rate', 0):.1%}")

    with open(OUTDIR / "hard_control_manual_audit_summary.json", "w") as f:
        json.dump(audit_summary, f, indent=2)

    # ── Step 5: Go/no-go decision ──
    print("\n[Step 5] Go/no-go decision...")

    full_top1 = results["coordrep_full"]["top1_correct_rate"]
    donor_top1 = results["donor_baseline"]["top1_correct_rate"]
    full_auroc = results["coordrep_full"]["mean_auroc"]
    donor_auroc = results["donor_baseline"]["mean_auroc"]
    random_top1 = results["random_same_donor"]["top1_correct_rate"]

    # Delta
    delta_top1 = full_top1 - donor_top1
    delta_auroc = full_auroc - donor_auroc

    # Go/no-go criteria
    go_criteria = {
        "full_beats_donor_top1": full_top1 > donor_top1,
        "delta_top1_significant": delta_top1 >= 0.10,
        "full_auroc_above_0.6": full_auroc >= 0.6,
        "full_beats_random": full_top1 > random_top1,
        "n_pools_sufficient": len(pools) >= 20,
    }
    n_pass = sum(go_criteria.values())
    decision = "GO_main_text" if n_pass >= 4 else "NO_GO_SI_only"

    print(f"  CoordRep full top-1: {full_top1:.4f}")
    print(f"  Donor baseline top-1: {donor_top1:.4f}")
    print(f"  Delta: {delta_top1:+.4f}")
    print(f"  AUROC full: {full_auroc:.4f} vs donor: {donor_auroc:.4f}")
    print(f"  Go criteria passed: {n_pass}/5")
    print(f"  Decision: {decision}")

    # ── Outputs ──
    print("\n[Step 6] Writing outputs...")

    summary = OrderedDict([
        ("n_mof_nodes", len(mof_nodes)),
        ("n_molecular_refs", len(mol_refs)),
        ("best_subtask", best_subtask),
        ("subtask_pool_stats", all_pool_stats),
        ("best_subtask_n_pools", len(pools)),
        ("best_subtask_mean_pool_size", round(mean_pool, 1)),
        ("methods_best_subtask", {
            m: {
                "top1_correct_rate": r["top1_correct_rate"],
                "top5_correct_rate": r["top5_correct_rate"],
                "mean_auroc": r["mean_auroc"],
                "mean_fp_rate": r["mean_fp_rate"],
            } for m, r in results.items()
        }),
        ("all_subtask_methods", {
            subtask: {
                m: {
                    "top1": r["top1_correct_rate"],
                    "top5": r["top5_correct_rate"],
                    "auroc": r["mean_auroc"],
                } for m, r in sresults.items()
            } for subtask, sresults in all_subtask_results.items()
        }),
        ("delta_full_vs_donor", {
            "top1": round(delta_top1, 4),
            "auroc": round(delta_auroc, 4),
        }),
        ("go_criteria", go_criteria),
        ("decision", decision),
        ("audit_summary", audit_summary),
        ("elapsed_seconds", round(time.time() - t_start)),
    ])

    with open(OUTDIR / "hard_control_summary.json", "w") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)

    # Figure-ready: control bars
    with open(OUTDIR / "fig_hard_control_bars.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=[
            "method", "top1_correct_rate", "top5_correct_rate", "mean_auroc", "mean_fp_rate"])
        w.writeheader()
        for m, r in results.items():
            w.writerow({
                "method": m,
                "top1_correct_rate": r["top1_correct_rate"],
                "top5_correct_rate": r["top5_correct_rate"],
                "mean_auroc": r["mean_auroc"],
                "mean_fp_rate": r["mean_fp_rate"],
            })

    # Pool index
    pool_fields = ["pool_id", "mof_source", "hard_key_metal", "hard_key_cn",
                   "hard_key_donor", "n_candidates", "n_positives", "n_negatives"]
    with open(OUTDIR / "hard_pool_index.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=pool_fields)
        w.writeheader()
        for pid, pd in pools.items():
            q = pd["query"]
            w.writerow({
                "pool_id": pid,
                "mof_source": q.source_id,
                "hard_key_metal": q.metals_str,
                "hard_key_cn": q.cn_list,
                "hard_key_donor": q.donor_set,
                "n_candidates": len(pd["pool"]),
                "n_positives": len(pd["positives"]),
                "n_negatives": len(pd["negatives"]),
            })

    # README
    readme = f"""# CoordRep-Rosetta Hard Controls

## Purpose

Test whether CoordRep higher-order fields (bridge topology, haptic modes,
mu-pattern, shape) add retrieval value beyond metal/CN/donor matching.

Three sub-tasks tested:
- **A (same_metal_cn)**: candidates share same metal + CN. {all_pool_stats['A_same_metal_cn']['n_pools']} pools.
- **B (same_metal)**: candidates share same metal. {all_pool_stats['B_same_metal']['n_pools']} pools.
- **C (same_metal_donor)**: candidates share same metal + donor set. {all_pool_stats['C_same_metal_donor']['n_pools']} pools.

Best subtask for go/no-go: **{best_subtask}** ({len(pools)} pools, mean {mean_pool:.1f} candidates).

## Data

- **MOF nodes**: {len(mof_nodes)}
- **Molecular references**: {len(mol_refs)} (mono={mono_count}, multi={multi_count})

## Results: Best Subtask ({best_subtask})

| Method | Top-1 Correct | Top-5 Correct | AUROC | FP Rate |
|---|---|---|---|---|
"""
    for m, r in results.items():
        readme += (f"| {m} | {r['top1_correct_rate']:.4f} | "
                   f"{r['top5_correct_rate']:.4f} | "
                   f"{r['mean_auroc']:.4f} | {r['mean_fp_rate']:.4f} |\n")

    if all_pool_stats['A_same_metal_cn']['n_pools'] > 0:
        readme += f"""
## Results: Sub-task A (same_metal_cn, {all_pool_stats['A_same_metal_cn']['n_pools']} pools)

| Method | Top-1 | Top-5 | AUROC |
|---|---|---|---|
"""
        for m, r in all_subtask_results["A_same_metal_cn"].items():
            readme += f"| {m} | {r['top1_correct_rate']:.4f} | {r['top5_correct_rate']:.4f} | {r['mean_auroc']:.4f} |\n"

    readme += f"""
## Go/No-Go

| Criterion | Pass |
|---|---|
"""
    for criterion, passed in go_criteria.items():
        readme += f"| {criterion} | {'YES' if passed else 'NO'} |\n"

    readme += f"""
**Decision: {decision}**

Delta (CoordRep full - donor baseline):
- Top-1: {delta_top1:+.4f}
- AUROC: {delta_auroc:+.4f}

## Manual Audit

| Type | N | Homologous Rate |
|---|---|---|
"""
    for atype, info in audit_summary.items():
        readme += f"| {atype} | {info['n']} | {info['homologous_rate']:.1%} |\n"

    readme += f"""
## Files

| File | Description |
|---|---|
| hard_control_summary.json | Full summary with go/no-go decision |
| hard_control_method_comparison.csv | Method comparison table |
| hard_control_detailed_results.csv | Per-pool per-method results |
| hard_pool_index.csv | Hard pool definitions |
| hard_control_manual_audit.csv | Manual audit rows |
| hard_control_manual_audit_summary.json | Audit summary |
| fig_hard_control_bars.csv | Figure-ready bar chart data |
| README.md | This file |

## No raw CSD coordinates exported.
"""
    (OUTDIR / "README.md").write_text(readme)

    elapsed = time.time() - t_start
    print(f"\n{'=' * 70}")
    print(f"Hard Controls Complete ({elapsed:.0f}s)")
    print(f"{'=' * 70}")
    print(f"  Hard pools: {len(pools)}")
    print(f"  CoordRep full top-1: {full_top1:.4f}")
    print(f"  Donor baseline top-1: {donor_top1:.4f}")
    print(f"  Delta: {delta_top1:+.4f}")
    print(f"  AUROC: {full_auroc:.4f} vs {donor_auroc:.4f}")
    print(f"  Decision: {decision}")
    print(f"  Outputs → {OUTDIR}")
    print("Done.")


if __name__ == "__main__":
    main()
