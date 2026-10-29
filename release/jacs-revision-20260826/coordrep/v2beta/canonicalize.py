"""Canonicalization for the CoordRep v2-beta molecular extensions.

The original v2-beta implementation used a fixed number of
Weisfeiler--Lehman refinement rounds and then relied on Python's stable sort.
That is not a canonical labelling algorithm: unresolved colour classes retain
the input order.  In particular, a chemically homogeneous Cu4 cycle produced
three different metal graphs under the 24 possible input permutations.

This module uses exact individualization--refinement for the (at most twelve)
metal vertices.  Every leaf is scored with a label-free code of the complete
metal/site/ligand incidence record and the lexicographically smallest code is
selected.  Sites and ligands are then ordered from label-free chemical and
incidence keys.  Raw source atom labels are retained only as provenance in
``site.meta``; they are never used to define a CoordRep identity.

The resulting guarantee is deliberately scoped: it is an exact canonical
labelling of the information present in ``MultiMetalRecord`` or
``HapticRecord``.  It cannot recover chemical distinctions that an upstream
adapter did not encode (for example, an adapter that emits ``[*]`` instead of
an atom-mapped ligand graph).
"""

from __future__ import annotations

import hashlib
from collections import Counter, defaultdict
from typing import Dict, Hashable, Iterable, List, Optional, Sequence, Tuple

from .core import (
    CoordinationSite,
    HapticIdentityKeys,
    HapticRecord,
    Ligand,
    MetalCenter,
    MetalEdge,
    MultiIdentityKeys,
    MultiMetalRecord,
)


def _short_hash(s: str, n: int = 12) -> str:
    return hashlib.sha256(s.encode()).hexdigest()[:n]


def _stable(value) -> str:
    """Return a comparable, deterministic scalar representation."""
    if value is None:
        return "<none>"
    if isinstance(value, float):
        return format(value, ".12g")
    if isinstance(value, (list, tuple)):
        return "(" + ",".join(_stable(v) for v in value) + ")"
    if isinstance(value, dict):
        return "{" + ",".join(
            f"{_stable(k)}:{_stable(v)}" for k, v in sorted(value.items(), key=lambda kv: _stable(kv[0]))
        ) + "}"
    return str(value)


def _metal_intrinsic(m: MetalCenter) -> Tuple[str, ...]:
    """Metal fields that are serialized, excluding arbitrary record labels."""
    return (
        m.element,
        m.ox_str,
        _stable(m.dcount),
        _stable(m.cn_site),
        _stable(m.cn_atom),
        _stable(m.eta_sum),
        m.local_shape_best,
        m.local_shape_class,
        _stable(m.local_shape_delta),
        m.local_shape_vals,
    )


def _donor_descriptors(site: CoordinationSite) -> Tuple[Tuple[str, str], ...]:
    """Label-free donor descriptors for a coordination site.

    An upstream adapter may provide ``donor_canonical_keys`` (for example,
    atom-graph refinement signatures).  When absent, the current record only
    contains an unordered element multiset, so that is the strongest invariant
    descriptor that can be used without inventing chemistry.
    """
    keys = site.meta.get("donor_canonical_keys")
    if not isinstance(keys, (list, tuple)) or len(keys) != len(site.donor_elements):
        keys = [""] * len(site.donor_elements)
    return tuple(sorted((str(e), _stable(k)) for e, k in zip(site.donor_elements, keys)))


def _site_intrinsic(site: CoordinationSite) -> Tuple:
    return (
        site.site_type,
        _donor_descriptors(site),
        _stable(site.eta),
        _stable(site.mu),
        site.mode,
        "centroid" if site.centroid_label else "",
    )


def _ligand_intrinsic(ligand: Ligand) -> Tuple[str, ...]:
    return (
        ligand.smiles,
        _stable(ligand.dent),
        _stable(ligand.charge),
    )


def _edge_intrinsic(edge: MetalEdge) -> Tuple[str, ...]:
    return (
        edge.relation,
        "" if edge.d_mm is None else f"{edge.d_mm:.2f}",
        edge.mm_bond,
    )


def _rank_signatures(signatures: Sequence[Hashable]) -> List[int]:
    """Replace arbitrary signatures by deterministic consecutive colours."""
    unique = sorted(set(signatures), key=_stable)
    rank = {signature: idx for idx, signature in enumerate(unique)}
    return [rank[signature] for signature in signatures]


def _validated_metal_index(metals: Sequence[MetalCenter]) -> Dict[str, int]:
    labels = [m.label for m in metals]
    if len(labels) != len(set(labels)):
        raise ValueError("metal labels must be unique before canonicalization")
    return {label: idx for idx, label in enumerate(labels)}


def _site_code(
    site: CoordinationSite,
    metal_position: Dict[str, int],
) -> Tuple:
    try:
        targets = tuple(sorted(metal_position[label] for label in site.target_metals))
    except KeyError as exc:
        raise ValueError(f"site {site.label!r} targets unknown metal {exc.args[0]!r}") from exc
    return (_site_intrinsic(site), targets)


def _record_code_for_metal_order(
    metals: Sequence[MetalCenter],
    edges: Sequence[MetalEdge],
    sites: Sequence[CoordinationSite],
    ligands: Sequence[Ligand],
    order: Sequence[int],
) -> Tuple:
    """Exact label-free code used to choose a metal canonical labelling."""
    if sorted(order) != list(range(len(metals))):
        raise ValueError("metal order is not a permutation")
    position = {metals[old_idx].label: new_idx for new_idx, old_idx in enumerate(order)}

    metal_code = tuple(_metal_intrinsic(metals[idx]) for idx in order)

    edge_code = []
    for edge in edges:
        if edge.m1 not in position or edge.m2 not in position:
            raise ValueError(f"metal edge targets unknown metal: {edge.m1!r}, {edge.m2!r}")
        a, b = sorted((position[edge.m1], position[edge.m2]))
        edge_code.append((a, b, _edge_intrinsic(edge)))

    sites_by_ligand: Dict[str, List[CoordinationSite]] = defaultdict(list)
    for site in sites:
        sites_by_ligand[site.ligand_label].append(site)

    ligand_labels = {ligand.label for ligand in ligands}
    ligand_code = []
    for ligand in ligands:
        attached = tuple(sorted(_site_code(site, position) for site in sites_by_ligand.get(ligand.label, [])))
        ligand_code.append((_ligand_intrinsic(ligand), attached))

    # Invalid/orphan sites remain visible to the code instead of disappearing.
    orphan_code = tuple(sorted(
        _site_code(site, position) for site in sites if site.ligand_label not in ligand_labels
    ))
    return (
        metal_code,
        tuple(sorted(edge_code)),
        tuple(sorted(ligand_code)),
        orphan_code,
    )


def _metal_refinement_context(
    metals: Sequence[MetalCenter],
    edges: Sequence[MetalEdge],
    sites: Sequence[CoordinationSite],
    ligands: Sequence[Ligand],
):
    metal_index = _validated_metal_index(metals)
    ligand_by_label = {ligand.label: ligand for ligand in ligands}
    sites_by_ligand: Dict[str, List[CoordinationSite]] = defaultdict(list)
    for site in sites:
        sites_by_ligand[site.ligand_label].append(site)

    incident_edges: List[List[Tuple[MetalEdge, int]]] = [[] for _ in metals]
    for edge in edges:
        if edge.m1 not in metal_index or edge.m2 not in metal_index:
            raise ValueError(f"metal edge targets unknown metal: {edge.m1!r}, {edge.m2!r}")
        i, j = metal_index[edge.m1], metal_index[edge.m2]
        incident_edges[i].append((edge, j))
        incident_edges[j].append((edge, i))

    incident_sites: List[List[CoordinationSite]] = [[] for _ in metals]
    for site in sites:
        for target in site.target_metals:
            if target not in metal_index:
                raise ValueError(f"site {site.label!r} targets unknown metal {target!r}")
            incident_sites[metal_index[target]].append(site)

    return metal_index, ligand_by_label, sites_by_ligand, incident_edges, incident_sites


def _refine_metal_colours(
    metals: Sequence[MetalCenter],
    edges: Sequence[MetalEdge],
    sites: Sequence[CoordinationSite],
    ligands: Sequence[Ligand],
    seeds: Sequence[Hashable],
) -> List[int]:
    """Equitable refinement of metal colours with site/ligand incidence."""
    (
        metal_index,
        ligand_by_label,
        sites_by_ligand,
        incident_edges,
        incident_sites,
    ) = _metal_refinement_context(metals, edges, sites, ligands)
    colours = _rank_signatures(list(seeds))

    for _ in range(max(2, len(metals) + 2)):
        fingerprints = []
        for idx, metal in enumerate(metals):
            edge_context = tuple(sorted(
                (_edge_intrinsic(edge), colours[other])
                for edge, other in incident_edges[idx]
            ))

            site_context = []
            for site in incident_sites[idx]:
                target_colours = tuple(sorted(colours[metal_index[t]] for t in site.target_metals))
                ligand = ligand_by_label.get(site.ligand_label)
                ligand_intrinsic = _ligand_intrinsic(ligand) if ligand is not None else ("<orphan>", "", "")
                ligand_site_context = tuple(sorted(
                    (
                        _site_intrinsic(other_site),
                        tuple(sorted(colours[metal_index[t]] for t in other_site.target_metals)),
                    )
                    for other_site in sites_by_ligand.get(site.ligand_label, [])
                ))
                site_context.append((
                    _site_intrinsic(site),
                    target_colours,
                    ligand_intrinsic,
                    ligand_site_context,
                ))

            fingerprints.append((
                colours[idx],
                _metal_intrinsic(metal),
                edge_context,
                tuple(sorted(site_context)),
            ))
        new_colours = _rank_signatures(fingerprints)
        if new_colours == colours:
            return colours
        colours = new_colours
    return colours


def _direct_swap_is_automorphism(
    metals: Sequence[MetalCenter],
    edges: Sequence[MetalEdge],
    sites: Sequence[CoordinationSite],
    ligands: Sequence[Ligand],
    i: int,
    j: int,
) -> bool:
    """Exact safe pruning test for true transposition automorphisms."""
    identity = list(range(len(metals)))
    swapped = list(identity)
    swapped[i], swapped[j] = swapped[j], swapped[i]
    return _record_code_for_metal_order(metals, edges, sites, ligands, identity) == (
        _record_code_for_metal_order(metals, edges, sites, ligands, swapped)
    )


def _exact_canonical_metal_indices(
    metals: Sequence[MetalCenter],
    edges: Sequence[MetalEdge],
    sites: Sequence[CoordinationSite],
    ligands: Sequence[Ligand],
) -> List[int]:
    """Return an exact canonical metal order by individualization--refinement."""
    n_metals = len(metals)
    if n_metals <= 1:
        return list(range(n_metals))

    root_seeds = [_metal_intrinsic(metal) for metal in metals]
    best_code = None
    best_order = None

    def visit(seeds: Sequence[Hashable]) -> None:
        nonlocal best_code, best_order
        colours = _refine_metal_colours(metals, edges, sites, ligands, seeds)
        cells: Dict[int, List[int]] = defaultdict(list)
        for idx, colour in enumerate(colours):
            cells[colour].append(idx)

        unresolved = [(colour, members) for colour, members in sorted(cells.items()) if len(members) > 1]
        if not unresolved:
            order = sorted(range(n_metals), key=lambda idx: colours[idx])
            code = _record_code_for_metal_order(metals, edges, sites, ligands, order)
            if best_code is None or code < best_code:
                best_code = code
                best_order = order
            return

        # The smallest unresolved cell controls branching.  True transposition
        # automorphisms are represented once; this avoids factorial work for a
        # complete set of indistinguishable metal twins without unsafe WL-only
        # pruning.
        _, cell = min(unresolved, key=lambda item: (len(item[1]), item[0]))
        representatives: List[int] = []
        for candidate in cell:
            if any(
                _direct_swap_is_automorphism(
                    metals, edges, sites, ligands, candidate, prior
                )
                for prior in representatives
            ):
                continue
            representatives.append(candidate)

        for candidate in representatives:
            individualized = [
                (colours[idx], "selected" if idx == candidate else "unselected")
                for idx in range(n_metals)
            ]
            visit(individualized)

    visit(root_seeds)
    if best_order is None:  # defensive; every finite search must reach a leaf
        raise RuntimeError("metal canonical-labelling search produced no leaf")
    return list(best_order)


def _wl_metal_signatures(
    metals: List[MetalCenter],
    edges: List[MetalEdge],
    sites: List[CoordinationSite],
    n_iter: int = 3,
) -> Dict[str, str]:
    """Compatibility helper returning stable refined metal signatures.

    The ``n_iter`` argument is retained for callers of the former private API;
    refinement now runs to stability.  These signatures may contain ties and
    are not themselves used as a canonical labelling.
    """
    del n_iter
    colours = _refine_metal_colours(
        metals, edges, sites, [], [_metal_intrinsic(metal) for metal in metals]
    )
    return {metal.label: f"c{colours[idx]}" for idx, metal in enumerate(metals)}


def _canonical_metal_order(
    metals: List[MetalCenter],
    edges: List[MetalEdge],
    sites: List[CoordinationSite],
    ligands: Optional[List[Ligand]] = None,
) -> List[MetalCenter]:
    order = _exact_canonical_metal_indices(metals, edges, sites, ligands or [])
    return [metals[idx] for idx in order]


def _relabel_metals(
    metals: List[MetalCenter],
    edges: List[MetalEdge],
    sites: List[CoordinationSite],
    ligands: List[Ligand],
) -> tuple:
    """Apply exact metal canonical labels and update every metal reference."""
    ordered = _canonical_metal_order(metals, edges, sites, ligands)
    mapping = {metal.label: f"M{idx + 1}" for idx, metal in enumerate(ordered)}

    for metal in ordered:
        metal.label = mapping[metal.label]
    for edge in edges:
        edge.m1 = mapping[edge.m1]
        edge.m2 = mapping[edge.m2]
        if edge.m1 > edge.m2:
            edge.m1, edge.m2 = edge.m2, edge.m1
    for site in sites:
        site.target_metals = sorted(mapping[target] for target in site.target_metals)
    return ordered, edges, sites, ligands


def _normalize_site_donors(site: CoordinationSite) -> None:
    """Remove source atom naming from the canonical record.

    The source labels remain in ``meta['source_donor_labels']`` for audit and
    back-mapping.  Canonical donor tokens are local to the site (C1, C2, ...
    after serialization), which is sufficient because a coordination site is
    defined as an unordered atom set in v2-beta.
    """
    site.meta.setdefault("source_donor_labels", list(site.donor_atoms))
    keys = site.meta.get("donor_canonical_keys")
    if not isinstance(keys, (list, tuple)) or len(keys) != len(site.donor_elements):
        keys = [""] * len(site.donor_elements)

    combined = sorted(
        zip(site.donor_elements, keys),
        key=lambda item: (str(item[0]), _stable(item[1])),
    )
    site.donor_elements = [str(element) for element, _ in combined]
    site.meta["donor_canonical_keys"] = [_stable(key) for _, key in combined]

    per_element = Counter()
    canonical_labels = []
    for element in site.donor_elements:
        per_element[element] += 1
        canonical_labels.append(str(per_element[element]))
    site.donor_atoms = canonical_labels
    if site.centroid_label:
        site.meta.setdefault("source_centroid_label", site.centroid_label)
        site.centroid_label = "site-centroid"


def _ligand_sort_key(ligand: Ligand, sites: List[CoordinationSite]) -> tuple:
    ligand_sites = tuple(sorted(
        site.site_signature() for site in sites if site.ligand_label == ligand.label
    ))
    return (_ligand_intrinsic(ligand), ligand_sites)


def _canonical_ligand_order(
    ligands: List[Ligand],
    sites: List[CoordinationSite],
) -> List[Ligand]:
    """Order ligands without altering which sites belong to each ligand."""
    ordered = sorted(ligands, key=lambda ligand: _ligand_sort_key(ligand, sites))
    mapping = {ligand.label: f"L{idx + 1}" for idx, ligand in enumerate(ordered)}
    for ligand in ordered:
        ligand.label = mapping[ligand.label]
    for site in sites:
        if site.ligand_label in mapping:
            site.ligand_label = mapping[site.ligand_label]
    return ordered


def _canonical_site_order(sites: List[CoordinationSite]) -> List[CoordinationSite]:
    for site in sites:
        _normalize_site_donors(site)
    return sorted(sites, key=lambda site: (site.site_signature(), site.ligand_label))


def _serialize_multi_content(rec: MultiMetalRecord) -> str:
    """Serialize without the ID block for stable hashing."""
    from .serialize import serialize_multi

    full = serialize_multi(rec)
    idx = full.find("[ID:")
    return full[:idx].strip() if idx >= 0 else full


def _serialize_haptic_content(rec: HapticRecord) -> str:
    """Serialize without the ID block for stable hashing."""
    from .serialize import serialize_haptic

    full = serialize_haptic(rec)
    idx = full.find("[ID:")
    return full[:idx].strip() if idx >= 0 else full


def _generate_multi_identity(
    rec: MultiMetalRecord,
    serialized: str,
) -> MultiIdentityKeys:
    del serialized
    content = _serialize_multi_content(rec)
    l0 = _short_hash(content, 16)

    shape_parts = [
        f"{metal.element}|{metal.local_shape_best}|CN{metal.cn_site}"
        for metal in rec.metals
    ]
    l1 = _short_hash(";".join(sorted(shape_parts)), 16)

    edge_parts = []
    for edge in rec.metal_edges:
        a, b = sorted([edge.m1, edge.m2])
        edge_parts.append(f"{a}-{b}|{edge.relation}")
    metal_parts = [metal.element for metal in rec.metals]
    l2_raw = ";".join(sorted(metal_parts)) + "|" + ";".join(sorted(edge_parts))
    l2 = _short_hash(l2_raw, 16)

    ligand_parts = sorted(ligand.smiles for ligand in rec.ligands)
    l3_raw = ";".join(sorted(metal_parts)) + "|" + ";".join(ligand_parts)
    l3 = _short_hash(l3_raw, 16)

    local_ids = {}
    for metal in rec.metals:
        local_sites = sorted(
            site.site_signature()
            for site in rec.sites
            if metal.label in site.target_metals
        )
        local_raw = (
            f"{metal.element}|{metal.ox_str}|{metal.local_shape_best}|{local_sites}"
        )
        local_ids[metal.label] = _short_hash(local_raw, 12)

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
    del serialized
    content = _serialize_haptic_content(rec)
    l0 = _short_hash(content, 16)

    metal = rec.metal
    shape_raw = (
        f"{metal.element}|{metal.ox_str}|{metal.local_shape_best}|CNsite:{metal.cn_site}"
    )
    l1 = _short_hash(shape_raw, 16)

    site_parts = sorted(site.site_signature() for site in rec.sites)
    l2 = _short_hash(f"{metal.element}|{site_parts}", 16)

    ligand_parts = sorted(ligand.smiles for ligand in rec.ligands)
    l3 = _short_hash(f"{metal.element}|{';'.join(ligand_parts)}", 16)

    return HapticIdentityKeys(
        L0_HapticState=l0,
        L1_HapticShape=l1,
        L2_SiteTopo=l2,
        L3_Connectivity=l3,
    )


def _canonicalize_multi_once(rec: MultiMetalRecord) -> str:
    metals, edges, sites, ligands = _relabel_metals(
        rec.metals, rec.metal_edges, rec.sites, rec.ligands
    )
    rec.metals = list(metals)
    rec.metal_edges = sorted(edges, key=lambda edge: edge.canonical_key())

    # Donor normalization must precede ligand/site key construction.
    for site in sites:
        _normalize_site_donors(site)
    rec.ligands = _canonical_ligand_order(ligands, sites)

    old_site_labels = {id(site): site.label for site in sites}
    rec.sites = _canonical_site_order(sites)
    site_mapping = {}
    for idx, site in enumerate(rec.sites):
        old_label = old_site_labels[id(site)]
        site.label = f"S{idx + 1}"
        site_mapping[old_label] = site.label

    for edge in rec.metal_edges:
        edge.bridges = sorted(site_mapping.get(label, label) for label in edge.bridges)
    for ligand in rec.ligands:
        ligand.sites = [site.label for site in rec.sites if site.ligand_label == ligand.label]
    for metal in rec.metals:
        metal.local_sites = [
            site.label for site in rec.sites if metal.label in site.target_metals
        ]
    return _serialize_multi_content(rec)


def canonicalize_multi(rec: MultiMetalRecord) -> MultiMetalRecord:
    """Canonicalize a multinuclear abstract record to a fixed point."""
    previous = ""
    for _ in range(5):
        current = _canonicalize_multi_once(rec)
        if current == previous:
            break
        previous = current

    from .serialize import serialize_multi

    rec.identity = _generate_multi_identity(rec, serialize_multi(rec))
    return rec


def _canonicalize_haptic_once(rec: HapticRecord) -> str:
    rec.metal.label = "M1"
    for site in rec.sites:
        site.target_metals = ["M1"]
        _normalize_site_donors(site)
    rec.ligands = _canonical_ligand_order(rec.ligands, rec.sites)

    rec.sites = _canonical_site_order(rec.sites)
    for idx, site in enumerate(rec.sites):
        site.label = f"S{idx + 1}"
    for ligand in rec.ligands:
        ligand.sites = [site.label for site in rec.sites if site.ligand_label == ligand.label]
    rec.metal.local_sites = [site.label for site in rec.sites]
    return _serialize_haptic_content(rec)


def canonicalize_haptic(rec: HapticRecord) -> HapticRecord:
    """Canonicalize a haptic abstract record to a fixed point."""
    previous = ""
    for _ in range(5):
        current = _canonicalize_haptic_once(rec)
        if current == previous:
            break
        previous = current

    from .serialize import serialize_haptic

    rec.identity = _generate_haptic_identity(rec, serialize_haptic(rec))
    return rec
