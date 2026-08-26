"""Strict nuisance-invariance tests for the v2-beta molecular extensions."""

from __future__ import annotations

import copy
import itertools
import random

from coordrep.v2beta.canonicalize import (
    _serialize_haptic_content,
    _serialize_multi_content,
    canonicalize_haptic,
    canonicalize_multi,
)
from coordrep.v2beta.core import (
    CoordinationSite,
    HapticRecord,
    Ligand,
    MetalCenter,
    MetalEdge,
    MultiMetalRecord,
)


def _bare_cu4_cycle() -> MultiMetalRecord:
    metals = [
        MetalCenter(label=f"X{idx}", element="Cu", oxidation=2, dcount=9)
        for idx in range(4)
    ]
    edges = [
        MetalEdge(
            m1=f"X{idx}",
            m2=f"X{(idx + 1) % 4}",
            relation="bridged",
            mm_bond="no",
        )
        for idx in range(4)
    ]
    return MultiMetalRecord(
        case_id="cu4_cycle",
        refcode="SYNTHETIC_CU4",
        metals=metals,
        metal_edges=edges,
        sites=[],
        ligands=[],
    )


def _symmetric_bridged_cu4_cycle() -> MultiMetalRecord:
    rec = _bare_cu4_cycle()
    rec.sites = []
    rec.ligands = []
    for idx in range(4):
        bridge_ligand = Ligand(
            label=f"LB{idx}", smiles="[O-]", dent=1, charge=-1
        )
        terminal_ligand = Ligand(
            label=f"LT{idx}", smiles="N", dent=1, charge=0
        )
        rec.ligands.extend([bridge_ligand, terminal_ligand])
        rec.sites.extend([
            CoordinationSite(
                label=f"SB{idx}",
                site_type="atom",
                ligand_label=bridge_ligand.label,
                donor_atoms=[f"arbitrary_O_{idx}"],
                donor_elements=["O"],
                eta=1,
                mu=2,
                target_metals=[f"X{idx}", f"X{(idx + 1) % 4}"],
                mode="mu2-O",
            ),
            CoordinationSite(
                label=f"ST{idx}",
                site_type="atom",
                ligand_label=terminal_ligand.label,
                donor_atoms=[f"arbitrary_N_{idx}"],
                donor_elements=["N"],
                eta=1,
                mu=1,
                target_metals=[f"X{idx}"],
            ),
        ])
    for metal in rec.metals:
        local = [site for site in rec.sites if metal.label in site.target_metals]
        metal.cn_site = len(local)
        metal.cn_atom = sum(site.eta for site in local)
        metal.eta_sum = metal.cn_atom
        metal.local_sites = [site.label for site in local]
    return rec


def _rename_and_shuffle_multi(
    record: MultiMetalRecord,
    seed: int,
) -> MultiMetalRecord:
    rng = random.Random(seed)
    rec = copy.deepcopy(record)
    rng.shuffle(rec.metals)
    rng.shuffle(rec.metal_edges)
    rng.shuffle(rec.sites)
    rng.shuffle(rec.ligands)

    metal_map = {
        metal.label: f"source_metal_{seed}_{idx}"
        for idx, metal in enumerate(rec.metals)
    }
    ligand_map = {
        ligand.label: f"source_ligand_{seed}_{idx}"
        for idx, ligand in enumerate(rec.ligands)
    }
    site_map = {
        site.label: f"source_site_{seed}_{idx}"
        for idx, site in enumerate(rec.sites)
    }

    for metal in rec.metals:
        old_label = metal.label
        metal.label = metal_map[old_label]
        metal.local_sites = [site_map.get(label, label) for label in metal.local_sites]
    for edge in rec.metal_edges:
        edge.m1 = metal_map[edge.m1]
        edge.m2 = metal_map[edge.m2]
        edge.bridges = [site_map.get(label, label) for label in edge.bridges]
    for ligand in rec.ligands:
        old_label = ligand.label
        ligand.label = ligand_map[old_label]
        ligand.sites = [site_map.get(label, label) for label in ligand.sites]
    for site_idx, site in enumerate(rec.sites):
        site.label = site_map[site.label]
        site.ligand_label = ligand_map[site.ligand_label]
        site.target_metals = [metal_map[label] for label in site.target_metals]

        keys = site.meta.get("donor_canonical_keys", [""] * len(site.donor_elements))
        combined = list(zip(site.donor_elements, keys))
        rng.shuffle(combined)
        site.donor_elements = [element for element, _ in combined]
        site.meta["donor_canonical_keys"] = [key for _, key in combined]
        site.donor_atoms = [
            f"renamed_atom_{seed}_{site_idx}_{atom_idx}"
            for atom_idx in range(len(combined))
        ]
        if site.centroid_label:
            site.centroid_label = f"renamed_centroid_{seed}_{site_idx}"
    return rec


def _ferrocene_like_record() -> HapticRecord:
    metal = MetalCenter(
        label="Fe_source",
        element="Fe",
        oxidation=2,
        dcount=6,
        cn_site=2,
        cn_atom=10,
        eta_sum=10,
    )
    ligands = [
        Ligand(label="ring_a", smiles="[cH-]1cccc1", dent=1, charge=-1),
        Ligand(label="ring_b", smiles="[cH-]1cccc1", dent=1, charge=-1),
    ]
    sites = [
        CoordinationSite(
            label="upper_ring",
            site_type="pi_fragment",
            ligand_label="ring_a",
            donor_atoms=[f"top_{idx}" for idx in range(5)],
            donor_elements=["C"] * 5,
            eta=5,
            mu=1,
            target_metals=["Fe_source"],
            mode="Cp",
            centroid_label="upper_centroid",
        ),
        CoordinationSite(
            label="lower_ring",
            site_type="pi_fragment",
            ligand_label="ring_b",
            donor_atoms=[f"bottom_{idx}" for idx in range(5)],
            donor_elements=["C"] * 5,
            eta=5,
            mu=1,
            target_metals=["Fe_source"],
            mode="Cp",
            centroid_label="lower_centroid",
        ),
    ]
    return HapticRecord(
        case_id="ferrocene_like",
        refcode="SYNTHETIC_FECP2",
        metal=metal,
        sites=sites,
        ligands=ligands,
    )


def _rename_rotate_reverse_haptic(record: HapticRecord, seed: int) -> HapticRecord:
    rng = random.Random(seed)
    rec = copy.deepcopy(record)
    rng.shuffle(rec.sites)
    rng.shuffle(rec.ligands)
    old_metal = rec.metal.label
    rec.metal.label = f"renamed_metal_{seed}"

    ligand_map = {
        ligand.label: f"renamed_ligand_{seed}_{idx}"
        for idx, ligand in enumerate(rec.ligands)
    }
    for ligand in rec.ligands:
        ligand.label = ligand_map[ligand.label]

    for site_idx, site in enumerate(rec.sites):
        site.label = f"renamed_site_{seed}_{site_idx}"
        site.ligand_label = ligand_map[site.ligand_label]
        site.target_metals = [
            rec.metal.label if target == old_metal else target
            for target in site.target_metals
        ]
        combined = list(zip(site.donor_atoms, site.donor_elements))
        shift = seed % len(combined)
        combined = combined[shift:] + combined[:shift]
        if seed % 2:
            combined.reverse()
        site.donor_elements = [element for _, element in combined]
        site.donor_atoms = [
            f"renamed_ring_atom_{seed}_{site_idx}_{idx}"
            for idx in range(len(combined))
        ]
        site.centroid_label = f"renamed_centroid_{seed}_{site_idx}"
    return rec


def test_cu4_24_metal_permutations_collapse_from_three_to_one():
    base = _bare_cu4_cycle()
    outputs = set()
    ids = set()
    for permutation in itertools.permutations(range(4)):
        rec = copy.deepcopy(base)
        old_labels = [f"X{idx}" for idx in range(4)]
        rec.metals = [rec.metals[idx] for idx in permutation]
        metal_map = {
            old_labels[old_idx]: f"input_{new_idx}"
            for new_idx, old_idx in enumerate(permutation)
        }
        for metal in rec.metals:
            metal.label = metal_map[metal.label]
        for edge in rec.metal_edges:
            edge.m1 = metal_map[edge.m1]
            edge.m2 = metal_map[edge.m2]
        canonicalize_multi(rec)
        outputs.add(_serialize_multi_content(rec))
        ids.add(rec.identity.L0_GlobalState)
    assert len(outputs) == 1
    assert len(ids) == 1


def test_symmetric_bridges_survive_full_container_shuffle_and_label_rename():
    base = canonicalize_multi(_symmetric_bridged_cu4_cycle())
    expected = _serialize_multi_content(base)
    expected_id = base.identity.L0_GlobalState
    for seed in range(40):
        variant = _rename_and_shuffle_multi(base, seed)
        canonicalize_multi(variant)
        assert _serialize_multi_content(variant) == expected
        assert variant.identity.L0_GlobalState == expected_id


def test_eta5_rotations_reversals_and_label_renames_collapse():
    base = canonicalize_haptic(_ferrocene_like_record())
    expected = _serialize_haptic_content(base)
    expected_id = base.identity.L0_HapticState
    for seed in range(20):
        variant = _rename_rotate_reverse_haptic(base, seed)
        canonicalize_haptic(variant)
        assert _serialize_haptic_content(variant) == expected
        assert variant.identity.L0_HapticState == expected_id


def test_nonisomorphic_cu4_path_and_cycle_do_not_collide():
    cycle = canonicalize_multi(_bare_cu4_cycle())
    path = _bare_cu4_cycle()
    path.metal_edges = path.metal_edges[:3]
    path = canonicalize_multi(path)
    assert _serialize_multi_content(path) != _serialize_multi_content(cycle)
    assert path.identity.L0_GlobalState != cycle.identity.L0_GlobalState


def test_twelve_identical_complete_graph_metals_do_not_factorially_expand():
    labels = [f"Z{idx}" for idx in range(12)]
    metals = [MetalCenter(label=label, element="Cu") for label in labels]
    edges = [
        MetalEdge(a, b, relation="contact", mm_bond="contact")
        for a, b in itertools.combinations(labels, 2)
    ]
    rec = MultiMetalRecord("cu12", "SYNTHETIC_CU12", metals, edges, [], [])
    expected = _serialize_multi_content(canonicalize_multi(rec))
    variant = _rename_and_shuffle_multi(rec, seed=12012)
    assert _serialize_multi_content(canonicalize_multi(variant)) == expected
