#!/usr/bin/env python3
"""Reproducible nuisance-invariance audit for the v2-beta upgrade."""

from __future__ import annotations

import copy
import itertools
import json
import random
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from coordrep.v2beta.canonicalize import (  # noqa: E402
    _serialize_haptic_content,
    _serialize_multi_content,
    canonicalize_haptic,
    canonicalize_multi,
)
from coordrep.v2beta.core import (  # noqa: E402
    HapticRecord,
    MetalCenter,
    MetalEdge,
    MultiMetalRecord,
)
from coordrep.v2beta.validate import validate_haptic, validate_multi  # noqa: E402
from scripts.generate_v2beta_extension import (  # noqa: E402
    build_all_haptic,
    build_all_multi,
)

N_VARIANTS = 100
OUTDIR = ROOT / "revision_results" / "canonical_invariance_upgrade"


def rename_shuffle_multi(record: MultiMetalRecord, seed: int) -> MultiMetalRecord:
    rng = random.Random(seed)
    rec = copy.deepcopy(record)
    rng.shuffle(rec.metals)
    rng.shuffle(rec.metal_edges)
    rng.shuffle(rec.sites)
    rng.shuffle(rec.ligands)

    metal_map = {
        metal.label: f"metal_source_{seed}_{idx}"
        for idx, metal in enumerate(rec.metals)
    }
    ligand_map = {
        ligand.label: f"ligand_source_{seed}_{idx}"
        for idx, ligand in enumerate(rec.ligands)
    }
    site_map = {
        site.label: f"site_source_{seed}_{idx}"
        for idx, site in enumerate(rec.sites)
    }

    for metal in rec.metals:
        old = metal.label
        metal.label = metal_map[old]
        metal.local_sites = [site_map.get(label, label) for label in metal.local_sites]
    for edge in rec.metal_edges:
        edge.m1 = metal_map[edge.m1]
        edge.m2 = metal_map[edge.m2]
        edge.bridges = [site_map.get(label, label) for label in edge.bridges]
    for ligand in rec.ligands:
        old = ligand.label
        ligand.label = ligand_map[old]
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
            f"donor_source_{seed}_{site_idx}_{atom_idx}"
            for atom_idx in range(len(combined))
        ]
        if site.centroid_label:
            site.centroid_label = f"centroid_source_{seed}_{site_idx}"
    return rec


def rename_rotate_haptic(record: HapticRecord, seed: int) -> HapticRecord:
    rng = random.Random(seed)
    rec = copy.deepcopy(record)
    rng.shuffle(rec.sites)
    rng.shuffle(rec.ligands)

    old_metal = rec.metal.label
    rec.metal.label = f"metal_source_{seed}"
    ligand_map = {
        ligand.label: f"ligand_source_{seed}_{idx}"
        for idx, ligand in enumerate(rec.ligands)
    }
    for ligand in rec.ligands:
        ligand.label = ligand_map[ligand.label]

    for site_idx, site in enumerate(rec.sites):
        site.label = f"site_source_{seed}_{site_idx}"
        site.ligand_label = ligand_map[site.ligand_label]
        site.target_metals = [
            rec.metal.label if target == old_metal else target
            for target in site.target_metals
        ]

        keys = site.meta.get("donor_canonical_keys", [""] * len(site.donor_elements))
        combined = list(zip(site.donor_elements, keys))
        if combined:
            shift = seed % len(combined)
            combined = combined[shift:] + combined[:shift]
            if seed % 2:
                combined.reverse()
        site.donor_elements = [element for element, _ in combined]
        site.meta["donor_canonical_keys"] = [key for _, key in combined]
        site.donor_atoms = [
            f"eta_source_{seed}_{site_idx}_{atom_idx}"
            for atom_idx in range(len(combined))
        ]
        if site.centroid_label:
            site.centroid_label = f"centroid_source_{seed}_{site_idx}"
    return rec


def bare_cu4_cycle() -> MultiMetalRecord:
    metals = [
        MetalCenter(label=f"X{idx}", element="Cu", oxidation=2, dcount=9)
        for idx in range(4)
    ]
    edges = [
        MetalEdge(f"X{idx}", f"X{(idx + 1) % 4}", "bridged", mm_bond="no")
        for idx in range(4)
    ]
    return MultiMetalRecord("cu4", "SYNTHETIC_CU4", metals, edges, [], [])


def permuted_cu4_records():
    base = bare_cu4_cycle()
    for permutation in itertools.permutations(range(4)):
        rec = copy.deepcopy(base)
        original = [f"X{idx}" for idx in range(4)]
        rec.metals = [rec.metals[idx] for idx in permutation]
        mapping = {
            original[old_idx]: f"source_{new_idx}"
            for new_idx, old_idx in enumerate(permutation)
        }
        for metal in rec.metals:
            metal.label = mapping[metal.label]
        for edge in rec.metal_edges:
            edge.m1 = mapping[edge.m1]
            edge.m2 = mapping[edge.m2]
        yield rec


def legacy_cu4_edge_code(record: MultiMetalRecord) -> tuple:
    """Reproduce the old stable-sort tie behaviour without importing it."""
    mapping = {metal.label: f"M{idx + 1}" for idx, metal in enumerate(record.metals)}
    return tuple(sorted(
        tuple(sorted((mapping[edge.m1], mapping[edge.m2])))
        for edge in record.metal_edges
    ))


def main() -> int:
    start = time.perf_counter()
    multi_records = build_all_multi()
    haptic_records = build_all_haptic()

    failures = []
    multi_trials = 0
    for record_idx, record in enumerate(multi_records):
        expected = _serialize_multi_content(record)
        expected_id = record.identity.L0_GlobalState
        for variant_idx in range(N_VARIANTS):
            seed = 100_000 * record_idx + variant_idx
            variant = canonicalize_multi(rename_shuffle_multi(record, seed))
            multi_trials += 1
            if (
                _serialize_multi_content(variant) != expected
                or variant.identity.L0_GlobalState != expected_id
            ):
                failures.append({
                    "record_type": "multi",
                    "case_id": record.case_id,
                    "seed": seed,
                })

    haptic_trials = 0
    for record_idx, record in enumerate(haptic_records):
        expected = _serialize_haptic_content(record)
        expected_id = record.identity.L0_HapticState
        for variant_idx in range(N_VARIANTS):
            seed = 200_000 * record_idx + variant_idx
            variant = canonicalize_haptic(rename_rotate_haptic(record, seed))
            haptic_trials += 1
            if (
                _serialize_haptic_content(variant) != expected
                or variant.identity.L0_HapticState != expected_id
            ):
                failures.append({
                    "record_type": "haptic",
                    "case_id": record.case_id,
                    "seed": seed,
                })

    cu4_variants = list(permuted_cu4_records())
    legacy_cu4_distinct = len({legacy_cu4_edge_code(record) for record in cu4_variants})
    upgraded_cu4_outputs = {
        _serialize_multi_content(canonicalize_multi(record))
        for record in cu4_variants
    }

    cycle = canonicalize_multi(bare_cu4_cycle())
    path = bare_cu4_cycle()
    path.metal_edges = path.metal_edges[:3]
    path = canonicalize_multi(path)

    multi_validation = [validate_multi(record).all_passed for record in multi_records]
    haptic_validation = [validate_haptic(record).all_passed for record in haptic_records]

    summary = {
        "algorithm": "exact individualization-refinement over metal vertices; label-free site/ligand incidence leaf code",
        "scope": "geometry-resolved abstract coordination record graph",
        "variants_per_curated_record": N_VARIANTS,
        "curated_multinuclear_records": len(multi_records),
        "curated_multinuclear_nuisance_trials": multi_trials,
        "curated_multinuclear_passed": multi_trials - sum(f["record_type"] == "multi" for f in failures),
        "curated_haptic_records": len(haptic_records),
        "curated_haptic_nuisance_trials": haptic_trials,
        "curated_haptic_passed": haptic_trials - sum(f["record_type"] == "haptic" for f in failures),
        "cu4_input_permutations": len(cu4_variants),
        "cu4_legacy_distinct_serializations": legacy_cu4_distinct,
        "cu4_upgraded_distinct_serializations": len(upgraded_cu4_outputs),
        "cycle_path_noncollision": cycle.identity.L0_GlobalState != path.identity.L0_GlobalState,
        "multi_validator_passed": sum(multi_validation),
        "multi_validator_total": len(multi_validation),
        "haptic_validator_passed": sum(haptic_validation),
        "haptic_validator_total": len(haptic_validation),
        "failures": failures,
        "all_passed": not failures
        and len(upgraded_cu4_outputs) == 1
        and legacy_cu4_distinct == 3
        and all(multi_validation)
        and all(haptic_validation),
        "elapsed_seconds": round(time.perf_counter() - start, 6),
        "claim_boundary": {
            "supported": "canonical identity of all fields encoded in the geometry-resolved molecular coordination record graph",
            "not_supported_without_adapter_upgrade": "full ligand chemical identity when the CSD adapter emits element-only or [*] placeholder ligands",
            "not_a_species_identifier": "different experimentally or computationally resolved geometries may intentionally yield different L0 states",
        },
    }

    OUTDIR.mkdir(parents=True, exist_ok=True)
    output = OUTDIR / "canonical_invariance_audit_summary.json"
    output.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0 if summary["all_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
