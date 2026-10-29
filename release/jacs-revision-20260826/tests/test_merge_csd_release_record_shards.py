from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import pytest

from coordrep.audit.records import record_digest
from coordrep.audit.metal_policy import METAL_POLICY_ID, METAL_POLICY_SHA256
from scripts import analyze_csd_release_audit as analyze
from scripts import merge_csd_release_record_shards as merge


def _write_jsonl(path: Path, rows):
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for row in rows:
            handle.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")


def _record(
    source_index: int,
    *,
    scope: str = "mono_eta1",
    scope_flags=(),
    external_center: bool = False,
    external_shared_donor: bool = False,
    donor_kind: str | None = None,
):
    scope_flags = tuple(sorted({*scope_flags, "metal_block_d"}))
    metal_site_count = 2 if scope.startswith("multi_") else 1
    donor_groups = []
    if donor_kind is not None:
        donor_groups.append({
            "kind": donor_kind,
            "hapticity": 5 if donor_kind == "pi_fragment" else 1,
            "bridge_degree": 1,
            "donor_elements": ["C"],
        })
    body = {
        "schema_version": "CoordRep-Record/3",
        "scope": scope,
        "scope_flags": list(scope_flags),
        "metal_sites": [{
            "element": "Fe",
            "record_level": "topology",
            "resolution": "partial",
            "canonical_status": "exact_label_free",
            "geometry_kind": "distance-angle-profile",
            "shape_best": "",
            "cn_site": 0,
            "cn_atom": 0,
            "issues": [],
        } for _ in range(metal_site_count)],
        "donor_groups": donor_groups,
        "incidences": [],
        "metal_relations": [],
        "external_center_relations": ([{
            "metal": "M1",
            "partner_element": "Xx",
            "partner_topology_key": "external-center",
            "bond_type": "Single",
            "distance": 2.5,
            "image_delta": [0, 0, 0],
            "image_delta_resolved": True,
        }] if external_center else []),
        "external_donor_relations": ([{
            "donor_group": "G1",
            "donor_ordinal": 0,
            "partner_element": "Xx",
            "partner_topology_key": "external-shared-donor",
            "bond_type": "Single",
            "distance": 2.1,
            "image_delta": [0, 0, 0],
            "image_delta_resolved": True,
        }] if external_shared_donor else []),
        "donor_group_relations": [],
        "donor_pair_relations": [],
    }
    return {
        **body,
        "source_index": source_index,
        "refcode": f"FAKE{source_index}",
        "expected_metal_sites": metal_site_count,
        "entry_resolution": "EMITTED_PARTIAL",
        "issue_codes": [],
        "record_checksum_sha256": record_digest(body),
    }


def _outcome(record):
    metal_site_count = record["expected_metal_sites"]
    return {
        "source_index": record["source_index"],
        "refcode": record["refcode"],
        "status": record["entry_resolution"],
        "scope": record["scope"],
        "scope_flags": record["scope_flags"],
        "metal_elements": ["Fe"],
        "metal_blocks": ["d"],
        "metal_sites": metal_site_count,
        "structural_metal_sites": metal_site_count,
        "audit_only_metal_sites": 0,
        "all_metal_sites_structural": True,
        "donor_groups": len(record["donor_groups"]),
        "incidences": len(record["incidences"]),
        "structural_record": True,
        "schema_valid": True,
        "record_checksum_sha256": record["record_checksum_sha256"],
        "issue_codes": [],
    }


def _make_two_shards(root: Path):
    source_fingerprint = {"adapter": "abc"}
    runtime = {"python": "test"}
    for shard_id in range(2):
        directory = root / f"shard_{shard_id:04d}_of_0002"
        directory.mkdir(parents=True)
        records = [_record(2 * shard_id), _record(2 * shard_id + 1)]
        outcomes = [_outcome(record) for record in records]
        _write_jsonl(directory / "RECORDS_INTERNAL.jsonl", records)
        _write_jsonl(directory / "ENTRY_OUTCOMES.jsonl", outcomes)
        _write_jsonl(directory / "FAILURES_INTERNAL.jsonl", [])
        summary = {
            "config": {
                "database": {"entry_count": 4},
                "database_sha256": "a" * 64,
                "runtime": runtime,
                "sources": source_fingerprint,
                "source_release": "test-release",
                "protocol_id": "test-protocol",
                "metal_policy_id": METAL_POLICY_ID,
                "metal_policy_sha256": METAL_POLICY_SHA256,
                "work_kind": "database_range",
                "start": 2 * shard_id,
                "stop": 2 * shard_id + 2,
                "shard_id": shard_id,
                "num_shards": 2,
            },
            "processed_to": 2 * shard_id + 2,
            "counters": {
                "inputs_attempted": 2,
                "target_entries": 2,
                "records_emitted": 2,
                "structural_records_emitted": 2,
                "audit_only_records_emitted": 0,
                "EMITTED_PARTIAL": 2,
            },
        }
        (directory / "SUMMARY.json").write_text(
            json.dumps(summary), encoding="utf-8"
        )


def _invoke(monkeypatch, root: Path, output: Path):
    census = root / "METAL_DOMAIN_CENSUS_PUBLIC.json"
    census.write_text(json.dumps({
        "database_entries": 4,
        "database_sha256": "a" * 64,
        "source_release": "test-release",
        "metal_policy_id": METAL_POLICY_ID,
        "metal_policy_sha256": METAL_POLICY_SHA256,
        "metal_containing_entries": 5,
        "metal_containing_3d_entries": 4,
        "metal_containing_without_3d_entries": 1,
    }), encoding="utf-8")
    monkeypatch.setattr(sys, "argv", [
        "merge",
        "--root", str(root),
        "--num-shards", "2",
        "--output-dir", str(output),
        "--expected-database-entries", "4",
        "--expected-target-entries", "4",
        "--metal-domain-census", str(census),
    ])
    merge.main()


def test_merge_proves_pairing_ranges_and_digest(tmp_path, monkeypatch):
    root = tmp_path / "shards"
    _make_two_shards(root)
    output = tmp_path / "merged"
    _invoke(monkeypatch, root, output)
    summary = json.loads(
        (output / "SUMMARY_PUBLIC_AGGREGATE.json").read_text(encoding="utf-8")
    )
    assert summary["target_entries"] == 4
    assert summary["metal_containing_entries"] == 5
    assert summary["metal_containing_without_3d_entries"] == 1
    assert summary["structural_record_coverage"] == 1.0
    assert summary["structural_metal_site_coverage"] == 1.0


def test_merge_rejects_tampered_record_body(tmp_path, monkeypatch):
    root = tmp_path / "shards"
    _make_two_shards(root)
    path = root / "shard_0000_of_0002" / "RECORDS_INTERNAL.jsonl"
    rows = [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines()]
    rows[0]["scope"] = "multi_eta1"
    _write_jsonl(path, rows)
    with pytest.raises(SystemExit, match="record digest mismatch"):
        _invoke(monkeypatch, root, tmp_path / "merged")


def test_merge_rejects_tampered_metal_domain_metadata(tmp_path, monkeypatch):
    root = tmp_path / "shards"
    _make_two_shards(root)
    path = root / "shard_0000_of_0002" / "ENTRY_OUTCOMES.jsonl"
    rows = [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines()]
    rows[0]["metal_elements"] = ["Al"]
    rows[0]["metal_blocks"] = ["p"]
    _write_jsonl(path, rows)
    with pytest.raises(SystemExit, match="record/outcome statistics mismatch"):
        _invoke(monkeypatch, root, tmp_path / "merged")


def test_claim_ready_analysis_is_hash_locked_census(tmp_path, monkeypatch):
    root = tmp_path / "shards"
    output = tmp_path / "merged"
    _make_two_shards(root)
    _invoke(monkeypatch, root, output)
    claim_output = tmp_path / "CLAIM_READY.json"
    monkeypatch.setattr(sys, "argv", [
        "analyze",
        "--records", str(output / "RECORDS_INTERNAL.jsonl"),
        "--outcomes", str(output / "ENTRY_OUTCOMES_INTERNAL.jsonl"),
        "--merged-summary", str(output / "SUMMARY_PUBLIC_AGGREGATE.json"),
        "--output", str(claim_output),
    ])
    analyze.main()
    claim = json.loads(claim_output.read_text(encoding="utf-8"))
    assert claim["statistical_frame"]["kind"] == "frozen_release_census"
    assert claim["coverage"]["entry_all_metal_sites_structural"]["fraction"] == 1.0
    assert "wilson_95_low" not in claim["coverage"]["entry_any_structural"]
    assert all(claim["accounting_checks"].values())


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    return digest.hexdigest()


def test_claim_ready_external_relation_partition_is_mutually_exclusive(
    tmp_path, monkeypatch
):
    records = [
        # A: one in-domain metal, no external relation, nonpolymeric -> strict.
        _record(
            0,
            scope="mono_haptic",
            scope_flags=("haptic",),
            donor_kind="pi_fragment",
        ),
        # B: one in-domain metal, direct external co-center -> external.
        _record(
            1,
            scope="mono_haptic",
            scope_flags=("haptic", "heterometal_co_center"),
            external_center=True,
            donor_kind="pi_fragment",
        ),
        # C: a pi candidate with an external shared donor is not haptic.
        _record(
            2,
            scope="mono_haptic",
            scope_flags=(
                "ambiguous_pi_candidate",
                "external_metal_shared_donor",
            ),
            external_shared_donor=True,
            donor_kind="candidate_atom_set",
        ),
        # D: polymeric takes precedence over an external shared donor.
        _record(
            3,
            scope="mono_haptic",
            scope_flags=(
                "external_metal_shared_donor",
                "haptic",
                "polymeric",
            ),
            external_shared_donor=True,
            donor_kind="pi_fragment",
        ),
        # E: multiple in-domain metal centers -> multiple.
        _record(4, scope="multi_eta1", scope_flags=("multimetal",)),
        # F: direct and shared external relations count once in the union.
        _record(
            5,
            scope_flags=(
                "external_metal_shared_donor",
                "heterometal_co_center",
            ),
            external_center=True,
            external_shared_donor=True,
        ),
    ]
    outcomes = [_outcome(record) for record in records]
    records_path = tmp_path / "RECORDS_INTERNAL.jsonl"
    outcomes_path = tmp_path / "ENTRY_OUTCOMES_INTERNAL.jsonl"
    merged_summary_path = tmp_path / "SUMMARY_PUBLIC_AGGREGATE.json"
    claim_output = tmp_path / "CLAIM_READY.json"
    _write_jsonl(records_path, records)
    _write_jsonl(outcomes_path, outcomes)
    total_sites = sum(record["expected_metal_sites"] for record in records)
    merged_summary = {
        "source_release": "test-release",
        "protocol_id": "test-protocol",
        "database_sha256": "a" * 64,
        "target_entries": len(records),
        "metal_policy_id": METAL_POLICY_ID,
        "metal_policy_sha256": METAL_POLICY_SHA256,
        "records_emitted": len(records),
        "failed_outcomes": 0,
        "structural_records_emitted": len(records),
        "entries_all_metal_sites_structural": len(records),
        "schema_valid_records": len(records),
        "metal_sites_accounted": total_sites,
        "structural_metal_sites": total_sites,
        "audit_only_metal_sites": 0,
        "failed_metal_sites": 0,
        "artifacts": {
            "records_internal_sha256": _sha256(records_path),
            "outcomes_internal_sha256": _sha256(outcomes_path),
        },
    }
    merged_summary_path.write_text(
        json.dumps(merged_summary), encoding="utf-8"
    )
    monkeypatch.setattr(sys, "argv", [
        "analyze",
        "--records", str(records_path),
        "--outcomes", str(outcomes_path),
        "--merged-summary", str(merged_summary_path),
        "--output", str(claim_output),
    ])
    analyze.main()
    claim = json.loads(claim_output.read_text(encoding="utf-8"))
    scope = claim["chemical_scope"]
    assert scope["strict_nonpolymeric_single_metal_center_entries"] == 1
    assert (
        scope["nonpolymeric_multiple_metal_centers_or_external_relation_entries"]
        == 4
    )
    assert scope["polymeric_entries"] == 1
    assert scope["multiple_metal_centers_or_external_relation_entries"] == 5
    assert scope["single_metal_center_entries"] == 5
    assert scope["multiple_metal_center_entries"] == 1
    assert scope["metal_block_entry_counts"] == {"d": 6}
    assert scope["external_metal_center_relation_entries"] == 2
    assert scope["external_metal_shared_donor_entries"] == 3
    assert scope["haptic_entries"] == 3
    assert scope["ambiguous_pi_candidate_entries"] == 1
    assert scope[
        "strict_nonpolymeric_single_metal_center_haptic_entries"
    ] == 1
    assert scope[
        "nonpolymeric_multiple_metal_centers_or_external_relation_haptic_entries"
    ] == 1
    assert scope["polymeric_haptic_entries"] == 1
    assert (
        scope["strict_nonpolymeric_single_metal_center_haptic_entries"]
        + scope[
            "nonpolymeric_multiple_metal_centers_or_external_relation_haptic_entries"
        ]
        + scope["polymeric_haptic_entries"]
        == scope["haptic_entries"]
    )
    assert (
        scope["strict_nonpolymeric_single_metal_center_entries"]
        + scope[
            "nonpolymeric_multiple_metal_centers_or_external_relation_entries"
        ]
        + scope["polymeric_entries"]
        == claim["emitted_entries"]
    )
    assert all(claim["accounting_checks"].values())
