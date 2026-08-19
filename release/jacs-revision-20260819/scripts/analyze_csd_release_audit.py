#!/usr/bin/env python3
"""Create claim-ready aggregate statistics from a merged CSD audit.

Only aggregate counts are written.  No CSD coordinates, atom labels, or
individual licensed records are included in the public summary.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Dict


BASE = Path(__file__).resolve().parent.parent
if str(BASE) not in sys.path:
    sys.path.insert(0, str(BASE))

from coordrep.audit.records import ENTRY_RESOLUTIONS, record_digest


FAILED_STATUSES = {"FAILED_INPUT", "FAILED_PIPELINE"}
OUTPUT_SCOPES = {"mono_eta1", "mono_haptic", "multi_eta1", "multi_haptic"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--records", type=Path, required=True)
    parser.add_argument("--outcomes", type=Path, required=True)
    parser.add_argument("--merged-summary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _atomic_json(path: Path, payload: Dict[str, Any]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    os.replace(temporary, path)


def _rate(successes: int, total: int) -> Dict[str, Any]:
    return {
        "numerator": successes,
        "denominator": total,
        "fraction": successes / total if total else None,
        "percent": 100.0 * successes / total if total else None,
    }


def main() -> None:
    args = parse_args()
    merged_summary = json.loads(args.merged_summary.read_text(encoding="utf-8"))
    expected_artifacts = merged_summary.get("artifacts", {})
    observed_records_sha256 = _sha256(args.records)
    observed_outcomes_sha256 = _sha256(args.outcomes)
    if observed_records_sha256 != expected_artifacts.get("records_internal_sha256"):
        raise SystemExit("records file does not match the merged-summary SHA-256")
    if observed_outcomes_sha256 != expected_artifacts.get("outcomes_internal_sha256"):
        raise SystemExit("outcomes file does not match the merged-summary SHA-256")
    expected_entries = int(merged_summary["transition_metal_entries"])
    expected_records = int(merged_summary["records_emitted"])

    outcome_counts = Counter()
    per_scope = defaultdict(Counter)
    emitted_refcodes = set()
    outcome_refcodes = set()
    source_indices = set()
    emitted_outcomes: Dict[str, Dict[str, Any]] = {}
    with args.outcomes.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            row = json.loads(line)
            refcode = str(row["refcode"])
            source_index = int(row["source_index"])
            if refcode in outcome_refcodes:
                raise SystemExit(f"duplicate outcome refcode: {refcode}")
            if source_index in source_indices:
                raise SystemExit(f"duplicate outcome source index: {source_index}")
            outcome_refcodes.add(refcode)
            source_indices.add(source_index)
            outcome_counts["entries"] += 1
            status = str(row["status"])
            if status not in ENTRY_RESOLUTIONS | FAILED_STATUSES:
                raise SystemExit(f"unknown terminal status for {refcode}: {status}")
            outcome_counts[f"status:{status}"] += 1
            sites = int(row.get("metal_sites", 0))
            structural_sites = int(row.get("structural_metal_sites", 0))
            audit_only_sites = int(row.get("audit_only_metal_sites", 0))
            failed_sites = int(row.get("failed_metal_sites", 0))
            if min(sites, structural_sites, audit_only_sites, failed_sites) < 0:
                raise SystemExit(f"negative site count for {refcode}")
            if structural_sites + audit_only_sites + failed_sites != sites:
                raise SystemExit(f"metal-site fields do not close for {refcode}")
            if status in FAILED_STATUSES:
                if sites < 1 or failed_sites != sites:
                    raise SystemExit(f"failed site accounting mismatch for {refcode}")
                outcome_counts["failed_entries"] += 1
                outcome_counts[f"failure:{row.get('issue_code', 'UNKNOWN')}"] += 1
                outcome_counts["metal_sites"] += sites
                outcome_counts["failed_metal_sites"] += failed_sites
                continue
            if row.get("schema_valid") is not True:
                raise SystemExit(f"emitted outcome is not schema-valid for {refcode}")
            if failed_sites:
                raise SystemExit(f"emitted outcome has failed metal sites for {refcode}")
            any_structural = bool(row["structural_record"])
            all_structural = bool(row["all_metal_sites_structural"])
            if any_structural != (structural_sites > 0):
                raise SystemExit(f"structural-record flag mismatch for {refcode}")
            if all_structural != (audit_only_sites == 0):
                raise SystemExit(f"all-sites flag mismatch for {refcode}")
            if (status == "EMITTED_AUDIT_ONLY") != (structural_sites == 0):
                raise SystemExit(f"audit-only terminal status mismatch for {refcode}")
            emitted_refcodes.add(refcode)
            emitted_outcomes[refcode] = row
            scope = str(row["scope"])
            if scope not in OUTPUT_SCOPES:
                raise SystemExit(f"unknown output scope for {refcode}: {scope}")
            flags = set(row.get("scope_flags", []))
            outcome_counts["emitted_entries"] += 1
            outcome_counts["entries_any_structural"] += int(any_structural)
            outcome_counts["entries_all_sites_structural"] += int(all_structural)
            outcome_counts["metal_sites"] += sites
            outcome_counts["structural_metal_sites"] += structural_sites
            outcome_counts["audit_only_metal_sites"] += audit_only_sites
            outcome_counts[f"scope:{scope}"] += 1
            one_dblock = scope.startswith("mono_")
            multiple_dblock = scope.startswith("multi_")
            is_polymeric = "polymeric" in flags
            has_external_metal_cocenter = "heterometal_co_center" in flags
            has_external_metal_shared_donor = (
                "external_metal_shared_donor" in flags
            )
            has_any_external_metal_relation = (
                has_external_metal_cocenter
                or has_external_metal_shared_donor
            )
            has_confirmed_haptic = "haptic" in flags
            is_strict_nonpolymeric_single = (
                one_dblock
                and not is_polymeric
                and not has_any_external_metal_relation
            )
            is_nonpolymeric_multiple_or_external = (
                not is_polymeric
                and (multiple_dblock or has_any_external_metal_relation)
            )
            outcome_counts["one_in_domain_dblock_center_entries"] += int(one_dblock)
            outcome_counts["multiple_in_domain_dblock_center_entries"] += int(
                multiple_dblock
            )
            outcome_counts["strict_nonpolymeric_single_metal_center_entries"] += int(
                is_strict_nonpolymeric_single
            )
            outcome_counts[
                "nonpolymeric_multiple_metal_centers_or_external_relation_entries"
            ] += int(is_nonpolymeric_multiple_or_external)
            outcome_counts[
                "multiple_metal_centers_or_external_relation_entries"
            ] += int(multiple_dblock or has_any_external_metal_relation)
            outcome_counts[
                "strict_nonpolymeric_single_metal_center_haptic_entries"
            ] += int(is_strict_nonpolymeric_single and has_confirmed_haptic)
            outcome_counts[
                "nonpolymeric_multiple_metal_centers_or_external_relation_haptic_entries"
            ] += int(
                is_nonpolymeric_multiple_or_external and has_confirmed_haptic
            )
            outcome_counts["polymeric_haptic_entries"] += int(
                is_polymeric and has_confirmed_haptic
            )
            for flag in flags:
                outcome_counts[f"flag:{flag}"] += 1
            for issue in row.get("issue_codes", []):
                outcome_counts[f"issue:{issue}"] += 1
            bucket = per_scope[scope]
            bucket["entries"] += 1
            bucket["entries_any_structural"] += int(any_structural)
            bucket["entries_all_sites_structural"] += int(all_structural)
            bucket["metal_sites"] += sites
            bucket["structural_metal_sites"] += structural_sites

    record_counts = Counter()
    seen_refcodes = set()
    with args.records.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            record = json.loads(line)
            refcode = record["refcode"]
            if refcode in seen_refcodes:
                raise SystemExit(f"duplicate record refcode: {refcode}")
            seen_refcodes.add(refcode)
            outcome = emitted_outcomes.get(refcode)
            if outcome is None:
                raise SystemExit(f"record has no emitted terminal outcome: {refcode}")
            if record.get("entry_resolution") != outcome.get("status"):
                raise SystemExit(f"record/outcome status mismatch for {refcode}")
            if record.get("record_checksum_sha256") != outcome.get(
                "record_checksum_sha256"
            ):
                raise SystemExit(f"record/outcome checksum mismatch for {refcode}")
            structural_body = {
                key: record[key]
                for key in (
                    "schema_version",
                    "scope",
                    "scope_flags",
                    "metal_sites",
                    "donor_groups",
                    "incidences",
                    "metal_relations",
                    "external_center_relations",
                    "external_donor_relations",
                    "donor_group_relations",
                    "donor_pair_relations",
                )
            }
            if record_digest(structural_body) != record.get("record_checksum_sha256"):
                raise SystemExit(f"record structural digest mismatch for {refcode}")
            structural_sites = sum(
                site.get("record_level") != "audit_only"
                for site in record["metal_sites"]
            )
            expected_outcome_fields = {
                "scope": record["scope"],
                "scope_flags": record["scope_flags"],
                "metal_sites": len(record["metal_sites"]),
                "structural_metal_sites": structural_sites,
                "audit_only_metal_sites": len(record["metal_sites"]) - structural_sites,
                "structural_record": structural_sites > 0,
                "all_metal_sites_structural": structural_sites == len(record["metal_sites"]),
                "donor_groups": len(record["donor_groups"]),
                "incidences": len(record["incidences"]),
            }
            if any(outcome.get(key) != value for key, value in expected_outcome_fields.items()):
                raise SystemExit(f"record/outcome field mismatch for {refcode}")
            record_counts["records"] += 1
            site_statuses = {
                site["canonical_status"] for site in record["metal_sites"]
            }
            record_counts["entries_exact_label_free"] += int(
                site_statuses == {"exact_label_free"}
            )
            record_counts["entries_source_order_noncanonical"] += int(
                "source_order_noncanonical" in site_statuses
            )
            has_bridge = False
            has_haptic = False
            has_pi_candidate = False
            has_metal_hydrogen = False
            has_nonzero_incidence_translation = False
            has_nonzero_metal_relation_translation = False
            has_nonzero_external_center_translation = False
            has_nonzero_external_donor_translation = False
            has_unresolved_translation_edge = False
            for site in record["metal_sites"]:
                record_counts[f"site_record_level:{site['record_level']}"] += 1
                record_counts[f"site_resolution:{site['resolution']}"] += 1
                record_counts[f"site_canonical:{site['canonical_status']}"] += 1
                record_counts[f"geometry_kind:{site['geometry_kind']}"] += 1
                if site.get("shape_best"):
                    record_counts[f"shape_best:{site['shape_best']}"] += 1
                record_counts[f"cn_site:{site['cn_site']}"] += 1
                record_counts[f"cn_atom:{site['cn_atom']}"] += 1
            for group in record["donor_groups"]:
                has_bridge = has_bridge or int(group["bridge_degree"]) > 1
                has_haptic = has_haptic or (
                    group["kind"] == "pi_fragment" and int(group["hapticity"]) > 1
                )
                has_pi_candidate = has_pi_candidate or (
                    group["kind"] == "candidate_atom_set"
                )
                has_metal_hydrogen = has_metal_hydrogen or "H" in group["donor_elements"]
            for incidence in record["incidences"]:
                has_unresolved_translation_edge = (
                    has_unresolved_translation_edge
                    or not bool(incidence.get("image_delta_resolved", True))
                )
                has_nonzero_incidence_translation = (
                    has_nonzero_incidence_translation
                    or (
                        bool(incidence.get("image_delta_resolved", True))
                        and any(
                            int(value) != 0 for value in incidence["image_delta"]
                        )
                    )
                )
            for metal_relation in record["metal_relations"]:
                has_unresolved_translation_edge = (
                    has_unresolved_translation_edge
                    or not bool(metal_relation.get("image_delta_resolved", True))
                )
                has_nonzero_metal_relation_translation = (
                    has_nonzero_metal_relation_translation
                    or (
                        bool(metal_relation.get("image_delta_resolved", True))
                        and any(
                            int(value) != 0
                            for value in metal_relation.get("image_delta", (0, 0, 0))
                        )
                    )
                )
            for external_relation in record.get("external_center_relations", []):
                has_unresolved_translation_edge = (
                    has_unresolved_translation_edge
                    or not bool(external_relation.get("image_delta_resolved", True))
                )
                has_nonzero_external_center_translation = (
                    has_nonzero_external_center_translation
                    or (
                        bool(external_relation.get("image_delta_resolved", True))
                        and any(
                            int(value) != 0
                            for value in external_relation.get("image_delta", (0, 0, 0))
                        )
                    )
                )
            for external_relation in record.get("external_donor_relations", []):
                has_unresolved_translation_edge = (
                    has_unresolved_translation_edge
                    or not bool(external_relation.get("image_delta_resolved", True))
                )
                has_nonzero_external_donor_translation = (
                    has_nonzero_external_donor_translation
                    or (
                        bool(external_relation.get("image_delta_resolved", True))
                        and any(
                            int(value) != 0
                            for value in external_relation.get("image_delta", (0, 0, 0))
                        )
                    )
                )
            record_counts["entries_with_bridge"] += int(has_bridge)
            record_counts["entries_with_haptic_site"] += int(has_haptic)
            record_counts["entries_with_haptic_site_nonambiguous"] += int(
                has_haptic and record.get("entry_resolution") != "EMITTED_AMBIGUOUS"
            )
            record_counts["entries_with_haptic_site_disorder_ambiguous"] += int(
                has_haptic and record.get("entry_resolution") == "EMITTED_AMBIGUOUS"
            )
            record_counts["entries_with_ambiguous_pi_candidate"] += int(
                has_pi_candidate
            )
            record_counts["entries_with_any_shared_donor"] += int(
                has_bridge or bool(record.get("external_donor_relations", []))
            )
            record_counts["entries_with_metal_bound_hydrogen"] += int(has_metal_hydrogen)
            record_counts["entries_with_nonzero_translation_incidence"] += int(
                has_nonzero_incidence_translation
            )
            record_counts["entries_with_nonzero_translation_edge"] += int(
                has_nonzero_incidence_translation
                or has_nonzero_metal_relation_translation
                or has_nonzero_external_center_translation
                or has_nonzero_external_donor_translation
            )
            record_counts["entries_with_unresolved_translation_edge"] += int(
                has_unresolved_translation_edge
            )
            record_counts["entries_with_direct_metal_bond"] += int(any(
                relation["relation"] == "direct_metal_bond"
                for relation in record["metal_relations"]
            ))
            record_counts["entries_with_external_metal_center"] += int(
                bool(record.get("external_center_relations", []))
            )
            record_counts["entries_with_external_metal_shared_donor"] += int(
                bool(record.get("external_donor_relations", []))
            )
            for relation in record.get("donor_pair_relations", []):
                record_counts[f"angle_class:{relation['angle_class']}"] += 1

    if seen_refcodes != emitted_refcodes:
        raise SystemExit("record and emitted-outcome refcode sets differ")

    total_entries = outcome_counts["entries"]
    total_sites = outcome_counts["metal_sites"]
    if total_entries != expected_entries:
        raise SystemExit(
            f"outcome rows are {total_entries}, merged summary requires {expected_entries}"
        )
    if record_counts["records"] != expected_records:
        raise SystemExit(
            f"record rows are {record_counts['records']}, merged summary requires "
            f"{expected_records}"
        )
    if outcome_counts["emitted_entries"] != expected_records:
        raise SystemExit("emitted outcome count does not match merged record count")
    merged_count_expectations = {
        "failed_outcomes": outcome_counts["failed_entries"],
        "structural_records_emitted": outcome_counts["entries_any_structural"],
        "entries_all_metal_sites_structural": outcome_counts[
            "entries_all_sites_structural"
        ],
        "schema_valid_records": outcome_counts["emitted_entries"],
        "metal_sites_accounted": total_sites,
        "structural_metal_sites": outcome_counts["structural_metal_sites"],
        "audit_only_metal_sites": outcome_counts["audit_only_metal_sites"],
        "failed_metal_sites": outcome_counts["failed_metal_sites"],
    }
    for field, observed in merged_count_expectations.items():
        if int(merged_summary.get(field, -1)) != observed:
            raise SystemExit(
                f"merged-summary field {field} is inconsistent: "
                f"{merged_summary.get(field)!r} != {observed}"
            )
    if (
        outcome_counts["structural_metal_sites"]
        + outcome_counts["audit_only_metal_sites"]
        + outcome_counts["failed_metal_sites"]
        != total_sites
    ):
        raise SystemExit("aggregate metal-site fields do not close")
    if (
        record_counts["entries_exact_label_free"]
        + record_counts["entries_source_order_noncanonical"]
        != record_counts["records"]
    ):
        raise SystemExit("canonical-status entry partition does not close")
    if (
        outcome_counts["one_in_domain_dblock_center_entries"]
        + outcome_counts["multiple_in_domain_dblock_center_entries"]
        != outcome_counts["emitted_entries"]
    ):
        raise SystemExit("one/multiple in-domain d-block-center partition does not close")
    if (
        outcome_counts["strict_nonpolymeric_single_metal_center_entries"]
        + outcome_counts[
            "nonpolymeric_multiple_metal_centers_or_external_relation_entries"
        ]
        + outcome_counts["flag:polymeric"]
        != outcome_counts["emitted_entries"]
    ):
        raise SystemExit("nonpolymeric mono/multinuclear and polymeric partition does not close")
    if (
        outcome_counts[
            "strict_nonpolymeric_single_metal_center_haptic_entries"
        ]
        + outcome_counts[
            "nonpolymeric_multiple_metal_centers_or_external_relation_haptic_entries"
        ]
        + outcome_counts["polymeric_haptic_entries"]
        != record_counts["entries_with_haptic_site"]
    ):
        raise SystemExit("confirmed-haptic main-class partition does not close")
    summary = {
        "statistical_frame": {
            "kind": "frozen_release_census",
            "source_release": merged_summary.get("source_release"),
            "protocol_id": merged_summary.get("protocol_id"),
            "note": (
                "Fractions are exact descriptive proportions for the frozen "
                "release; no binomial confidence interval or population "
                "generalization is implied."
            ),
        },
        "provenance": {
            "merged_summary_sha256": _sha256(args.merged_summary),
            "analysis_script_sha256": _sha256(Path(__file__).resolve()),
            "records_internal_sha256": observed_records_sha256,
            "outcomes_internal_sha256": observed_outcomes_sha256,
            "database_sha256": merged_summary.get("database_sha256"),
        },
        "definitions": {
            "entry_any_structural": (
                "At least one metal site has a native donor incidence or a "
                "native in-domain or external metal-center relation."
            ),
            "entry_all_metal_sites_structural": (
                "Every transition-metal site in the entry is structural rather "
                "than audit-only."
            ),
            "scope_coverage": (
                "Conditional completeness within output-assigned emitted scopes; "
                "not an independently stratified source-population coverage rate."
            ),
            "bridged_entries": (
                "Emitted records containing a donor group shared by more than "
                "one distinct serialized metal site."
            ),
            "haptic_entries": (
                "Emitted records containing a reconstructed pi_fragment with "
                "hapticity greater than one."
            ),
            "haptic_main_class_partition": (
                "Confirmed haptic entries, identified by the haptic scope flag, "
                "are partitioned among strict nonpolymeric single-center, "
                "nonpolymeric multiple/external-relation, and polymeric records. "
                "Ambiguous pi candidates are excluded."
            ),
            "strict_single_and_multiple_metal_center_entries": (
                "A strict nonpolymeric single-metal-center entry has exactly one "
                "in-domain d-block center and no direct or shared-donor relation "
                "to an external metal. "
                "The multiple-center category includes either multiple in-domain "
                "d-block centers or a direct or shared-donor relation to an "
                "external metal. The two "
                "nonpolymeric classes and the polymeric class are mutually exclusive."
            ),
            "nonzero_translation_edge_entries": (
                "Emitted records containing at least one resolved nonzero image "
                "delta on a metal-donor, in-domain metal-metal, or external "
                "metal-center edge."
            ),
        },
        "transition_metal_entries": total_entries,
        "emitted_entries": outcome_counts["emitted_entries"],
        "failed_entries": outcome_counts["failed_entries"],
        "coverage": {
            "entry_any_structural": _rate(
                outcome_counts["entries_any_structural"], total_entries
            ),
            "entry_all_metal_sites_structural": _rate(
                outcome_counts["entries_all_sites_structural"], total_entries
            ),
            "metal_site_structural": _rate(
                outcome_counts["structural_metal_sites"], total_sites
            ),
        },
        "chemical_scope": {
            "one_in_domain_dblock_center_entries": outcome_counts[
                "one_in_domain_dblock_center_entries"
            ],
            "multiple_in_domain_dblock_center_entries": outcome_counts[
                "multiple_in_domain_dblock_center_entries"
            ],
            "strict_nonpolymeric_single_metal_center_entries": outcome_counts[
                "strict_nonpolymeric_single_metal_center_entries"
            ],
            "nonpolymeric_multiple_metal_centers_or_external_relation_entries": outcome_counts[
                "nonpolymeric_multiple_metal_centers_or_external_relation_entries"
            ],
            "multiple_metal_centers_or_external_relation_entries": outcome_counts[
                "multiple_metal_centers_or_external_relation_entries"
            ],
            "haptic_entries": record_counts["entries_with_haptic_site"],
            "haptic_entries_nonambiguous": record_counts[
                "entries_with_haptic_site_nonambiguous"
            ],
            "haptic_entries_disorder_ambiguous": record_counts[
                "entries_with_haptic_site_disorder_ambiguous"
            ],
            "strict_nonpolymeric_single_metal_center_haptic_entries": outcome_counts[
                "strict_nonpolymeric_single_metal_center_haptic_entries"
            ],
            "nonpolymeric_multiple_metal_centers_or_external_relation_haptic_entries": outcome_counts[
                "nonpolymeric_multiple_metal_centers_or_external_relation_haptic_entries"
            ],
            "polymeric_haptic_entries": outcome_counts[
                "polymeric_haptic_entries"
            ],
            "ambiguous_pi_candidate_entries": record_counts[
                "entries_with_ambiguous_pi_candidate"
            ],
            "bridged_entries": record_counts["entries_with_bridge"],
            "any_in_domain_or_external_shared_donor_entries": record_counts[
                "entries_with_any_shared_donor"
            ],
            "direct_metal_bond_entries": record_counts["entries_with_direct_metal_bond"],
            "external_metal_center_relation_entries": record_counts[
                "entries_with_external_metal_center"
            ],
            "external_metal_shared_donor_entries": record_counts[
                "entries_with_external_metal_shared_donor"
            ],
            "metal_bound_hydrogen_entries": record_counts[
                "entries_with_metal_bound_hydrogen"
            ],
            "heterometal_co_center_entries": outcome_counts[
                "flag:heterometal_co_center"
            ],
            "polymeric_entries": outcome_counts["flag:polymeric"],
            "nonzero_translation_edge_entries": record_counts[
                "entries_with_nonzero_translation_edge"
            ],
            "nonzero_translation_metal_donor_incidence_entries": record_counts[
                "entries_with_nonzero_translation_incidence"
            ],
            "unresolved_translation_edge_entries": record_counts[
                "entries_with_unresolved_translation_edge"
            ],
            "disordered_entries": outcome_counts["flag:disorder"],
            "high_nuclearity_entries": outcome_counts["flag:high_nuclearity"],
            "high_cn_entries": outcome_counts["flag:high_cn"],
        },
        "canonical_accounting": {
            "entries_exact_label_free": record_counts["entries_exact_label_free"],
            "entries_source_order_noncanonical": record_counts[
                "entries_source_order_noncanonical"
            ],
            "site_status_counts": {
                key.split(":", 1)[1]: value
                for key, value in sorted(record_counts.items())
                if key.startswith("site_canonical:")
            },
        },
        "status_counts": {
            key.split(":", 1)[1]: value
            for key, value in sorted(outcome_counts.items())
            if key.startswith("status:")
        },
        "scope_counts": {
            key.split(":", 1)[1]: value
            for key, value in sorted(outcome_counts.items())
            if key.startswith("scope:")
        },
        "scope_coverage": {
            scope: {
                "entry_any_structural": _rate(
                    counts["entries_any_structural"], counts["entries"]
                ),
                "entry_all_metal_sites_structural": _rate(
                    counts["entries_all_sites_structural"], counts["entries"]
                ),
                "metal_site_structural": _rate(
                    counts["structural_metal_sites"], counts["metal_sites"]
                ),
            }
            for scope, counts in sorted(per_scope.items())
        },
        "site_record_level_counts": {
            key.split(":", 1)[1]: value
            for key, value in sorted(record_counts.items())
            if key.startswith("site_record_level:")
        },
        "site_resolution_counts": {
            key.split(":", 1)[1]: value
            for key, value in sorted(record_counts.items())
            if key.startswith("site_resolution:")
        },
        "geometry_kind_counts": {
            key.split(":", 1)[1]: value
            for key, value in sorted(record_counts.items())
            if key.startswith("geometry_kind:")
        },
        "donor_pair_angle_class_counts": {
            key.split(":", 1)[1]: value
            for key, value in sorted(record_counts.items())
            if key.startswith("angle_class:")
        },
        "issue_counts": {
            key.split(":", 1)[1]: value
            for key, value in sorted(outcome_counts.items())
            if key.startswith("issue:")
        },
        "failure_issue_counts": {
            key.split(":", 1)[1]: value
            for key, value in sorted(outcome_counts.items())
            if key.startswith("failure:")
        },
        "accounting_checks": {
            "outcome_rows_equal_frozen_tm_denominator": total_entries == expected_entries,
            "record_rows_equal_emitted_outcomes": (
                record_counts["records"] == outcome_counts["emitted_entries"]
            ),
            "metal_site_fields_closed": (
                total_sites
                == outcome_counts["structural_metal_sites"]
                + outcome_counts["audit_only_metal_sites"]
                + outcome_counts["failed_metal_sites"]
            ),
            "canonical_entry_partition_closed": (
                record_counts["records"]
                == record_counts["entries_exact_label_free"]
                + record_counts["entries_source_order_noncanonical"]
            ),
            "one_multiple_dblock_center_partition_closed": (
                outcome_counts["emitted_entries"]
                == outcome_counts["one_in_domain_dblock_center_entries"]
                + outcome_counts["multiple_in_domain_dblock_center_entries"]
            ),
            "nonpolymeric_mono_multi_and_polymeric_partition_closed": (
                outcome_counts["emitted_entries"]
                == outcome_counts[
                    "strict_nonpolymeric_single_metal_center_entries"
                ]
                + outcome_counts[
                    "nonpolymeric_multiple_metal_centers_or_external_relation_entries"
                ]
                + outcome_counts["flag:polymeric"]
            ),
            "confirmed_haptic_main_class_partition_closed": (
                record_counts["entries_with_haptic_site"]
                == outcome_counts[
                    "strict_nonpolymeric_single_metal_center_haptic_entries"
                ]
                + outcome_counts[
                    "nonpolymeric_multiple_metal_centers_or_external_relation_haptic_entries"
                ]
                + outcome_counts["polymeric_haptic_entries"]
            ),
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    _atomic_json(args.output, summary)
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
