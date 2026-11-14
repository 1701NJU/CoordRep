#!/usr/bin/env python
"""Relation-aware CoordStatePairs OOF benchmark.

This is a deliberately small wrapper around ``run_coordstatepairs_oof.py``.
It leaves the frozen 48,057-record cohort, attributed-2D-graph folds, model
architecture, optimizer, and evaluation code unchanged.  The only addition is
an invariant CoordRep donor-relation block.  The paired
``hybrid_mask_relation`` control keeps the complete input width but replaces
all relation columns by their train-partition means after standardization.

The relation fingerprint intentionally discards ligand labels and donor-site
indices.  Each serialized cis/trans relation contributes one count to

    relation_rel=<cis/trans>|donor_pair=<sorted elements>|same_ligand=<0/1>

so relabelling ligands, relabelling donor sites, or reordering relation tokens
does not change the features.  Ligand-payload identities are not included in
this first, conservative experiment.
"""

from __future__ import annotations

from collections import Counter
from typing import Dict, Iterable, Tuple

import run_coordstatepairs_oof as benchmark
from run_downstream_graph_benchmark import record_to_field_dict as _base_fields


RELATION_PREFIX = "relation_"


def _relation_endpoint(token: str) -> Tuple[str, str]:
    """Return (ligand label, donor element) from ``L1:N:2``."""
    fields = str(token).strip().split(":")
    if len(fields) < 3 or not fields[0] or not fields[1]:
        raise ValueError(f"Malformed relation endpoint: {token!r}")
    return fields[0], fields[1]


def _iter_relations(record: dict) -> Iterable[Tuple[str, str]]:
    for relation_type, field in (("cis", "cis_pairs"), ("trans", "trans_pairs")):
        for token in record.get(field) or []:
            yield relation_type, str(token)


def relation_features(record: dict) -> Dict[str, float]:
    """Build a ligand-permutation- and traversal-invariant relation block."""
    counts: Counter[str] = Counter()
    parse_failures = 0
    parsed_relations = 0
    for relation_type, token in _iter_relations(record):
        try:
            left, right = token.split("--", 1)
            left_ligand, left_element = _relation_endpoint(left)
            right_ligand, right_element = _relation_endpoint(right)
        except (TypeError, ValueError):
            parse_failures += 1
            continue
        donor_pair = "~".join(sorted((left_element, right_element)))
        same_ligand = int(left_ligand == right_ligand)
        key = (
            f"{RELATION_PREFIX}rel={relation_type}"
            f"|donor_pair={donor_pair}|same_ligand={same_ligand}"
        )
        counts[key] += 1
        parsed_relations += 1

    # These audit/missingness columns are namespaced so the masked control
    # removes them together with the chemically resolved relation counts.
    source_missing = not (
        "cis_pairs" in record or "trans_pairs" in record
    )
    features: Dict[str, float] = {key: float(value) for key, value in counts.items()}
    features[f"{RELATION_PREFIX}parsed_count"] = float(parsed_relations)
    features[f"{RELATION_PREFIX}parse_failure_count"] = float(parse_failures)
    features[f"{RELATION_PREFIX}source_missing"] = float(source_missing)
    return features


def record_to_relation_aware_field_dict(record: dict) -> Dict[str, float]:
    features = _base_fields(record)
    features.update(relation_features(record))
    return features


def _invariance_self_test() -> None:
    first = {
        "cis_pairs": ["L1:N:1--L2:O:2", "L1:N:2--L1:N:1"],
        "trans_pairs": ["L1:N:1--L3:Cl:1"],
    }
    relabelled_and_reordered = {
        "cis_pairs": ["L7:N:9--L7:N:4", "L4:O:8--L7:N:2"],
        "trans_pairs": ["L9:Cl:6--L7:N:5"],
    }
    if relation_features(first) != relation_features(relabelled_and_reordered):
        raise RuntimeError("Relation fingerprint invariance self-test failed")


def main() -> None:
    _invariance_self_test()
    # load_raw_cohort resolves this name in the imported benchmark module.
    benchmark.record_to_field_dict = record_to_relation_aware_field_dict
    benchmark.main()


if __name__ == "__main__":
    main()
