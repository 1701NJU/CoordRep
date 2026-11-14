#!/usr/bin/env python3
"""
csd_external_audit.py
=====================
Task 1: CSD coverage / filtering waterfall.

Scans CSD entries, applies the CoordRep filtering pipeline, attempts
CoordRep generation, and reports aggregate statistics.

Usage:
    python scripts/csd_external_audit.py [--max_entries N] [--out DIR]

License: no raw CSD coordinates or CIF data are exported.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
import traceback
from pathlib import Path

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from ccdc.io import EntryReader

from coordrep.encode import encode_molecule
from coordrep.core import CoordRepConfig
from coordrep.canonical.canonicalize import canonicalize_complex
from coordrep.identity import extract_identity_keys
from coordrep_tools.csd_adapter import (
    filter_csd_entry, csd_entry_to_raw_molecule, refcode_family,
)
from coordrep_tools.csd_filtering_stats import WaterfallTracker
from coordrep_tools.validate import is_valid_coordrep


def _process_one_entry(entry, tracker, config, retained):
    """Process one CSD entry through the full waterfall. Appends to retained if passes."""
    from coordrep_tools.csd_adapter import (
        _count_metals, _get_donor_neighbors, csd_entry_to_raw_molecule, refcode_family,
    )

    # Stage: has_3d
    if not entry.has_3d_structure:
        tracker.rejection_reasons["no_3d_structure"] += 1
        return
    tracker.record_waterfall_stage("has_3d")

    mol_csd = entry.molecule
    if mol_csd is None:
        tracker.rejection_reasons["no_molecule"] += 1
        return

    # Stage: transition_metal
    metals, n_metals = _count_metals(mol_csd)
    if n_metals == 0:
        tracker.rejection_reasons["no_transition_metal"] += 1
        return
    tracker.record_waterfall_stage("transition_metal")

    metal_elem = metals[0].atomic_symbol

    # Stage: mononuclear
    if n_metals > 1:
        tracker.rejection_reasons["multinuclear"] += 1
        return
    tracker.record_waterfall_stage("mononuclear")

    # Stage: hapticity — check bond types for pi/delocalized
    metal_atom = metals[0]
    has_haptic = False
    for b in metal_atom.bonds:
        bt = str(b.bond_type) if b.bond_type else ""
        if "pi" in bt.lower() or "delocalized" in bt.lower():
            has_haptic = True
            break
    if has_haptic:
        tracker.rejection_reasons["hapticity"] += 1
        return
    tracker.record_waterfall_stage("eta1")

    # Stage: CN range
    donors = _get_donor_neighbors(mol_csd, metal_atom)
    cn = len(donors)
    if cn < 2 or cn > 6:
        tracker.rejection_reasons[f"cn_out_of_range_{cn}"] += 1
        return
    tracker.record_waterfall_stage("cn_ok")

    # Stage: disorder
    if entry.has_disorder:
        tracker.rejection_reasons["disorder"] += 1
        return
    tracker.record_waterfall_stage("no_disorder")

    # Stage: polymeric
    if entry.is_polymeric:
        tracker.rejection_reasons["polymeric"] += 1
        return
    tracker.record_waterfall_stage("no_polymer")

    # Convert to RawMolecule
    raw_mol = csd_entry_to_raw_molecule(entry)
    if raw_mol is None:
        tracker.rejection_reasons["raw_mol_conversion_failed"] += 1
        return
    tracker.record_waterfall_stage("donor_ok")

    # Encode to CoordComplex
    cc = encode_molecule(raw_mol, config)
    if cc.metal.element == "?":
        tracker.rejection_reasons["encode_no_metal"] += 1
        return

    # Check ligand SMILES validity
    has_valid_smiles = (all(lig.smiles and len(lig.smiles) > 0
                           for lig in cc.ligands) and len(cc.ligands) > 0)
    if not has_valid_smiles:
        tracker.rejection_reasons["invalid_smiles"] += 1
        return
    tracker.record_waterfall_stage("valid_smiles")

    # Canonicalize + serialize
    cc_canon = canonicalize_complex(cc)
    coordrep_str = cc_canon.to_string()

    if not is_valid_coordrep(coordrep_str, strict=True):
        tracker.rejection_reasons["invalid_coordrep"] += 1
        return
    tracker.record_waterfall_stage("valid_coordrep")

    # Identity keys
    try:
        keys = extract_identity_keys(coordrep_str)
    except Exception:
        keys = None

    is_boundary = False
    if keys and keys.binned_shape:
        is_boundary = keys.binned_shape.is_boundary

    if is_boundary:
        tracker.record_boundary()

    tracker.record_filter_result(True, "", metal=metal_elem, cn=cn)

    retained.append({
        "refcode": entry.identifier,
        "family": refcode_family(entry.identifier),
        "metal": metal_elem,
        "cn": cn,
        "coordrep": coordrep_str,
        "L0": keys.L0_StateKey if keys else "",
        "L1": keys.L1_ShapeID if keys else "",
        "L2": keys.L2_TopoID if keys else "",
        "L3": keys.L3_ConnID if keys else "",
        "is_boundary": is_boundary,
        "best_shape": keys.best_shape if keys else "",
    })


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--max_entries", type=int, default=0,
                        help="Max CSD entries to scan (0 = all)")
    parser.add_argument("--out", default="revision_results/csd_external")
    parser.add_argument("--batch_log", type=int, default=50000,
                        help="Log progress every N entries")
    args = parser.parse_args()

    os.makedirs(args.out, exist_ok=True)
    config = CoordRepConfig.default()

    print("Opening CSD …")
    reader = EntryReader("CSD")
    total = len(reader)
    limit = args.max_entries if args.max_entries > 0 else total
    print(f"  CSD has {total:,} entries; scanning up to {limit:,}")

    tracker = WaterfallTracker()

    # Retained entries: store only refcode + CoordRep string + identity keys
    retained = []  # list of dicts (refcode, metal, cn, coordrep_str, L0..L3, is_boundary)

    t0 = time.time()

    from coordrep_tools.csd_adapter import _count_metals, _get_donor_neighbors

    for idx in range(min(limit, total)):
        tracker.record_scan()

        if idx > 0 and idx % args.batch_log == 0:
            elapsed = time.time() - t0
            rate = idx / elapsed
            eta = (limit - idx) / rate if rate > 0 else 0
            print(f"  [{idx:,}/{limit:,}] retained={len(retained):,} "
                  f"({elapsed:.0f}s, {rate:.0f}/s, ETA {eta:.0f}s)")

        # Wrap entire per-entry processing — CSD API can segfault
        # on corrupt entries; we use subprocess-level protection in
        # the outer driver if needed, but at minimum catch exceptions.
        try:
            entry = reader[idx]
        except Exception:
            continue

        try:
            result = _process_one_entry(
                entry, tracker, config, retained)
        except Exception:
            tracker.rejection_reasons["unexpected_error"] += 1

    elapsed = time.time() - t0
    print(f"\nDone: scanned {tracker.n_scanned:,}, retained {len(retained):,} "
          f"({elapsed:.0f}s)")

    # ── Write outputs ─────────────────────────────────────────
    tracker.write_all(args.out)

    # Save retained refcode-level summary (no raw coords)
    retained_path = os.path.join(args.out, "csd_retained_entries.jsonl")
    with open(retained_path, "w") as f:
        for r in retained:
            f.write(json.dumps(r) + "\n")
    print(f"  Retained entries: {retained_path}")

    # Print summary
    s = tracker.summary_dict()
    print(f"\n{'='*60}")
    print("CSD EXTERNAL AUDIT SUMMARY")
    print(f"{'='*60}")
    for k, v in s.items():
        if k != "top_rejection_reasons":
            print(f"  {k}: {v}")
    print("  Top rejection reasons:")
    for k, v in s["top_rejection_reasons"].items():
        print(f"    {k}: {v}")
    print()


if __name__ == "__main__":
    main()
