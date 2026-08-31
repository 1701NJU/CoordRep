#!/usr/bin/env python3
"""
run_full_csd_scan.py
====================
Full CSD scan → retained entries → PathFinder analysis.

Run with 1701 env (has CSD Python API):
  nohup /data/miniconda3/envs/1701/bin/python scripts/run_full_csd_scan.py > logs/full_csd_scan.log 2>&1 &

Estimated runtime: 3-5 hours on full 1.4M CSD entries.
"""

from __future__ import annotations

import json
import os
import sys
import time
import traceback
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from ccdc.io import EntryReader

from coordrep.encode import encode_molecule
from coordrep.core import CoordRepConfig
from coordrep.canonical.canonicalize import canonicalize_complex
from coordrep.identity import extract_identity_keys
from coordrep_tools.csd_adapter import (
    csd_entry_to_raw_molecule, refcode_family,
    _count_metals, _get_donor_neighbors, _has_hapticity,
)
from coordrep_tools.validate import is_valid_coordrep

BASE = Path(__file__).parent.parent
OUT_DIR = BASE / "revision_results/csd_pathfinder_full"
RETAINED_OUT = OUT_DIR / "full_csd_retained_entries.jsonl"

os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(BASE / "logs", exist_ok=True)


def main():
    config = CoordRepConfig.default()

    print("Opening CSD …")
    reader = EntryReader("CSD")
    total = len(reader)
    print(f"  CSD has {total:,} entries")

    retained = []
    waterfall = {
        "total_scanned": 0,
        "has_3d": 0,
        "transition_metal": 0,
        "mononuclear": 0,
        "eta1": 0,
        "cn_ok": 0,
        "no_disorder": 0,
        "no_polymer": 0,
        "donor_ok": 0,
        "valid_smiles": 0,
        "valid_coordrep": 0,
    }
    rejection_reasons = {}
    n_boundary = 0

    t0 = time.time()
    batch_log = 50000

    for idx in range(total):
        waterfall["total_scanned"] += 1

        if idx > 0 and idx % batch_log == 0:
            elapsed = time.time() - t0
            rate = idx / elapsed
            eta = (total - idx) / rate if rate > 0 else 0
            print(f"  [{idx:,}/{total:,}] retained={len(retained):,} "
                  f"({elapsed:.0f}s, {rate:.0f}/s, ETA {eta:.0f}s)")
            sys.stdout.flush()

        try:
            entry = reader[idx]
        except Exception:
            continue

        try:
            # Waterfall stages
            if not entry.has_3d_structure:
                rejection_reasons["no_3d"] = rejection_reasons.get("no_3d", 0) + 1
                continue
            waterfall["has_3d"] += 1

            mol_csd = entry.molecule
            if mol_csd is None:
                continue

            metals, n_metals = _count_metals(mol_csd)
            if n_metals == 0:
                rejection_reasons["no_TM"] = rejection_reasons.get("no_TM", 0) + 1
                continue
            waterfall["transition_metal"] += 1

            if n_metals > 1:
                rejection_reasons["multinuclear"] = rejection_reasons.get("multinuclear", 0) + 1
                continue
            waterfall["mononuclear"] += 1

            metal_atom = metals[0]
            if _has_hapticity(mol_csd, metal_atom):
                rejection_reasons["hapticity"] = rejection_reasons.get("hapticity", 0) + 1
                continue
            waterfall["eta1"] += 1

            donors = _get_donor_neighbors(mol_csd, metal_atom)
            cn = len(donors)
            if cn < 2 or cn > 6:
                rejection_reasons[f"cn_{cn}"] = rejection_reasons.get(f"cn_{cn}", 0) + 1
                continue
            waterfall["cn_ok"] += 1

            if entry.has_disorder:
                rejection_reasons["disorder"] = rejection_reasons.get("disorder", 0) + 1
                continue
            waterfall["no_disorder"] += 1

            if entry.is_polymeric:
                rejection_reasons["polymeric"] = rejection_reasons.get("polymeric", 0) + 1
                continue
            waterfall["no_polymer"] += 1

            raw_mol = csd_entry_to_raw_molecule(entry)
            if raw_mol is None:
                rejection_reasons["raw_mol_fail"] = rejection_reasons.get("raw_mol_fail", 0) + 1
                continue
            waterfall["donor_ok"] += 1

            cc = encode_molecule(raw_mol, config)
            if cc.metal.element == "?":
                continue
            if not all(lig.smiles and len(lig.smiles) > 0 for lig in cc.ligands) or len(cc.ligands) == 0:
                rejection_reasons["invalid_smiles"] = rejection_reasons.get("invalid_smiles", 0) + 1
                continue
            waterfall["valid_smiles"] += 1

            cc_canon = canonicalize_complex(cc)
            coordrep_str = cc_canon.to_string()

            if not is_valid_coordrep(coordrep_str, strict=True):
                rejection_reasons["invalid_coordrep"] = rejection_reasons.get("invalid_coordrep", 0) + 1
                continue
            waterfall["valid_coordrep"] += 1

            try:
                keys = extract_identity_keys(coordrep_str)
            except Exception:
                keys = None

            is_boundary = False
            if keys and keys.binned_shape:
                is_boundary = keys.binned_shape.is_boundary

            if is_boundary:
                n_boundary += 1

            retained.append({
                "refcode": entry.identifier,
                "family": refcode_family(entry.identifier),
                "metal": metal_atom.atomic_symbol,
                "cn": cn,
                "coordrep": coordrep_str,
                "L0": keys.L0_StateKey if keys else "",
                "L1": keys.L1_ShapeID if keys else "",
                "L2": keys.L2_TopoID if keys else "",
                "L3": keys.L3_ConnID if keys else "",
                "is_boundary": is_boundary,
                "best_shape": keys.best_shape if keys else "",
            })

        except Exception:
            rejection_reasons["error"] = rejection_reasons.get("error", 0) + 1

    elapsed = time.time() - t0
    print(f"\nDone: scanned {total:,}, retained {len(retained):,} ({elapsed:.0f}s)")

    # Write retained entries
    with open(RETAINED_OUT, "w") as f:
        for r in retained:
            f.write(json.dumps(r) + "\n")
    print(f"  Written: {RETAINED_OUT}")

    # Write summary
    summary = {
        "csd_release": "2024.3",
        "scan_date": time.strftime("%Y-%m-%d %H:%M:%S"),
        "total_scanned": total,
        "retained": len(retained),
        "retention_rate": round(len(retained) / total, 4),
        "n_boundary": n_boundary,
        "boundary_fraction": round(n_boundary / max(len(retained), 1), 4),
        "waterfall": waterfall,
        "top_rejections": dict(sorted(rejection_reasons.items(), key=lambda x: -x[1])[:15]),
        "elapsed_seconds": round(elapsed, 1),
    }
    with open(OUT_DIR / "full_csd_scan_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"  Summary: {json.dumps(summary, indent=2)}")

    # Now run PathFinder analysis on the new data
    print("\n" + "=" * 70)
    print("Running PathFinder analysis on full-CSD data…")
    print("=" * 70)
    os.system(f"{sys.executable} {BASE}/scripts/run_csd_pathfinder.py --use-existing")


if __name__ == "__main__":
    main()
