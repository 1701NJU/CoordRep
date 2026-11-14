#!/usr/bin/env python3
"""
scope_table_gen.py
==================
Task 5: Generate scope table / exclusion boundary summary.

Produces markdown and CSV tables for the main text and SI.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

SUPPORTED = [
    ("Mononuclear transition-metal complexes", "3d/4d/5d metals, single metal center"),
    ("Atom-resolved eta-1 donors", "Each donor atom tracked by element and rank"),
    ("Monodentate ligands", "Single-donor ligands (e.g., Cl, NH3, H2O)"),
    ("Bidentate eta-1 chelating ligands", "e.g., en, bpy, acac, oxalate"),
    ("Multidentate eta-1 chelating ligands", "dent >= 3 (e.g., EDTA, tpy, porphyrin)"),
    ("Coordination number 2-14", "Full CN range from linear to high-coordinate"),
    ("cis/trans, fac/mer, donor-relation tokens", "Stereochemical constraint grammar"),
    ("CoordRep-ID L0-L3", "Hierarchical identity: StateKey, ShapeID, TopoID, ConnID"),
]

OUT_OF_SCOPE = [
    ("Multinuclear complexes", "Multiple metal centers, bridging ligands"),
    ("MOFs / coordination polymers", "Periodic / extended lattice structures"),
    ("Polyoxometalates", "Large cluster metal-oxide frameworks"),
    ("eta-n organometallics", "Metallocenes, allyl, arene complexes (eta > 1)"),
    ("Positional disorder / partial occupancy", "Crystallographic disorder"),
    ("Ambiguous donor assignment", "Weak interactions below distance threshold"),
]

FUTURE = [
    ("CoordRep-Multi", "Multi-metal graph + bridging ligand relations"),
    ("CoordRep-Periodic", "Metal node / linker / periodic topology"),
    ("CoordRep-Haptic", "eta-n fragment token + centroid / slippage geometry"),
    ("CoordRep-Occ", "Occupancy-aware disorder representation"),
]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="revision_results/scope")
    args = parser.parse_args()

    os.makedirs(args.out, exist_ok=True)

    # ── Markdown table ────────────────────────────────────
    p = os.path.join(args.out, "coordrep_v1_scope_table.md")
    with open(p, "w") as f:
        f.write("# CoordRep v1 Scope Table\n\n")
        f.write("## Supported in CoordRep v1\n\n")
        f.write("| Feature | Details |\n|---------|--------|\n")
        for feat, detail in SUPPORTED:
            f.write(f"| {feat} | {detail} |\n")

        f.write("\n## Explicitly Out of v1 Scope\n\n")
        f.write("| Exclusion | Reason |\n|-----------|--------|\n")
        for feat, detail in OUT_OF_SCOPE:
            f.write(f"| {feat} | {detail} |\n")

        f.write("\n## Future Extension Route\n\n")
        f.write("| Extension | Description |\n|-----------|------------|\n")
        for feat, detail in FUTURE:
            f.write(f"| {feat} | {detail} |\n")

        f.write("\n---\n")
        f.write("*CoordRep v1 focuses on mononuclear eta-1 coordination snapshots, ")
        f.write("including monodentate and multidentate chelating ligands.*\n")
    print(f"  {p}")

    # ── CSV table ─────────────────────────────────────────
    p = os.path.join(args.out, "coordrep_v1_scope_table.csv")
    with open(p, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["category", "feature", "details"])
        for feat, detail in SUPPORTED:
            w.writerow(["supported_v1", feat, detail])
        for feat, detail in OUT_OF_SCOPE:
            w.writerow(["out_of_scope_v1", feat, detail])
        for feat, detail in FUTURE:
            w.writerow(["future_extension", feat, detail])
    print(f"  {p}")

    # ── scope_summary.json ────────────────────────────────
    summary = {
        "n_supported_features": len(SUPPORTED),
        "n_excluded_features": len(OUT_OF_SCOPE),
        "n_future_extensions": len(FUTURE),
        "supported": [f for f, _ in SUPPORTED],
        "excluded": [f for f, _ in OUT_OF_SCOPE],
        "future": [f for f, _ in FUTURE],
    }
    p = os.path.join(args.out, "scope_summary.json")
    with open(p, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"  {p}")

    print(f"\n{'='*60}")
    print("SCOPE TABLE GENERATED")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
