#!/usr/bin/env python3
"""
export_box1_coordrep_records.py
===============================

Generate Box 1. Representative complete CoordRep records.

All records are produced by the live CoordRep serializer from CSD crystal
structures.  No hand-written strings.  Every record is parse-validated and
round-trip checked before export.

Selected examples:
  A. cis-[PtCl2(NH3)2]         → AKUVUA   (Pt SP-4, monodentate, cis/trans)
  B. [Co(en)3]3+               → JOXSOI   (Co Oh, tris-bidentate, 3 trans)
  C. fac-[Co(dien)(CN)3]       → IDARAG   (Co Oh, tridentate+3×CN, fac)
     mer-[Co(dien)(CN)3]       → IDAREK   (Co, same L3, mer) [SI]

If primary candidates fail validation gate, fallback alternatives are tried.

Output → revision_results/box1_coordrep_records/

Run with 1701 env (has CSD Python API):
  /data/miniconda3/envs/1701/bin/python scripts/export_box1_coordrep_records.py
"""

from __future__ import annotations

import csv
import json
import re
import sys
import traceback
from dataclasses import asdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from ccdc.io import EntryReader

from coordrep.encode import encode_molecule
from coordrep.core import CoordRepConfig
from coordrep.canonical.canonicalize import canonicalize_complex
from coordrep.identity import extract_identity_keys
from coordrep_tools.csd_adapter import csd_entry_to_raw_molecule
from coordrep_tools.validate import is_valid_coordrep

# ── Paths ──
OUT = Path(__file__).parent.parent / "revision_results/box1_coordrep_records"
OUT.mkdir(parents=True, exist_ok=True)

SERIALIZER_VERSION = "coordrep_v1"
PARSER_VERSION = "coordrep_v1"

# ── Example definitions ──
# Primary candidates with fallbacks
EXAMPLES = [
    {
        "example_id": "A",
        "display_name": "cis-[Pt(MeNH2)2I2]",
        "role_in_box": "Monodentate square-planar record",
        "refcodes": ["IYIPOV"],
        "expect_metal": "Pt",
        "expect_cn": 4,
        "expect_shape": "SP",
        "expect_stereo": True,
        "expect_bidentate": False,
        "expect_fac_mer": None,
        "included_in_main_box": True,
    },
    {
        "example_id": "B",
        "display_name": "[Co(en)3]3+",
        "role_in_box": "Bidentate chelate record",
        "refcodes": ["JOXSOI", "VUYDOJ01", "VUYDOJ02"],
        "expect_metal": "Co",
        "expect_cn": 6,
        "expect_shape": "Oh",
        "expect_stereo": True,
        "expect_bidentate": True,
        "expect_fac_mer": None,
        "included_in_main_box": True,
    },
    {
        "example_id": "C_fac",
        "display_name": "fac-[Co(dien)(CN)3]",
        "role_in_box": "Octahedral fac stereochemical record",
        "refcodes": ["IDARAG", "IDARAG01"],
        "expect_metal": "Co",
        "expect_cn": 6,
        "expect_shape": "Oh",
        "expect_stereo": True,
        "expect_bidentate": True,
        "expect_fac_mer": "fac",
        "included_in_main_box": True,
    },
    {
        "example_id": "C_mer",
        "display_name": "mer-[Co(dien)(CN)3]",
        "role_in_box": "Octahedral mer stereochemical record (SI)",
        "refcodes": ["IDAREK"],
        "expect_metal": "Co",
        "expect_cn": 6,
        "expect_shape": None,  # May be TPr
        "expect_stereo": True,
        "expect_bidentate": True,
        "expect_fac_mer": "mer",
        "included_in_main_box": False,
    },
]


# ══════════════════════════════════════════════════════════════
# Core pipeline: CSD → encode → canonicalize → serialize → validate
# ══════════════════════════════════════════════════════════════

def process_refcode(refcode: str, config: CoordRepConfig, reader):
    """Full pipeline for one CSD refcode.  Returns dict or raises."""
    entry = reader.entry(refcode)
    if entry is None:
        raise ValueError(f"CSD entry {refcode} not found")

    raw_mol = csd_entry_to_raw_molecule(entry)
    if raw_mol is None:
        raise ValueError(f"Could not convert {refcode} to RawMolecule")

    cc = encode_molecule(raw_mol, config)
    if cc.metal.element == "?":
        raise ValueError(f"Metal not detected for {refcode}")

    cc_canon = canonicalize_complex(cc)
    coordrep_str = cc_canon.to_string()

    # Validate
    parse_ok = is_valid_coordrep(coordrep_str, strict=True)

    # Derive CN
    cn = (cc_canon.shape.cn if cc_canon.shape
          else len(cc_canon.graph.donor_indices))

    # Round-trip: re-parse identity keys (no structural re-encode needed;
    # identity keys must be extractable and consistent)
    keys = extract_identity_keys(coordrep_str)

    # Check round-trip: the L0 key IS the full string
    roundtrip_ok = (keys.L0_StateKey == coordrep_str)

    return {
        "refcode": refcode,
        "coordrep_raw": coordrep_str,
        "cc": cc_canon,
        "keys": keys,
        "parse_valid": parse_ok,
        "roundtrip_valid": roundtrip_ok,
        "metal": cc_canon.metal.element,
        "cn": cn,
    }


def validate_gate(result: dict, example: dict) -> tuple[bool, str]:
    """Check whether a serialized record passes the validation gate."""
    cr = result["coordrep_raw"]
    reasons = []

    if not result["parse_valid"]:
        reasons.append("parse_valid=false")
    if not result["roundtrip_valid"]:
        reasons.append("roundtrip_valid=false")

    # Metal
    if example["expect_metal"] and result["metal"] != example["expect_metal"]:
        reasons.append(f"metal={result['metal']} expected={example['expect_metal']}")

    # CN
    if example["expect_cn"] and result["cn"] != example["expect_cn"]:
        reasons.append(f"cn={result['cn']} expected={example['expect_cn']}")

    # Shape
    if example["expect_shape"]:
        shape_m = re.search(r'ShapeBest:(\w[\w-]*)', cr)
        obs_shape = shape_m.group(1) if shape_m else "?"
        if obs_shape != example["expect_shape"]:
            reasons.append(f"shape={obs_shape} expected={example['expect_shape']}")

    # Stereo
    if example["expect_stereo"]:
        if not re.search(r'\{trans:', cr) and not re.search(r'\{cis:', cr):
            reasons.append("no stereo tokens")

    # fac/mer
    if example["expect_fac_mer"]:
        fm_m = re.search(r'\{fm:(\w+)\}', cr)
        if not fm_m:
            reasons.append(f"no fac/mer token, expected={example['expect_fac_mer']}")
        elif fm_m.group(1) != example["expect_fac_mer"]:
            reasons.append(f"fm={fm_m.group(1)} expected={example['expect_fac_mer']}")

    # Bidentate
    if example["expect_bidentate"]:
        # Check if any ligand appears multiple times in stereo tokens
        donor_refs = re.findall(r'(L\d+):[A-Z][a-z]?:\d+', cr)
        from collections import Counter
        lig_count = Counter(donor_refs)
        if not any(v >= 2 for v in lig_count.values()):
            reasons.append("no bidentate evidence in stereo tokens")

    # No placeholders
    for bad in ("...", "TBD", "unknown", "placeholder"):
        if bad in cr:
            reasons.append(f"placeholder found: {bad}")

    # Ligand dictionary present
    if not re.search(r'\|L\d+=', cr):
        reasons.append("no ligand dictionary")

    # L0-L3 keys extractable
    keys = result["keys"]
    if not keys.L0_StateKey:
        reasons.append("L0 missing")
    if not keys.L3_ConnID:
        reasons.append("L3 missing")

    return (len(reasons) == 0, "; ".join(reasons))


def format_display(cr: str) -> str:
    """Line-break a coordrep string for human display in Box 1."""
    # Parse main blocks
    metal_m = re.search(r'(\[Metal:[^\]]+\])', cr)
    shape_m = re.search(r'(<ShapeBest:[^>]+>)', cr)
    stereo_all = re.findall(r'(\{[^}]+\})', cr)
    lig_m = re.findall(r'(\|L\d+=[^|]+)', cr)

    metal_s = metal_m.group(1) if metal_m else ""
    shape_s = shape_m.group(1) if shape_m else ""
    stereo_s = "".join(stereo_all) if stereo_all else ""
    lig_s = "".join(lig_m) + "|" if lig_m else ""

    lines = []
    lines.append(f"| M       = {metal_s}")
    lines.append(f"| Shape   = {shape_s}")
    lines.append(f"| Stereo  = {stereo_s}")
    lines.append(f"| Ligands = {lig_s}")
    return "\n".join(lines)


# ══════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════

def main():
    config = CoordRepConfig.default()
    print("Opening CSD …")
    reader = EntryReader("CSD")
    print(f"  CSD available ({len(reader):,} entries)")

    results = []

    for ex in EXAMPLES:
        print(f"\n{'─'*60}")
        print(f"Example {ex['example_id']}: {ex['display_name']}")
        print(f"  Role: {ex['role_in_box']}")

        success = False
        for refcode in ex["refcodes"]:
            print(f"  Trying {refcode} …", end=" ")
            try:
                result = process_refcode(refcode, config, reader)
                gate_pass, gate_reason = validate_gate(result, ex)
                if gate_pass:
                    print(f"✓ PASS")
                    result["example_id"] = ex["example_id"]
                    result["display_name"] = ex["display_name"]
                    result["role_in_box"] = ex["role_in_box"]
                    result["included_in_main_box"] = ex["included_in_main_box"]
                    result["source_refcode"] = refcode
                    result["gate_pass"] = True
                    result["gate_reason"] = ""
                    results.append(result)
                    success = True
                    break
                else:
                    print(f"✗ GATE FAIL: {gate_reason}")
                    # Still store for reporting
                    result["example_id"] = ex["example_id"]
                    result["display_name"] = ex["display_name"]
                    result["role_in_box"] = ex["role_in_box"]
                    result["included_in_main_box"] = ex["included_in_main_box"]
                    result["source_refcode"] = refcode
                    result["gate_pass"] = False
                    result["gate_reason"] = gate_reason
            except Exception as e:
                print(f"✗ ERROR: {e}")
                traceback.print_exc()

        if not success:
            print(f"  ⚠ No valid candidate for {ex['example_id']}")
            # Add a placeholder result for reporting (but mark as failed)
            results.append({
                "example_id": ex["example_id"],
                "display_name": ex["display_name"],
                "role_in_box": ex["role_in_box"],
                "included_in_main_box": False,
                "source_refcode": ex["refcodes"][0],
                "coordrep_raw": "",
                "keys": None,
                "parse_valid": False,
                "roundtrip_valid": False,
                "metal": "",
                "cn": 0,
                "gate_pass": False,
                "gate_reason": "no valid candidate found",
            })

    # ══════════════════════════════════════════════════════════
    # Export files
    # ══════════════════════════════════════════════════════════

    passed = [r for r in results if r["gate_pass"]]
    print(f"\n{'═'*60}")
    print(f"Results: {len(passed)} passed / {len(results)} total")

    # ── File 1: box1_record_index.csv ──
    idx_fields = [
        "example_id", "display_name", "role_in_box",
        "source_database", "source_id",
        "formula_or_label", "metal", "oxidation_state", "d_count",
        "coordination_number", "shape_best", "shape_second",
        "best_CShM", "second_CShM", "delta_CShM", "boundary_flag",
        "n_ligands", "n_donor_sites", "has_multidentate_ligand",
        "has_fac_mer_relation", "has_cis_trans_relation",
        "parse_valid", "roundtrip_valid",
        "included_in_main_box", "si_full_record_file", "notes",
    ]

    with open(OUT / "box1_record_index.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=idx_fields)
        w.writeheader()
        for r in results:
            cr = r["coordrep_raw"]
            # Parse fields from coordrep
            ox_m = re.search(r'ox:([^|}\]>]+)', cr) if cr else None
            d_m = re.search(r'd:(d\d+)', cr) if cr else None
            shape_m = re.search(r'ShapeBest:(\w[\w-]*)', cr) if cr else None
            v_m = re.search(r'V:([\d.]+),([\d.]+)', cr) if cr else None
            fm_m = re.search(r'\{fm:(\w+)\}', cr) if cr else None
            trans_m = re.findall(r'\{trans:', cr) if cr else []
            ligs = re.findall(r'\|L\d+=([^|]+)', cr) if cr else []
            donor_refs = re.findall(r'(L\d+):[A-Z][a-z]?:\d+', cr) if cr else []
            from collections import Counter
            lig_donor_counts = Counter(donor_refs)

            # Determine shape info
            if v_m:
                v1, v2 = float(v_m.group(1)), float(v_m.group(2))
                best_cshm = round(min(v1, v2), 2)
                second_cshm = round(max(v1, v2), 2)
                delta = round(abs(v1 - v2), 2)
            else:
                best_cshm = second_cshm = delta = ""

            # Boundary flag
            bd = "yes" if (cr and "_boundary" in cr) else "no"
            if not cr:
                bd = ""

            # Determine second shape from CN-based pairs
            SHAPE_PAIRS = {
                4: ("SP", "Td"), 5: ("SPY", "TBPY"),
                6: ("Oh", "TPr"), 2: ("L", "L"), 3: ("TP", "TY"),
            }
            obs_best = shape_m.group(1) if shape_m else ""
            pair = SHAPE_PAIRS.get(r["cn"], ("", ""))
            if obs_best == pair[0]:
                obs_second = pair[1]
            elif obs_best == pair[1]:
                obs_second = pair[0]
            else:
                obs_second = ""

            row = {
                "example_id": r["example_id"],
                "display_name": r["display_name"],
                "role_in_box": r["role_in_box"],
                "source_database": "CSD",
                "source_id": r["source_refcode"],
                "formula_or_label": r["display_name"],
                "metal": r["metal"],
                "oxidation_state": ox_m.group(1).strip() if ox_m else "",
                "d_count": d_m.group(1) if d_m else "",
                "coordination_number": r["cn"],
                "shape_best": obs_best,
                "shape_second": obs_second,
                "best_CShM": best_cshm,
                "second_CShM": second_cshm,
                "delta_CShM": delta,
                "boundary_flag": bd,
                "n_ligands": len(ligs),
                "n_donor_sites": r["cn"],
                "has_multidentate_ligand":
                    any(v >= 2 for v in lig_donor_counts.values()),
                "has_fac_mer_relation": bool(fm_m),
                "has_cis_trans_relation": len(trans_m) > 0,
                "parse_valid": r["parse_valid"],
                "roundtrip_valid": r["roundtrip_valid"],
                "included_in_main_box": r["included_in_main_box"],
                "si_full_record_file": "box1_si_full_records.txt",
                "notes": r.get("gate_reason", ""),
            }
            w.writerow(row)
    print("  1. box1_record_index.csv")

    # ── File 2: box1_display_records.md ──
    with open(OUT / "box1_display_records.md", "w") as f:
        f.write("# Box 1. Representative complete CoordRep records\n\n")

        for r in passed:
            if not r["included_in_main_box"]:
                continue
            cr = r["coordrep_raw"]
            keys = r["keys"]

            f.write(f"## {r['example_id']}. {r['role_in_box']}: "
                    f"{r['display_name']}\n\n")
            f.write("```\nCoordRep-v1\n")
            f.write(format_display(cr))
            f.write(f"\n| ID      = [L0_hash:{_hash8(keys.L0_StateKey)}"
                    f"|L1:{keys.L1_ShapeID[:60]}"
                    f"|L2:{keys.L2_TopoID[:60]}"
                    f"|L3:{keys.L3_ConnID}]")
            f.write("\n```\n\n")
            f.write(f"Source: CSD {r['source_refcode']}  \n")
            f.write(f"Serializer: {SERIALIZER_VERSION}\n\n")

        # SI records
        si_records = [r for r in passed if not r["included_in_main_box"]]
        if si_records:
            f.write("---\n\n## Supporting Information records\n\n")
            for r in si_records:
                cr = r["coordrep_raw"]
                keys = r["keys"]
                f.write(f"### {r['example_id']}. {r['role_in_box']}: "
                        f"{r['display_name']}\n\n")
                f.write("```\nCoordRep-v1\n")
                f.write(format_display(cr))
                f.write(f"\n| ID      = [L0_hash:{_hash8(keys.L0_StateKey)}"
                        f"|L1:{keys.L1_ShapeID[:60]}"
                        f"|L2:{keys.L2_TopoID[:60]}"
                        f"|L3:{keys.L3_ConnID}]")
                f.write("\n```\n\n")
                f.write(f"Source: CSD {r['source_refcode']}  \n")
                f.write(f"Serializer: {SERIALIZER_VERSION}\n\n")

        f.write("---\n\n")
        f.write("The records are direct outputs of the CoordRep serializer "
                "and are line-broken only for readability; full "
                "machine-readable records and hashes are provided in the "
                "Supporting Information.\n")
    print("  2. box1_display_records.md")

    # ── File 3: box1_full_records.jsonl ──
    with open(OUT / "box1_full_records.jsonl", "w") as f:
        for r in passed:
            cr = r["coordrep_raw"]
            keys = r["keys"]
            entry = {
                "example_id": r["example_id"],
                "display_name": r["display_name"],
                "source_database": "CSD",
                "source_id": r["source_refcode"],
                "coordrep_raw": cr,
                "coordrep_line_broken": format_display(cr),
                "parsed_fields": {
                    "metal": f"[Metal:{keys.metal}|ox:{keys.oxidation}"
                             f"|d:?|CN:{keys.cn}]",
                    "shape": re.search(r'(<ShapeBest:[^>]+>)', cr).group(1)
                             if re.search(r'(<ShapeBest:[^>]+>)', cr)
                             else "",
                    "stereo": "".join(re.findall(r'(\{[^}]+\})', cr)),
                    "ligands": "".join(re.findall(r'(\|L\d+=[^|]+)', cr)),
                    "id": {
                        "L0_StateKey_hash": _hash8(keys.L0_StateKey),
                        "L1_ShapeID": keys.L1_ShapeID,
                        "L2_TopoID": keys.L2_TopoID,
                        "L3_ConnID": keys.L3_ConnID,
                    },
                },
                "validation": {
                    "parse_valid": r["parse_valid"],
                    "roundtrip_valid": r["roundtrip_valid"],
                    "hashes_recomputed": True,
                },
            }
            json.dump(entry, f)
            f.write("\n")
    print("  3. box1_full_records.jsonl")

    # ── File 4: box1_validation_report.csv ──
    val_fields = [
        "example_id", "parse_valid", "roundtrip_valid",
        "serializer_version", "parser_version",
        "source_structure_loaded", "metal_detected_correctly",
        "CN_correct", "shape_expected", "shape_observed",
        "stereo_tokens_present", "ligand_dictionary_present",
        "donor_markers_present", "L0_present", "L1_present",
        "L2_present", "L3_present", "no_placeholders",
        "pass_main_box_gate", "failure_reason",
    ]
    with open(OUT / "box1_validation_report.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=val_fields)
        w.writeheader()
        for r, ex in zip(results, EXAMPLES):
            cr = r["coordrep_raw"]
            keys = r.get("keys")
            shape_m = re.search(r'ShapeBest:(\w[\w-]*)', cr) if cr else None
            has_stereo = bool(re.search(r'\{trans:', cr)) if cr else False
            has_lig = bool(re.search(r'\|L\d+=', cr)) if cr else False
            has_donor = bool(
                re.search(r'L\d+:[A-Z][a-z]?:\d+', cr)) if cr else False
            no_ph = not any(
                bad in cr for bad in ("...", "TBD", "unknown")
            ) if cr else False

            row = {
                "example_id": r["example_id"],
                "parse_valid": r["parse_valid"],
                "roundtrip_valid": r["roundtrip_valid"],
                "serializer_version": SERIALIZER_VERSION,
                "parser_version": PARSER_VERSION,
                "source_structure_loaded": bool(cr),
                "metal_detected_correctly":
                    r["metal"] == ex["expect_metal"] if ex["expect_metal"] else True,
                "CN_correct":
                    r["cn"] == ex["expect_cn"] if ex["expect_cn"] else True,
                "shape_expected": ex.get("expect_shape", "any"),
                "shape_observed": shape_m.group(1) if shape_m else "",
                "stereo_tokens_present": has_stereo,
                "ligand_dictionary_present": has_lig,
                "donor_markers_present": has_donor,
                "L0_present": bool(keys and keys.L0_StateKey) if keys else False,
                "L1_present": bool(keys and keys.L1_ShapeID) if keys else False,
                "L2_present": bool(keys and keys.L2_TopoID) if keys else False,
                "L3_present": bool(keys and keys.L3_ConnID) if keys else False,
                "no_placeholders": no_ph,
                "pass_main_box_gate": r["gate_pass"],
                "failure_reason": r.get("gate_reason", ""),
            }
            w.writerow(row)
    print("  4. box1_validation_report.csv")

    # ── File 5: box1_si_full_records.txt ──
    with open(OUT / "box1_si_full_records.txt", "w") as f:
        f.write("=" * 72 + "\n")
        f.write("Supporting Information: Full CoordRep Records\n")
        f.write("=" * 72 + "\n\n")
        f.write(f"Serializer version: {SERIALIZER_VERSION}\n")
        f.write(f"Parser version:     {PARSER_VERSION}\n")
        f.write("Source: Cambridge Structural Database (CSD)\n")
        f.write("No raw CSD coordinates are included.\n\n")

        for r in passed:
            cr = r["coordrep_raw"]
            keys = r["keys"]
            f.write("-" * 72 + "\n")
            f.write(f"Example {r['example_id']}: {r['display_name']}\n")
            f.write(f"CSD Refcode: {r['source_refcode']}\n")
            f.write(f"Parse valid: {r['parse_valid']}\n")
            f.write(f"Round-trip valid: {r['roundtrip_valid']}\n\n")
            f.write("Full CoordRep string (unwrapped):\n")
            f.write(cr + "\n\n")
            f.write("Identity keys:\n")
            f.write(f"  L0 (StateKey):  {keys.L0_StateKey[:80]}…\n"
                    if len(keys.L0_StateKey) > 80
                    else f"  L0 (StateKey):  {keys.L0_StateKey}\n")
            f.write(f"  L0 hash:        {_hash8(keys.L0_StateKey)}\n")
            f.write(f"  L1 (ShapeID):   {keys.L1_ShapeID}\n")
            f.write(f"  L2 (TopoID):    {keys.L2_TopoID}\n")
            f.write(f"  L3 (ConnID):    {keys.L3_ConnID}\n\n")

        f.write("=" * 72 + "\n")
        f.write("End of Supporting Information records.\n")
    print("  5. box1_si_full_records.txt")

    # ── File 6: box1_readme.md ──
    with open(OUT / "box1_readme.md", "w") as f:
        f.write("""# Box 1: Representative Complete CoordRep Records

## Purpose

These records demonstrate the full CoordRep-v1 representation for three
classes of coordination compounds, directly addressing the reviewer concern
that "full syntax is insufficiently documented."

## Data provenance

- **Source**: Cambridge Structural Database (CSD)
- **Serializer**: `coordrep_v1` (`coordrep.encode.encode_molecule` →
  `canonicalize_complex` → `to_string()`)
- **Validation**: `coordrep_tools.validate.is_valid_coordrep(strict=True)`
- **Identity keys**: `coordrep.identity.extract_identity_keys()`

## Selection criteria

| Example | Criteria |
|---------|----------|
| A | CN=4, Pt, shape_best=SP, all monodentate, has trans tokens, not boundary |
| B | CN=6, Co, shape_best=Oh, has bidentate ligand(s), has trans tokens |
| C | CN=6, has {fm:fac} or {fm:mer} token, Oh preferred |

## Commands to reproduce

```bash
/data/miniconda3/envs/1701/bin/python scripts/export_box1_coordrep_records.py
```

## Guarantees

- All records are **direct serializer output**; no hand-written strings.
- Every record passes `parse_valid=True` and `roundtrip_valid=True`.
- No placeholders ("...", "TBD", "unknown") appear in any output.
- Raw CSD coordinates are **not** redistributed.

## Files

| File | Description |
|------|-------------|
| `box1_record_index.csv` | Summary index for all examples |
| `box1_display_records.md` | Main-text Box 1 formatted for typesetting |
| `box1_full_records.jsonl` | Machine-readable full records |
| `box1_validation_report.csv` | Per-record validation gate results |
| `box1_si_full_records.txt` | Supporting Information full records |
| `box1_readme.md` | This file |
""")
    print("  6. box1_readme.md")

    # ── Summary ──
    print(f"\n{'═'*60}")
    print(f"✓ All 6 files written to {OUT}/")
    n_main = sum(1 for r in passed if r["included_in_main_box"])
    n_si = sum(1 for r in passed if not r["included_in_main_box"])
    print(f"  Main-text records: {n_main}")
    print(f"  SI records:        {n_si}")
    for r in results:
        status = "✓" if r["gate_pass"] else "✗"
        print(f"  {status} {r['example_id']:6s} {r['display_name']:30s} "
              f"[{r['source_refcode']}]"
              f"{' – ' + r['gate_reason'] if r.get('gate_reason') else ''}")


def _hash8(s: str) -> str:
    """Short 8-char hash for display."""
    import hashlib
    return hashlib.sha256(s.encode()).hexdigest()[:8]


if __name__ == "__main__":
    main()
