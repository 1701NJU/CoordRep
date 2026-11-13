"""Compare independently supplied complete donor references; never feeds the encoder.

Reference vectors and generated outputs may be licensed. Keep them local.
The positional tolerance is fixed at the reported 0.005 angstrom.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.optimize import linear_sum_assignment

POSITION_TOLERANCE_A = 0.005


def compare_images(expected, observed):
    """Maximum one-to-one compatible matching; count both omissions and extras."""
    n, m = len(expected), len(observed)
    compatible = np.zeros((n, m), dtype=bool)
    distances = np.full((n, m), np.inf)
    for i, e in enumerate(expected):
        for j, o in enumerate(observed):
            if e["element"] != o["element"]:
                continue
            if "relative_translation" in e and list(e["relative_translation"]) != list(o.get("relative_translation", [])):
                continue
            d = float(np.linalg.norm(np.asarray(e["vector_A"], dtype=float) - np.asarray(o["vector_A"], dtype=float)))
            distances[i, j] = d
            compatible[i, j] = d <= POSITION_TOLERANCE_A
    if n and m:
        # Large incompatibility cost maximizes compatible matches before distance.
        cost = np.where(compatible, distances, 1e6)
        rr, cc = linear_sum_assignment(cost)
        matches = [(int(i), int(j), float(distances[i, j])) for i, j in zip(rr, cc) if compatible[i, j]]
    else:
        matches = []
    ei, oi = {i for i, _, _ in matches}, {j for _, j, _ in matches}
    missing = [i for i in range(n) if i not in ei]
    extra = [j for j in range(m) if j not in oi]
    return {"complete_set_equal": not missing and not extra,
            "expected_count": n, "observed_count": m, "matched": len(matches),
            "missing_reference_indices0": missing, "additional_observed_indices0": extra,
            "matched_positions": matches, "tolerance_A": POSITION_TOLERANCE_A}


def compare_reference(output, reference):
    route = output["route"]
    if route == "xyz":
        candidates = [x for x in output["sites"] if x["metal_index0"] == reference["metal_index0"]]
        if len(candidates) != 1: raise ValueError("Reference must identify exactly one original metal index")
        if "selected_donor_indices0" not in candidates[0]:
            raise ValueError("Donor perception failed; absent output is not an empty chemical reference")
        actual = set(candidates[0]["selected_donor_indices0"])
        expected = set(reference.get("donor_indices0", reference.get("expected_donor_indices0", [])))
        if "donor_indices0" not in reference and "expected_donor_indices0" not in reference:
            raise ValueError("A complete expected donor-index set is required")
        return {"complete_set_equal": actual == expected, "expected_count": len(expected),
                "observed_count": len(actual), "matched": len(actual & expected),
                "missing_donor_indices0": sorted(expected - actual),
                "additional_donor_indices0": sorted(actual - expected)}
    if route == "cif":
        candidates = [x for x in output["outcomes"]
                      if x["target"]["label"] == reference["target_label"]
                      and x["condition"] == reference["condition"]]
    elif route == "native":
        candidates = [x for x in output["native_physical_donor_spheres"]
                      if x["source_metal_label"] == reference["target_label"]]
    else: raise ValueError("Unsupported output route")
    if len(candidates) != 1:
        raise ValueError("Reference must uniquely identify target and condition; do not choose the best matching output")
    observed = candidates[0]["donor_images"]
    if route == "native" and not reference.get("include_metal_neighbours", False):
        observed = [x for x in observed if not x["is_metal"]]
    result = compare_images(reference["donor_images"], observed)
    result["typed_incidence_count"] = len(candidates[0]["record"]["incidences"]) if route == "cif" else None
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output-record", required=True, type=Path)
    p.add_argument("--references", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    a = p.parse_args()
    if a.output.exists(): raise FileExistsError("Choose a new comparison output path")
    record = json.loads(a.output_record.read_text(encoding="utf-8"))
    refs = json.loads(a.references.read_text(encoding="utf-8"))
    if isinstance(refs, dict): refs = refs.get("references", refs.get("cases"))
    if refs is None: raise ValueError("Reference document requires references or cases array")
    refs = [r for r in refs if "case_id" not in r or r["case_id"] == record["case_id"]]
    if not refs: raise ValueError("No reference matches this case_id")
    results=[]
    for ref in refs:
        row={"reference_id":ref.get("reference_id", ref.get("case_id"))}
        try: row.update(status="COMPARED", **compare_reference(record,ref))
        except Exception as exc: row.update(status="REFERENCE_OR_MAPPING_UNRESOLVED",error=type(exc).__name__+": "+str(exc))
        results.append(row)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    with a.output.open("x",encoding="utf-8") as f:
        json.dump({"comparisons":results,"references_are_not_encoder_inputs":True},f,indent=2)
        f.write("\n")


if __name__ == "__main__": main()
