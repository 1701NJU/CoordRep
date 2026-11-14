"""Verify public derived tables and hashes, without loading any structure file."""
from pathlib import Path
import csv
import hashlib
import json
import re


ROOT = Path(__file__).resolve().parent


def read_csv(name):
    with (ROOT / name).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def total(rows, key):
    return sum(int(row[key]) for row in rows)


def check(condition, label):
    if not condition:
        raise ValueError(label)


def main():
    register = read_csv("case_register_84.csv")
    fe = read_csv("molecular44_execution_and_CN.csv")
    molecular = read_csv("Table_S16_molecular_donor_references.csv")
    crystal = read_csv("Table_S17_crystallographic_references.csv")
    disorder = read_csv("Table_S18_disorder_conditional_references.csv")
    references = read_csv("references_S1_S16.csv")
    conditional = read_csv("conditional_source_inputs.csv")
    public_indices = json.loads((ROOT / "molecular_public_reference_indices.json").read_text(encoding="utf-8"))
    summary = json.loads((ROOT / "validation_summary.json").read_text(encoding="utf-8"))
    check(len(register) == len({r["case_id"] for r in register}) == 84, "84 unique cases")
    check(all(r["retained_in_fixed_panel"] == "True" for r in register), "No fixed case removed")
    check(sum(r["route"] == "xyz" for r in register) == 44, "44 XYZ inputs")
    check(sum(r["route"] == "native" for r in register) == 40, "40 native-source inputs")
    check(all(r["processing_and_fresh_replay"] == "COMPLETE" for r in register), "All case outputs replayed")
    for row in register:
        if row["route"] == "xyz":
            check(bool(re.fullmatch(r"[0-9a-f]{64}", row["input_sha256"])), "XYZ SHA-256")
            check(row["input_basename"] == Path(row["input_basename"]).name, "No input directory")
        else:
            check(not row["input_basename"] and not row["input_sha256"], "No licensed native input copy")
    check(len(fe) == 44 and total(fe, "local_Fe_records") == 47, "44 structures and 47 Fe records")
    check(sum(r["single_center_CN_equal"] == "True" for r in fe) == 40, "40 CN agreements")
    check([r["case_id"] for r in fe if r["single_center_CN_equal"] == "False"] == ["28"], "Case 28 CN difference retained")
    check(sum(r["single_center_CN_equal"] == "" for r in fe) == 3, "Three two-Fe cases")
    check(len(molecular) == 21 and total(molecular, "expected_donor_count") == 103, "21/103 molecular reference endpoint")
    check(total(molecular, "retained_donor_count") == 103, "103 molecular donors retained")
    check(total(molecular, "missing_donors") == total(molecular, "additional_donors") == 0, "Exact molecular sets")
    check(public_indices["index_base"] == 0 and len(public_indices["cases"]) == 21, "Public XYZ references")
    check(sum(len(r["expected_donor_indices0"]) for r in public_indices["cases"]) == 103, "103 public reference indices")
    check({r["case_id"] for r in public_indices["cases"]} == {r["case_id"] for r in molecular}, "Reference index case binding")
    case_map = {r["case_id"]: r for r in register}
    for row in public_indices["cases"]:
        check(row["input_sha256"] == case_map[row["case_id"]]["input_sha256"], "Reference input SHA binding")
        check(len(row["expected_donor_indices0"]) == len(set(row["expected_donor_indices0"])), "Distinct donor indices")
        check(row["metal_index0"] not in row["expected_donor_indices0"], "Metal excluded from donor set")
    check(len(crystal) == 12 and total(crystal, "qualified_metal_targets") == 19, "12 entries and 19 crystal targets")
    check(total(crystal, "complete_sets_equal") == 17, "17 exact crystal sets")
    check(total(crystal, "assigned_nonmetal_contacts") == total(crystal, "retained_nonmetal_contacts") == 107, "107 crystal contacts")
    check(total(crystal, "missing_nonmetal_donors") == 0 and total(crystal, "additional_nonmetal_donors") == 4, "Cu extras retained")
    check([r["refcode"] for r in crystal if int(r["additional_nonmetal_donors"])] == ["NUXNIG"], "Only Cu case has extras")
    check(total(crystal, "directed_metal_metal_relations_expected") == total(crystal, "directed_metal_metal_relations_retained") == 6, "Six directed Au-Tl relations")
    check(len(disorder) == 2 and total(disorder, "qualified_local_alternatives") == total(disorder, "complete_local_sets") == 16, "16 disorder alternatives")
    check(total(disorder, "expected_physical_donor_images") == total(disorder, "matched_physical_donor_images") == 72, "72 donor images")
    check(total(disorder, "missing_donors") == total(disorder, "additional_donors") == 0, "Exact conditional donor sets")
    check(total(disorder, "separately_retained_unscored_candidates") == 2, "Two unscored candidates retained")
    check(total(conditional, "expected_local_configurations") == 18 and total(conditional, "qualified_reference_configurations") == 16, "18 emitted, 16 reference-qualified")
    check([r["reference"] for r in references] == [f"S{i}" for i in range(1, 17)], "Sequential SI references")
    check(all(r["reference"] in {x["reference"] for x in references} for r in molecular + crystal + disorder), "All table references resolve")
    check(summary["fixed_panel"]["structures"] == 84, "Summary panel")
    check(summary["molecular_complete_references"]["assigned_donors"] == 103, "Summary molecular counts")
    check(summary["crystal_complete_references"]["exact_complete_sets"] == 17, "Summary crystal counts")
    check(summary["disorder_complete_references"]["retained_donor_images"] == 72, "Summary disorder counts")
    stability = summary["stability_controls"]
    check(stability["Fe_centre_level_trials"] == stability["Fe_centre_level_trials_passed"] == 470, "470 stability trials")
    check(stability["additional_case13_centre_level_trials"] == stability["additional_case13_centre_level_trials_passed"] == 2000, "2000 case 13 trials")
    check(stability["conditional_CIF_input_control_runs"] == 24 and stability["conditional_local_record_comparisons"] == 216, "24 runs versus 216 comparisons")
    check(not summary["scope"]["prospective_or_human_blind_benchmark"], "No blind benchmark claim")
    check(not summary["scope"]["whole_CSD_chemical_accuracy_estimated"], "No whole-CSD accuracy extrapolation")
    for line in (ROOT / "SHA256SUMS.txt").read_text(encoding="utf-8").splitlines():
        digest, name = line.split("  ", 1)
        path = ROOT / name
        check(path.resolve().parent == ROOT, "Checksum target stays in package")
        check(hashlib.sha256(path.read_bytes()).hexdigest() == digest, "Checksum mismatch: " + name)
    print(json.dumps({"status": "PASS", "fixed_cases": 84, "molecular_complete_sets": 21,
                      "crystal_exact_sets": "17/19", "conditional_complete_sets": 16,
                      "conditional_donor_images": 72, "restricted_structure_inputs_loaded": False}))


if __name__ == "__main__":
    main()
