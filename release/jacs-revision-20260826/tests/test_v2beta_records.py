"""Test v2-beta Box 1 records (Examples 4–5) validate correctly."""

import json
from pathlib import Path

import pytest

RESULTS = Path(__file__).resolve().parent.parent / "revision_results"
V2BETA_DIR = RESULTS / "box1_v2beta_examples"


@pytest.fixture(scope="module")
def v2beta_records():
    path = V2BETA_DIR / "box1_v2beta_full_records.jsonl"
    if not path.exists():
        pytest.skip("v2beta records JSONL not found")
    records = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line:
                records.append(json.loads(line))
    return records


@pytest.fixture(scope="module")
def validation_report():
    path = V2BETA_DIR / "box1_v2beta_validation_report.csv"
    if not path.exists():
        pytest.skip("v2beta validation report not found")
    import csv
    with open(path, encoding="utf-8") as f:
        return list(csv.DictReader(f))


class TestV2BetaRecordsExist:
    def test_records_file_exists(self):
        assert (V2BETA_DIR / "box1_v2beta_full_records.jsonl").exists()

    def test_validation_file_exists(self):
        assert (V2BETA_DIR / "box1_v2beta_validation_report.csv").exists()

    def test_identity_keys_file_exists(self):
        assert (V2BETA_DIR / "box1_v2beta_identity_keys.csv").exists()

    def test_display_file_exists(self):
        assert (V2BETA_DIR / "box1_v2beta_display_records.txt").exists()

    def test_selection_notes_exist(self):
        assert (V2BETA_DIR / "box1_v2beta_selection_notes.md").exists()


class TestExample4Multinuclear:
    def test_record_present(self, v2beta_records):
        ex4 = [r for r in v2beta_records if r["example_id"] == "example_4"]
        assert len(ex4) == 1

    def test_refcode(self, v2beta_records):
        ex4 = [r for r in v2beta_records if r["example_id"] == "example_4"][0]
        assert ex4["refcode"] == "ETOREN"

    def test_record_type(self, v2beta_records):
        ex4 = [r for r in v2beta_records if r["example_id"] == "example_4"][0]
        assert ex4["record_type"] == "molecular_multinuclear_graph"

    def test_two_metals(self, v2beta_records):
        ex4 = [r for r in v2beta_records if r["example_id"] == "example_4"][0]
        assert len(ex4["metals"]) == 2

    def test_metals_are_cu(self, v2beta_records):
        ex4 = [r for r in v2beta_records if r["example_id"] == "example_4"][0]
        elements = {m["element"] for m in ex4["metals"]}
        assert elements == {"Cu"}

    def test_has_bridges(self, v2beta_records):
        ex4 = [r for r in v2beta_records if r["example_id"] == "example_4"][0]
        bridges = [s for s in ex4["sites"] if s["mu"] > 1]
        assert len(bridges) >= 2, "Expected at least 2 bridge sites"

    def test_bridge_not_duplicated(self, v2beta_records):
        ex4 = [r for r in v2beta_records if r["example_id"] == "example_4"][0]
        bridges = [s for s in ex4["sites"] if s["mu"] > 1]
        for b in bridges:
            assert len(b["target_metals"]) == b["mu"], (
                f"Bridge {b['label']}: mu={b['mu']} but "
                f"target_metals={b['target_metals']}"
            )

    def test_local_sphere_consistency(self, v2beta_records):
        ex4 = [r for r in v2beta_records if r["example_id"] == "example_4"][0]
        for m in ex4["metals"]:
            local = [s for s in ex4["sites"]
                     if m["label"] in s["target_metals"]]
            assert len(local) == m["cn_site"], (
                f"{m['label']}: cn_site={m['cn_site']} but "
                f"found {len(local)} local sites"
            )

    def test_identity_keys_present(self, v2beta_records):
        ex4 = [r for r in v2beta_records if r["example_id"] == "example_4"][0]
        identity = ex4["identity"]
        for key in ["L0_GlobalState", "L1_GlobalShape",
                     "L2_MetalGraphTopo", "L3_Connectivity"]:
            assert identity.get(key), f"Missing identity key: {key}"

    def test_no_mm_bond_claimed(self, v2beta_records):
        ex4 = [r for r in v2beta_records if r["example_id"] == "example_4"][0]
        for edge in ex4.get("metal_edges", []):
            assert edge["mm_bond"] == "no", (
                "ETOREN should not claim a direct Cu-Cu bond"
            )


class TestExample5Haptic:
    def test_record_present(self, v2beta_records):
        ex5 = [r for r in v2beta_records if r["example_id"] == "example_5"]
        assert len(ex5) == 1

    def test_refcode(self, v2beta_records):
        ex5 = [r for r in v2beta_records if r["example_id"] == "example_5"][0]
        assert ex5["refcode"] == "FEROCE01"

    def test_record_type(self, v2beta_records):
        ex5 = [r for r in v2beta_records if r["example_id"] == "example_5"][0]
        assert ex5["record_type"] == "haptic_pi_site_object"

    def test_single_metal_fe(self, v2beta_records):
        ex5 = [r for r in v2beta_records if r["example_id"] == "example_5"][0]
        assert ex5["metal"]["element"] == "Fe"

    def test_cn_site_2(self, v2beta_records):
        ex5 = [r for r in v2beta_records if r["example_id"] == "example_5"][0]
        assert ex5["metal"]["cn_site"] == 2, (
            "Ferrocene should have CN_site=2 (two Cp sites)"
        )

    def test_eta_sum_10(self, v2beta_records):
        ex5 = [r for r in v2beta_records if r["example_id"] == "example_5"][0]
        assert ex5["metal"]["eta_sum"] == 10, (
            "Ferrocene should have eta_sum=10 (5+5)"
        )

    def test_two_pi_sites(self, v2beta_records):
        ex5 = [r for r in v2beta_records if r["example_id"] == "example_5"][0]
        pi_sites = [s for s in ex5["sites"] if s["site_type"] == "pi_fragment"]
        assert len(pi_sites) == 2

    def test_each_site_eta5(self, v2beta_records):
        ex5 = [r for r in v2beta_records if r["example_id"] == "example_5"][0]
        for s in ex5["sites"]:
            assert s["eta"] == 5, f"Site {s['label']} should be eta=5"

    def test_site_atoms_contiguous_5(self, v2beta_records):
        ex5 = [r for r in v2beta_records if r["example_id"] == "example_5"][0]
        for s in ex5["sites"]:
            assert len(s["donor_atoms"]) == 5, (
                f"Site {s['label']}: expected 5 donor atoms for eta5-Cp"
            )

    def test_not_expanded_as_eta1(self, v2beta_records):
        ex5 = [r for r in v2beta_records if r["example_id"] == "example_5"][0]
        # Should have exactly 2 sites (not 10 eta1 sites)
        assert len(ex5["sites"]) == 2, (
            f"Expected 2 sites for ferrocene, got {len(ex5['sites'])}. "
            "Cp rings must not be expanded as individual eta1 donors."
        )

    def test_identity_keys_present(self, v2beta_records):
        ex5 = [r for r in v2beta_records if r["example_id"] == "example_5"][0]
        identity = ex5["identity"]
        for key in ["L0_HapticState", "L1_HapticShape",
                     "L2_SiteTopo", "L3_Connectivity"]:
            assert identity.get(key), f"Missing identity key: {key}"


class TestValidationReport:
    def test_all_checks_pass(self, validation_report):
        for row in validation_report:
            status = row["status"]
            # Status is either "true" for boolean checks, or a
            # descriptive string for label checks (e.g. scope_label)
            if status in ("true", "false"):
                assert status == "true", (
                    f"{row['example_id']} / {row['check_name']}: "
                    f"expected true, got {status}"
                )
            else:
                # Non-boolean status (e.g. scope label value) — just
                # verify it's not empty
                assert len(status) > 0, (
                    f"{row['example_id']} / {row['check_name']}: "
                    f"empty status value"
                )
