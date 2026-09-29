"""Pure aggregate-check tests; no CCDC installation or CIF is required."""

import unittest

from compare_public_summary import extract_headline


class PublicSummaryComparisonTests(unittest.TestCase):
    def test_extracts_primary_headline_without_site_rows(self) -> None:
        batch = {
            "status": "full_zip_run_requires_claim_review",
            "selected_entries": 15906,
            "start_index": 0,
            "totals": {"emitted_sites": 4},
            "unique_metric_ids": 2,
            "unique_stereo_ids": 3,
            "metric_to_stereo_multiplicity": {"1": 1, "2": 1},
            "state_bearing_entries_by_flag": {"Sohncke": 1, "not flagged Sohncke": 1},
            "paired_entries_by_flag": {"not flagged Sohncke": 1},
            "state_bearing_families_by_flag": {"Sohncke": 1, "not flagged Sohncke": 1},
            "paired_families_by_flag": {"Sohncke": 0, "not flagged Sohncke": 1},
        }
        observed = extract_headline(batch)
        self.assertEqual(observed["emitted_site_instances"], 4)
        self.assertEqual(observed["metric_ids_with_multiple_stereo_ids"], 1)
        self.assertEqual(observed["entries_with_within_entry_pair"], 1)
        self.assertEqual(
            observed["spacegroup_blind_join"]["not flagged Sohncke"]["paired_fraction_percent"],
            100.0,
        )

    def test_rejects_partial_output(self) -> None:
        with self.assertRaises(ValueError):
            extract_headline({"status": "partial_smoke", "selected_entries": 10})


if __name__ == "__main__":
    unittest.main()
