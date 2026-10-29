"""Optional public-MOF-archive/CCDC integration test.

Set COORDREP_MOF_ZIP to the locally downloaded pinned collection archive.
No CIF content or site rows are written by the tested command.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


@unittest.skipUnless(os.environ.get("COORDREP_MOF_ZIP"), "COORDREP_MOF_ZIP is unset")
class AbadugIntegrationTest(unittest.TestCase):
    def test_public_zip_one_cif_commitments(self) -> None:
        archive = os.environ["COORDREP_MOF_ZIP"]
        result = subprocess.run(
            [sys.executable, "-B", str(ROOT / "run_abadug.py"), "--mof-zip", archive],
            cwd=str(ROOT),
            check=True,
            capture_output=True,
            text=True,
        )
        payload = json.loads(result.stdout)
        self.assertEqual(
            payload["status"],
            "PASS_ABADUG_public_ZIP_one_CIF_matches_frozen_v3_commitments",
        )
        self.assertTrue(all(payload["checks"].values()))
        self.assertEqual(payload["result"]["emitted_sites"], 6)

    def test_public_zip_ten_entry_batch_smoke(self) -> None:
        archive = os.environ["COORDREP_MOF_ZIP"]
        result = subprocess.run(
            [sys.executable, "-B", str(ROOT / "run_batch.py"), "--mof-zip", archive, "--limit", "10"],
            cwd=str(ROOT),
            check=True,
            capture_output=True,
            text=True,
        )
        payload = json.loads(result.stdout)
        self.assertEqual(payload["status"], "partial_smoke")
        self.assertEqual(payload["selected_entries"], 10)
        self.assertEqual(payload["totals"]["emitted_sites"], 46)
        self.assertEqual(payload["unique_metric_ids"], 17)
        self.assertEqual(payload["unique_stereo_ids"], 25)
        self.assertEqual(payload["private_frozen_QA"], "not_requested")
        parallel = subprocess.run(
            [
                sys.executable, "-B", str(ROOT / "run_batch.py"),
                "--mof-zip", archive, "--limit", "10", "--workers", "2",
            ],
            cwd=str(ROOT),
            check=True,
            capture_output=True,
            text=True,
        )
        parallel_payload = json.loads(parallel.stdout)
        self.assertEqual(parallel_payload["totals"], payload["totals"])
        self.assertEqual(parallel_payload["unique_metric_ids"], payload["unique_metric_ids"])
        self.assertEqual(parallel_payload["unique_stereo_ids"], payload["unique_stereo_ids"])
        self.assertEqual(
            parallel_payload["metric_to_stereo_multiplicity"],
            payload["metric_to_stereo_multiplicity"],
        )


if __name__ == "__main__":
    unittest.main()
