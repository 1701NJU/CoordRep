"""Execute every manifest row without dropping failures or substituting inputs.

Generated files may contain licensed source information. Keep the output local.
Chemical reference labels in the public register are not passed to the encoder.
"""
import argparse
import csv
import json
import os
from pathlib import Path
import subprocess
import sys
from run_input import sha

HERE = Path(__file__).resolve().parent


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--manifest", required=True, type=Path)
    p.add_argument("--input-directory", required=True, type=Path)
    p.add_argument("--output-directory", required=True, type=Path)
    p.add_argument("--database", default="CSD")
    p.add_argument("--source-release", default="user-local-CSD")
    p.add_argument("--timeout-seconds", type=int, default=120)
    a = p.parse_args()
    if a.timeout_seconds != 120:
        p.error("This panel protocol uses a fixed 120-second per-case budget")
    a.output_directory.mkdir(parents=True, exist_ok=False)
    with a.manifest.open(encoding="utf-8-sig", newline="") as f:
        rows = list(csv.DictReader(f))
    results = []
    for i, row in enumerate(rows):
        case_id, route = row["case_id"], row["route"]
        if not case_id or Path(case_id).name != case_id or "/" in case_id or "\\" in case_id:
            raise ValueError("case_id must be a nonempty filename-safe identifier")
        result = {"case_id": case_id, "route": route, "manifest_row0": i}
        output = a.output_directory / (str(i).zfill(3) + "_" + case_id + ".json")
        args = [sys.executable, "-B", str(HERE / "run_input.py"), route,
                "--case-id", case_id, "--output", str(output)]
        try:
            if route in ("xyz", "cif"):
                name = row["input_basename"]
                if not name or Path(name).name != name or "/" in name or "\\" in name:
                    raise ValueError("input_basename must identify one file in the input directory")
                source = a.input_directory / name
                expected = row.get("input_sha256", "").strip().lower()
                actual = sha(source)
                if not expected or actual != expected:
                    raise ValueError("Input checksum missing or different; input was not substituted")
                args += ["--input", str(source)]
                if route == "xyz": args += ["--metal-element", "Fe"]
                elif row.get("source_block"): args += ["--block", row["source_block"]]
            elif route == "native":
                args += ["--refcode", row["refcode"], "--database", a.database,
                         "--source-index", str(i), "--source-release", a.source_release]
            else:
                raise ValueError("Unknown manifest route")
            completed = subprocess.run(args, timeout=120, capture_output=True,
                text=True, encoding="utf-8", errors="replace",
                env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"})
            result.update(returncode=completed.returncode, stdout=completed.stdout,
                          stderr=completed.stderr,
                          status="ENCODER_RETURNED" if completed.returncode == 0 else "EXECUTION_FAILED_RETAINED")
            if output.exists():
                result["output_sha256"] = sha(output)
                # Successful execution is not scored as chemical agreement.
                result["chemical_agreement_scored"] = False
        except subprocess.TimeoutExpired:
            result["status"] = "TIMEOUT_RETAINED"
        except Exception as exc:
            result.update(status="INPUT_OR_EXECUTION_ISSUE_RETAINED",
                          error=type(exc).__name__ + ": " + str(exc))
        results.append(result)
        with (a.output_directory / ("receipt_" + str(i).zfill(3) + ".json")).open("x", encoding="utf-8") as f:
            json.dump(result, f, indent=2)
    with (a.output_directory / "panel_execution.json").open("x", encoding="utf-8") as f:
        json.dump({"schema": "CoordRep-public-panel-execution/1", "rows_attempted": len(rows),
            "rows_retained": len(results), "per_case_timeout_seconds": 120,
            "chemical_agreement_is_a_separate_comparison": True, "results": results}, f, indent=2)
    print(json.dumps({"rows_attempted": len(rows), "rows_retained": len(results)}))


if __name__ == "__main__": main()
