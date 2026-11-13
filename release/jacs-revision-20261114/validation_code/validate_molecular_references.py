"""Recompute the public molecular reference subset from user-downloaded XYZs.

The published index mappings apply only to exact-checksum external molecular
inputs. They do not supply or stand in for restricted CSD crystal references.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

from run_input import HERE, sha, write_new
from compare_local_sets import compare_reference


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--references",required=True,type=Path)
    p.add_argument("--input-directory",required=True,type=Path)
    p.add_argument("--output-directory",required=True,type=Path)
    a=p.parse_args()
    a.output_directory.mkdir(parents=True,exist_ok=False)
    document=json.loads(a.references.read_text(encoding="utf-8-sig"))
    if document.get("index_base") != 0: raise ValueError("Only declared zero-based references are accepted")
    results=[]
    for i,ref in enumerate(document["cases"]):
        row={"case_id":ref["case_id"],"reference_is_encoder_input":False}
        try:
            name=ref["input_basename"]
            if not name or Path(name).name != name or "/" in name or "\\" in name:
                raise ValueError("Input name must be a basename")
            path=a.input_directory/name
            if sha(path)!=ref["input_sha256"].lower():
                raise ValueError("Input checksum differs; reference atom indices cannot be transferred")
            output=a.output_directory/(str(i).zfill(3)+"_encoded.json")
            command=[sys.executable,"-B",str(HERE/"run_input.py"),"xyz","--input",str(path),
                     "--metal-element","Fe","--case-id",ref["case_id"],"--output",str(output)]
            process=subprocess.run(command,capture_output=True,text=True,encoding="utf-8",errors="replace",
                timeout=120,env={**os.environ,"PYTHONDONTWRITEBYTECODE":"1"})
            row.update(returncode=process.returncode,stdout=process.stdout,stderr=process.stderr)
            if process.returncode: raise RuntimeError("Encoder failed; see captured stderr")
            encoded=json.loads(output.read_text(encoding="utf-8"))
            row.update(status="COMPARED",**compare_reference(encoded,ref))
        except subprocess.TimeoutExpired:
            row.update(status="TIMEOUT_RETAINED",complete_set_equal=False)
        except Exception as exc:
            row.update(status="INPUT_OR_COMPARISON_ISSUE_RETAINED",complete_set_equal=False,
                       error=type(exc).__name__+": "+str(exc))
        results.append(row)
    matched=sum(r.get("complete_set_equal",False) for r in results)
    summary={"schema":"CoordRep-public-molecular-reference-replay/1",
             "reference_count":len(results),"complete_donor_set_matches":matched,
             "all_cases_retained":True,"references_are_not_encoder_inputs":True,
             "status":"PASS" if matched==len(results) else "MISMATCH_OR_ISSUE_RETAINED",
             "results":results}
    write_new(a.output_directory/"molecular_validation.json",summary)
    print(json.dumps({k:v for k,v in summary.items() if k!="results"},indent=2))
    if matched!=len(results):raise SystemExit(1)


if __name__=="__main__":main()
