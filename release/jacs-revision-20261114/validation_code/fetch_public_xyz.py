"""Download external molecular inputs from the recorded immutable source pin.

No CSD files are downloaded. Third-party terms remain applicable; the CoordRep
MIT license does not relicense the external coordinates. Files are not replaced.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import ssl
from urllib.parse import urlparse
from urllib.request import urlopen, Request

PIN="839c5024787724df5b13850f41fe631b80017074"
HERE=Path(__file__).resolve().parent


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--directory",required=True,type=Path)
    p.add_argument("--ca-file",type=Path,help="Optional trusted PEM CA bundle; TLS verification remains enabled")
    a=p.parse_args()
    a.directory.mkdir(parents=True,exist_ok=False)
    with (HERE/"public_xyz_sources_44.csv").open(encoding="utf-8",newline="") as f:
        rows=list(csv.DictReader(f))
    if len(rows)!=44 or len({r["case_id"] for r in rows})!=44:
        raise ValueError("The pinned source manifest must retain all44 cases")
    results=[]
    for row in rows:
        receipt={"case_id":row["case_id"],"source_url":row["source_url"]}
        try:
            url=urlparse(row["source_url"])
            if url.scheme!="https" or url.netloc!="raw.githubusercontent.com" or "/"+PIN+"/" not in url.path:
                raise ValueError("Only the recorded immutable source is accepted")
            if Path(row["input_basename"]).name!=row["input_basename"]:
                raise ValueError("Expected a basename")
            context=ssl.create_default_context(cafile=str(a.ca_file)) if a.ca_file else ssl.create_default_context()
            with urlopen(Request(row["source_url"],headers={"User-Agent":"CoordRep-validation-input-fetch/1"}),timeout=60,context=context) as response:
                data=response.read()
            digest=hashlib.sha256(data).hexdigest()
            if digest!=row["input_sha256"]:
                raise ValueError("Downloaded bytes differ from the recorded input checksum")
            with (a.directory/row["input_basename"]).open("xb") as f:f.write(data)
            receipt.update(status="DOWNLOADED_CHECKSUM_MATCH",sha256=digest)
        except Exception as exc:
            receipt.update(status="DOWNLOAD_ISSUE_RETAINED",error=type(exc).__name__+": "+str(exc))
        results.append(receipt)
    with (a.directory/"download_receipts.json").open("x",encoding="utf-8") as f:
        json.dump({"upstream_commit":PIN,"results":results},f,indent=2)
    complete=sum(r["status"]=="DOWNLOADED_CHECKSUM_MATCH" for r in results)
    print(json.dumps({"downloaded_and_verified":complete,"expected":44}))
    if complete!=44:raise SystemExit(1)


if __name__=="__main__":main()
