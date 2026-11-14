"""Fresh-process replay of an input command, keeping every generated file."""
import argparse,json,os,subprocess,sys
from pathlib import Path

HERE=Path(__file__).resolve().parent

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--directory",type=Path,required=True)
    p.add_argument("encoder_arguments",nargs=argparse.REMAINDER)
    a=p.parse_args()
    args=a.encoder_arguments
    if args and args[0]=="--":args=args[1:]
    if not args or "--output" in args: p.error("Provide route arguments after --, without --output")
    a.directory.mkdir(parents=True,exist_ok=False)
    payloads=[];receipts=[]
    for i in (1,2):
        output=a.directory/f"run{i}.json"
        command=[sys.executable,"-B",str(HERE/"run_input.py"),*args,"--output",str(output)]
        result=subprocess.run(command,encoding="utf-8",errors="replace",capture_output=True,
            env={**os.environ,"PYTHONDONTWRITEBYTECODE":"1"})
        receipts.append({"run":i,"returncode":result.returncode,"stdout":result.stdout,"stderr":result.stderr})
        if result.returncode==0:payloads.append(json.loads(output.read_text(encoding="utf-8")))
    passed=len(payloads)==2 and payloads[0]==payloads[1]
    (a.directory/"replay.json").write_text(json.dumps({"status":"PASS" if passed else "FAIL_RETAINED",
        "fresh_processes":2,"complete_payload_equal":passed,"runs":receipts},indent=2),encoding="utf-8")
    if not passed:raise SystemExit(1)
    print("PASS: both fresh-process payloads agree")

if __name__=="__main__":main()
