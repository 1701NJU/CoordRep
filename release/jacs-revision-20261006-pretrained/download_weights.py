"""Download and verify pretrained release assets without Git LFS or credentials."""
from pathlib import Path, PurePosixPath
import argparse
import hashlib
import json
import shutil
import tempfile
import urllib.request
import zipfile

HERE = Path(__file__).resolve().parent

def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def download(output):
    destination = Path(output).resolve()
    if destination.exists():
        raise FileExistsError(f"Refusing to overwrite existing output: {destination}")
    assets = json.loads((HERE / "download_manifest.json").read_text(encoding="utf-8"))["assets"]
    manifest = json.loads((HERE / "manifest.json").read_text(encoding="utf-8"))
    models = {model["key"]: model for model in manifest["models"]}
    if len(assets) != len(models) or {asset["model_key"] for asset in assets} != set(models):
        raise ValueError("Download assets do not cover exactly the released models")
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="coordrep-download-", dir=destination.parent) as temporary:
        scratch = Path(temporary)
        extracted = scratch / "extracted"
        extracted.mkdir()
        for info in assets:
            archive = scratch / info["file"]
            print(f"Downloading {info['model_key']}", flush=True)
            request = urllib.request.Request(info["url"], headers={"User-Agent": "CoordRep-checkpoint-downloader"})
            with urllib.request.urlopen(request, timeout=60) as source, archive.open("wb") as target:
                shutil.copyfileobj(source, target)
            if archive.stat().st_size != info["bytes"] or sha256(archive) != info["sha256"]:
                raise ValueError("Downloaded archive failed size/SHA-256 validation")
            expected = {row["path"] for row in models[info["model_key"]]["files"]}
            expected.update({"manifest.json", "LICENSE"})
            with zipfile.ZipFile(archive) as bundle:
                names = [member.filename for member in bundle.infolist()]
                if len(names) != len(set(names)) or set(names) != expected:
                    raise ValueError("Archive membership does not match the release manifest")
                for member in bundle.infolist():
                    path = PurePosixPath(member.filename)
                    if path.is_absolute() or ".." in path.parts or "\\" in member.filename:
                        raise ValueError("Unsafe archive member path")
                    if (member.external_attr >> 16) & 0o170000 == 0o120000:
                        raise ValueError("Symlink archive members are not supported")
                recovered = json.loads(bundle.read("manifest.json").decode("utf-8"))
                if recovered != manifest:
                    raise ValueError("Downloaded model manifest does not match the Git manifest")
                if bundle.read("LICENSE") != (HERE / "LICENSE").read_bytes():
                    raise ValueError("Bundle license does not match the Git license")
                bundle.extractall(extracted)
        for model in manifest["models"]:
            for row in model["files"]:
                path = extracted / row["path"]
                if path.stat().st_size != row["bytes"] or sha256(path) != row["sha256"]:
                    raise ValueError(f"Bundle member failed validation: {row['path']}")
        if (extracted / "LICENSE").read_bytes() != (HERE / "LICENSE").read_bytes():
            raise ValueError("Bundle license does not match the Git license")
        extracted.rename(destination)
    print(json.dumps({"downloaded": str(destination), "models": len(models),
                      "verified_archive_sha256": [asset["sha256"] for asset in assets]}, indent=2))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="downloaded")
    download(parser.parse_args().output)
