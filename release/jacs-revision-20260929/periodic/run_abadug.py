#!/usr/bin/env python3
"""Public, no-row-output periodic-v3 worked example for ABADUG.

Reads one P1 CIF from the user's CCDC MOF Collection ZIP, uses their local
CCDC Python API to obtain the explicit molecular bond graph, reconstructs
translation-labelled quotient donor edges, and applies the exact frozen
CR-PLS/2 MID/SID canonicalizer. No CIF or per-site rows are written.

This is a deliberately narrow ABADUG demonstration, not the full v2/v3
eligibility waterfall or a whole-collection census runner.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import zipfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np

from periodic_v3 import (
    PeriodicDonorImage,
    canonical_entry_multiset_id,
    canonical_periodic_groups,
    canonicalize_fractional_site,
)


ARCHIVE_SHA256 = "ABF9F2A7288FD83E2A59F19B717C202414E744A1B2C8C1F9EEC1F909554346FA"
CANONICALIZER_SHA256 = "292AFD4DB89D1DD552F9A03433B207FD648BB30E9DED717BFB8DBCE17D4694AC"
QUOTIENT_SHA256 = "8C9DE81E1EBABDECE083BD6393786357BDC35CE3D6AD3142AAB205E2ECB19EEA"
TARGET_REFCODE = "ABADUG"
TRANSITION_METALS = {
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "La", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
}
DONOR_ELEMENTS = {"N", "O", "S", "P", "Cl", "Br", "I", "F", "C", "Se", "Te", "As"}

# Commitments computed by the original frozen v3 run. They are hashes of
# derived identity payloads, not a copy of source coordinates or site rows.
EXPECTED = {
    "periodic_atom_orbits": 282,
    "emitted_sites": 6,
    "sites_with_nonzero_translation": 3,
    "entry_mid_multiset_sha256": "fd1969008e6147f3077418b22588d2084bfb4767452f2afce82ef71e7fadbe61",
    "entry_sid_multiset_sha256": "b8d4823c52c4435e0be8052177f78b50e7187e48fc528ebf31b4552a9bb67f8d",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest().upper()


def _load_target_cif(archive: zipfile.ZipFile) -> Tuple[str, Dict[str, str]]:
    details = archive.read("CSD_MOF_Collection/Framework details.csv")
    records = csv.DictReader(io.StringIO(details.decode("utf-8-sig")))
    selected = [
        row for row in records if row["CSD refcode"].strip().upper() == TARGET_REFCODE
    ]
    if len(selected) != 1:
        raise ValueError(f"Expected one Framework details row for {TARGET_REFCODE}")
    row = selected[0]
    stem = row["CIF filename"].strip()
    matches = [name for name in archive.namelist() if Path(name).stem.lower() == stem.lower()]
    if len(matches) != 1:
        raise ValueError(f"Expected one CIF member for {stem}, found {len(matches)}")
    if row["Unreliable chemistry"].strip().lower() not in {"", "-", "0", "false", "no"}:
        raise ValueError("ABADUG has an unexpected unreliable-chemistry flag")
    return archive.read(matches[0]).decode("utf-8"), row


def _edge_key(donor_group: int, donor_shift: Tuple[int, ...], metal_shift: Tuple[int, ...]):
    return (donor_group,) + tuple(
        int(donor_shift[index]) - int(metal_shift[index]) for index in range(3)
    )


def calculate(crystal: Any) -> Dict[str, Any]:
    molecule = crystal.molecule
    if molecule is None:
        raise ValueError("CCDC returned a crystal with no molecule")
    cell = np.asarray(crystal.fractional_to_orthogonal.rotation, dtype=float)
    if cell.shape != (3, 3) or not np.all(np.isfinite(cell)) or abs(np.linalg.det(cell)) < 1e-10:
        raise ValueError("Invalid crystallographic cell transform")
    atom_info, atom_to_group, groups = canonical_periodic_groups(molecule, cell)
    metal_groups = sorted(
        group_id for group_id, group in groups.items()
        if group["element"] in TRANSITION_METALS
    )
    edges: Dict[int, Dict[Tuple[int, int, int, int], Dict[str, Any]]] = defaultdict(dict)
    for bond in molecule.bonds:
        left, right = bond.atoms
        orientations = []
        if left.atomic_symbol in TRANSITION_METALS and right.atomic_symbol in DONOR_ELEMENTS:
            orientations.append((left, right))
        if right.atomic_symbol in TRANSITION_METALS and left.atomic_symbol in DONOR_ELEMENTS:
            orientations.append((right, left))
        bond_type = str(bond.bond_type or "")
        if orientations and ("pi" in bond_type.lower() or "deloc" in bond_type.lower()):
            raise ValueError("ABADUG includes a haptic/delocalized metal-donor bond")
        for metal_atom, donor_atom in orientations:
            metal_info = atom_info[metal_atom.index]
            donor_info = atom_info[donor_atom.index]
            metal_group = atom_to_group[metal_atom.index]
            donor_group = atom_to_group[donor_atom.index]
            key = _edge_key(
                donor_group, donor_info["translation"], metal_info["translation"]
            )
            observation = np.asarray(donor_info["cartesian"]) - np.asarray(metal_info["cartesian"])
            record = edges[metal_group].setdefault(
                key, {"element": donor_atom.atomic_symbol, "observations": []}
            )
            if record["element"] != donor_atom.atomic_symbol:
                raise ValueError("Conflicting donor element for one quotient edge")
            record["observations"].append(observation)

    metric_ids: List[str] = []
    stereo_ids: List[str] = []
    cn_counts: Counter[int] = Counter()
    translated_sites = 0
    for metal_group in metal_groups:
        keys = sorted(edges.get(metal_group, {}))
        cn = len(keys)
        if cn not in {2, 3, 4, 5, 6}:
            raise ValueError(f"Unexpected ABADUG coordination number {cn}")
        metal = groups[metal_group]
        metal_fractional = np.asarray(metal["wrapped_fractional"], dtype=float)
        donors = []
        has_translation = False
        for key in keys:
            donor_group = groups[key[0]]
            donor_fractional = np.asarray(donor_group["wrapped_fractional"], dtype=float)
            vector = (donor_fractional + np.asarray(key[1:], dtype=float) - metal_fractional) @ cell
            distance = float(np.linalg.norm(vector))
            if not math.isfinite(distance) or distance < 0.60 or distance > 3.50:
                raise ValueError("ABADUG donor edge violates frozen distance gate")
            for observation in edges[metal_group][key]["observations"]:
                if float(np.linalg.norm(observation - vector)) > 0.001:
                    raise ValueError("Periodic quotient vector witness failed")
            donors.append(
                PeriodicDonorImage(
                    element=edges[metal_group][key]["element"],
                    fractional_base=tuple(float(value) for value in donor_fractional),
                    translation=tuple(int(value) for value in key[1:]),
                    orbit=int(key[0]),
                    bond_role="sigma",
                )
            )
            has_translation |= any(value != 0 for value in key[1:])
        state = canonicalize_fractional_site(
            metal=str(metal["element"]),
            metal_fractional_base=tuple(float(value) for value in metal_fractional),
            metal_translation=(0, 0, 0),
            donors=donors,
            cell_matrix=cell,
            decimals=6,
        )
        metric_ids.append(state.metric_id)
        stereo_ids.append(state.stereo_id)
        cn_counts[cn] += 1
        translated_sites += int(has_translation)
    return {
        "periodic_atom_orbits": len(groups),
        "emitted_sites": len(metric_ids),
        "sites_with_nonzero_translation": translated_sites,
        "coordination_numbers": {str(key): value for key, value in sorted(cn_counts.items())},
        "entry_mid_multiset_sha256": canonical_entry_multiset_id(metric_ids),
        "entry_sid_multiset_sha256": canonical_entry_multiset_id(
            stereo_ids, version="CR-PERIODIC-STEREO-ATLAS/1"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mof-zip", required=True, type=Path)
    args = parser.parse_args()
    if not args.mof_zip.is_file():
        parser.error(f"MOF archive not found: {args.mof_zip}")
    archive_hash = sha256(args.mof_zip)
    if archive_hash != ARCHIVE_SHA256:
        raise ValueError("MOF archive SHA-256 differs from frozen v3 protocol")
    package_root = Path(__file__).resolve().parent / "periodic_v3"
    if sha256(package_root / "canonical_periodic_site.py") != CANONICALIZER_SHA256:
        raise ValueError("Vendored canonicalizer differs from frozen source")
    if sha256(package_root / "canonical_csd_quotient.py") != QUOTIENT_SHA256:
        raise ValueError("Vendored quotient builder differs from frozen source")
    try:
        import ccdc
        from ccdc.crystal import Crystal
    except ImportError as exc:
        raise SystemExit("The CCDC Python API is required in the selected Python environment") from exc
    with zipfile.ZipFile(args.mof_zip) as archive:
        cif_text, _ = _load_target_cif(archive)
    crystal = Crystal.from_string(cif_text, format="cif")
    observed = calculate(crystal)
    checks = {key: observed[key] == value for key, value in EXPECTED.items()}
    if not all(checks.values()):
        raise AssertionError({"checks": checks, "observed": observed})
    print(json.dumps({
        "status": "PASS_ABADUG_public_ZIP_one_CIF_matches_frozen_v3_commitments",
        "refcode": TARGET_REFCODE,
        "archive_sha256": archive_hash,
        "ccdc_api": ccdc.__version__,
        "result": observed,
        "checks": checks,
        "scope": "single-entry smoke, not full-collection reproduction",
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
