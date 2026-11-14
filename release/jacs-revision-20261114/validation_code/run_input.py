"""Configurable wrappers around the unchanged, versioned input implementations.

Outputs can contain coordinates and detailed mappings from the user's private
inputs. Keep such output local; the public release contains no deposited CIFs.
"""
from __future__ import annotations
import argparse
import dataclasses
import hashlib
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.dont_write_bytecode = True
sys.path[:0] = [str(HERE / "frozen_xyz3_runtime"), str(HERE)]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_new(path, value):
    path = Path(path)
    if path.exists():
        raise FileExistsError("Use a new output path; previous results are retained")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, ensure_ascii=True, allow_nan=False,
                  default=lambda x: x.item() if hasattr(x, "item") else str(x))
        handle.write("\n")


def encode_xyz(path, elements=None, case_id=None):
    from coordrep import __version__, encode_molecule, CoordRepConfig
    from coordrep.encode import _load_molecule
    from coordrep.io.tmqm_reader import TRANSITION_METALS
    from coordrep.graph.neighbors import perceive_neighbors
    from coordrep.graph.build_graph import GraphBuilder
    path = Path(path)
    mol = _load_molecule(path)
    eligible = set(elements) if elements else TRANSITION_METALS
    sites = []
    for index, atom in enumerate(mol.atoms):
        if atom.element not in eligible:
            continue
        site = {"metal_index0": index, "element": atom.element}
        try:
            neighbours, trace = perceive_neighbors(mol, index)
            site.update(perception_trace=trace,
                        selected_donor_indices0=[n.atom_idx for n in neighbours],
                        atom_contact_cn=len(neighbours))
            cc = encode_molecule(mol, config=CoordRepConfig.default(), metal_index=index)
            site["encoded_record"] = cc.to_dict()
            site["record_issues"] = [dataclasses.asdict(i) for i in cc.issues]
            site["ligand_payload_provenance"] = [l.payload_provenance for l in cc.ligands]
            graph = GraphBuilder().build(mol, metal_index=index)
            site["graph_metadata_metal_index0"] = graph.metal_idx
            site["graph_selected_donor_indices0"] = graph.donor_indices
            site["graph_edges"] = sorted([sorted(e) for e in graph.molecular_graph.edges])
            site["status"] = "ENCODED" if cc.graph is not None else "CONTROLLED_NO_COORDINATION_RECORD"
            canonical = cc.canonicalize()
            site["canonical_record"] = canonical.to_dict()
            site["coordrep_string"] = canonical.to_string()
            site["canonicalization_status"] = "COMPLETE"
        except Exception as exc:
            site.setdefault("status", "SITE_EXECUTION_FAILED")
            if site.get("canonicalization_status") != "COMPLETE":
                site["canonicalization_status"] = "FAILED_RETAINED"
            site["error"] = type(exc).__name__ + ": " + str(exc)
        sites.append(site)
    return {"schema": "CoordRep-public-input-run/1", "route": "xyz",
            "case_id": case_id or path.stem, "adapter_version": __version__,
            "input_sha256": sha(path), "source_name": path.name,
            "reference_fields_supplied_to_encoder": False, "BO_supplied": False,
            "local_metal_count": len(sites), "sites": sites,
            "scope": "LOCAL_CENTRE_RECORDS_NOT_WHOLE_MULTIMETAL_ASSEMBLY_IDENTITY"}


def encode_cif(path, block=None, case_id=None):
    from coordrep_cifcfg.adapter import encode_cif as frozen_encode_cif
    path = Path(path)
    result = frozen_encode_cif(path, block)
    result.update(schema="CoordRep-public-input-run/1", route="cif",
                  case_id=case_id or path.stem, input_sha256=sha(path),
                  source_name=path.name)
    return result


def native_spheres(entry, adapter):
    """Author-owned audit export; the returned coordinates must remain local."""
    import numpy as np
    def vector(atom, name):
        try:
            x = getattr(atom, name)
            return None if x is None else [float(x.x), float(x.y), float(x.z)]
        except Exception:
            return None
    atoms, bonds, adjacency = adapter._collect_graph(entry.molecule)
    metals = [i for i, a in enumerate(atoms) if adapter._is_metal_atom(a)]
    colours = adapter._wl_colours(atoms, adjacency)
    drafts, _, external = adapter._build_group_drafts(atoms, bonds, adjacency, metals)
    keyed = [(adapter._metal_sort_key(i, atoms, colours, drafts, tuple(external)), i) for i in metals]
    order = [i for _, i in sorted(keyed, key=lambda pair: (pair[0], pair[1]))]
    mapping = {i: "M" + str(j + 1) for j, i in enumerate(order)}
    spheres = []
    for mi in order:
        centre = atoms[mi]
        mc = vector(centre, "coordinates")
        neighbours = []
        for bond in bonds:
            if mi not in (bond.a, bond.b):
                continue
            j = bond.b if bond.a == mi else bond.a
            atom = atoms[j]
            dc = vector(atom, "coordinates")
            v = None if mc is None or dc is None else [dc[k] - mc[k] for k in range(3)]
            delta = adapter._relative_image_delta(centre, atom)
            neighbours.append({"source_atom_index0": j, "source_atom_label": str(atom.label),
                "element": str(atom.atomic_symbol), "is_metal": j in mapping,
                "target_metal_record_label": mapping.get(j), "bond_type": bond.bond_type,
                "vector_A": v, "distance_A": None if v is None else float(np.linalg.norm(v)),
                "relative_translation": delta, "translation_resolved": delta is not None,
                "occupancy": getattr(atom, "occupancy", None)})
        spheres.append({"record_metal_label": mapping[mi], "source_metal_label": str(centre.label),
                        "element": str(centre.atomic_symbol), "donor_images": neighbours})
    return spheres


def encode_native(refcode, database="CSD", source_index=0, source_release="user-local-CSD", case_id=None):
    from ccdc.io import EntryReader
    import coordrep.audit.csd_adapter as adapter
    reader = EntryReader(database)
    entry = reader.entry(refcode)
    record = adapter.adapt_csd_entry(entry, source_index, source_release=source_release)
    return {"schema": "CoordRep-public-input-run/1", "route": "native",
            "case_id": case_id or refcode, "refcode": refcode,
            "record": record.to_dict(), "native_physical_donor_spheres": native_spheres(entry, adapter),
            "reference_fields_supplied_to_encoder": False,
            "scope": "SOURCE_DERIVED_RELATIONS_CHEMICAL_REFERENCE_COMPARISON_IS_SEPARATE"}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="route", required=True)
    for route in ("xyz", "cif", "native"):
        q = sub.add_parser(route)
        q.add_argument("--output", required=True, type=Path)
        q.add_argument("--case-id")
        if route in ("xyz", "cif"):
            q.add_argument("--input", required=True, type=Path)
        if route == "xyz": q.add_argument("--metal-element", action="append")
        if route == "cif": q.add_argument("--block")
        if route == "native":
            q.add_argument("--refcode", required=True)
            q.add_argument("--database", default="CSD")
            q.add_argument("--source-index", type=int, default=0)
            q.add_argument("--source-release", default="user-local-CSD")
    args = p.parse_args()
    if args.route == "xyz": value = encode_xyz(args.input, args.metal_element, args.case_id)
    elif args.route == "cif": value = encode_cif(args.input, args.block, args.case_id)
    else: value = encode_native(args.refcode, args.database, args.source_index, args.source_release, args.case_id)
    write_new(args.output, value)
    print(json.dumps({"route": args.route, "case_id": value["case_id"],
                      "output_sha256": sha(args.output), "output_saved": True}))


if __name__ == "__main__": main()
