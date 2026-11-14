"""Source-coordinate-only local configurations; no reference bond tables read."""
from dataclasses import dataclass
from collections import Counter
from itertools import product, combinations
from pathlib import Path
import hashlib
import json
import math
import re
import numpy as np
import gemmi

from coordrep.audit.metal_policy import ALL_METAL_SET
from coordrep.audit.records import (DonorGroup, CoordinationIncidence, DonorPairRelation,
    MetalSiteRecord, EntryRecord, validate_entry_record)
from coordrep.graph.neighbors import (DONOR_ELEMENTS, covalent_radius,
    COVALENT_RADIUS_FACTOR, RADIUS_BOUNDARY_TOLERANCE_A, GABRIEL_MARGIN_TOLERANCE_A2)
from coordrep.geometry.shape import ShapeCalculator
from frozen_canonical_periodic_site import canonicalize_cartesian_site

VERSION = "cifcfg1.20261113"
POSITION_DEDUP_DECIMALS = 10
MAX_LOCAL_CONDITIONS = 4096
ALLOWED_COLUMNS = (
    "label", "type_symbol", "fract_x", "fract_y", "fract_z", "occupancy",
    "disorder_assembly", "disorder_group",
)


def clean(raw):
    if raw is None:
        return None
    token = str(raw)
    if token in (".", "?"):
        return token
    return gemmi.cif.as_string(token)


def num(raw):
    x = gemmi.cif.as_number(str(raw))
    if not math.isfinite(x):
        raise ValueError("Missing or invalid required numerical CIF value")
    return float(x)


def canonical_json(value):
    return json.dumps(value, ensure_ascii=True, sort_keys=True, separators=(",", ":"), allow_nan=False)


def digest(value):
    return hashlib.sha256(canonical_json(value).encode()).hexdigest()


@dataclass
class Snapshot:
    block: str
    atoms: list
    matrix: np.ndarray
    operations: list
    issues: list


def read_snapshot(path, block_name=None):
    doc = gemmi.cif.read_file(str(path))
    if block_name is None:
        if len(doc) != 1:
            raise ValueError("Specify a block for a multiple-block source CIF")
        block = doc[0]
    else:
        block = doc.find_block(block_name)
        if block is None:
            raise ValueError("Requested source block is absent")
    cell = gemmi.UnitCell(*[num(block.find_value("_cell_" + k)) for k in
        ("length_a", "length_b", "length_c", "angle_alpha", "angle_beta", "angle_gamma")])
    matrix = np.column_stack([list(cell.orthogonalize(gemmi.Fractional(*v))) for v in np.eye(3)])
    if np.linalg.det(matrix) <= 0:
        raise ValueError("A finite right-handed cell is required")
    columns = {k: list(block.find_values("_atom_site_" + k)) for k in ALLOWED_COLUMNS}
    n = len(columns["label"])
    if not n or any(len(columns[k]) != n for k in ("type_symbol", "fract_x", "fract_y", "fract_z")):
        raise ValueError("Complete labelled fractional atom coordinates are required")
    for k, values in columns.items():
        if values and len(values) != n:
            raise ValueError("Inconsistent atom-site column length: " + k)
    atoms, issues = [], []
    for i in range(n):
        get = lambda k, default=None: clean(columns[k][i]) if columns[k] else default
        symbol = get("type_symbol")
        match = re.match(r"^([A-Z][a-z]?)", symbol or "")
        if match is None or gemmi.Element(match.group(1)).atomic_number <= 0:
            raise ValueError("Unknown source element")
        occ_text = get("occupancy")
        occupancy = None if occ_text in (None, "", ".", "?") else num(occ_text)
        if occupancy is not None and not 0 <= occupancy <= 1.000001:
            raise ValueError("Occupancy outside zero to one")
        assembly, group = get("disorder_assembly", "."), get("disorder_group", ".")
        if group not in (".", "?", "", None) and assembly in (".", "?", "", None):
            assembly = "__unspecified_single_assembly__"
            issues.append("DISORDER_ASSEMBLY_UNSPECIFIED")
        atoms.append({"source_index0": i, "label": get("label"), "element": match.group(1),
            "fractional": [num(columns["fract_" + k][i]) for k in "xyz"],
            "occupancy": occupancy, "occupancy_source_token": occ_text,
            "assembly": assembly, "group": group})
    if len({a["label"] for a in atoms}) != n:
        raise ValueError("Source atom labels must be unique within a block")
    op_strings = list(block.find_values("_space_group_symop_operation_xyz")) or list(
        block.find_values("_symmetry_equiv_pos_as_xyz"))
    if op_strings:
        operations = [gemmi.Op(clean(s)) for s in op_strings]
    else:
        name = clean(block.find_value("_space_group_name_H-M_alt") or block.find_value("_symmetry_space_group_name_H-M"))
        if not name:
            raise ValueError("Source symmetry operations or a space-group name are required")
        operations = list(gemmi.SpaceGroup(name).operations())
        issues.append("SYMMETRY_OPERATIONS_FROM_SOURCE_SPACE_GROUP_NAME")
    # Deduplicate equivalent source operations without using their source order.
    operations = [gemmi.Op(s) for s in sorted({op.triplet() for op in operations})]
    return Snapshot(block.name, atoms, matrix, operations, sorted(set(issues)))


def expanded_local_images(snapshot, target, radius):
    """All physical images within radius, including the full special-site orbit."""
    center = np.asarray(target["fractional"], dtype=float)
    bound = radius * np.linalg.norm(np.linalg.inv(snapshot.matrix), axis=1)
    found = {}
    for atom in snapshot.atoms:
        if atom["element"] == "H" or atom["occupancy"] == 0:
            continue
        for op in snapshot.operations:
            f = np.asarray(op.apply_to_xyz(atom["fractional"]), dtype=float)
            wrapped = f - np.floor(f + 1e-12)
            # A quotient-graph orbit is modulo lattice translations, not modulo
            # the full space group. Rotated wrapped positions remain distinct.
            wrapped_key = tuple(round(float(x), POSITION_DEDUP_DECIMALS) % 1.0 for x in wrapped)
            orbit = "source-position:" + str(atom["source_index0"]) + ":" + canonical_json(wrapped_key)
            ranges = [range(math.ceil(center[k] - wrapped[k] - bound[k] - 1e-10),
                            math.floor(center[k] - wrapped[k] + bound[k] + 1e-10) + 1) for k in range(3)]
            for shift in product(*ranges):
                image = wrapped + np.asarray(shift, dtype=float)
                vector = snapshot.matrix @ (image - center)
                distance = float(np.linalg.norm(vector))
                if distance > radius + 1e-10 or distance <= 1e-10:
                    continue
                key = (atom["source_index0"], *(round(float(x), POSITION_DEDUP_DECIMALS) for x in image))
                row = {**atom, "wrapped_fractional": wrapped.tolist(),
                       "fractional_image": image.tolist(), "translation": list(shift),
                       "relative_translation": (np.asarray(shift,dtype=int)-np.floor(center+1e-12).astype(int)).tolist(),
                       "vector_A": vector.tolist(), "distance_A": distance,
                       "symmetry_operation": op.triplet(), "orbit": orbit}
                if key not in found or row["symmetry_operation"] < found[key]["symmetry_operation"]:
                    found[key] = row
    return sorted(found.values(), key=lambda r: (r["element"], tuple(np.round(r["vector_A"], 10)), r["label"]))


def group_kind(atom):
    group = atom["group"]
    if group in (None, "", ".", "?"):
        return "unassigned_partial" if atom["occupancy"] is None or atom["occupancy"] < 0.999999 else "common"
    return "special_position" if str(group).startswith("-") else "positive"


def local_conditions(snapshot, target, images):
    # Only points within a candidate's distance can block its Gabriel ball.
    potential = [p for p in images if p["element"] in DONOR_ELEMENTS and p["distance_A"] <
        COVALENT_RADIUS_FACTOR * (covalent_radius(target["element"]) + covalent_radius(p["element"]))]
    max_distance = max((p["distance_A"] for p in potential), default=0.0)
    relevant = [p for p in images if p["distance_A"] <= max_distance + 1e-10]
    choices = {}
    for atom in relevant + [target]:
        if group_kind(atom) == "positive":
            choices.setdefault(atom["assembly"], set()).add(str(atom["group"]))
    if group_kind(target) == "positive":
        choices[target["assembly"]] = {str(target["group"])}
    keys = sorted(choices)
    count = math.prod(len(choices[k]) for k in keys) if keys else 1
    if count > MAX_LOCAL_CONDITIONS:
        raise ValueError("Local condition product exceeds fixed budget; do not silently prune")
    conditions = [dict(zip(keys, values)) for values in product(*(sorted(choices[k]) for k in keys))] if keys else [{}]
    return conditions, max_distance


def select_context(target, images, condition, max_distance):
    selected, unresolved = [], []
    if group_kind(target) == "special_position":
        unresolved.append({"label":target["label"],"assembly":target["assembly"],
                           "group":target["group"],"occupancy":target["occupancy"],
                           "reason":"target_special_position_disorder"})
    for atom in images:
        kind = group_kind(atom)
        if kind == "positive":
            if atom["assembly"] in condition and str(atom["group"]) != condition[atom["assembly"]]:
                continue
        if kind in ("special_position", "unassigned_partial") and atom["distance_A"] <= max_distance + 1e-10:
            unresolved.append({"label": atom["label"], "assembly": atom["assembly"],
                               "group": atom["group"], "occupancy": atom["occupancy"], "reason": kind})
        selected.append(atom)
    # A partly occupied target is conditioned as present; it is not its own blocker.
    return selected, unresolved


def geometric_contacts(target, context):
    accepted, trace = [], []
    for i, atom in enumerate(context):
        if atom["element"] not in DONOR_ELEMENTS:
            continue
        cap = COVALENT_RADIUS_FACTOR * (covalent_radius(target["element"]) + covalent_radius(atom["element"]))
        d = atom["distance_A"]
        if d > cap + RADIUS_BOUNDARY_TOLERANCE_A:
            continue
        row = {"source_label": atom["label"], "element": atom["element"], "vector_A": atom["vector_A"],
               "distance_A": d, "radius_cap_A": cap, "blockers": [], "accepted": False}
        if abs(d - cap) <= RADIUS_BOUNDARY_TOLERANCE_A:
            row["reason"] = "NUMERICAL_RADIUS_BOUNDARY_UNRESOLVED"
        else:
            v = np.asarray(atom["vector_A"])
            for j, blocker in enumerate(context):
                if i == j:
                    continue
                w = np.asarray(blocker["vector_A"])
                margin = float(w @ w + (w-v) @ (w-v) - d*d)
                if margin < -GABRIEL_MARGIN_TOLERANCE_A2:
                    row["blockers"].append({"source_label": blocker["label"], "element": blocker["element"],
                                            "vector_A": blocker["vector_A"], "margin_A2": margin})
            row["accepted"] = not row["blockers"]
            row["reason"] = "EMPTY_BALL_VISIBLE_CANDIDATE" if row["accepted"] else "GEOMETRICALLY_OCCLUDED_CANDIDATE_NOT_PROVEN_NONBOND"
        trace.append(row)
        if row["accepted"]:
            accepted.append(atom)
    return accepted, trace


def make_record(snapshot, target, donors, condition, unresolved):
    donors = sorted(donors, key=lambda x: (x["element"], tuple(np.round(x["vector_A"], 10)), x["label"]))
    groups, incidences, pairs = [], [], []
    for i, atom in enumerate(donors):
        label = "G%04d" % (i+1)
        groups.append(DonorGroup(label, "atom", (atom["element"],),
            (digest({"element": atom["element"], "provenance": "coordinate_inferred_atom"}),),
            ("M1",), periodic_orbit_key=atom["orbit"],
            evidence="CIF_coordinate_radius_Gabriel", bond_types=("coordinate_inferred",)))
        incidences.append(CoordinationIncidence("M1", label, 0, atom["element"], "sigma_candidate",
            "coordinate_inferred", tuple(int(x) for x in atom["relative_translation"]), True, round(atom["distance_A"], 4)))
    vectors = np.asarray([d["vector_A"] for d in donors], dtype=float).reshape((-1,3))
    radial = tuple(sorted(round(d["distance_A"],4) for d in donors))
    angular = []
    for i,j in combinations(range(len(donors)),2):
        a,b=vectors[i],vectors[j]
        cosine=float(np.clip(a@b/(np.linalg.norm(a)*np.linalg.norm(b)),-1,1))
        angle=math.degrees(math.acos(cosine)); angular.append(round(cosine,6))
        kind="trans" if angle>=150 else "cis" if 60<=angle<=120 else "other"
        pairs.append(DonorPairRelation("M1", groups[i].label,groups[j].label,round(angle,4),kind,
                     round(float(np.linalg.norm(a-b)),4),"CIF_coordinate_geometry"))
    issues=["COORDINATE_INFERRED_CONTACTS", "LOCAL_CONDITION_NOT_GLOBAL_DISORDER_ENSEMBLE"]
    if unresolved: issues.append("UNRESOLVED_LOCAL_DISORDER_CONTEXT")
    shape_best=""; shape_values=(); geometry_kind="distance-angle-profile"
    if 2<=len(donors)<=6:
        try:
            shape=ShapeCalculator(rounding_decimals=4).compute(vectors,np.zeros(3))
            shape_values=tuple(float(x) for x in shape.values_rounded)
            if shape_values:
                shape_best=str(shape.ref_shapes[int(np.argmin(shape_values))]); geometry_kind="CShM"
        except Exception as exc:
            issues.append("CSHM_UNAVAILABLE:"+type(exc).__name__)
    site=MetalSiteRecord("M1",target["element"],tuple(g.label for g in groups),len(groups),len(donors),len(donors),
        0,0,0,tuple(sorted(Counter(d["element"] for d in donors).items())),geometry_kind,
        shape_best,shape_values,radial,tuple(sorted(angular)),"local_state" if donors else "audit_only",
        "ambiguous" if unresolved else "partial", "source_order_noncanonical",tuple(issues))
    record=EntryRecord(0,"", "conditional_periodic_local",("coordinate_input",),1,(site,),tuple(groups),
        tuple(incidences),(),(),(),(),tuple(pairs),"EMITTED_AMBIGUOUS" if unresolved else "EMITTED_PARTIAL",
        tuple(issues),source_database="CIF",protocol_id=VERSION).with_digest()
    errors=validate_entry_record(record)
    if errors: raise RuntimeError("Typed record validation failed: "+",".join(errors))
    keys=None
    if 2<=len(donors)<=6:
        keys=dict(canonicalize_cartesian_site(metal=target["element"], donor_elements=[d["element"] for d in donors],
            donor_vectors=vectors,donor_orbits=[d["orbit"] for d in donors],
            bond_roles=["coordinate_inferred_atom"]*len(donors)).as_dict())
        # Returned permutation indices are audit positions, not invariant identity.
    return {"record":record.to_dict(),"schema_validation_errors":errors,"local_metric_stereo":keys,
            "donor_images":donors,"condition":condition,"condition_probability":None,
            "target_presence_condition":{"source_label":target["label"],"assumed_present_for_local_record":True,
                                         "source_occupancy":target["occupancy"]},
            "unresolved_local_context":unresolved,"global_disorder_model_resolved":False,
            "input_edge_evidence":"COORDINATE_INFERRED_NOT_AUTHOR_GEOMETRY_TABLE"}


def encode_snapshot(snapshot):
    targets=[a for a in snapshot.atoms if a["element"] in ALL_METAL_SET and a["occupancy"]!=0]
    outcomes=[]
    supported=[a["element"] for a in snapshot.atoms if a["element"] in DONOR_ELEMENTS]
    for target in targets:
        radius=max((COVALENT_RADIUS_FACTOR*(covalent_radius(target["element"])+covalent_radius(e)) for e in supported),default=0.0)
        images=expanded_local_images(snapshot,target,radius)
        conditions,max_distance=local_conditions(snapshot,target,images)
        for condition in conditions:
            context,unresolved=select_context(target,images,condition,max_distance)
            donors,trace=geometric_contacts(target,context)
            result=make_record(snapshot,target,donors,condition,unresolved)
            result.update({"target":target,"block":snapshot.block,"source_issues":snapshot.issues,
                           "candidate_trace":trace,"expanded_local_heavy_image_count":len(images),
                           "selected_local_context_count":len(context),"selected_local_context":context,
                           "search_radius_A":radius})
            outcomes.append(result)
    return {"adapter_version":VERSION,"original_ASU_metal_positions":len(targets),"local_condition_records":len(outcomes),
            "source_atom_count":len(snapshot.atoms),"source_block":snapshot.block,"outcomes":outcomes,
            "source_cell_matrix_A":snapshot.matrix.tolist(),
            "reference_tables_read":False,"global_disorder_ensemble_claimed":False}


def encode_cif(path,block_name=None):
    return encode_snapshot(read_snapshot(path,block_name))
