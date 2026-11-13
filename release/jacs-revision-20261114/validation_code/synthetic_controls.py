"""Source-independent constructed controls; no deposited structure is bundled.

Creates only temporary, invented coordinates. A pass is a software control, not
an 84-case chemical-validation result. No CCDC installation is required.
"""
import argparse
import itertools
import json
from pathlib import Path
import tempfile
import numpy as np

from run_input import HERE, encode_xyz, encode_cif, write_new
from compare_local_sets import compare_images
from frozen_canonical_periodic_site import canonicalize_cartesian_site


def xyz_text(elements, positions):
    return str(len(elements)) + "\nConstructed software control, not a compound\n" + "\n".join(
        e + " " + " ".join("%.12f" % float(c) for c in v)
        for e, v in zip(elements, positions)) + "\n"


def cif_text(rows):
    return """data_constructed
_cell_length_a 10
_cell_length_b 10
_cell_length_c 10
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_space_group_name_H-M_alt 'P 1'
loop_
_space_group_symop_operation_xyz
'x,y,z'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
_atom_site_disorder_assembly
_atom_site_disorder_group
""" + "\n".join(rows) + "\n"


def run_controls():
    results = []
    def check(name, condition):
        results.append({"control": name, "passed": bool(condition)})
        if not condition: raise AssertionError(name)

    with tempfile.TemporaryDirectory(prefix="coordrep_constructed_") as d:
        root = Path(d)
        elements = ["Fe", "N", "O", "S", "Cl"]
        xyz = np.array([[0., 0., 0.], [2., 0., 0.], [0., 2.1, 0.],
                        [0., 0., 2.2], [-1.4, -1.5, -1.6]])
        first = root / "sphere.xyz"
        first.write_text(xyz_text(elements, xyz), encoding="ascii")
        encoded = encode_xyz(first, ["Fe"], "constructed")
        site = encoded["sites"][0]
        check("XYZ four-donor sphere encoded", site.get("atom_contact_cn") == 4 and
              site.get("canonicalization_status") == "COMPLETE")
        key = site["coordrep_string"]
        angle = .371
        rotation = np.array([[np.cos(angle), -np.sin(angle), 0.],
                             [np.sin(angle), np.cos(angle), 0.], [0., 0., 1.]])
        permutation = [4, 2, 0, 3, 1]
        second = root / "reexpressed.xyz"
        changed = (xyz @ rotation.T + np.array([3.4, -2.1, .7]))[permutation]
        second.write_text(xyz_text([elements[i] for i in permutation], changed), encoding="ascii")
        other = encode_xyz(second, ["Fe"], "constructed")["sites"][0]
        check("XYZ rotation translation and atom-order invariance", other.get("coordrep_string") == key)

        # Direct neighbour tests isolate the fixed admission/screening policy.
        from coordrep.encode import _load_molecule
        from coordrep.graph.neighbors import perceive_neighbors
        silicon = root / "silicon.xyz"
        silicon.write_text(xyz_text(["Fe", "Si"], [[0, 0, 0], [2.3, 0, 0]]), encoding="ascii")
        donors, _ = perceive_neighbors(_load_molecule(silicon), 0)
        check("Si is an eligible donor", [x.atom_idx for x in donors] == [1])
        blocked = root / "blocked.xyz"
        blocked.write_text(xyz_text(["Fe", "C", "B"], [[0,0,0], [2,0,0], [1,0,0]]), encoding="ascii")
        donors, _ = perceive_neighbors(_load_molecule(blocked), 0)
        check("Original non-donor heavy atom screens an occluded candidate", 1 not in [x.atom_idx for x in donors])
        hydrogen = root / "hydrogen.xyz"
        hydrogen.write_text(xyz_text(["Fe", "C", "H"], [[0,0,0], [2,0,0], [1,0,0]]), encoding="ascii")
        donors, _ = perceive_neighbors(_load_molecule(hydrogen), 0)
        check("H is excluded as a donor and blocker", [x.atom_idx for x in donors] == [1])

        periodic = root / "periodic.cif"
        periodic.write_text(cif_text([
            "Na1 Na .95 .5 .5 1 . .", "F1 F .19 .5 .5 1 . .",
            "F2 F .71 .5 .5 1 . .", "F3 F .95 .74 .5 1 . .",
            "F4 F .95 .26 .5 1 . .", "F5 F .95 .5 .74 1 . .",
            "F6 F .95 .5 .26 1 . ."]), encoding="ascii")
        out = encode_cif(periodic, case_id="constructed_periodic")
        sphere = out["outcomes"][0]
        check("Periodic image reconstruction retains six donors", len(out["outcomes"]) == 1 and
              len(sphere["donor_images"]) == 6)
        check("Periodic nonzero image is retained", any(any(x["relative_translation"]) for x in sphere["donor_images"]))
        check("Periodic typed incidences retained", len(sphere["record"]["incidences"]) == 6)

        disorder = root / "disorder.cif"
        disorder.write_text(cif_text([
            "Zn1 Zn .5 .5 .5 1 . .", "N1 N .69 .5 .5 1 . .", "N2 N .31 .5 .5 1 . .",
            "Cl1 Cl .5 .72 .5 .5 A 1", "Cl2 Cl .5 .5 .72 .5 A 1",
            "Cl3 Cl .5 .28 .5 .5 A 2", "Cl4 Cl .5 .5 .28 .5 A 2"]), encoding="ascii")
        out = encode_cif(disorder, case_id="constructed_disorder")
        check("Alternative disorder conditions are not unioned", len(out["outcomes"]) == 2 and
              all(len(s["donor_images"]) == 4 for s in out["outcomes"]))
        check("Source disorder selectors preserved", {tuple(sorted(s["condition"].items())) for s in out["outcomes"]}
              == {(("A", "1"),), (("A", "2"),)})
        check("No unsupported conformer probability assigned", all(s["condition_probability"] is None for s in out["outcomes"]))

        vec = np.array([[-1.5,-1.7,-1.2], [2.,.3,.2], [-.4,2.2,.5], [-.7,-.5,2.1]])
        labels = ["Cl", "N", "O", "S"]
        def identity(v, labels):
            return canonicalize_cartesian_site(metal="Zn", donor_elements=labels, donor_vectors=v,
                                               donor_orbits=list(range(4)), bond_roles=["sigma"]*4)
        base = identity(vec, labels)
        mirror = identity(vec @ np.diag([-1.,1.,1.]), labels)
        check("Reflection preserves MID", mirror.metric_id == base.metric_id)
        check("Reflection changes SID for this labelled chiral frame", mirror.stereo_id != base.stereo_id)
        check("All 24 donor orders preserve MID and SID", all(
            (lambda k: k.metric_id == base.metric_id and k.stereo_id == base.stereo_id)(
                identity(vec[list(p)], [labels[i] for i in p]))
            for p in itertools.permutations(range(4))))

        reference = [{"element":"N", "vector_A":[2.,0.,0.], "relative_translation":[1,0,0]}]
        check("Complete donor-image set matches", compare_images(reference, reference)["complete_set_equal"])
        check("Missing donor is a retained mismatch", not compare_images(reference, [])["complete_set_equal"])
        check("Additional donor is a retained mismatch", not compare_images(reference, reference*2)["complete_set_equal"])
        near = [{"element":"N", "vector_A":[2.004,0.,0.], "relative_translation":[1,0,0]}]
        far = [{"element":"N", "vector_A":[2.006,0.,0.], "relative_translation":[1,0,0]}]
        check("Fixed 0.005 A positional tolerance", compare_images(reference, near)["complete_set_equal"] and
              not compare_images(reference, far)["complete_set_equal"])
    return {"schema":"CoordRep-constructed-software-controls/1", "status":"PASS",
            "checks_passed":len(results), "deposited_compounds_tested":0,
            "licensed_CSD_input_required":False, "checks":results}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output",type=Path)
    a=p.parse_args()
    result=run_controls()
    if a.output: write_new(a.output,result)
    print(json.dumps(result,indent=2))


if __name__=="__main__":main()
