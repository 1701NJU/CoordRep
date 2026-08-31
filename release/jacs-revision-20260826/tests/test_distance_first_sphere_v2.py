"""Chemical boundary tests for the experimental distance-completion layer."""

from coordrep.audit.distance_first_sphere_v2 import (
    COVALENT_RADII_ANGSTROM,
    DEFAULT_POLICY,
    ContactObservation,
    GeometryAtom,
    NativeEdgeKey,
    classify_distance_first_sphere,
    enumerate_molecular_observations,
    enumerate_periodic_observations,
    summarize_distance_results,
)


def obs(
    metal_key,
    metal,
    block,
    donor_key,
    donor,
    distance,
    **kwargs,
):
    return ContactObservation(
        metal_key=metal_key,
        metal_element=metal,
        metal_block=block,
        donor_key=donor_key,
        donor_element=donor,
        distance_A=distance,
        **kwargs,
    )


def test_native_edges_are_separate_and_immutable():
    result = classify_distance_first_sphere(
        [
            obs("Fe1", "Fe", "d", "O_native", "O", 2.02, native_edge=True),
            obs("Fe1", "Fe", "d", "N_added", "N", 2.05),
        ]
    )
    assert len(result.native_edges) == 1
    assert result.native_edges[0].status == "native"
    assert result.native_edges[0].reason == "NATIVE_EDGE_IMMUTABLE"
    assert len(result.accepted) == 1
    assert result.accepted[0].donor_key == "N_added"


def test_strong_distance_layer_covers_s_p_d_and_f_centres():
    result = classify_distance_first_sphere(
        [
            obs("Li1", "Li", "s", "O1", "O", 2.00),
            obs("Al1", "Al", "p", "O2", "O", 1.90),
            obs("Ni1", "Ni", "d", "N1", "N", 2.00),
            obs("U1", "U", "f", "O3", "O", 2.40),
        ]
    )
    assert {(item.metal_block, item.status) for item in result.inferred_contacts} == {
        ("s", "accepted"),
        ("p", "accepted"),
        ("d", "accepted"),
        ("f", "accepted"),
    }


def test_long_ionic_contact_requires_first_shell_gap():
    with_gap = classify_distance_first_sphere(
        [
            obs("Li1", "Li", "s", "O_first", "O", 2.45),
            obs("Li1", "Li", "s", "O_second", "O", 3.25),
        ]
    )
    first = next(item for item in with_gap.inferred_contacts if item.donor_key == "O_first")
    assert first.status == "accepted"
    assert first.reason == "IONIC_EXTENSION_WITH_FIRST_SHELL_GAP"

    without_gap = classify_distance_first_sphere(
        [obs("Li1", "Li", "s", "O_first", "O", 2.45)]
    )
    assert without_gap.inferred_contacts[0].status == "candidate"
    assert without_gap.inferred_contacts[0].reason == "WITHIN_CANDIDATE_BAND_ONLY"


def test_possible_haptic_contact_is_grouped_but_eta_is_not_inferred():
    observations = [
        obs(
            "Fe1",
            "Fe",
            "d",
            f"C{i}",
            "C",
            2.03 + i * 0.01,
            pi_component_key="ring-A",
            ligand_component_key="Cp-A",
        )
        for i in range(5)
    ]
    result = classify_distance_first_sphere(observations)
    assert len(result.collective_site_candidates) == 1
    site = result.collective_site_candidates[0]
    assert site.hapticity is None
    assert len(site.member_keys) == 5
    assert {item.status for item in result.inferred_contacts} == {"candidate"}
    assert {item.reason for item in result.inferred_contacts} == {
        "HAPTIC_ASSIGNMENT_WITHHELD"
    }


def test_periodic_images_of_one_wrapped_orbit_remain_distinct():
    result = classify_distance_first_sphere(
        [
            obs(
                "Cu1", "Cu", "d", "Cl_a", "Cl", 2.30,
                periodic=True, donor_orbit_key="Cl_orbit", image_delta=(0, 0, 0),
            ),
            obs(
                "Cu1", "Cu", "d", "Cl_b", "Cl", 2.31,
                periodic=True, donor_orbit_key="Cl_orbit", image_delta=(0, 0, 1),
            ),
        ]
    )
    assert len(result.accepted) == 2
    assert {item.image_delta for item in result.accepted} == {(0, 0, 0), (0, 0, 1)}


def test_unresolved_periodicity_and_disorder_never_auto_accept():
    result = classify_distance_first_sphere(
        [
            obs(
                "Ba1", "Ba", "s", "O_unknown_image", "O", 2.75,
                periodic=True, image_delta=None,
            ),
            obs(
                "Co1", "Co", "d", "N_alt", "N", 2.00,
                disorder_state="incompatible",
            ),
            obs(
                "Li1", "Li", "s", "O_half", "O", 2.00,
                disorder_state="compatible", donor_occupancy=0.5,
            ),
        ]
    )
    decisions = {item.donor_key: item for item in result.inferred_contacts}
    assert decisions["O_unknown_image"].status == "candidate"
    assert decisions["O_unknown_image"].reason == "PBC_IMAGE_UNRESOLVED"
    assert decisions["N_alt"].status == "rejected"
    assert decisions["N_alt"].reason == "DISORDER_ALTERNATIVES_INCOMPATIBLE"
    assert decisions["O_half"].status == "candidate"
    assert decisions["O_half"].reason == "DISORDER_ADJUDICATION_REQUIRED"


def test_unsafe_donor_and_metal_metal_proximity_are_candidates_only():
    result = classify_distance_first_sphere(
        [
            obs("Ni1", "Ni", "d", "C_close", "C", 1.90),
            obs("Ni1", "Ni", "d", "Co_close", "Co", 2.40, donor_is_metal=True),
        ]
    )
    decisions = {item.donor_key: item for item in result.inferred_contacts}
    assert decisions["C_close"].status == "candidate"
    assert decisions["C_close"].reason == "DONOR_CHEMISTRY_ADJUDICATION_REQUIRED"
    assert decisions["Co_close"].status == "candidate"
    assert decisions["Co_close"].reason == "METAL_METAL_INFERENCE_WITHHELD"


def test_second_shell_and_missing_radius_are_protocol_rejections():
    result = classify_distance_first_sphere(
        [
            obs("Ni1", "Ni", "d", "O_first", "O", 2.00),
            obs("Ni1", "Ni", "d", "O_second", "O", 2.80),
            obs("Og1", "Og", "p", "O1", "O", 2.00),
        ]
    )
    decisions = {item.donor_key + ":" + item.metal_key: item for item in result.inferred_contacts}
    assert decisions["O_first:Ni1"].status == "accepted"
    assert decisions["O_second:Ni1"].status == "rejected"
    assert decisions["O1:Og1"].reason == "FROZEN_RADIUS_UNAVAILABLE"
    assert COVALENT_RADII_ANGSTROM["Og"] is None


def test_release_summary_is_stratified_and_policy_digest_is_stable():
    result = classify_distance_first_sphere(
        [
            obs("Al1", "Al", "p", "O1", "O", 1.90),
            obs("U1", "U", "f", "C1", "C", 2.50),
        ]
    )
    summary = summarize_distance_results([result])
    assert summary["distance_observations"] == 2
    assert summary["by_block_and_status"]["p:accepted"] == 1
    assert summary["by_block_and_status"]["f:candidate"] == 1
    assert result.policy_digest == DEFAULT_POLICY.digest()
    assert len(result.policy_digest) == 64


def test_molecular_enumerator_retains_native_edges_beyond_distance_decisions():
    atoms = [
        GeometryAtom("Al1", "Al", (0.0, 0.0, 0.0), is_metal=True, metal_block="p"),
        GeometryAtom("O1", "O", (1.9, 0.0, 0.0)),
        GeometryAtom("N1", "N", (0.0, 2.0, 0.0)),
    ]
    observations = enumerate_molecular_observations(
        atoms, native_edges=[NativeEdgeKey("Al1", "O1")]
    )
    result = classify_distance_first_sphere(observations)
    assert [item.donor_key for item in result.native_edges] == ["O1"]
    assert [item.donor_key for item in result.accepted] == ["N1"]


def test_periodic_enumerator_emits_two_images_of_one_orbit():
    atoms = [
        GeometryAtom(
            "Cu1", "Cu", (0.0, 0.0, 0.0), is_metal=True, metal_block="d",
            orbit_key="Cu_orbit",
        ),
        GeometryAtom("Cl1", "Cl", (0.5, 0.0, 0.0), orbit_key="Cl_orbit"),
    ]
    observations = enumerate_periodic_observations(
        atoms,
        cell_vectors_A=((4.0, 0.0, 0.0), (0.0, 10.0, 0.0), (0.0, 0.0, 10.0)),
        native_edges=[NativeEdgeKey("Cu1", "Cl1", (0, 0, 0))],
    )
    selected = [
        item for item in observations
        if item.metal_key == "Cu1" and item.donor_key == "Cl1" and item.distance_A <= 2.01
    ]
    assert {item.image_delta for item in selected} == {(-1, 0, 0), (0, 0, 0)}
    result = classify_distance_first_sphere(selected)
    assert len(result.native_edges) == 1
    assert len(result.accepted) == 1
    assert result.accepted[0].image_delta == (-1, 0, 0)
