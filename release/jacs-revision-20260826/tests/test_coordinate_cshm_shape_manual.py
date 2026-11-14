"""Independent numerical anchors for the coordinate-space CShM repair."""

from itertools import permutations

import numpy as np
import pytest

from coordrep.geometry.shape import (
    IDEAL_GEOMETRIES,
    ShapeCalculator,
    compute_coordinate_cshm,
)


# Worked examples printed in the SHAPE 2.1 manual.
ABOXIY = np.array([
    [2.7159, 8.9642, 15.0153],
    [3.9023, 7.5659, 14.2563],
    [3.9912, 8.8145, 16.4883],
    [1.0864, 9.0325, 13.3242],
    [1.4893, 10.6356, 16.0313],
])

ACACPT = np.array([
    [0.0, 0.0, 0.0],
    [0.6294, 1.3760, -1.2703],
    [-0.6294, -1.3760, 1.2703],
    [-1.7613, 0.9024, 0.3486],
    [1.7613, -0.9024, -0.3486],
])


@pytest.mark.parametrize(
    "structure,expected_tetrahedral,expected_square_planar",
    [
        (ABOXIY, 31.374675123038607, 0.9695734707444207),
        (ACACPT, 33.43969489575136, 0.15954234362706737),
    ],
)
def test_shape_2_1_manual_cn4_examples(
    structure,
    expected_tetrahedral,
    expected_square_planar,
):
    metal = structure[0]
    donors = structure[1:]
    tetrahedral, _ = compute_coordinate_cshm(
        donors, metal, IDEAL_GEOMETRIES[4]["Td"]
    )
    square_planar, _ = compute_coordinate_cshm(
        donors, metal, IDEAL_GEOMETRIES[4]["SP"]
    )

    assert tetrahedral == pytest.approx(expected_tetrahedral, abs=5e-9)
    assert square_planar == pytest.approx(expected_square_planar, abs=5e-9)


@pytest.mark.parametrize("shape_name", ["Td", "SP"])
def test_all_24_donor_orders_recover_the_ideal_cn4_shape(shape_name):
    metal = np.array([3.1, -2.7, 5.4])
    template = IDEAL_GEOMETRIES[4][shape_name]

    # A nontrivial proper rotation, scale, and translation ensure this is not
    # merely a comparison in the template's native coordinate frame.
    axis = np.array([1.0, 2.0, -0.5])
    axis /= np.linalg.norm(axis)
    angle = 0.73
    cross = np.array([
        [0.0, -axis[2], axis[1]],
        [axis[2], 0.0, -axis[0]],
        [-axis[1], axis[0], 0.0],
    ])
    rotation = (
        np.eye(3) * np.cos(angle)
        + (1.0 - np.cos(angle)) * np.outer(axis, axis)
        + np.sin(angle) * cross
    )
    donors = metal + 2.35 * (template @ rotation)

    calculator = ShapeCalculator(rounding_decimals=8)
    for order in permutations(range(4)):
        result = calculator.compute(donors[list(order)], metal)
        assert result.best_shape == shape_name
        assert result.best_cshm == pytest.approx(0.0, abs=1e-20)


@pytest.mark.parametrize(
    "coordination_number,shape_name",
    [
        (2, "L"),
        (3, "TP"),
        (4, "Td"),
        (4, "SP"),
        (5, "TBP"),
        (5, "SPY"),
        (6, "Oh"),
        (6, "TPr"),
    ],
)
def test_each_cosym_reference_is_a_zero_of_its_own_score(
    coordination_number,
    shape_name,
):
    template = IDEAL_GEOMETRIES[coordination_number][shape_name]
    score, _ = compute_coordinate_cshm(
        donor_coords=template,
        metal_coord=np.zeros(3),
        template_donors=template,
    )
    assert score == pytest.approx(0.0, abs=1e-20)


def test_nonideal_shape_is_invariant_to_translation_rotation_scale_and_order():
    baseline = ShapeCalculator(rounding_decimals=12).compute(
        ABOXIY[1:], ABOXIY[0]
    )

    axis = np.array([-0.4, 0.7, 1.2])
    axis /= np.linalg.norm(axis)
    angle = 1.17
    cross = np.array([
        [0.0, -axis[2], axis[1]],
        [axis[2], 0.0, -axis[0]],
        [-axis[1], axis[0], 0.0],
    ])
    rotation = (
        np.eye(3) * np.cos(angle)
        + (1.0 - np.cos(angle)) * np.outer(axis, axis)
        + np.sin(angle) * cross
    )
    translation = np.array([11.0, -4.5, 8.25])
    transformed = 3.7 * (ABOXIY @ rotation) + translation
    donor_order = [2, 0, 3, 1]
    transformed_result = ShapeCalculator(rounding_decimals=12).compute(
        transformed[1:][donor_order], transformed[0]
    )

    assert transformed_result.ref_shapes == baseline.ref_shapes
    assert transformed_result.values == pytest.approx(
        baseline.values, abs=2e-12
    )
    assert transformed_result.best_shape == baseline.best_shape
    assert transformed_result.delta_to_second == pytest.approx(
        baseline.delta_to_second, abs=2e-12
    )


def test_kabsch_guard_rejects_an_improper_reflection_for_generic_geometry():
    template = np.array([
        [1.2, 0.1, -0.2],
        [-0.4, 1.7, 0.3],
        [0.2, -0.6, 2.1],
        [-1.4, -0.8, -0.5],
    ])
    metal = np.zeros(3)
    self_score, _ = compute_coordinate_cshm(template, metal, template)

    reflection = np.diag([-1.0, 1.0, 1.0])
    reflected = template @ reflection
    reflected_score, _ = compute_coordinate_cshm(
        reflected, metal, template
    )

    assert self_score == pytest.approx(0.0, abs=1e-20)
    # An implementation that silently permits det(R)=-1 would return zero.
    assert reflected_score > 0.1


@pytest.mark.parametrize(
    "coordination_number,observed_name,target_name,expected",
    [
        (4, "Td", "SP", 33.333333333323246),
        (4, "SP", "Td", 33.33333333332325),
        (5, "TBP", "SPY", 5.383811324174554),
        (5, "SPY", "TBP", 5.383811324174555),
        (6, "Oh", "TPr", 16.736755360769074),
        (6, "TPr", "Oh", 16.736755360769077),
    ],
)
def test_official_cosym_cross_reference_values(
    coordination_number,
    observed_name,
    target_name,
    expected,
):
    observed = IDEAL_GEOMETRIES[coordination_number][observed_name]
    target = IDEAL_GEOMETRIES[coordination_number][target_name]
    score, _ = compute_coordinate_cshm(observed, np.zeros(3), target)
    assert score == pytest.approx(expected, abs=5e-10)
