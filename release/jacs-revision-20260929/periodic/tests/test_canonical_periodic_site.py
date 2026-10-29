"""Regression tests for the CR-PLS/2 periodic local-state canonicalizer."""

from __future__ import annotations

import itertools
import math
import random
import sys
import unittest
from pathlib import Path
from typing import Hashable, Iterable, List, Sequence, Tuple

import numpy as np


SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR.parent))

from periodic_v3 import (  # noqa: E402
    PeriodicDonorImage,
    canonical_entry_multiset_id,
    canonicalize_cartesian_site,
    canonicalize_fractional_site,
)


def _split_fractional_image(
    value: Sequence[float],
) -> Tuple[Tuple[float, float, float], Tuple[int, int, int]]:
    image = np.asarray(value, dtype=float)
    translation = np.floor(image).astype(int)
    base = image - translation
    return tuple(float(item) for item in base), tuple(
        int(item) for item in translation
    )


def _make_donors(
    totals: Sequence[Sequence[float]],
    elements: Sequence[str],
    orbits: Sequence[Hashable],
) -> List[PeriodicDonorImage]:
    result = []
    for total, element, orbit in zip(totals, elements, orbits):
        base, translation = _split_fractional_image(total)
        result.append(
            PeriodicDonorImage(
                element=element,
                fractional_base=base,
                translation=translation,
                orbit=orbit,
            )
        )
    return result


class PeriodicCanonicalUpgradeTests(unittest.TestCase):
    def setUp(self) -> None:
        self.cell = np.array(
            [[9.2, 0.0, 0.0], [1.1, 10.3, 0.0], [0.7, 1.5, 11.7]],
            dtype=float,
        )
        self.metal_total = np.array([0.75, 0.75, -0.3485], dtype=float)
        self.cartesian_vectors = np.array(
            [
                [2.062149, 0.0, 0.0],
                [0.0, 2.062149, 0.0],
                [-2.062149, 0.0, 0.0],
                [0.0, -2.062149, 0.0],
                [0.0, 0.0, -2.695543],
                [0.0, 0.0, 2.767657],
            ],
            dtype=float,
        )
        self.elements = ["N", "N", "N", "N", "Cl", "Cl"]
        # The two axial chlorides are translated images of one quotient orbit.
        self.orbits = ["n1", "n2", "n3", "n4", "cl", "cl"]
        relative_fractional = self.cartesian_vectors @ np.linalg.inv(self.cell)
        self.donor_totals = self.metal_total[None, :] + relative_fractional
        self.metal_base, self.metal_translation = _split_fractional_image(
            self.metal_total
        )
        self.donors = _make_donors(
            self.donor_totals, self.elements, self.orbits
        )
        self.reference = canonicalize_fractional_site(
            metal="Cu",
            metal_fractional_base=self.metal_base,
            metal_translation=self.metal_translation,
            donors=self.donors,
            cell_matrix=self.cell,
        )

    def assert_same(self, result) -> None:
        self.assertEqual(result.metric_id, self.reference.metric_id)
        self.assertEqual(result.stereo_id, self.reference.stereo_id)
        self.assertEqual(result.token, self.reference.token)

    def test_all_input_orders_for_cn4_and_random_orders_for_cn6(self) -> None:
        cn4_reference = canonicalize_cartesian_site(
            metal="Pt",
            donor_elements=["N", "Cl", "P", "S"],
            donor_vectors=self.cartesian_vectors[:4],
            donor_orbits=[10, 20, 30, 40],
        )
        for permutation in itertools.permutations(range(4)):
            result = canonicalize_cartesian_site(
                metal="Pt",
                donor_elements=[
                    ["N", "Cl", "P", "S"][index]
                    for index in permutation
                ],
                donor_vectors=self.cartesian_vectors[list(permutation[:4])],
                donor_orbits=[[10, 20, 30, 40][index] for index in permutation],
            )
            self.assertEqual(result.metric_id, cn4_reference.metric_id)
            self.assertEqual(result.stereo_id, cn4_reference.stereo_id)

        rng = random.Random(20260803)
        for _ in range(50):
            permutation = list(range(6))
            rng.shuffle(permutation)
            self.assert_same(
                canonicalize_fractional_site(
                    metal="Cu",
                    metal_fractional_base=self.metal_base,
                    metal_translation=self.metal_translation,
                    donors=[self.donors[index] for index in permutation],
                    cell_matrix=self.cell,
                )
            )

    def test_global_origin_and_common_image_shifts(self) -> None:
        rng = np.random.default_rng(20260803)
        for _ in range(30):
            origin = rng.uniform(-2.0, 2.0, size=3)
            common_image = rng.integers(-4, 5, size=3)
            metal_total = self.metal_total + origin + common_image
            donor_totals = self.donor_totals + origin + common_image
            metal_base, metal_translation = _split_fractional_image(metal_total)
            donors = _make_donors(donor_totals, self.elements, self.orbits)
            self.assert_same(
                canonicalize_fractional_site(
                    metal="Cu",
                    metal_fractional_base=metal_base,
                    metal_translation=metal_translation,
                    donors=donors,
                    cell_matrix=self.cell,
                )
            )

    def test_independent_representative_gauge_and_orbit_renaming(self) -> None:
        rng = np.random.default_rng(20260804)
        for trial in range(30):
            metal_gauge = rng.integers(-5, 6, size=3)
            metal_base = np.asarray(self.metal_base) + metal_gauge
            metal_translation = np.asarray(self.metal_translation) - metal_gauge
            orbit_rename = {
                value: f"trial{trial}_orbit{index * 17 + 3}"
                for index, value in enumerate(dict.fromkeys(self.orbits))
            }
            donors = []
            for donor in self.donors:
                gauge = rng.integers(-5, 6, size=3)
                donors.append(
                    PeriodicDonorImage(
                        element=donor.element,
                        fractional_base=tuple(
                            np.asarray(donor.fractional_base) + gauge
                        ),
                        translation=tuple(
                            int(value)
                            for value in np.asarray(donor.translation) - gauge
                        ),
                        orbit=orbit_rename[donor.orbit],
                    )
                )
            self.assert_same(
                canonicalize_fractional_site(
                    metal="Cu",
                    metal_fractional_base=metal_base,
                    metal_translation=metal_translation,
                    donors=donors,
                    cell_matrix=self.cell,
                )
            )

    def test_unimodular_cell_basis_transformations(self) -> None:
        transforms = [
            np.eye(3, dtype=int),
            np.array([[1, 1, 0], [0, 1, 0], [0, 0, 1]], dtype=int),
            np.array([[1, 0, 0], [0, 1, 1], [0, 0, 1]], dtype=int),
            np.array([[0, 1, 0], [1, 0, 0], [0, 0, 1]], dtype=int),
            np.array([[0, 1, 0], [-1, 0, 0], [0, 0, 1]], dtype=int),
            np.array([[-1, 0, 0], [0, -1, 0], [0, 0, 1]], dtype=int),
            np.array([[1, 1, 0], [1, 0, 0], [0, 0, 1]], dtype=int),
            np.array([[1, 0, 1], [0, 1, 0], [0, 0, 1]], dtype=int),
            np.array([[1, 0, 0], [1, 1, 0], [0, 1, 1]], dtype=int),
            np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]], dtype=int),
        ]
        for transform in transforms:
            with self.subTest(transform=transform.tolist()):
                determinant = int(round(float(np.linalg.det(transform))))
                self.assertIn(determinant, {-1, 1})
                inverse = np.linalg.inv(transform)
                transformed_cell = transform @ self.cell
                transformed_metal = self.metal_total @ inverse
                transformed_donors = self.donor_totals @ inverse
                metal_base, metal_translation = _split_fractional_image(
                    transformed_metal
                )
                donors = _make_donors(
                    transformed_donors,
                    self.elements,
                    [101, 203, 307, 409, 503, 503],
                )
                self.assert_same(
                    canonicalize_fractional_site(
                        metal="Cu",
                        metal_fractional_base=metal_base,
                        metal_translation=metal_translation,
                        donors=donors,
                        cell_matrix=transformed_cell,
                    )
                )

    def test_proper_rotation_and_mirror_layer_boundary(self) -> None:
        vectors = np.array(
            [
                [1.1, 0.2, 0.3],
                [-0.4, 1.3, 0.1],
                [0.2, -0.5, 1.4],
                [-1.0, -0.7, -0.9],
            ]
        )
        elements = ["N", "O", "P", "S"]
        orbits = ["a", "b", "c", "d"]
        reference = canonicalize_cartesian_site(
            metal="Co",
            donor_elements=elements,
            donor_vectors=vectors,
            donor_orbits=orbits,
        )
        angle = 0.731
        rotation = np.array(
            [
                [math.cos(angle), -math.sin(angle), 0.0],
                [math.sin(angle), math.cos(angle), 0.0],
                [0.0, 0.0, 1.0],
            ]
        )
        rotated = canonicalize_cartesian_site(
            metal="Co",
            donor_elements=elements,
            donor_vectors=vectors @ rotation,
            donor_orbits=orbits,
        )
        self.assertEqual(rotated.metric_id, reference.metric_id)
        self.assertEqual(rotated.stereo_id, reference.stereo_id)

        reflection = np.diag([-1.0, 1.0, 1.0])
        mirrored = canonicalize_cartesian_site(
            metal="Co",
            donor_elements=elements,
            donor_vectors=vectors @ reflection,
            donor_orbits=orbits,
        )
        self.assertEqual(mirrored.metric_id, reference.metric_id)
        self.assertNotEqual(mirrored.stereo_id, reference.stereo_id)

    def test_periodic_orbit_partition_is_chemical_information(self) -> None:
        shared_orbit = canonicalize_cartesian_site(
            metal="Cu",
            donor_elements=self.elements,
            donor_vectors=self.cartesian_vectors,
            donor_orbits=self.orbits,
        )
        distinct_orbits = canonicalize_cartesian_site(
            metal="Cu",
            donor_elements=self.elements,
            donor_vectors=self.cartesian_vectors,
            donor_orbits=["n1", "n2", "n3", "n4", "cl1", "cl2"],
        )
        self.assertNotEqual(shared_orbit.metric_id, distinct_orbits.metric_id)
        self.assertNotEqual(shared_orbit.stereo_id, distinct_orbits.stereo_id)

    def test_a_genuinely_different_periodic_image_is_not_merged(self) -> None:
        altered = list(self.donors)
        donor = altered[-1]
        shifted_translation = list(donor.translation)
        shifted_translation[2] += 1
        altered[-1] = PeriodicDonorImage(
            element=donor.element,
            fractional_base=donor.fractional_base,
            translation=tuple(shifted_translation),
            orbit=donor.orbit,
        )
        changed = canonicalize_fractional_site(
            metal="Cu",
            metal_fractional_base=self.metal_base,
            metal_translation=self.metal_translation,
            donors=altered,
            cell_matrix=self.cell,
        )
        self.assertNotEqual(changed.metric_id, self.reference.metric_id)

    def test_entry_multiset_preserves_multiplicity_not_order(self) -> None:
        values = ["b", "a", "b", "c"]
        self.assertEqual(
            canonical_entry_multiset_id(values),
            canonical_entry_multiset_id(reversed(values)),
        )
        self.assertNotEqual(
            canonical_entry_multiset_id(values),
            canonical_entry_multiset_id(["a", "b", "c"]),
        )


if __name__ == "__main__":
    unittest.main(verbosity=2)
