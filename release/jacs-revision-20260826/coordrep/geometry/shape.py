"""Continuous shape measures for coordination environments.

The publication-facing score in this module is a coordinate-space continuous
shape measure (CShM).  The central atom is kept as a distinguished, fixed
vertex while every donor-to-template assignment is enumerated.  For each
assignment, the experimental and reference point sets are centered, the
reference is aligned with a proper Kabsch rotation, and an optimal uniform
scale is fitted:

    S = 100 * min ||X - s Y R||^2 / ||X||^2

This definition reproduces the worked SHAPE 2.1 examples used by the
regression tests in ``tests/test_coordinate_cshm_shape_manual.py``.

The original CoordRep donor-distance-matrix scores are retained under
explicit ``legacy_distance_matrix`` names solely for backward-compatible
forensic comparisons.  They are not CShM and are not used by
``ShapeCalculator``.
"""

from dataclasses import dataclass
from itertools import permutations
from typing import Dict, List, Tuple

import numpy as np
from scipy.spatial.distance import cdist


def _normalize_template(coords: np.ndarray) -> np.ndarray:
    """Scale donor coordinates about the fixed central atom at the origin.

    Translation is deliberately not removed here.  In particular, the donor
    centroid of a square pyramid does not coincide with its central atom.
    Coordinate-space CShM performs the required translation after the central
    atom has been included in the point set.
    """

    coords = np.asarray(coords, dtype=float)
    avg_radius = np.linalg.norm(coords, axis=1).mean()
    if avg_radius > 1e-12:
        coords = coords / avg_radius
    return coords


# Reference coordinates below follow CoSyM 0.10.9
# ``cosymlib/shape/ideal_structures_center.yaml``.  CoordRep retains its
# historical compact labels in serialized output; their CoSyM/SHAPE labels are
# L-2, TP-3, T-4, SP-4, TBPY-5, SPY-5, OC-6, and TPR-6, respectively.

# CN=2
_LINEAR_2 = _normalize_template(np.array([
    [0, 0, 1],
    [0, 0, -1],
]))

# CN=3
_TRIGONAL_PLANAR_3 = _normalize_template(np.array([
    [1, 0, 0],
    [-0.5, np.sqrt(3) / 2, 0],
    [-0.5, -np.sqrt(3) / 2, 0],
]))

# CN=4
_TETRAHEDRAL_4 = _normalize_template(np.array([
    [0, 0.912870929175, -0.645497224368],
    [0, -0.912870929175, -0.645497224368],
    [0.912870929175, 0, 0.645497224368],
    [-0.912870929175, 0, 0.645497224368],
]))

_SQUARE_PLANAR_4 = _normalize_template(np.array([
    [1.118033988750, 0, 0],
    [0, 1.118033988750, 0],
    [-1.118033988750, 0, 0],
    [0, -1.118033988750, 0],
]))

# CN=5
_TRIGONAL_BIPYRAMIDAL_5 = _normalize_template(np.array([
    [0, 0, -1.095445115010],
    [1.095445115010, 0, 0],
    [-0.547722557505, 0.948683298051, 0],
    [-0.547722557505, -0.948683298051, 0],
    [0, 0, 1.095445115010],
]))

_SQUARE_PYRAMIDAL_5 = _normalize_template(np.array([
    [0, 0, 1.095445115010],
    [1.060660171780, 0, -0.273861278753],
    [0, 1.060660171780, -0.273861278753],
    [-1.060660171780, 0, -0.273861278753],
    [0, -1.060660171780, -0.273861278753],
]))

# CN=6
_OCTAHEDRAL_6 = _normalize_template(np.array([
    [0, 0, -1.080123449735],
    [1.080123449735, 0, 0],
    [0, 1.080123449735, 0],
    [-1.080123449735, 0, 0],
    [0, -1.080123449735, 0],
    [0, 0, 1.080123449735],
]))

_TRIGONAL_PRISMATIC_6 = _normalize_template(np.array([
    [0.816496580928, 0, -0.707106781187],
    [-0.408248290464, 0.707106781187, -0.707106781187],
    [-0.408248290464, -0.707106781187, -0.707106781187],
    [0.816496580928, 0, 0.707106781187],
    [-0.408248290464, 0.707106781187, 0.707106781187],
    [-0.408248290464, -0.707106781187, 0.707106781187],
]))


IDEAL_GEOMETRIES = {
    2: {"L": _LINEAR_2},
    3: {"TP": _TRIGONAL_PLANAR_3},
    4: {"Td": _TETRAHEDRAL_4, "SP": _SQUARE_PLANAR_4},
    5: {"TBP": _TRIGONAL_BIPYRAMIDAL_5, "SPY": _SQUARE_PYRAMIDAL_5},
    6: {"Oh": _OCTAHEDRAL_6, "TPr": _TRIGONAL_PRISMATIC_6},
}


_PERMUTATION_INDEX_CACHE: Dict[int, np.ndarray] = {}
_REFERENCE_BATCH_CACHE: Dict[
    Tuple[int, bytes], Tuple[np.ndarray, float]
] = {}


def _permutation_indices(n: int) -> np.ndarray:
    """Return all lexicographically ordered permutations for ``n`` donors."""

    cached = _PERMUTATION_INDEX_CACHE.get(n)
    if cached is None:
        cached = np.asarray(list(permutations(range(n))), dtype=np.int64)
        _PERMUTATION_INDEX_CACHE[n] = cached
    return cached


def _centered_reference_batch(
    template: np.ndarray,
) -> Tuple[np.ndarray, float]:
    """Cache all centered, central-atom-fixed reference permutations."""

    key = (len(template), template.tobytes())
    cached = _REFERENCE_BATCH_CACHE.get(key)
    if cached is None:
        assignments = _permutation_indices(len(template))
        reference_donors = template[assignments]
        reference = np.concatenate(
            (
                np.zeros((len(assignments), 1, 3), dtype=float),
                reference_donors,
            ),
            axis=1,
        )
        centered = reference - reference.mean(axis=1, keepdims=True)
        denominator = float(np.einsum("ni,ni->", centered[0], centered[0]))
        cached = (centered, denominator)
        _REFERENCE_BATCH_CACHE[key] = cached
    return cached


def compute_coordinate_cshm(
    donor_coords: np.ndarray,
    metal_coord: np.ndarray,
    template_donors: np.ndarray,
) -> Tuple[float, np.ndarray]:
    """Compute the globally minimized coordinate-space CShM.

    Parameters
    ----------
    donor_coords
        Experimental donor coordinates with shape ``(CN, 3)``.
    metal_coord
        Experimental central-atom coordinate with shape ``(3,)``.
    template_donors
        Reference donor coordinates with the reference central atom fixed at
        the origin.

    Returns
    -------
    score, assignment
        The unrounded CShM and the donor-template permutation attaining the
        minimum.  The first lexicographic assignment is retained on exact
        ties, making the diagnostic deterministic.
    """

    donors = np.asarray(donor_coords, dtype=float)
    metal = np.asarray(metal_coord, dtype=float)
    template = np.asarray(template_donors, dtype=float)

    if donors.ndim != 2 or donors.shape[1:] != (3,):
        raise ValueError("donor_coords must have shape (CN, 3)")
    if metal.shape != (3,):
        raise ValueError("metal_coord must have shape (3,)")
    if template.shape != donors.shape:
        return float("inf"), np.arange(len(donors), dtype=np.int64)

    cn = len(donors)
    assignments = _permutation_indices(cn)

    # The central atom is the first and fixed point in both point sets.
    experimental = np.vstack((metal, donors))
    experimental_centered = experimental - experimental.mean(axis=0)
    experimental_denominator = float(
        np.einsum("ni,ni->", experimental_centered, experimental_centered)
    )
    if experimental_denominator < 1e-12:
        return float("inf"), np.arange(cn, dtype=np.int64)

    reference_centered, reference_denominator = _centered_reference_batch(
        template
    )

    # Closed-form batched proper Kabsch residual.  If H = U S V^T, the
    # maximum proper-rotation inner product is sum(S), except that the
    # smallest singular value is subtracted when det(UV^T) is negative.
    # det(H) has the same sign for full-rank H; for rank-deficient H the
    # smallest singular value is zero and the choice is immaterial.
    covariance = np.einsum(
        "pni,nj->pij", reference_centered, experimental_centered
    )
    singular_values = np.linalg.svd(covariance, compute_uv=False)
    maximum_inner_products = singular_values.sum(axis=1)
    improper = np.linalg.det(covariance) < 0.0
    maximum_inner_products[improper] -= (
        2.0 * singular_values[improper, -1]
    )
    residual_sums = (
        experimental_denominator
        - maximum_inner_products**2 / reference_denominator
    )
    residual_sums = np.maximum(residual_sums, 0.0)
    residual_sums[
        residual_sums <= 1e-12 * experimental_denominator
    ] = 0.0
    scores = 100.0 * residual_sums / experimental_denominator

    best_index = int(np.argmin(scores))
    return float(scores[best_index]), assignments[best_index].copy()


def compute_atom_signature(sorted_distances: np.ndarray) -> int:
    """Return the frozen deterministic donor signature used by CoordRep 1.1.1."""

    discretized = np.round(sorted_distances * 10000).astype(np.int64)
    signature = 0
    for distance in discretized[:6]:
        signature = signature * 100003 + int(distance)
    return signature % (10**15)


def compute_legacy_distance_matrix_score_greedy(
    donor_vectors: np.ndarray,
    template_donors: np.ndarray,
    atom_signatures: np.ndarray,
) -> Tuple[float, np.ndarray]:
    """Reproduce the frozen greedy donor-distance-matrix score.

    This function is retained only for forensic comparison with CoordRep
    1.1.1.  It is not a coordinate-space CShM.
    """

    donors = np.asarray(donor_vectors, dtype=float)
    template = np.asarray(template_donors, dtype=float)
    n = len(donors)
    if n != len(template):
        return float("inf"), np.arange(n, dtype=np.int64)

    donor_distances = cdist(donors, donors, metric="euclidean")
    template_distances = cdist(template, template, metric="euclidean")
    donor_features = np.asarray(
        [np.sort(donor_distances[index]) for index in range(n)]
    )
    template_features = np.asarray(
        [np.sort(template_distances[index]) for index in range(n)]
    )

    signature_order = np.argsort(atom_signatures)
    assigned = set()
    assignment = np.zeros(n, dtype=np.int64)
    for donor_index in signature_order:
        best_template_index = -1
        best_cost = float("inf")
        for template_index in range(n):
            if template_index in assigned:
                continue
            cost = float(
                np.sum(
                    (
                        donor_features[donor_index]
                        - template_features[template_index]
                    )
                    ** 2
                )
            )
            cost += 1e-10 * template_index
            if cost < best_cost:
                best_cost = cost
                best_template_index = template_index
        assignment[donor_index] = best_template_index
        assigned.add(best_template_index)

    permuted_template_distances = template_distances[assignment][:, assignment]
    numerator = float(
        np.sum((donor_distances - permuted_template_distances) ** 2)
    )
    denominator = float(np.sum(donor_distances**2))
    if denominator < 1e-12:
        return float("inf"), assignment
    return 100.0 * numerator / denominator, assignment


def compute_legacy_distance_matrix_score_exact(
    donor_vectors: np.ndarray,
    template_donors: np.ndarray,
) -> Tuple[float, np.ndarray]:
    """Globally minimize the legacy donor-distance score over assignments.

    This exact enumeration diagnoses the frozen greedy-assignment defect but
    remains a donor-distance-matrix score, not a CShM.
    """

    donors = np.asarray(donor_vectors, dtype=float)
    template = np.asarray(template_donors, dtype=float)
    n = len(donors)
    if n != len(template):
        return float("inf"), np.arange(n, dtype=np.int64)

    donor_distances = cdist(donors, donors, metric="euclidean")
    template_distances = cdist(template, template, metric="euclidean")
    denominator = float(np.sum(donor_distances**2))
    if denominator < 1e-12:
        return float("inf"), np.arange(n, dtype=np.int64)

    assignments = _permutation_indices(n)
    permuted = template_distances[assignments[:, :, None], assignments[:, None, :]]
    differences = donor_distances[None, :, :] - permuted
    numerators = np.einsum("pij,pij->p", differences, differences)
    best_index = int(np.argmin(numerators))
    score = 100.0 * float(numerators[best_index]) / denominator
    return score, assignments[best_index].copy()


def compute_cshm_with_tiebreak(
    P: np.ndarray,
    Q_template: np.ndarray,
    atom_signatures: np.ndarray,
) -> Tuple[float, np.ndarray]:
    """Deprecated compatibility wrapper for the frozen legacy score."""

    return compute_legacy_distance_matrix_score_greedy(
        P, Q_template, atom_signatures
    )


@dataclass
class ShapeResult:
    """Continuous-shape calculation result."""

    cn: int
    ref_shapes: List[str]
    values: np.ndarray
    values_rounded: np.ndarray
    assignments: Dict[str, np.ndarray]
    best_shape: str = ""
    best_cshm: float = 0.0
    delta_to_second: float = 0.0
    shape_class: str = ""

    def to_canonical_string(self) -> str:
        if not self.ref_shapes:
            return ""

        if self.delta_to_second < 0.5:
            delta_bin = 0
        elif self.delta_to_second < 2:
            delta_bin = 1
        elif self.delta_to_second < 5:
            delta_bin = 2
        else:
            delta_bin = 3
        return f"<ShapeBest:{self.best_shape}|Delta:{delta_bin}>"

    def to_dict(self) -> dict:
        return {
            "cn": self.cn,
            "ref_shapes": self.ref_shapes,
            "values_rounded": self.values_rounded.tolist(),
            "best_shape": self.best_shape,
            "best_cshm": round(self.best_cshm, 2),
            "delta_to_second": round(self.delta_to_second, 2),
        }


class ShapeCalculator:
    """Assignment-, translation-, rotation-, and scale-invariant CShM."""

    def __init__(self, rounding_decimals: int = 2):
        self.rounding_decimals = rounding_decimals

    def _compute_atom_signatures(
        self, donor_coords: np.ndarray, metal_coord: np.ndarray
    ) -> np.ndarray:
        """Compute frozen signatures for backward-compatible diagnostics."""

        n = len(donor_coords)
        signatures = np.zeros(n, dtype=np.int64)
        distances_to_metal = np.linalg.norm(
            donor_coords - metal_coord, axis=1
        )
        donor_distances = cdist(donor_coords, donor_coords, metric="euclidean")
        for index in range(n):
            feature = np.concatenate(
                ([distances_to_metal[index]], np.sort(donor_distances[index]))
            )
            signatures[index] = compute_atom_signature(feature)
        return signatures

    def _canonicalize_donor_order(
        self,
        donor_coords: np.ndarray,
        metal_coord: np.ndarray,
        signatures: np.ndarray,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Retain the frozen deterministic donor order for serialized diagnostics."""

        if len(donor_coords) <= 1:
            return donor_coords, signatures
        sorted_indices = np.argsort(signatures)
        return donor_coords[sorted_indices], signatures[sorted_indices]

    def compute(
        self, donor_coords: np.ndarray, metal_coord: np.ndarray
    ) -> ShapeResult:
        """Compute coordinate-space CShM against all references for the CN."""

        donors = np.asarray(donor_coords, dtype=float)
        metal = np.asarray(metal_coord, dtype=float)
        cn = len(donors)

        if cn not in IDEAL_GEOMETRIES:
            return ShapeResult(
                cn=cn,
                ref_shapes=[],
                values=np.array([]),
                values_rounded=np.array([]),
                assignments={},
                best_shape=f"CN{cn}",
                best_cshm=0.0,
                delta_to_second=0.0,
                shape_class="unknown",
            )

        templates = IDEAL_GEOMETRIES[cn]
        ref_shapes = list(templates.keys())
        values = []
        assignments = {}
        for shape_name in ref_shapes:
            cshm, assignment = compute_coordinate_cshm(
                donors, metal, templates[shape_name]
            )
            values.append(cshm)
            assignments[shape_name] = assignment

        values_array = np.asarray(values, dtype=float)
        values_rounded = np.round(values_array, self.rounding_decimals)
        sorted_indices = np.argsort(values_array)
        best_index = int(sorted_indices[0])
        best_shape = ref_shapes[best_index]
        best_cshm = float(values_array[best_index])

        if len(values_array) > 1:
            second_index = int(sorted_indices[1])
            delta_to_second = float(
                values_array[second_index] - values_array[best_index]
            )
        else:
            delta_to_second = float("inf")

        if best_cshm < 1.0:
            shape_class = "ideal"
        elif best_cshm < 5.0:
            shape_class = "distorted"
        else:
            shape_class = "irregular"
        if delta_to_second < 1.0:
            shape_class = "ambiguous"

        return ShapeResult(
            cn=cn,
            ref_shapes=ref_shapes,
            values=values_array,
            values_rounded=values_rounded,
            assignments=assignments,
            best_shape=best_shape,
            best_cshm=best_cshm,
            delta_to_second=delta_to_second,
            shape_class=shape_class,
        )
