"""Canonical identity for a periodic first-sphere coordination state.

The integer translation label attached to a quotient-graph edge is an
excellent audit witness, but it is a coordinate *gauge*: its three integers
change after an equivalent choice of unit-cell basis or atom representative.
Consequently, raw ``(quotient_atom_id, tau_x, tau_y, tau_z)`` tuples must not
be used as a cross-CIF identity key.

This module instead canonicalizes the physical local coordination object:

* donor attributes (element, bond role, and an optional chemical site tag),
* the equivalence partition saying which edges are images of the same
  quotient donor orbit, and
* the metal-centred Cartesian vectors reconstructed from the periodic edges.

The labelled Gram matrix is invariant to translation, rotation, atom order,
periodic-image gauge, and an exact unimodular change of lattice basis.  A
separate oriented-volume layer retains stereochemical handedness under proper
rotations.  Coordination numbers are deliberately limited to 2--6, so an
exhaustive donor permutation search (at most 720 permutations) is both simple
and independently auditable.

The resulting object is a canonical *local-state* identifier.  It is not a
canonical identifier for a complete MOF framework, pore, guest population, or
crystallographic species.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import math
import re
from collections import Counter
from dataclasses import dataclass
from typing import Any, Hashable, Iterable, Mapping, Optional, Sequence, Tuple

import numpy as np


_ELEMENT = re.compile(r"^[A-Z][a-z]?$" )


@dataclass(frozen=True)
class PeriodicDonorImage:
    """One physical donor image participating in a metal first sphere.

    ``fractional_base + translation`` is the selected physical image.  The
    base need not be wrapped; accepting an arbitrary representative makes the
    image-gauge invariance explicit and testable.  ``orbit`` may be any
    hashable identifier.  Only equality among orbit identifiers is used, so
    renaming them cannot affect the canonical result.
    """

    element: str
    fractional_base: Tuple[float, float, float]
    translation: Tuple[int, int, int]
    orbit: Hashable
    bond_role: str = "sigma"
    site_tag: str = ""

    @property
    def fractional_image(self) -> np.ndarray:
        return np.asarray(self.fractional_base, dtype=float) + np.asarray(
            self.translation, dtype=float
        )

    @property
    def attributes(self) -> Tuple[str, str, str]:
        return (str(self.element), str(self.bond_role), str(self.site_tag))


@dataclass(frozen=True)
class CanonicalPeriodicSite:
    """Canonical payload and identity keys for one local metal state."""

    metric_payload: str
    stereo_payload: str
    metric_id: str
    stereo_id: str
    token: str
    canonical_metric_order: Tuple[int, ...]
    canonical_stereo_order: Tuple[int, ...]
    decimals: int

    def as_dict(self) -> Mapping[str, Any]:
        return {
            "metric_payload": self.metric_payload,
            "stereo_payload": self.stereo_payload,
            "metric_id": self.metric_id,
            "stereo_id": self.stereo_id,
            "token": self.token,
            "canonical_metric_order": list(self.canonical_metric_order),
            "canonical_stereo_order": list(self.canonical_stereo_order),
            "decimals": self.decimals,
        }


def _validate_element(value: str, field: str) -> str:
    text = str(value)
    if _ELEMENT.fullmatch(text) is None:
        raise ValueError(f"Invalid {field} element: {value!r}")
    return text


def _validate_vectors(vectors: Sequence[Sequence[float]]) -> np.ndarray:
    array = np.asarray(vectors, dtype=float)
    if array.ndim != 2 or array.shape[1:] != (3,):
        raise ValueError("Donor vectors must have shape (CN, 3)")
    if array.shape[0] not in {2, 3, 4, 5, 6}:
        raise ValueError("Canonical periodic sites support CN 2 through 6")
    if not np.all(np.isfinite(array)):
        raise ValueError("Donor vectors contain nonfinite values")
    if np.any(np.linalg.norm(array, axis=1) <= 1.0e-12):
        raise ValueError("A donor vector is degenerate at the metal centre")
    return array


def _quantize(value: float, scale: int) -> int:
    if not math.isfinite(float(value)):
        raise ValueError("Cannot quantize a nonfinite invariant")
    return int(np.rint(float(value) * scale))


def _restricted_growth_partition(
    orbit_labels: Sequence[Hashable], permutation: Sequence[int]
) -> Tuple[int, ...]:
    """Relabel an orbit partition by first occurrence.

    This removes all dependence on the actual quotient-group identifiers.
    For example, ``[91, 17, 91, 23]`` becomes ``[0, 1, 0, 2]``.
    """

    mapping: dict[Hashable, int] = {}
    result = []
    for input_index in permutation:
        label = orbit_labels[int(input_index)]
        if label not in mapping:
            mapping[label] = len(mapping)
        result.append(mapping[label])
    return tuple(result)


def _flatten_upper(matrix: np.ndarray, permutation: Sequence[int]) -> Tuple[int, ...]:
    result = []
    for row in range(len(permutation)):
        for column in range(row, len(permutation)):
            result.append(int(matrix[permutation[row], permutation[column]]))
    return tuple(result)


def _flatten_oriented_volumes(
    tensor: np.ndarray, permutation: Sequence[int]
) -> Tuple[int, ...]:
    result = []
    for left in range(len(permutation)):
        for middle in range(left + 1, len(permutation)):
            for right in range(middle + 1, len(permutation)):
                result.append(
                    int(
                        tensor[
                            permutation[left],
                            permutation[middle],
                            permutation[right],
                        ]
                    )
                )
    return tuple(result)


def _payload(
    *,
    metal: str,
    attributes: Tuple[Tuple[str, str, str], ...],
    partition: Tuple[int, ...],
    gram: Tuple[int, ...],
    decimals: int,
    oriented_volumes: Optional[Tuple[int, ...]] = None,
) -> str:
    value: dict[str, Any] = {
        "v": "CR-PLS/2",
        "metal": metal,
        "cn": len(attributes),
        "donors": [list(item) for item in attributes],
        "orbit_partition": list(partition),
        "gram_A2_q": list(gram),
        "quantization_decimals": int(decimals),
    }
    if oriented_volumes is not None:
        value["oriented_volume_A3_q"] = list(oriented_volumes)
        value["frame_group"] = "SO(3)"
    else:
        value["frame_group"] = "O(3)"
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def _sha256(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def canonicalize_cartesian_site(
    *,
    metal: str,
    donor_elements: Sequence[str],
    donor_vectors: Sequence[Sequence[float]],
    donor_orbits: Sequence[Hashable],
    bond_roles: Optional[Sequence[str]] = None,
    site_tags: Optional[Sequence[str]] = None,
    decimals: int = 6,
) -> CanonicalPeriodicSite:
    """Canonicalize a labelled periodic first-sphere point set.

    The metric ID identifies configurations up to the full orthogonal group
    O(3), and therefore intentionally merges mirror images.  The stereo ID
    additionally includes signed triple products and identifies configurations
    up to proper rotations SO(3).  It is the appropriate key when the input
    Cartesian convention is known to be right-handed.
    """

    metal_text = _validate_element(metal, "metal")
    vectors = _validate_vectors(donor_vectors)
    cn = int(vectors.shape[0])
    if len(donor_elements) != cn or len(donor_orbits) != cn:
        raise ValueError("Donor elements, vectors, and orbits must have equal length")
    if bond_roles is None:
        bond_roles = ("sigma",) * cn
    if site_tags is None:
        site_tags = ("",) * cn
    if len(bond_roles) != cn or len(site_tags) != cn:
        raise ValueError("Bond roles and site tags must have length CN")
    if not isinstance(decimals, int) or decimals < 0 or decimals > 9:
        raise ValueError("decimals must be an integer between 0 and 9")

    attributes = tuple(
        (
            _validate_element(donor_elements[index], "donor"),
            str(bond_roles[index]),
            str(site_tags[index]),
        )
        for index in range(cn)
    )
    try:
        for orbit in donor_orbits:
            hash(orbit)
    except TypeError as exc:
        raise ValueError("Every donor orbit identifier must be hashable") from exc

    scale = 10**decimals
    gram_float = vectors @ vectors.T
    gram_quantized = np.rint(gram_float * scale).astype(object)
    # Build the complete ordered triple-product tensor once; permutation
    # candidates then only index this small integer array.
    oriented_quantized = np.empty((cn, cn, cn), dtype=object)
    for left in range(cn):
        for middle in range(cn):
            cross = np.cross(vectors[left], vectors[middle])
            for right in range(cn):
                oriented_quantized[left, middle, right] = _quantize(
                    float(np.dot(cross, vectors[right])), scale
                )

    metric_best: Optional[Tuple[Any, ...]] = None
    stereo_best: Optional[Tuple[Any, ...]] = None
    metric_order: Optional[Tuple[int, ...]] = None
    stereo_order: Optional[Tuple[int, ...]] = None
    for permutation in itertools.permutations(range(cn)):
        donor_key = tuple(attributes[index] for index in permutation)
        partition_key = _restricted_growth_partition(donor_orbits, permutation)
        gram_key = _flatten_upper(gram_quantized, permutation)
        shared = (donor_key, partition_key, gram_key)
        if metric_best is None or shared < metric_best:
            metric_best = shared
            metric_order = tuple(permutation)
        volume_key = _flatten_oriented_volumes(oriented_quantized, permutation)
        stereo = shared + (volume_key,)
        if stereo_best is None or stereo < stereo_best:
            stereo_best = stereo
            stereo_order = tuple(permutation)

    assert metric_best is not None and metric_order is not None
    assert stereo_best is not None and stereo_order is not None
    metric_payload = _payload(
        metal=metal_text,
        attributes=metric_best[0],
        partition=metric_best[1],
        gram=metric_best[2],
        decimals=decimals,
    )
    stereo_payload = _payload(
        metal=metal_text,
        attributes=stereo_best[0],
        partition=stereo_best[1],
        gram=stereo_best[2],
        decimals=decimals,
        oriented_volumes=stereo_best[3],
    )
    metric_id = _sha256(metric_payload)
    stereo_id = _sha256(stereo_payload)
    donor_counts = Counter(str(element) for element in donor_elements)
    donor_formula = ",".join(
        f"{element}:{donor_counts[element]}" for element in sorted(donor_counts)
    )
    orbit_pattern = ".".join(str(value) for value in stereo_best[1])
    token = (
        f"CR-PLS/2|M={metal_text}|CN={cn}|D={donor_formula}"
        f"|OR={orbit_pattern}|MID={metric_id[:20]}|SID={stereo_id[:20]}"
    )
    return CanonicalPeriodicSite(
        metric_payload=metric_payload,
        stereo_payload=stereo_payload,
        metric_id=metric_id,
        stereo_id=stereo_id,
        token=token,
        canonical_metric_order=metric_order,
        canonical_stereo_order=stereo_order,
        decimals=decimals,
    )


def canonicalize_fractional_site(
    *,
    metal: str,
    metal_fractional_base: Sequence[float],
    metal_translation: Sequence[int],
    donors: Sequence[PeriodicDonorImage],
    cell_matrix: Sequence[Sequence[float]],
    decimals: int = 6,
) -> CanonicalPeriodicSite:
    """Canonicalize a periodic site from fractional image data.

    Row-vector convention is used throughout: ``cartesian = fractional @ A``.
    Under a unimodular basis change ``A' = U @ A``, coordinates transform as
    ``f' = f @ inv(U)`` and produce the identical physical vectors.
    """

    matrix = np.asarray(cell_matrix, dtype=float)
    if matrix.shape != (3, 3) or not np.all(np.isfinite(matrix)):
        raise ValueError("cell_matrix must be a finite 3 by 3 matrix")
    determinant = float(np.linalg.det(matrix))
    if not math.isfinite(determinant) or abs(determinant) <= 1.0e-10:
        raise ValueError("cell_matrix is singular")
    metal_base = np.asarray(metal_fractional_base, dtype=float)
    metal_shift = np.asarray(metal_translation, dtype=float)
    if metal_base.shape != (3,) or metal_shift.shape != (3,):
        raise ValueError("Metal fractional base and translation must be length 3")
    metal_image = metal_base + metal_shift
    vectors = [
        (donor.fractional_image - metal_image) @ matrix for donor in donors
    ]
    return canonicalize_cartesian_site(
        metal=metal,
        donor_elements=[donor.element for donor in donors],
        donor_vectors=vectors,
        donor_orbits=[donor.orbit for donor in donors],
        bond_roles=[donor.bond_role for donor in donors],
        site_tags=[donor.site_tag for donor in donors],
        decimals=decimals,
    )


def canonical_entry_multiset_id(
    site_ids: Iterable[str], *, version: str = "CR-PERIODIC-ATLAS/1"
) -> str:
    """Hash a sorted multiset of local-state IDs, preserving multiplicity.

    This is useful for an entry-level atlas checksum.  It is intentionally not
    called a framework ID because first-sphere multisets omit linker and pore
    topology.
    """

    values = sorted(str(value) for value in site_ids)
    payload = json.dumps(
        {"v": str(version), "site_ids": values},
        sort_keys=True,
        separators=(",", ":"),
    )
    return _sha256(payload)
