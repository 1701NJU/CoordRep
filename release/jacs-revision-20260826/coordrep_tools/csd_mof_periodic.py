"""Canonical CoordRep-compatible payload for periodic local metal sites.

This is intentionally a *local-state* adapter, not a serializer for a complete
MOF framework.  It carries only fields already defined by the molecular
CoordRep coordination-state layer: metal, coordination number, donor-element
multiset, best reference shape, and corrected coordinate-space CShM.

Format
------
``CR-PLS/1|M=Cu|CN=6|D=Cl:2,N:4|SH=Oh|S=1.8748``

The grammar is deliberately small, deterministic, and coordinate-free so
that it can be round-tripped without redistributing source CIF coordinates.
"""

from __future__ import annotations

import re
from collections import Counter
from dataclasses import dataclass
from typing import Iterable, Tuple


_TOKEN_PATTERN = re.compile(
    r"^CR-PLS/1\|M=([A-Z][a-z]?)\|CN=([2-6])\|D=([^|]+)"
    r"\|SH=([A-Za-z0-9]+)\|S=([0-9]+\.[0-9]{4})$"
)


@dataclass(frozen=True)
class PeriodicLocalState:
    """Coordinate-free periodic first-sphere state payload."""

    metal: str
    cn: int
    donors: Tuple[Tuple[str, int], ...]
    best_shape: str
    best_cshm: float

    @classmethod
    def from_fields(
        cls,
        metal: str,
        cn: int,
        donor_elements: Iterable[str],
        best_shape: str,
        best_cshm: float,
    ) -> "PeriodicLocalState":
        counts = Counter(str(element) for element in donor_elements)
        donors = tuple(sorted(counts.items()))
        state = cls(
            metal=str(metal),
            cn=int(cn),
            donors=donors,
            best_shape=str(best_shape),
            best_cshm=round(float(best_cshm), 4),
        )
        state.validate()
        return state

    def validate(self) -> None:
        if not re.fullmatch(r"[A-Z][a-z]?", self.metal):
            raise ValueError(f"Invalid metal element: {self.metal!r}")
        if self.cn not in {2, 3, 4, 5, 6}:
            raise ValueError(f"Unsupported coordination number: {self.cn}")
        if sum(count for _, count in self.donors) != self.cn:
            raise ValueError("Donor multiplicities do not sum to CN")
        if not self.donors:
            raise ValueError("At least one donor element is required")
        previous = ""
        for element, count in self.donors:
            if not re.fullmatch(r"[A-Z][a-z]?", element):
                raise ValueError(f"Invalid donor element: {element!r}")
            if element <= previous:
                raise ValueError("Donor elements must be unique and sorted")
            if int(count) <= 0:
                raise ValueError("Donor multiplicities must be positive")
            previous = element
        if not re.fullmatch(r"[A-Za-z0-9]+", self.best_shape):
            raise ValueError(f"Invalid shape label: {self.best_shape!r}")
        if self.best_cshm < 0.0:
            raise ValueError("CShM must be nonnegative")

    def to_token(self) -> str:
        self.validate()
        donor_text = ",".join(
            f"{element}:{count}" for element, count in self.donors
        )
        return (
            f"CR-PLS/1|M={self.metal}|CN={self.cn}|D={donor_text}"
            f"|SH={self.best_shape}|S={self.best_cshm:.4f}"
        )

    @classmethod
    def from_token(cls, token: str) -> "PeriodicLocalState":
        match = _TOKEN_PATTERN.fullmatch(str(token))
        if match is None:
            raise ValueError("Invalid CR-PLS/1 token")
        metal, cn_text, donor_text, shape, cshm_text = match.groups()
        donors = []
        for item in donor_text.split(","):
            parts = item.split(":")
            if len(parts) != 2:
                raise ValueError("Invalid donor field")
            donors.append((parts[0], int(parts[1])))
        state = cls(
            metal=metal,
            cn=int(cn_text),
            donors=tuple(donors),
            best_shape=shape,
            best_cshm=float(cshm_text),
        )
        state.validate()
        return state
