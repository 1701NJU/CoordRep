"""Small, lossless parser for the CoordRep stereochemical constraint block."""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Tuple


_TRANS_RE = re.compile(r"\{trans:([^}]+)\}")
_CIS_RE = re.compile(r"\{cis:([^}]+)\}")
_FM_RE = re.compile(r"\{fm:([^}]+)\}")


@dataclass(frozen=True)
class ParsedConstraintBlock:
    """Exact string payloads from a serialized constraint block."""

    trans: Tuple[str, ...] = ()
    cis: Tuple[str, ...] = ()
    fac_mer: Tuple[str, ...] = ()

    def to_string(self) -> str:
        """Return the canonical serializer ordering for this block."""
        return (
            "".join(f"{{trans:{item}}}" for item in self.trans)
            + "".join(f"{{cis:{item}}}" for item in self.cis)
            + "".join(f"{{fm:{item}}}" for item in self.fac_mer)
        )


def parse_constraint_block(coordrep_string: str) -> ParsedConstraintBlock:
    """Parse all trans/cis/fac-mer relations without discarding site labels."""
    return ParsedConstraintBlock(
        trans=tuple(_TRANS_RE.findall(coordrep_string)),
        cis=tuple(_CIS_RE.findall(coordrep_string)),
        fac_mer=tuple(_FM_RE.findall(coordrep_string)),
    )
