"""Frozen element-domain policy for the all-metal CSD structure audit.

The symbols were obtained once from ``ccdc.molecule.Atom.is_metal`` under
CCDC Python API 3.6.0 and are stored explicitly here.  Production results do
not depend on a runtime definition that could drift between API releases.
"""

from __future__ import annotations

import hashlib
import json
from typing import Iterable, Tuple


METAL_POLICY_ID = "ccdc-3.6.0-is_metal-frozen-20260825-v1"

ALL_METALS: Tuple[str, ...] = tuple(
    "Li Be Na Mg Al K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge Rb Sr Y Zr Nb Mo "
    "Tc Ru Rh Pd Ag Cd In Sn Sb Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm "
    "Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po Fr Ra Ac Th Pa U Np Pu Am Cm "
    "Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og".split()
)
ALL_METAL_SET = frozenset(ALL_METALS)

S_BLOCK_METALS = frozenset("Li Be Na Mg K Ca Rb Sr Cs Ba Fr Ra".split())
D_BLOCK_METALS = frozenset(
    "Sc Ti V Cr Mn Fe Co Ni Cu Zn Y Zr Nb Mo Tc Ru Rh Pd Ag Cd Hf Ta W Re Os "
    "Ir Pt Au Hg Rf Db Sg Bh Hs Mt Ds Rg Cn".split()
)
F_BLOCK_METALS = frozenset(
    "La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Ac Th Pa U Np Pu Am Cm Bk "
    "Cf Es Fm Md No Lr".split()
)
P_BLOCK_METALS = frozenset(ALL_METAL_SET - S_BLOCK_METALS - D_BLOCK_METALS - F_BLOCK_METALS)

METAL_POLICY_SHA256 = hashlib.sha256(
    json.dumps(
        {
            "policy_id": METAL_POLICY_ID,
            "symbols": ALL_METALS,
            "s": sorted(S_BLOCK_METALS),
            "p": sorted(P_BLOCK_METALS),
            "d": sorted(D_BLOCK_METALS),
            "f": sorted(F_BLOCK_METALS),
        },
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
).hexdigest()


def is_target_metal_symbol(symbol: object) -> bool:
    """Return whether *symbol* belongs to the frozen all-metal domain."""

    return str(symbol) in ALL_METAL_SET


def metal_block(symbol: object) -> str:
    """Return the periodic-table block used for stratified audit summaries."""

    text = str(symbol)
    if text in S_BLOCK_METALS:
        return "s"
    if text in P_BLOCK_METALS:
        return "p"
    if text in D_BLOCK_METALS:
        return "d"
    if text in F_BLOCK_METALS:
        return "f"
    return "outside"


def metal_blocks(symbols: Iterable[object]) -> Tuple[str, ...]:
    """Return sorted unique block labels for a collection of symbols."""

    order = {"s": 0, "p": 1, "d": 2, "f": 3, "outside": 4}
    return tuple(sorted({metal_block(symbol) for symbol in symbols}, key=order.__getitem__))


if len(ALL_METALS) != len(ALL_METAL_SET):  # pragma: no cover - import-time guard
    raise RuntimeError("the frozen all-metal policy contains duplicate symbols")
if S_BLOCK_METALS | P_BLOCK_METALS | D_BLOCK_METALS | F_BLOCK_METALS != ALL_METAL_SET:
    raise RuntimeError("metal block sets do not close to the frozen all-metal domain")
if any(
    left & right
    for index, left in enumerate((S_BLOCK_METALS, P_BLOCK_METALS, D_BLOCK_METALS, F_BLOCK_METALS))
    for right in (S_BLOCK_METALS, P_BLOCK_METALS, D_BLOCK_METALS, F_BLOCK_METALS)[index + 1 :]
):  # pragma: no cover - import-time guard
    raise RuntimeError("metal block sets overlap")
