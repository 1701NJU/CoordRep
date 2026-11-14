#!/usr/bin/env python3
"""Print the CCDC metal classification for the periodic table.

This is a provenance probe only.  The resulting symbols are copied into the
frozen audit policy rather than queried dynamically during production runs.
"""

from __future__ import annotations

import json

import ccdc
from ccdc.molecule import Atom


ELEMENTS = (
    "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn "
    "Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce "
    "Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn "
    "Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl "
    "Mc Lv Ts Og"
).split()


def main() -> None:
    rows = []
    for atomic_number, symbol in enumerate(ELEMENTS, start=1):
        try:
            atom = Atom(atomic_symbol=symbol)
        except TypeError:
            atom = Atom(symbol)
        try:
            is_metal = bool(atom.is_metal)
            error = None
        except Exception as exc:  # pragma: no cover - diagnostic probe
            is_metal = None
            error = f"{type(exc).__name__}: {exc}"
        rows.append({
            "atomic_number": atomic_number,
            "symbol": symbol,
            "is_metal": is_metal,
            "error": error,
        })
    payload = {
        "ccdc_version": str(getattr(ccdc, "__version__", "unknown")),
        "metal_symbols": [row["symbol"] for row in rows if row["is_metal"]],
        "rows": rows,
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
