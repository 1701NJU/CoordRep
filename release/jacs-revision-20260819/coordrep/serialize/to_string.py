"""
Serialization module.

Converts a CoordComplex into its canonical CoordRep string using the
bucketed shape format:
    <ShapeBest:Oh|Class:ideal|Delta:2|V:0.23,8.10>
"""

from ..core import CoordComplex


def serialize_complex(cc: CoordComplex) -> str:
    """
    Serialize a CoordComplex to a CoordRep string.

    Format:
        [Metal:Fe|ox:+2|d:d6|CN:6]
        <ShapeBest:Oh|Class:ideal|Delta:2|V:...>
        {trans:L1:N:1--L2:Cl:1}
        |L1=...|L2=...|
    """
    parts = []

    # 1. Metal block
    parts.append(_serialize_metal(cc))

    # 2. Shape block (bucketed format)
    shape_str = _serialize_shape(cc)
    if shape_str:
        parts.append(shape_str)

    # 3. Constraints block
    constraints_str = _serialize_constraints(cc)
    if constraints_str:
        parts.append(constraints_str)

    # 4. Ligand dictionary
    parts.append(_serialize_ligands(cc))

    return ''.join(parts)


def _serialize_metal(cc: CoordComplex) -> str:
    """Serialize the metal block."""
    m = cc.metal

    items = [f"Metal:{m.element}"]

    if m.oxidation is not None:
        ox_str = f"+{m.oxidation}" if m.oxidation >= 0 else str(m.oxidation)
        items.append(f"ox:{ox_str}")

    if m.dcount is not None:
        items.append(f"d:d{m.dcount}")

    cn = len(cc.graph.donor_indices) if cc.graph else 0
    items.append(f"CN:{cn}")

    return '[' + '|'.join(items) + ']'


def _serialize_shape(cc: CoordComplex) -> str:
    """
    Serialize the shape block using bucketed format.

    Output: <ShapeBest:Oh|Class:ideal|Delta:2|V:0.23,8.10>
    """
    s = cc.shape

    if s is None:
        cn = len(cc.graph.donor_indices) if cc.graph else 0
        return (
            f"<ShapeStatus:unsupported|Reason:missing-shape-object|CN:{cn}>"
        )

    if len(s.ref_shapes) == 0:
        return (
            f"<ShapeStatus:unsupported|Reason:no-reference-implementation"
            f"|CN:{s.cn}>"
        )

    values = s.values_rounded if s.values_rounded is not None else s.values

    if values is None or len(values) == 0:
        return (
            f"<ShapeStatus:unsupported|Reason:empty-cshm-vector|CN:{s.cn}>"
        )

    min_idx = 0
    min_val = float('inf')
    for i, v in enumerate(values):
        if v < min_val:
            min_val = v
            min_idx = i

    best_shape = s.ref_shapes[min_idx]
    best_cshm = values[min_idx]

    sorted_vals = sorted(values)
    if len(sorted_vals) > 1:
        delta = sorted_vals[1] - sorted_vals[0]
    else:
        delta = 0

    # Delta bucketing
    if delta < 0.5:
        delta_bin = 0
    elif delta < 2:
        delta_bin = 1
    elif delta < 5:
        delta_bin = 2
    else:
        delta_bin = 3

    # Shape classification
    if best_cshm < 1:
        shape_class = "ideal"
    elif best_cshm < 3:
        shape_class = "good"
    elif best_cshm < 8:
        shape_class = "dist"
    else:
        shape_class = "irreg"

    vals_str = ','.join(f"{v:.2f}" for v in values)

    return f"<ShapeBest:{best_shape}|Class:{shape_class}|Delta:{delta_bin}|V:{vals_str}>"


def _serialize_constraints(cc: CoordComplex) -> str:
    """Serialize the stereochemical constraints block."""
    c = cc.constraints
    parts = []

    for s1, s2 in c.trans_pairs:
        parts.append(f"{{trans:{s1}--{s2}}}")

    for s1, s2 in c.cis_pairs:
        parts.append(f"{{cis:{s1}--{s2}}}")

    if 'fac_mer' in c.notes:
        parts.append(f"{{fm:{c.notes['fac_mer']}}}")

    return ''.join(parts)


def _serialize_ligands(cc: CoordComplex) -> str:
    """Serialize the ligand dictionary."""
    if not cc.ligands:
        return '||'

    items = []
    for lig in cc.ligands:
        provenance = (lig.payload_provenance or "SMILES").upper()
        if provenance not in {"SMILES", "FORMULA"}:
            provenance = "FORMULA"
        items.append(f"{lig.lig_id}={provenance}:{lig.smiles}")

    return '|' + '|'.join(items) + '|'
