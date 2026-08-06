"""
CoordRep Core Data Structures

Defines the canonical data model for coordination complex representation:
CoordComplex, MetalState, ShapeVector, LigandModule, ConstraintSet, etc.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Any
import numpy as np


@dataclass
class MetalState:
    """Metal center state."""
    element: str
    oxidation: Optional[int] = None
    spin: Optional[str] = None             # "LS" / "HS" / None
    dcount: Optional[int] = None           # d-electron count
    phys_vec: Optional[np.ndarray] = None  # [row, group, EN, radius, ...]

    def to_dict(self) -> dict:
        return {
            'element': self.element,
            'oxidation': self.oxidation,
            'spin': self.spin,
            'dcount': self.dcount,
        }


@dataclass
class ShapeVector:
    """
    Continuous Shape Measures (CShM) vector.

    Encodes geometry as continuous coordinates on a shape manifold
    rather than discrete labels. E.g. Zn is not "Td or SP" but
    "(Td_dist=0.2, SP_dist=12.5)".
    """
    cn: int
    ref_shapes: List[str]                  # e.g., ["Td", "SP"]
    values: np.ndarray                     # raw CShM values
    values_rounded: Optional[np.ndarray] = None

    def round(self, decimals: int = 2) -> 'ShapeVector':
        """Round to specified decimal places for canonical serialization."""
        self.values_rounded = np.round(self.values, decimals)
        return self

    def to_dict(self) -> dict:
        return {
            'cn': self.cn,
            'ref_shapes': self.ref_shapes,
            'values': self.values.tolist() if self.values is not None else None,
            'values_rounded': self.values_rounded.tolist() if self.values_rounded is not None else None,
        }

    def __repr__(self):
        if self.values_rounded is not None:
            parts = [f"{s}={v:.2f}" for s, v in zip(self.ref_shapes, self.values_rounded)]
        elif self.values is not None:
            parts = [f"{s}={v:.2f}" for s, v in zip(self.ref_shapes, self.values)]
        else:
            parts = []
        return f"<Shape CN={self.cn}: {', '.join(parts)}>"


@dataclass
class LigandModule:
    """
    Ligand module with attachment points.

    Each ligand is an independent, reusable module identified by
    canonical SMILES and donor-atom markers.
    """
    lig_id: str                            # L1, L2, ...
    smiles: str                            # RDKit canonical SMILES
    attach_atoms: List[int]                # donor atom indices within ligand graph
    donor_elements: List[str]              # donor atom elements, e.g. ['N', 'N'] for en
    dent: int                              # denticity: 1 / 2 / 3 ...
    eta: Optional[int] = None              # hapticity (if applicable)
    charge: Optional[int] = None
    meta: Dict[str, Any] = field(default_factory=dict)  # MW, etc.
    payload_provenance: str = "SMILES"       # "SMILES" / "FORMULA"
    connectivity_status: str = "unspecified"

    def get_sort_key(self) -> tuple:
        """
        Deterministic ligand sort key for canonicalization.

        Priority (descending):
            1. denticity (desc)
            2. donor atomic-number multiset (lex)
            3. hapticity (desc, None last)
            4. |charge| (desc, None last)
            5. molecular weight (desc)
            6. SMILES (lex)
        """
        ATOMIC_NUMBERS = {
            'H': 1, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'P': 15, 'S': 16, 'Cl': 17,
            'Br': 35, 'I': 53
        }

        donor_z = sorted([ATOMIC_NUMBERS.get(e, 0) for e in self.donor_elements], reverse=True)

        return (
            -self.dent,
            tuple(donor_z),
            -(self.eta or 0),
            -(abs(self.charge) if self.charge else 0),
            -self.meta.get('mw', 0),
            self.payload_provenance,
            self.smiles,
        )

    def to_dict(self) -> dict:
        return {
            'lig_id': self.lig_id,
            'smiles': self.smiles,
            'dent': self.dent,
            'eta': self.eta,
            'charge': self.charge,
            'donor_elements': self.donor_elements,
            'payload_provenance': self.payload_provenance,
            'connectivity_status': self.connectivity_status,
        }


@dataclass
class DonorSite:
    """A coordination site used in constraint expressions."""
    lig_id: str          # L1, L2, ...
    donor_element: str   # N, O, Cl, ...
    donor_rank: int      # rank among same-element donors in this ligand (1-based)

    def __str__(self):
        return f"{self.lig_id}:{self.donor_element}:{self.donor_rank}"

    def __lt__(self, other):
        return str(self) < str(other)

    def __eq__(self, other):
        return str(self) == str(other)

    def __hash__(self):
        return hash(str(self))


@dataclass
class ConstraintSet:
    """
    Relative configuration constraints (trans/cis/fac/mer).

    Expressed purely in terms of relational invariants between
    donor sites, independent of slot numbering.
    """
    trans_pairs: List[Tuple[DonorSite, DonorSite]] = field(default_factory=list)
    cis_pairs: List[Tuple[DonorSite, DonorSite]] = field(default_factory=list)
    notes: Dict[str, Any] = field(default_factory=dict)

    def canonicalize(self) -> 'ConstraintSet':
        """Canonicalize constraints by sorting pairs."""
        def sort_pair(pair):
            a, b = pair
            return (min(a, b), max(a, b))

        self.trans_pairs = sorted([sort_pair(p) for p in self.trans_pairs])
        self.cis_pairs = sorted([sort_pair(p) for p in self.cis_pairs])
        return self

    def to_dict(self) -> dict:
        return {
            'trans_pairs': [(str(a), str(b)) for a, b in self.trans_pairs],
            'cis_pairs': [(str(a), str(b)) for a, b in self.cis_pairs],
            'notes': self.notes,
        }


@dataclass
class AssemblyGraph:
    """Metal--ligand assembly graph."""
    metal_idx: int
    donor_indices: List[int]
    lig_assignments: Dict[int, str]        # donor_idx -> lig_id
    bond_orders: Dict[int, float]          # donor_idx -> M-L bond order

    def to_dict(self) -> dict:
        return {
            'metal_idx': self.metal_idx,
            'donor_indices': self.donor_indices,
            'n_donors': len(self.donor_indices),
        }


@dataclass
class CoordRepIssue:
    """Encoding issue or warning."""
    code: str
    message: str
    severity: str       # "error" / "warning" / "info"
    context: Dict = field(default_factory=dict)

    def __str__(self):
        return f"[{self.severity.upper()}] {self.code}: {self.message}"


class IssueCodes:
    NO_METAL = "NO_METAL"
    MULTI_METAL = "MULTI_METAL"
    NO_DONORS = "NO_DONORS"

    SHAPE_COMPUTATION_FAILED = "SHAPE_COMPUTATION_FAILED"
    SHAPE_REFERENCE_UNSUPPORTED = "SHAPE_REFERENCE_UNSUPPORTED"
    RDKIT_SANITIZE_FAILED = "RDKIT_SANITIZE_FAILED"
    FORMULA_FALLBACK = "FORMULA_FALLBACK"

    LOW_BOND_ORDER = "LOW_BOND_ORDER"
    UNUSUAL_GEOMETRY = "UNUSUAL_GEOMETRY"
    AMBIGUOUS_OXIDATION = "AMBIGUOUS_OXIDATION"

    HIGH_COORDINATION = "HIGH_COORDINATION"


@dataclass
class CoordComplex:
    """
    Central CoordRep object.

    Guarantees: the same structure produces a bit-identical canonical
    string under arbitrary atom reordering, coordinate rotation/translation,
    and minor numerical noise (after rounding).
    """
    metal: MetalState
    ligands: List[LigandModule]
    graph: AssemblyGraph
    shape: Optional[ShapeVector]
    constraints: ConstraintSet

    source_id: str = ""
    is_canonical: bool = False
    issues: List[CoordRepIssue] = field(default_factory=list)

    def canonicalize(self) -> 'CoordComplex':
        """Return a canonicalized copy of this CoordComplex."""
        from .canonical.canonicalize import canonicalize_complex
        return canonicalize_complex(self)

    def to_string(self) -> str:
        """
        Serialize to a CoordRep string.

        Format:
            [Metal:Fe|ox:+2|d:d6|CN:6]
            <ShapeBest:Oh|Class:ideal|Delta:3|V:0.23,8.10>
            {trans:L1:N:1--L2:Cl:1}
            |L1=NCCN|L2=Cl|
        """
        from .serialize.to_string import serialize_complex
        return serialize_complex(self)

    def to_features(self) -> dict:
        """
        Extract ML feature vectors.

        Returns:
            dict with keys "tokens", "geom", "metal", "ligands".
        """
        from .features.to_features import extract_features
        return extract_features(self)

    def validate(self) -> List[CoordRepIssue]:
        """Validate coordination constraints for consistency."""
        from .validate.checks import validate_complex
        return validate_complex(self)

    def to_dict(self) -> dict:
        """Convert to a JSON-serializable dictionary."""
        return {
            'source_id': self.source_id,
            'metal': self.metal.to_dict(),
            'ligands': [lig.to_dict() for lig in self.ligands],
            'shape': self.shape.to_dict() if self.shape else None,
            'constraints': self.constraints.to_dict(),
            'is_canonical': self.is_canonical,
            'issues': [str(i) for i in self.issues],
        }


@dataclass
class CoordRepConfig:
    """Encoding configuration."""
    bo_threshold_mode: str = "adaptive"
    bo_threshold_fixed: float = 0.3

    rounding_decimals: int = 2
    cn4_ref_shapes: List[str] = field(default_factory=lambda: ["Td", "SP"])
    cn6_ref_shapes: List[str] = field(default_factory=lambda: ["Oh", "TP"])

    trans_angle_range: Tuple[float, float] = (165.0, 195.0)
    cis_angle_range: Tuple[float, float] = (75.0, 105.0)

    allow_multimetal: bool = False
    strict_rdkit_sanitize: bool = False

    strict_invariance_cn_max: int = 6
    high_cn_policy: str = "equivalence_class"

    @classmethod
    def default(cls) -> 'CoordRepConfig':
        return cls()
