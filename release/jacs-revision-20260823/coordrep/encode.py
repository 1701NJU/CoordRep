"""
Main Encoding API

Converts raw coordination complex structures (XYZ, CIF) into canonical
CoordRep objects with continuous shape measures, ligand decomposition,
and stereochemical constraints.
"""

from pathlib import Path
from typing import List, Optional
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed

from .core import (
    CoordComplex, MetalState, AssemblyGraph, ConstraintSet,
    CoordRepConfig, CoordRepIssue, IssueCodes, ShapeVector, LigandModule
)
from .io.tmqm_reader import TMQMReader, RawMolecule, Atom, TRANSITION_METALS
from .graph.neighbors import get_coordination_neighbors, get_donor_indices, DONOR_ELEMENTS
from .graph.build_graph import GraphBuilder, METAL_INFO
from .graph.ligand_module import LigandExtractor
from .graph.donor_sites import build_all_donor_sites
from .geometry.shape import ShapeCalculator, ShapeResult
from .geometry.rel_config import ConstraintDetector


def encode(xyz_path: str = None,
          bo_path: str = None,
          molecule: RawMolecule = None,
          config: CoordRepConfig = None) -> CoordComplex:
    """
    Encode XYZ + BO files (or a RawMolecule) into a CoordComplex.

    Args:
        xyz_path: Path to the XYZ coordinate file.
        bo_path: Path to the Wiberg bond-order file (optional).
        molecule: Pre-loaded RawMolecule (overrides xyz/bo paths).
        config: CoordRepConfig instance for tuning thresholds.

    Returns:
        CoordComplex ready for canonicalization and serialization.
    """
    config = config or CoordRepConfig.default()

    if molecule is None:
        molecule = _load_molecule(xyz_path, bo_path)

    return encode_molecule(molecule, config)


def encode_molecule(mol: RawMolecule,
                   config: CoordRepConfig = None) -> CoordComplex:
    """
    Encode a RawMolecule into a CoordComplex.

    Pipeline:
        1. Identify metal center
        2. Determine coordination neighbors (H filtered)
        3. Build molecular graph
        4. Extract ligand modules
        5. Compute assignment-invariant shape vector (CShM)
        6. Build donor-site mapping
        7. Detect stereochemical constraints
        8. Assemble CoordComplex
    """
    config = config or CoordRepConfig.default()
    all_issues = []

    # 1. Find metal center
    metal_idx = _find_metal(mol)
    if metal_idx is None:
        return _create_failed_complex(mol.mol_id, "No transition metal found")

    metal_element = mol.atoms[metal_idx].element

    # 2. Determine coordination neighbors (H excluded)
    neighbors = get_coordination_neighbors(
        mol, metal_idx,
        bo_threshold=config.bo_threshold_fixed,
        distance_cutoff=2.8
    )

    if not neighbors:
        return _create_failed_complex(mol.mol_id, "No donors found (H excluded)")

    donor_indices = get_donor_indices(neighbors)
    cn = len(donor_indices)

    # 3. Build molecular graph
    graph_builder = GraphBuilder(config)
    try:
        graph_result = graph_builder.build(mol)
        mol_graph = graph_result.molecular_graph
    except Exception as e:
        mol_graph = None
        all_issues.append(CoordRepIssue(
            code="GRAPH_BUILD_FAILED",
            message=str(e),
            severity="warning"
        ))

    # 4. Extract ligands
    ligand_extractor = LigandExtractor(config.strict_rdkit_sanitize)
    ligands, lig_issues = ligand_extractor.extract_ligands(
        mol,
        mol_graph,
        metal_idx,
        donor_indices
    )
    all_issues.extend(lig_issues)

    # 5. Compute assignment-invariant shape vector
    donor_coords = np.array([mol.atoms[i].coords for i in donor_indices])
    metal_coord = mol.atoms[metal_idx].coords

    shape_calc = ShapeCalculator(rounding_decimals=config.rounding_decimals)
    try:
        shape_result = shape_calc.compute(donor_coords, metal_coord)

        shape = ShapeVector(
            cn=shape_result.cn,
            ref_shapes=shape_result.ref_shapes,
            values=shape_result.values,
            values_rounded=shape_result.values_rounded
        )
    except Exception as e:
        shape = None
        all_issues.append(CoordRepIssue(
            code=IssueCodes.SHAPE_COMPUTATION_FAILED,
            message=str(e),
            severity="warning"
        ))

    # 6. Build donor-site mapping
    donor_sites = build_all_donor_sites(mol, ligands, metal_idx, mol_graph)

    # 7. Detect stereochemical constraints
    constraint_detector = ConstraintDetector(config)
    constraints = constraint_detector.detect_constraints_with_sites(
        mol, metal_idx, donor_indices, ligands, donor_sites
    )

    # 8. Create metal state
    metal_state = _create_metal_state(metal_element, cn)

    # 9. Build AssemblyGraph
    donor_bond_orders = {n.atom_idx: n.bond_order for n in neighbors}
    assembly_graph = AssemblyGraph(
        metal_idx=metal_idx,
        donor_indices=donor_indices,
        lig_assignments={d: lig.lig_id for lig in ligands for d in lig.attach_atoms if d in donor_indices},
        bond_orders=donor_bond_orders
    )

    # 10. Assemble CoordComplex
    cc = CoordComplex(
        metal=metal_state,
        ligands=ligands,
        graph=assembly_graph,
        shape=shape,
        constraints=constraints,
        source_id=mol.mol_id,
        is_canonical=False,
        issues=all_issues
    )

    return cc


def _find_metal(mol: RawMolecule) -> Optional[int]:
    """Find the transition-metal atom index."""
    for i, atom in enumerate(mol.atoms):
        if atom.element in TRANSITION_METALS:
            return i
    return None


def _create_metal_state(element: str, cn: int) -> MetalState:
    """Create a MetalState with estimated oxidation and d-electron count."""
    info = METAL_INFO.get(element, {})

    oxidation = _estimate_oxidation(element, cn)

    dcount = None
    if oxidation is not None:
        group = info.get('group', 0)
        if 3 <= group <= 12:
            dcount = max(0, min(10, group - oxidation))

    phys_vec = np.array([
        info.get('row', 0),
        info.get('group', 0),
        info.get('en', 0),
    ])

    return MetalState(
        element=element,
        oxidation=oxidation,
        spin=None,
        dcount=dcount,
        phys_vec=phys_vec
    )


def _estimate_oxidation(element: str, cn: int) -> Optional[int]:
    """Estimate oxidation state from common (element, CN) patterns."""
    common_ox = {
        'Fe': {4: 2, 5: 2, 6: 2},
        'Co': {4: 2, 5: 2, 6: 3},
        'Ni': {4: 2, 5: 2, 6: 2},
        'Cu': {4: 2, 5: 2, 6: 2},
        'Zn': {4: 2, 5: 2, 6: 2},
        'Pd': {4: 2, 5: 2, 6: 4},
        'Pt': {4: 2, 5: 2, 6: 4},
        'Ru': {5: 2, 6: 2},
        'Rh': {4: 1, 5: 3, 6: 3},
    }

    if element in common_ox:
        return common_ox[element].get(cn, 2)
    return None


def batch_encode(molecules: List[RawMolecule],
                config: CoordRepConfig = None,
                n_jobs: int = 1,
                show_progress: bool = True) -> List[CoordComplex]:
    """
    Encode a list of RawMolecule objects into CoordComplex objects.

    Args:
        molecules: List of RawMolecule instances.
        config: CoordRepConfig for encoding parameters.
        n_jobs: Number of parallel workers (currently sequential only).
        show_progress: Show a tqdm progress bar if available.

    Returns:
        List of CoordComplex objects (failed conversions return
        a placeholder with issues).
    """
    config = config or CoordRepConfig.default()

    results = []
    iterator = molecules
    if show_progress:
        try:
            from tqdm import tqdm
            iterator = tqdm(molecules, desc="Encoding")
        except ImportError:
            pass

    for mol in iterator:
        try:
            cc = encode_molecule(mol, config)
            results.append(cc)
        except Exception as e:
            results.append(_create_failed_complex(mol.mol_id, str(e)))

    return results


def _load_molecule(xyz_path: str, bo_path: str = None) -> RawMolecule:
    """Load a molecule from XYZ (+ optional BO) files."""
    xyz_path = Path(xyz_path)

    atoms = []
    mol_id = xyz_path.stem

    with open(xyz_path, 'r') as f:
        lines = f.readlines()

    n_atoms = int(lines[0].strip())

    for i, line in enumerate(lines[2:n_atoms + 2]):
        parts = line.strip().split()
        if len(parts) >= 4:
            atoms.append(Atom(
                index=i,
                element=parts[0],
                x=float(parts[1]),
                y=float(parts[2]),
                z=float(parts[3])
            ))

    bond_orders = None
    if bo_path:
        bo_path = Path(bo_path)
        if bo_path.exists():
            rows = []
            with open(bo_path, 'r') as f:
                for line in f:
                    if line.strip() and not line.startswith('#'):
                        try:
                            values = [float(x) for x in line.split()]
                            if values:
                                rows.append(values)
                        except ValueError:
                            continue
            if rows:
                bond_orders = np.array(rows)

    return RawMolecule(
        mol_id=mol_id,
        atoms=atoms,
        bond_orders=bond_orders
    )


def _create_failed_complex(mol_id: str, error_msg: str) -> CoordComplex:
    """Create a placeholder CoordComplex for failed encoding."""
    return CoordComplex(
        metal=MetalState(element="?", oxidation=None),
        ligands=[],
        graph=None,
        shape=None,
        constraints=ConstraintSet(),
        source_id=mol_id,
        is_canonical=False,
        issues=[CoordRepIssue(
            code="ENCODING_FAILED",
            message=error_msg,
            severity="error"
        )]
    )
