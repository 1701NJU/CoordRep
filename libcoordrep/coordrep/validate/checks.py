"""
Validation Checks

"""

from typing import List

from ..core import CoordComplex, CoordRepIssue, IssueCodes


def validate_complex(cc: CoordComplex) -> List[CoordRepIssue]:
    """
 CoordComplex
    
    Returns:
    """
    issues = []
    
    # 1. Basic checks
    issues.extend(_check_basic_structure(cc))
    
    issues.extend(_check_geometry_consistency(cc))
    
    issues.extend(_check_ligand_consistency(cc))
    
    issues.extend(_check_constraint_consistency(cc))
    
    return issues


def _check_basic_structure(cc: CoordComplex) -> List[CoordRepIssue]:
 """"""
    issues = []
    
    if cc.metal is None:
        issues.append(CoordRepIssue(
            code=IssueCodes.NO_METAL,
            message="Metal state is missing",
            severity="error"
        ))
    
    if cc.graph:
        cn = len(cc.graph.donor_indices)
        if cn == 0:
            issues.append(CoordRepIssue(
                code=IssueCodes.NO_DONORS,
                message="No donor atoms found",
                severity="error"
            ))
        elif cn > 12:
            issues.append(CoordRepIssue(
                code=IssueCodes.UNUSUAL_GEOMETRY,
                message=f"Unusually high coordination number: {cn}",
                severity="warning"
            ))
    
    return issues


def _check_geometry_consistency(cc: CoordComplex) -> List[CoordRepIssue]:
 """"""
    issues = []
    
    if cc.shape is None:
        issues.append(CoordRepIssue(
            code=IssueCodes.SHAPE_COMPUTATION_FAILED,
            message="Shape vector is missing",
            severity="warning"
        ))
        return issues
    
    expected_cn = cc.shape.cn
    if cc.graph:
        actual_cn = len(cc.graph.donor_indices)
        if expected_cn != actual_cn:
            issues.append(CoordRepIssue(
                code=IssueCodes.UNUSUAL_GEOMETRY,
                message=f"CN mismatch: shape CN={expected_cn}, actual CN={actual_cn}",
                severity="warning"
            ))
    
    if cc.shape.values is not None:
        min_cshm = min(cc.shape.values)
        if min_cshm > 10:
            issues.append(CoordRepIssue(
                code=IssueCodes.UNUSUAL_GEOMETRY,
                message=f"Very distorted geometry: min CShM = {min_cshm:.2f}",
                severity="warning",
                context={'min_cshm': min_cshm}
            ))
    
    return issues


def _check_ligand_consistency(cc: CoordComplex) -> List[CoordRepIssue]:
 """"""
    issues = []
    
    for lig in cc.ligands:
        if lig.dent != len(lig.donor_elements):
            issues.append(CoordRepIssue(
                code="LIGAND_INCONSISTENCY",
                message=f"{lig.lig_id}: dent={lig.dent} but {len(lig.donor_elements)} donor elements",
                severity="warning",
                context={'lig_id': lig.lig_id}
            ))
        
        if not lig.smiles:
            issues.append(CoordRepIssue(
                code=IssueCodes.RDKIT_SANITIZE_FAILED,
                message=f"{lig.lig_id}: Empty SMILES",
                severity="warning"
            ))
    
    return issues


def _check_constraint_consistency(cc: CoordComplex) -> List[CoordRepIssue]:
 """"""
    issues = []
    
    seen_pairs = set()
    for s1, s2 in cc.constraints.trans_pairs:
        pair_key = tuple(sorted([str(s1), str(s2)]))
        if pair_key in seen_pairs:
            issues.append(CoordRepIssue(
                code="DUPLICATE_CONSTRAINT",
                message=f"Duplicate trans pair: {s1}--{s2}",
                severity="warning"
            ))
        seen_pairs.add(pair_key)
    
    trans_set = set()
    for s1, s2 in cc.constraints.trans_pairs:
        trans_set.add(tuple(sorted([str(s1), str(s2)])))
    
    for s1, s2 in cc.constraints.cis_pairs:
        pair_key = tuple(sorted([str(s1), str(s2)]))
        if pair_key in trans_set:
            issues.append(CoordRepIssue(
                code="CONSTRAINT_CONTRADICTION",
                message=f"Pair is both trans and cis: {s1}--{s2}",
                severity="error"
            ))
    
    return issues


def validate_batch(complexes: List[CoordComplex]) -> dict:
    """
 Statistics
    
    Returns:
 Statistics
    """
    stats = {
        'total': len(complexes),
        'valid': 0,
        'with_warnings': 0,
        'with_errors': 0,
        'error_counts': {},
    }
    
    for cc in complexes:
        issues = validate_complex(cc)
        
        has_error = any(i.severity == 'error' for i in issues)
        has_warning = any(i.severity == 'warning' for i in issues)
        
        if not issues:
            stats['valid'] += 1
        elif has_error:
            stats['with_errors'] += 1
        elif has_warning:
            stats['with_warnings'] += 1
        
        for issue in issues:
            stats['error_counts'][issue.code] = stats['error_counts'].get(issue.code, 0) + 1
    
    return stats





