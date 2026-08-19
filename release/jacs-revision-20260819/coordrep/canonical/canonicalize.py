"""
Canonicalization engine.

Produces a unique, deterministic CoordRep string for any given
coordination complex regardless of input atom ordering:
    1. Sort ligands by a composite key
    2. Relabel ligand IDs (L1, L2, ...)
    3. Enumerate equivalent permutations and pick lexicographic minimum
    4. Canonicalize constraints
    5. Round shape vector
"""

from typing import List, Dict, Set
from copy import deepcopy
from itertools import permutations, product
from math import factorial
from collections import defaultdict

from ..core import CoordComplex, LigandModule, ConstraintSet, DonorSite


def canonicalize_complex(cc: CoordComplex) -> CoordComplex:
    """
    Canonicalize a CoordComplex.

    Steps:
        1. Standardize ligand SMILES (already canonical via RDKit)
        2. Sort ligands and relabel IDs
        3. Equivalent relabeling (permute duplicate ligands, pick minimum)
        4. Rewrite constraint symbols
        5. Round shape vector
        6. Mark as canonical
    """
    if cc.is_canonical:
        return cc

    cc = deepcopy(cc)

    # Step 1: Sort ligands
    sorted_ligands = sorted(cc.ligands, key=lambda lig: lig.get_sort_key())

    # Step 2: Relabel ligand IDs
    old_to_new_id = {}
    for i, lig in enumerate(sorted_ligands):
        old_id = lig.lig_id
        new_id = f"L{i + 1}"
        old_to_new_id[old_id] = new_id
        lig.lig_id = new_id

    cc.ligands = sorted_ligands

    # Step 3: Rewrite constraint symbols
    cc.constraints = _rewrite_constraints(cc.constraints, old_to_new_id)

    # Step 4: Equivalent relabeling (permute duplicates, pick lex-minimum)
    cc = _apply_equivalent_relabeling(cc)

    # Step 5: Round shape vector
    if cc.shape is not None:
        cc.shape.round(decimals=2)

    # Step 6: Mark as canonical
    cc.is_canonical = True

    return cc


def _rewrite_constraints(constraints: ConstraintSet,
                        id_map: Dict[str, str]) -> ConstraintSet:
    """Rewrite constraints with new ligand IDs."""
    new_trans = []
    for s1, s2 in constraints.trans_pairs:
        new_s1 = DonorSite(
            lig_id=id_map.get(s1.lig_id, s1.lig_id),
            donor_element=s1.donor_element,
            donor_rank=s1.donor_rank
        )
        new_s2 = DonorSite(
            lig_id=id_map.get(s2.lig_id, s2.lig_id),
            donor_element=s2.donor_element,
            donor_rank=s2.donor_rank
        )
        new_trans.append((new_s1, new_s2))

    new_cis = []
    for s1, s2 in constraints.cis_pairs:
        new_s1 = DonorSite(
            lig_id=id_map.get(s1.lig_id, s1.lig_id),
            donor_element=s1.donor_element,
            donor_rank=s1.donor_rank
        )
        new_s2 = DonorSite(
            lig_id=id_map.get(s2.lig_id, s2.lig_id),
            donor_element=s2.donor_element,
            donor_rank=s2.donor_rank
        )
        new_cis.append((new_s1, new_s2))

    return ConstraintSet(
        trans_pairs=new_trans,
        cis_pairs=new_cis,
        notes=constraints.notes
    ).canonicalize()


def _apply_equivalent_relabeling(cc: CoordComplex) -> CoordComplex:
    """
    Equivalent relabeling for duplicate ligands.

    Groups chemically equivalent ligand records and enumerates the Cartesian
    product of their within-group permutations.  Evaluating groups one at a
    time is not exact when the global lexicographic minimum requires two or
    more duplicate groups to be permuted simultaneously.

    The validated molecular core has CN <= 6, so the worst case is 6! = 720
    candidates.  Larger searches are rejected explicitly rather than silently
    truncating a duplicate group and returning an input-order-dependent string.
    """
    smiles_groups = defaultdict(list)
    for lig in cc.ligands:
        equivalence_key = (
            lig.payload_provenance,
            lig.smiles,
            tuple(sorted(zip(lig.donor_elements, lig.attach_atoms))),
            lig.dent,
            lig.eta,
            lig.charge,
            lig.connectivity_status,
        )
        smiles_groups[equivalence_key].append(lig.lig_id)

    duplicate_groups = [ids for ids in smiles_groups.values() if len(ids) > 1]

    if not duplicate_groups:
        return cc

    candidate_count = 1
    for group in duplicate_groups:
        candidate_count *= factorial(len(group))
    max_exact_candidates = 100_000
    if candidate_count > max_exact_candidates:
        raise ValueError(
            "exact equivalent-ligand canonicalization would require "
            f"{candidate_count} candidates (limit {max_exact_candidates}); "
            "the input lies outside the validated CN<=6 core"
        )

    best_cc = cc
    best_string = cc.to_string()
    permutation_families = [tuple(permutations(group)) for group in duplicate_groups]
    for permutation_choice in product(*permutation_families):
        perm_map = {}
        for group, perm in zip(duplicate_groups, permutation_choice):
            perm_map.update(dict(zip(group, perm)))
        candidate = _apply_permutation(cc, perm_map)
        candidate_string = candidate.to_string()
        if candidate_string < best_string:
            best_cc = candidate
            best_string = candidate_string

    return best_cc


def _apply_permutation(cc: CoordComplex, perm_map: Dict[str, str]) -> CoordComplex:
    """Apply a ligand-ID permutation."""
    cc = deepcopy(cc)

    for lig in cc.ligands:
        if lig.lig_id in perm_map:
            lig.lig_id = perm_map[lig.lig_id]

    cc.ligands = sorted(cc.ligands, key=lambda l: l.lig_id)

    new_trans = []
    for s1, s2 in cc.constraints.trans_pairs:
        new_s1 = DonorSite(
            lig_id=perm_map.get(s1.lig_id, s1.lig_id),
            donor_element=s1.donor_element,
            donor_rank=s1.donor_rank
        )
        new_s2 = DonorSite(
            lig_id=perm_map.get(s2.lig_id, s2.lig_id),
            donor_element=s2.donor_element,
            donor_rank=s2.donor_rank
        )
        new_trans.append((new_s1, new_s2))

    new_cis = []
    for s1, s2 in cc.constraints.cis_pairs:
        new_s1 = DonorSite(
            lig_id=perm_map.get(s1.lig_id, s1.lig_id),
            donor_element=s1.donor_element,
            donor_rank=s1.donor_rank
        )
        new_s2 = DonorSite(
            lig_id=perm_map.get(s2.lig_id, s2.lig_id),
            donor_element=s2.donor_element,
            donor_rank=s2.donor_rank
        )
        new_cis.append((new_s1, new_s2))

    cc.constraints = ConstraintSet(
        trans_pairs=new_trans,
        cis_pairs=new_cis,
        notes=cc.constraints.notes
    ).canonicalize()

    return cc


def verify_canonical_invariance(cc: CoordComplex,
                               other: CoordComplex) -> bool:
    """Verify that two CoordComplex objects produce the same canonical string."""
    cc_can = cc.canonicalize()
    other_can = other.canonicalize()

    return cc_can.to_string() == other_can.to_string()
