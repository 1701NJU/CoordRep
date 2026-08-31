#!/usr/bin/env python3
"""
CoordRep

 is_valid_coordrep CoordRep
"""

import re
from typing import Tuple, List, Optional


def validate_brackets(s: str) -> Tuple[bool, Optional[str]]:
    """
    Validate bracket matching

    Returns:
        (is_valid, error_message)
    """
    stack = []
    pairs = {'(': ')', '[': ']', '{': '}', '<': '>'}

    for i, char in enumerate(s):
        if char in pairs:
            stack.append((char, i))
        elif char in pairs.values():
            if not stack:
                return False, f"Unmatched closing bracket '{char}' at position {i}"
            open_char, _ = stack.pop()
            if pairs[open_char] != char:
                return False, f"Mismatched brackets: expected '{pairs[open_char]}', got '{char}' at position {i}"

    if stack:
        unmatched = stack[0]
        return False, f"Unmatched opening bracket '{unmatched[0]}' at position {unmatched[1]}"

    return True, None


def validate_metal_token(s: str) -> Tuple[bool, Optional[str]]:
    """
 Metal token

 Valid format: [Metal:La|CN:6] [Metal:Fe|ox:+2|d:d6|CN:8]
    """
    metal_pattern = r'\[Metal:[A-Z][a-z]?\|(?:ox:[^\|]+\|)?(?:d:[^\|]+\|)?CN:\d+\]'
    if not re.search(metal_pattern, s):
        # May not have metal token; can still be valid depending on context
        return True, None
    return True, None


def validate_shape_token(s: str) -> Tuple[bool, Optional[str]]:
    """
 Shape token

    Valid format: <ShapeBest:Oh|Class:dist|Delta:3|V:4.99,12.35>
    """
    shape_pattern = r'<Shape[^>]*>'
    matches = re.findall(shape_pattern, s)
    # Shape token Optional
    return True, None


def validate_ligand_tokens(s: str) -> Tuple[bool, Optional[str]]:
    """
 Ligand token

    Valid format: |L1=SMILES|L2=SMILES|...
    """
    # Ligand Starts with |L\d+=
    ligand_pattern = r'\|L\d+='
    matches = re.findall(ligand_pattern, s)

    # Check if ligand indices are consecutive
    if matches:
        numbers = sorted(int(re.search(r'\d+', m).group()) for m in matches)
        expected = list(range(1, len(numbers) + 1))
        # Allow non-starting from 1, but must be consecutive
        if numbers != list(range(numbers[0], numbers[-1] + 1)):
            return False, f"Non-consecutive ligand indices: {numbers}"

    return True, None


def validate_stereo_constraints(s: str) -> Tuple[bool, Optional[str]]:
    """

 Valid format: {trans:L1:N:1--L2:O:1} {cis:...}
    """
    stereo_pattern = r'\{(trans|cis|mer|fac):[^}]+\}'
    matches = re.findall(stereo_pattern, s)
    # Stereo constraints are optional
    return True, None


def validate_charge_tokens(s: str) -> Tuple[bool, Optional[str]]:
    """
 token
    """
    # Perform basic validation here
    invalid_charge = re.findall(r'[\+\-]\d{3,}', s)  # Should not have values like +100
    if invalid_charge:
        return False, f"Invalid charge format: {invalid_charge}"
    return True, None


def validate_structure(s: str) -> Tuple[bool, List[str]]:
    """
    Full structural validation

    Returns:
        (is_valid, list_of_errors)
    """
    errors = []

    valid, err = validate_brackets(s)
    if not valid:
        errors.append(err)

    # 2. Metal token
    valid, err = validate_metal_token(s)
    if not valid:
        errors.append(err)

    # 3. Shape token
    valid, err = validate_shape_token(s)
    if not valid:
        errors.append(err)

    # 4. Ligand tokens
    valid, err = validate_ligand_tokens(s)
    if not valid:
        errors.append(err)

    # 5. Stereo constraints
    valid, err = validate_stereo_constraints(s)
    if not valid:
        errors.append(err)

    # 6. Charge tokens
    valid, err = validate_charge_tokens(s)
    if not valid:
        errors.append(err)

    return len(errors) == 0, errors


def is_valid_coordrep(s: str, strict: bool = False) -> bool:
    """
 CoordRep

    Args:
 s: CoordRep
 strict:

    Returns:
 bool:
    """
    if not s or not isinstance(s, str):
        return False

    # Basic checks
    s = s.strip()
    if len(s) < 10:
        return False

    # Bracket matching (most important)
    valid, _ = validate_brackets(s)
    if not valid:
        return False

    # Must contain Metal token or at least resemble CoordRep
    has_metal = '[Metal:' in s or re.search(r'\[[A-Z][a-z]?\]', s)
    has_ligand = '|L' in s or re.search(r'\|L\d+=', s)

    if not has_metal and not has_ligand:
        return False

    if strict:
        valid, errors = validate_structure(s)
        return valid

    return True


def get_validation_errors(s: str) -> List[str]:
    """
    Get all validation errors

    Args:
 s: CoordRep

    Returns:
    """
    if not s or not isinstance(s, str):
        return ["Empty or invalid input"]

    _, errors = validate_structure(s)
    return errors


def classify_error_type(s: str) -> Optional[str]:
    """
    Classify error type

    Returns:
 : 'bracket', 'charge', 'delimiter', 'truncation', 'other', None (if valid)
    """
    if is_valid_coordrep(s):
        return None

    # Check brackets
    valid, err = validate_brackets(s)
    if not valid:
        return 'bracket'

    # Check for truncation (incomplete ending)
    if s.endswith('|') or s.endswith(':') or s.endswith('='):
        return 'truncation'

    # Check delimiters
    if '||' in s or ';;' in s or ',,' in s:
        return 'delimiter'

    # Check charges
    valid, err = validate_charge_tokens(s)
    if not valid:
        return 'charge'

    return 'other'
