"""
tokenizer_variants.py
=====================
Factorized vs composite metal-block tokenization variants.

Default (after fix): FACTORIZED
  tokenizer.py _tokenize_metal now always emits [Fe] ;ox=+2 ;d=6 ;CN=6
  regardless of whether the input uses pipe or semicolon separators.

Legacy (composite):
  The entire metal block [Metal:Fe|ox:+2|d:d6|CN:6] was treated as one token
  due to a separator mismatch (serializer used | while tokenizer split on ;).
  Kept here only for ablation comparison.
"""

from __future__ import annotations

import re
from typing import List, Tuple


def retokenize_metal_factorized(coordrep_str: str) -> str:
    """
    Convert a CoordRep string from composite metal block format
    to factorized metal block format.

    [Metal:Fe|ox:+2|d:d6|CN:6] → [Fe;ox=+2;d=6;CN=6]

    The factorized format uses ';' as internal separator so
    _tokenize_metal splits it into separate tokens.
    """
    m = re.match(r'(\[Metal:([A-Z][a-z]?)\|([^\]]*)\])(.*)', coordrep_str, re.DOTALL)
    if not m:
        return coordrep_str

    metal = m.group(2)
    fields_str = m.group(3)  # e.g. "ox:+2|d:d6|CN:6"
    rest = m.group(4)

    # Parse fields
    parts = [f for f in fields_str.split('|') if f]
    factorized_parts = []
    for part in parts:
        if part.startswith('ox:'):
            val = part[3:]  # +2
            factorized_parts.append(f';ox={val}')
        elif part.startswith('d:d'):
            val = part[2:]  # d6
            factorized_parts.append(f';d={val[1:]}')  # ;d=6
        elif part.startswith('d:'):
            val = part[2:]
            factorized_parts.append(f';d={val}')
        elif part.startswith('CN:'):
            val = part[3:]
            factorized_parts.append(f';CN={val}')

    factorized_metal = f'[{metal}{"".join(factorized_parts)}]'
    return factorized_metal + rest


def retokenize_metal_composite(coordrep_str: str) -> str:
    """
    Ensure a CoordRep string uses composite metal block format.
    If already composite, return as-is.
    If factorized [Fe;ox=+2;d=6;CN=6], convert to [Metal:Fe|ox:+2|d:d6|CN:6].
    """
    # Check if already composite
    if coordrep_str.startswith('[Metal:'):
        return coordrep_str

    m = re.match(r'\[([A-Z][a-z]?)([^\]]*)\](.*)', coordrep_str, re.DOTALL)
    if not m:
        return coordrep_str

    metal = m.group(1)
    fields_str = m.group(2)  # e.g. ";ox=+2;d=6;CN=6"
    rest = m.group(3)

    parts = [f for f in fields_str.split(';') if f]
    composite_parts = []
    for part in parts:
        if part.startswith('ox='):
            composite_parts.append(f'ox:{part[3:]}')
        elif part.startswith('d='):
            composite_parts.append(f'd:d{part[2:]}')
        elif part.startswith('CN='):
            composite_parts.append(f'CN:{part[3:]}')
        elif part.startswith('row='):
            composite_parts.append(f'row:{part[4:]}')

    fields = '|'.join([f'Metal:{metal}'] + composite_parts)
    return f'[{fields}]' + rest


def tokenize_metal_block_factorized(block: str) -> List[str]:
    """
    Tokenize a factorized metal block into separate tokens.
    Input: [Fe;ox=+2;d=6;CN=6]
    Output: ['[Fe]', ';ox=+2', ';d=6', ';CN=6']
    """
    inner = block[1:-1]  # Fe;ox=+2;d=6;CN=6
    parts = inner.split(';')
    tokens = [f'[{parts[0]}]']
    for part in parts[1:]:
        if part:
            tokens.append(f';{part}')
    return tokens


def tokenize_metal_block_composite(block: str) -> List[str]:
    """
    Tokenize a composite metal block as a single token.
    Input: [Metal:Fe|ox:+2|d:d6|CN:6]
    Output: ['[Metal:Fe|ox:+2|d:d6|CN:6]']
    """
    return [block]


def get_factorized_tokens(coordrep_str: str) -> List[str]:
    """
    Re-tokenize a composite CoordRep string using factorized metal tokens.
    Returns full token list.
    """
    factorized = retokenize_metal_factorized(coordrep_str)

    # Now manually tokenize using the factorized scheme
    # The factorized string has [Fe;ox=+2;d=6;CN=6] which
    # the standard _tokenize_metal will correctly split on ';'
    return factorized
