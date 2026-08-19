"""
    (Coordinated Masking)

""

1. Block-wise Masking: Mask Shape Constraint
2. Random Masking: 15% Mask
3. Span Masking: Mask Token
"""

import random
import numpy as np
from typing import List, Tuple, Dict
from dataclasses import dataclass


@dataclass
class MaskingConfig:
    """"""
    mlm_probability: float = 0.15

    block_mask_probability: float = 0.30

    span_mask_probability: float = 0.20
    max_span_length: int = 5

    protected_tokens: List[str] = None

    def __post_init__(self):
        if self.protected_tokens is None:
            self.protected_tokens = ["[CLS]", "[SEP]", "[PAD]"]


class CoordinatedMasker:
    """

    1. 30% Mask <Shape> {Constraint}
    → [Metal] |Ligands|
    2. 20% Span Masking
    →
    3. 15% MLM
    """

    def __init__(self, config: MaskingConfig = None):
        self.config = config or MaskingConfig()

    def create_mlm_inputs(
        self,
        input_ids: List[int],
        tokens: List[str],
        mask_token_id: int
    ) -> Tuple[List[int], List[int]]:
        """
        MLM

        Returns:
        masked_ids:
    labels: Token-100
        """
        masked_ids = input_ids.copy()
        labels = [-100] * len(input_ids) # -100

        if random.random() < self.config.block_mask_probability:
            masked_ids, labels = self._block_mask(
                input_ids, tokens, mask_token_id
            )
        elif random.random() < self.config.span_mask_probability:
            masked_ids, labels = self._span_mask(
                input_ids, tokens, mask_token_id
            )
        else:
            masked_ids, labels = self._random_mask(
                input_ids, tokens, mask_token_id
            )

        return masked_ids, labels

    def _block_mask(
        self,
        input_ids: List[int],
        tokens: List[str],
        mask_token_id: int
    ) -> Tuple[List[int], List[int]]:
        """
        Block-wise Masking: Mask Shape Constraint

    [Metal] |Ligands|
        """
        masked_ids = input_ids.copy()
        labels = [-100] * len(input_ids)

        shape_start = None
        shape_end = None
        for i, tok in enumerate(tokens):
            if tok == "<Shape:":
                shape_start = i
            elif tok == ">" and shape_start is not None:
                shape_end = i + 1
                break

        constraint_ranges = []
        in_constraint = False
        start = None
        for i, tok in enumerate(tokens):
            if tok in ["{trans:", "{cis:"]:
                in_constraint = True
                start = i
            elif tok == "}" and in_constraint:
                constraint_ranges.append((start, i + 1))
                in_constraint = False

        if shape_start is not None and shape_end is not None:
            for i in range(shape_start, shape_end):
                if tokens[i] not in self.config.protected_tokens:
                    labels[i] = input_ids[i]
                    masked_ids[i] = mask_token_id

        for start, end in constraint_ranges:
            for i in range(start, end):
                if tokens[i] not in self.config.protected_tokens:
                    labels[i] = input_ids[i]
                    masked_ids[i] = mask_token_id

        return masked_ids, labels

    def _span_mask(
        self,
        input_ids: List[int],
        tokens: List[str],
        mask_token_id: int
    ) -> Tuple[List[int], List[int]]:
        """
        Span Masking: Mask Token
        """
        masked_ids = input_ids.copy()
        labels = [-100] * len(input_ids)

        valid_positions = [
            i for i, tok in enumerate(tokens)
            if tok not in self.config.protected_tokens
        ]

        if not valid_positions:
            return masked_ids, labels

        num_to_mask = max(1, int(len(valid_positions) * self.config.mlm_probability))

        masked = set()
        while len(masked) < num_to_mask:
            start = random.choice(valid_positions)
            span_len = random.randint(1, self.config.max_span_length)

            for i in range(start, min(start + span_len, len(tokens))):
                if i in valid_positions and i not in masked:
                    masked.add(i)

        for i in masked:
            labels[i] = input_ids[i]
            masked_ids[i] = mask_token_id

        return masked_ids, labels

    def _random_mask(
        self,
        input_ids: List[int],
        tokens: List[str],
        mask_token_id: int
    ) -> Tuple[List[int], List[int]]:
        """
        Random Masking (15%)
        """
        masked_ids = input_ids.copy()
        labels = [-100] * len(input_ids)

        for i, tok in enumerate(tokens):
            if tok in self.config.protected_tokens:
                continue

            if random.random() < self.config.mlm_probability:
                labels[i] = input_ids[i]

                # 80% MASK, 10% random, 10% original
                r = random.random()
                if r < 0.8:
                    masked_ids[i] = mask_token_id
                elif r < 0.9:
                    masked_ids[i] = random.randint(5, 100)  # random token
                # else: keep original

        return masked_ids, labels


class LigandShuffler:
    """
    (Ligand Shuffle Augmentation)

    1. N!
    2.
    3.
    """

    def shuffle(self, coordrep_string: str) -> str:
        """

        : [Fe...]<Shape...>{trans:L1:N:1--L2:N:1}(L1@1,2)(L2@3,4)|L1=NCCN||L2=Cl|
    : [Fe...]<Shape...>{trans:L2:N:1--L1:N:1}(L2@1,2)(L1@3,4)|L2=NCCN||L1=Cl|
        """
        import re

        ligand_pattern = re.compile(r'\|L(\d+)=([^|]+)\|')
        ligands = ligand_pattern.findall(coordrep_string)

        if len(ligands) < 2:
            return coordrep_string

        old_ids = [f"L{lid}" for lid, _ in ligands]
        new_ids = old_ids.copy()
        random.shuffle(new_ids)

        id_mapping = {old: new for old, new in zip(old_ids, new_ids)}

        result = coordrep_string

        temp_mapping = {old: f"__TEMP_{i}__" for i, old in enumerate(old_ids)}
        for old, temp in temp_mapping.items():
            result = result.replace(old, temp)

        for temp, new in zip(temp_mapping.values(), new_ids):
            result = result.replace(temp, new)

        return result

    def augment_batch(self, strings: List[str], n_augments: int = 2) -> List[str]:
        """
        """
        augmented = []
        for s in strings:
            augmented.append(s) #
            for _ in range(n_augments):
                augmented.append(self.shuffle(s))
        return augmented


if __name__ == "__main__":
    print("=== Coordinated Masking Test ===")

    masker = CoordinatedMasker()

    tokens = ["[CLS]", "[Fe]", ";ox=+2", "<Shape:", "Oh", "=", "V_025", ">",
              "{trans:", "L1:N:1", "--", "L2:N:1", "}", "[SEP]"]
    input_ids = list(range(len(tokens)))

    masked_ids, labels = masker._block_mask(input_ids, tokens, mask_token_id=4)

    print("Tokens:", tokens)
    print("Original IDs:", input_ids)
    print("Masked IDs:", masked_ids)
    print("Labels:", labels)
    print()

    print("=== Ligand Shuffle Test ===")

    shuffler = LigandShuffler()

    test_input = "[Fe;ox=+2]{trans:L1:N:1--L2:N:1}(L1@1,2)(L2@3,4)|L1=NCCN||L2=Cl|"
    print("Original:", test_input)

    for i in range(3):
        shuffled = shuffler.shuffle(test_input)
        print(f"Shuffle {i+1}:", shuffled)
