#!/usr/bin/env python3
"""
Tool B: / (Spellchecker)

/ CoordRep
"""

import json
import re
from typing import List, Tuple, Dict, Optional
from pathlib import Path
from collections import defaultdict

from .infer import mlm_topk, load_model_and_tokenizer, get_token_probs_at_position
from .validate import is_valid_coordrep, validate_brackets, classify_error_type


STRUCTURE_TOKENS = ['(', ')', '[', ']', '{', '}', '<', '>', '|', ';', ',', '+', '-', '=', ':']


class CoordRepRepairer:
    """
 CoordRep
    """

    def __init__(self, checkpoint_path: str, tokenizer_path: str = None, device: str = "cuda"):
        """
        Args:
 checkpoint_path: checkpoint
 tokenizer_path: tokenizer
 device:
        """
        self.checkpoint_path = checkpoint_path
        self.tokenizer_path = tokenizer_path
        self.device = device

        # Preload model and tokenizer
        self.model, self.tokenizer = load_model_and_tokenizer(checkpoint_path, tokenizer_path, device)

    def tokenize(self, s: str) -> List[str]:
        """Tokenize the string"""
        return self.tokenizer.tokenize(s)

    def detokenize(self, tokens: List[str]) -> str:
        """Merge tokens back into string"""
        # Simple concatenation (CoordRep typically has no spaces)
        return ''.join(tokens)

    def find_candidate_positions(self, tokens: List[str]) -> List[int]:
        """
        Identify positions that may need repair

        Prioritize structure-related tokens
        """
        candidates = []

        for i, tok in enumerate(tokens):
            if tok in STRUCTURE_TOKENS:
                candidates.append(i)
            elif tok.startswith('[UNK'):
                candidates.append(i)

        return candidates

    def try_repair_position(self, tokens: List[str], position: int, k: int = 5) -> Optional[str]:
        """
 Attempt repair

        Args:
 tokens: token
 position:
 k:

        Returns:
 token
        """
        # Create masked version
        masked_tokens = tokens.copy()
        original_token = masked_tokens[position]
        masked_tokens[position] = '[MASK]'

        # Get predictions
        masked_str = self.detokenize(masked_tokens)
        preds = mlm_topk(masked_str, [position], k, self.checkpoint_path, self.tokenizer_path, self.device)

        if not preds or not preds[0]:
            return None

        # Check if each candidate makes the sequence valid
        for tok, prob in preds[0]:
            test_tokens = tokens.copy()
            test_tokens[position] = tok
            test_str = self.detokenize(test_tokens)

            if is_valid_coordrep(test_str):
                return tok

        return None

    def repair(self, coordrep_corrupted: str, max_iters: int = 5, k: int = 5) -> Dict:
        """
        Repair CoordRep sequence

        Args:
 coordrep_corrupted:
 max_iters:
 k:

        Returns:
            {
                'repaired': str,
                'was_valid_before': bool,
                'is_valid_after': bool,
                'edits': [{'pos': int, 'from': str, 'to': str, 'prob': float}],
                'failure_reason': str | None
            }
        """
        result = {
            'repaired': coordrep_corrupted,
            'was_valid_before': is_valid_coordrep(coordrep_corrupted),
            'is_valid_after': False,
            'edits': [],
            'failure_reason': None,
        }

        # Return directly if already valid
        if result['was_valid_before']:
            result['is_valid_after'] = True
            return result

        # Tokenize
        try:
            tokens = self.tokenize(coordrep_corrupted)
        except Exception as e:
            result['failure_reason'] = f'tokenization_error: {str(e)}'
            return result

        # Identify candidate positions
        candidates = self.find_candidate_positions(tokens)

        if not candidates:
            # No obvious candidates, scan all positions
            candidates = list(range(len(tokens)))

        # Iterative repair
        current_tokens = tokens.copy()
        edits = []

        for iteration in range(max_iters):
            # Check if repaired
            current_str = self.detokenize(current_tokens)
            if is_valid_coordrep(current_str):
                result['repaired'] = current_str
                result['is_valid_after'] = True
                result['edits'] = edits
                return result

            # Attempt repair at each candidate position
            repaired_this_iter = False

            for pos in candidates:
                if pos >= len(current_tokens):
                    continue

                # Attempt substitution
                masked_tokens = current_tokens.copy()
                original = masked_tokens[pos]
                masked_tokens[pos] = '[MASK]'

                masked_str = self.detokenize(masked_tokens)
                preds = mlm_topk(masked_str, [pos], k, self.checkpoint_path, self.tokenizer_path, self.device)

                if not preds or not preds[0]:
                    continue

                # Check each candidate
                for tok, prob in preds[0]:
                    if tok == original:
                        continue # token

                    test_tokens = current_tokens.copy()
                    test_tokens[pos] = tok
                    test_str = self.detokenize(test_tokens)

                    # Check for improvement (closer to valid or already valid)
                    if is_valid_coordrep(test_str):
                        current_tokens = test_tokens
                        edits.append({
                            'pos': pos,
                            'from': original,
                            'to': tok,
                            'prob': prob,
                        })
                        repaired_this_iter = True
                        break

            if not repaired_this_iter:
                break

        # Final check
        result['repaired'] = self.detokenize(current_tokens)
        result['is_valid_after'] = is_valid_coordrep(result['repaired'])
        result['edits'] = edits

        if not result['is_valid_after']:
            result['failure_reason'] = 'max_iterations_reached'

        return result


def repair_coordrep(
    coordrep_corrupted: str,
    checkpoint_path: str,
    tokenizer_path: str = None,
    max_iters: int = 5,
    k: int = 5
) -> Dict:
    """
 Repair CoordRep sequence
    """
    repairer = CoordRepRepairer(checkpoint_path, tokenizer_path)
    return repairer.repair(coordrep_corrupted, max_iters, k)


def evaluate_structure_recovery(
    samples: List[dict],
    checkpoint_path: str,
    tokenizer_path: str = None,
    device: str = "cuda",
    k: int = 5
) -> Dict:
    """
    Evaluate structure token recovery performance

    Args:
 samples:
 checkpoint_path: checkpoint
 tokenizer_path: tokenizer
 device:
        k: top-k

    Returns:
    """
    results = {
        'total_positions': 0,
        'correct_top1': 0,
        'correct_top5': 0,
        'by_cn': defaultdict(lambda: {'total': 0, 'correct_top1': 0, 'correct_top5': 0}),
        'predictions': [],
    }

    for sample in samples:
        coordrep_masked = sample['coordrep_masked']
        mask_positions = sample['mask_positions']
        target_tokens = sample.get('target_tokens', [sample.get('target_token')])
        meta = sample.get('meta', {})

        # Ensure target_tokens is a list
        if not isinstance(target_tokens, list):
            target_tokens = [target_tokens]

        # Predict
        preds = mlm_topk(coordrep_masked, mask_positions, k, checkpoint_path, tokenizer_path, device)

        for i, (pos, target) in enumerate(zip(mask_positions, target_tokens)):
            if i >= len(preds):
                continue

            pred_topk = preds[i]
            pred_tokens = [tok for tok, prob in pred_topk]

            results['total_positions'] += 1

            is_top1 = pred_tokens[0] == target if pred_tokens else False
            is_top5 = target in pred_tokens

            if is_top1:
                results['correct_top1'] += 1
            if is_top5:
                results['correct_top5'] += 1

            # Stratify by CN
            cn = meta.get('cn')
            if cn is not None:
                results['by_cn'][cn]['total'] += 1
                if is_top1:
                    results['by_cn'][cn]['correct_top1'] += 1
                if is_top5:
                    results['by_cn'][cn]['correct_top5'] += 1

    # Compute accuracy
    n = results['total_positions']
    results['top1_acc'] = results['correct_top1'] / n if n > 0 else 0
    results['top5_acc'] = results['correct_top5'] / n if n > 0 else 0

    return results


def evaluate_repair_on_synthetic(
    corrupted_samples: List[dict],
    checkpoint_path: str,
    tokenizer_path: str = None,
    max_iters: int = 5,
    k: int = 5
) -> Dict:
    """
 Evaluate repair

    Args:
 corrupted_samples: 'coordrep_clean' 'coordrep_corrupted'
 checkpoint_path: checkpoint
 tokenizer_path: tokenizer
 max_iters:
        k: top-k

    Returns:
    """
    repairer = CoordRepRepairer(checkpoint_path, tokenizer_path)

    results = {
        'total': 0,
        'valid_before': 0,
        'valid_after': 0,
 'exact_match': 0, #
        'by_corruption_type': defaultdict(lambda: {
            'total': 0,
            'valid_before': 0,
            'valid_after': 0,
            'exact_match': 0
        }),
        'error_types': defaultdict(int),
        'repairs': [],
    }

    for sample in corrupted_samples:
        coordrep_clean = sample['coordrep_clean']
        coordrep_corrupted = sample['coordrep_corrupted']
        corruption_type = sample.get('corruption_type', 'unknown')

        results['total'] += 1

        repair_result = repairer.repair(coordrep_corrupted, max_iters, k)

        # Statistics
        if repair_result['was_valid_before']:
            results['valid_before'] += 1
            results['by_corruption_type'][corruption_type]['valid_before'] += 1

        if repair_result['is_valid_after']:
            results['valid_after'] += 1
            results['by_corruption_type'][corruption_type]['valid_after'] += 1

        if repair_result['repaired'] == coordrep_clean:
            results['exact_match'] += 1
            results['by_corruption_type'][corruption_type]['exact_match'] += 1

        results['by_corruption_type'][corruption_type]['total'] += 1

        if repair_result['failure_reason']:
            results['error_types'][repair_result['failure_reason']] += 1
        elif not repair_result['is_valid_after']:
            error_type = classify_error_type(repair_result['repaired'])
            results['error_types'][error_type or 'unknown'] += 1

        results['repairs'].append({
            'id': sample.get('id', ''),
            'corruption_type': corruption_type,
            'was_valid_before': repair_result['was_valid_before'],
            'is_valid_after': repair_result['is_valid_after'],
            'exact_match': repair_result['repaired'] == coordrep_clean,
            'n_edits': len(repair_result['edits']),
        })

    # Compute ratios
    n = results['total']
    results['valid_rate_before'] = results['valid_before'] / n if n > 0 else 0
    results['valid_rate_after'] = results['valid_after'] / n if n > 0 else 0
    results['exact_match_rate'] = results['exact_match'] / n if n > 0 else 0
    results['uplift'] = results['valid_rate_after'] - results['valid_rate_before']

    return results
