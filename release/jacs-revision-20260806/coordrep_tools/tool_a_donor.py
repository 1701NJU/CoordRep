#!/usr/bin/env python3
"""
Tool A: Donor (Design Assistant)

        " + " Top-k donor atom
"""

import json
import re
from typing import List, Tuple, Dict, Optional
from pathlib import Path
from collections import defaultdict

from .infer import mlm_topk, load_model_and_tokenizer
from .baselines import RandomBaseline, FrequencyBaseline, ConditionalFrequencyBaseline, DONOR_ATOMS


class DonorPredictor:
    """
        Donor atom Predict
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

        # Preload model
        load_model_and_tokenizer(checkpoint_path, tokenizer_path, device)

    def predict(self, coordrep_masked: str, mask_positions: List[int], k: int = 5) -> List[List[Tuple[str, float]]]:
        """
        Predict donor atoms

        Args:
        coordrep_masked: [MASK]
        mask_positions: mask
            k: top-k

        Returns:
        mask top-k Predict
        """
        return mlm_topk(
            coordrep_masked,
            mask_positions,
            k,
            self.checkpoint_path,
            self.tokenizer_path,
            self.device
        )

    def predict_single(self, coordrep_masked: str, mask_position: int, k: int = 5) -> List[Tuple[str, float]]:
        """Predict single position"""
        results = self.predict(coordrep_masked, [mask_position], k)
        return results[0] if results else []

    def filter_donor_predictions(self, predictions: List[Tuple[str, float]],
                                  allowed_donors: List[str] = None) -> List[Tuple[str, float]]:
        """
        Filter predictions, keep only donor atoms

        Args:
        predictions: Predict
        allowed_donors: donor

        Returns:
        Predict
        """
        if allowed_donors is None:
            allowed_donors = DONOR_ATOMS

        filtered = [(tok, prob) for tok, prob in predictions if tok in allowed_donors]

        # Renormalize
        if filtered:
            total = sum(p for _, p in filtered)
            if total > 0:
                filtered = [(tok, prob / total) for tok, prob in filtered]

        return filtered


def evaluate_donor_prediction(
    samples: List[dict],
    checkpoint_path: str,
    tokenizer_path: str = None,
    device: str = "cuda",
    k: int = 5
) -> Dict:
    """
    Evaluate donor prediction performance

    Args:
        samples: 'coordrep_masked', 'mask_positions', 'target_token', 'meta'
        checkpoint_path: checkpoint
        tokenizer_path: tokenizer
        device:
        k: top-k

    Returns:
    """
    predictor = DonorPredictor(checkpoint_path, tokenizer_path, device)

    results = {
        'total': 0,
        'correct_top1': 0,
        'correct_top5': 0,
        'by_donor': defaultdict(lambda: {'total': 0, 'correct_top1': 0, 'correct_top5': 0}),
        'by_cn': defaultdict(lambda: {'total': 0, 'correct_top1': 0, 'correct_top5': 0}),
        'by_metal': defaultdict(lambda: {'total': 0, 'correct_top1': 0, 'correct_top5': 0}),
        'predictions': [],
    }

    for sample in samples:
        coordrep_masked = sample['coordrep_masked']
        mask_positions = sample['mask_positions']
        target_token = sample['target_token']
        meta = sample.get('meta', {})

        # Predict
        preds = predictor.predict(coordrep_masked, mask_positions, k)

        if not preds or not preds[0]:
            continue

        pred_topk = preds[0]
        pred_tokens = [tok for tok, prob in pred_topk]

        # Statistics
        results['total'] += 1

        is_top1 = pred_tokens[0] == target_token if pred_tokens else False
        is_top5 = target_token in pred_tokens

        if is_top1:
            results['correct_top1'] += 1
        if is_top5:
            results['correct_top5'] += 1

        # Stratify by donor
        donor_key = target_token
        results['by_donor'][donor_key]['total'] += 1
        if is_top1:
            results['by_donor'][donor_key]['correct_top1'] += 1
        if is_top5:
            results['by_donor'][donor_key]['correct_top5'] += 1

        # Stratify by CN
        cn = meta.get('cn')
        if cn is not None:
            results['by_cn'][cn]['total'] += 1
            if is_top1:
                results['by_cn'][cn]['correct_top1'] += 1
            if is_top5:
                results['by_cn'][cn]['correct_top5'] += 1

        # Stratify by metal
        metal = meta.get('metal')
        if metal:
            results['by_metal'][metal]['total'] += 1
            if is_top1:
                results['by_metal'][metal]['correct_top1'] += 1
            if is_top5:
                results['by_metal'][metal]['correct_top5'] += 1

        # Record prediction details (for case cards)
        results['predictions'].append({
            'id': sample.get('id', ''),
            'target': target_token,
            'pred_top5': pred_topk[:5],
            'correct_top1': is_top1,
            'correct_top5': is_top5,
            'rank': pred_tokens.index(target_token) + 1 if target_token in pred_tokens else -1,
            'meta': meta,
        })

    # Compute accuracy
    n = results['total']
    results['top1_acc'] = results['correct_top1'] / n if n > 0 else 0
    results['top5_acc'] = results['correct_top5'] / n if n > 0 else 0

    return results


def compare_with_baselines(
    samples: List[dict],
    checkpoint_path: str,
    tokenizer_path: str = None,
    device: str = "cuda",
    training_samples: List[dict] = None
) -> Dict:
    """
    Compare with baselines

    Args:
        samples:
        checkpoint_path: checkpoint
        tokenizer_path: tokenizer
        device:
        training_samples: frequency baselines

    Returns:
    """
    model_results = evaluate_donor_prediction(samples, checkpoint_path, tokenizer_path, device)

    # Random baseline
    random_bl = RandomBaseline(candidate_set=DONOR_ATOMS)
    random_results = random_bl.get_accuracy(samples)

    # Global frequency baseline
    freq_bl = FrequencyBaseline()
    if training_samples:
        freq_bl.fit(training_samples)
    else:
        freq_bl.fit(samples) #
    freq_results = freq_bl.get_accuracy(samples)

    # Conditional frequency baseline
    cond_freq_bl = ConditionalFrequencyBaseline()
    if training_samples:
        cond_freq_bl.fit(training_samples)
    else:
        cond_freq_bl.fit(samples)
    cond_freq_results = cond_freq_bl.get_accuracy(samples)

    return {
        'model': {
            'top1_acc': model_results['top1_acc'],
            'top5_acc': model_results['top5_acc'],
            'n': model_results['total'],
        },
        'random': random_results,
        'global_frequency': freq_results,
        'conditional_frequency': cond_freq_results,
    }


def generate_casecards(
    predictions: List[dict],
    n_cards: int = 10,
    include_failures: bool = True
) -> List[dict]:
    """
    Generate case cards for visualization

    Args:
        predictions: Predict
        n_cards:
        include_failures: Predict

    Returns:
        Case cards
    """
    # Select representative examples
    successes = [p for p in predictions if p['correct_top1']]
    failures = [p for p in predictions if not p['correct_top5']]
    partial = [p for p in predictions if p['correct_top5'] and not p['correct_top1']]

    cards = []

    if successes:
        successes_sorted = sorted(successes, key=lambda x: x['pred_top5'][0][1] if x['pred_top5'] else 0, reverse=True)
        cards.extend(successes_sorted[:n_cards // 3])

    if partial:
        cards.extend(partial[:n_cards // 3])

    if include_failures and failures:
        cards.extend(failures[:n_cards // 3])

    # Format
    formatted_cards = []
    for card in cards[:n_cards]:
        formatted_cards.append({
            'id': card['id'],
            'metal': card['meta'].get('metal', 'N/A'),
            'cn': card['meta'].get('cn', 'N/A'),
            'geom': card['meta'].get('geom', 'N/A'),
            'target_token': card['target'],
            'predictions': [
                {'token': tok, 'prob': round(prob, 4)}
                for tok, prob in card['pred_top5'][:5]
            ],
            'target_rank': card['rank'],
            'status': 'correct_top1' if card['correct_top1'] else ('correct_top5' if card['correct_top5'] else 'failed'),
        })

    return formatted_cards
