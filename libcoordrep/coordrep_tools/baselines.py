#!/usr/bin/env python3
"""
Baseline

 RandomFrequencyConditional Frequency baseline
"""

import random
import numpy as np
from typing import List, Tuple, Dict, Optional
from collections import Counter, defaultdict


# Standard donor atom set
DONOR_ATOMS = ['C', 'N', 'O', 'S', 'P', 'F', 'Cl', 'Br', 'I', 'Se', 'Te', 'As', 'Si', 'B', 'H']

# Standard structure tokens
STRUCTURE_TOKENS = ['L', 'Td', 'SP', 'Oh', 'TP', 'TBP', 'SPY', 'TPr', '(', ')', '[', ']', '{', '}', '|', ';', ',', '+', '-']


class RandomBaseline:
    """
 baseline
    """
    
    def __init__(self, candidate_set: List[str] = None, seed: int = 42):
        """
        Args:
 candidate_set: token DONOR_ATOMS
 seed:
        """
        self.candidate_set = candidate_set or DONOR_ATOMS
        self.seed = seed
        random.seed(seed)
        np.random.seed(seed)
    
    def predict(self, context: dict = None, k: int = 5) -> List[Tuple[str, float]]:
        """
 Predict top-k
        
        Args:
 context: baseline
            k: top-k
        
        Returns:
            [(token, prob), ...]
        """
        # Uniform probability
        prob = 1.0 / len(self.candidate_set)
        
        # Random ordering
        shuffled = self.candidate_set.copy()
        random.shuffle(shuffled)
        
        return [(tok, prob) for tok in shuffled[:k]]
    
    def predict_batch(self, contexts: List[dict], k: int = 5) -> List[List[Tuple[str, float]]]:
 """Predict"""
        return [self.predict(ctx, k) for ctx in contexts]
    
    def get_accuracy(self, samples: List[dict], target_field: str = 'target_token') -> Dict:
        """
        Compute accuracy
        
        Args:
 samples:
 target_field:
        
        Returns:
            {'top1_acc': float, 'top5_acc': float, 'n': int}
        """
        correct_top1 = 0
        correct_top5 = 0
        n = len(samples)
        
        for sample in samples:
            target = sample[target_field]
            pred = self.predict(sample, k=5)
            pred_tokens = [t for t, p in pred]
            
            if pred_tokens[0] == target:
                correct_top1 += 1
            if target in pred_tokens:
                correct_top5 += 1
        
        return {
            'top1_acc': correct_top1 / n if n > 0 else 0,
            'top5_acc': correct_top5 / n if n > 0 else 0,
            'n': n,
            'baseline_type': 'random'
        }


class FrequencyBaseline:
    """
 baseline token
    """
    
    def __init__(self, training_samples: List[dict] = None, target_field: str = 'target_token'):
        """
        Args:
 training_samples:
 target_field:
        """
        self.target_field = target_field
        self.freq = Counter()
        
        if training_samples:
            self.fit(training_samples)
    
    def fit(self, samples: List[dict]):
 """"""
        self.freq = Counter()
        for sample in samples:
            if self.target_field in sample:
                self.freq[sample[self.target_field]] += 1
        
        # Normalize
        total = sum(self.freq.values())
        self.probs = {k: v / total for k, v in self.freq.items()}
        
        # Sort by frequency
        self.ranked = sorted(self.probs.items(), key=lambda x: -x[1])
    
    def predict(self, context: dict = None, k: int = 5) -> List[Tuple[str, float]]:
        """
 Predict top-k
        """
        if not self.ranked:
            return []
        return self.ranked[:k]
    
    def predict_batch(self, contexts: List[dict], k: int = 5) -> List[List[Tuple[str, float]]]:
 """Predict"""
        return [self.predict(ctx, k) for ctx in contexts]
    
    def get_accuracy(self, samples: List[dict]) -> Dict:
        """Compute accuracy"""
        correct_top1 = 0
        correct_top5 = 0
        n = len(samples)
        
        for sample in samples:
            target = sample[self.target_field]
            pred = self.predict(sample, k=5)
            pred_tokens = [t for t, p in pred]
            
            if pred_tokens and pred_tokens[0] == target:
                correct_top1 += 1
            if target in pred_tokens:
                correct_top5 += 1
        
        return {
            'top1_acc': correct_top1 / n if n > 0 else 0,
            'top5_acc': correct_top5 / n if n > 0 else 0,
            'n': n,
            'baseline_type': 'global_frequency'
        }


class ConditionalFrequencyBaseline:
    """
 baseline (metal, cn) Predict
    """
    
    def __init__(self, training_samples: List[dict] = None, 
                 target_field: str = 'target_token',
                 condition_fields: List[str] = None):
        """
        Args:
 training_samples:
 target_field:
 condition_fields: ['metal', 'cn']
        """
        self.target_field = target_field
        self.condition_fields = condition_fields or ['metal', 'cn']
        self.cond_freq = defaultdict(Counter)
        self.global_freq = Counter()
        
        if training_samples:
            self.fit(training_samples)
    
    def _get_condition_key(self, sample: dict) -> tuple:
 """ key"""
        meta = sample.get('meta', sample)
        values = []
        for field in self.condition_fields:
            val = meta.get(field, 'unknown')
            values.append(val)
        return tuple(values)
    
    def fit(self, samples: List[dict]):
 """"""
        self.cond_freq = defaultdict(Counter)
        self.global_freq = Counter()
        
        for sample in samples:
            if self.target_field not in sample:
                continue
            
            target = sample[self.target_field]
            cond_key = self._get_condition_key(sample)
            
            self.cond_freq[cond_key][target] += 1
            self.global_freq[target] += 1
        
        # Precompute sorted order for each condition
        self.cond_ranked = {}
        for cond_key, freq in self.cond_freq.items():
            total = sum(freq.values())
            probs = [(k, v / total) for k, v in freq.items()]
            self.cond_ranked[cond_key] = sorted(probs, key=lambda x: -x[1])
        
        # Global fallback
        total = sum(self.global_freq.values())
        self.global_ranked = sorted(
            [(k, v / total) for k, v in self.global_freq.items()],
            key=lambda x: -x[1]
        )
    
    def predict(self, context: dict, k: int = 5) -> List[Tuple[str, float]]:
        """
 Predict top-k
        """
        cond_key = self._get_condition_key(context)
        
        # Try conditional frequency
        if cond_key in self.cond_ranked:
            return self.cond_ranked[cond_key][:k]
        
        # Fall back to global frequency
        return self.global_ranked[:k]
    
    def predict_batch(self, contexts: List[dict], k: int = 5) -> List[List[Tuple[str, float]]]:
 """Predict"""
        return [self.predict(ctx, k) for ctx in contexts]
    
    def get_accuracy(self, samples: List[dict]) -> Dict:
        """Compute accuracy"""
        correct_top1 = 0
        correct_top5 = 0
        n = len(samples)
        
        for sample in samples:
            target = sample[self.target_field]
            pred = self.predict(sample, k=5)
            pred_tokens = [t for t, p in pred]
            
            if pred_tokens and pred_tokens[0] == target:
                correct_top1 += 1
            if target in pred_tokens:
                correct_top5 += 1
        
        return {
            'top1_acc': correct_top1 / n if n > 0 else 0,
            'top5_acc': correct_top5 / n if n > 0 else 0,
            'n': n,
            'baseline_type': 'conditional_frequency'
        }


class KNNBaseline:
    """
 kNN baseline
 Vote
    """
    
    def __init__(self, training_samples: List[dict] = None,
                 target_field: str = 'target_token',
                 k_neighbors: int = 10):
        """
        Args:
 training_samples:
 target_field:
 k_neighbors:
        """
        self.target_field = target_field
        self.k_neighbors = k_neighbors
        self.samples = []
        self.features = []
        
        if training_samples:
            self.fit(training_samples)
    
    def _extract_features(self, sample: dict) -> np.ndarray:
        """Extract simple features"""
        meta = sample.get('meta', sample)
        
        cn = meta.get('cn', 0)
        metal = meta.get('metal', '')
        geom = meta.get('geom', '')
        
        # Simple hash
        metal_hash = hash(metal) % 100
        geom_hash = hash(geom) % 100
        
        return np.array([cn, metal_hash, geom_hash], dtype=float)
    
    def fit(self, samples: List[dict]):
        """Build index"""
        self.samples = []
        self.features = []
        
        for sample in samples:
            if self.target_field in sample:
                self.samples.append(sample)
                self.features.append(self._extract_features(sample))
        
        if self.features:
            self.features = np.array(self.features)
    
    def predict(self, context: dict, k: int = 5) -> List[Tuple[str, float]]:
        """kNN VotePredict"""
        if len(self.features) == 0:
            return []
        
        query = self._extract_features(context)
        
        # Compute distances
        distances = np.linalg.norm(self.features - query, axis=1)
        
        # Find nearest neighbors
        nn_indices = np.argsort(distances)[:self.k_neighbors]
        
        # Vote
        votes = Counter()
        for idx in nn_indices:
            target = self.samples[idx][self.target_field]
            votes[target] += 1
        
        # Convert to probabilities
        total = sum(votes.values())
        ranked = sorted([(t, c / total) for t, c in votes.items()], key=lambda x: -x[1])
        
        return ranked[:k]
    
    def get_accuracy(self, samples: List[dict]) -> Dict:
        """Compute accuracy"""
        correct_top1 = 0
        correct_top5 = 0
        n = len(samples)
        
        for sample in samples:
            target = sample[self.target_field]
            pred = self.predict(sample, k=5)
            pred_tokens = [t for t, p in pred]
            
            if pred_tokens and pred_tokens[0] == target:
                correct_top1 += 1
            if target in pred_tokens:
                correct_top5 += 1
        
        return {
            'top1_acc': correct_top1 / n if n > 0 else 0,
            'top5_acc': correct_top5 / n if n > 0 else 0,
            'n': n,
            'baseline_type': 'knn'
        }
