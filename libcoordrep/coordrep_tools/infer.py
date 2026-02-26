#!/usr/bin/env python3
"""
MLM

 mlm_topk Tool A/B baselines
"""

import sys
import json
import torch
import numpy as np
from pathlib import Path
from typing import List, Tuple, Optional, Union

sys.path.insert(0, str(Path(__file__).parent.parent))

from brain.model import CoordRepForMLM, CoordRepModelConfig
from brain.tokenizer import CoordRepTokenizer, TokenizerConfig


# Global cache (avoid redundant loading)
_MODEL_CACHE = {}
_TOKENIZER_CACHE = {}


def load_model_and_tokenizer(
    checkpoint_path: str,
    tokenizer_path: str = None,
    device: str = "cuda"
) -> Tuple:
    """
    Load model and tokenizer (with caching)
    
    Args:
 checkpoint_path: checkpoint
 tokenizer_path: tokenizer JSON None checkpoint
 device: (cuda/cpu)
    
    Returns:
        (model, tokenizer)
    """
    global _MODEL_CACHE, _TOKENIZER_CACHE
    
    cache_key = f"{checkpoint_path}_{device}"
    
    if cache_key in _MODEL_CACHE:
        return _MODEL_CACHE[cache_key], _TOKENIZER_CACHE[cache_key]
    
    # Infer tokenizer path
    if tokenizer_path is None:
        ckpt_dir = Path(checkpoint_path).parent
        tokenizer_path = str(ckpt_dir / "tokenizer.json")
    
    # Load tokenizer
    with open(tokenizer_path) as f:
        data = json.load(f)
    
    config_dict = data.get('config', {})
    config = TokenizerConfig(**config_dict)
    tokenizer = CoordRepTokenizer(config)
    tokenizer.token2id = data['token2id']
    tokenizer.id2token = {v: k for k, v in tokenizer.token2id.items()}
    
    # Load model
    model_config = CoordRepModelConfig.small(max_length=768)
 model_config.vocab_size = 10000 # checkpoint
    
    model = CoordRepForMLM(model_config)
    checkpoint = torch.load(checkpoint_path, map_location='cpu', weights_only=False)
    model.load_state_dict(checkpoint['model_state_dict'])
    model.eval()
    model = model.to(device)
    
    # Cache
    _MODEL_CACHE[cache_key] = model
    _TOKENIZER_CACHE[cache_key] = tokenizer
    
    return model, tokenizer


def mlm_topk(
    coordrep_masked: str,
    mask_positions: List[int],
    k: int,
    checkpoint_path: str,
    tokenizer_path: str = None,
    device: str = "cuda",
    return_all_probs: bool = False
) -> List[List[Tuple[str, float]]]:
    """
 masked MLM mask position top-k Predict
    
    Args:
 coordrep_masked: [MASK] CoordRep
 mask_positions: mask token (0-indexed, tokenize )
 k: top-k
 checkpoint_path: checkpoint
 tokenizer_path: tokenizer
 device:
 return_all_probs: token
    
    Returns:
 List[List[Tuple[str, float]]]: mask position top-k (token, prob)
    """
    # Set deterministic mode
    torch.manual_seed(42)
    np.random.seed(42)
    
    model, tokenizer = load_model_and_tokenizer(checkpoint_path, tokenizer_path, device)
    
    # Tokenize
    tokens = tokenizer.tokenize(coordrep_masked)
    
    # Convert to token IDs
    token_ids = []
    for tok in tokens:
        if tok == '[MASK]':
            token_ids.append(tokenizer.token2id.get('[MASK]', 3))
        elif tok in tokenizer.token2id:
            token_ids.append(tokenizer.token2id[tok])
        else:
            token_ids.append(tokenizer.token2id.get('[UNK]', 1))
    
    # Limit sequence length
    max_len = 512
    if len(token_ids) > max_len:
        token_ids = token_ids[:max_len]
    
    # Model inference
    results = []
    
    with torch.no_grad():
        input_tensor = torch.tensor([token_ids], device=device)
        attention_mask = torch.ones_like(input_tensor)
        
        logits, _ = model(input_tensor, attention_mask=attention_mask)
        
        for pos in mask_positions:
            if pos >= logits.shape[1]:
                results.append([])
                continue
            
            pos_logits = logits[0, pos, :]
            probs = torch.softmax(pos_logits, dim=-1)
            
            # Top-k
            topk_probs, topk_ids = torch.topk(probs, min(k, probs.shape[0]))
            
            topk_results = []
            for prob, tid in zip(topk_probs.cpu().numpy(), topk_ids.cpu().numpy()):
                tok = tokenizer.id2token.get(int(tid), f'[UNK:{tid}]')
                topk_results.append((tok, float(prob)))
            
            results.append(topk_results)
    
    return results


def mlm_topk_batch(
    samples: List[dict],
    k: int,
    checkpoint_path: str,
    tokenizer_path: str = None,
    device: str = "cuda",
    batch_size: int = 32
) -> List[List[List[Tuple[str, float]]]]:
    """
    Batch MLM inference
    
    Args:
 samples: 'coordrep_masked' 'mask_positions'
        k: top-k
 checkpoint_path: checkpoint
 tokenizer_path: tokenizer
 device:
 batch_size:
    
    Returns:
 mlm_topk
    """
    results = []
    for sample in samples:
        result = mlm_topk(
            sample['coordrep_masked'],
            sample['mask_positions'],
            k,
            checkpoint_path,
            tokenizer_path,
            device
        )
        results.append(result)
    return results


def get_token_probs_at_position(
    coordrep: str,
    position: int,
    candidate_tokens: List[str],
    checkpoint_path: str,
    tokenizer_path: str = None,
    device: str = "cuda"
) -> dict:
    """
    Get candidate token probabilities at specified position
    
    Args:
 coordrep: CoordRep
 position: token
 candidate_tokens: token
 checkpoint_path: checkpoint
 tokenizer_path: tokenizer
 device:
    
    Returns:
 {token: prob}
    """
    model, tokenizer = load_model_and_tokenizer(checkpoint_path, tokenizer_path, device)
    
    # Tokenize
    tokens = tokenizer.tokenize(coordrep)
    
    # Create masked version
    masked_tokens = tokens.copy()
    if position < len(masked_tokens):
        masked_tokens[position] = '[MASK]'
    
    token_ids = []
    for tok in masked_tokens:
        if tok == '[MASK]':
            token_ids.append(tokenizer.token2id.get('[MASK]', 3))
        elif tok in tokenizer.token2id:
            token_ids.append(tokenizer.token2id[tok])
        else:
            token_ids.append(tokenizer.token2id.get('[UNK]', 1))
    
    max_len = 512
    if len(token_ids) > max_len:
        token_ids = token_ids[:max_len]
    
    with torch.no_grad():
        input_tensor = torch.tensor([token_ids], device=device)
        attention_mask = torch.ones_like(input_tensor)
        
        logits, _ = model(input_tensor, attention_mask=attention_mask)
        
        if position >= logits.shape[1]:
            return {tok: 0.0 for tok in candidate_tokens}
        
        pos_logits = logits[0, position, :]
        probs = torch.softmax(pos_logits, dim=-1)
        
        results = {}
        for tok in candidate_tokens:
            if tok in tokenizer.token2id:
                tid = tokenizer.token2id[tok]
                results[tok] = float(probs[tid].cpu())
            else:
                results[tok] = 0.0
        
        return results
