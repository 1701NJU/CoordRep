#!/usr/bin/env python3
"""
Tool B

 Fig5d
"""

import sys
import json
import argparse
import pandas as pd
from pathlib import Path
from typing import Dict
from collections import defaultdict

sys.path.insert(0, str(Path(__file__).parent.parent))

from coordrep_tools.tool_b_repair import (
    CoordRepRepairer,
    evaluate_structure_recovery,
    evaluate_repair_on_synthetic
)
from coordrep_tools.validate import is_valid_coordrep, classify_error_type


def create_struct_eval_dataset(
    source_path: str,
    output_path: str,
    n_samples: int = 5000
) -> None:
    """
 structure token
    """
    print(f"Loading source data from {source_path}...")
    
    with open(source_path) as f:
        source_samples = [json.loads(line) for line in f]
    
    struct_samples = [s for s in source_samples if s.get('task') == 'mask_structure']
    
    print(f"Found {len(struct_samples)} structure mask samples")
    
    eval_samples = []
    for i, s in enumerate(struct_samples[:n_samples]):
        eval_samples.append({
            'id': s.get('id', f'struct_{i}'),
            'coordrep_masked': ''.join(s.get('masked_tokens', [])),
            'mask_positions': [s.get('mask_pos', 0)],
            'target_tokens': [s.get('y_true', {}).get('structure_token')],
            'meta': s.get('group', {})
        })
    
    # Save
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_path, 'w') as f:
        for sample in eval_samples:
            f.write(json.dumps(sample) + '\n')
    
    print(f"Saved {len(eval_samples)} samples to {output_path}")


def run_structure_recovery_eval(
    eval_data_path: str,
    checkpoint_path: str,
    tokenizer_path: str,
    output_dir: str,
    device: str = "cuda"
) -> Dict:
 """ token """
    
    print("=" * 50)
    print("Tool B: Structure Token Recovery Evaluation")
    print("=" * 50)
    
    # Load data
    print(f"\nLoading evaluation data from {eval_data_path}...")
    with open(eval_data_path) as f:
        samples = [json.loads(line) for line in f]
    
    print(f"Loaded {len(samples)} samples")
    
    print("\nRunning evaluation...")
    results = evaluate_structure_recovery(
        samples, checkpoint_path, tokenizer_path, device
    )
    
    print(f"\n=== Results ===")
    print(f"Top-1 Accuracy: {results['top1_acc']*100:.1f}%")
    print(f"Top-5 Accuracy: {results['top5_acc']*100:.1f}%")
    print(f"Total positions: {results['total_positions']}")
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    cn_data = []
    for cn, stats in sorted(results['by_cn'].items()):
        if stats['total'] > 0:
            cn_data.append({
                'cn': cn,
                'n_samples': stats['total'],
                'top1_acc': round(stats['correct_top1'] / stats['total'], 4),
                'top5_acc': round(stats['correct_top5'] / stats['total'], 4),
            })
    
    cn_df = pd.DataFrame(cn_data)
    cn_df.to_csv(output_dir / 'fig5d_cn_stratified.csv', index=False)
    print(f"\nSaved: {output_dir / 'fig5d_cn_stratified.csv'}")
    
    return results


def run_repair_eval(
    corrupted_data_path: str,
    checkpoint_path: str,
    tokenizer_path: str,
    output_dir: str,
    device: str = "cuda"
) -> Dict:
 """"""
    
    print("=" * 50)
    print("Tool B: Repair Evaluation on Synthetic Corruption")
    print("=" * 50)
    
    # Load data
    print(f"\nLoading corrupted data from {corrupted_data_path}...")
    with open(corrupted_data_path) as f:
        samples = [json.loads(line) for line in f]
    
    print(f"Loaded {len(samples)} samples")
    
    print("\nRunning repair evaluation...")
    results = evaluate_repair_on_synthetic(
        samples, checkpoint_path, tokenizer_path
    )
    
    print(f"\n=== Results ===")
    print(f"Valid rate before: {results['valid_rate_before']*100:.1f}%")
    print(f"Valid rate after:  {results['valid_rate_after']*100:.1f}%")
    print(f"Uplift: +{results['uplift']*100:.1f}%")
    print(f"Exact match rate: {results['exact_match_rate']*100:.1f}%")
    
    output_dir = Path(output_dir)
    
    # 1. fig5d_validity_uplift.csv
    uplift_data = [{
        'metric': 'overall',
        'valid_rate_before': round(results['valid_rate_before'], 4),
        'valid_rate_after': round(results['valid_rate_after'], 4),
        'uplift': round(results['uplift'], 4),
        'n_samples': results['total'],
    }]
    
    for ctype, stats in results['by_corruption_type'].items():
        if stats['total'] > 0:
            uplift_data.append({
                'metric': ctype,
                'valid_rate_before': round(stats['valid_before'] / stats['total'], 4),
                'valid_rate_after': round(stats['valid_after'] / stats['total'], 4),
                'uplift': round((stats['valid_after'] - stats['valid_before']) / stats['total'], 4),
                'n_samples': stats['total'],
            })
    
    uplift_df = pd.DataFrame(uplift_data)
    uplift_df.to_csv(output_dir / 'fig5d_validity_uplift.csv', index=False)
    print(f"\nSaved: {output_dir / 'fig5d_validity_uplift.csv'}")
    
    # 2. fig5d_error_types.csv
    error_data = []
    for error_type, count in results['error_types'].items():
        error_data.append({
            'error_type': error_type,
            'count': count,
            'rate': round(count / results['total'], 4) if results['total'] > 0 else 0,
        })
    
    error_df = pd.DataFrame(error_data)
    error_df = error_df.sort_values('count', ascending=False)
    error_df.to_csv(output_dir / 'fig5d_error_types.csv', index=False)
    print(f"Saved: {output_dir / 'fig5d_error_types.csv'}")
    
    return results


def main():
    parser = argparse.ArgumentParser(description='Run Tool B Evaluation')
    parser.add_argument('--struct_data', type=str,
                        default='/data/CoordRep/CoordSMILES/libcoordrep/figures/fig5/data/fig5_tasks.jsonl',
                        help='Structure eval data path')
    parser.add_argument('--corrupted_data', type=str,
                        default=None,
                        help='Corrupted data path for repair eval')
    parser.add_argument('--checkpoint', type=str,
                        default='/data/CoordRep/CoordSMILES/libcoordrep/checkpoints/pretrain_v3/best_model.pt',
                        help='Model checkpoint path')
    parser.add_argument('--tokenizer', type=str,
                        default='/data/CoordRep/CoordSMILES/libcoordrep/checkpoints/pretrain_v3/tokenizer.json',
                        help='Tokenizer path')
    parser.add_argument('--output_dir', type=str,
                        default='/data/CoordRep/CoordSMILES/libcoordrep/outputs/fig5',
                        help='Output directory')
    parser.add_argument('--device', type=str, default='cuda',
                        help='Device')
    parser.add_argument('--create_eval_data', action='store_true',
                        help='Create evaluation data from source')
    parser.add_argument('--create_corrupted', action='store_true',
                        help='Create synthetic corrupted data')
    parser.add_argument('--skip_repair', action='store_true',
                        help='Skip repair evaluation')
    
    args = parser.parse_args()
    
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    struct_eval_path = args.struct_data
    if args.create_eval_data:
        struct_eval_path = str(output_dir / 'struct_eval.jsonl')
        create_struct_eval_dataset(args.struct_data, struct_eval_path)
    
    run_structure_recovery_eval(
        struct_eval_path,
        args.checkpoint,
        args.tokenizer,
        args.output_dir,
        args.device
    )
    
    if not args.skip_repair:
        corrupted_path = args.corrupted_data
        
        if args.create_corrupted or corrupted_path is None:
            from build_synth_corruption import generate_corrupted_samples
            
            print("\n" + "=" * 50)
            print("Creating synthetic corrupted data...")
            
            clean_data_path = '/data/CoordRep/CoordSMILES/libcoordrep/figures/fig5/data/geometry_cloze_test.jsonl'
            if Path(clean_data_path).exists():
                with open(clean_data_path) as f:
                    clean_samples = [json.loads(line) for line in f][:5000]
                
                corrupted_samples = generate_corrupted_samples(clean_samples, n_output=2000)
                
                corrupted_path = str(output_dir / 'synth_corrupted.jsonl')
                with open(corrupted_path, 'w') as f:
                    for s in corrupted_samples:
                        f.write(json.dumps(s) + '\n')
                
                print(f"Saved corrupted data to {corrupted_path}")
        
        if corrupted_path and Path(corrupted_path).exists():
            run_repair_eval(
                corrupted_path,
                args.checkpoint,
                args.tokenizer,
                args.output_dir,
                args.device
            )
        else:
            print("\nSkipping repair evaluation: no corrupted data available")
    
    print("\n=== All Evaluations Complete ===")


if __name__ == '__main__':
    main()
