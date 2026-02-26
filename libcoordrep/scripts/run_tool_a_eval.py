#!/usr/bin/env python3
"""
Tool A

 Donor Fig5b
"""

import sys
import json
import argparse
import pandas as pd
from pathlib import Path
from typing import Dict
from collections import defaultdict

sys.path.insert(0, str(Path(__file__).parent.parent))

from coordrep_tools.tool_a_donor import (
    DonorPredictor, 
    evaluate_donor_prediction, 
    compare_with_baselines,
    generate_casecards
)
from coordrep_tools.baselines import DONOR_ATOMS


def create_donor_eval_dataset(
    source_path: str,
    output_path: str,
    n_samples: int = 5000
) -> None:
    """
 donor
    """
    print(f"Loading source data from {source_path}...")
    
    with open(source_path) as f:
        source_samples = [json.loads(line) for line in f]
    
    donor_samples = [s for s in source_samples if s.get('task') == 'mask_donor']
    
    print(f"Found {len(donor_samples)} donor mask samples")
    
    eval_samples = []
    for i, s in enumerate(donor_samples[:n_samples]):
        eval_samples.append({
            'id': s.get('id', f'donor_{i}'),
            'coordrep_masked': ''.join(s.get('masked_tokens', [])),
            'mask_positions': [s.get('mask_pos', 0)],
            'target_token': s.get('y_true', {}).get('donor_atom'),
            'meta': s.get('group', {})
        })
    
    # Save
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_path, 'w') as f:
        for sample in eval_samples:
            f.write(json.dumps(sample) + '\n')
    
    print(f"Saved {len(eval_samples)} samples to {output_path}")


def run_evaluation(
    eval_data_path: str,
    checkpoint_path: str,
    tokenizer_path: str,
    output_dir: str,
    device: str = "cuda"
) -> Dict:
 """"""
    
    print("=" * 50)
    print("Tool A: Donor Prediction Evaluation")
    print("=" * 50)
    
    # Load data
    print(f"\nLoading evaluation data from {eval_data_path}...")
    with open(eval_data_path) as f:
        samples = [json.loads(line) for line in f]
    
    print(f"Loaded {len(samples)} samples")
    
    print("\nRunning model evaluation...")
    model_results = evaluate_donor_prediction(
        samples, checkpoint_path, tokenizer_path, device
    )
    
    print(f"\n=== Model Results ===")
    print(f"Top-1 Accuracy: {model_results['top1_acc']*100:.1f}%")
    print(f"Top-5 Accuracy: {model_results['top5_acc']*100:.1f}%")
    print(f"Total samples: {model_results['total']}")
    
    print("\nComparing with baselines...")
    comparison = compare_with_baselines(
        samples, checkpoint_path, tokenizer_path, device,
 training_samples=samples # baseline
    )
    
    print(f"\n=== Comparison ===")
    for method, res in comparison.items():
        print(f"  {method}: Top-1 {res['top1_acc']*100:.1f}%, Top-5 {res['top5_acc']*100:.1f}%")
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 1. fig5b_leaderboard.csv
    leaderboard_data = []
    for method, res in comparison.items():
        leaderboard_data.append({
            'method': method,
            'top1_acc': round(res['top1_acc'], 4),
            'top5_acc': round(res['top5_acc'], 4),
            'n': res['n'],
        })
    
    leaderboard_df = pd.DataFrame(leaderboard_data)
    leaderboard_df.to_csv(output_dir / 'fig5b_leaderboard.csv', index=False)
    print(f"\nSaved: {output_dir / 'fig5b_leaderboard.csv'}")
    
    # 2. fig5b_casecards.json
    casecards = generate_casecards(model_results['predictions'], n_cards=10)
    with open(output_dir / 'fig5b_casecards.json', 'w') as f:
        json.dump(casecards, f, indent=2)
    print(f"Saved: {output_dir / 'fig5b_casecards.json'}")
    
    # 3. fig5c_gain_vs_random.csv
    random_baseline = 1.0 / len(DONOR_ATOMS)
    gain_data = []
    for donor, stats in model_results['by_donor'].items():
        if stats['total'] > 0:
            top1_acc = stats['correct_top1'] / stats['total']
            gain_data.append({
                'donor_atom': donor,
                'top1_acc': round(top1_acc, 4),
                'n_samples': stats['total'],
                'random_baseline': round(random_baseline, 4),
                'x_random': round(top1_acc / random_baseline, 2) if random_baseline > 0 else 0,
            })
    
    gain_df = pd.DataFrame(gain_data)
    gain_df = gain_df.sort_values('n_samples', ascending=False)
    gain_df.to_csv(output_dir / 'fig5c_gain_vs_random.csv', index=False)
    print(f"Saved: {output_dir / 'fig5c_gain_vs_random.csv'}")
    
    print("\n=== Evaluation Complete ===")
    
    return {
        'model': model_results,
        'comparison': comparison,
    }


def main():
    parser = argparse.ArgumentParser(description='Run Tool A Evaluation')
    parser.add_argument('--eval_data', type=str, 
                        default='/data/CoordRep/CoordSMILES/libcoordrep/figures/fig5/data/fig5_tasks.jsonl',
                        help='Evaluation data path (or source to create from)')
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
    
    args = parser.parse_args()
    
    if args.create_eval_data:
        donor_eval_path = Path(args.output_dir) / 'donor_eval.jsonl'
        create_donor_eval_dataset(args.eval_data, str(donor_eval_path))
        args.eval_data = str(donor_eval_path)
    
    run_evaluation(
        args.eval_data,
        args.checkpoint,
        args.tokenizer,
        args.output_dir,
        args.device
    )


if __name__ == '__main__':
    main()
