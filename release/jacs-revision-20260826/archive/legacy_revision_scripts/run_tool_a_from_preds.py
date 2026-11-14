#!/usr/bin/env python3
"""
Tool A - Predict Fig5

 fig5_preds.jsonl donor Predict
"""

import sys
import json
import pandas as pd
from pathlib import Path
from collections import defaultdict, Counter

sys.path.insert(0, str(Path(__file__).parent.parent))

from coordrep_tools.baselines import DONOR_ATOMS


def main():
    preds_path = '/data/CoordRep/CoordSMILES/libcoordrep/figures/fig5/data/fig5_preds.jsonl'
    output_dir = Path('/data/CoordRep/CoordSMILES/libcoordrep/outputs/fig5')
    output_dir.mkdir(parents=True, exist_ok=True)

    print("Loading predictions...")
    with open(preds_path) as f:
        all_preds = [json.loads(line) for line in f]

    donor_preds = [p for p in all_preds if p['task'] == 'mask_donor']
    print(f"Found {len(donor_preds)} donor predictions")

    # Statistics
    total = len(donor_preds)
    correct_top1 = sum(1 for p in donor_preds if p['correct_top1'])
    correct_top5 = sum(1 for p in donor_preds if p['correct_top5'])

    print(f"\n=== Model Results ===")
    print(f"Top-1 Accuracy: {correct_top1/total*100:.1f}%")
    print(f"Top-5 Accuracy: {correct_top5/total*100:.1f}%")
    print(f"Total: {total}")

    # Stratify by donor
    by_donor = defaultdict(lambda: {'total': 0, 'correct_top1': 0, 'correct_top5': 0})
    for p in donor_preds:
        donor = p['y_true'].get('donor_atom')
        if donor:
            by_donor[donor]['total'] += 1
            if p['correct_top1']:
                by_donor[donor]['correct_top1'] += 1
            if p['correct_top5']:
                by_donor[donor]['correct_top5'] += 1

    # Stratify by CN
    by_cn = defaultdict(lambda: {'total': 0, 'correct_top1': 0, 'correct_top5': 0})
    for p in donor_preds:
        cn = p['group'].get('cn')
        if cn is not None:
            by_cn[cn]['total'] += 1
            if p['correct_top1']:
                by_cn[cn]['correct_top1'] += 1
            if p['correct_top5']:
                by_cn[cn]['correct_top5'] += 1

    # Random baseline
    random_top1 = 1.0 / len(DONOR_ATOMS)
    random_top5 = min(5, len(DONOR_ATOMS)) / len(DONOR_ATOMS)

    # Global frequency baseline
    donor_counts = Counter(p['y_true'].get('donor_atom') for p in donor_preds)
    total_donor = sum(donor_counts.values())
    freq_ranked = [d for d, c in donor_counts.most_common()]

    freq_correct_top1 = sum(1 for p in donor_preds if p['y_true'].get('donor_atom') == freq_ranked[0])
    freq_correct_top5 = sum(1 for p in donor_preds if p['y_true'].get('donor_atom') in freq_ranked[:5])

    metal_donor_freq = defaultdict(Counter)
    for p in donor_preds:
        metal = p['group'].get('metal_element')
        donor = p['y_true'].get('donor_atom')
        if metal and donor:
            metal_donor_freq[metal][donor] += 1

    cond_correct_top1 = 0
    cond_correct_top5 = 0
    for p in donor_preds:
        metal = p['group'].get('metal_element')
        donor = p['y_true'].get('donor_atom')
        if metal in metal_donor_freq:
            ranked = [d for d, c in metal_donor_freq[metal].most_common()]
            if ranked and ranked[0] == donor:
                cond_correct_top1 += 1
            if donor in ranked[:5]:
                cond_correct_top5 += 1

    print(f"\n=== Baselines ===")
    print(f"Random: Top-1 {random_top1*100:.1f}%, Top-5 {random_top5*100:.1f}%")
    print(f"Global Freq: Top-1 {freq_correct_top1/total*100:.1f}%, Top-5 {freq_correct_top5/total*100:.1f}%")
    print(f"Cond Freq: Top-1 {cond_correct_top1/total*100:.1f}%, Top-5 {cond_correct_top5/total*100:.1f}%")

    leaderboard_data = [
        {'method': 'model', 'top1_acc': round(correct_top1/total, 4), 'top5_acc': round(correct_top5/total, 4), 'n': total},
        {'method': 'random', 'top1_acc': round(random_top1, 4), 'top5_acc': round(random_top5, 4), 'n': total},
        {'method': 'global_frequency', 'top1_acc': round(freq_correct_top1/total, 4), 'top5_acc': round(freq_correct_top5/total, 4), 'n': total},
        {'method': 'conditional_frequency', 'top1_acc': round(cond_correct_top1/total, 4), 'top5_acc': round(cond_correct_top5/total, 4), 'n': total},
    ]
    pd.DataFrame(leaderboard_data).to_csv(output_dir / 'fig5b_leaderboard.csv', index=False)
    print(f"\nSaved: {output_dir / 'fig5b_leaderboard.csv'}")

    casecards = []

    successes = [p for p in donor_preds if p['correct_top1']]
    successes_sorted = sorted(successes, key=lambda x: x['pred']['top1_prob'], reverse=True)

    for p in successes_sorted[:4]:
        casecards.append({
            'id': p['id'],
            'metal': p['group'].get('metal_element', 'N/A'),
            'cn': p['group'].get('cn', 'N/A'),
            'target_token': p['y_true'].get('donor_atom'),
            'predictions': [{'token': t, 'prob': round(pr, 4)} for t, pr in p['pred']['topk'][:5]],
            'target_rank': 1,
            'status': 'correct_top1',
        })

    partial = [p for p in donor_preds if p['correct_top5'] and not p['correct_top1']]
    for p in partial[:3]:
        target = p['y_true'].get('donor_atom')
        rank = -1
        for i, (t, pr) in enumerate(p['pred']['topk'][:5]):
            if t == target:
                rank = i + 1
                break
        casecards.append({
            'id': p['id'],
            'metal': p['group'].get('metal_element', 'N/A'),
            'cn': p['group'].get('cn', 'N/A'),
            'target_token': target,
            'predictions': [{'token': t, 'prob': round(pr, 4)} for t, pr in p['pred']['topk'][:5]],
            'target_rank': rank,
            'status': 'correct_top5',
        })

    failures = [p for p in donor_preds if not p['correct_top5']]
    for p in failures[:3]:
        casecards.append({
            'id': p['id'],
            'metal': p['group'].get('metal_element', 'N/A'),
            'cn': p['group'].get('cn', 'N/A'),
            'target_token': p['y_true'].get('donor_atom'),
            'predictions': [{'token': t, 'prob': round(pr, 4)} for t, pr in p['pred']['topk'][:5]],
            'target_rank': -1,
            'status': 'failed',
        })

    with open(output_dir / 'fig5b_casecards.json', 'w') as f:
        json.dump(casecards, f, indent=2)
    print(f"Saved: {output_dir / 'fig5b_casecards.json'}")

    gain_data = []
    for donor, stats in sorted(by_donor.items(), key=lambda x: -x[1]['total']):
        if stats['total'] > 0:
            top1_acc = stats['correct_top1'] / stats['total']
            gain_data.append({
                'donor_atom': donor,
                'top1_acc': round(top1_acc, 4),
                'top5_acc': round(stats['correct_top5'] / stats['total'], 4),
                'n_samples': stats['total'],
                'random_baseline': round(random_top1, 4),
                'x_random': round(top1_acc / random_top1, 2) if random_top1 > 0 else 0,
            })

    pd.DataFrame(gain_data).to_csv(output_dir / 'fig5c_gain_vs_random.csv', index=False)
    print(f"Saved: {output_dir / 'fig5c_gain_vs_random.csv'}")

    cn_data = []
    for cn, stats in sorted(by_cn.items()):
        if stats['total'] > 0:
            cn_data.append({
                'cn': cn,
                'n_samples': stats['total'],
                'top1_acc': round(stats['correct_top1'] / stats['total'], 4),
                'top5_acc': round(stats['correct_top5'] / stats['total'], 4),
            })

    pd.DataFrame(cn_data).to_csv(output_dir / 'fig5_donor_by_cn.csv', index=False)
    print(f"Saved: {output_dir / 'fig5_donor_by_cn.csv'}")

    print("\n=== Tool A Evaluation Complete ===")


if __name__ == '__main__':
    main()
