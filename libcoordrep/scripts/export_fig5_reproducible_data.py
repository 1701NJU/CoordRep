#!/usr/bin/env python3
"""
Export complete reproducible data for all Fig5 B/C/D/E panels

Each panel's data is exported to separate CSV files with complete column names and descriptions
"""

import sys
import json
import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict, Counter

sys.path.insert(0, str(Path(__file__).parent.parent))


def load_predictions(preds_path: str):
    """Load prediction data"""
    with open(preds_path) as f:
        all_preds = [json.loads(line) for line in f]
    return all_preds


def export_panel_b_data(all_preds: list, output_dir: Path):
    """
 Panel B: ΔTop1
    
    - panel_b_delta_histogram.csv: Histogram bin data
 - panel_b_delta_by_metal_cn.csv: (metal, CN)
 - panel_b_top_gains.csv: Top 20
 - panel_b_casecards.csv: CondFreq
    """
    print("\n=== Exporting Panel B Data ===")
    
    donor_preds = [p for p in all_preds if p['task'] == 'mask_donor']
    
    # Compute conditional frequency baseline
    cond_freq = defaultdict(Counter)
    for p in donor_preds:
        metal = p['group'].get('metal_element')
        cn = p['group'].get('cn')
        donor = p['y_true'].get('donor_atom')
        if metal and cn is not None and donor:
            cond_freq[(metal, cn)][donor] += 1
    
    cond_ranked = {}
    for key, counter in cond_freq.items():
        cond_ranked[key] = [d for d, c in counter.most_common()]
    
    # Group statistics by (metal, CN)
    groups = defaultdict(lambda: {
        'total': 0, 
        'model_top1': 0, 
        'condfreq_top1': 0,
        'model_top5': 0,
        'condfreq_top5': 0,
    })
    
 corrections = [] # CondFreq
    
    for p in donor_preds:
        metal = p['group'].get('metal_element')
        cn = p['group'].get('cn')
        donor = p['y_true'].get('donor_atom')
        
        if not metal or cn is None:
            continue
        
        key = (metal, cn)
        groups[key]['total'] += 1
        
        # Model predictions
        model_top1_correct = p['correct_top1']
        model_top5_correct = p['correct_top5']
        
        if model_top1_correct:
            groups[key]['model_top1'] += 1
        if model_top5_correct:
            groups[key]['model_top5'] += 1
        
        # CondFreq predictions
        condfreq_pred = cond_ranked.get(key, [''])[0]
        condfreq_top5 = cond_ranked.get(key, [])[:5]
        
        condfreq_top1_correct = condfreq_pred == donor
        condfreq_top5_correct = donor in condfreq_top5
        
        if condfreq_top1_correct:
            groups[key]['condfreq_top1'] += 1
        if condfreq_top5_correct:
            groups[key]['condfreq_top5'] += 1
        
        if model_top1_correct and not condfreq_top1_correct:
            model_topk = p['pred']['topk'][:5]
            corrections.append({
                'sample_id': p['id'],
                'metal': metal,
                'cn': cn,
                'target_donor': donor,
                'condfreq_prediction': condfreq_pred,
                'model_prediction': model_topk[0][0] if model_topk else '',
                'model_confidence': model_topk[0][1] if model_topk else 0,
                'model_top2': model_topk[1][0] if len(model_topk) > 1 else '',
                'model_top2_prob': model_topk[1][1] if len(model_topk) > 1 else 0,
                'model_top3': model_topk[2][0] if len(model_topk) > 2 else '',
                'model_top3_prob': model_topk[2][1] if len(model_topk) > 2 else 0,
            })
    
    # 1. Export detailed data for each (metal, CN) group
    metal_cn_data = []
    for (metal, cn), stats in groups.items():
 if stats['total'] >= 10: # 10
            model_top1_acc = stats['model_top1'] / stats['total']
            condfreq_top1_acc = stats['condfreq_top1'] / stats['total']
            model_top5_acc = stats['model_top5'] / stats['total']
            condfreq_top5_acc = stats['condfreq_top5'] / stats['total']
            
            metal_cn_data.append({
                'metal': metal,
                'cn': cn,
                'n_samples': stats['total'],
                'model_top1_correct': stats['model_top1'],
                'model_top1_acc': round(model_top1_acc, 6),
                'condfreq_top1_correct': stats['condfreq_top1'],
                'condfreq_top1_acc': round(condfreq_top1_acc, 6),
                'delta_top1': round(model_top1_acc - condfreq_top1_acc, 6),
                'delta_top1_pp': round((model_top1_acc - condfreq_top1_acc) * 100, 2),
                'model_top5_acc': round(model_top5_acc, 6),
                'condfreq_top5_acc': round(condfreq_top5_acc, 6),
                'delta_top5': round(model_top5_acc - condfreq_top5_acc, 6),
            })
    
    df_metal_cn = pd.DataFrame(metal_cn_data)
    df_metal_cn = df_metal_cn.sort_values('delta_top1', ascending=False)
    df_metal_cn.to_csv(output_dir / 'panel_b_delta_by_metal_cn.csv', index=False)
    print(f"  Saved: panel_b_delta_by_metal_cn.csv ({len(df_metal_cn)} rows)")
    
    # 2. Export histogram bin data
    deltas = df_metal_cn['delta_top1_pp'].values
    bins = np.linspace(-40, 45, 18)
    hist_counts, hist_edges = np.histogram(deltas, bins=bins)
    
    hist_data = []
    for i in range(len(hist_counts)):
        bin_center = (hist_edges[i] + hist_edges[i+1]) / 2
        hist_data.append({
            'bin_left': round(hist_edges[i], 2),
            'bin_right': round(hist_edges[i+1], 2),
            'bin_center': round(bin_center, 2),
            'count': int(hist_counts[i]),
            'is_positive': bin_center >= 0,
        })
    
    df_hist = pd.DataFrame(hist_data)
    df_hist.to_csv(output_dir / 'panel_b_delta_histogram.csv', index=False)
    print(f"  Saved: panel_b_delta_histogram.csv ({len(df_hist)} bins)")
    
    # 3. Export top gains
    df_top_gains = df_metal_cn.head(20).copy()
    df_top_gains.to_csv(output_dir / 'panel_b_top_gains.csv', index=False)
    print(f"  Saved: panel_b_top_gains.csv (20 rows)")
    
    # 4. Export correction cases
    corrections_sorted = sorted(corrections, key=lambda x: x['model_confidence'], reverse=True)
    df_corrections = pd.DataFrame(corrections_sorted[:30])  # Top 30
    df_corrections.to_csv(output_dir / 'panel_b_casecards.csv', index=False)
    print(f"  Saved: panel_b_casecards.csv ({len(df_corrections)} cases)")
    
    # 5. Export statistical summary
    summary = {
        'total_groups': len(df_metal_cn),
        'model_wins': int((df_metal_cn['delta_top1'] > 0).sum()),
        'condfreq_wins': int((df_metal_cn['delta_top1'] < 0).sum()),
        'ties': int((df_metal_cn['delta_top1'] == 0).sum()),
        'mean_delta_top1_pp': round(df_metal_cn['delta_top1_pp'].mean(), 2),
        'median_delta_top1_pp': round(df_metal_cn['delta_top1_pp'].median(), 2),
        'std_delta_top1_pp': round(df_metal_cn['delta_top1_pp'].std(), 2),
        'total_corrections': len(corrections),
    }
    df_summary = pd.DataFrame([summary])
    df_summary.to_csv(output_dir / 'panel_b_summary.csv', index=False)
    print(f"  Saved: panel_b_summary.csv")
    
    return df_metal_cn


def export_panel_c_data(all_preds: list, output_dir: Path):
    """
 Panel C: Confidence-Coverage + Donor-type Δ
    
    - panel_c_confidence_coverage.csv: Confidence-coverage curve data
    - panel_c_donor_type_delta.csv: Delta data by donor type
    """
    print("\n=== Exporting Panel C Data ===")
    
    donor_preds = [p for p in all_preds if p['task'] == 'mask_donor']
    
    # Compute conditional frequency baseline
    cond_freq = defaultdict(Counter)
    for p in donor_preds:
        metal = p['group'].get('metal_element')
        cn = p['group'].get('cn')
        donor = p['y_true'].get('donor_atom')
        if metal and cn is not None and donor:
            cond_freq[(metal, cn)][donor] += 1
    
    cond_ranked = {}
    for key, counter in cond_freq.items():
        cond_ranked[key] = [d for d, c in counter.most_common()]
    
    data = []
    for p in donor_preds:
        if not p['pred']['topk']:
            continue
        
        confidence = p['pred']['topk'][0][1]
        correct = p['correct_top1']
        data.append((confidence, correct))
    
    # Sort by confidence
    data.sort(key=lambda x: x[0], reverse=True)
    
    # Compute at different thresholds
    thresholds = list(np.arange(0.05, 1.0, 0.05))
    thresholds.extend([0.85, 0.90, 0.95, 0.99])
    thresholds = sorted(set(thresholds))
    
    cc_data = []
    total = len(data)
    
    for threshold in thresholds:
        filtered = [(conf, corr) for conf, corr in data if conf >= threshold]
        
        if not filtered:
            continue
        
        coverage = len(filtered) / total
        accuracy = sum(corr for _, corr in filtered) / len(filtered)
        n_correct = sum(corr for _, corr in filtered)
        
        cc_data.append({
            'threshold': round(threshold, 3),
            'n_samples': len(filtered),
            'n_correct': n_correct,
            'coverage': round(coverage, 6),
            'coverage_pct': round(coverage * 100, 2),
            'accuracy': round(accuracy, 6),
            'accuracy_pct': round(accuracy * 100, 2),
        })
    
    df_cc = pd.DataFrame(cc_data)
    df_cc.to_csv(output_dir / 'panel_c_confidence_coverage.csv', index=False)
    print(f"  Saved: panel_c_confidence_coverage.csv ({len(df_cc)} thresholds)")
    
    by_donor = defaultdict(lambda: {
        'total': 0, 
        'model_top1': 0, 
        'condfreq_top1': 0,
        'model_top5': 0,
    })
    
    for p in donor_preds:
        metal = p['group'].get('metal_element')
        cn = p['group'].get('cn')
        donor = p['y_true'].get('donor_atom')
        
        if not donor:
            continue
        
        by_donor[donor]['total'] += 1
        
        if p['correct_top1']:
            by_donor[donor]['model_top1'] += 1
        if p['correct_top5']:
            by_donor[donor]['model_top5'] += 1
        
        key = (metal, cn)
        if key in cond_ranked:
            if cond_ranked[key][0] == donor:
                by_donor[donor]['condfreq_top1'] += 1
    
 random_baseline = 1.0 / 15 # 15 donor atoms
    
    donor_data = []
    for donor, stats in by_donor.items():
        if stats['total'] < 10:
            continue
        
        model_acc = stats['model_top1'] / stats['total']
        condfreq_acc = stats['condfreq_top1'] / stats['total']
        model_top5_acc = stats['model_top5'] / stats['total']
        
        donor_data.append({
            'donor_atom': donor,
            'n_samples': stats['total'],
            'model_top1_correct': stats['model_top1'],
            'model_top1_acc': round(model_acc, 6),
            'model_top1_acc_pct': round(model_acc * 100, 2),
            'condfreq_top1_correct': stats['condfreq_top1'],
            'condfreq_top1_acc': round(condfreq_acc, 6),
            'condfreq_top1_acc_pct': round(condfreq_acc * 100, 2),
            'delta_top1': round(model_acc - condfreq_acc, 6),
            'delta_top1_pp': round((model_acc - condfreq_acc) * 100, 2),
            'model_top5_acc': round(model_top5_acc, 6),
            'random_baseline': round(random_baseline, 6),
            'x_random': round(model_acc / random_baseline, 2) if random_baseline > 0 else 0,
        })
    
    df_donor = pd.DataFrame(donor_data)
    df_donor = df_donor.sort_values('n_samples', ascending=False)
    df_donor.to_csv(output_dir / 'panel_c_donor_type_delta.csv', index=False)
    print(f"  Saved: panel_c_donor_type_delta.csv ({len(df_donor)} donor types)")
    
    return df_cc, df_donor


def export_panel_d_data(output_dir: Path):
    """
 Panel D: Validity Uplift
    
 - panel_d_validity_uplift.csv:
    """
    print("\n=== Exporting Panel D Data ===")
    
    # Read existing uplift data
    uplift_path = output_dir / 'fig5d_validity_uplift.csv'
    repair_results_path = output_dir / 'fig5d_repair_results.json'
    
    if uplift_path.exists():
        df_uplift = pd.read_csv(uplift_path)
        
        # Reformat to ensure clear column names
        df_formatted = df_uplift.copy()
        df_formatted = df_formatted.rename(columns={
            'corruption_type': 'corruption_type',
            'n_samples': 'n_samples',
            'valid_before_corruption': 'valid_before_corruption_rate',
            'valid_after_corruption': 'valid_after_corruption_rate',
            'valid_after_repair': 'valid_after_repair_rate',
            'uplift': 'repair_uplift',
        })
        
        # Add percentage columns
        if 'valid_after_corruption_rate' in df_formatted.columns:
            df_formatted['valid_after_corruption_pct'] = (df_formatted['valid_after_corruption_rate'] * 100).round(2)
            df_formatted['valid_after_repair_pct'] = (df_formatted['valid_after_repair_rate'] * 100).round(2)
            df_formatted['repair_uplift_pp'] = (df_formatted['repair_uplift'] * 100).round(2)
        
        df_formatted.to_csv(output_dir / 'panel_d_validity_uplift.csv', index=False)
        print(f"  Saved: panel_d_validity_uplift.csv ({len(df_formatted)} rows)")
        
        return df_formatted
    else:
        print("  Warning: fig5d_validity_uplift.csv not found")
        return None


def export_panel_e_data(all_preds: list, output_dir: Path):
    """
 Panel E: CN
    
 - panel_e_cn_stratified.csv: Stratify by CN
    """
    print("\n=== Exporting Panel E Data ===")
    
    struct_preds = [p for p in all_preds if p['task'] == 'mask_structure']
    
    # Group statistics by CN
    by_cn = defaultdict(lambda: {
        'total': 0, 
        'correct_top1': 0, 
        'correct_top5': 0,
    })
    
    for p in struct_preds:
        cn = p['group'].get('cn')
        if cn is None:
            continue
        
        by_cn[cn]['total'] += 1
        if p['correct_top1']:
            by_cn[cn]['correct_top1'] += 1
        if p['correct_top5']:
            by_cn[cn]['correct_top5'] += 1
    
    cn_data = []
    for cn in sorted(by_cn.keys()):
        stats = by_cn[cn]
 if stats['total'] >= 10: # 10
            top1_acc = stats['correct_top1'] / stats['total']
            top5_acc = stats['correct_top5'] / stats['total']
            
            cn_data.append({
                'cn': cn,
                'n_samples': stats['total'],
                'top1_correct': stats['correct_top1'],
                'top1_acc': round(top1_acc, 6),
                'top1_acc_pct': round(top1_acc * 100, 2),
                'top5_correct': stats['correct_top5'],
                'top5_acc': round(top5_acc, 6),
                'top5_acc_pct': round(top5_acc * 100, 2),
            })
    
    df_cn = pd.DataFrame(cn_data)
    df_cn.to_csv(output_dir / 'panel_e_cn_stratified.csv', index=False)
    print(f"  Saved: panel_e_cn_stratified.csv ({len(df_cn)} CN values)")
    
    # Export overall statistics
    total = sum(s['total'] for s in by_cn.values())
    correct_top1 = sum(s['correct_top1'] for s in by_cn.values())
    correct_top5 = sum(s['correct_top5'] for s in by_cn.values())
    
    summary = {
        'total_samples': total,
        'top1_correct': correct_top1,
        'top1_acc': round(correct_top1 / total, 6) if total > 0 else 0,
        'top1_acc_pct': round(correct_top1 / total * 100, 2) if total > 0 else 0,
        'top5_correct': correct_top5,
        'top5_acc': round(correct_top5 / total, 6) if total > 0 else 0,
        'top5_acc_pct': round(correct_top5 / total * 100, 2) if total > 0 else 0,
    }
    df_summary = pd.DataFrame([summary])
    df_summary.to_csv(output_dir / 'panel_e_summary.csv', index=False)
    print(f"  Saved: panel_e_summary.csv")
    
    return df_cn


def main():
    preds_path = '/data/CoordRep/CoordSMILES/libcoordrep/figures/fig5/data/fig5_preds.jsonl'
    output_dir = Path('/data/CoordRep/CoordSMILES/libcoordrep/outputs/fig5/reproducible_data')
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 60)
    print("Exporting Fig5 Reproducible Data")
    print("=" * 60)
    
    print(f"\nLoading predictions from {preds_path}...")
    all_preds = load_predictions(preds_path)
    print(f"Loaded {len(all_preds)} predictions")
    
    # Export data for each panel
    export_panel_b_data(all_preds, output_dir)
    export_panel_c_data(all_preds, output_dir)
 export_panel_d_data(output_dir.parent) #
    export_panel_e_data(all_preds, output_dir)
    
    import shutil
    src_uplift = output_dir.parent / 'fig5d_validity_uplift.csv'
    if src_uplift.exists():
        shutil.copy(src_uplift, output_dir / 'panel_d_validity_uplift.csv')
        print(f"\n  Copied: panel_d_validity_uplift.csv")
    
    print("\n" + "=" * 60)
    print("All Reproducible Data Exported!")
    print("=" * 60)
    
    # Print file list
    print(f"\nOutput directory: {output_dir}")
    print("\nGenerated files:")
    for f in sorted(output_dir.glob('*.csv')):
        size = f.stat().st_size
        print(f"  {f.name} ({size:,} bytes)")


if __name__ == '__main__':
    main()
