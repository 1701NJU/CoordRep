#!/usr/bin/env python3
"""
Generate Fig5 B/C/D panel plots

Fig5B: ΔTop1 + Top gains + Case Cards
Fig5C: Confidence-Coverage + Donor-type Δ
Fig5D: Valid Parse Uplift + CN
"""

import sys
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 11
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['xtick.labelsize'] = 9
plt.rcParams['ytick.labelsize'] = 9

# Color scheme
COLORS = {
    'model': '#2E86AB',
    'baseline': '#A23B72',
    'positive': '#28A745',
    'negative': '#DC3545',
    'neutral': '#6C757D',
}


def plot_fig5b(output_dir: Path):
    """
 Fig5B: ΔTop1 + Top gains
    """
    print("\n=== Plotting Fig5B ===")
    
    delta_path = output_dir / "fig5_delta_top1_by_metal_cn.csv"
    gains_path = output_dir / "fig5_delta_top1_gains.csv"
    
    df = pd.read_csv(delta_path)
    gains = pd.read_csv(gains_path)
    
    fig = plt.figure(figsize=(7.2, 3.2))
    
    ax = fig.add_axes([0.10, 0.18, 0.52, 0.72])
    
    deltas = df['delta_top1'].values * 100  # Convert to percentage points
    
    # Color by sign: positive=blue, negative=red
    bins = np.linspace(-40, 45, 18)
    n, bins_out, patches = ax.hist(deltas, bins=bins, edgecolor='white', linewidth=0.5)
    
    # Color by value
    for patch, left, right in zip(patches, bins_out[:-1], bins_out[1:]):
        center = (left + right) / 2
        if center >= 0:
            patch.set_facecolor(COLORS['positive'])
            patch.set_alpha(0.7)
        else:
            patch.set_facecolor(COLORS['negative'])
            patch.set_alpha(0.7)
    
    ax.axvline(0, color='black', linewidth=1.5, linestyle='-')
    ax.axvline(4.26, color=COLORS['model'], linewidth=2, linestyle='--', label='Mean +4.26pp')
    
    ax.set_xlabel(r'$\Delta$Top-1 (Model − CondFreq), percentage points', fontsize=11)
    ax.set_ylabel('Count of (metal, CN) groups', fontsize=11)
    ax.set_title('Beyond Conditional Lookup', fontsize=12, fontweight='bold')
    ax.legend(loc='upper left', frameon=False)
    
    # Add statistical annotations
    n_pos = (df['delta_top1'] > 0).sum()
    n_neg = (df['delta_top1'] < 0).sum()
    ax.text(0.97, 0.95, f'Model wins: {n_pos}/179 (67%)\nCondFreq wins: {n_neg}/179 (17%)',
            transform=ax.transAxes, ha='right', va='top', fontsize=9,
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    ax2 = fig.add_axes([0.68, 0.18, 0.30, 0.72])
    ax2.axis('off')
    
    top5 = gains.head(5)
    cell_text = []
    for _, r in top5.iterrows():
        cell_text.append([
            r['metal'], 
            str(int(r['cn'])), 
            f"+{r['delta_top1']*100:.1f}"
        ])
    
    tbl = ax2.table(
        cellText=cell_text,
        colLabels=['Metal', 'CN', r'$\Delta$pp'],
        loc='center',
        cellLoc='center',
        colColours=['#E8E8E8'] * 3,
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(10)
    tbl.scale(1.0, 1.4)
    
    # Set table style
    for (row, col), cell in tbl.get_celld().items():
        if row == 0:
            cell.set_text_props(fontweight='bold')
        cell.set_edgecolor('#CCCCCC')
    
    ax2.set_title('Top Gains', fontsize=11, fontweight='bold', pad=10)
    
    fig.savefig(output_dir / "fig5b_delta_distribution.pdf", bbox_inches='tight', dpi=150)
    fig.savefig(output_dir / "fig5b_delta_distribution.svg", bbox_inches='tight')
    fig.savefig(output_dir / "fig5b_delta_distribution.png", bbox_inches='tight', dpi=150)
    plt.close(fig)
    
    print(f"Saved: {output_dir / 'fig5b_delta_distribution.pdf'}")


def plot_fig5b_casecards(output_dir: Path):
    """
 Fig5B Case Cards: CondFreq
    """
    print("\n=== Plotting Fig5B Case Cards ===")
    
    case_path = output_dir / "fig5b_correction_casecards.json"
    with open(case_path) as f:
        cases = json.load(f)
    
    # Select top 5 representative cases
    cases = cases[:5]
    
    fig = plt.figure(figsize=(7.2, 2.0))
    ax = fig.add_axes([0.02, 0.05, 0.96, 0.85])
    ax.axis('off')
    
    rows = []
    for c in cases:
        metal = c.get('metal', '')
        cn = c.get('cn', '')
        target = c.get('target_donor', c.get('target', ''))
        cond = c.get('condfreq_prediction', c.get('condfreq_top1', ''))
        model = c.get('model_prediction', c.get('model_top1', ''))
        conf = c.get('model_confidence', 0)
        
        rows.append([
            f"{metal} CN={cn}",
            target,
            f"{cond} (wrong)",
            f"{model} (correct, {conf:.0%})"
        ])
    
    tbl = ax.table(
        cellText=rows,
        colLabels=['Group', 'Target', 'CondFreq', 'Model'],
        loc='center',
        cellLoc='center',
        colColours=['#E8E8E8'] * 4,
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(10)
    tbl.scale(1.0, 1.5)
    
    for (row, col), cell in tbl.get_celld().items():
        if row == 0:
            cell.set_text_props(fontweight='bold')
        cell.set_edgecolor('#CCCCCC')
 if col == 2 and row > 0: # CondFreq
            cell.set_text_props(color=COLORS['negative'])
 if col == 3 and row > 0: # Model
            cell.set_text_props(color=COLORS['positive'])
    
    ax.set_title('Model Corrects CondFreq Mistakes (1,899 cases total)', 
                 fontsize=11, fontweight='bold', y=1.02)
    
    fig.savefig(output_dir / "fig5b_casecards.pdf", bbox_inches='tight', dpi=150)
    fig.savefig(output_dir / "fig5b_casecards.svg", bbox_inches='tight')
    plt.close(fig)
    
    print(f"Saved: {output_dir / 'fig5b_casecards.pdf'}")


def plot_fig5c(output_dir: Path):
    """
 Fig5C: Confidence-Coverage + Donor-type Δ
    """
    print("\n=== Plotting Fig5C ===")
    
    cc_path = output_dir / "fig5_confidence_coverage.csv"
    type_path = output_dir / "fig5c_gain_with_delta.csv"
    
    cc = pd.read_csv(cc_path)
    tp = pd.read_csv(type_path)
    
    fig = plt.figure(figsize=(7.2, 4.5))
    
    ax1 = fig.add_axes([0.12, 0.58, 0.83, 0.36])
    
    coverages = cc['coverage'].values * 100
    accuracies = cc['accuracy'].values * 100
    
    ax1.plot(coverages, accuracies, 'o-', color=COLORS['model'], 
             linewidth=2, markersize=4, label='Model')
    ax1.axhline(y=85.7, color=COLORS['neutral'], linestyle='--', 
                alpha=0.7, label='Overall Top-1 (85.7%)')
    
    # Annotate key points (threshold=0.85)
    idx_85 = (cc['threshold'] - 0.85).abs().idxmin()
    ax1.scatter([coverages[idx_85]], [accuracies[idx_85]], 
                s=100, c=COLORS['positive'], zorder=5, edgecolor='white', linewidth=2)
    ax1.annotate(f'τ=0.85\n{accuracies[idx_85]:.1f}%@{coverages[idx_85]:.0f}%',
                 xy=(coverages[idx_85], accuracies[idx_85]),
                 xytext=(coverages[idx_85]-15, accuracies[idx_85]-3),
                 fontsize=9, ha='right',
                 arrowprops=dict(arrowstyle='->', color='gray', lw=0.8))
    
    ax1.set_xlabel('Coverage (%)', fontsize=11)
    ax1.set_ylabel('Top-1 Accuracy (%)', fontsize=11)
    ax1.set_title('Human-in-the-Loop Deployment', fontsize=12, fontweight='bold')
    ax1.set_xlim(55, 105)
    ax1.set_ylim(82, 100)
    ax1.legend(loc='lower left', frameon=False)
    ax1.grid(True, alpha=0.3)
    
    ax2 = fig.add_axes([0.12, 0.12, 0.83, 0.36])
    
    tp_sorted = tp.sort_values('delta_top1', ascending=True)
    
    colors = [COLORS['positive'] if d > 0 else COLORS['negative'] 
              for d in tp_sorted['delta_top1'].values]
    
    bars = ax2.barh(tp_sorted['donor'].astype(str), 
                    tp_sorted['delta_top1'].values * 100, 
                    color=colors, alpha=0.8, edgecolor='white')
    
    ax2.axvline(0, color='black', linewidth=1)
    ax2.set_xlabel(r'$\Delta$Top-1 vs CondFreq (pp)', fontsize=11)
    ax2.set_title('Where Model Improves Over Lookup', fontsize=12, fontweight='bold')
    
    # Add value labels
    for bar, val in zip(bars, tp_sorted['delta_top1'].values * 100):
        if val >= 0:
            ax2.text(val + 1, bar.get_y() + bar.get_height()/2, 
                     f'+{val:.1f}', va='center', fontsize=8)
        else:
            ax2.text(val - 1, bar.get_y() + bar.get_height()/2, 
                     f'{val:.1f}', va='center', ha='right', fontsize=8)
    
    fig.savefig(output_dir / "fig5c_confidence_and_type.pdf", bbox_inches='tight', dpi=150)
    fig.savefig(output_dir / "fig5c_confidence_and_type.svg", bbox_inches='tight')
    fig.savefig(output_dir / "fig5c_confidence_and_type.png", bbox_inches='tight', dpi=150)
    plt.close(fig)
    
    print(f"Saved: {output_dir / 'fig5c_confidence_and_type.pdf'}")


def plot_fig5d(output_dir: Path):
    """
 Fig5D: Valid Parse Uplift + CN
    """
    print("\n=== Plotting Fig5D ===")
    
    uplift_path = output_dir / "fig5d_validity_uplift.csv"
    cn_path = output_dir / "fig5d_cn_stratified.csv"
    
    up = pd.read_csv(uplift_path)
    cn = pd.read_csv(cn_path)
    
    fig = plt.figure(figsize=(7.2, 3.5))
    
    ax1 = fig.add_axes([0.10, 0.18, 0.52, 0.72])
    
    up_types = up[up['corruption_type'] != 'overall'].copy()
    
    x = np.arange(len(up_types))
    w = 0.35
    
    bars1 = ax1.bar(x - w/2, up_types['valid_after_corruption'].values * 100, 
                    width=w, label='Before Repair', color=COLORS['negative'], alpha=0.7)
    bars2 = ax1.bar(x + w/2, up_types['valid_after_repair'].values * 100, 
                    width=w, label='After Repair', color=COLORS['positive'], alpha=0.7)
    
    # Simplify labels
    labels = up_types['corruption_type'].str.replace('_', '\n').values
    ax1.set_xticks(x)
    ax1.set_xticklabels(labels, rotation=0, ha='center', fontsize=8)
    ax1.set_ylabel('Valid Parse Rate (%)', fontsize=11)
    ax1.set_title('Automated Repair (Validity Uplift)', fontsize=12, fontweight='bold')
    ax1.legend(loc='upper right', frameon=False, fontsize=9)
    ax1.set_ylim(0, 115)
    
    # Add uplift annotations
    for i, (b1, b2) in enumerate(zip(bars1, bars2)):
        uplift = up_types.iloc[i]['uplift'] * 100
        if uplift > 5:
            ax1.annotate(f'+{uplift:.0f}pp', 
                         xy=(b2.get_x() + b2.get_width()/2, b2.get_height()),
                         xytext=(0, 3), textcoords='offset points',
                         ha='center', fontsize=7, color=COLORS['positive'])
    
    ax2 = fig.add_axes([0.70, 0.18, 0.28, 0.72])
    
    cn_main = cn[(cn['cn'] >= 4) & (cn['cn'] <= 8) & (cn['n_samples'] >= 100)]
    
    bars = ax2.bar(cn_main['cn'].astype(str), cn_main['top1_acc'].values * 100,
                   color=COLORS['model'], alpha=0.8, edgecolor='white')
    
    ax2.set_ylim(90, 102)
    ax2.set_xlabel('Coordination Number', fontsize=10)
    ax2.set_ylabel('Recovery (%)', fontsize=10)
    ax2.set_title('Stable Across CN', fontsize=11, fontweight='bold')
    ax2.axhline(y=97.8, color=COLORS['neutral'], linestyle='--', alpha=0.7)
    
    # Add value labels
    for bar in bars:
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3,
                 f'{bar.get_height():.1f}', ha='center', fontsize=8)
    
    fig.savefig(output_dir / "fig5d_uplift_and_cn.pdf", bbox_inches='tight', dpi=150)
    fig.savefig(output_dir / "fig5d_uplift_and_cn.svg", bbox_inches='tight')
    fig.savefig(output_dir / "fig5d_uplift_and_cn.png", bbox_inches='tight', dpi=150)
    plt.close(fig)
    
    print(f"Saved: {output_dir / 'fig5d_uplift_and_cn.pdf'}")


def plot_fig5_combined(output_dir: Path):
    """
 Combine all panels Fig5
    """
    print("\n=== Plotting Combined Fig5 ===")
    
    # Load all data
    delta_df = pd.read_csv(output_dir / "fig5_delta_top1_by_metal_cn.csv")
    gains_df = pd.read_csv(output_dir / "fig5_delta_top1_gains.csv")
    cc_df = pd.read_csv(output_dir / "fig5_confidence_coverage.csv")
    type_df = pd.read_csv(output_dir / "fig5c_gain_with_delta.csv")
    uplift_df = pd.read_csv(output_dir / "fig5d_validity_uplift.csv")
    cn_df = pd.read_csv(output_dir / "fig5d_cn_stratified.csv")
    
    fig = plt.figure(figsize=(10, 8))
    
    ax_b = fig.add_axes([0.08, 0.55, 0.40, 0.38])
    deltas = delta_df['delta_top1'].values * 100
    bins = np.linspace(-40, 45, 18)
    n, bins_out, patches = ax_b.hist(deltas, bins=bins, edgecolor='white', linewidth=0.5)
    for patch, left, right in zip(patches, bins_out[:-1], bins_out[1:]):
        center = (left + right) / 2
        patch.set_facecolor(COLORS['positive'] if center >= 0 else COLORS['negative'])
        patch.set_alpha(0.7)
    ax_b.axvline(0, color='black', linewidth=1.5)
    ax_b.axvline(4.26, color=COLORS['model'], linewidth=2, linestyle='--')
    ax_b.set_xlabel(r'$\Delta$Top-1 (pp)', fontsize=10)
    ax_b.set_ylabel('Count', fontsize=10)
    ax_b.set_title('(B) Beyond Conditional Lookup', fontsize=11, fontweight='bold')
    ax_b.text(0.97, 0.95, f'Model wins: 67%\nMean: +4.26pp',
              transform=ax_b.transAxes, ha='right', va='top', fontsize=8,
              bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    ax_c1 = fig.add_axes([0.58, 0.55, 0.38, 0.38])
    coverages = cc_df['coverage'].values * 100
    accuracies = cc_df['accuracy'].values * 100
    ax_c1.plot(coverages, accuracies, 'o-', color=COLORS['model'], linewidth=2, markersize=3)
    ax_c1.axhline(y=85.7, color=COLORS['neutral'], linestyle='--', alpha=0.7)
    idx_85 = (cc_df['threshold'] - 0.85).abs().idxmin()
    ax_c1.scatter([coverages[idx_85]], [accuracies[idx_85]], s=80, c=COLORS['positive'], 
                  zorder=5, edgecolor='white', linewidth=2)
    ax_c1.set_xlabel('Coverage (%)', fontsize=10)
    ax_c1.set_ylabel('Accuracy (%)', fontsize=10)
    ax_c1.set_title('(C) Confidence-Coverage Trade-off', fontsize=11, fontweight='bold')
    ax_c1.set_xlim(55, 105)
    ax_c1.set_ylim(82, 100)
    ax_c1.grid(True, alpha=0.3)
    ax_c1.text(0.03, 0.03, f'95.6%@74.3% (τ=0.85)',
               transform=ax_c1.transAxes, fontsize=8)
    
    ax_d1 = fig.add_axes([0.08, 0.08, 0.40, 0.38])
    up_types = uplift_df[uplift_df['corruption_type'] != 'overall'].copy()
    x = np.arange(len(up_types))
    w = 0.35
    ax_d1.bar(x - w/2, up_types['valid_after_corruption'].values * 100, 
              width=w, label='Before', color=COLORS['negative'], alpha=0.7)
    ax_d1.bar(x + w/2, up_types['valid_after_repair'].values * 100, 
              width=w, label='After', color=COLORS['positive'], alpha=0.7)
    labels = ['bracket\nswap', 'delimiter\nswap', 'missing\nbracket', 
              'missing\ncharge', 'truncation']
    ax_d1.set_xticks(x)
    ax_d1.set_xticklabels(labels, fontsize=7)
    ax_d1.set_ylabel('Valid Rate (%)', fontsize=10)
    ax_d1.set_title('(D) Repair Uplift (+26.8pp overall)', fontsize=11, fontweight='bold')
    ax_d1.legend(loc='upper right', frameon=False, fontsize=8)
    ax_d1.set_ylim(0, 115)
    
    ax_d2 = fig.add_axes([0.58, 0.08, 0.38, 0.38])
    cn_main = cn_df[(cn_df['cn'] >= 4) & (cn_df['cn'] <= 8) & (cn_df['n_samples'] >= 100)]
    bars = ax_d2.bar(cn_main['cn'].astype(str), cn_main['top1_acc'].values * 100,
                     color=COLORS['model'], alpha=0.8, edgecolor='white')
    ax_d2.set_ylim(90, 102)
    ax_d2.set_xlabel('Coordination Number', fontsize=10)
    ax_d2.set_ylabel('Recovery (%)', fontsize=10)
    ax_d2.set_title('(E) Structure Recovery by CN', fontsize=11, fontweight='bold')
    ax_d2.axhline(y=97.8, color=COLORS['neutral'], linestyle='--', alpha=0.7)
    for bar in bars:
        ax_d2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3,
                   f'{bar.get_height():.1f}', ha='center', fontsize=7)
    
    fig.savefig(output_dir / "Fig5_tools_combined.pdf", bbox_inches='tight', dpi=150)
    fig.savefig(output_dir / "Fig5_tools_combined.svg", bbox_inches='tight')
    fig.savefig(output_dir / "Fig5_tools_combined.png", bbox_inches='tight', dpi=150)
    plt.close(fig)
    
    print(f"Saved: {output_dir / 'Fig5_tools_combined.pdf'}")


def main():
    output_dir = Path('/data/CoordRep/CoordSMILES/libcoordrep/outputs/fig5')
    
    plot_fig5b(output_dir)
    plot_fig5b_casecards(output_dir)
    plot_fig5c(output_dir)
    plot_fig5d(output_dir)
    
    plot_fig5_combined(output_dir)
    
    print("\n" + "=" * 50)
    print("All Fig5 Panels Generated!")
    print("=" * 50)


if __name__ == '__main__':
    main()
