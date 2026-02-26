#!/usr/bin/env python3
"""
Reproduce Fig5 B/C/D/E from CSV files

Only requires CSV files in this directory:
- panel_b_delta_histogram.csv
- panel_b_top_gains.csv
- panel_b_summary.csv
- panel_c_confidence_coverage.csv
- panel_c_donor_type_delta.csv
- panel_d_validity_uplift.csv
- panel_e_cn_stratified.csv

Run: python plot_fig5_from_csv.py
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Set plot style
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 11
plt.rcParams['axes.titlesize'] = 12

# Color scheme
COLORS = {
    'model': '#2E86AB',
    'baseline': '#A23B72',
    'positive': '#28A745',
    'negative': '#DC3545',
    'neutral': '#6C757D',
}


def plot_panel_b(data_dir: Path, output_dir: Path):
 """Panel B: ΔTop1 + Top gains"""
    hist = pd.read_csv(data_dir / 'panel_b_delta_histogram.csv')
    gains = pd.read_csv(data_dir / 'panel_b_top_gains.csv')
    summary = pd.read_csv(data_dir / 'panel_b_summary.csv').iloc[0]
    
    fig = plt.figure(figsize=(7.2, 3.2))
    
    ax = fig.add_axes([0.10, 0.18, 0.52, 0.72])
    
    colors = [COLORS['positive'] if p else COLORS['negative'] 
              for p in hist['is_positive']]
    
    ax.bar(hist['bin_center'], hist['count'], width=4.5, 
           color=colors, alpha=0.7, edgecolor='white')
    ax.axvline(0, color='black', linewidth=1.5)
    ax.axvline(summary['mean_delta_top1_pp'], color=COLORS['model'], 
               linewidth=2, linestyle='--', label=f'Mean +{summary["mean_delta_top1_pp"]:.2f}pp')
    
    ax.set_xlabel(r'$\Delta$Top-1 (Model − CondFreq), pp', fontsize=11)
    ax.set_ylabel('Count of (metal, CN) groups', fontsize=11)
    ax.set_title('(B) Beyond Conditional Lookup', fontsize=12, fontweight='bold')
    ax.legend(loc='upper left', frameon=False)
    
    ax.text(0.97, 0.95, 
            f'Model wins: {summary["model_wins"]}/{summary["total_groups"]} ({summary["model_wins"]/summary["total_groups"]*100:.0f}%)',
            transform=ax.transAxes, ha='right', va='top', fontsize=9,
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    ax2 = fig.add_axes([0.68, 0.18, 0.30, 0.72])
    ax2.axis('off')
    
    top5 = gains.head(5)
    cell_text = [[r['metal'], str(int(r['cn'])), f"+{r['delta_top1_pp']:.1f}"] 
                 for _, r in top5.iterrows()]
    
    tbl = ax2.table(cellText=cell_text, colLabels=['Metal', 'CN', 'Δpp'],
                    loc='center', cellLoc='center', colColours=['#E8E8E8']*3)
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(10)
    tbl.scale(1.0, 1.4)
    ax2.set_title('Top Gains', fontsize=11, fontweight='bold', pad=10)
    
    fig.savefig(output_dir / 'panel_b.pdf', bbox_inches='tight', dpi=150)
    fig.savefig(output_dir / 'panel_b.png', bbox_inches='tight', dpi=150)
    plt.close(fig)
    print(f"Saved: panel_b.pdf/png")


def plot_panel_c(data_dir: Path, output_dir: Path):
    """Panel C: Confidence-Coverage + Donor Δ"""
    cc = pd.read_csv(data_dir / 'panel_c_confidence_coverage.csv')
    donor = pd.read_csv(data_dir / 'panel_c_donor_type_delta.csv')
    
    # Deduplicate
    cc = cc.drop_duplicates(subset=['threshold'])
    
    fig = plt.figure(figsize=(7.2, 4.5))
    
    ax1 = fig.add_axes([0.12, 0.58, 0.83, 0.36])
    
    ax1.plot(cc['coverage_pct'], cc['accuracy_pct'], 'o-', 
             color=COLORS['model'], linewidth=2, markersize=4)
    ax1.axhline(y=85.71, color=COLORS['neutral'], linestyle='--', 
                alpha=0.7, label='Overall (85.7%)')
    
    # Annotate tau=0.85
    row_85 = cc[cc['threshold'] == 0.85].iloc[0]
    ax1.scatter([row_85['coverage_pct']], [row_85['accuracy_pct']], 
                s=100, c=COLORS['positive'], zorder=5, edgecolor='white', linewidth=2)
    ax1.annotate(f'τ=0.85\n{row_85["accuracy_pct"]:.1f}%@{row_85["coverage_pct"]:.0f}%',
                 xy=(row_85['coverage_pct'], row_85['accuracy_pct']),
                 xytext=(row_85['coverage_pct']-15, row_85['accuracy_pct']-3),
                 fontsize=9, ha='right')
    
    ax1.set_xlabel('Coverage (%)', fontsize=11)
    ax1.set_ylabel('Top-1 Accuracy (%)', fontsize=11)
    ax1.set_title('(C) Confidence-Coverage Trade-off', fontsize=12, fontweight='bold')
    ax1.set_xlim(55, 105)
    ax1.set_ylim(82, 100)
    ax1.legend(loc='lower left', frameon=False)
    ax1.grid(True, alpha=0.3)
    
    ax2 = fig.add_axes([0.12, 0.12, 0.83, 0.36])
    
    donor_sorted = donor.sort_values('delta_top1_pp', ascending=True)
    colors = [COLORS['positive'] if d > 0 else COLORS['negative'] 
              for d in donor_sorted['delta_top1_pp']]
    
    bars = ax2.barh(donor_sorted['donor_atom'], donor_sorted['delta_top1_pp'], 
                    color=colors, alpha=0.8)
    ax2.axvline(0, color='black', linewidth=1)
    ax2.set_xlabel(r'$\Delta$Top-1 vs CondFreq (pp)', fontsize=11)
    ax2.set_title('Improvement by Donor Type', fontsize=11, fontweight='bold')
    
    for bar, val in zip(bars, donor_sorted['delta_top1_pp']):
        x_pos = val + 1 if val >= 0 else val - 1
        ha = 'left' if val >= 0 else 'right'
        ax2.text(x_pos, bar.get_y() + bar.get_height()/2, 
                 f'{val:+.1f}', va='center', ha=ha, fontsize=8)
    
    fig.savefig(output_dir / 'panel_c.pdf', bbox_inches='tight', dpi=150)
    fig.savefig(output_dir / 'panel_c.png', bbox_inches='tight', dpi=150)
    plt.close(fig)
    print(f"Saved: panel_c.pdf/png")


def plot_panel_d(data_dir: Path, output_dir: Path):
    """Panel D: Validity Uplift"""
    up = pd.read_csv(data_dir / 'panel_d_validity_uplift.csv')
    
    fig = plt.figure(figsize=(6, 3.5))
    ax = fig.add_axes([0.12, 0.20, 0.85, 0.70])
    
    # Filter out overall
    up_types = up[up['corruption_type'] != 'overall'].copy()
    
    x = np.arange(len(up_types))
    w = 0.35
    
    bars1 = ax.bar(x - w/2, up_types['valid_after_corruption'] * 100, width=w,
                   label='Before Repair', color=COLORS['negative'], alpha=0.7)
    bars2 = ax.bar(x + w/2, up_types['valid_after_repair'] * 100, width=w,
                   label='After Repair', color=COLORS['positive'], alpha=0.7)
    
    # Simplify labels
    labels = [s.replace('_', '\n') for s in up_types['corruption_type']]
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel('Valid Parse Rate (%)', fontsize=11)
    ax.set_title('(D) Repair Validity Uplift (+26.8pp overall)', fontsize=12, fontweight='bold')
    ax.legend(loc='upper right', frameon=False, fontsize=9)
    ax.set_ylim(0, 115)
    
    # Add uplift annotations
    for i, (b1, b2) in enumerate(zip(bars1, bars2)):
        uplift = up_types.iloc[i]['uplift'] * 100
        if uplift > 5:
            ax.annotate(f'+{uplift:.0f}pp', 
                        xy=(b2.get_x() + b2.get_width()/2, b2.get_height()),
                        xytext=(0, 3), textcoords='offset points',
                        ha='center', fontsize=7, color=COLORS['positive'])
    
    fig.savefig(output_dir / 'panel_d.pdf', bbox_inches='tight', dpi=150)
    fig.savefig(output_dir / 'panel_d.png', bbox_inches='tight', dpi=150)
    plt.close(fig)
    print(f"Saved: panel_d.pdf/png")


def plot_panel_e(data_dir: Path, output_dir: Path):
    """Panel E: CN Stratified Recovery"""
    cn = pd.read_csv(data_dir / 'panel_e_cn_stratified.csv')
    
    fig = plt.figure(figsize=(5, 3.5))
    ax = fig.add_axes([0.15, 0.18, 0.80, 0.72])
    
    cn_main = cn[(cn['cn'] >= 4) & (cn['cn'] <= 8) & (cn['n_samples'] >= 100)]
    
    bars = ax.bar(cn_main['cn'].astype(str), cn_main['top1_acc_pct'],
                  color=COLORS['model'], alpha=0.8, edgecolor='white')
    
    ax.axhline(y=97.79, color=COLORS['neutral'], linestyle='--', alpha=0.7,
               label='Overall (97.8%)')
    ax.set_ylim(90, 102)
    ax.set_xlabel('Coordination Number', fontsize=11)
    ax.set_ylabel('Top-1 Recovery (%)', fontsize=11)
    ax.set_title('(E) Structure Recovery by CN', fontsize=12, fontweight='bold')
    ax.legend(loc='lower right', frameon=False, fontsize=9)
    
    for bar in bars:
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3,
                f'{bar.get_height():.1f}', ha='center', fontsize=8)
    
    fig.savefig(output_dir / 'panel_e.pdf', bbox_inches='tight', dpi=150)
    fig.savefig(output_dir / 'panel_e.png', bbox_inches='tight', dpi=150)
    plt.close(fig)
    print(f"Saved: panel_e.pdf/png")


def plot_combined(data_dir: Path, output_dir: Path):
    """Combine all panels"""
    # Load all data
    hist = pd.read_csv(data_dir / 'panel_b_delta_histogram.csv')
    summary = pd.read_csv(data_dir / 'panel_b_summary.csv').iloc[0]
    cc = pd.read_csv(data_dir / 'panel_c_confidence_coverage.csv').drop_duplicates(subset=['threshold'])
    up = pd.read_csv(data_dir / 'panel_d_validity_uplift.csv')
    cn = pd.read_csv(data_dir / 'panel_e_cn_stratified.csv')
    
    fig = plt.figure(figsize=(10, 8))
    
    ax_b = fig.add_axes([0.08, 0.55, 0.40, 0.38])
    colors_b = [COLORS['positive'] if p else COLORS['negative'] for p in hist['is_positive']]
    ax_b.bar(hist['bin_center'], hist['count'], width=4.5, color=colors_b, alpha=0.7, edgecolor='white')
    ax_b.axvline(0, color='black', linewidth=1.5)
    ax_b.axvline(summary['mean_delta_top1_pp'], color=COLORS['model'], linewidth=2, linestyle='--')
    ax_b.set_xlabel(r'$\Delta$Top-1 (pp)', fontsize=10)
    ax_b.set_ylabel('Count', fontsize=10)
    ax_b.set_title('(B) Beyond Conditional Lookup', fontsize=11, fontweight='bold')
    ax_b.text(0.97, 0.95, f'Model wins: 67%\nMean: +4.26pp',
              transform=ax_b.transAxes, ha='right', va='top', fontsize=8,
              bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    ax_c = fig.add_axes([0.58, 0.55, 0.38, 0.38])
    ax_c.plot(cc['coverage_pct'], cc['accuracy_pct'], 'o-', color=COLORS['model'], linewidth=2, markersize=3)
    ax_c.axhline(y=85.71, color=COLORS['neutral'], linestyle='--', alpha=0.7)
    row_85 = cc[cc['threshold'] == 0.85].iloc[0]
    ax_c.scatter([row_85['coverage_pct']], [row_85['accuracy_pct']], s=80, c=COLORS['positive'], 
                 zorder=5, edgecolor='white', linewidth=2)
    ax_c.set_xlabel('Coverage (%)', fontsize=10)
    ax_c.set_ylabel('Accuracy (%)', fontsize=10)
    ax_c.set_title('(C) Confidence-Coverage', fontsize=11, fontweight='bold')
    ax_c.set_xlim(55, 105)
    ax_c.set_ylim(82, 100)
    ax_c.grid(True, alpha=0.3)
    ax_c.text(0.03, 0.03, '95.6%@74.3% (τ=0.85)', transform=ax_c.transAxes, fontsize=8)
    
    ax_d = fig.add_axes([0.08, 0.08, 0.40, 0.38])
    up_types = up[up['corruption_type'] != 'overall']
    x = np.arange(len(up_types))
    w = 0.35
    ax_d.bar(x - w/2, up_types['valid_after_corruption'] * 100, width=w,
             label='Before', color=COLORS['negative'], alpha=0.7)
    ax_d.bar(x + w/2, up_types['valid_after_repair'] * 100, width=w,
             label='After', color=COLORS['positive'], alpha=0.7)
    labels = ['bracket\nswap', 'delimiter\nswap', 'missing\nbracket', 'missing\ncharge', 'truncation']
    ax_d.set_xticks(x)
    ax_d.set_xticklabels(labels, fontsize=7)
    ax_d.set_ylabel('Valid Rate (%)', fontsize=10)
    ax_d.set_title('(D) Repair Uplift (+26.8pp)', fontsize=11, fontweight='bold')
    ax_d.legend(loc='upper right', frameon=False, fontsize=8)
    ax_d.set_ylim(0, 115)
    
    ax_e = fig.add_axes([0.58, 0.08, 0.38, 0.38])
    cn_main = cn[(cn['cn'] >= 4) & (cn['cn'] <= 8) & (cn['n_samples'] >= 100)]
    bars = ax_e.bar(cn_main['cn'].astype(str), cn_main['top1_acc_pct'],
                    color=COLORS['model'], alpha=0.8, edgecolor='white')
    ax_e.axhline(y=97.79, color=COLORS['neutral'], linestyle='--', alpha=0.7)
    ax_e.set_ylim(90, 102)
    ax_e.set_xlabel('Coordination Number', fontsize=10)
    ax_e.set_ylabel('Recovery (%)', fontsize=10)
    ax_e.set_title('(E) Structure Recovery by CN', fontsize=11, fontweight='bold')
    for bar in bars:
        ax_e.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3,
                  f'{bar.get_height():.1f}', ha='center', fontsize=7)
    
    fig.savefig(output_dir / 'Fig5_combined.pdf', bbox_inches='tight', dpi=150)
    fig.savefig(output_dir / 'Fig5_combined.png', bbox_inches='tight', dpi=150)
    plt.close(fig)
    print(f"Saved: Fig5_combined.pdf/png")


def main():
    data_dir = Path(__file__).parent
    output_dir = data_dir / 'plots'
    output_dir.mkdir(exist_ok=True)
    
    print("=" * 50)
    print("Reproducing Fig5 from CSV data")
    print("=" * 50)
    
    plot_panel_b(data_dir, output_dir)
    plot_panel_c(data_dir, output_dir)
    plot_panel_d(data_dir, output_dir)
    plot_panel_e(data_dir, output_dir)
    plot_combined(data_dir, output_dir)
    
    print("\n" + "=" * 50)
    print(f"All plots saved to: {output_dir}")
    print("=" * 50)


if __name__ == '__main__':
    main()
