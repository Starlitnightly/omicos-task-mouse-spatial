#!/usr/bin/env python
"""Figure 1: Whole-body AD genetic risk atlas for Cell paper."""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import pandas as pd
import numpy as np
from scipy import stats

plt.rcParams.update({'font.family': 'Arial', 'font.size': 11})

# ---------- Load data ----------
print("Loading spot data...")
spot = pd.read_csv('results/all_organs_spot_AD_pvalues.csv.gz')
spot_ctrl1 = spot[spot['sample'] == 'CTRL_1'].copy()
print(f"  CTRL_1 spots: {len(spot_ctrl1):,}")

print("Loading cauchy data...")
cauchy = pd.read_csv('results/all_organs_cauchy_AD.csv.gz')

# ---------- Setup ----------
organs = sorted(spot_ctrl1['organ'].unique())
n_organs = len(organs)
print(f"  {n_organs} organs")

# 16 distinct colors
tab20 = plt.cm.tab20(np.linspace(0, 1, 20))
organ_colors = {org: tab20[i] for i, org in enumerate(organs)}

fig, axes = plt.subplots(2, 2, figsize=(16, 14))

# ========== Panel A: Organ annotation map ==========
ax = axes[0, 0]
for org in organs:
    mask = spot_ctrl1['organ'] == org
    ax.scatter(spot_ctrl1.loc[mask, 'x'], spot_ctrl1.loc[mask, 'y'],
               c=[organ_colors[org]], s=0.1, label=org, rasterized=True)
ax.legend(fontsize=8, markerscale=10, loc='center left', bbox_to_anchor=(1.02, 0.5),
          frameon=False, handletextpad=0.3)
ax.set_title('Organ annotation', fontsize=14, fontweight='bold')
ax.set_aspect('equal')
ax.set_xticks([])
ax.set_yticks([])
for spine in ax.spines.values():
    spine.set_visible(False)
ax.text(-0.05, 1.05, 'A', transform=ax.transAxes, fontsize=20, fontweight='bold', va='top')

# ========== Panel B: AD risk map ==========
ax = axes[0, 1]
vmax = np.percentile(spot_ctrl1['logp'], 99)
sc = ax.scatter(spot_ctrl1['x'], spot_ctrl1['y'], c=spot_ctrl1['logp'],
                cmap='RdYlBu_r', s=0.1, vmin=0, vmax=vmax, rasterized=True)
cbar = fig.colorbar(sc, ax=ax, shrink=0.6, pad=0.02)
cbar.set_label(r'$-\log_{10}(p)$', fontsize=12)
ax.set_title('AD genetic risk', fontsize=14, fontweight='bold')
ax.set_aspect('equal')
ax.set_xticks([])
ax.set_yticks([])
for spine in ax.spines.values():
    spine.set_visible(False)
ax.text(-0.05, 1.05, 'B', transform=ax.transAxes, fontsize=20, fontweight='bold', va='top')

# ========== Panel C: Organ-level enrichment bar chart ==========
ax = axes[1, 0]
cauchy_ctrl1 = cauchy[cauchy['sample'] == 'CTRL_1']
organ_min_p = cauchy_ctrl1.groupby('organ')['p_cauchy'].min().reset_index()
organ_min_p['neg_logp'] = -np.log10(organ_min_p['p_cauchy'])
organ_min_p = organ_min_p.sort_values('neg_logp', ascending=True)

bonf_thresh = -np.log10(0.05 / n_organs)
colors_bar = ['#d62728' if v > bonf_thresh else '#888888' for v in organ_min_p['neg_logp']]

ax.barh(range(len(organ_min_p)), organ_min_p['neg_logp'], color=colors_bar, edgecolor='none')
ax.set_yticks(range(len(organ_min_p)))
ax.set_yticklabels(organ_min_p['organ'], fontsize=10)
ax.axvline(bonf_thresh, color='black', linestyle='--', linewidth=1, label=f'Bonferroni (p=0.05/{n_organs})')
ax.legend(fontsize=9, loc='lower right')
ax.set_xlabel(r'$-\log_{10}$(Cauchy $p$)', fontsize=12)
ax.set_title('Organ-level AD enrichment', fontsize=14, fontweight='bold')
ax.text(-0.05, 1.05, 'C', transform=ax.transAxes, fontsize=20, fontweight='bold', va='top')

# ========== Panel D: Replicate concordance ==========
ax = axes[1, 1]
cauchy_ctrl2 = cauchy[cauchy['sample'] == 'CTRL_2']
organ_min_p1 = cauchy_ctrl1.groupby('organ')['p_cauchy'].min()
organ_min_p2 = cauchy_ctrl2.groupby('organ')['p_cauchy'].min()

# Align on shared organs
shared = sorted(set(organ_min_p1.index) & set(organ_min_p2.index))
x_vals = -np.log10(organ_min_p1[shared].values)
y_vals = -np.log10(organ_min_p2[shared].values)

ax.scatter(x_vals, y_vals, s=40, c='#1f77b4', edgecolors='black', linewidths=0.5, zorder=3)

# Label points
for i, org in enumerate(shared):
    ax.annotate(org, (x_vals[i], y_vals[i]), fontsize=7, ha='left', va='bottom',
                xytext=(4, 4), textcoords='offset points')

# Identity line
lim_max = max(max(x_vals), max(y_vals)) * 1.1
ax.plot([0, lim_max], [0, lim_max], 'k--', linewidth=0.8, alpha=0.5)
ax.set_xlim(0, lim_max)
ax.set_ylim(0, lim_max)

r, pval = stats.pearsonr(x_vals, y_vals)
ax.text(0.05, 0.95, f'Pearson r = {r:.3f}', transform=ax.transAxes,
        fontsize=11, va='top', ha='left',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='gray', alpha=0.8))

ax.set_xlabel(r'CTRL_1: $-\log_{10}$(Cauchy $p$)', fontsize=12)
ax.set_ylabel(r'CTRL_2: $-\log_{10}$(Cauchy $p$)', fontsize=12)
ax.set_title('Replicate concordance', fontsize=14, fontweight='bold')
ax.set_aspect('equal')
ax.text(-0.05, 1.05, 'D', transform=ax.transAxes, fontsize=20, fontweight='bold', va='top')

# ---------- Save ----------
plt.tight_layout()
for fmt in ['png', 'pdf']:
    outpath = f'figures/paper_v2/Figure1.{fmt}'
    fig.savefig(outpath, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Saved: {outpath}")

plt.close()

# Report file sizes
import os
for fmt in ['png', 'pdf']:
    path = f'figures/paper_v2/Figure1.{fmt}'
    size_mb = os.path.getsize(path) / 1e6
    print(f"  {path}: {size_mb:.1f} MB")
