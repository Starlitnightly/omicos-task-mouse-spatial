#!/usr/bin/env python
"""Figure 2: Tissue-resident macrophages universally carry AD risk — Cell paper."""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd
import numpy as np
import os

plt.rcParams.update({'font.family': 'Arial', 'font.size': 11})

# ---------- Load data ----------
print("Loading data...")
cauchy = pd.read_csv('results/all_organs_cauchy_AD.csv.gz')
age = pd.read_csv('results/age_all_cauchy_AD.csv.gz')

# CTRL_1 only for cauchy
cauchy_ctrl1 = cauchy[cauchy['sample_name'].str.contains('CTRL1')].copy()

# ---------- Figure setup ----------
fig = plt.figure(figsize=(16, 16), constrained_layout=True)
gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.08, wspace=0.08,
                       height_ratios=[1, 1.2])

# ========== Panel A: Brain cell type barplot ==========
ax_a = fig.add_subplot(gs[0, 0])

# Use age data for brain cell types — pick 18m as representative (has named cell types)
brain_age = age[age['organ'].str.contains('Brain', na=False)].copy()

# For Panel A, use 18m data which has well-annotated cell types
brain_18m = brain_age[brain_age['age'] == '18m'].copy()
brain_18m['neg_logp'] = -np.log10(brain_18m['p_cauchy'].clip(lower=1e-300))
brain_18m = brain_18m.sort_values('neg_logp', ascending=True)

# Myeloid detection
def is_myeloid(name):
    name_lower = str(name).lower()
    return any(kw in name_lower for kw in ['macrophage', 'microgl', 'myeloid', 'meninges'])

brain_18m['is_myeloid'] = brain_18m['annotation'].apply(is_myeloid)
colors_a = ['#d62728' if m else 'steelblue' for m in brain_18m['is_myeloid']]

n_annot = len(brain_18m)
bonf_a = -np.log10(0.05 / n_annot)

ax_a.barh(range(n_annot), brain_18m['neg_logp'].values, color=colors_a, edgecolor='none')
ax_a.set_yticks(range(n_annot))
ax_a.set_yticklabels(brain_18m['annotation'].values, fontsize=9)
ax_a.axvline(bonf_a, color='black', linestyle='--', linewidth=1,
             label=f'Bonferroni (p=0.05/{n_annot})')
ax_a.legend(fontsize=8, loc='lower right')
ax_a.set_xlabel(r'$-\log_{10}$(Cauchy $p$)', fontsize=12)
ax_a.set_title('Brain cell types', fontsize=14, fontweight='bold')

# Legend patches for myeloid vs other
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor='#d62728', label='Myeloid'),
                   Patch(facecolor='steelblue', label='Other')]
leg = ax_a.legend(handles=legend_elements, fontsize=8, loc='lower right',
                  title='', frameon=True)
# Add threshold line to legend manually
ax_a.axvline(bonf_a, color='black', linestyle='--', linewidth=1)

ax_a.text(-0.05, 1.05, 'A', transform=ax_a.transAxes, fontsize=20,
          fontweight='bold', va='top')

# ========== Panel B: Microglia vs non-myeloid across ages ==========
ax_b = fig.add_subplot(gs[0, 1])

# Normalize age names
age_map = {
    '03_months': 3, '06_months': 6, '12_months': 12, '16_months': 16,
    '18m': 18, '21m': 21, '21_months': 21, '23_months': 23, '24m': 24,
    '30m': 30, '3m': 3, '1m': 1, '6m': 6, '12m': 12,
}
brain_age['age_num'] = brain_age['age'].map(age_map)
brain_age['is_myeloid'] = brain_age['annotation'].apply(is_myeloid)

# For each age: min p for myeloid, median p for non-myeloid
ages_sorted = sorted(brain_age['age_num'].dropna().unique())

microglia_vals = []
nonmyeloid_vals = []
valid_ages = []

for a in ages_sorted:
    sub = brain_age[brain_age['age_num'] == a]
    myeloid_sub = sub[sub['is_myeloid']]
    nonmyeloid_sub = sub[~sub['is_myeloid']]
    if len(myeloid_sub) == 0:
        continue
    # Min p for myeloid (best signal)
    mic_p = myeloid_sub['p_cauchy'].min()
    # Median p for non-myeloid
    nonm_p = nonmyeloid_sub['p_cauchy'].median() if len(nonmyeloid_sub) > 0 else 1.0
    microglia_vals.append(-np.log10(max(mic_p, 1e-300)))
    nonmyeloid_vals.append(-np.log10(max(nonm_p, 1e-300)))
    valid_ages.append(int(a))

x = np.arange(len(valid_ages))
width = 0.35

ax_b.bar(x - width/2, microglia_vals, width, color='#d62728', label='Microglia/macrophage (best)')
ax_b.bar(x + width/2, nonmyeloid_vals, width, color='steelblue', label='Non-myeloid (median)')

# Bonferroni — use typical brain annotation count
typical_n = 18
bonf_b = -np.log10(0.05 / typical_n)
ax_b.axhline(bonf_b, color='black', linestyle='--', linewidth=1,
             label=f'Bonferroni threshold')

ax_b.set_xticks(x)
ax_b.set_xticklabels([f'{a}m' for a in valid_ages], fontsize=10)
ax_b.set_xlabel('Age (months)', fontsize=12)
ax_b.set_ylabel(r'$-\log_{10}$(Cauchy $p$)', fontsize=12)
ax_b.set_title('Microglia vs non-myeloid across ages', fontsize=14, fontweight='bold')
ax_b.legend(fontsize=8, loc='upper left', frameon=True)

ax_b.text(-0.05, 1.05, 'B', transform=ax_b.transAxes, fontsize=20,
          fontweight='bold', va='top')

# ========== Panel C: Top 25 cell type annotations across organs ==========
ax_c = fig.add_subplot(gs[1, :])

# For each organ, get all annotations with -log10(p_cauchy)
cauchy_ctrl1['neg_logp'] = -np.log10(cauchy_ctrl1['p_cauchy'].clip(lower=1e-300))
cauchy_ctrl1['label'] = cauchy_ctrl1['organ'] + ': ' + cauchy_ctrl1['annotation'].str.strip()

# Sort and take top 25
top25 = cauchy_ctrl1.nlargest(25, 'neg_logp').sort_values('neg_logp', ascending=True)

# Color by organ
organs_unique = sorted(top25['organ'].unique())
tab20 = plt.cm.tab20(np.linspace(0, 1, 20))
organ_colors = {org: tab20[i] for i, org in enumerate(sorted(cauchy_ctrl1['organ'].unique()))}

colors_c = [organ_colors[org] for org in top25['organ']]

n_total = len(cauchy_ctrl1)
bonf_c = -np.log10(0.05 / n_total)

ax_c.barh(range(len(top25)), top25['neg_logp'].values, color=colors_c, edgecolor='none')
ax_c.set_yticks(range(len(top25)))
ax_c.set_yticklabels(top25['label'].values, fontsize=9)
ax_c.axvline(bonf_c, color='black', linestyle='--', linewidth=1,
             label=f'Bonferroni (p=0.05/{n_total})')

# Organ legend
legend_patches = [Patch(facecolor=organ_colors[org], label=org) for org in organs_unique]
ax_c.legend(handles=legend_patches, fontsize=8, loc='lower right',
            ncol=2, frameon=True, title='Organ')

ax_c.set_xlabel(r'$-\log_{10}$(Cauchy $p$)', fontsize=12)
ax_c.set_title('Top 25 cell type annotations across organs', fontsize=14, fontweight='bold')

ax_c.text(-0.02, 1.05, 'C', transform=ax_c.transAxes, fontsize=20,
          fontweight='bold', va='top')

# ---------- Save ----------
os.makedirs('figures/paper_v2', exist_ok=True)
for fmt in ['png', 'pdf']:
    outpath = f'figures/paper_v2/Figure2.{fmt}'
    fig.savefig(outpath, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Saved: {outpath}")

plt.close()

# Report file sizes
for fmt in ['png', 'pdf']:
    path = f'figures/paper_v2/Figure2.{fmt}'
    size_mb = os.path.getsize(path) / 1e6
    print(f"  {path}: {size_mb:.1f} MB")
