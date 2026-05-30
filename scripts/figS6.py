#!/usr/bin/env python
"""Figure S6: Sensitivity and replicate analysis"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 11

BASE = os.environ.get('BASE', '.')
OUT = f'{BASE}/figures/paper_v2'

cauchy = pd.read_csv(f'{BASE}/results/all_organs_cauchy_AD.csv.gz')

# ---- Panel A: Replicate barplot ----
# For each organ, min p_cauchy for CTRL_1 and CTRL_2
rep_data = cauchy.groupby(['organ', 'sample'])['p_cauchy'].min().reset_index()
rep_pivot = rep_data.pivot(index='organ', columns='sample', values='p_cauchy')
# Some organs may not have CTRL_2
for col in ['CTRL_1', 'CTRL_2']:
    if col not in rep_pivot.columns:
        rep_pivot[col] = np.nan

organs = sorted(rep_pivot.index)
rep_pivot = rep_pivot.loc[organs]

# ---- Panel B: Replicate scatter ----
# Per annotation, -log10(p_cauchy) for CTRL_1 vs CTRL_2
ctrl1 = cauchy[cauchy['sample'] == 'CTRL_1'][['annotation', 'organ', 'p_cauchy']].rename(columns={'p_cauchy': 'p1'})
ctrl2 = cauchy[cauchy['sample'] == 'CTRL_2'][['annotation', 'organ', 'p_cauchy']].rename(columns={'p_cauchy': 'p2'})
merged = ctrl1.merge(ctrl2, on=['annotation', 'organ'], how='inner')
merged['logp1'] = -np.log10(merged['p1'].clip(1e-300))
merged['logp2'] = -np.log10(merged['p2'].clip(1e-300))

# ---- Panel C: Volcano-style annotation plot ----
ctrl1_all = cauchy[cauchy['sample'] == 'CTRL_1'].copy()
ctrl1_all['logp_cauchy'] = -np.log10(ctrl1_all['p_cauchy'].clip(1e-300))
ctrl1_all['logp_median'] = -np.log10(ctrl1_all['p_median'].clip(1e-300))

bonf = -np.log10(0.05 / len(ctrl1_all))

# ---- PLOT ----
fig = plt.figure(figsize=(16, 16))

# Panel A: top-left
ax1 = fig.add_axes([0.08, 0.58, 0.4, 0.38])
x_pos = np.arange(len(organs))
width = 0.35
v1 = -np.log10(rep_pivot['CTRL_1'].values.clip(1e-300))
v2 = -np.log10(rep_pivot['CTRL_2'].values.clip(1e-300))
ax1.bar(x_pos - width/2, v1, width, color='#4C72B0', label='CTRL_1', edgecolor='white', linewidth=0.5)
ax1.bar(x_pos + width/2, v2, width, color='#DD8452', label='CTRL_2', edgecolor='white', linewidth=0.5)
ax1.set_xticks(x_pos)
ax1.set_xticklabels(organs, rotation=45, ha='right', fontsize=9)
ax1.set_ylabel('-log10(min p_cauchy)', fontsize=11)
ax1.set_title('Replicate concordance (min p per organ)', fontsize=12, fontweight='bold')
ax1.legend(fontsize=10)
ax1.text(-0.08, 1.02, 'A', transform=ax1.transAxes, fontsize=16, fontweight='bold')

# Panel B: top-right
ax2 = fig.add_axes([0.58, 0.58, 0.36, 0.38])
# Color by organ
organ_list = sorted(merged['organ'].unique())
colors = plt.cm.tab20(np.linspace(0, 1, len(organ_list)))
organ_color = {o: colors[i] for i, o in enumerate(organ_list)}

for o in organ_list:
    sub = merged[merged['organ'] == o]
    ax2.scatter(sub['logp1'], sub['logp2'], c=[organ_color[o]], label=o, s=30, alpha=0.8, edgecolors='white', linewidth=0.3)

# Pearson r
valid = merged.dropna()
if len(valid) > 2:
    r, pval = pearsonr(valid['logp1'], valid['logp2'])
    ax2.text(0.05, 0.95, f'r = {r:.3f}', transform=ax2.transAxes, fontsize=11, va='top')

# Label top 5
top5 = merged.nlargest(5, 'logp1')
for _, row in top5.iterrows():
    ax2.annotate(row['annotation'], (row['logp1'], row['logp2']), fontsize=7,
                 xytext=(5, 5), textcoords='offset points')

lim = max(merged['logp1'].max(), merged['logp2'].max()) * 1.1
ax2.plot([0, lim], [0, lim], 'k--', alpha=0.3, linewidth=0.8)
ax2.set_xlabel('-log10(p_cauchy) CTRL_1', fontsize=11)
ax2.set_ylabel('-log10(p_cauchy) CTRL_2', fontsize=11)
ax2.set_title('Annotation-level replicate correlation', fontsize=12, fontweight='bold')
ax2.legend(fontsize=6, ncol=2, loc='lower right', framealpha=0.7)
ax2.text(-0.08, 1.02, 'B', transform=ax2.transAxes, fontsize=16, fontweight='bold')

# Panel C: bottom wide
ax3 = fig.add_axes([0.08, 0.08, 0.85, 0.38])
for o in organ_list:
    sub = ctrl1_all[ctrl1_all['organ'] == o]
    ax3.scatter(sub['logp_median'], sub['logp_cauchy'], c=[organ_color[o]], label=o, s=40, alpha=0.8, edgecolors='white', linewidth=0.3)

ax3.axhline(bonf, color='red', linestyle='--', linewidth=1, alpha=0.7, label=f'Bonferroni (p={0.05/len(ctrl1_all):.1e})')
ax3.set_xlabel('-log10(p_median)', fontsize=11)
ax3.set_ylabel('-log10(p_cauchy)', fontsize=11)
ax3.set_title('Annotation enrichment: Cauchy vs Median significance (CTRL_1)', fontsize=12, fontweight='bold')

# Label top 5
top5c = ctrl1_all.nlargest(5, 'logp_cauchy')
for _, row in top5c.iterrows():
    ax3.annotate(f"{row['annotation']}\n({row['organ']})", (row['logp_median'], row['logp_cauchy']),
                 fontsize=7, xytext=(8, 4), textcoords='offset points',
                 arrowprops=dict(arrowstyle='->', color='gray', lw=0.5))

ax3.legend(fontsize=7, ncol=3, loc='upper left', framealpha=0.7)
ax3.text(-0.04, 1.02, 'C', transform=ax3.transAxes, fontsize=16, fontweight='bold')

plt.savefig(f'{OUT}/FigureS6.png', dpi=200, bbox_inches='tight')
plt.savefig(f'{OUT}/FigureS6.pdf', dpi=200, bbox_inches='tight')
plt.close()
print("FigureS6 saved")
for ext in ['png', 'pdf']:
    sz = os.path.getsize(f'{OUT}/FigureS6.{ext}')
    print(f"  FigureS6.{ext}: {sz/1024:.0f} KB")
