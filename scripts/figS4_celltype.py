#!/usr/bin/env python
"""Figure S4: Extended cell-type analysis (3 panels)"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
import seaborn as sns
import scanpy as sc
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 11

BASEDIR = ''
OUTDIR = f'{BASEDIR}figures/paper_v2'

short_names = {
    'Bone Marrow': 'BM', 'Brain': 'Brain', 'Brown Fat': 'BrFat',
    'Colon': 'Colon', 'Heart': 'Heart', 'Kidney': 'Kidney',
    'Liver': 'Liver', 'Lung': 'Lung', 'Lymph Node': 'LN',
    'Muscle': 'Muscle', 'Pancreas': 'Panc', 'Skin': 'Skin',
    'Small Intestine': 'SI', 'Spleen': 'Spleen', 'Stomach': 'Stom',
    'Thymus': 'Thymus'
}

# Load cauchy data
cauchy = pd.read_csv(f'{BASEDIR}results/all_organs_cauchy_AD.csv.gz')
ctrl1 = cauchy[cauchy['sample'] == 'CTRL_1'].copy()
ctrl2 = cauchy[cauchy['sample'] == 'CTRL_2'].copy()

ctrl1['logp'] = -np.log10(ctrl1['p_cauchy'].clip(1e-20))
ctrl2['logp'] = -np.log10(ctrl2['p_cauchy'].clip(1e-20))

fig = plt.figure(figsize=(16, 16))
gs = gridspec.GridSpec(2, 2, height_ratios=[1, 0.8], hspace=0.35, wspace=0.3,
                       left=0.08, right=0.95, top=0.95, bottom=0.05)

# Panel A: Heatmap
print("Panel A: Annotation-level heatmap...")
ax_a = fig.add_subplot(gs[0, 0])

pivot = ctrl1.pivot_table(index='annotation', columns='organ', values='logp', aggfunc='first')
pivot = pivot.loc[pivot.max(axis=1).sort_values(ascending=False).index]

if len(pivot) > 30:
    pivot = pivot.head(30)

pivot.index = [a[:30] + '...' if len(a) > 30 else a for a in pivot.index]
pivot.columns = [short_names.get(c, c) for c in pivot.columns]

sns.heatmap(pivot, ax=ax_a, cmap='YlOrRd', annot=True, fmt='.1f',
            annot_kws={'fontsize': 5}, cbar_kws={'label': '-log$_{10}$(Cauchy p)', 'shrink': 0.8},
            linewidths=0.5, linecolor='white')
ax_a.set_title('Annotation-level AD enrichment', fontsize=12)
ax_a.tick_params(axis='y', labelsize=7)
ax_a.tick_params(axis='x', labelsize=8, rotation=45)
ax_a.set_ylabel('')
ax_a.set_xlabel('')
ax_a.text(-0.05, 1.03, 'A', fontsize=20, fontweight='bold',
          va='top', transform=ax_a.transAxes)

# Panel B: Replicate correlation at annotation level
print("Panel B: Replicate correlation...")
ax_b = fig.add_subplot(gs[0, 1])

merged = ctrl1[['annotation', 'organ', 'logp']].merge(
    ctrl2[['annotation', 'organ', 'logp']],
    on=['annotation', 'organ'], suffixes=('_c1', '_c2'), how='inner')

organs_unique = sorted(merged['organ'].unique())
cmap_organs = plt.cm.get_cmap('tab20', len(organs_unique))
organ_colors = {o: cmap_organs(i) for i, o in enumerate(organs_unique)}

for org in organs_unique:
    mask = merged['organ'] == org
    ax_b.scatter(merged.loc[mask, 'logp_c1'], merged.loc[mask, 'logp_c2'],
                 c=[organ_colors[org]], s=20, alpha=0.7,
                 label=short_names.get(org, org), edgecolors='none')

r, pval = stats.pearsonr(merged['logp_c1'], merged['logp_c2'])
maxval = max(merged['logp_c1'].max(), merged['logp_c2'].max()) * 1.1
ax_b.plot([0, maxval], [0, maxval], 'k--', alpha=0.5, linewidth=1)
ax_b.set_xlabel('-log$_{10}$(Cauchy p) CTRL1', fontsize=11)
ax_b.set_ylabel('-log$_{10}$(Cauchy p) CTRL2', fontsize=11)
ax_b.set_title(f'Annotation replicate correlation (r = {r:.3f})', fontsize=12)
ax_b.legend(fontsize=6, ncol=2, markerscale=2, loc='upper left')
ax_b.text(-0.05, 1.03, 'B', fontsize=20, fontweight='bold',
          va='top', transform=ax_b.transAxes)

# Panel C: Top 5 annotations spatial maps
# Highlight the specific annotation region within each organ
print("Panel C: Top 5 spatial maps...")
gs_c = gs[1, :].subgridspec(1, 5, wspace=0.15)

top5 = ctrl1.nlargest(5, 'logp')
print("  Top 5 annotations:")
for _, row in top5.iterrows():
    print(f"    {row['annotation']} ({row['organ']}): p={row['p_cauchy']:.2e}")

# Load h5ad obs for needed organs to get annotation column
needed_organs = top5['organ'].unique()
organ_obs_cache = {}
for org in needed_organs:
    org_file = org.replace(' ', '_')
    f = f'{BASEDIR}data/st/per_organ/{org_file}_CTRL1.h5ad'
    adata = sc.read_h5ad(f, backed='r')
    obs = adata.obs[['x_plotting', 'y_plotting', 'annotation']].copy()
    organ_obs_cache[org] = obs
    adata.file.close()
    print(f"  Loaded {org}: {len(obs)} spots")

for idx, (_, row) in enumerate(top5.iterrows()):
    ax = fig.add_subplot(gs_c[idx])
    organ = row['organ']
    annotation = row['annotation']

    obs = organ_obs_cache[organ]
    highlight = obs['annotation'] == annotation
    bg = ~highlight

    # Plot background in gray
    if bg.sum() > 0:
        obs_bg = obs[bg]
        if len(obs_bg) > 30000:
            obs_bg = obs_bg.sample(30000, random_state=42)
        ax.scatter(obs_bg['x_plotting'], obs_bg['y_plotting'],
                   c='lightgray', s=0.3, rasterized=True, alpha=0.3)

    # Plot highlighted annotation in red
    obs_fg = obs[highlight]
    if len(obs_fg) > 0:
        ax.scatter(obs_fg['x_plotting'], obs_fg['y_plotting'],
                   c='red', s=0.5, rasterized=True, alpha=0.7)

    title = f'{annotation}\n({organ})\np={row["p_cauchy"]:.2e}'
    ax.set_title(title, fontsize=8)
    ax.set_aspect('equal')
    ax.axis('off')

    if idx == 0:
        ax.text(-0.05, 1.18, 'C', fontsize=20, fontweight='bold',
                va='top', transform=ax.transAxes)

for fmt in ['png', 'pdf']:
    fig.savefig(f'{OUTDIR}/FigureS4.{fmt}', dpi=300, bbox_inches='tight')
    print(f"Saved FigureS4.{fmt}")

import os
for fmt in ['png', 'pdf']:
    sz = os.path.getsize(f'{OUTDIR}/FigureS4.{fmt}')
    print(f"FigureS4.{fmt}: {sz/1024/1024:.1f} MB")
plt.close()
print("Done!")
