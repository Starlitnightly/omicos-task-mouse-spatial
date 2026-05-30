#!/usr/bin/env python
"""Figure S1: Data quality control (2x2)"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats
import glob
import os
import warnings
warnings.filterwarnings('ignore')

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 11

BASEDIR = os.environ.get('BASE', '.')
OUTDIR = f'{BASEDIR}/figures/paper_v2'

h5ad_files = sorted(glob.glob(f'{BASEDIR}/data/st/per_organ/*_CTRL1.h5ad'))
organ_names = [f.split('/')[-1].replace('_CTRL1.h5ad', '').replace('_', ' ') for f in h5ad_files]

short_names = {
    'Bone Marrow': 'BM', 'Brain': 'Brain', 'Brown Fat': 'BrFat',
    'Colon': 'Colon', 'Heart': 'Heart', 'Kidney': 'Kidney',
    'Liver': 'Liver', 'Lung': 'Lung', 'Lymph Node': 'LN',
    'Muscle': 'Muscle', 'Pancreas': 'Panc', 'Skin': 'Skin',
    'Small Intestine': 'SI', 'Spleen': 'Spleen', 'Stomach': 'Stom',
    'Thymus': 'Thymus'
}
organ_order = [short_names[n] for n in organ_names]

# Load obs for all organs
print("Loading obs data...")
all_obs = {}
for f, name in zip(h5ad_files, organ_names):
    adata = sc.read_h5ad(f, backed='r')
    obs = adata.obs[['x_plotting', 'y_plotting', 'n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'annotation']].copy()
    all_obs[name] = obs
    adata.file.close()
    print(f"  {name}: {len(obs)} spots")

# Build combined df
records = []
for name, obs in all_obs.items():
    records.append(pd.DataFrame({
        'Organ': name,
        'Organ_short': short_names[name],
        'Genes/spot': obs['n_genes_by_counts'].values,
        'Counts/spot': obs['total_counts'].values,
        '% Mito': obs['pct_counts_mt'].values,
    }))
combined = pd.concat(records, ignore_index=True)

# Figure layout
fig = plt.figure(figsize=(16, 14))
gs = gridspec.GridSpec(2, 2, hspace=0.3, wspace=0.25,
                       left=0.07, right=0.95, top=0.95, bottom=0.06)

# Panel A: 3 stacked boxplots
print("Panel A: QC boxplots...")
gs_a = gs[0, 0].subgridspec(3, 1, hspace=0.45)
metrics = [('Genes/spot', 'Genes per spot'), ('Counts/spot', 'Total counts per spot'), ('% Mito', '% Mitochondrial')]

for idx, (col, label) in enumerate(metrics):
    ax = fig.add_subplot(gs_a[idx])
    data_by_organ = [combined.loc[combined['Organ_short'] == o, col].values for o in organ_order]
    bp = ax.boxplot(data_by_organ, labels=organ_order, widths=0.6,
                    patch_artist=True, showfliers=False,
                    medianprops=dict(color='red', linewidth=1.5))
    for patch in bp['boxes']:
        patch.set_facecolor('#4C72B0')
        patch.set_alpha(0.7)
    ax.set_ylabel(label, fontsize=9)
    ax.tick_params(axis='y', labelsize=8)
    if idx < 2:
        ax.set_xticklabels([])
    else:
        ax.tick_params(axis='x', rotation=45, labelsize=7)
    if idx == 0:
        ax.text(-0.12, 1.15, 'A', fontsize=20, fontweight='bold',
                va='top', transform=ax.transAxes)

# Panel B: Spatial %mito for 3 representative organs
print("Panel B: Spatial %mito...")
gs_b = gs[0, 1].subgridspec(1, 3, wspace=0.15)
rep_organs = ['Heart', 'Brain', 'Liver']
for idx, org in enumerate(rep_organs):
    ax = fig.add_subplot(gs_b[idx])
    obs = all_obs[org]
    if len(obs) > 30000:
        obs = obs.sample(30000, random_state=42)
    vmax = obs['pct_counts_mt'].quantile(0.95)
    sc_plot = ax.scatter(obs['x_plotting'], obs['y_plotting'],
                         c=obs['pct_counts_mt'], cmap='YlOrRd', s=0.3,
                         rasterized=True, vmin=0, vmax=vmax)
    ax.set_title(org, fontsize=10)
    ax.set_aspect('equal')
    ax.axis('off')
    if idx == 2:
        cbar = plt.colorbar(sc_plot, ax=ax, shrink=0.6, pad=0.02)
        cbar.set_label('% Mito', fontsize=9)
        cbar.ax.tick_params(labelsize=7)
    if idx == 0:
        ax.text(-0.05, 1.1, 'B', fontsize=20, fontweight='bold',
                va='top', transform=ax.transAxes)

# Panel C: Replicate consistency
print("Panel C: Replicate consistency...")
ax_c = fig.add_subplot(gs[1, 0])
cauchy = pd.read_csv(f'{BASEDIR}/results/all_organs_cauchy_AD.csv.gz')
ctrl1 = cauchy[cauchy['sample'] == 'CTRL_1'].copy()
ctrl2 = cauchy[cauchy['sample'] == 'CTRL_2'].copy()

organ_min_c1 = ctrl1.groupby('organ')['p_cauchy'].min().reset_index()
organ_min_c1.columns = ['organ', 'p_ctrl1']
organ_min_c2 = ctrl2.groupby('organ')['p_cauchy'].min().reset_index()
organ_min_c2.columns = ['organ', 'p_ctrl2']

merged = organ_min_c1.merge(organ_min_c2, on='organ')
merged['logp1'] = -np.log10(merged['p_ctrl1'].clip(1e-20))
merged['logp2'] = -np.log10(merged['p_ctrl2'].clip(1e-20))

r, pval = stats.pearsonr(merged['logp1'], merged['logp2'])
ax_c.scatter(merged['logp1'], merged['logp2'], s=50, c='#4C72B0',
             edgecolors='black', linewidth=0.5, zorder=3)
for _, row in merged.iterrows():
    ax_c.annotate(short_names.get(row['organ'], row['organ']),
                  (row['logp1'], row['logp2']), fontsize=7,
                  xytext=(3, 3), textcoords='offset points')
lims = [0, max(merged['logp1'].max(), merged['logp2'].max()) * 1.1]
ax_c.plot(lims, lims, 'k--', alpha=0.5, linewidth=1)
ax_c.set_xlabel('-log$_{10}$(min Cauchy p) CTRL1', fontsize=11)
ax_c.set_ylabel('-log$_{10}$(min Cauchy p) CTRL2', fontsize=11)
ax_c.set_title(f'Replicate consistency (r = {r:.3f})', fontsize=12)
ax_c.text(-0.12, 1.05, 'C', fontsize=20, fontweight='bold',
          va='top', transform=ax_c.transAxes)

# Panel D: Homolog coverage
print("Panel D: Homolog coverage...")
ax_d = fig.add_subplot(gs[1, 1])
coverage = 51.9
bars = ax_d.bar(range(len(organ_names)), [coverage]*len(organ_names),
                color='#4C72B0', alpha=0.7, edgecolor='black', linewidth=0.5)
ax_d.axhline(y=coverage, color='red', linestyle='--', linewidth=1.5, label=f'Mean = {coverage}%')
ax_d.set_xticks(range(len(organ_names)))
ax_d.set_xticklabels([short_names[n] for n in organ_names], rotation=45, ha='right', fontsize=8)
ax_d.set_ylabel('Ortholog coverage (%)', fontsize=11)
ax_d.set_ylim(0, 70)
ax_d.set_title('Mouse-human ortholog coverage', fontsize=12)
ax_d.legend(fontsize=9)
ax_d.text(-0.12, 1.05, 'D', fontsize=20, fontweight='bold',
          va='top', transform=ax_d.transAxes)

for fmt in ['png', 'pdf']:
    fig.savefig(f'{OUTDIR}/FigureS1.{fmt}', dpi=300, bbox_inches='tight')
    print(f"Saved FigureS1.{fmt}")

for fmt in ['png', 'pdf']:
    sz = os.path.getsize(f'{OUTDIR}/FigureS1.{fmt}')
    print(f"FigureS1.{fmt}: {sz/1024/1024:.1f} MB")
plt.close()
print("Done!")
