#!/usr/bin/env python
"""Figure S2: gsMap technical validation (2x2)"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats
import umap
import warnings
warnings.filterwarnings('ignore')

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 11

BASEDIR = '.'
OUTDIR = f'{BASEDIR}/figures/paper_v2'

fig = plt.figure(figsize=(16, 14))
gs = gridspec.GridSpec(2, 2, hspace=0.35, wspace=0.3,
                       left=0.08, right=0.95, top=0.95, bottom=0.08)

# Top row: split into left (Brain UMAP) and right (Bone Marrow UMAP) for Panel A
gs_top = gs[0, :].subgridspec(1, 2, wspace=0.3)

# Panel A: GVAE UMAP for Brain and Bone Marrow
print("Panel A: GVAE UMAP...")
organs_umap = ['Brain', 'Bone_Marrow']
organ_labels = ['Brain', 'Bone Marrow']

for idx, (org, label) in enumerate(zip(organs_umap, organ_labels)):
    ax = fig.add_subplot(gs_top[idx])
    latent_f = f'{BASEDIR}/models/gsmap_output/{org}_CTRL1/find_latent_representations/{org}_CTRL1_add_latent.h5ad'
    adata = sc.read_h5ad(latent_f, backed='r')
    latent = adata.obsm['latent_GVAE'][:].copy()
    annotations = adata.obs['annotation'].values.copy()
    adata.file.close()

    n = len(latent)
    if n > 20000:
        idx_sub = np.random.RandomState(42).choice(n, 20000, replace=False)
    else:
        idx_sub = np.arange(n)

    latent_sub = latent[idx_sub]
    ann_sub = annotations[idx_sub]

    print(f"  Computing UMAP for {label}...")
    reducer = umap.UMAP(n_components=2, random_state=42, n_neighbors=30, min_dist=0.3)
    embedding = reducer.fit_transform(latent_sub)

    unique_ann = sorted(set(ann_sub))
    cmap = plt.cm.get_cmap('tab20', len(unique_ann))
    colors = {a: cmap(i) for i, a in enumerate(unique_ann)}

    for a in unique_ann:
        mask = ann_sub == a
        ax.scatter(embedding[mask, 0], embedding[mask, 1], c=[colors[a]],
                   s=0.5, alpha=0.5, label=a, rasterized=True)
    ax.set_title(f'{label} GVAE latent UMAP', fontsize=12)
    ax.set_xlabel('UMAP1', fontsize=10)
    ax.set_ylabel('UMAP2', fontsize=10)
    ax.legend(fontsize=5, ncol=2, markerscale=5, loc='best',
              framealpha=0.8, handletextpad=0.1, columnspacing=0.5)
    if idx == 0:
        ax.text(-0.08, 1.05, 'A', fontsize=20, fontweight='bold',
                va='top', transform=ax.transAxes)

# Load spot data (needed for panels B, C, D)
print("Loading spot data...")
spot_data = pd.read_csv(f'{BASEDIR}/results/all_organs_spot_AD_pvalues.csv.gz')
ctrl1_spots = spot_data[spot_data['sample'] == 'CTRL_1'].copy()

short_names = {
    'Bone Marrow': 'BM', 'Brain': 'Brain', 'Brown Fat': 'BrFat',
    'Colon': 'Colon', 'Heart': 'Heart', 'Kidney': 'Kidney',
    'Liver': 'Liver', 'Lung': 'Lung', 'Lymph Node': 'LN',
    'Muscle': 'Muscle', 'Pancreas': 'Panc', 'Skin': 'Skin',
    'Small Intestine': 'SI', 'Spleen': 'Spleen', 'Stomach': 'Stom',
    'Thymus': 'Thymus'
}

# Panel B: GSS distribution
print("Panel B: GSS distribution...")
ax_b = fig.add_subplot(gs[1, 0])
import pyarrow.feather as feather

organs_gss = ['Brain', 'Bone_Marrow', 'Heart', 'Liver']
colors_gss = ['#4C72B0', '#DD8452', '#55A868', '#C44E52']

for org, color in zip(organs_gss, colors_gss):
    gss_f = f'{BASEDIR}/models/gsmap_output/{org}_CTRL1/latent_to_gene/{org}_CTRL1_gene_marker_score.feather'
    df = feather.read_feather(gss_f)
    spot_cols = df.columns[1:]
    sample_cols = np.random.RandomState(42).choice(spot_cols, min(1000, len(spot_cols)), replace=False)
    mean_gss = df[sample_cols].mean(axis=1).values
    mean_gss = mean_gss[mean_gss > 0]
    ax_b.hist(mean_gss, bins=50, alpha=0.5, color=color, label=org.replace('_', ' '),
              density=True, histtype='stepfilled', edgecolor=color, linewidth=0.5)

ax_b.set_xlabel('Mean Gene Specificity Score', fontsize=11)
ax_b.set_ylabel('Density', fontsize=11)
ax_b.set_title('Gene Specificity Score distribution', fontsize=12)
ax_b.legend(fontsize=9)
ax_b.set_xlim(0, ax_b.get_xlim()[1])
ax_b.text(-0.08, 1.05, 'B', fontsize=20, fontweight='bold',
          va='top', transform=ax_b.transAxes)

# Panel C: QQ plot
print("Panel C: QQ plot...")
ax_c = fig.add_subplot(gs[1, 1])

organs_list = sorted(ctrl1_spots['organ'].unique())
cmap_organs = plt.cm.get_cmap('tab20', len(organs_list))

for i, org in enumerate(organs_list):
    pvals = ctrl1_spots.loc[ctrl1_spots['organ'] == org, 'p'].values
    pvals = pvals[np.isfinite(pvals) & (pvals > 0)]
    pvals_sorted = np.sort(pvals)
    n = len(pvals_sorted)
    expected = -np.log10(np.arange(1, n+1) / (n+1))
    observed = -np.log10(pvals_sorted)[::-1]
    if n > 2000:
        idx_sub = np.linspace(0, n-1, 2000, dtype=int)
        expected = expected[idx_sub]
        observed = observed[idx_sub]
    ax_c.scatter(expected, observed, s=1, alpha=0.6,
                 color=cmap_organs(i), label=short_names.get(org, org), rasterized=True)

maxval = max(ax_c.get_xlim()[1], ax_c.get_ylim()[1])
ax_c.plot([0, maxval], [0, maxval], 'k--', linewidth=1, alpha=0.5)
ax_c.set_xlabel('Expected -log$_{10}$(p)', fontsize=11)
ax_c.set_ylabel('Observed -log$_{10}$(p)', fontsize=11)
ax_c.set_title('QQ plot by organ', fontsize=12)
ax_c.legend(fontsize=6, ncol=2, markerscale=5, loc='upper left')
ax_c.text(-0.08, 1.05, 'C', fontsize=20, fontweight='bold',
          va='top', transform=ax_c.transAxes)

# Now add Panel D as an inset or extra panel - size bias check
# We need to fit it in. Let's overlay it as text + inset in Panel C area
# Actually, let's restructure: make it a 2x3 grid with Panel D in position [1,2]
# But we already created the figure. Let's add Panel D as an inset axes.
# Better approach: add a small inset inside Panel C's area
print("Panel D: Size bias check (inset)...")

# Compute per-organ stats
organ_stats = ctrl1_spots.groupby('organ').agg(
    n_spots=('spot', 'count'),
    median_logp=('logp', 'median')
).reset_index()

# Add as inset in the bottom-right area
ax_d = fig.add_axes([0.62, 0.12, 0.30, 0.22])
ax_d.scatter(organ_stats['n_spots'], organ_stats['median_logp'],
             s=50, c='#4C72B0', edgecolors='black', linewidth=0.5, zorder=3)
for _, row in organ_stats.iterrows():
    ax_d.annotate(short_names.get(row['organ'], row['organ']),
                  (row['n_spots'], row['median_logp']),
                  fontsize=6, xytext=(3, 3), textcoords='offset points')

r_size, p_size = stats.pearsonr(organ_stats['n_spots'], organ_stats['median_logp'])
ax_d.set_xlabel('Number of spots', fontsize=9)
ax_d.set_ylabel('Median -log$_{10}$(p)', fontsize=9)
ax_d.set_title(f'Size bias (r = {r_size:.3f}, p = {p_size:.3f})', fontsize=9)
ax_d.tick_params(labelsize=7)
ax_d.text(-0.08, 1.12, 'D', fontsize=20, fontweight='bold',
          va='top', transform=ax_d.transAxes)

for fmt in ['png', 'pdf']:
    fig.savefig(f'{OUTDIR}/FigureS2.{fmt}', dpi=300, bbox_inches='tight')
    print(f"Saved FigureS2.{fmt}")

import os
for fmt in ['png', 'pdf']:
    sz = os.path.getsize(f'{OUTDIR}/FigureS2.{fmt}')
    print(f"FigureS2.{fmt}: {sz/1024/1024:.1f} MB")
plt.close()
print("Done!")
