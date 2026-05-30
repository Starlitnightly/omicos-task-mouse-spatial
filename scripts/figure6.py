#!/usr/bin/env python
"""Figure 6: AD risk gene landscape for Cell paper."""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import numpy as np
import pandas as pd
import scanpy as sc
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Set Arial font
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 11

BASE = Path('')
OUT = BASE / 'figures' / 'paper_v2'

# ── Data loading ──
ORGAN_ORDER = ['Brain', 'Bone_Marrow', 'Thymus', 'Spleen', 'Lymph_Node',
               'Lung', 'Heart', 'Liver', 'Kidney', 'Pancreas',
               'Stomach', 'Small_Intestine', 'Colon', 'Muscle', 'Skin', 'Brown_Fat']

AD_RISK_GENES = ['APOE', 'BIN1', 'CLU', 'TREM2', 'PICALM', 'CR1', 'CD33',
                 'ABCA7', 'SORL1', 'APP', 'PSEN1', 'ADAM10', 'CD2AP', 'INPP5D',
                 'SPI1', 'PLCG2', 'ABI3', 'MEF2C', 'PTK2B']

def load_all_organs():
    """Load Gene_Diagnostic_Info for all CTRL1 organs."""
    all_data = {}
    for organ in ORGAN_ORDER:
        csv_path = (BASE / 'models' / 'gsmap_output' / f'{organ}_CTRL1' /
                    'report' / 'AD' / f'{organ}_CTRL1_AD_Gene_Diagnostic_Info.csv')
        if csv_path.exists():
            df = pd.read_csv(csv_path)
            all_data[organ] = df
    return all_data

all_organs = load_all_organs()

# ══════════════════════════════════════════════
# Figure 6
# ══════════════════════════════════════════════
fig = plt.figure(figsize=(16, 16))

# Layout: top row has 2 panels, bottom row has 1 wide panel
gs = fig.add_gridspec(2, 2, height_ratios=[1, 0.8], hspace=0.35, wspace=0.35,
                      left=0.08, right=0.95, top=0.95, bottom=0.05)

# ── Panel A: Risk gene expression heatmap across organs ──
ax_a = fig.add_subplot(gs[0, 0])

# Build max PCC matrix: gene x organ
pcc_matrix = pd.DataFrame(index=AD_RISK_GENES, columns=ORGAN_ORDER, dtype=float)
for organ, df in all_organs.items():
    for gene in AD_RISK_GENES:
        gene_data = df[df['Gene'] == gene]
        if len(gene_data) > 0:
            pcc_matrix.loc[gene, organ] = gene_data['PCC'].max()

# Filter to organs present
organs_present = [o for o in ORGAN_ORDER if o in all_organs]
pcc_sub = pcc_matrix[organs_present].astype(float)

# Pretty organ names for display
organ_labels = [o.replace('_', ' ') for o in organs_present]

im = ax_a.imshow(pcc_sub.values, aspect='auto', cmap='YlOrRd',
                 interpolation='nearest')
# Gray for NaN
ax_a.set_facecolor('#d0d0d0')
im.set_clim(vmin=np.nanmin(pcc_sub.values), vmax=np.nanmax(pcc_sub.values))

# Annotate
for i in range(pcc_sub.shape[0]):
    for j in range(pcc_sub.shape[1]):
        val = pcc_sub.values[i, j]
        if not np.isnan(val):
            color = 'white' if val > 0.5 * np.nanmax(pcc_sub.values) else 'black'
            ax_a.text(j, i, f'{val:.2f}', ha='center', va='center',
                     fontsize=5.5, color=color)

ax_a.set_xticks(range(len(organs_present)))
ax_a.set_xticklabels(organ_labels, rotation=45, ha='right', fontsize=8)
ax_a.set_yticks(range(len(AD_RISK_GENES)))
ax_a.set_yticklabels(AD_RISK_GENES, fontsize=8)
ax_a.set_title('AD risk genes: PCC with AD enrichment across organs', fontsize=11, fontweight='bold')
plt.colorbar(im, ax=ax_a, shrink=0.6, label='PCC')
ax_a.text(-0.15, 1.05, 'A', transform=ax_a.transAxes, fontsize=20, fontweight='bold', va='top')

# ── Panel B: Top 20 genes by PCC (brain) ──
ax_b = fig.add_subplot(gs[0, 1])

brain_df = all_organs['Brain']
# For each gene, get the row with max PCC
brain_top = brain_df.loc[brain_df.groupby('Gene')['PCC'].idxmax()]
brain_top = brain_top.nlargest(20, 'PCC')

# Color by annotation
annotations = brain_top['Annotation'].unique()
cmap_annot = plt.cm.Set2(np.linspace(0, 1, len(annotations)))
annot_colors = {a: cmap_annot[i] for i, a in enumerate(annotations)}

colors = [annot_colors[a] for a in brain_top['Annotation']]
brain_top_sorted = brain_top.sort_values('PCC', ascending=True)
colors_sorted = [annot_colors[a] for a in brain_top_sorted['Annotation']]

bars = ax_b.barh(range(len(brain_top_sorted)), brain_top_sorted['PCC'].values,
                 color=colors_sorted, edgecolor='gray', linewidth=0.5)
ax_b.set_yticks(range(len(brain_top_sorted)))
ax_b.set_yticklabels(brain_top_sorted['Gene'].values, fontsize=8)
ax_b.set_xlabel('PCC', fontsize=10)
ax_b.set_title('Brain: top 20 AD-correlated genes', fontsize=11, fontweight='bold')

# Legend for annotations
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor=annot_colors[a], edgecolor='gray', label=a)
                   for a in annotations]
ax_b.legend(handles=legend_elements, loc='lower right', fontsize=6,
            title='Annotation', title_fontsize=7, framealpha=0.9)
ax_b.text(-0.15, 1.05, 'B', transform=ax_b.transAxes, fontsize=20, fontweight='bold', va='top')

# ── Panel C: Spatial expression of AD risk genes in brain ──
gs_bottom = gs[1, :].subgridspec(1, 6, wspace=0.3)

# Load brain h5ad
print("Loading Brain h5ad...")
adata = sc.read_h5ad(BASE / 'data' / 'st' / 'per_organ' / 'Brain_CTRL1.h5ad')

genes_mouse = ['Apoe', 'Bin1', 'Picalm', 'Trem2', 'App', 'Sorl1']
genes_human = ['APOE', 'BIN1', 'PICALM', 'TREM2', 'APP', 'SORL1']

x = adata.obs['x_plotting'].values.astype(float)
y = adata.obs['y_plotting'].values.astype(float)

for idx, (gm, gh) in enumerate(zip(genes_mouse, genes_human)):
    ax = fig.add_subplot(gs_bottom[idx])

    if gm in adata.var_names:
        gene_idx = list(adata.var_names).index(gm)
        expr_raw = adata.X[:, gene_idx]
        import scipy.sparse
        if scipy.sparse.issparse(expr_raw):
            expr = np.asarray(expr_raw.todense()).flatten()
        else:
            expr = np.array(expr_raw).flatten()
        expr_log = np.log1p(expr)
    else:
        expr_log = np.zeros(len(x))

    sc_plot = ax.scatter(x, y, c=expr_log, cmap='viridis', s=0.3,
                         rasterized=True, edgecolors='none')
    ax.set_title(gh, fontsize=10, fontweight='bold')
    ax.set_aspect('equal')
    ax.axis('off')
    plt.colorbar(sc_plot, ax=ax, shrink=0.6, pad=0.02,
                 label='log1p(expr)' if idx == 5 else '')

# Panel C label
fig.text(0.03, 0.38, 'C', fontsize=20, fontweight='bold', va='top')

print("Saving Figure 6...")
fig.savefig(OUT / 'Figure6.png', dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(OUT / 'Figure6.pdf', bbox_inches='tight', facecolor='white')
plt.close()
print("Figure 6 done.")
