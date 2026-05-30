#!/usr/bin/env python
"""Figure S5: AD risk gene extended analysis"""
import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
import scipy.cluster.hierarchy as sch
import scanpy as sc

# Style
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 11

BASE = ''
OUT = f'{BASE}figures/paper_v2'
os.makedirs(OUT, exist_ok=True)

# ---- Panel A: Gene PCC heatmap ----
organs_dirs = [
    'Bone_Marrow_CTRL1', 'Brain_CTRL1', 'Brown_Fat_CTRL1', 'Colon_CTRL1',
    'Heart_CTRL1', 'Kidney_CTRL1', 'Liver_CTRL1', 'Lung_CTRL1',
    'Lymph_Node_CTRL1', 'Muscle_CTRL1', 'Pancreas_CTRL1', 'Skin_CTRL1',
    'Small_Intestine_CTRL1', 'Spleen_CTRL1', 'Stomach_CTRL1', 'Thymus_CTRL1'
]

all_gene_data = []
for od in organs_dirs:
    organ_name = od.replace('_CTRL1', '')
    csv_path = f'{BASE}models/gsmap_output/{od}/report/AD/{od}_AD_Gene_Diagnostic_Info.csv'
    if not os.path.exists(csv_path):
        print(f"WARNING: missing {csv_path}")
        continue
    df = pd.read_csv(csv_path)
    df['organ'] = organ_name
    all_gene_data.append(df)

gene_df = pd.concat(all_gene_data, ignore_index=True)

# Max PCC per gene per organ
pcc_pivot = gene_df.groupby(['Gene', 'organ'])['PCC'].max().reset_index()
pcc_matrix = pcc_pivot.pivot(index='Gene', columns='organ', values='PCC').fillna(0)

# Top 50 by mean PCC
mean_pcc = pcc_matrix.mean(axis=1).sort_values(ascending=False)
top50 = mean_pcc.head(50).index
heatmap_data = pcc_matrix.loc[top50]

# Cluster rows
if len(heatmap_data) > 1:
    linkage = sch.linkage(heatmap_data.values, method='ward', metric='euclidean')
    dendro = sch.dendrogram(linkage, no_plot=True)
    row_order = dendro['leaves']
    heatmap_data = heatmap_data.iloc[row_order]

# Clean organ names for display
organ_display = [c.replace('_', ' ') for c in heatmap_data.columns]

# ---- Panel B: Top 10 genes spatial expression ----
top10_genes = mean_pcc.head(10).index.tolist()
print(f"Top 10 genes: {top10_genes}")

# Load Brain h5ad
adata = sc.read_h5ad(f'{BASE}data/st/per_organ/Brain_CTRL1.h5ad')
x = adata.obs['x_plotting'].values
y = adata.obs['y_plotting'].values

# Build case-insensitive lookup: uppercase -> actual var_name
var_upper_map = {v.upper(): v for v in adata.var_names}

# Get expression for each gene
gene_exprs = {}
for g in top10_genes:
    actual_name = var_upper_map.get(g.upper())
    if actual_name is not None:
        expr = np.array(adata[:, actual_name].X.todense()).flatten() if hasattr(adata.X, 'todense') else adata[:, actual_name].X.flatten()
        gene_exprs[g] = np.log1p(expr)
    else:
        print(f"Gene {g} not found in Brain h5ad")
        gene_exprs[g] = np.zeros(len(x))

# ---- PLOT ----
fig = plt.figure(figsize=(16, 20))

# Panel A
ax_heat = fig.add_axes([0.08, 0.52, 0.85, 0.44])
im = ax_heat.imshow(heatmap_data.values, aspect='auto', cmap='YlOrRd', interpolation='nearest')
ax_heat.set_yticks(range(len(heatmap_data)))
ax_heat.set_yticklabels(heatmap_data.index, fontsize=8)
ax_heat.set_xticks(range(len(organ_display)))
ax_heat.set_xticklabels(organ_display, rotation=45, ha='right', fontsize=9)
ax_heat.set_title('Gene PCC across organs', fontsize=13, fontweight='bold', pad=10)
cb = fig.colorbar(im, ax=ax_heat, fraction=0.02, pad=0.02)
cb.set_label('PCC', fontsize=10)
ax_heat.text(-0.04, 1.02, 'A', transform=ax_heat.transAxes, fontsize=16, fontweight='bold', va='bottom')

# Panel B: 2x5 grid
genes_to_plot = [g for g in top10_genes if g in gene_exprs][:10]
for i, g in enumerate(genes_to_plot):
    row, col = divmod(i, 5)
    ax = fig.add_axes([0.05 + col * 0.19, 0.26 - row * 0.22, 0.17, 0.20])
    sc_plot = ax.scatter(x, y, c=gene_exprs[g], s=0.3, cmap='viridis', rasterized=True, vmin=0)
    ax.set_title(g, fontsize=10, fontweight='bold')
    ax.set_aspect('equal')
    ax.axis('off')
    if i == 0:
        ax.text(-0.1, 1.15, 'B', transform=ax.transAxes, fontsize=16, fontweight='bold', va='bottom')

# Colorbar for Panel B
cax = fig.add_axes([0.92, 0.06, 0.015, 0.16])
cb2 = fig.colorbar(sc_plot, cax=cax)
cb2.set_label('log1p(expression)', fontsize=9)

plt.savefig(f'{OUT}/FigureS5.png', dpi=200, bbox_inches='tight')
plt.savefig(f'{OUT}/FigureS5.pdf', dpi=200, bbox_inches='tight')
plt.close()
print(f"FigureS5 saved")
# Report sizes
for ext in ['png', 'pdf']:
    sz = os.path.getsize(f'{OUT}/FigureS5.{ext}')
    print(f"  FigureS5.{ext}: {sz/1024:.0f} KB")
