#!/usr/bin/env python
"""
Figure 3: Two distinct AD risk pathways — myeloid and gut epithelial APP processing
2×2 grid for Cell paper
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.colors import Normalize
import glob
import os
import warnings
warnings.filterwarnings('ignore')

# ---------- Config ----------
BASE = os.environ.get('BASE', '.')
OUT_DIR = os.path.join(BASE, 'figures/paper_v2')
os.makedirs(OUT_DIR, exist_ok=True)

# Font setup
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 11
plt.rcParams['axes.titlesize'] = 13
plt.rcParams['axes.titleweight'] = 'bold'

MYELOID_GENES = ['TREM2', 'CD33', 'SPI1', 'INPP5D', 'PLCG2', 'ABI3', 'MS4A6A']
APP_GENES = ['APP', 'SORL1', 'PSEN1', 'ADAM10', 'CD2AP', 'PICALM', 'CR1']

# ---------- Load data ----------
cauchy = pd.read_csv(os.path.join(BASE, 'results/all_organs_cauchy_AD.csv.gz'))
cauchy_ctrl1 = cauchy[cauchy['sample'] == 'CTRL_1'].copy()

# Load Gene Diagnostic Info for all organs
gdi_frames = []
for f in sorted(glob.glob(os.path.join(BASE, 'models/gsmap_output/*_CTRL1/report/AD/*_Gene_Diagnostic_Info.csv'))):
    organ = f.split('gsmap_output/')[1].split('_CTRL1')[0].replace('_', ' ')
    df = pd.read_csv(f)
    df['organ'] = organ
    gdi_frames.append(df)
gdi = pd.concat(gdi_frames, ignore_index=True)

# ---------- Figure ----------
fig, axes = plt.subplots(2, 2, figsize=(16, 14))
fig.subplots_adjust(hspace=0.35, wspace=0.35, top=0.94, bottom=0.08, left=0.10, right=0.96)

# ======== PANEL A: Gut immune vs epithelial enrichment ========
ax = axes[0, 0]
gut_organs = ['Colon', 'Small Intestine', 'Stomach']
gut = cauchy_ctrl1[cauchy_ctrl1['organ'].isin(gut_organs)].copy()

# Classify annotations
immune_keywords = ['Myeloid', 'Immune', 'Lymph', 'macrophage', 'Dendri', 'Other']
epithelial_keywords = ['Tissue', 'Epithelium', 'Goblet', 'Secretory', 'Enterocyte',
                       'Foveolar', 'Gland', 'Colon', 'Small Intestine', 'Chief', 'Parietal']

def classify_annot(name):
    name_lower = name.lower()
    # 'Other' is ambiguous but in spatial data often contains immune/stromal
    for kw in ['myeloid', 'immune', 'lymph', 'macrophage', 'dendri']:
        if kw in name_lower:
            return 'Immune/Myeloid'
    if 'other' in name_lower:
        return 'Immune/Myeloid'
    return 'Epithelial/Tissue'

gut['category'] = gut['annotation'].apply(classify_annot)
gut['neg_log10_p'] = -np.log10(gut['p_cauchy'].clip(lower=1e-300))

# For each gut organ, get min p_cauchy (max -log10p) per category
grouped = gut.groupby(['organ', 'category'])['neg_log10_p'].max().reset_index()

x_pos = np.arange(len(gut_organs))
width = 0.35
colors = {'Immune/Myeloid': '#D32F2F', 'Epithelial/Tissue': '#388E3C'}

for i, cat in enumerate(['Immune/Myeloid', 'Epithelial/Tissue']):
    vals = []
    for org in gut_organs:
        sub = grouped[(grouped['organ'] == org) & (grouped['category'] == cat)]
        vals.append(sub['neg_log10_p'].values[0] if len(sub) > 0 else 0)
    bars = ax.bar(x_pos + (i - 0.5) * width, vals, width, label=cat,
                  color=colors[cat], edgecolor='white', linewidth=0.5)

ax.axhline(-np.log10(0.05), color='grey', ls='--', lw=1, label='p = 0.05')
ax.axhline(-np.log10(0.05 / 86), color='grey', ls=':', lw=1, label='Bonferroni')
ax.set_xticks(x_pos)
ax.set_xticklabels(gut_organs, fontsize=11)
ax.set_ylabel('$-\\log_{10}$(Cauchy $p$)', fontsize=11)
ax.set_title('Gut: immune vs epithelial AD enrichment', fontsize=13, fontweight='bold')
ax.legend(fontsize=9, loc='upper right')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# ======== PANEL B: Two-pathway heatmap ========
ax = axes[0, 1]

all_organs = sorted(gdi['organ'].unique())
all_pathway_genes = MYELOID_GENES + APP_GENES

# For each organ, get max PCC per gene across annotations
heatmap_data = []
for gene in all_pathway_genes:
    row = []
    for organ in all_organs:
        sub = gdi[(gdi['organ'] == organ) & (gdi['Gene'] == gene)]
        if len(sub) > 0:
            row.append(sub['PCC'].max())
        else:
            row.append(np.nan)
    heatmap_data.append(row)

heatmap_arr = np.array(heatmap_data)

im = ax.imshow(heatmap_arr, aspect='auto', cmap='YlOrRd', interpolation='nearest')
ax.set_xticks(range(len(all_organs)))
ax.set_xticklabels([o.replace(' ', '\n') for o in all_organs], fontsize=8, rotation=45, ha='right')
ax.set_yticks(range(len(all_pathway_genes)))
ax.set_yticklabels(all_pathway_genes, fontsize=9)

# Horizontal line separating pathways
ax.axhline(len(MYELOID_GENES) - 0.5, color='black', lw=2)

# Add colored brackets on the left side using y-axis tick label colors
for idx in range(len(MYELOID_GENES)):
    ax.get_yticklabels()[idx].set_color('#D32F2F')
for idx in range(len(MYELOID_GENES), len(all_pathway_genes)):
    ax.get_yticklabels()[idx].set_color('#388E3C')

cbar = fig.colorbar(im, ax=ax, shrink=0.6, pad=0.02)
cbar.set_label('Max PCC', fontsize=10)
ax.set_title('Two AD risk gene pathways', fontsize=13, fontweight='bold')

# ======== PANEL C: Gut APP pathway by annotation ========
ax = axes[1, 0]

# Use all gut organs' annotations for richer view
gut_organs_underscore = ['Colon', 'Small Intestine', 'Stomach']
gut_gdi = gdi[gdi['organ'].isin(gut_organs_underscore)].copy()
gut_gdi['organ_annot'] = gut_gdi['organ'] + ': ' + gut_gdi['Annotation']

gut_annots = sorted(gut_gdi['organ_annot'].unique())

c_data = []
for gene in APP_GENES:
    row = []
    for annot in gut_annots:
        sub = gut_gdi[(gut_gdi['Gene'] == gene) & (gut_gdi['organ_annot'] == annot)]
        if len(sub) > 0:
            row.append(sub['PCC'].max())
        else:
            row.append(np.nan)
    c_data.append(row)

c_arr = np.array(c_data)

im2 = ax.imshow(c_arr, aspect='auto', cmap='YlOrRd', interpolation='nearest')
ax.set_xticks(range(len(gut_annots)))
ax.set_xticklabels([a.replace(': ', ':\n') for a in gut_annots], fontsize=8, rotation=45, ha='right')
ax.set_yticks(range(len(APP_GENES)))
ax.set_yticklabels(APP_GENES, fontsize=10)

# Annotate PCC values
for i in range(len(APP_GENES)):
    for j in range(len(gut_annots)):
        val = c_arr[i, j]
        if not np.isnan(val):
            color = 'white' if val > 0.3 else 'black'
            ax.text(j, i, f'{val:.2f}', ha='center', va='center', fontsize=7, color=color)

cbar2 = fig.colorbar(im2, ax=ax, shrink=0.6, pad=0.02)
cbar2.set_label('PCC', fontsize=10)
ax.set_title('Gut organs: APP pathway genes by region', fontsize=13, fontweight='bold')

# ======== PANEL D: Myeloid vs APP pathway per organ ========
ax = axes[1, 1]

organ_myeloid_mean = []
organ_app_mean = []
organ_labels = []

for organ in all_organs:
    organ_data = gdi[gdi['organ'] == organ]
    # Mean of max PCC per gene
    myeloid_vals = []
    for g in MYELOID_GENES:
        sub = organ_data[organ_data['Gene'] == g]
        if len(sub) > 0:
            myeloid_vals.append(sub['PCC'].max())
    app_vals = []
    for g in APP_GENES:
        sub = organ_data[organ_data['Gene'] == g]
        if len(sub) > 0:
            app_vals.append(sub['PCC'].max())

    if myeloid_vals and app_vals:
        organ_myeloid_mean.append(np.mean(myeloid_vals))
        organ_app_mean.append(np.mean(app_vals))
        organ_labels.append(organ)

mx = np.array(organ_myeloid_mean)
ay = np.array(organ_app_mean)

# Diagonal
lim_min = min(mx.min(), ay.min()) - 0.01
lim_max = max(mx.max(), ay.max()) + 0.01
ax.plot([lim_min, lim_max], [lim_min, lim_max], 'k--', lw=1, alpha=0.5, zorder=0)

# Scatter
for i, label in enumerate(organ_labels):
    is_above = ay[i] > mx[i]  # APP > Myeloid → above diagonal
    if label == 'Colon':
        color = '#388E3C'
        marker = '*'
        ms = 200
        zorder = 5
    elif is_above:
        color = '#66BB6A'
        marker = 'o'
        ms = 60
        zorder = 3
    else:
        color = '#E57373'
        marker = 'o'
        ms = 60
        zorder = 3
    ax.scatter(mx[i], ay[i], c=color, s=ms, marker=marker, edgecolors='black',
               linewidths=0.5, zorder=zorder)
    # Label offset
    offset_x, offset_y = 0.003, 0.003
    if label == 'Colon':
        offset_y = 0.006
    ax.annotate(label, (mx[i], ay[i]), fontsize=8,
                xytext=(offset_x, offset_y), textcoords='offset fontsize',
                ha='left', va='bottom')

ax.set_xlabel('Mean myeloid pathway PCC', fontsize=11)
ax.set_ylabel('Mean APP pathway PCC', fontsize=11)
ax.set_title('Myeloid vs APP pathway per organ', fontsize=13, fontweight='bold')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Shade regions
ax.fill_between([lim_min, lim_max], [lim_min, lim_max], lim_max,
                alpha=0.05, color='green', label='APP dominant')
ax.fill_between([lim_min, lim_max], lim_min, [lim_min, lim_max],
                alpha=0.05, color='red', label='Myeloid dominant')
ax.legend(fontsize=8, loc='lower right')

# ======== Panel labels ========
for i, (label, axi) in enumerate(zip(['A', 'B', 'C', 'D'], axes.flat)):
    axi.text(-0.08, 1.08, label, transform=axi.transAxes,
             fontsize=20, fontweight='bold', va='top', ha='right')

# ======== Save ========
fig.savefig(os.path.join(OUT_DIR, 'Figure3.png'), dpi=300, bbox_inches='tight',
            facecolor='white')
fig.savefig(os.path.join(OUT_DIR, 'Figure3.pdf'), bbox_inches='tight',
            facecolor='white')
plt.close()
print('Figure 3 saved successfully.')
print(f'  PNG: {OUT_DIR}/Figure3.png')
print(f'  PDF: {OUT_DIR}/Figure3.pdf')
