#!/usr/bin/env python
"""Figure 7: Peripheral-brain axis for Cell paper."""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import pearsonr
import warnings
warnings.filterwarnings('ignore')

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 11

BASE = Path('')
OUT = BASE / 'figures' / 'paper_v2'

ORGAN_ORDER = ['Brain', 'Bone_Marrow', 'Thymus', 'Spleen', 'Lymph_Node',
               'Lung', 'Heart', 'Liver', 'Kidney', 'Pancreas',
               'Stomach', 'Small_Intestine', 'Colon', 'Muscle', 'Skin', 'Brown_Fat']

def load_all_organs():
    all_data = {}
    for organ in ORGAN_ORDER:
        csv_path = (BASE / 'models' / 'gsmap_output' / f'{organ}_CTRL1' /
                    'report' / 'AD' / f'{organ}_CTRL1_AD_Gene_Diagnostic_Info.csv')
        if csv_path.exists():
            df = pd.read_csv(csv_path)
            all_data[organ] = df
    return all_data

all_organs = load_all_organs()
organs_present = [o for o in ORGAN_ORDER if o in all_organs]

# Build gene x organ max PCC pivot
records = []
for organ, df in all_organs.items():
    max_pcc = df.groupby('Gene')['PCC'].max().reset_index()
    max_pcc['Organ'] = organ
    records.append(max_pcc)
gene_organ_df = pd.concat(records, ignore_index=True)
pivot = gene_organ_df.pivot_table(index='Gene', columns='Organ', values='PCC', aggfunc='max')

# ══════════════════════════════════════════════
# Figure 7
# ══════════════════════════════════════════════
fig, axes = plt.subplots(2, 2, figsize=(16, 14))
fig.subplots_adjust(hspace=0.4, wspace=0.35, left=0.08, right=0.95, top=0.95, bottom=0.06)

# ── Panel A: Cross-organ correlation ──
ax_a = axes[0, 0]

# Filter genes present in at least 8 organs
gene_counts = pivot[organs_present].notna().sum(axis=1)
genes_8plus = gene_counts[gene_counts >= 8].index
pivot_filt = pivot.loc[genes_8plus, organs_present].copy()

# Compute organ-organ Pearson correlation using DataFrame.corr()
n_org = len(organs_present)
corr_mat = pivot_filt[organs_present].corr(method='pearson', min_periods=10)

# Lower triangle mask
mask = np.triu(np.ones((n_org, n_org), dtype=bool), k=1)
corr_vals = corr_mat.values.astype(float).copy()
corr_vals[mask] = np.nan

im_a = ax_a.imshow(corr_vals, cmap='RdBu_r', vmin=-0.5, vmax=1, aspect='auto')
ax_a.set_facecolor('#e0e0e0')

organ_labels = [o.replace('_', ' ') for o in organs_present]
ax_a.set_xticks(range(n_org))
ax_a.set_xticklabels(organ_labels, rotation=45, ha='right', fontsize=7)
ax_a.set_yticks(range(n_org))
ax_a.set_yticklabels(organ_labels, fontsize=7)

# Annotate lower triangle
for i in range(n_org):
    for j in range(n_org):
        if not mask[i, j] and not np.isnan(corr_vals[i, j]):
            val = corr_vals[i, j]
            color = 'white' if abs(val) > 0.5 else 'black'
            ax_a.text(j, i, f'{val:.2f}', ha='center', va='center',
                     fontsize=5, color=color)

plt.colorbar(im_a, ax=ax_a, shrink=0.6, label='Pearson r')
ax_a.set_title('Cross-organ correlation of AD risk gene profiles', fontsize=10, fontweight='bold')
ax_a.text(-0.12, 1.05, 'A', transform=ax_a.transAxes, fontsize=20, fontweight='bold', va='top')

# ── Panel B: Brain-peripheral shared genes matrix ──
ax_b = axes[0, 1]

# Top 100 genes per organ by PCC
top100 = {}
for organ, df in all_organs.items():
    t = df.loc[df.groupby('Gene')['PCC'].idxmax()].nlargest(100, 'PCC')
    top100[organ] = set(t['Gene'].values)

brain_top100 = top100['Brain']
peripheral_organs = [o for o in organs_present if o != 'Brain']

# Find shared genes
shared_genes_set = set()
for org in peripheral_organs:
    shared_genes_set.update(brain_top100 & top100[org])

shared_genes = sorted(shared_genes_set)

# Build matrix
share_mat = pd.DataFrame(0, index=shared_genes,
                         columns=[o.replace('_', ' ') for o in peripheral_organs])
for org in peripheral_organs:
    overlap = brain_top100 & top100[org]
    for g in overlap:
        share_mat.loc[g, org.replace('_', ' ')] = 1

# Sort genes by total sharing
share_mat['total'] = share_mat.sum(axis=1)
share_mat = share_mat.sort_values('total', ascending=False)
total_col = share_mat.pop('total')

# Show top 30 for readability
n_show = min(30, len(share_mat))
share_show = share_mat.iloc[:n_show]

im_b = ax_b.imshow(share_show.values, cmap='YlGnBu', aspect='auto',
                   interpolation='nearest', vmin=0, vmax=1)
ax_b.set_xticks(range(len(peripheral_organs)))
ax_b.set_xticklabels(share_show.columns, rotation=45, ha='right', fontsize=6)
ax_b.set_yticks(range(n_show))
ax_b.set_yticklabels(share_show.index, fontsize=6)
ax_b.set_title('AD risk genes shared between brain and peripheral organs',
               fontsize=10, fontweight='bold')
ax_b.text(-0.12, 1.05, 'B', transform=ax_b.transAxes, fontsize=20, fontweight='bold', va='top')

# ── Panel C: Number of organs sharing each brain gene ──
ax_c = axes[1, 0]

# Count how many peripheral organs share each brain top100 gene
gene_organ_count = {}
for g in brain_top100:
    count = sum(1 for org in peripheral_organs if g in top100[org])
    if count > 0:
        gene_organ_count[g] = count

count_df = pd.DataFrame({'Gene': list(gene_organ_count.keys()),
                          'Count': list(gene_organ_count.values())})
count_df = count_df.sort_values('Count', ascending=False).head(20)
count_df = count_df.sort_values('Count', ascending=True)

colors_c = []
for c in count_df['Count']:
    if c >= 5:
        colors_c.append('#d32f2f')
    elif c >= 3:
        colors_c.append('#1976d2')
    else:
        colors_c.append('#9e9e9e')

ax_c.barh(range(len(count_df)), count_df['Count'].values, color=colors_c,
          edgecolor='gray', linewidth=0.5)
ax_c.set_yticks(range(len(count_df)))
ax_c.set_yticklabels(count_df['Gene'].values, fontsize=8)
ax_c.set_xlabel('Number of peripheral organs', fontsize=10)
ax_c.set_title('Brain AD genes in peripheral organs', fontsize=10, fontweight='bold')

# Legend
from matplotlib.patches import Patch
legend_c = [Patch(facecolor='#d32f2f', label='>=5 organs'),
            Patch(facecolor='#1976d2', label='>=3 organs'),
            Patch(facecolor='#9e9e9e', label='<3 organs')]
ax_c.legend(handles=legend_c, loc='lower right', fontsize=8, framealpha=0.9)
ax_c.text(-0.12, 1.05, 'C', transform=ax_c.transAxes, fontsize=20, fontweight='bold', va='top')

# ── Panel D: Top 20 universally AD-correlated genes ──
ax_d = axes[1, 1]

# Mean PCC across organs
pivot_all = pivot[organs_present].copy()
mean_pcc = pivot_all.mean(axis=1, skipna=True)
top20_genes = mean_pcc.nlargest(20).index.tolist()

heatmap_data = pivot_all.loc[top20_genes, organs_present].astype(float)

im_d = ax_d.imshow(heatmap_data.values, cmap='YlOrRd', aspect='auto')
ax_d.set_facecolor('#d0d0d0')

# Annotate
for i in range(heatmap_data.shape[0]):
    for j in range(heatmap_data.shape[1]):
        val = heatmap_data.values[i, j]
        if not np.isnan(val):
            color = 'white' if val > 0.5 * np.nanmax(heatmap_data.values) else 'black'
            ax_d.text(j, i, f'{val:.2f}', ha='center', va='center',
                     fontsize=5, color=color)

organ_labels_d = [o.replace('_', ' ') for o in organs_present]
ax_d.set_xticks(range(len(organs_present)))
ax_d.set_xticklabels(organ_labels_d, rotation=45, ha='right', fontsize=7)
ax_d.set_yticks(range(len(top20_genes)))
ax_d.set_yticklabels(top20_genes, fontsize=8)
plt.colorbar(im_d, ax=ax_d, shrink=0.6, label='PCC')
ax_d.set_title('Top 20 universally AD-correlated genes', fontsize=10, fontweight='bold')
ax_d.text(-0.12, 1.05, 'D', transform=ax_d.transAxes, fontsize=20, fontweight='bold', va='top')

print("Saving Figure 7...")
fig.savefig(OUT / 'Figure7.png', dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(OUT / 'Figure7.pdf', bbox_inches='tight', facecolor='white')
plt.close()
print("Figure 7 done.")
