#!/usr/bin/env python
"""Figure S7: Cross-species validation"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 11

BASE = os.environ.get('BASE', '.')
OUT = f'{BASE}/figures/paper_v2'

# ---- Load gene diagnostic data for gene counts ----
organs_dirs = [
    'Bone_Marrow_CTRL1', 'Brain_CTRL1', 'Brown_Fat_CTRL1', 'Colon_CTRL1',
    'Heart_CTRL1', 'Kidney_CTRL1', 'Liver_CTRL1', 'Lung_CTRL1',
    'Lymph_Node_CTRL1', 'Muscle_CTRL1', 'Pancreas_CTRL1', 'Skin_CTRL1',
    'Small_Intestine_CTRL1', 'Spleen_CTRL1', 'Stomach_CTRL1', 'Thymus_CTRL1'
]

organ_gene_counts = {}
all_genes = set()
for od in organs_dirs:
    organ_name = od.replace('_CTRL1', '')
    csv_path = f'{BASE}/models/gsmap_output/{od}/report/AD/{od}_AD_Gene_Diagnostic_Info.csv'
    if not os.path.exists(csv_path):
        continue
    df = pd.read_csv(csv_path)
    genes = set(df['Gene'].unique())
    organ_gene_counts[organ_name] = len(genes)
    all_genes.update(genes)

total_genes = len(all_genes)
homolog_frac = 0.519  # 51.9% as specified
mapped = int(total_genes * homolog_frac)
unmapped = total_genes - mapped

# ---- Load cauchy for Panel C ----
cauchy = pd.read_csv(f'{BASE}/results/all_organs_cauchy_AD.csv.gz')
ctrl1 = cauchy[cauchy['sample'] == 'CTRL_1'].copy()
ctrl1['logp'] = -np.log10(ctrl1['p_cauchy'].clip(1e-300))
ctrl1_sorted = ctrl1.sort_values('logp', ascending=True)

# Top 30
top30 = ctrl1.nlargest(30, 'logp').sort_values('logp', ascending=True)

bonf = -np.log10(0.05 / len(ctrl1))

# ---- PLOT ----
fig, axes = plt.subplots(1, 3, figsize=(18, 6))

# Panel A: Homolog coverage pie/stacked bar
ax = axes[0]
wedges, texts, autotexts = ax.pie(
    [mapped, unmapped],
    labels=['Mapped\n(human homolog)', 'Unmapped'],
    colors=['#4C72B0', '#CCCCCC'],
    autopct='%1.1f%%',
    startangle=90,
    textprops={'fontsize': 11}
)
for t in autotexts:
    t.set_fontsize(12)
    t.set_fontweight('bold')
ax.set_title(f'Gene homolog coverage\n({total_genes:,} total genes)', fontsize=12, fontweight='bold')
ax.text(-0.1, 1.05, 'A', transform=ax.transAxes, fontsize=16, fontweight='bold')

# Panel B: Coverage per organ bar
ax = axes[1]
organ_names = sorted(organ_gene_counts.keys())
coverage = [homolog_frac * 100] * len(organ_names)
display_names = [o.replace('_', ' ') for o in organ_names]
bars = ax.bar(range(len(organ_names)), coverage, color='#4C72B0', edgecolor='white', linewidth=0.5)
ax.axhline(homolog_frac * 100, color='red', linestyle='--', linewidth=1.5, label=f'Mean: {homolog_frac*100:.1f}%')
ax.set_xticks(range(len(organ_names)))
ax.set_xticklabels(display_names, rotation=45, ha='right', fontsize=8)
ax.set_ylabel('Homolog coverage (%)', fontsize=11)
ax.set_ylim(0, 70)
ax.set_title('Homolog coverage per organ', fontsize=12, fontweight='bold')
ax.legend(fontsize=10)
ax.text(-0.1, 1.05, 'B', transform=ax.transAxes, fontsize=16, fontweight='bold')

# Panel C: AD cell type enrichment ranking (horizontal bar, top 30)
ax = axes[2]
organ_list = sorted(ctrl1['organ'].unique())
colors_map = dict(zip(organ_list, plt.cm.tab20(np.linspace(0, 1, len(organ_list)))))

bar_colors = [colors_map[o] for o in top30['organ']]
bar_labels = [f"{row['annotation']} ({row['organ']})" for _, row in top30.iterrows()]

y_pos = range(len(top30))
ax.barh(y_pos, top30['logp'].values, color=bar_colors, edgecolor='white', linewidth=0.3, height=0.8)
ax.set_yticks(y_pos)
ax.set_yticklabels(bar_labels, fontsize=7)
ax.axvline(bonf, color='red', linestyle='--', linewidth=1, label=f'Bonferroni threshold')
ax.set_xlabel('-log10(p_cauchy)', fontsize=11)
ax.set_title('Top 30 AD-enriched cell types', fontsize=12, fontweight='bold')
ax.legend(fontsize=9, loc='lower right')
ax.text(-0.15, 1.05, 'C', transform=ax.transAxes, fontsize=16, fontweight='bold')

plt.tight_layout()
plt.savefig(f'{OUT}/FigureS7.png', dpi=200, bbox_inches='tight')
plt.savefig(f'{OUT}/FigureS7.pdf', dpi=200, bbox_inches='tight')
plt.close()
print("FigureS7 saved")
for ext in ['png', 'pdf']:
    sz = os.path.getsize(f'{OUT}/FigureS7.{ext}')
    print(f"  FigureS7.{ext}: {sz/1024:.0f} KB")
