#!/usr/bin/env python
"""Figure S8: Age progression across organs"""
import os, re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 11

BASE = ''
OUT = f'{BASE}figures/paper_v2'

def parse_age_months(age_str):
    age_str = str(age_str)
    m = re.match(r'(\d+)m$', age_str)
    if m: return int(m.group(1))
    m = re.match(r'(\d+)_months$', age_str)
    if m: return int(m.group(1))
    return None

def fix_organ(name):
    return re.sub(r'_\d{2}$', '', str(name))

# Load data
age_df = pd.read_csv(f'{BASE}results/age_all_cauchy_AD.csv.gz')
print(f"Age data shape: {age_df.shape}")
print(f"Organs: {age_df['organ'].unique()}")
print(f"Ages: {age_df['age'].unique()}")

# Fix organ names and parse ages
age_df['organ_clean'] = age_df['organ'].apply(fix_organ)
age_df['age_months'] = age_df['age'].apply(parse_age_months)

# Drop rows where age couldn't be parsed
n_before = len(age_df)
age_df = age_df.dropna(subset=['age_months'])
age_df['age_months'] = age_df['age_months'].astype(int)
print(f"Dropped {n_before - len(age_df)} rows with unparseable age")

organs_order = ['BAT', 'Bone_Marrow', 'Brain', 'Colon', 'Heart', 'Ileum',
                'Kidney', 'Liver', 'Lung', 'Muscle', 'Pancreas', 'Skin',
                'Spleen', 'Stomach', 'Thymus']

# Bonferroni: use total annotations across all ages/organs
n_tests = len(age_df)
bonf = -np.log10(0.05 / n_tests)

fig, axes = plt.subplots(5, 3, figsize=(24, 40))
axes = axes.flatten()

for idx, organ in enumerate(organs_order):
    ax = axes[idx]
    sub = age_df[age_df['organ_clean'] == organ].copy()

    if len(sub) == 0:
        ax.set_title(organ.replace('_', ' '), fontsize=13, fontweight='bold')
        ax.text(0.5, 0.5, 'No data', transform=ax.transAxes, ha='center', va='center', fontsize=14, color='gray')
        ax.set_xlabel('Age (months)')
        ax.set_ylabel('-log10(p_cauchy)')
        continue

    # For each age, find the annotation with min p_cauchy (most significant)
    # First, deduplicate: for same organ-age from different samples, take min p_cauchy per annotation
    sub_agg = sub.groupby(['age_months', 'annotation']).agg({'p_cauchy': 'min'}).reset_index()

    # For each age, find the best (min p) annotation
    best_per_age = sub_agg.loc[sub_agg.groupby('age_months')['p_cauchy'].idxmin()]
    best_per_age = best_per_age.sort_values('age_months')

    ages = best_per_age['age_months'].values
    logp = -np.log10(best_per_age['p_cauchy'].values.clip(1e-300))
    labels = best_per_age['annotation'].values

    bar_colors = plt.cm.viridis(np.linspace(0.3, 0.9, len(ages)))
    bars = ax.bar(range(len(ages)), logp, color=bar_colors, edgecolor='white', linewidth=0.5)
    ax.set_xticks(range(len(ages)))
    ax.set_xticklabels([str(a) for a in ages], fontsize=10)

    # Label best cell type above each bar
    for i, (b, lbl) in enumerate(zip(bars, labels)):
        # Truncate long labels
        short_lbl = lbl if len(lbl) < 20 else lbl[:18] + '..'
        ax.text(b.get_x() + b.get_width()/2, b.get_height() + 0.1,
                short_lbl, ha='center', va='bottom', fontsize=7, rotation=45)

    ax.axhline(bonf, color='red', linestyle='--', linewidth=1, alpha=0.7)
    ax.set_xlabel('Age (months)', fontsize=11)
    ax.set_ylabel('-log10(p_cauchy)', fontsize=11)
    ax.set_title(organ.replace('_', ' '), fontsize=13, fontweight='bold')

    # Add some padding at top for labels
    if len(logp) > 0:
        ax.set_ylim(0, max(logp) * 1.4)

plt.suptitle('AD enrichment across ages by organ', fontsize=16, fontweight='bold', y=0.995)
plt.tight_layout(rect=[0, 0, 1, 0.99])

plt.savefig(f'{OUT}/FigureS8.png', dpi=150, bbox_inches='tight')
plt.savefig(f'{OUT}/FigureS8.pdf', dpi=150, bbox_inches='tight')
plt.close()
print("FigureS8 saved")
for ext in ['png', 'pdf']:
    sz = os.path.getsize(f'{OUT}/FigureS8.{ext}')
    print(f"  FigureS8.{ext}: {sz/1024:.0f} KB")
