#!/usr/bin/env python
"""
Generate Figure 4 and Figure 5 for AD GWAS spatial mapping paper.
"""

import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

# Style
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 11
plt.rcParams['axes.titlesize'] = 14
plt.rcParams['axes.titleweight'] = 'bold'

BASEDIR = '.'
OUTDIR = f'{BASEDIR}/figures/paper_v2'

# ── Data loading and cleaning ──────────────────────────────────────────

def fix_organ(name):
    """Remove trailing _## patterns from PanSci organ names."""
    name = re.sub(r'_\d{2}$', '', str(name))
    name = re.sub(r'_\d{2}_months$', '', str(name))
    return name

def parse_age(age_str, organ_raw=None):
    """Normalize age strings to numeric months."""
    age_str = str(age_str).strip()
    # Try ##_months format
    m = re.match(r'^(\d+)_months$', age_str)
    if m:
        return int(m.group(1))
    # Try ##m format
    m = re.match(r'^(\d+)m$', age_str)
    if m:
        return int(m.group(1))
    # If age is 'months', extract from organ name (PanSci)
    if age_str == 'months' and organ_raw is not None:
        m = re.search(r'_(\d+)$', str(organ_raw))
        if m:
            return int(m.group(1))
    return None

def load_trait(trait):
    """Load and clean a trait's age-resolved data."""
    path = f'{BASEDIR}/results/age_all_cauchy_{trait}.csv.gz'
    df = pd.read_csv(path)
    if 'trait' not in df.columns:
        df['trait'] = trait
    # Parse age from age column + organ_raw
    df['organ_raw'] = df['organ'].copy()
    df['age_months'] = df.apply(lambda r: parse_age(r['age'], r['organ_raw']), axis=1)
    df['organ'] = df['organ_raw'].apply(fix_organ)
    df = df.dropna(subset=['age_months'])
    df['age_months'] = df['age_months'].astype(int)
    return df

ad = load_trait('AD')
pd_df = load_trait('PD')
als = load_trait('ALS')

print(f"AD: {len(ad)} rows, organs={sorted(ad['organ'].unique())}, ages={sorted(ad['age_months'].unique())}")
print(f"PD: {len(pd_df)} rows, organs={sorted(pd_df['organ'].unique())}, ages={sorted(pd_df['age_months'].unique())}")
print(f"ALS: {len(als)} rows, organs={sorted(als['organ'].unique())}, ages={sorted(als['age_months'].unique())}")

# Anatomical organ order
ORGAN_ORDER = [
    'Brain', 'Spleen', 'Thymus', 'Bone_Marrow',
    'Lung', 'Heart', 'Liver', 'Kidney', 'Pancreas',
    'Stomach', 'Ileum', 'Colon',
    'Muscle', 'Skin', 'BAT'
]

def organ_sort_key(organ):
    if organ in ORGAN_ORDER:
        return ORGAN_ORDER.index(organ)
    return len(ORGAN_ORDER)

# ── Figure 4 ───────────────────────────────────────────────────────────

def make_figure4():
    # Panel A & B data: organ × age heatmap
    grouped = ad.groupby(['organ', 'age_months'])['p_cauchy'].min().reset_index()
    grouped['neg_log10_p'] = -np.log10(grouped['p_cauchy'].clip(lower=1e-300))
    # For duplicate organ-age (from different datasets), take min p
    grouped = grouped.groupby(['organ', 'age_months']).agg({'p_cauchy': 'min', 'neg_log10_p': 'max'}).reset_index()

    organs = sorted(grouped['organ'].unique(), key=organ_sort_key)
    ages = sorted(grouped['age_months'].unique())

    # Build matrix
    mat = pd.DataFrame(np.nan, index=organs, columns=ages)
    pmat = pd.DataFrame(1.0, index=organs, columns=ages)
    for _, row in grouped.iterrows():
        mat.loc[row['organ'], row['age_months']] = row['neg_log10_p']
        pmat.loc[row['organ'], row['age_months']] = row['p_cauchy']

    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    fig.subplots_adjust(wspace=0.35)

    # Panel A: Heatmap
    ax = axes[0]
    # Mask non-significant (p > 0.05) as 0
    mask = pmat > 0.05
    mat_plot = mat.copy()
    mat_plot[mask] = 0

    # Annotation: show values, but mark NS
    annot_arr = mat.copy().astype(object)
    for o in organs:
        for a in ages:
            v = mat.loc[o, a]
            p = pmat.loc[o, a]
            if pd.isna(v):
                annot_arr.loc[o, a] = ''
            elif p > 0.05:
                annot_arr.loc[o, a] = 'NS'
            else:
                annot_arr.loc[o, a] = f'{v:.1f}'

    # Replace NaN with 0 for plotting
    mat_plot = mat_plot.fillna(0)

    sns.heatmap(mat_plot, ax=ax, cmap='YlOrRd', annot=annot_arr, fmt='',
                linewidths=0.5, linecolor='white', cbar_kws={'label': '-log10(p)', 'shrink': 0.7},
                vmin=0, vmax=mat.max().max())
    ax.set_xlabel('Age (months)', fontweight='bold')
    ax.set_ylabel('')
    ax.set_title('AD enrichment across organs and ages')
    ax.set_xticklabels([str(int(a)) for a in ages], rotation=0)
    ax.set_yticklabels(organs, rotation=0)

    # Panel label
    ax.text(-0.12, 1.05, 'A', transform=ax.transAxes, fontsize=20, fontweight='bold', va='top')

    # Panel B: Trajectory line plot
    ax2 = axes[1]
    cmap_lines = plt.cm.get_cmap('tab20', len(organs))
    for i, organ in enumerate(organs):
        sub = grouped[grouped['organ'] == organ].sort_values('age_months')
        if len(sub) > 0:
            ax2.plot(sub['age_months'], sub['neg_log10_p'], 'o-', color=cmap_lines(i),
                     label=organ, linewidth=1.5, markersize=4)

    # Bonferroni threshold
    n_tests = ad['annotation'].nunique()
    bonf = -np.log10(0.05 / n_tests)
    ax2.axhline(bonf, color='red', linestyle='--', linewidth=1, alpha=0.7, label=f'Bonferroni (p={0.05/n_tests:.1e})')
    ax2.axhline(-np.log10(0.05), color='gray', linestyle=':', linewidth=1, alpha=0.5, label='p=0.05')

    ax2.set_xlabel('Age (months)', fontweight='bold')
    ax2.set_ylabel('-log10(min p_cauchy)', fontweight='bold')
    ax2.set_title('AD risk trajectory across lifespan')
    ax2.set_xticks(ages)
    ax2.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=8, frameon=False)

    ax2.text(-0.12, 1.05, 'B', transform=ax2.transAxes, fontsize=20, fontweight='bold', va='top')

    fig.suptitle('Figure 4: AD risk is stable across the mouse lifespan', fontsize=16, fontweight='bold', y=1.02)
    fig.tight_layout()

    for ext in ['png', 'pdf']:
        fig.savefig(f'{OUTDIR}/Figure4.{ext}', dpi=300, bbox_inches='tight')
    plt.close(fig)
    print("Figure 4 saved.")

# ── Figure 5 ───────────────────────────────────────────────────────────

def make_figure5():
    all_data = pd.concat([ad, pd_df, als], ignore_index=True)

    # Panel A: Organ-level comparison
    # For each trait+organ, min p_cauchy across all ages and annotations
    organ_trait = all_data.groupby(['trait', 'organ'])['p_cauchy'].min().reset_index()
    organ_trait['neg_log10_p'] = -np.log10(organ_trait['p_cauchy'].clip(lower=1e-300))

    # Pivot
    pivot_a = organ_trait.pivot(index='organ', columns='trait', values='neg_log10_p').fillna(0)
    # Only keep organs present in at least AD
    ad_organs = set(ad['organ'].unique())
    pivot_a = pivot_a[pivot_a.index.isin(ad_organs)]
    for t in ['AD', 'PD', 'ALS']:
        if t not in pivot_a.columns:
            pivot_a[t] = 0
    pivot_a = pivot_a.sort_values('AD', ascending=True)

    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.subplots_adjust(hspace=0.35, wspace=0.35)

    # Panel A: Grouped horizontal bars
    ax = axes[0, 0]
    y = np.arange(len(pivot_a))
    h = 0.25
    ax.barh(y - h, pivot_a['AD'], height=h, color='#d62728', label='AD', zorder=3)
    ax.barh(y, pivot_a['PD'], height=h, color='#1f77b4', label='PD', zorder=3)
    ax.barh(y + h, pivot_a['ALS'], height=h, color='#2ca02c', label='ALS', zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(pivot_a.index)
    ax.set_xlabel('-log10(min p_cauchy)', fontweight='bold')
    ax.set_title('Organ-level enrichment: AD vs PD vs ALS')
    # Bonferroni
    n_all = all_data['annotation'].nunique()
    bonf = -np.log10(0.05 / n_all)
    ax.axvline(bonf, color='red', linestyle='--', linewidth=1, alpha=0.7, label='Bonferroni')
    ax.axvline(-np.log10(0.05), color='gray', linestyle=':', linewidth=1, alpha=0.5, label='p=0.05')
    ax.legend(fontsize=9, loc='lower right')
    ax.grid(axis='x', alpha=0.3)
    ax.text(-0.15, 1.05, 'A', transform=ax.transAxes, fontsize=20, fontweight='bold', va='top')

    # Panel B: Myeloid enrichment heatmap
    ax = axes[0, 1]
    myeloid_kw = ['myeloid', 'macro', 'micro', 'mono', 'microglia']
    is_myeloid = all_data['annotation'].str.lower().apply(
        lambda x: any(k in x for k in myeloid_kw))
    myeloid_data = all_data[is_myeloid].copy()
    my_ot = myeloid_data.groupby(['trait', 'organ'])['p_cauchy'].min().reset_index()
    my_ot['neg_log10_p'] = -np.log10(my_ot['p_cauchy'].clip(lower=1e-300))
    pivot_b = my_ot.pivot(index='organ', columns='trait', values='neg_log10_p').fillna(0)
    pivot_b = pivot_b[pivot_b.index.isin(ad_organs)]
    for t in ['AD', 'PD', 'ALS']:
        if t not in pivot_b.columns:
            pivot_b[t] = 0
    pivot_b = pivot_b[['AD', 'ALS', 'PD']]
    pivot_b = pivot_b.sort_values('AD', ascending=False)

    # Only keep organs with some signal
    pivot_b = pivot_b[(pivot_b > 0).any(axis=1)]

    sns.heatmap(pivot_b, ax=ax, cmap='YlOrRd', annot=True, fmt='.1f',
                linewidths=0.5, linecolor='white', cbar_kws={'label': '-log10(p)', 'shrink': 0.7})
    ax.set_title('Myeloid enrichment: AD >> ALS > PD')
    ax.set_ylabel('')
    ax.text(-0.15, 1.05, 'B', transform=ax.transAxes, fontsize=20, fontweight='bold', va='top')

    # Panel C: Brain cell types AD vs PD
    ax = axes[1, 0]
    brain_ad = ad[ad['organ'] == 'Brain'].groupby('annotation')['p_cauchy'].min().reset_index()
    brain_ad['neg_log10_p'] = -np.log10(brain_ad['p_cauchy'].clip(lower=1e-300))
    brain_ad['trait'] = 'AD'

    brain_pd = pd_df[pd_df['organ'] == 'Brain'].groupby('annotation')['p_cauchy'].min().reset_index()
    brain_pd['neg_log10_p'] = -np.log10(brain_pd['p_cauchy'].clip(lower=1e-300))
    brain_pd['trait'] = 'PD'

    brain_both = pd.merge(brain_ad[['annotation', 'neg_log10_p']],
                          brain_pd[['annotation', 'neg_log10_p']],
                          on='annotation', suffixes=('_AD', '_PD'), how='outer').fillna(0)

    # Classify cell types
    def classify_brain(name):
        nl = name.lower()
        if any(k in nl for k in myeloid_kw):
            return 'Myeloid'
        if 'oligodendro' in nl:
            return 'Oligodendrocyte'
        return 'Other'

    brain_both['cell_class'] = brain_both['annotation'].apply(classify_brain)
    brain_both = brain_both.sort_values('neg_log10_p_AD', ascending=True)

    # Top 20 for readability
    brain_both_top = brain_both.nlargest(20, 'neg_log10_p_AD').sort_values('neg_log10_p_AD', ascending=True)

    y = np.arange(len(brain_both_top))
    colors_map = {'Myeloid': '#d62728', 'Oligodendrocyte': '#9467bd', 'Other': '#7f7f7f'}
    colors_ad = [colors_map[c] for c in brain_both_top['cell_class']]
    colors_pd = ['#1f77b4'] * len(brain_both_top)

    ax.barh(y - 0.15, brain_both_top['neg_log10_p_AD'], height=0.3, color=colors_ad, label='AD', zorder=3)
    ax.barh(y + 0.15, brain_both_top['neg_log10_p_PD'], height=0.3, color=colors_pd, alpha=0.7, label='PD', zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(brain_both_top['annotation'], fontsize=8)
    ax.set_xlabel('-log10(min p_cauchy)', fontweight='bold')
    ax.set_title('Brain: AD vs PD cell types')
    # Custom legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#d62728', label='AD - Myeloid'),
        Patch(facecolor='#9467bd', label='AD - Oligodendrocyte'),
        Patch(facecolor='#7f7f7f', label='AD - Other'),
        Patch(facecolor='#1f77b4', alpha=0.7, label='PD'),
    ]
    ax.legend(handles=legend_elements, fontsize=8, loc='lower right')
    ax.grid(axis='x', alpha=0.3)
    ax.text(-0.15, 1.05, 'C', transform=ax.transAxes, fontsize=20, fontweight='bold', va='top')

    # Panel D: Gut comparison AD vs PD
    ax = axes[1, 1]
    gut_organs = ['Colon', 'Ileum', 'Stomach']

    immune_kw = ['myeloid', 'lymphoid', 'immune', 'macro', 'mono', 'micro', 'b cell', 't cell',
                 'plasma', 'dendritic', 'neutrophil', 'hematopoietic', 'leukocyte']
    epithelial_kw = ['epithelial', 'enterocyte', 'goblet', 'paneth', 'tuft', 'brush',
                     'enteroendocrine', 'stem cell', 'crypt', 'mucous', 'parietal', 'chief',
                     'secretory', 'intestinal']

    def classify_gut(name):
        nl = name.lower()
        if any(k in nl for k in immune_kw):
            return 'Immune'
        if any(k in nl for k in epithelial_kw):
            return 'Epithelial'
        return 'Other'

    gut_data = all_data[all_data['organ'].isin(gut_organs) & all_data['trait'].isin(['AD', 'PD'])].copy()
    gut_data['cell_class'] = gut_data['annotation'].apply(classify_gut)
    gut_data = gut_data[gut_data['cell_class'].isin(['Immune', 'Epithelial'])]

    gut_summary = gut_data.groupby(['organ', 'trait', 'cell_class'])['p_cauchy'].min().reset_index()
    gut_summary['neg_log10_p'] = -np.log10(gut_summary['p_cauchy'].clip(lower=1e-300))

    # Create grouped bars: organ × cell_class, comparing AD vs PD
    categories = []
    ad_vals = []
    pd_vals = []
    for organ in gut_organs:
        for cc in ['Immune', 'Epithelial']:
            label = f'{organ}\n{cc}'
            categories.append(label)
            ad_v = gut_summary[(gut_summary['organ'] == organ) &
                               (gut_summary['trait'] == 'AD') &
                               (gut_summary['cell_class'] == cc)]['neg_log10_p']
            pd_v = gut_summary[(gut_summary['organ'] == organ) &
                               (gut_summary['trait'] == 'PD') &
                               (gut_summary['cell_class'] == cc)]['neg_log10_p']
            ad_vals.append(ad_v.values[0] if len(ad_v) > 0 else 0)
            pd_vals.append(pd_v.values[0] if len(pd_v) > 0 else 0)

    x = np.arange(len(categories))
    w = 0.35
    ax.bar(x - w/2, ad_vals, w, color='#d62728', label='AD', zorder=3)
    ax.bar(x + w/2, pd_vals, w, color='#1f77b4', label='PD', zorder=3)
    ax.set_xticks(x)
    ax.set_xticklabels(categories, fontsize=9)
    ax.set_ylabel('-log10(min p_cauchy)', fontweight='bold')
    ax.set_title('Gut: AD vs PD enrichment')
    ax.axhline(-np.log10(0.05), color='gray', linestyle=':', linewidth=1, alpha=0.5)
    ax.legend(fontsize=9)
    ax.grid(axis='y', alpha=0.3)
    ax.text(-0.15, 1.05, 'D', transform=ax.transAxes, fontsize=20, fontweight='bold', va='top')

    fig.suptitle('Figure 5: Multi-trait comparison', fontsize=16, fontweight='bold', y=1.02)
    fig.tight_layout()

    for ext in ['png', 'pdf']:
        fig.savefig(f'{OUTDIR}/Figure5.{ext}', dpi=300, bbox_inches='tight')
    plt.close(fig)
    print("Figure 5 saved.")

# ── Execute ────────────────────────────────────────────────────────────

make_figure4()
make_figure5()

# Report sizes
import os
for f in ['Figure4.png', 'Figure4.pdf', 'Figure5.png', 'Figure5.pdf']:
    path = f'{OUTDIR}/{f}'
    size = os.path.getsize(path)
    print(f"  {f}: {size/1024:.0f} KB")
