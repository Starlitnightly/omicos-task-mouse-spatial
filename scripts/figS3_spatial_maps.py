#!/usr/bin/env python
"""Figure S3: Per-organ spatial risk maps (4x4 grid)"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings('ignore')

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 11

BASEDIR = ''
OUTDIR = f'{BASEDIR}figures/paper_v2'

# Load data
print("Loading spot data...")
spot_data = pd.read_csv(f'{BASEDIR}results/all_organs_spot_AD_pvalues.csv.gz')
ctrl1 = spot_data[spot_data['sample'] == 'CTRL_1'].copy()
print(f"  CTRL1 spots: {len(ctrl1)}")

cauchy = pd.read_csv(f'{BASEDIR}results/all_organs_cauchy_AD.csv.gz')
cauchy_c1 = cauchy[cauchy['sample'] == 'CTRL_1']
# Min cauchy p per organ
organ_p = cauchy_c1.groupby('organ')['p_cauchy'].min().to_dict()

organs = sorted(ctrl1['organ'].unique())
print(f"  Organs: {len(organs)}")

# Global color range
vmin = 0
vmax = ctrl1['logp'].quantile(0.99)

fig, axes = plt.subplots(4, 4, figsize=(20, 20))
axes_flat = axes.flatten()

for idx, org in enumerate(organs):
    ax = axes_flat[idx]
    organ_data = ctrl1[ctrl1['organ'] == org]

    # Subsample if too many
    if len(organ_data) > 50000:
        organ_data = organ_data.sample(50000, random_state=42)

    sc = ax.scatter(organ_data['x'], organ_data['y'],
                    c=organ_data['logp'], cmap='RdYlBu_r',
                    s=0.5, rasterized=True, vmin=vmin, vmax=vmax)
    cauchy_p = organ_p.get(org, np.nan)
    ax.set_title(f'{org}\n(Cauchy p = {cauchy_p:.2e})', fontsize=10)
    ax.set_aspect('equal')
    ax.axis('off')

# Hide unused axes
for idx in range(len(organs), 16):
    axes_flat[idx].set_visible(False)

# Shared colorbar
cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
sm = plt.cm.ScalarMappable(cmap='RdYlBu_r', norm=mcolors.Normalize(vmin=vmin, vmax=vmax))
sm.set_array([])
cbar = fig.colorbar(sm, cax=cbar_ax)
cbar.set_label('-log$_{10}$(p)', fontsize=13)

plt.subplots_adjust(left=0.02, right=0.90, top=0.96, bottom=0.02, wspace=0.05, hspace=0.15)
fig.suptitle('Spatial AD risk maps per organ (CTRL1)', fontsize=16, fontweight='bold', y=0.99)

for fmt in ['png', 'pdf']:
    fig.savefig(f'{OUTDIR}/FigureS3.{fmt}', dpi=300, bbox_inches='tight')
    print(f"Saved FigureS3.{fmt}")

import os
for fmt in ['png', 'pdf']:
    sz = os.path.getsize(f'{OUTDIR}/FigureS3.{fmt}')
    print(f"FigureS3.{fmt}: {sz/1024/1024:.1f} MB")
plt.close()
print("Done!")
