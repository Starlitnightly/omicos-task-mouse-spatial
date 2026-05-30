"""Aggregate gsMap results from all completed organs into combined CSV files."""
import os, glob
import numpy as np
import pandas as pd
import anndata as ad

BASE = os.environ.get('BASE', '')
WORKDIR = f'{BASE}models/gsmap_output' if BASE.endswith('/') or BASE == '' else f'{BASE}/models/gsmap_output'
ST_DIR = f'{BASE}data/st/per_organ' if BASE.endswith('/') or BASE == '' else f'{BASE}/data/st/per_organ'
RESULTS_DIR = f'{BASE}results' if BASE.endswith('/') or BASE == '' else f'{BASE}/results'

# Collect spatial LDSC spot-level results
all_spot = []
all_cauchy = []

for ldsc_file in sorted(glob.glob(f'{WORKDIR}/*/spatial_ldsc/*_AD.csv.gz')):
    sample_name = ldsc_file.split('/')[-3]
    h5ad_path = f'{ST_DIR}/{sample_name}.h5ad'

    if not os.path.exists(h5ad_path):
        continue

    df = pd.read_csv(ldsc_file)
    adata = ad.read_h5ad(h5ad_path, backed='r')
    organ = str(adata.obs['Organ_Full_Name'].iloc[0])
    sample = str(adata.obs['Sample'].iloc[0])

    df['organ'] = organ
    df['sample'] = sample
    df['sample_name'] = sample_name
    df['x'] = adata.obsm['spatial'][:len(df), 0]
    df['y'] = adata.obsm['spatial'][:len(df), 1]
    all_spot.append(df)

    # Cauchy results
    cauchy_file = f'{WORKDIR}/{sample_name}/cauchy_combination/{sample_name}_AD.Cauchy.csv.gz'
    if os.path.exists(cauchy_file):
        dfc = pd.read_csv(cauchy_file)
        dfc['organ'] = organ
        dfc['sample'] = sample
        dfc['sample_name'] = sample_name
        all_cauchy.append(dfc)

if all_spot:
    spot_df = pd.concat(all_spot, ignore_index=True)
    # Add logp column
    spot_df['logp'] = -np.log10(spot_df['p'].clip(lower=1e-300))
    spot_df.to_csv(f'{RESULTS_DIR}/all_organs_spot_AD_pvalues.csv.gz', index=False, compression='gzip')
    print(f"Spot results: {len(spot_df)} spots, {spot_df['organ'].nunique()} organs")
    print(spot_df.groupby('organ').size().sort_values(ascending=False))

if all_cauchy:
    cauchy_df = pd.concat(all_cauchy, ignore_index=True)
    cauchy_df.to_csv(f'{RESULTS_DIR}/all_organs_cauchy_AD.csv.gz', index=False, compression='gzip')
    print(f"\nCauchy results: {len(cauchy_df)} annotations")
    print(cauchy_df.sort_values('p_cauchy').head(20))
