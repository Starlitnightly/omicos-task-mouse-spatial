"""Merge the 32 original 6-week-old Array-seq pseudo-spatial samples (16 organs x 2 CTRL
replicates) into a single h5ad, attaching per-spot p-values from 15 traits' gsMap output.

Source files:
  - h5ad:  data/st/per_organ/{Organ}_CTRL{1|2}.h5ad
  - gsMap: models/gsmap_output/{Organ}_CTRL{1|2}/spatial_ldsc/{sample}_{TRAIT}.csv.gz

Output: data/age_merged/original_st_all_organs_all_traits.h5ad
"""
import os, glob
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy import sparse

ST_DIR    = 'data/st/per_organ'
LDSC_ROOT = 'models/gsmap_output'
OUT_DIR   = 'data/age_merged'
OUT_FILE  = f'{OUT_DIR}/original_st_all_organs_all_traits.h5ad'

os.makedirs(OUT_DIR, exist_ok=True)

TRAITS = ['AD','PD','ALS','MS','FTD','SCZ','MDD','ADHD',
          'T2D','RA','CRP','BMD','GERD','IBS','Height']

# Enumerate samples — every {Organ}_CTRL{1|2}.h5ad
h5ads = sorted(glob.glob(f'{ST_DIR}/*.h5ad'))
samples = []
for p in h5ads:
    s = os.path.basename(p).replace('.h5ad','')
    ldsc_csv = f'{LDSC_ROOT}/{s}/spatial_ldsc/{s}_AD.csv.gz'
    if not os.path.exists(ldsc_csv):
        print(f'   SKIP {s}: missing AD gsMap output')
        continue
    # Parse organ from {Organ}_CTRL{R}
    parts = s.split('_CTRL')
    organ = parts[0]
    rep   = f'CTRL{parts[1]}' if len(parts) > 1 else 'unknown'
    samples.append((s, organ, rep, p))
print(f'[1/4] Found {len(samples)} samples with both h5ad and gsMap output')

def load_one(sample, organ, rep, path):
    a = sc.read_h5ad(path)
    a.obs['sample']    = sample
    a.obs['organ']     = organ
    a.obs['replicate'] = rep

    # Trait p/z per spot — sample obs index should match LDSC 'spot' column
    ldsc_dir = f'{LDSC_ROOT}/{sample}/spatial_ldsc'
    spot_index = a.obs.index.astype(str)
    for trait in TRAITS:
        csv = f'{ldsc_dir}/{sample}_{trait}.csv.gz'
        if not os.path.exists(csv):
            a.obs[f'p_{trait}'] = np.nan
            a.obs[f'z_{trait}'] = np.nan
            continue
        df = pd.read_csv(csv)
        m = df.set_index(df['spot'].astype(str))
        a.obs[f'p_{trait}'] = m['p'].reindex(spot_index).values
        a.obs[f'z_{trait}'] = m['z'].reindex(spot_index).values

    if not sparse.issparse(a.X):
        a.X = sparse.csr_matrix(a.X)
    return a

adatas = []
for i, (s, o, r, p) in enumerate(samples, 1):
    try:
        a = load_one(s, o, r, p)
        adatas.append(a)
        print(f'   [{i:>2}/{len(samples)}] {s:<28} shape={a.shape}')
    except Exception as e:
        print(f'   FAIL {s}: {type(e).__name__}: {e}')

print(f'[2/4] Concatenating {len(adatas)} AnnData (inner-join on genes)...')
merged = ad.concat(adatas, join='inner', label='sample_id',
                   keys=[a.obs['sample'].iloc[0] for a in adatas],
                   index_unique='__')
spatial_all = np.vstack([a.obsm['spatial'] for a in adatas])
merged.obsm['spatial'] = spatial_all
print(f'   Merged shape: {merged.shape}')
print(f'   X sparse: {sparse.issparse(merged.X)}')

# Categoricals
for c in ['organ','replicate','sample']:
    if c in merged.obs.columns:
        merged.obs[c] = merged.obs[c].astype('category')

print(f'[3/4] Writing to {OUT_FILE}...')
merged.write_h5ad(OUT_FILE, compression='gzip')
sz_gb = os.path.getsize(OUT_FILE) / 1e9
print(f'   wrote {sz_gb:.2f} GB')

print('\n=== SUMMARY ===')
print(f'Total spots: {merged.n_obs:,}')
print(f'Common genes (inner-join): {merged.n_vars:,}')
print(f'Samples: {merged.obs["sample"].nunique()}')
print(f'Organs: {sorted(merged.obs["organ"].unique())}')
print(f'\nPer-organ spot counts (CTRL1 + CTRL2):')
print(merged.obs.groupby(['organ','replicate'], observed=True).size().unstack(fill_value=0))

# trait sanity
print(f'\n=== Per-trait median -log10(p) (lower = no enrichment, ~0.3 = noise floor) ===')
for tr in TRAITS:
    col = f'p_{tr}'
    if col in merged.obs.columns:
        med = -np.log10(merged.obs[col].clip(lower=1e-30)).median()
        print(f'  {tr:>6}: {med:.3f}')

# Per-organ AD median
print(f'\n=== AD per-organ median -log10(p) ===')
adv = merged.obs[['organ','p_AD']].copy()
adv['neglogp'] = -np.log10(adv['p_AD'].clip(lower=1e-30))
print(adv.groupby('organ', observed=True)['neglogp'].median().sort_values(ascending=False))
print(f'\n[4/4] Done.')
