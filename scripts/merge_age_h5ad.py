"""Merge 126 age × organ CMAP-mapped pseudo-spatial samples into a single h5ad.

For each sample:
  - Load h5ad from models/cmap_output/{organ}/{sample}.h5ad
  - Read 15 traits' per-spot p-values from models/gsmap_age_output/{sample}/spatial_ldsc/
  - Attach trait p/z to obs as `p_{TRAIT}`, `z_{TRAIT}`
  - Tag with sample/organ/age in obs

Concatenate via sc.concat(..., join='inner') to keep common-gene set; X kept sparse.

Output: data/age_merged/age_all_organs_all_traits.h5ad
"""
import os, glob, sys
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy import sparse

CMAP_DIR  = 'models/cmap_output'
LDSC_ROOT = 'models/gsmap_age_output'
OUT_DIR   = 'data/age_merged'
OUT_FILE  = f'{OUT_DIR}/age_all_organs_all_traits.h5ad'

os.makedirs(OUT_DIR, exist_ok=True)

TRAITS = ['AD','PD','ALS','MS','FTD','SCZ','MDD','ADHD',
          'T2D','RA','CRP','BMD','GERD','IBS','Height']

# Enumerate samples with full trait coverage
sample_dirs = sorted(os.listdir(LDSC_ROOT))
samples = []
for s in sample_dirs:
    ad_csv = f'{LDSC_ROOT}/{s}/spatial_ldsc/{s}_AD.csv.gz'
    if os.path.exists(ad_csv):
        samples.append(s)
print(f'[1/4] Found {len(samples)} samples with gsMap AD output')

# Map sample -> CMAP h5ad path
def find_cmap(sample):
    organ = sample.split('_')[0]
    if organ in ('Bone','BAT'):  # underscore organs
        if sample.startswith('Bone_Marrow'):
            organ = 'Bone_Marrow'
    # Re-derive organ from filename match
    for o in os.listdir(CMAP_DIR):
        if not os.path.isdir(f'{CMAP_DIR}/{o}'): continue
        p = f'{CMAP_DIR}/{o}/{sample}.h5ad'
        if os.path.exists(p):
            return o, p
    return None, None

# Pre-flight: confirm all samples map to a CMAP file
pairs = []
missing = []
for s in samples:
    organ, path = find_cmap(s)
    if path is None:
        missing.append(s)
    else:
        pairs.append((s, organ, path))
print(f'[2/4] {len(pairs)}/{len(samples)} samples have a matching CMAP h5ad')
if missing:
    print(f'   missing: {missing[:5]}{"..." if len(missing)>5 else ""}')

# Load each, attach trait p/z, append
def load_one(sample, organ, path):
    a = sc.read_h5ad(path)
    a.obs['sample'] = sample
    a.obs['organ']  = organ
    # Parse age — accept both "3m" and "03_months"
    age_str = sample.replace(f'{organ}_', '')
    age_m = None
    if age_str.endswith('m'):
        try: age_m = int(age_str[:-1])
        except: pass
    elif '_months' in age_str:
        try: age_m = int(age_str.split('_')[0])
        except: pass
    a.obs['age_months'] = age_m
    a.obs['age_label']  = age_str

    # Attach trait p/z
    ldsc_dir = f'{LDSC_ROOT}/{sample}/spatial_ldsc'
    spot_index = a.obs.index.astype(str)
    for trait in TRAITS:
        csv = f'{ldsc_dir}/{sample}_{trait}.csv.gz'
        if not os.path.exists(csv):
            a.obs[f'p_{trait}'] = np.nan
            a.obs[f'z_{trait}'] = np.nan
            continue
        df = pd.read_csv(csv)
        # spot column should match obs index
        m = df.set_index(df['spot'].astype(str))
        a.obs[f'p_{trait}'] = m['p'].reindex(spot_index).values
        a.obs[f'z_{trait}'] = m['z'].reindex(spot_index).values

    # Sparsify X if dense
    if not sparse.issparse(a.X):
        a.X = sparse.csr_matrix(a.X)
    return a

merged = None
adatas = []
for i, (s, o, p) in enumerate(pairs, 1):
    try:
        a = load_one(s, o, p)
        adatas.append(a)
        if i % 20 == 0 or i == len(pairs):
            print(f'   [{i}/{len(pairs)}] loaded {s}  shape={a.shape}')
    except Exception as e:
        print(f'   FAIL {s}: {type(e).__name__}: {e}')

print(f'[3/4] Concatenating {len(adatas)} AnnData objects (inner-join on genes)...')
merged = ad.concat(adatas, join='inner', label='sample_id',
                   keys=[a.obs['sample'].iloc[0] for a in adatas],
                   index_unique='__')
# anndata's concat drops obsm by default unless aligned shapes. We need spatial coords.
# Rebuild obsm['spatial'] by stacking each sample's spatial obsm.
spatial_all = np.vstack([a.obsm['spatial'] for a in adatas])
merged.obsm['spatial'] = spatial_all
# Original sample-local obs.index may have been suffixed by index_unique='__' — keep mapping
print(f'   Merged shape: {merged.shape}')
print(f'   X sparse: {sparse.issparse(merged.X)}')
print(f'   obs cols: {list(merged.obs.columns)[:10]}... ({len(merged.obs.columns)} total)')
print(f'   obsm: {list(merged.obsm.keys())}')

# Final dtype clean-up: categoricals for organ/age/sample/cellType to save space
for c in ['organ','age_label','sample','cellType','annotation','Age_group']:
    if c in merged.obs.columns:
        merged.obs[c] = merged.obs[c].astype('category')

print(f'[4/4] Writing to {OUT_FILE}...')
merged.write_h5ad(OUT_FILE, compression='gzip')
sz_gb = os.path.getsize(OUT_FILE) / 1e9
print(f'   wrote {sz_gb:.2f} GB')

# Print a quick summary
print('\n=== SUMMARY ===')
print(f'Total spots: {merged.n_obs:,}')
print(f'Common genes (inner-join): {merged.n_vars:,}')
print(f'Samples: {merged.obs["sample"].nunique()}')
print(f'Organs: {sorted(merged.obs["organ"].unique())}')
print(f'Ages (months): {sorted(merged.obs["age_months"].dropna().unique().tolist())}')
print(f'Traits: {TRAITS}')
print(f'\nPer-organ × age spot counts:')
print(merged.obs.groupby(['organ','age_months']).size().unstack(fill_value=0))
