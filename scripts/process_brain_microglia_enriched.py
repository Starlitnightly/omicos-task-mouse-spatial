#!/usr/bin/env python
"""
Merge all brain data sources with enriched microglia, run CMAP per age.

Sources:
1. PanSci Brain (3m/6m/21m) - full brain, ~30K cells, ~338 microglia
2. TMS Brain_Myeloid (3m/18m/24m) - 13K pure microglia+macrophage
3. TMS Brain_NonMyeloid (3m/18m/24m) - 7K non-myeloid brain cells
4. GSE179358 (6m/12m/18m/24m) - 29K pure hippocampal microglia

Target ages: 3m, 6m, 12m, 18m, 21m, 24m
"""

import anndata as ad, scanpy as sc, numpy as np, pandas as pd
import os, gc, sys, time, warnings
warnings.filterwarnings('ignore')
import logging; logging.getLogger('harmonypy').setLevel(logging.WARNING)

sys.path.insert(0, '../26_CMAP')
import torch, CMAP_py

BASE = ''
CMAP_DIR = f'{BASE}models/cmap_output/Brain'
MG_DIR = f'{BASE}data/microglia_aging'
os.makedirs(CMAP_DIR, exist_ok=True)
np.random.seed(42)

# ── Load all data sources ──
print("Loading data sources...")

# 1. PanSci Brain
pansci = ad.read_h5ad(f'{BASE}data/pansci/per_organ_subset/Brain_WT_subset.h5ad')
pansci.obs['cellType'] = pansci.obs['Main_cell_type'].values
pansci.obs['age_std'] = pansci.obs['Age_group_standardized'].map({
    '03_months': '3m', '06_months': '6m', '21_months': '21m'})
pansci.obs['source'] = 'PanSci'
print(f"  PanSci: {pansci.shape}, ages={sorted(pansci.obs['age_std'].unique())}")

# 2. TMS Brain_Myeloid
myeloid = ad.read_h5ad(f'{BASE}data/tms/per_organ/TMS_Brain_Myeloid_facs.h5ad')
myeloid.obs['cellType'] = myeloid.obs['cell_ontology_class'].values
myeloid.obs['age_std'] = myeloid.obs['age'].values
myeloid.obs['source'] = 'TMS_Myeloid'
print(f"  TMS Myeloid: {myeloid.shape}, ages={sorted(myeloid.obs['age_std'].unique())}")

# 3. TMS Brain_NonMyeloid
nonmyeloid = ad.read_h5ad(f'{BASE}data/tms/per_organ/TMS_Brain_NonMyeloid_facs.h5ad')
nonmyeloid.obs['cellType'] = nonmyeloid.obs['cell_ontology_class'].values
nonmyeloid.obs['age_std'] = nonmyeloid.obs['age'].values
nonmyeloid.obs['source'] = 'TMS_NonMyeloid'
print(f"  TMS NonMyeloid: {nonmyeloid.shape}, ages={sorted(nonmyeloid.obs['age_std'].unique())}")

# 4. GSE179358 hippocampal microglia
print("  Loading GSE179358...")
mg_adatas = {}
for age, prefix in [('6m','GSM5416205_6mth'), ('12m','GSM5416206_12mth'),
                     ('18m','GSM5416207_18mth'), ('24m','GSM5416208_24mth')]:
    a = sc.read_mtx(f'{MG_DIR}/{prefix}_matrix.mtx.gz').T
    gn = pd.read_csv(f'{MG_DIR}/{prefix}_genes.tsv.gz', header=None, sep='\t')
    bc = pd.read_csv(f'{MG_DIR}/{prefix}_barcodes.tsv.gz', header=None, sep='\t')
    a.var_names = gn[1].values if gn.shape[1] > 1 else gn[0].values
    a.obs_names = [f'GSE179358_{age}_{x}' for x in bc[0].values]
    a.var_names_make_unique()
    a.obs['cellType'] = 'microglia_hippocampal'
    a.obs['age_std'] = age
    a.obs['source'] = 'GSE179358'
    mg_adatas[age] = a
    print(f"    {age}: {a.shape[0]} cells")

# ── Load ST ──
print("Loading ST Brain...")
st = ad.read_h5ad(f'{BASE}data/st/per_organ/Brain_CTRL1.h5ad')
if st.shape[0] > 15000:
    st = st[np.random.choice(st.shape[0], 15000, replace=False)].copy()
st_genes = set(st.var_names)

# ── Find shared genes across ALL sources ──
all_var = [set(pansci.var_names), set(myeloid.var_names), set(nonmyeloid.var_names), st_genes]
for a in mg_adatas.values():
    all_var.append(set(a.var_names))
shared = sorted(set.intersection(*all_var))
print(f"Shared genes across all sources + ST: {len(shared)}")

# Subset all to shared genes
pansci = pansci[:, shared].copy()
myeloid = myeloid[:, shared].copy()
nonmyeloid = nonmyeloid[:, shared].copy()
for k in mg_adatas:
    mg_adatas[k] = mg_adatas[k][:, shared].copy()

# ── Process each age ──
all_ages = ['3m', '6m', '12m', '18m', '21m', '24m']

for age in all_ages:
    print(f"\n{'='*50}")
    print(f"Age: {age}")
    print(f"{'='*50}")

    parts = []

    # PanSci non-microglia (provides diverse brain cell types)
    p_mask = pansci.obs['age_std'] == age
    if p_mask.sum() > 0:
        p_sub = pansci[p_mask]
        non_mg = p_sub[p_sub.obs['cellType'] != 'Microglia']
        if non_mg.shape[0] > 5000:
            non_mg = non_mg[np.random.choice(non_mg.shape[0], 5000, replace=False)]
        parts.append(non_mg.copy())
        print(f"  PanSci non-microglia: {non_mg.shape[0]}")

    # TMS NonMyeloid
    n_mask = nonmyeloid.obs['age_std'] == age
    if n_mask.sum() > 0:
        n_sub = nonmyeloid[n_mask]
        if n_sub.shape[0] > 2000:
            n_sub = n_sub[np.random.choice(n_sub.shape[0], 2000, replace=False)]
        parts.append(n_sub.copy())
        print(f"  TMS NonMyeloid: {n_sub.shape[0]}")

    # TMS Brain_Myeloid
    m_mask = myeloid.obs['age_std'] == age
    if m_mask.sum() > 0:
        m_sub = myeloid[m_mask]
        if m_sub.shape[0] > 3000:
            m_sub = m_sub[np.random.choice(m_sub.shape[0], 3000, replace=False)]
        parts.append(m_sub.copy())
        print(f"  TMS Myeloid: {m_sub.shape[0]}")

    # GSE179358 hippocampal microglia
    if age in mg_adatas:
        mg_sub = mg_adatas[age]
        if mg_sub.shape[0] > 3000:
            mg_sub = mg_sub[np.random.choice(mg_sub.shape[0], 3000, replace=False)]
        parts.append(mg_sub.copy())
        print(f"  GSE179358 microglia: {mg_sub.shape[0]}")

    # PanSci microglia (small but add them)
    if p_mask.sum() > 0:
        p_mg = pansci[p_mask & (pansci.obs['cellType'] == 'Microglia')]
        if p_mg.shape[0] > 0:
            parts.append(p_mg.copy())
            print(f"  PanSci microglia: {p_mg.shape[0]}")

    if not parts:
        print(f"  [SKIP] No data")
        continue

    merged = ad.concat(parts)
    merged.obs_names_make_unique()
    merged.obs['annotation'] = merged.obs['cellType'].values
    print(f"  Total: {merged.shape[0]} cells, {merged.obs['cellType'].nunique()} types")

    mg_count = merged.obs['cellType'].str.contains('icroglia|acrophage', case=False).sum()
    print(f"  Microglia+macrophage: {mg_count} ({mg_count/len(merged)*100:.1f}%)")

    # Run CMAP
    shared_cmap = sorted(set(merged.var_names) & set(st.var_names))
    print(f"  CMAP: {merged.shape[0]} cells, {st.shape[0]} spots, {len(shared_cmap)} genes")

    t0 = time.time()
    CMAP_py.run_cmap_fast(merged, st, shared_cmap,
        spatial_data_type='square', cluster_key='annotation', svm_threshold=0.5,
        num_epochs=2000, n_near_spot=5, dis_cut=None, radius=1/6,
        device='cuda' if torch.cuda.is_available() else 'cpu')
    print(f"  CMAP done in {time.time()-t0:.0f}s")

    coords = merged.obsm['spatial_CMAP']
    mapped = np.isfinite(coords[:, 0])
    print(f"  Mapped: {mapped.sum()}/{len(merged)} ({mapped.mean():.1%})")

    if mapped.sum() < 30:
        print(f"  [SKIP] Too few")
        continue

    a = merged[mapped].copy()
    a.obsm['spatial'] = a.obsm['spatial_CMAP']
    a.layers['count'] = a.X.copy()
    out_path = f'{CMAP_DIR}/Brain_{age}.h5ad'
    a.write_h5ad(out_path)
    print(f"  Saved: {out_path} ({a.shape[0]} cells)")

    del merged, a, parts; gc.collect()

print("\n=== Brain microglia-enriched CMAP complete ===")
print("Now run gsMap on Brain_3m/6m/12m/18m/21m/24m")
