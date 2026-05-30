#!/usr/bin/env python
"""Process remaining ★ gaps: TMS Droplet 21m + other missing ages.
Serial CMAP to prevent OOM."""

import anndata as ad, numpy as np, pandas as pd, scanpy as sc
import os, gc, sys, time, warnings
warnings.filterwarnings('ignore')
import logging; logging.getLogger('harmonypy').setLevel(logging.WARNING)
sys.path.insert(0, '../26_CMAP')
import torch, CMAP_py

BASE = ''
TMS = f'{BASE}data/tms/per_organ'
CMAP_DIR = f'{BASE}models/cmap_output'
ST_DIR = f'{BASE}data/st/per_organ'
np.random.seed(42)

st_organ_map = {
    'Heart':'Heart','Kidney':'Kidney','Lung':'Lung','Liver':'Liver',
    'BAT':'Brown_Fat','Muscle':'Muscle','Colon':'Colon','Skin':'Skin',
    'Pancreas':'Pancreas',
}

def run_cmap_and_save(adata_sc, organ, age_str, ct_col='cell_ontology_class'):
    """Run CMAP on single-cell data and save per-age h5ad"""
    out_dir = f'{CMAP_DIR}/{organ}'
    out_path = f'{out_dir}/{organ}_{age_str}.h5ad'
    os.makedirs(out_dir, exist_ok=True)

    if os.path.exists(out_path):
        print(f"    [EXISTS] {out_path}")
        return

    if adata_sc.shape[0] < 30:
        print(f"    [SKIP] Too few cells: {adata_sc.shape[0]}")
        return

    st_name = st_organ_map.get(organ, organ)
    st = ad.read_h5ad(f'{ST_DIR}/{st_name}_CTRL1.h5ad')
    if st.shape[0] > 15000:
        st = st[np.random.choice(st.shape[0], 15000, replace=False)].copy()

    # Prepare SC data
    adata_sc.obs['annotation'] = adata_sc.obs[ct_col].astype(str).values
    adata_sc.obs['cellType'] = adata_sc.obs[ct_col].astype(str).values
    shared = sorted(set(adata_sc.var_names) & set(st.var_names))

    if len(shared) < 3000:
        print(f"    [SKIP] Too few shared genes: {len(shared)}")
        return

    print(f"    CMAP: {adata_sc.shape[0]} cells, {st.shape[0]} spots, {len(shared)} genes")

    t0 = time.time()
    CMAP_py.run_cmap_fast(adata_sc, st, shared,
        spatial_data_type='square', cluster_key='annotation', svm_threshold=0.5,
        num_epochs=2000, n_near_spot=5, dis_cut=None, radius=1/6,
        device='cuda' if torch.cuda.is_available() else 'cpu')
    print(f"    CMAP done in {time.time()-t0:.0f}s")

    coords = adata_sc.obsm['spatial_CMAP']
    mapped = np.isfinite(coords[:, 0])
    print(f"    Mapped: {mapped.sum()}/{len(adata_sc)} ({mapped.mean():.1%})")

    if mapped.sum() < 30:
        print(f"    [SKIP] Too few mapped")
        del st; gc.collect()
        return

    a = adata_sc[mapped].copy()
    a.obsm['spatial'] = a.obsm['spatial_CMAP']
    a.layers['count'] = a.X.copy()
    a.write_h5ad(out_path)
    print(f"    Saved: {out_path} ({a.shape[0]} cells)")
    del st, a; gc.collect()


# ═══ Process TMS Droplet 21m gaps ═══
tms_droplet_gaps = [
    # (TMS file, organ, target ages)
    ('TMS_Fat_droplet.h5ad', 'BAT', ['21m', '30m']),
    ('TMS_Heart_Aorta_droplet.h5ad', 'Heart', ['21m']),  # already have from FACS but droplet has more cells
    ('TMS_Kidney_droplet.h5ad', 'Kidney', ['21m']),
    ('TMS_Liver_droplet.h5ad', 'Liver', ['1m', '21m', '30m']),
    ('TMS_Lung_droplet.h5ad', 'Lung', ['21m']),
    ('TMS_Limb_Muscle_droplet.h5ad', 'Muscle', ['1m', '21m', '30m']),
    ('TMS_Skin_droplet.h5ad', 'Skin', ['18m', '21m']),
    ('TMS_Pancreas_facs.h5ad', 'Pancreas', ['3m']),
    ('TMS_Large_Intestine_droplet.h5ad', 'Colon', ['30m']),
]

print("=" * 60)
print("Processing TMS Droplet/FACS gaps")
print("=" * 60)

for fname, organ, target_ages in tms_droplet_gaps:
    path = f'{TMS}/{fname}'
    if not os.path.exists(path):
        print(f"\n[SKIP] {fname} not found")
        continue

    a = ad.read_h5ad(path)
    print(f"\n--- {fname} -> {organ} ---")
    print(f"  Total: {a.shape[0]}, ages: {sorted(a.obs['age'].unique())}")

    for age in target_ages:
        if age not in a.obs['age'].values:
            print(f"  [SKIP] {age} not in dataset")
            continue

        sub = a[a.obs['age'] == age].copy()
        if sub.shape[0] > 10000:
            sub = sub[np.random.choice(sub.shape[0], 10000, replace=False)].copy()

        print(f"\n  {organ}_{age}: {sub.shape[0]} cells, {sub.obs['cell_ontology_class'].nunique()} types")
        run_cmap_and_save(sub, organ, age)
        del sub; gc.collect()

    del a; gc.collect()

print("\n=== All TMS Droplet gaps processed ===")
