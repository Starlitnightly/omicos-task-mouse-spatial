#!/usr/bin/env python
"""Process MCA2 (GSE153562) data for missing organ-age combinations.
Loads dge.txt.gz files, merges replicates, runs CMAP, saves per-age h5ad."""

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import os, gc, sys, glob, gzip, time, warnings
warnings.filterwarnings('ignore')
import logging; logging.getLogger('harmonypy').setLevel(logging.WARNING)

BASE = os.environ.get('BASE', '')
MCA2_DIR = f'{BASE}/data/mca2/extracted'
ST_DIR = f'{BASE}/data/st/per_organ'
CMAP_DIR = f'{BASE}/models/cmap_output'

sys.path.insert(0, '../26_CMAP')
import torch, CMAP_py

np.random.seed(42)

st_organ_map = {
    'Stomach': 'Stomach', 'Intestine': 'Small_Intestine',
    'Spleen': 'Spleen', 'Thymus': 'Thymus',
}
# Map to our organ naming
organ_out_map = {
    'Stomach': 'Stomach', 'Intestine': 'Ileum',
    'Spleen': 'Spleen', 'Thymus': 'Thymus',
}
age_map = {
    'OneYear': '12m', 'EighteenMonths': '18m', 'TwoYears': '24m',
    'Adult': '6w',  # skip adult (6-8 weeks, we already have PanSci 3m)
}

# Ages we already have - skip these
existing_cmap = {}
for organ_out in ['Stomach', 'Ileum', 'Spleen', 'Thymus']:
    d = f'{CMAP_DIR}/{organ_out}'
    if os.path.exists(d):
        existing = [f.replace(f'{organ_out}_','').replace('.h5ad','')
                    for f in os.listdir(d) if f.endswith('.h5ad') and 'cmap_mapped' not in f]
        existing_cmap[organ_out] = set(existing)
    else:
        existing_cmap[organ_out] = set()

print("Existing CMAP ages:")
for k,v in existing_cmap.items():
    print(f"  {k}: {sorted(v)}")

def load_dge(filepath):
    """Load a single dge.txt.gz file as anndata (cells x genes)"""
    df = pd.read_csv(filepath, sep='\t', index_col=0, compression='gzip')
    # df is genes x cells, transpose
    adata = ad.AnnData(X=df.T.values.astype(np.float32),
                       obs=pd.DataFrame(index=df.columns),
                       var=pd.DataFrame(index=df.index))
    return adata

def process_organ(organ_name):
    """Process all ages for one MCA2 organ"""
    organ_out = organ_out_map[organ_name]
    st_name = st_organ_map[organ_name]

    # Find all files for this organ
    files = sorted(glob.glob(f'{MCA2_DIR}/*{organ_name}*_dge.txt.gz'))
    print(f"\n{'='*50}")
    print(f"Processing {organ_name} -> {organ_out}")
    print(f"Found {len(files)} files")

    # Group by age
    age_files = {}
    for f in files:
        bn = os.path.basename(f)
        # Parse age from filename
        for age_key, age_std in age_map.items():
            if age_key in bn:
                if age_std not in age_files:
                    age_files[age_std] = []
                age_files[age_std].append(f)
                break

    print(f"Age groups: {', '.join(f'{k}({len(v)})' for k,v in sorted(age_files.items()))}")

    # Load ST for gene intersection
    st_path = f'{ST_DIR}/{st_name}_CTRL1.h5ad'
    adata_st = ad.read_h5ad(st_path)
    if adata_st.shape[0] > 15000:
        np.random.seed(42)
        adata_st = adata_st[np.random.choice(adata_st.shape[0], 15000, replace=False)].copy()
    st_genes = set(adata_st.var_names)

    out_dir = f'{CMAP_DIR}/{organ_out}'
    os.makedirs(out_dir, exist_ok=True)

    for age_std, flist in sorted(age_files.items()):
        if age_std == '6w':
            print(f"  [SKIP] {age_std} (already have young data from PanSci)")
            continue

        if age_std in existing_cmap.get(organ_out, set()):
            print(f"  [EXISTS] {organ_out}_{age_std}.h5ad already exists")
            continue

        print(f"\n  --- {organ_out} {age_std} ({len(flist)} replicates) ---")

        # Load and merge replicates
        adatas = []
        for f in flist:
            a = load_dge(f)
            # Make obs_names unique per file
            prefix = os.path.basename(f).split('_dge')[0]
            a.obs_names = [f'{prefix}_{x}' for x in a.obs_names]
            adatas.append(a)
            print(f"    Loaded {os.path.basename(f)}: {a.shape}")

        adata_sc = ad.concat(adatas)
        adata_sc.var_names_make_unique()
        print(f"    Merged: {adata_sc.shape}")

        # Filter genes
        shared = sorted(set(adata_sc.var_names) & st_genes)
        print(f"    Shared genes with ST: {len(shared)}")
        if len(shared) < 3000:
            print(f"    [SKIP] Too few shared genes")
            continue

        adata_sc = adata_sc[:, shared].copy()

        # Subsample if too many cells
        if adata_sc.shape[0] > 10000:
            idx = np.random.choice(adata_sc.shape[0], 10000, replace=False)
            adata_sc = adata_sc[idx].copy()

        # MCA2 has no cell type annotations - run basic clustering
        print(f"    Running basic clustering...")
        sc.pp.normalize_total(adata_sc, target_sum=1e4)
        sc.pp.log1p(adata_sc)
        adata_sc.layers['lognorm'] = adata_sc.X.copy()

        # Store raw counts for CMAP
        # Reload raw counts
        adatas_raw = []
        for f in flist:
            a = load_dge(f)
            prefix = os.path.basename(f).split('_dge')[0]
            a.obs_names = [f'{prefix}_{x}' for x in a.obs_names]
            adatas_raw.append(a)
        adata_raw = ad.concat(adatas_raw)
        adata_raw.var_names_make_unique()
        adata_raw = adata_raw[:, shared].copy()
        if adata_raw.shape[0] > 10000:
            adata_raw = adata_raw[adata_sc.obs_names].copy()

        # Use raw counts for CMAP
        adata_sc.X = adata_raw.X.copy()
        del adata_raw, adatas_raw; gc.collect()

        # Quick clustering for annotation
        sc.pp.highly_variable_genes(adata_sc, n_top_genes=2000, flavor='seurat_v3',
                                     layer='lognorm' if 'lognorm' in adata_sc.layers else None)
        sc.tl.pca(adata_sc, n_comps=30, use_highly_variable=True)
        sc.pp.neighbors(adata_sc, n_pcs=30)
        sc.tl.leiden(adata_sc, resolution=0.5)
        adata_sc.obs['annotation'] = adata_sc.obs['leiden'].astype(str)
        adata_sc.obs['cellType'] = adata_sc.obs['leiden'].astype(str)
        adata_sc.obs['Age_group'] = age_std
        del adata_sc.layers['lognorm']

        print(f"    Clusters: {adata_sc.obs['annotation'].nunique()}")
        print(f"    Cells: {adata_sc.shape[0]}")

        # Run CMAP
        shared_cmap = sorted(set(adata_sc.var_names) & set(adata_st.var_names))
        print(f"    CMAP: SC={adata_sc.shape[0]}, ST={adata_st.shape[0]}, shared={len(shared_cmap)}")

        t0 = time.time()
        try:
            CMAP_py.run_cmap_fast(adata_sc, adata_st, shared_cmap,
                spatial_data_type='square', cluster_key='annotation', svm_threshold=0.5,
                num_epochs=2000, n_near_spot=5, dis_cut=None, radius=1/6,
                device='cuda' if torch.cuda.is_available() else 'cpu')
            print(f"    CMAP done in {time.time()-t0:.0f}s")
        except Exception as e:
            print(f"    CMAP FAILED: {e}")
            del adata_sc; gc.collect()
            continue

        coords = adata_sc.obsm['spatial_CMAP']
        mapped = np.isfinite(coords[:, 0])
        print(f"    Mapped: {mapped.sum()}/{len(mapped)} ({mapped.mean():.1%})")

        if mapped.sum() < 30:
            print(f"    [SKIP] Too few mapped cells")
            del adata_sc; gc.collect()
            continue

        # Save per-age h5ad
        a = adata_sc[mapped].copy()
        a.obsm['spatial'] = a.obsm['spatial_CMAP']
        a.layers['count'] = a.X.copy()
        out_path = f'{out_dir}/{organ_out}_{age_std}.h5ad'
        a.write_h5ad(out_path)
        print(f"    Saved: {out_path} ({a.shape[0]} cells)")

        del adata_sc, a; gc.collect()

    del adata_st; gc.collect()


# Main
if __name__ == '__main__':
    for organ in ['Stomach', 'Intestine', 'Spleen', 'Thymus']:
        process_organ(organ)

    print("\n=== MCA2 processing complete ===")
