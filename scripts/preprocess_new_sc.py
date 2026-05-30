#!/usr/bin/env python
"""
Preprocess new single-cell datasets to fill age gaps:
1. TMS FACS/Droplet: Heart, Kidney, Lung, Liver, BAT, Muscle, Colon, Skin, Brain (18m/24m etc.)
2. MCA2 (GSE153562): Stomach, Intestine, Spleen, Thymus (12m/18m/24m)
3. GSE210669: Intestine aging (24m)

For each: load -> filter ages -> subsample 10K/age -> subset to shared genes with ST -> save
"""

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import os, gc, sys, glob, warnings
warnings.filterwarnings('ignore')

BASE = ''
ST_DIR = f'{BASE}data/st/per_organ'
TMS_DIR = f'{BASE}data/tms/per_organ'
MCA2_DIR = f'{BASE}data/mca2'
OUT_DIR = f'{BASE}data/pansci/per_organ_subset'  # reuse same output dir
CMAP_DIR = f'{BASE}models/cmap_output'

np.random.seed(42)

# ST organ name mapping for gene intersection
st_organ_map = {
    'Heart': 'Heart', 'Kidney': 'Kidney', 'Lung': 'Lung', 'Liver': 'Liver',
    'BAT': 'Brown_Fat', 'Muscle': 'Muscle', 'Colon': 'Colon',
    'Skin': 'Skin', 'Brain': 'Brain', 'Stomach': 'Stomach',
    'Ileum': 'Small_Intestine', 'Spleen': 'Spleen', 'Thymus': 'Thymus',
    'Bone_Marrow': 'Bone_Marrow', 'Pancreas': 'Pancreas',
}

def get_st_genes(organ):
    st_name = st_organ_map.get(organ, organ)
    st = ad.read_h5ad(f'{ST_DIR}/{st_name}_CTRL1.h5ad', backed='r')
    return set(st.var_names)

def subsample_per_age(adata, age_col, max_per_age=10000):
    subsets = []
    for age in sorted(adata.obs[age_col].unique()):
        idx = adata.obs_names[adata.obs[age_col] == age]
        n = min(max_per_age, len(idx))
        sampled = np.random.choice(idx, n, replace=False)
        subsets.extend(sampled)
    return adata[subsets].copy()

# ═══════════════════════════════════════════════
# PART 1: Process TMS FACS/Droplet new organs
# ═══════════════════════════════════════════════

tms_files = {
    # organ: (file, type, target_organ_for_ST)
    'Heart': ('TMS_Heart_Aorta_droplet.h5ad', 'droplet', 'Heart'),
    'Kidney': ('TMS_Kidney_droplet.h5ad', 'droplet', 'Kidney'),
    'Lung': ('TMS_Lung_droplet.h5ad', 'droplet', 'Lung'),
    'Liver': ('TMS_Liver_facs.h5ad', 'facs', 'Liver'),
    'BAT': ('TMS_BAT_facs.h5ad', 'facs', 'BAT'),
    'Muscle': ('TMS_Limb_Muscle_facs.h5ad', 'facs', 'Muscle'),
    'Colon': ('TMS_Large_Intestine_facs.h5ad', 'facs', 'Colon'),
    'Skin': ('TMS_Skin_facs.h5ad', 'facs', 'Skin'),
    'Brain_Myeloid': ('TMS_Brain_Myeloid_facs.h5ad', 'facs', 'Brain'),
    'Brain_NonMyeloid': ('TMS_Brain_NonMyeloid_facs.h5ad', 'facs', 'Brain'),
}

# Ages we already have per organ (to avoid duplicate processing)
existing_ages = {
    'Heart': {'03_months','06_months','12_months','16_months','23_months'},
    'Kidney': {'03_months','06_months','12_months','16_months','23_months'},
    'Lung': {'03_months','06_months','12_months','16_months','23_months'},
    'Liver': {'03_months','06_months','12_months','16_months','23_months'},
    'BAT': {'03_months','06_months','12_months','16_months','23_months'},
    'Muscle': {'03_months','06_months','12_months','16_months','23_months'},
    'Colon': {'03_months','06_months','12_months','16_months','23_months'},
    'Skin': {'21m'},
    'Brain': {'03_months','06_months','21_months'},
}

# TMS age -> standardized month string
tms_age_to_std = {'1m': '1m', '3m': '3m', '18m': '18m', '21m': '21m', '24m': '24m', '30m': '30m'}

def process_tms():
    print("\n" + "="*60)
    print("PART 1: Processing TMS FACS/Droplet data")
    print("="*60)

    for organ_label, (fname, dtype, st_organ) in tms_files.items():
        path = f'{TMS_DIR}/{fname}'
        if not os.path.exists(path):
            print(f"\n[SKIP] {organ_label}: {fname} not found")
            continue

        print(f"\n--- {organ_label} ({fname}) ---")
        adata = ad.read_h5ad(path)
        print(f"  Loaded: {adata.shape}")

        # Get ages
        ages = sorted(adata.obs['age'].unique())
        print(f"  Ages: {ages}")

        # Filter to ages we don't already have
        new_ages = [a for a in ages if a not in existing_ages.get(st_organ, set())]
        if not new_ages:
            print(f"  [SKIP] All ages already covered")
            continue
        print(f"  New ages to process: {new_ages}")

        adata = adata[adata.obs['age'].isin(new_ages)].copy()

        # Get shared genes with ST
        st_genes = get_st_genes(st_organ)
        shared = sorted(set(adata.var_names) & st_genes)
        print(f"  Shared genes with ST: {len(shared)}")

        if len(shared) < 5000:
            print(f"  [WARN] Too few shared genes, skipping")
            continue

        adata = adata[:, shared].copy()

        # Subsample
        adata = subsample_per_age(adata, 'age', max_per_age=10000)

        # Add annotation column
        adata.obs['annotation'] = adata.obs['cell_ontology_class'].values
        adata.obs['Age_group'] = adata.obs['age'].values
        adata.obs['cellType'] = adata.obs['cell_ontology_class'].values

        # For Brain, merge Myeloid and NonMyeloid later
        if 'Brain' in organ_label:
            out_name = f'Brain_TMS_{organ_label.split("_")[1]}_subset.h5ad'
        else:
            out_name = f'{st_organ}_TMS_subset.h5ad'

        out_path = f'{OUT_DIR}/{out_name}'
        adata.write_h5ad(out_path)
        print(f"  Saved: {out_path} ({adata.shape[0]} cells, {os.path.getsize(out_path)/1e6:.0f}MB)")

        # Now run CMAP + save per-age h5ad for each new age
        process_cmap_ages(adata, st_organ, new_ages)

        del adata; gc.collect()


def process_cmap_ages(adata_sc, st_organ, ages):
    """Run CMAP mapping and save per-age h5ad files"""
    import sys
    sys.path.insert(0, '../26_CMAP')
    import torch, CMAP_py
    import logging; logging.getLogger('harmonypy').setLevel(logging.WARNING)

    st_name = st_organ_map.get(st_organ, st_organ)
    st_path = f'{ST_DIR}/{st_name}_CTRL1.h5ad'
    adata_st = ad.read_h5ad(st_path)

    # Subsample ST if too large
    if adata_st.shape[0] > 25000:
        np.random.seed(42)
        adata_st = adata_st[np.random.choice(adata_st.shape[0], 25000, replace=False)].copy()

    shared = sorted(set(adata_sc.var_names) & set(adata_st.var_names))
    print(f"  CMAP: SC={adata_sc.shape[0]}, ST={adata_st.shape[0]}, shared={len(shared)}")

    # Run CMAP
    import time
    t0 = time.time()
    CMAP_py.run_cmap_fast(adata_sc, adata_st, shared,
        spatial_data_type='square', cluster_key='annotation', svm_threshold=0.5,
        num_epochs=2000, n_near_spot=5, dis_cut=None, radius=1/6,
        device='cuda' if torch.cuda.is_available() else 'cpu')
    print(f"  CMAP done in {time.time()-t0:.0f}s")

    coords = adata_sc.obsm['spatial_CMAP']
    mapped = np.isfinite(coords[:, 0])
    print(f"  Mapped: {mapped.sum()}/{len(mapped)} ({mapped.mean():.1%})")

    # Save per-age h5ad
    out_dir = f'{CMAP_DIR}/{st_organ}'
    os.makedirs(out_dir, exist_ok=True)

    age_col = 'Age_group' if 'Age_group' in adata_sc.obs.columns else 'age'
    for age in sorted(adata_sc.obs.loc[mapped, age_col].unique()):
        am = mapped & (adata_sc.obs[age_col] == age).values
        if am.sum() < 30:
            print(f"  [SKIP] {age}: only {am.sum()} cells mapped")
            continue
        a = adata_sc[am].copy()
        a.obsm['spatial'] = a.obsm['spatial_CMAP']
        a.layers['count'] = a.X.copy()
        a.obs['annotation'] = a.obs['cellType']
        ac = str(age).replace(' ', '_').replace('/', '_')
        out_path = f'{out_dir}/{st_organ}_{ac}.h5ad'

        # Skip if already exists
        if os.path.exists(out_path):
            print(f"  [EXISTS] {out_path}")
            continue

        a.write_h5ad(out_path)
        print(f"  Saved {st_organ}_{ac}: {a.shape[0]} cells")

    del adata_st; gc.collect()


# ═══════════════════════════════════════════════
# PART 2: Process MCA2 (GSE153562)
# ═══════════════════════════════════════════════

def process_mca2():
    print("\n" + "="*60)
    print("PART 2: Processing MCA2 (GSE153562)")
    print("="*60)

    tar_path = f'{MCA2_DIR}/GSE153562_RAW.tar'
    if not os.path.exists(tar_path):
        print(f"[SKIP] {tar_path} not found")
        return

    # Extract tar
    extract_dir = f'{MCA2_DIR}/extracted'
    if not os.path.exists(extract_dir):
        os.makedirs(extract_dir, exist_ok=True)
        import tarfile
        print("  Extracting tar...")
        with tarfile.open(tar_path, 'r') as tar:
            tar.extractall(extract_dir)
        print("  Done")

    # List extracted files
    files = sorted(glob.glob(f'{extract_dir}/*.txt.gz') + glob.glob(f'{extract_dir}/*.csv.gz'))
    print(f"  Found {len(files)} files")
    for f in files[:10]:
        print(f"    {os.path.basename(f)}")
    if len(files) > 10:
        print(f"    ... and {len(files)-10} more")

    # MCA2 uses Microwell-seq, files are typically named by sample
    # We need to identify which files correspond to which organ/age
    # Parse sample names from filenames
    organ_age_files = {}
    for f in files:
        bn = os.path.basename(f).replace('.txt.gz', '').replace('.csv.gz', '')
        # Try to parse organ and age from filename
        # Typical: GSMxxxxxx_Organ_Age.txt.gz
        print(f"  Parsing: {bn}")

    # If files are just GSM IDs, we need the sample metadata
    # Let's check the first file format
    if files:
        import gzip
        with gzip.open(files[0], 'rt') as fh:
            header = fh.readline().strip()
            line1 = fh.readline().strip()
            print(f"  First file header: {header[:100]}")
            print(f"  First file line1: {line1[:100]}")

    # We'll need to map GSM IDs to organs/ages
    # For now, print what we have and continue


# ═══════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════

if __name__ == '__main__':
    mode = sys.argv[1] if len(sys.argv) > 1 else 'tms'

    if mode == 'tms':
        process_tms()
    elif mode == 'mca2':
        process_mca2()
    elif mode == 'all':
        process_tms()
        process_mca2()
    else:
        print(f"Usage: {sys.argv[0]} [tms|mca2|all]")
