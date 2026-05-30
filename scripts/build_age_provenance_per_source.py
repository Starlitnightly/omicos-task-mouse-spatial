"""Build per-(organ, age, source) cell count table — not just one primary source per sample.

For each of the 126 samples:
  - If obs.source exists (Brain), count cells per source value.
  - Otherwise use the inferred primary source for the whole sample.
"""
import os
import scanpy as sc
import pandas as pd

LDSC_ROOT = 'models/gsmap_age_output'
CMAP_DIR  = 'models/cmap_output'

def find_cmap(sample):
    for o in os.listdir(CMAP_DIR):
        p = f'{CMAP_DIR}/{o}/{sample}.h5ad'
        if os.path.exists(p):
            return o, p
    return None, None

def infer_primary(sample, organ, obs_cols, obs_head):
    if 'Sub_cell_type' in obs_cols and 'Lineage' in obs_cols:
        return 'PanSci'
    if 'PCR_sample_name' in obs_cols or 'Ligation_barcode' in obs_cols:
        return 'EasySci'
    if 'method' in obs_cols and 'mouse.id' in obs_cols:
        return 'TMS'
    if sample == 'Brain_21m':
        return 'TMS'
    age_label = sample.replace(f'{organ}_', '')
    age_m = None
    if age_label.endswith('m'):
        try: age_m = int(age_label[:-1])
        except: pass
    elif '_months' in age_label:
        try: age_m = int(age_label.split('_')[0])
        except: pass
    if organ in ('Ileum','Stomach','Spleen','Thymus','Colon') and age_m in (6,12,18,21,24):
        return 'MCA2'
    return 'unknown'

# Map fine-grained TMS labels to single 'TMS' bucket; MCA2 / PanSci / EasySci / Hammond stay as-is
def normalize_src(s):
    s = str(s)
    if s.startswith('TMS'):    return 'TMS'
    if s == 'GSE179358':       return 'Hammond_microglia'
    return s

samples = sorted([s for s in os.listdir(LDSC_ROOT)
                  if os.path.exists(f'{LDSC_ROOT}/{s}/spatial_ldsc/{s}_AD.csv.gz')])

rows = []
for s in samples:
    organ, path = find_cmap(s)
    if not path: continue
    a = sc.read_h5ad(path, backed='r')
    age_label = s.replace(f'{organ}_', '')
    age_m = None
    if age_label.endswith('m'):
        try: age_m = int(age_label[:-1])
        except: pass
    elif '_months' in age_label:
        try: age_m = int(age_label.split('_')[0])
        except: pass

    # If obs.source is present, count per source — otherwise use primary inferred
    if 'source' in a.obs.columns:
        vc = a.obs['source'].astype(str).value_counts()
        for src_raw, n in vc.items():
            rows.append({'sample': s, 'organ': organ, 'age_label': age_label,
                         'age_m': age_m, 'source': normalize_src(src_raw),
                         'n_cells': int(n)})
    else:
        primary = infer_primary(s, organ, list(a.obs.columns), a.obs.head(50))
        rows.append({'sample': s, 'organ': organ, 'age_label': age_label,
                     'age_m': age_m, 'source': primary, 'n_cells': a.shape[0]})

df = pd.DataFrame(rows)
# Collapse rows that map to the same source within a sample (e.g. TMS_Myeloid + TMS_NonMyeloid -> TMS)
df = (df.groupby(['sample','organ','age_label','age_m','source'], observed=True)['n_cells']
        .sum().reset_index())

OUT = 'results/age_sample_provenance_per_source.csv'
df.to_csv(OUT, index=False)
print(f'wrote {OUT}')
print(f'\nTotal rows: {len(df)}  (was 126 single-primary; now multi-source per sample)')
print(f'Total cells: {df["n_cells"].sum():,}')
print(f'\nPer source:')
print(df.groupby('source')['n_cells'].sum().sort_values(ascending=False))

# How many samples have multi-source breakdown?
mc = df.groupby('sample')['source'].nunique()
multi = (mc > 1).sum()
print(f'\nSamples with >1 source: {multi} / {df["sample"].nunique()}')
print(df[df["sample"].isin(mc[mc>1].index)].sort_values(['sample','source']).head(20).to_string())
