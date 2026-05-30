"""Build curated provenance table: 126 (organ x age) age samples -> source dataset / paper.

Sources identified by:
  - obs columns + first-row metadata (PanSci/EasySci/TMS markers)
  - explicit `source` column when present (Brain)
  - preprocess_new_sc.py logic (MCA2 covered Stomach/Intestine/Spleen/Thymus 12-24m)
  - memory note (Brain_21m = TMS microglia 18+24m proxy)
"""
import os, json
import scanpy as sc
import pandas as pd

BASE = os.environ.get('BASE', '../')
DATA_BASE = os.environ.get('DATA_BASE', '')

LDSC_ROOT = os.path.join(DATA_BASE, 'models/gsmap_age_output')
CMAP_DIR  = os.path.join(DATA_BASE, 'models/cmap_output')

samples = sorted([s for s in os.listdir(LDSC_ROOT)
                  if os.path.exists(f'{LDSC_ROOT}/{s}/spatial_ldsc/{s}_AD.csv.gz')])

def find_cmap(sample):
    for o in os.listdir(CMAP_DIR):
        p = f'{CMAP_DIR}/{o}/{sample}.h5ad'
        if os.path.exists(p):
            return o, p
    return None, None

def primary_source(sample, organ, obs_cols, obs_head):
    # Explicit source column (Brain mostly)
    if 'source' in obs_cols:
        srcs = sorted(set(obs_head['source'].dropna().astype(str).unique()))
        # Map mixed sources to a "Brain atlas" label
        # PanSci/TMS_NonMyeloid/TMS_Myeloid/MCA2/GSE179358
        if any('TMS' in s for s in srcs):
            return 'TMS'
        if 'PanSci' in srcs:
            return 'PanSci'
        if 'MCA2' in srcs:
            return 'MCA2'
        if 'GSE179358' in srcs:
            return 'Hammond_microglia'
        return srcs[0]

    # PanSci signature
    if 'Sub_cell_type' in obs_cols and 'Lineage' in obs_cols:
        return 'PanSci'
    # EasySci (PanSci brain)
    if 'PCR_sample_name' in obs_cols or 'Ligation_barcode' in obs_cols:
        return 'EasySci'
    # TMS signature
    if 'method' in obs_cols and 'mouse.id' in obs_cols:
        m = str(obs_head['method'].dropna().iloc[0]) if 'method' in obs_head.columns else ''
        return 'TMS'  # collapse facs/droplet under TMS for the figure

    # Brain_21m special case (memorised: microglia proxy from TMS 18m + 24m)
    if sample == 'Brain_21m':
        return 'TMS'

    # Hand-curate unknowns based on preprocess_new_sc.py:
    # MCA2 covered Stomach / Ileum (intestine) / Spleen / Thymus at 12-24m
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

rows = []
for s in samples:
    organ, path = find_cmap(s)
    if not path:
        rows.append({'sample': s, 'organ': '?', 'source': 'NO_CMAP', 'age_label': '?', 'age_m': None, 'n_cells': 0})
        continue
    a = sc.read_h5ad(path, backed='r')
    src = primary_source(s, organ, list(a.obs.columns), a.obs.head(50))
    age_label = s.replace(f'{organ}_', '')
    age_m = None
    if age_label.endswith('m'):
        try: age_m = int(age_label[:-1])
        except: pass
    elif '_months' in age_label:
        try: age_m = int(age_label.split('_')[0])
        except: pass
    rows.append({'sample': s, 'organ': organ, 'source': src,
                 'age_label': age_label, 'age_m': age_m, 'n_cells': a.shape[0]})

df = pd.DataFrame(rows)
df.to_csv(os.path.join(DATA_BASE, 'results/age_sample_provenance.csv'), index=False)
print(f'Total: {len(df)}')
print(f'\nSource counts:')
print(df['source'].value_counts())
print(f'\nOrgans:')
print(df['organ'].value_counts())
print(f'\nAges (months):')
print(df['age_m'].value_counts().sort_index())
print(f'\nCross-tab: organ × source (cells)')
ct = df.groupby(['organ','source'], observed=True)['n_cells'].sum().unstack(fill_value=0)
print(ct)
