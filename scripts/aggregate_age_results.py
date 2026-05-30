#!/usr/bin/env python
"""Aggregate gsMap age-resolved results across all organs and age groups."""

import os
import glob
import pandas as pd
import numpy as np

BASE = os.environ.get('BASE', '.')
AGE_DIR = f'{BASE}/models/gsmap_age_output'
OUT_DIR = f'{BASE}/results'
os.makedirs(OUT_DIR, exist_ok=True)

# ── 1. Aggregate Cauchy combination p-values ──
print("=== Aggregating Cauchy combination p-values ===")
cauchy_files = sorted(glob.glob(f'{AGE_DIR}/*/cauchy_combination/*_AD.Cauchy.csv.gz'))
print(f"Found {len(cauchy_files)} Cauchy result files")

all_cauchy = []
for f in cauchy_files:
    sample = os.path.basename(os.path.dirname(os.path.dirname(f)))
    df = pd.read_csv(f)
    # Parse organ and age from sample name
    # Formats: Organ_03_months, Organ_3m, etc.
    parts = sample.rsplit('_', 2)
    if parts[-1] == 'months':
        organ = '_'.join(parts[:-2])
        age = parts[-2] + '_months'
    elif parts[-1].endswith('m'):
        organ = '_'.join(parts[:-1])
        age = parts[-1]
    else:
        organ = sample
        age = 'unknown'
    df['sample'] = sample
    df['organ'] = organ
    df['age'] = age
    all_cauchy.append(df)

if all_cauchy:
    cauchy_df = pd.concat(all_cauchy, ignore_index=True)
    cauchy_df.to_csv(f'{OUT_DIR}/age_all_cauchy_AD.csv.gz', index=False, compression='gzip')
    print(f"Saved: {OUT_DIR}/age_all_cauchy_AD.csv.gz")
    print(f"  Organs: {cauchy_df['organ'].nunique()}")
    print(f"  Age groups: {sorted(cauchy_df['age'].unique())}")
    print(f"  Total annotation entries: {len(cauchy_df)}")
else:
    print("No Cauchy files found!")

# ── 2. Aggregate per-spot p-values ──
print("\n=== Aggregating per-spot p-values ===")
spot_files = sorted(glob.glob(f'{AGE_DIR}/*/spatial_ldsc/*_AD.csv.gz'))
print(f"Found {len(spot_files)} spatial LDSC result files")

all_spots = []
for f in spot_files:
    sample = os.path.basename(os.path.dirname(os.path.dirname(f)))
    df = pd.read_csv(f)
    parts = sample.rsplit('_', 2)
    if parts[-1] == 'months':
        organ = '_'.join(parts[:-2])
        age = parts[-2] + '_months'
    elif parts[-1].endswith('m'):
        organ = '_'.join(parts[:-1])
        age = parts[-1]
    else:
        organ = sample
        age = 'unknown'
    df['sample'] = sample
    df['organ'] = organ
    df['age'] = age
    all_spots.append(df)

if all_spots:
    spots_df = pd.concat(all_spots, ignore_index=True)
    spots_df.to_csv(f'{OUT_DIR}/age_all_spots_AD_pvalues.csv.gz', index=False, compression='gzip')
    print(f"Saved: {OUT_DIR}/age_all_spots_AD_pvalues.csv.gz")
    print(f"  Total spots: {len(spots_df)}")
else:
    print("No spot-level files found!")

# ── 3. Summary table: organ x age Cauchy p-values ──
if all_cauchy:
    print("\n=== Summary: min Cauchy p-value per organ x age ===")
    # Use the annotation-level Cauchy p-value (min across annotations per organ)
    p_col = [c for c in cauchy_df.columns if 'p' in c.lower() and 'cauchy' in c.lower()]
    if not p_col:
        p_col = [c for c in cauchy_df.columns if 'p' in c.lower()]
    if p_col:
        pcol = p_col[0]
        summary = cauchy_df.groupby(['organ', 'age'])[pcol].min().unstack()
        summary.to_csv(f'{OUT_DIR}/age_organ_summary.csv')
        print(f"Saved: {OUT_DIR}/age_organ_summary.csv")
        print(summary.to_string())

print("\n=== Done ===")
