#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate Supplementary Table S2.csv
Family × Phenotype Category contingency table analysis with Haberman adjusted
standardized residuals (z-scores) and BH-FDR correction.

Data source: EMS_M2_Cleaned_1425plants.csv
Output: Supplementary Table S2.csv
"""

import pandas as pd
import numpy as np
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests
import os

# ── 1. Load data ───────────────────────────────────────────────────────────
input_path = r"C:\Users\...\03_Outputs_and_Figures\Cleaned_Data\EMS_M2_Cleaned_1425plants.csv"
output_path = r"C:\Users\...\03_Outputs_and_Figures\Supplementary Table S2.csv"

df = pd.read_csv(input_path)
print(f"Loaded: {input_path}")
print(f"Dimensions: {df.shape[0]} rows x {df.shape[1]} cols")
print()

# ── 2. Extract phenotype categories ────────────────────────────────────────
# Phenotype columns (AllCategories versions contain semicolon-separated standardized categories)
pheno_cols = [
    "Phenotype_day30_AllCategories",
    "Phenotype_day60_AllCategories",
    "Phenotype_day90_AllCategories",
    "Phenotype_day120_AllCategories",
]

# Build a long-format dataframe: one row per (Family, Plant, Phenotype_Category)
records = []
for _, row in df.iterrows():
    fam = row["Family"]
    plant_id = row["Plant_ID"]
    # Collect all unique phenotype categories across the 4 time points for this plant
    cats_for_plant = set()
    for col in pheno_cols:
        val = row[col]
        if pd.notna(val) and str(val).strip() != "":
            # Split semicolon-separated categories
            parts = str(val).split(";")
            for p in parts:
                p = p.strip()
                if p != "":
                    cats_for_plant.add(p)
    # Emit one record per unique (Family, Plant, Category)
    for cat in cats_for_plant:
        records.append({"Family": fam, "Plant_ID": plant_id, "Phenotype_Category": cat})

long_df = pd.DataFrame(records)
print(f"Long-format records (Family × Plant × Category): {len(long_df)}")
print(f"Unique families: {long_df['Family'].nunique()}")
print(f"Unique phenotype categories: {long_df['Phenotype_Category'].nunique()}")
print()

# ── 3. Construct Observed count matrix O_ij ────────────────────────────────
# Rows = Family (98), Columns = Phenotype Category
# Count: number of PLANTS in family i that exhibit category j (binary per plant)
obs_matrix = long_df.groupby(["Family", "Phenotype_Category"]).size().unstack(fill_value=0)

# Ensure all 98 families are present (some may have zero records if no phenotypes)
all_families = sorted(df["Family"].unique())
obs_matrix = obs_matrix.reindex(index=all_families, fill_value=0)

print(f"Contingency table dimensions: {obs_matrix.shape[0]} families x {obs_matrix.shape[1]} categories")
print()

# ── 4. Calculate Adjusted Standardized Residuals (Haberman z) ─────────────
N = obs_matrix.values.sum()  # grand total
R = obs_matrix.sum(axis=1).values  # row totals (per family)
C = obs_matrix.sum(axis=0).values  # column totals (per category)

print(f"Grand total N = {N}")
print(f"Row totals range: {R.min()} - {R.max()}")
print(f"Column totals range: {C.min()} - {C.max()}")
print()

# Expected counts E_ij = (R_i * C_j) / N
R_mat = R[:, np.newaxis]  # column vector
C_mat = C[np.newaxis, :]  # row vector
E = (R_mat * C_mat) / N

# Adjusted standardized residuals (Haberman)
# z_ij = (O_ij - E_ij) / sqrt(E_ij * (1 - R_i/N) * (1 - C_j/N))
denom = np.sqrt(E * (1 - R_mat / N) * (1 - C_mat / N))
# Avoid division by zero
denom = np.where(denom == 0, np.nan, denom)
z = (obs_matrix.values - E) / denom

# ── 5. Raw P-values (two-tailed) ───────────────────────────────────────────
# P_ij = 2 * (1 - Phi(|z_ij|))
abs_z = np.abs(z)
p_raw = 2 * (1 - norm.cdf(abs_z))
# Handle NaN (set p-value to 1.0 where z is NaN)
p_raw = np.where(np.isnan(p_raw), 1.0, p_raw)

# ── 6. BH-FDR Correction ───────────────────────────────────────────────────
# Flatten all p-values, apply BH, reshape back
p_flat = p_raw.flatten()
_, q_flat, _, _ = multipletests(p_flat, method="fdr_bh", alpha=0.05)
q = q_flat.reshape(p_raw.shape)

# ── 7. Build output dataframe ──────────────────────────────────────────────
families = obs_matrix.index.tolist()
categories = obs_matrix.columns.tolist()

rows_out = []
for i, fam in enumerate(families):
    for j, cat in enumerate(categories):
        o_val = int(obs_matrix.values[i, j])
        e_val = float(E[i, j])
        z_val = float(z[i, j]) if not np.isnan(z[i, j]) else float("nan")
        p_val = float(p_raw[i, j])
        q_val = float(q[i, j])
        rows_out.append({
            "Family": fam,
            "Phenotype_Category": cat,
            "Observed_Count": o_val,
            "Expected_Count": round(e_val, 4),
            "Standardized_Residual_z": round(z_val, 4) if not np.isnan(z_val) else float("nan"),
            "Raw_P_value": p_val,
            "FDR_q_value": q_val,
        })

result_df = pd.DataFrame(rows_out)

# Sort: significant enrichment hotspots first (z > 0 and q < 0.05), then by q_value
result_df["is_hotspot"] = (result_df["Standardized_Residual_z"] > 0) & (result_df["FDR_q_value"] < 0.05)
result_df = result_df.sort_values(["is_hotspot", "FDR_q_value"], ascending=[False, True]).drop(columns=["is_hotspot"]).reset_index(drop=True)

# Save
result_df.to_csv(output_path, index=False)
print(f"Saved to: {output_path}")
print(f"Total cells in matrix: {len(result_df)}")
print()

# ── 8. Aggregate Summary ───────────────────────────────────────────────────
hotspots = result_df[(result_df["Standardized_Residual_z"] > 0) & (result_df["FDR_q_value"] < 0.05)]
print("=" * 70)
print("AGGREGATE SUMMARY")
print("=" * 70)
print(f"Total enrichment hotspots (z > 0, q < 0.05): {len(hotspots)}")
print(f"Matrix dimensions: {len(families)} families x {len(categories)} categories = {len(families)*len(categories)} cells")
print()

# Top 15 hotspots
print("TOP 15 ENRICHMENT HOTSPOTS:")
print(hotspots.head(15).to_string(index=False))
print()

# ── 9. Family 80 Targeted Extraction ──────────────────────────────────────
print("=" * 70)
print("TARGETED EXTRACTION: FAMILY 80 — CHLOROSIS/PALLOR")
print("=" * 70)
# Find chlorosis-related categories
chlorosis_keywords = ["chlorosis", "pallor", "albinism", "chim"]
fam80_rows = result_df[result_df["Family"] == 80]
chlorosis_rows = fam80_rows[fam80_rows["Phenotype_Category"].str.lower().str.contains("|".join(chlorosis_keywords), na=False)]

if len(chlorosis_rows) > 0:
    print(f"Found {len(chlorosis_rows)} chlorosis-related row(s) for Family 80:")
    print(chlorosis_rows.to_string(index=False))
else:
    print("No chlorosis-related rows found for Family 80.")
    print("\nAll Family 80 rows with non-zero observed counts:")
    fam80_nonzero = fam80_rows[fam80_rows["Observed_Count"] > 0]
    print(fam80_nonzero.to_string(index=False))

# Also check the specific category name
print()
print("All unique phenotype categories containing 'chlor' or 'pallor' or 'albin':")
all_cats = sorted(result_df["Phenotype_Category"].unique())
chlor_cats = [c for c in all_cats if any(k in c.lower() for k in ["chlor", "pallor", "albin"])]
for c in chlor_cats:
    print(f"  - {c}")
