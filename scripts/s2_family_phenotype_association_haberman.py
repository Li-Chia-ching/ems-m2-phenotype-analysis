#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate Supplementary Table S2.csv
Family × Phenotype Category contingency table analysis with Haberman adjusted
standardized residuals (z-scores) and BH-FDR correction.
"""

import pandas as pd
import numpy as np
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests
import os

# ── 1. 使用相对路径 ───────────────────────────────────────────────────────
INPUT_DIR = os.path.join(".", "03_Outputs_and_Figures", "Cleaned_Data")
OUTPUT_DIR = os.path.join(".", "03_Outputs_and_Figures")

input_path = os.path.join(INPUT_DIR, "EMS_M2_Cleaned_1425plants.csv")
output_path = os.path.join(OUTPUT_DIR, "Supplementary_Table_S2.csv")

df = pd.read_csv(input_path)
print(f"Loaded: {input_path}")
print(f"Dimensions: {df.shape[0]} rows x {df.shape[1]} cols\n")

# ── 2. Extract phenotype categories ────────────────────────────────────────
pheno_cols = [
    "Phenotype_day30_AllCategories", "Phenotype_day60_AllCategories",
    "Phenotype_day90_AllCategories", "Phenotype_day120_AllCategories",
]

records = []
for _, row in df.iterrows():
    fam = row["Family"]
    plant_id = row["Plant_ID"]
    cats_for_plant = set()
    for col in pheno_cols:
        val = row[col]
        if pd.notna(val) and str(val).strip() != "":
            parts = str(val).split(";")
            for p in parts:
                p = p.strip()
                if p != "":
                    cats_for_plant.add(p)
    for cat in cats_for_plant:
        records.append({"Family": fam, "Plant_ID": plant_id, "Phenotype_Category": cat})

long_df = pd.DataFrame(records)
print(f"Long-format records: {len(long_df)}")
print(f"Unique families: {long_df['Family'].nunique()}\n")

# ── 3. Construct Observed count matrix O_ij ────────────────────────────────
obs_matrix = long_df.groupby(["Family", "Phenotype_Category"]).size().unstack(fill_value=0)
all_families = sorted(df["Family"].unique())
obs_matrix = obs_matrix.reindex(index=all_families, fill_value=0)

# ── 4. Calculate Adjusted Standardized Residuals (Haberman z) ─────────────
N = obs_matrix.values.sum()
R = obs_matrix.sum(axis=1).values
C = obs_matrix.sum(axis=0).values

R_mat = R[:, np.newaxis]
C_mat = C[np.newaxis, :]
E = (R_mat * C_mat) / N

denom = np.sqrt(E * (1 - R_mat / N) * (1 - C_mat / N))
denom = np.where(denom == 0, np.nan, denom)
z = (obs_matrix.values - E) / denom

# ── 5. Raw P-values & 6. BH-FDR Correction ────────────────────────────────
abs_z = np.abs(z)
p_raw = 2 * (1 - norm.cdf(abs_z))
p_raw = np.where(np.isnan(p_raw), 1.0, p_raw)

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
            "Family": fam, "Phenotype_Category": cat, "Observed_Count": o_val,
            "Expected_Count": round(e_val, 4),
            "Standardized_Residual_z": round(z_val, 4) if not np.isnan(z_val) else float("nan"),
            "Raw_P_value": p_val, "FDR_q_value": q_val,
        })

result_df = pd.DataFrame(rows_out)
result_df["is_hotspot"] = (result_df["Standardized_Residual_z"] > 0) & (result_df["FDR_q_value"] < 0.05)
result_df = result_df.sort_values(["is_hotspot", "FDR_q_value"], ascending=[False, True]).drop(columns=["is_hotspot"]).reset_index(drop=True)

os.makedirs(OUTPUT_DIR, exist_ok=True)
result_df.to_csv(output_path, index=False)
print(f"Saved to: {output_path}\n")

# ── 8. Aggregate Summary ───────────────────────────────────────────────────
hotspots = result_df[(result_df["Standardized_Residual_z"] > 0) & (result_df["FDR_q_value"] < 0.05)]
print("TOP 15 ENRICHMENT HOTSPOTS:")
print(hotspots.head(15).to_string(index=False), "\n")
