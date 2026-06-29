#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate Supplementary Table S1.csv
Family-level mutation burden analysis with Binomial Exact Test + BH-FDR correction.
"""

import pandas as pd
import numpy as np
from scipy.stats import binomtest
from statsmodels.stats.multitest import multipletests
import os

# ── 1. 使用相对路径 (消除隐私风险并提高可移植性) ──────────────────
# 假设运行目录为项目根目录
INPUT_DIR = os.path.join(".", "03_Outputs_and_Figures", "Cleaned_Data")
OUTPUT_DIR = os.path.join(".", "03_Outputs_and_Figures")

input_path = os.path.join(INPUT_DIR, "EMS_M2_Cleaned_1425plants.csv")
output_path = os.path.join(OUTPUT_DIR, "Supplementary_Table_S1.csv")

df = pd.read_csv(input_path)
print(f"Loaded: {input_path}")
print(f"Dimensions: {df.shape[0]} rows x {df.shape[1]} cols")
print(f"Columns: {df.columns.tolist()}\n")

# ── 2. Define 'Mutant' ────────────────────────────────────────────────────
pheno_cols = [c for c in df.columns if c.startswith("Phenotype_day")]
print(f"Phenotype columns detected: {pheno_cols}\n")

def is_mutant(row):
    for col in pheno_cols:
        val = row[col]
        if pd.notna(val) and str(val).strip() != "":
            return True
    return False

df["is_mutant"] = df.apply(is_mutant, axis=1)

# ── 3. Calculate Global Background ────────────────────────────────────────
N_total = len(df)
K_total = df["is_mutant"].sum()
P_global = K_total / N_total

print("=" * 60)
print("GLOBAL BACKGROUND")
print("=" * 60)
print(f"Total population (N): {N_total}")
print(f"Total mutant plants (K): {K_total}")
print(f"Global mutation rate (P_global): {P_global:.6f} ({P_global*100:.2f}%)\n")

# ── 4. Family-Level Analysis ──────────────────────────────────────────────
family_stats = df.groupby("Family").agg(
    Total_Plants=("is_mutant", "count"),
    Mutant_Count=("is_mutant", "sum")
).reset_index()

family_stats["Expected_Mutant_Count"] = (family_stats["Total_Plants"] * P_global).round(2)

# ── 5. Binomial Exact Test ────────────────────────────────────────────────
p_values = []
for _, row in family_stats.iterrows():
    k = int(row["Mutant_Count"])
    n = int(row["Total_Plants"])
    result = binomtest(k, n, p=P_global, alternative="greater")
    p_values.append(result.pvalue)

family_stats["P_value"] = p_values

# ── 6. BH-FDR Correction ──────────────────────────────────────────────────
_, q_values, _, _ = multipletests(family_stats["P_value"].values, method="fdr_bh", alpha=0.05)
family_stats["q_value"] = q_values

# ── 7. Format and sort ────────────────────────────────────────────────────
family_stats_display = family_stats.copy()
family_stats_display["P_value"] = family_stats_display["P_value"].apply(lambda x: f"{x:.6e}" if x < 0.001 else f"{x:.6f}")
family_stats_display["q_value"] = family_stats_display["q_value"].apply(lambda x: f"{x:.6e}" if x < 0.001 else f"{x:.6f}")

family_stats = family_stats.sort_values("q_value", ascending=True).reset_index(drop=True)
family_stats_display = family_stats_display.sort_values("q_value", ascending=True).reset_index(drop=True)

output_cols = ["Family", "Total_Plants", "Mutant_Count", "Expected_Mutant_Count", "P_value", "q_value"]
family_stats_out = family_stats[output_cols]

# ── 8. Save to CSV ────────────────────────────────────────────────────────
os.makedirs(OUTPUT_DIR, exist_ok=True)
family_stats_out.to_csv(output_path, index=False)
print(f"Saved to: {output_path}")
print(f"Total families: {len(family_stats_out)}\n")

# ── 9. Display Summary ────────────────────────────────────────────────────
print("=" * 80)
print("FIRST 10 ROWS (sorted by q_value ascending)")
print("=" * 80)
print(family_stats_display[output_cols].head(10).to_string(index=False), "\n")

print("=" * 60)
print("SUMMARY")
print("=" * 60)
print(f"Families with q < 0.05: {(family_stats['q_value'] < 0.05).sum()}")
print(f"Families with q < 0.01: {(family_stats['q_value'] < 0.01).sum()}")
print(f"Families with q < 0.001: {(family_stats['q_value'] < 0.001).sum()}")
