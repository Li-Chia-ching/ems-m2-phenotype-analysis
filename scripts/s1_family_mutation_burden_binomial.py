#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate Supplementary Table S1.csv
Family-level mutation burden analysis with Binomial Exact Test + BH-FDR correction.

Data source: EMS_M2_Cleaned_1425plants.csv
Output: Supplementary Table S1.csv
"""

import pandas as pd
import numpy as np
from scipy.stats import binomtest
from statsmodels.stats.multitest import multipletests

# ── 1. Load data ──────────────────────────────────────────────────────────
input_path = r"C:\Users\...\03_Outputs_and_Figures\Cleaned_Data\EMS_M2_Cleaned_1425plants.csv"
output_path = r"C:\Users\...\03_Outputs_and_Figures\Supplementary Table S1.csv"

df = pd.read_csv(input_path)
print(f"Loaded: {input_path}")
print(f"Dimensions: {df.shape[0]} rows x {df.shape[1]} cols")
print(f"Columns: {df.columns.tolist()}")
print()

# ── 2. Define 'Mutant' ────────────────────────────────────────────────────
# Identify phenotype columns (Phenotype_day30, day60, day90, day120 — both raw and AllCategories)
pheno_cols = [c for c in df.columns if c.startswith("Phenotype_day")]
print(f"Phenotype columns detected: {pheno_cols}")
print()

# A plant is 'Mutant' if ANY Phenotype_day* column has non-empty, non-NA content
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
print(f"Global mutation rate (P_global): {P_global:.6f} ({P_global*100:.2f}%)")
print()

# ── 4. Family-Level Analysis ──────────────────────────────────────────────
family_stats = df.groupby("Family").agg(
    Total_Plants=("is_mutant", "count"),
    Mutant_Count=("is_mutant", "sum")
).reset_index()

# Calculate expected mutants
family_stats["Expected_Mutant_Count"] = (family_stats["Total_Plants"] * P_global).round(2)

# ── 5. Binomial Exact Test (one-sided: observed >= expected) ──────────────
# H0: p = P_global; H1: p > P_global
p_values = []
for _, row in family_stats.iterrows():
    k = int(row["Mutant_Count"])
    n = int(row["Total_Plants"])
    result = binomtest(k, n, p=P_global, alternative="greater")
    p_values.append(result.pvalue)

family_stats["P_value"] = p_values

# ── 6. BH-FDR Correction ──────────────────────────────────────────────────
_, q_values, _, _ = multipletests(
    family_stats["P_value"].values,
    method="fdr_bh",
    alpha=0.05
)
family_stats["q_value"] = q_values

# ── 7. Format and sort ────────────────────────────────────────────────────
# Round P_value and q_value for display (keep full precision in file)
family_stats_display = family_stats.copy()
family_stats_display["P_value"] = family_stats_display["P_value"].apply(
    lambda x: f"{x:.6e}" if x < 0.001 else f"{x:.6f}"
)
family_stats_display["q_value"] = family_stats_display["q_value"].apply(
    lambda x: f"{x:.6e}" if x < 0.001 else f"{x:.6f}"
)

# Sort by q_value ascending
family_stats = family_stats.sort_values("q_value", ascending=True).reset_index(drop=True)
family_stats_display = family_stats_display.sort_values("q_value", ascending=True).reset_index(drop=True)

# Reorder columns for output
output_cols = ["Family", "Total_Plants", "Mutant_Count", "Expected_Mutant_Count", "P_value", "q_value"]
family_stats_out = family_stats[output_cols]

# ── 8. Save to CSV ────────────────────────────────────────────────────────
family_stats_out.to_csv(output_path, index=False)
print(f"Saved to: {output_path}")
print(f"Total families: {len(family_stats_out)}")
print()

# ── 9. Display first 10 rows ──────────────────────────────────────────────
print("=" * 80)
print("FIRST 10 ROWS (sorted by q_value ascending)")
print("=" * 80)
# Use display format for console
display_df = family_stats_display[output_cols].head(10)
print(display_df.to_string(index=False))
print()

# Summary statistics
print("=" * 60)
print("SUMMARY")
print("=" * 60)
sig_005 = (family_stats["q_value"] < 0.05).sum()
sig_001 = (family_stats["q_value"] < 0.01).sum()
sig_0001 = (family_stats["q_value"] < 0.001).sum()
print(f"Families with q < 0.05: {sig_005}")
print(f"Families with q < 0.01: {sig_001}")
print(f"Families with q < 0.001: {sig_0001}")
print(f"Family size range: {family_stats['Total_Plants'].min()}-{family_stats['Total_Plants'].max()}")
print(f"Family size median: {family_stats['Total_Plants'].median()}")
