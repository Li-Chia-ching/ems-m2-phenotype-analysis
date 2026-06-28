---

# EMS-induced M<sub>2</sub> Phenotyping Analysis Pipeline

## Overview

This repository contains a reproducible computational pipeline for the analysis of EMS-induced mutagenesis in *Medicago polymorpha* M<sub>2</sub> populations. The workflow integrates family-level mutation burden estimation, phenotype–family association analysis, and global phenotypic spectrum characterization.

The pipeline was developed for downstream publication and Supplementary Table/Figure generation.

---

## Data Input

All analyses are based on a curated M<sub>2</sub> phenotypic dataset:

* `EMS_M2_Cleaned_1425plants.csv`
* **Population size:** 1,425 M<sub>2</sub> plants
* **Structured by:**
* Family ID (*n* = 98 families)
* Multi-timepoint phenotype scoring (day 30–120)
* Standardized categorical phenotype annotations



> **Note:** Raw data are not included in this repository.

---

## Analytical Workflow

### 1. Family-Level Mutation Burden Analysis (Supplementary Table S1)

Mutation burden was defined as the presence of at least one scored phenotype across four developmental timepoints.

A global mutation probability ($P_{\text{global}}$) was estimated from the full population. Family-level enrichment of mutant frequency was tested using a one-sided binomial exact test:

$$H_0: p_{\text{family}} = P_{\text{global}}$$

$$H_1: p_{\text{family}} > P_{\text{global}}$$

Multiple testing correction was performed using the Benjamini–Hochberg false discovery rate (FDR).

**Outputs:**

* Family-level mutation burden
* Expected vs. observed mutant counts
* P-values and FDR-adjusted q-values

### 2. Phenotype–Family Association Analysis (Supplementary Table S2)

To quantify associations between genetic families and discrete phenotypic classes, a contingency table (Family × Phenotype Category) was constructed.

Adjusted standardized residuals (Haberman residuals) were calculated as:

$$z_{ij} = \frac{O_{ij} - E_{ij}}{\sqrt{E_{ij} \left(1 - \frac{R_i}{N}\right) \left(1 - \frac{C_j}{N}\right)}}$$

Where:

* $O_{ij}$: observed count
* $E_{ij}$: expected count under independence
* $R_i$: row marginal total
* $C_j$: column marginal total
* $N$: total sample size

Two-sided P-values were derived from the standard normal approximation. FDR correction (Benjamini–Hochberg) was applied across all cells.

**Outputs:**

* Observed/expected counts
* Standardized residuals (z-scores)
* Raw and FDR-adjusted significance values

### 3. Global Phenotypic Spectrum Analysis

Phenotypic categories were aggregated across all individuals to compute the global frequency ($GF$):

$$GF_j = \frac{n_j}{N_{\text{total}}}$$

Where $n_j$ is the number of plants exhibiting phenotype $j$. This analysis provides a population-level overview of EMS-induced phenotypic diversity.

### 4. Trait-Level Quantitative Analysis and Visualization

Continuous trait responses to EMS concentration were modeled using:

* **Log-logistic dose–response (LL.4)** for survival traits
* **Quadratic regression** for mutation rate optimization
* **Pearson correlation** for trait association structure
* **Faceted trend analysis** across developmental traits

All model fits were accompanied by bootstrap resampling (n = 1000) to estimate confidence intervals for key parameters (e.g., LD50, optimal EMS concentration).

---

## Statistical Design Justification

Statistical analyses were designed to account for the hierarchical structure of EMS-induced mutagenesis data and the non-independence inherent in family-based segregation.

* **Mutation Burden:** At the population level, mutation burden was evaluated using a binomial framework under a global mutation probability estimated from the full M2 cohort, ensuring a consistent null expectation across genetically structured families.
* **Phenotype–Family Association:** Contingency analysis was performed using adjusted standardized residuals (Haberman correction) to control for marginal effects of heterogeneous family sizes and unequal phenotype prevalence. This approach provides a variance-stabilized measure of deviation from independence in sparse count matrices.
* **Multiple Testing:** To mitigate inflation of Type I error arising from multiple testing across families and phenotype categories, all p-values were corrected using the Benjamini–Hochberg false discovery rate procedure.
* **Dose–Response Modeling:** Continuous trait responses to EMS dosage were modeled using parametric (LL.4 log-logistic) and non-linear quadratic regression frameworks selected based on biological interpretability and empirical goodness-of-fit, with parameter uncertainty quantified via non-parametric bootstrap resampling (1,000 iterations).

Collectively, this analytical design prioritizes robustness to unbalanced sampling, distributional heterogeneity, and multiple inference while maintaining interpretability for genetic and physiological inference in mutagenized plant populations.

---

## Software and Dependencies

### Python

* `pandas`
* `numpy`
* `scipy`
* `statsmodels`

### R

* `ggplot2`
* `dplyr`
* `tidyr`
* `patchwork`
* `viridis`
* `scico`
* `stringr`

---

## Outputs

The pipeline generates the following assets:

### Supplementary Tables

* **Supplementary Table S1:** Family mutation burden analysis
* **Supplementary Table S2:** Family × phenotype association matrix

### Figures

* **Fig. 1:** Dose–response and mutation optimization
* **Fig. 2:** Global phenotypic spectrum
* **Fig. 3:** Family segregation and trait correlation
* **Fig. 4:** Trait response dynamics

> *All figures are exported in publication-quality PDF format (300 DPI equivalent vector output).*

---

## Reproducibility

To ensure complete reproducibility of the results:

* Analyses were performed using fixed random seeds where applicable.
* Bootstrapping iterations were set to *n* = 1000.
* FDR correction utilized the standard Benjamini–Hochberg method.
* Missing values were handled using pairwise deletion or exclusion depending on specific model assumptions.

---

## Notes

This repository is intended as a companion analytical framework for EMS mutagenesis phenotyping studies and can be readily adapted to other plant genetic populations with categorical trait scoring systems.
