# =============================================================================
# EMS Mutagenesis — Master Plotting Script (V1.0)
# Refactored for new manuscript storyline: 4 figures, no heatmap
#
# Figure Hierarchy (NEW):
#   Fig. 1: Dose-Response (p1a | p1b)
#   Fig. 2: Global Bar Chart (p2a alone — for Illustrator compositing)
#   Fig. 3: Segregation & Correlation (p2b | p3a)
#   Fig. 4: Trait Trends (p3b alone)
#
# Data Sources:
#   - 01_Raw_Data/Archived/EMS_M2_MP2603V2_Sheet1.csv (dose-response)
#   - 03_Outputs_and_Figures/Cleaned_Data/EMS_M2_Cleaned_1425plants.csv (M2 phenotypes)
#
# Aesthetics: scico + viridis + Okabe-Ito palettes, sans-serif fonts, 300 DPI
# =============================================================================

# ── 0. Package Setup ─────────────────────────────────────────────────────────
# Ensure user library is in path (for scico installed outside system library)
user_lib <- "C:/Users/.../R/win-library/4.6"
if (dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths()))

packages <- c("dplyr", "tidyr", "ggplot2", "readr", "broom", "forcats",
              "stringr", "scales", "patchwork", "viridis", "scico", "readxl", "tools")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
    library(pkg, character.only = TRUE)
  }
}
library(conflicted)
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

# ── 1. Global Configuration ──────────────────────────────────────────────────
base_font_size <- 9
geom_text_size <- base_font_size / .pt  # ggplot2 internal conversion

# Okabe-Ito colorblind-friendly palette
okabe_ito <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7", "#000000"
)

# Publication theme — minimalist, crisp axes, no panel border
pub_theme <- theme_bw(base_size = base_font_size, base_family = "sans") +
  theme(
    text = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray95", linewidth = 0.3),
    strip.background = element_rect(fill = "transparent", color = NA),
    strip.text = element_text(face = "bold", size = base_font_size),
    legend.position = "right",
    legend.key.size = unit(0.35, "cm"),
    legend.title = element_text(size = base_font_size - 1),
    legend.text = element_text(size = base_font_size - 2),
    plot.title = element_text(face = "bold", size = base_font_size + 1),
    plot.subtitle = element_text(size = base_font_size - 1, color = "gray30"),
    plot.margin = margin(8, 10, 8, 8)
  )

# Paths
project_root <- "C:/Users/..."
out_dir <- file.path(project_root, "03_Outputs_and_Figures")
cleaned_data_path <- file.path(out_dir, "Cleaned_Data", "EMS_M2_Cleaned_1425plants.csv")
sheet1_path <- file.path(project_root, "01_Raw_Data", "Archived", "EMS_M2_MP2603V2_Sheet1.csv")

# ── 2. Data Loading ──────────────────────────────────────────────────────────
# 2A. Dose-response data (Sheet1)
dose_df <- read_csv(sheet1_path, show_col_types = FALSE)
names(dose_df) <- tolower(gsub("\\s+", "_", trimws(names(dose_df))))

dose_data <- dose_df %>%
  select(ems_concentration, matches("_(mean|sd)$")) %>%
  filter(!is.na(ems_concentration)) %>%
  rename_with(~"survival_mean", matches("seedling_survival_rate_mean")) %>%
  rename_with(~"survival_sd", matches("seedling_survival_rate_sd"))

# 2B. M2 phenotype data (Cleaned)
m2_df <- read_csv(cleaned_data_path, show_col_types = FALSE)
N_total <- nrow(m2_df)  # 1425

# =============================================================================
# FIGURE 1: Dose-Response & M1 Survival (p1a | p1b)
# =============================================================================

# ── 1A: LL.4 Log-Logistic Model Fit ──────────────────────────────────────────
ld_data <- dose_data %>%
  filter(ems_concentration > 0, !is.na(survival_mean)) %>%
  mutate(ems_pct = ems_concentration * 100)

# Fit LL.4 using nls() (drc package not available for R 4.6.0)
nls_fit <- tryCatch({
  nls(survival_mean ~ c + (d - c) / (1 + exp(b * (log(ems_concentration) - log(e)))),
      data = ld_data,
      start = list(c = 0, d = 0.95, e = 0.006, b = 5),
      control = nls.control(maxiter = 500, warnOnly = TRUE))
}, error = function(e) NULL)

# Extract parameters and compute LD values
ld_params <- data.frame(Parameter = character(), Estimate = numeric(),
                        Lower = numeric(), Upper = numeric(), stringsAsFactors = FALSE)

if (!is.null(nls_fit)) {
  coefs <- coef(nls_fit)
  c_val <- coefs["c"]; d_val <- coefs["d"]; e_val <- coefs["e"]; b_val <- coefs["b"]

  inverse_LL4 <- function(p, c_val, d_val, e_val, b_val) {
    ratio <- (d_val - c_val) / (p - c_val) - 1
    if (ratio <= 0) return(NA)
    e_val * ratio^(1 / b_val)
  }

  ld50 <- e_val
  ld10 <- inverse_LL4(d_val - 0.10 * (d_val - c_val), c_val, d_val, e_val, b_val)
  ld90 <- inverse_LL4(d_val - 0.90 * (d_val - c_val), c_val, d_val, e_val, b_val)

  # Bootstrap CIs for LD50
  set.seed(2026)
  n_boot <- 1000
  boot_ld50 <- numeric(n_boot)
  boot_ld10 <- numeric(n_boot)
  boot_ld90 <- numeric(n_boot)

  for (i in seq_len(n_boot)) {
    boot_idx <- sample(seq_len(nrow(ld_data)), replace = TRUE)
    boot_data <- ld_data[boot_idx, ]
    boot_fit <- tryCatch({
      nls(survival_mean ~ c + (d - c) / (1 + exp(b * (log(ems_concentration) - log(e)))),
          data = boot_data,
          start = list(c = c_val, d = d_val, e = e_val, b = b_val),
          control = nls.control(maxiter = 200, warnOnly = TRUE))
    }, error = function(e) NULL)

    if (!is.null(boot_fit)) {
      bc <- coef(boot_fit)
      boot_ld50[i] <- bc["e"]
      boot_ld10[i] <- tryCatch(inverse_LL4(bc["d"] - 0.10 * (bc["d"] - bc["c"]), bc["c"], bc["d"], bc["e"], bc["b"]), error = function(e) NA)
      boot_ld90[i] <- tryCatch(inverse_LL4(bc["d"] - 0.90 * (bc["d"] - bc["c"]), bc["c"], bc["d"], bc["e"], bc["b"]), error = function(e) NA)
    } else {
      boot_ld50[i] <- NA; boot_ld10[i] <- NA; boot_ld90[i] <- NA
    }
  }

  ld_params <- data.frame(
    Parameter = c("LD10", "LD50", "LD90"),
    Estimate = c(ld10, ld50, ld90),
    Lower = c(quantile(boot_ld10, 0.025, na.rm = TRUE),
              quantile(boot_ld50, 0.025, na.rm = TRUE),
              quantile(boot_ld90, 0.025, na.rm = TRUE)),
    Upper = c(quantile(boot_ld10, 0.975, na.rm = TRUE),
              quantile(boot_ld50, 0.975, na.rm = TRUE),
              quantile(boot_ld90, 0.975, na.rm = TRUE))
  )
}

# Prediction data for smooth curve
pred_x <- seq(min(ld_data$ems_concentration), max(ld_data$ems_concentration), length.out = 200)
if (!is.null(nls_fit)) {
  pred_y <- predict(nls_fit, newdata = data.frame(ems_concentration = pred_x))
  pred_df <- data.frame(ems_concentration = pred_x, ems_pct = pred_x * 100,
                        Prediction = pred_y)
} else {
  pred_df <- data.frame(ems_concentration = pred_x, ems_pct = pred_x * 100,
                        Prediction = NA)
}

# Plot 1A — Okabe-Ito blue line, vermillion points
p1a <- ggplot() +
  geom_line(data = pred_df, aes(x = ems_pct, y = Prediction),
            color = okabe_ito[5], linewidth = 0.9) +
  geom_point(data = ld_data, aes(x = ems_pct, y = survival_mean),
             size = 2.5, color = okabe_ito[6], alpha = 0.9) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray60", linewidth = 0.4) +
  labs(
    title = "A  LL.4 Log-Logistic Dose-Response",
    x = "EMS Concentration (% v/v)",
    y = "M1 Seedling Survival Rate"
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +
  pub_theme +
  theme(legend.position = "none")

# Annotate LD50 if available
if (!is.null(nls_fit) && !is.na(ld50)) {
  p1a <- p1a +
    annotate("point", x = ld50 * 100, y = 0.5, color = okabe_ito[3], size = 2.5, shape = 18) +
    annotate("text", x = ld50 * 100, y = 0.58,
             label = sprintf("LD[50] == %.2f*\"%%\"", ld50 * 100),
             parse = TRUE, size = geom_text_size, color = okabe_ito[3], hjust = 0)
}

# ── 1B: Mutation Rate Quadratic Optimization ────────────────────────────────
mut_data <- dose_data %>%
  select(ems_concentration, m2_mutation_rate_mean) %>%
  filter(!is.na(m2_mutation_rate_mean)) %>%
  mutate(ems_pct = ems_concentration * 100)

model_mut <- lm(m2_mutation_rate_mean ~ ems_concentration + I(ems_concentration^2), data = mut_data)
coefs <- coef(model_mut)
a_coef <- coefs["I(ems_concentration^2)"]
b_coef <- coefs["ems_concentration"]
opt_conc <- if (!is.na(a_coef) && a_coef < 0) -b_coef / (2 * a_coef) else NA

# Bootstrap CI for optimal concentration
set.seed(2026)
boot_opt <- replicate(1000, {
  boot_idx <- sample(seq_len(nrow(mut_data)), replace = TRUE)
  boot_fit <- tryCatch(lm(m2_mutation_rate_mean ~ ems_concentration + I(ems_concentration^2),
                          data = mut_data[boot_idx, ]), error = function(e) NULL)
  if (!is.null(boot_fit)) {
    bc <- coef(boot_fit)
    if (!is.na(bc[3]) && bc[3] < 0) -bc[2] / (2 * bc[3]) else NA
  } else NA
})
opt_ci <- quantile(boot_opt, c(0.025, 0.975), na.rm = TRUE)

# Prediction data
pred_mut <- data.frame(ems_concentration = seq(min(mut_data$ems_concentration),
                                                max(mut_data$ems_concentration), length.out = 100))
pred_mut$ems_pct <- pred_mut$ems_concentration * 100
pred_mut$fitted <- predict(model_mut, newdata = pred_mut)
pred_mut$se <- predict(model_mut, newdata = pred_mut, se.fit = TRUE)$se.fit
pred_mut$lower <- pred_mut$fitted - 1.96 * pred_mut$se
pred_mut$upper <- pred_mut$fitted + 1.96 * pred_mut$se

# Plot 1B — muted low-alpha ribbon, Okabe-Ito line/points
p1b <- ggplot() +
  geom_ribbon(data = pred_mut, aes(x = ems_pct, ymin = lower, ymax = upper),
              fill = "gray75", alpha = 0.35) +
  geom_line(data = pred_mut, aes(x = ems_pct, y = fitted),
            color = okabe_ito[5], linewidth = 0.9) +
  geom_point(data = mut_data, aes(x = ems_pct, y = m2_mutation_rate_mean),
             size = 2.5, color = okabe_ito[6]) +
  labs(
    title = "B  Quadratic Mutation Rate Optimization",
    x = "EMS Concentration (% v/v)",
    y = expression(M[2]~mutation~rate)
  ) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +
  pub_theme

if (!is.na(opt_conc)) {
  p1b <- p1b +
    geom_vline(xintercept = opt_conc * 100, linetype = "dashed", color = okabe_ito[3], linewidth = 0.6) +
    annotate("text", x = opt_conc * 100, y = max(mut_data$m2_mutation_rate_mean, na.rm = TRUE) * 0.85,
             label = sprintf("C[opt] == %.2f*\"%%\"", opt_conc * 100),
             parse = TRUE, size = geom_text_size, color = okabe_ito[3], hjust = -0.1)
}

# Export plot data
write.csv(ld_data %>% mutate(ems_pct = ems_concentration * 100),
          file.path(out_dir, "Fig1A_PlotData.csv"), row.names = FALSE)
write.csv(mut_data, file.path(out_dir, "Fig1B_PlotData.csv"), row.names = FALSE)

message("=> Fig. 1 plots generated")

# =============================================================================
# FIGURE 2 (Part): Global Phenotypic Mutation Spectrum — Bar Chart (p2a)
# Exported alone for Adobe Illustrator compositing with physical photos
# =============================================================================

# Extract all phenotype categories from Day 90 (richest data)
global_cats <- m2_df %>%
  filter(!is.na(Phenotype_day90_AllCategories), Phenotype_day90_AllCategories != "") %>%
  mutate(cat_list = str_split(Phenotype_day90_AllCategories, ";")) %>%
  unnest(cat_list) %>%
  mutate(category = str_trim(cat_list)) %>%
  filter(category != "")

# Global frequencies
global_freq <- global_cats %>%
  count(category, name = "n_mutants") %>%
  mutate(
    GF = n_mutants / N_total,
    GF_pct = GF * 100
  ) %>%
  arrange(desc(GF)) %>%
  mutate(category = fct_reorder(category, GF))

# Clean category names for display
clean_names <- c(
  "Dwarfism_StuntedGrowth" = "Dwarfism",
  "Flowering_Early" = "Early Flowering",
  "Multifoliate" = "Multifoliate",
  "Chlorosis_Pallor" = "Chlorosis",
  "Branching_Clustered" = "Clustered Branch",
  "Leaf_Morphology_Small" = "Small Leaves",
  "Stem_Elongated" = "Elongated Stem",
  "Stem_Thin" = "Thin Stem",
  "Pigmentation_Anthocyanin" = "Anthocyanin",
  "Leaf_Morphology_Narrow" = "Narrow Leaves",
  "Chimerism" = "Chimerism",
  "Leaf_Morphology_Folded" = "Folded Leaves",
  "Senescence_Early" = "Early Senescence",
  "Lethal" = "Lethal",
  "Root_StrongSystem" = "Strong Root",
  "GrowthHabit_Prostrate" = "Prostrate",
  "GrowthHabit_Ascending" = "Ascending",
  "Chlorosis_Albinism" = "Albinism",
  "GrowthHabit_Compact" = "Compact Growth",
  "Branching_Reduced" = "Reduced Branch",
  "Disease_Susceptible" = "Disease Susc.",
  "Pigmentation_GreenPetiole" = "Green Petiole",
  "Leaf_Morphology_DeepVein" = "Deep Vein",
  "Leaf_Morphology_ShallowVein" = "Shallow Vein"
)

global_freq <- global_freq %>%
  mutate(category_label = clean_names[as.character(category)],
         category_label = ifelse(is.na(category_label), as.character(category), category_label),
         category_label = fct_reorder(category_label, GF))

# Family-level frequencies (top 10 families by mutant count)
family_sizes <- m2_df %>% count(Family, name = "family_N")

family_burden <- global_cats %>%
  group_by(Family) %>%
  summarise(total_mutants = n(), .groups = "drop") %>%
  left_join(family_sizes, by = "Family") %>%
  arrange(desc(total_mutants)) %>%
  head(10)

family_freq <- global_cats %>%
  filter(Family %in% family_burden$Family) %>%
  group_by(Family, category) %>%
  summarise(n = n(), .groups = "drop") %>%
  left_join(family_sizes, by = "Family") %>%
  mutate(
    FF = n / family_N,
    FF_pct = FF * 100,
    Family_label = paste0("Fam ", Family, " (n=", family_N, ")"),
    category_label = clean_names[as.character(category)],
    category_label = ifelse(is.na(category_label), as.character(category), category_label)
  ) %>%
  filter(category_label %in% levels(global_freq$category_label)) %>%
  mutate(category_label = factor(category_label, levels = levels(global_freq$category_label)))

# ── p2a: Global bar chart (Y = category_label, X = GF_pct) ──
# Title has NO panel letter — will be panel D in Adobe Illustrator
p2a <- ggplot(global_freq, aes(x = GF_pct, y = category_label)) +
  geom_col(fill = okabe_ito[5], width = 0.7, alpha = 0.9) +
  geom_text(aes(label = sprintf("%.2f%%", GF_pct)),
            hjust = -0.15, size = 2.2, family = "sans", color = "gray30") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.20))) +
  labs(
    title = "Global Phenotypic Mutation Spectrum",
    subtitle = expression(italic(GF)[j] == frac(n[j], 1425)),
    x = "Global Frequency (%)", y = NULL
  ) +
  pub_theme +
  theme(panel.grid.major.y = element_line(color = "gray92", linewidth = 0.3))

# Export plot data
write.csv(global_freq, file.path(out_dir, "Fig2A_Global_BarChart_Data.csv"), row.names = FALSE)
write.csv(family_freq, file.path(out_dir, "Fig2_Family_Frequency.csv"), row.names = FALSE)

message("=> Fig. 2 plot generated")

# =============================================================================
# FIGURE 3: Segregation & Correlation (p2b | p3a)
# p2b = Panel A (Family Bubble Chart), p3a = Panel B (Correlation Matrix)
# =============================================================================

# ── p2b: Family-level bubble chart (Panel A of new Fig 3) ──
p2b <- ggplot(family_freq, aes(x = Family_label, y = category_label)) +
  geom_point(aes(size = FF_pct, color = FF_pct), shape = 16, alpha = 0.85) +
  scale_color_viridis_c(option = "magma", direction = -1, name = "Family Freq (%)") +
  scale_size_continuous(range = c(1, 8), breaks = c(5, 10, 20, 40), name = "Family Freq (%)") +
  guides(color = guide_legend(), size = guide_legend()) +
  labs(
    title = "A  Family-Level Mutation Frequencies",
    subtitle = expression(italic(FF)[ij] == frac(n[ij], N[i]) ~ " for Top 10 families"),
    x = NULL, y = NULL
  ) +
  pub_theme +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 40, hjust = 1, size = 7),
        panel.grid.major.y = element_line(color = "gray92", linewidth = 0.3))

# ── p3a: Correlation Matrix (Panel B of new Fig 3) ──
corr_cols <- dose_data %>%
  select(ends_with("_mean"), ems_concentration) %>%
  select(where(~ is.numeric(.) && sum(!is.na(.)) > 2 && var(., na.rm = TRUE) > 0))

raw_names <- names(corr_cols)
clean_corr_names <- tools::toTitleCase(gsub("_mean|_", " ", raw_names))
clean_corr_names <- gsub("Ems Concentration", "EMS Conc.", clean_corr_names)
clean_corr_names <- gsub("M2 Mutation Rate", "M2 Mut. Rate", clean_corr_names)
clean_corr_names <- gsub("Seedling Survival Rate", "Survival Rate", clean_corr_names)
clean_corr_names <- gsub("Germination Rate", "Germ. Rate", clean_corr_names)
clean_corr_names <- gsub("Harvested Pods Plant", "Pods/Plant", clean_corr_names)
clean_corr_names <- gsub("Harvested Seeds Plant Calc", "Seeds/Plant", clean_corr_names)
clean_corr_names <- gsub("Seeds Pod", "Seeds/Pod", clean_corr_names)
clean_corr_names <- gsub("Dead Embryo Ratio", "Dead Embryo", clean_corr_names)
clean_corr_names <- gsub("Fertile Mutation Rate", "Fertile Mut.", clean_corr_names)
clean_corr_names <- gsub("Flowering Rate", "Flower Rate", clean_corr_names)
clean_corr_names <- gsub("Pod Set Rate", "Pod Set", clean_corr_names)
names(corr_cols) <- clean_corr_names

cormat <- cor(corr_cols, use = "pairwise.complete.obs", method = "pearson")
cormat_df <- as.data.frame(as.table(cormat))
names(cormat_df) <- c("Var1", "Var2", "Value")
cormat_df$Var1 <- factor(cormat_df$Var1, levels = clean_corr_names)
cormat_df$Var2 <- factor(cormat_df$Var2, levels = rev(clean_corr_names))
cormat_df$Label <- ifelse(is.na(cormat_df$Value), "", sprintf("%.2f", cormat_df$Value))
cormat_df$Text_Color <- ifelse(!is.na(cormat_df$Value) & abs(cormat_df$Value) > 0.7, "white", "black")

p3a <- ggplot(cormat_df, aes(x = Var1, y = Var2, fill = Value)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = Label, color = Text_Color), size = geom_text_size * 0.9,
            family = "sans", na.rm = TRUE) +
  scale_color_identity() +
  scale_fill_scico(palette = "vik", limits = c(-1, 1), na.value = "gray85",
                    name = "Pearson\nr") +
  labs(title = "B  Trait Correlation Matrix", x = NULL, y = NULL) +
  theme_minimal(base_size = base_font_size, base_family = "sans") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 7),
    axis.text.y = element_text(color = "black", size = 7),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.position = "right",
    legend.key.height = unit(0.8, "cm"),
    legend.key.width = unit(0.25, "cm"),
    plot.title = element_text(face = "bold", size = base_font_size + 1)
  )

# Export plot data
write.csv(cormat_df, file.path(out_dir, "Fig3A_Correlation_Matrix.csv"), row.names = FALSE)

message("=> Fig. 3 plots generated")

# =============================================================================
# FIGURE 4: Trait Trends (p3b alone — Panel A)
# =============================================================================

# ── p3b: Faceted Multi-trait Trends (Panel A of new Fig 4) ──
ems_range <- diff(range(dose_data$ems_concentration, na.rm = TRUE))
err_bar_width <- ifelse(ems_range > 0, ems_range / 30, 0.005)

plot_data_long <- dose_data %>%
  pivot_longer(cols = -ems_concentration, names_to = c("trait_raw", ".value"),
               names_pattern = "(.*)_(mean|sd)$") %>%
  filter(!is.na(mean)) %>%
  mutate(
    ems_pct = ems_concentration * 100,
    trait_clean = tools::toTitleCase(gsub("_", " ", trait_raw)),
    trait_clean = gsub("Seedling Survival Rate", "Survival Rate", trait_clean),
    trait_clean = gsub("Germination Rate", "Germ. Rate", trait_clean),
    trait_clean = gsub("Harvested Pods Plant", "Pods/Plant", trait_clean),
    trait_clean = gsub("Harvested Seeds Plant Calc", "Seeds/Plant", trait_clean),
    trait_clean = gsub("Seeds Pod", "Seeds/Pod", trait_clean),
    trait_clean = gsub("Dead Embryo Ratio", "Dead Embryo", trait_clean),
    trait_clean = gsub("Fertile Mutation Rate", "Fertile Mut.", trait_clean),
    trait_clean = gsub("Flowering Rate", "Flower Rate", trait_clean),
    trait_clean = gsub("Pod Set Rate", "Pod Set", trait_clean),
    trait_clean = gsub("M2 Mutation Rate", "M2 Mut. Rate", trait_clean),
    trait_clean = gsub("M2 Germination Rate", "M2 Germ. Rate", trait_clean)
  )

p3b <- ggplot(plot_data_long, aes(x = ems_pct, y = mean)) +
  geom_line(color = okabe_ito[5], linewidth = 0.5) +
  geom_point(size = 1.5, color = okabe_ito[6]) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                width = err_bar_width * 100, color = "gray50", linewidth = 0.3, na.rm = TRUE) +
  facet_wrap(~trait_clean, scales = "free_y", ncol = 3) +
  labs(
    title = "A  Phenotypic Trait Responses to EMS Treatment",
    x = "EMS Concentration (% v/v)",
    y = "Measured Value"
  ) +
  pub_theme +
  theme(axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6),
        strip.background = element_rect(fill = "transparent", color = NA),
        strip.text = element_text(face = "bold", size = base_font_size - 1))

# Export plot data
write.csv(plot_data_long, file.path(out_dir, "Fig3B_Trait_Trends.csv"), row.names = FALSE)

message("=> Fig. 4 plot generated")

# =============================================================================
# FINAL PDF EXPORTS
# =============================================================================

# Export 1: Fig 1 — Dose-Response (p1a | p1b)
fig1_final <- (p1a | p1b) +
  plot_annotation(
    title = "EMS Dose-Response Analysis in Medicago polymorpha",
    subtitle = "Four-parameter log-logistic model (LL.4) and quadratic mutation rate optimization",
    theme = theme(plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
                  plot.subtitle = element_text(size = 8, hjust = 0.5, color = "gray30"))
  )

ggsave(file.path(out_dir, "Fig1_DoseResponse.pdf"),
       plot = fig1_final, width = 10, height = 5, units = "in", dpi = 300, bg = "white")

# Export 2: Fig 2 (Part) — Global Bar Chart alone (for Illustrator compositing)
ggsave(file.path(out_dir, "Fig2_Bottom_GlobalBarChart.pdf"),
       plot = p2a, width = 8, height = 4, units = "in", dpi = 300, bg = "white")

# Export 3: Fig 3 — Segregation & Correlation (p2b | p3a)
fig3_final <- (p2b | p3a) +
  plot_layout(widths = c(1.2, 1)) +
  plot_annotation(
    title = "Phenotype Segregation and Trait Correlation Analysis",
    subtitle = sprintf("M2 population: N = %d plants, %d families", N_total, length(unique(m2_df$Family))),
    theme = theme(plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
                  plot.subtitle = element_text(size = 8, hjust = 0.5, color = "gray30"))
  )

ggsave(file.path(out_dir, "Fig3_Segregation_Correlation.pdf"),
       plot = fig3_final, width = 12, height = 5, units = "in", dpi = 300, bg = "white")

# Export 4: Fig 4 — Trait Trends (p3b alone)
ggsave(file.path(out_dir, "Fig4_TraitTrends.pdf"),
       plot = p3b, width = 8, height = 6, units = "in", dpi = 300, bg = "white")

message("=== All 4 figures exported as PDF successfully ===")
