# ============================================================
# Figure 2: Main SOLAR result for Horvath residual acceleration
#
# Expected folder structure:
#
# Current working directory/
# ├── Output/
# │   ├── Output_from_M1_summary.csv
# │   ├── Output_from_M1_factor_strength.csv
# │   ├── Output_from_M1_latent_scores_and_fitted.csv
# │   ├── Output_from_M1_optimization_trace.csv
# │   └── Output_from_M1_top_CpG_importance.csv
# └── Age methylation data/
#     └── X_cpg_names.rds
#
# All plots are saved as PNG files inside Output/
# ============================================================

# ------------------------------------------------------------
# Packages
# ------------------------------------------------------------
setwd("U:/SOLAR/Case study")
pkg_needed <- c(
  "ggplot2",
  "dplyr",
  "tidyr",
  "readr",
  "stringr",
  "scales",
  "forcats",
  "ggrepel",
  "patchwork"
)

for (pkg in pkg_needed) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(scales)
library(forcats)
library(ggrepel)
library(patchwork)

# ------------------------------------------------------------
# Directories
# ------------------------------------------------------------

out_dir <- "Output"

if (!dir.exists(out_dir)) {
  stop("Output folder not found.")
}

# ------------------------------------------------------------
# Read SOLAR outputs
# ------------------------------------------------------------

summary_tbl <- read_csv(
  file.path(out_dir, "Output_from_M1_summary.csv"),
  show_col_types = FALSE
)

factor_tbl <- read_csv(
  file.path(out_dir, "Output_from_M1_factor_strength.csv"),
  show_col_types = FALSE
)

score_tbl <- read_csv(
  file.path(out_dir, "Output_from_M1_latent_scores_and_fitted.csv"),
  show_col_types = FALSE
)

trace_tbl <- read_csv(
  file.path(out_dir, "Output_from_M1_optimization_trace.csv"),
  show_col_types = FALSE
)

top_cpg_tbl <- read_csv(
  file.path(out_dir, "Output_from_M1_top_CpG_importance.csv"),
  show_col_types = FALSE
)

top_cpg_tbl <- top_cpg_tbl %>%
  dplyr::mutate(
    relative_importance = importance / max(importance, na.rm = TRUE)
  )

# ------------------------------------------------------------
# Optional CpG names
# ------------------------------------------------------------

cpg_name_file <- file.path("Age methylation data", "X_cpg_names.rds")

if (file.exists(cpg_name_file)) {
  
  cpg_names <- readRDS(cpg_name_file)
  
  if (length(cpg_names) >= max(top_cpg_tbl$CpG_index, na.rm = TRUE)) {
    top_cpg_tbl <- top_cpg_tbl %>%
      dplyr::mutate(CpG = cpg_names[CpG_index])
  } else {
    warning("X_cpg_names.rds exists but length is smaller than max CpG_index.")
    top_cpg_tbl <- top_cpg_tbl %>%
      dplyr::mutate(CpG = paste0("CpG_", CpG_index))
  }
  
} else {
  
  top_cpg_tbl <- top_cpg_tbl %>%
    dplyr::mutate(CpG = paste0("CpG_", CpG_index))
}

# ------------------------------------------------------------
# Extract key quantities
# ------------------------------------------------------------

q_hat <- summary_tbl$q_hat[1]
rmse_val <- summary_tbl$rmse_in_sample[1]
corr_val <- summary_tbl$corr_in_sample[1]
r2_val <- summary_tbl$R2_in_sample[1]
runtime_val <- summary_tbl$total_runtime_sec[1]

score_tbl <- score_tbl %>%
  dplyr::mutate(
    Age_group = factor(
      Age_months,
      levels = sort(unique(Age_months)),
      labels = paste0(sort(unique(Age_months)), " mo")
    )
  )

factor_cols <- grep("^H_factor_", names(score_tbl), value = TRUE)

# ------------------------------------------------------------
# Common plotting theme
# ------------------------------------------------------------

theme_main <- theme_classic(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12.5, color = "black"),
    axis.line = element_line(linewidth = 0.45),
    axis.ticks = element_line(linewidth = 0.45),
    legend.title = element_text(size = 12.5),
    legend.text = element_text(size = 11.5),
    legend.position = "right",
    plot.margin = margin(10, 12, 10, 12)
  )

theme_no_legend <- theme_main +
  theme(legend.position = "none")

col_main <- "#2CA25F"
col_secondary <- "#2B7BDA"
col_dark <- "#303030"
col_gray <- "gray70"

# ============================================================
# Plot 1: Selected rank and factor strength
# ============================================================
factor_long <- factor_tbl %>%
  dplyr::mutate(
    factor_label = paste0("Factor ", factor),
    factor_label = factor(factor_label, levels = paste0("Factor ", factor))
  )

p_factor_strength <- ggplot(
  factor_long,
  aes(
    x = factor_label,
    y = d_squared_abs_beta_prop
  )
) +
  geom_col(
    width = 0.68,
    fill = col_main,
    alpha = 0.85
  ) +
  geom_line(
    aes(y = d_squared_prop, group = 1),
    linewidth = 0.75,
    color = "black"
  ) +
  geom_point(
    aes(y = d_squared_prop),
    size = 3.1,
    color = "black"
  ) +
  annotate(
    "label",
    x = 2.8,
    y = 0.86,
    label = "Bars: supervised contribution\nLine: methylation factor strength",
    hjust = 0,
    vjust = 1,
    size = 3.7,
    linewidth = 0.35,
    label.r = unit(0.10, "lines"),
    label.padding = unit(0.22, "lines"),
    fill = "white",
    color = "black"
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 0.95),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title = "Factor strength (Horvath)",
    x = NULL,
    y = "Relative contribution"
  ) +
  theme_main +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 15.5,
      hjust = 0.5
    ),
    legend.position = "none",
    axis.text.x = element_text(size = 10.5),
    axis.text.y = element_text(size = 10.5),
    axis.title.y = element_text(size = 12.0),
    panel.border = element_blank()
  )

p_factor_strength

ggsave(
  filename = file.path(out_dir, "Plot_M1_factor_strength.png"),
  plot = p_factor_strength,
  width = 5.1,
  height = 4.2,
  dpi = 450,
  bg = "white"
)

# ============================================================
# Plot 2: Latent score map, Factors 1 and 2
# Colored by Horvath residual acceleration
# ============================================================

p_score_continuous <- ggplot(
  score_tbl,
  aes(
    x = H_factor_01,
    y = H_factor_02,
    color = y_observed
  )
) +
  geom_point(
    size = 2.1,
    alpha = 0.82
  ) +
  scale_color_gradient2(
    low = "#2B7BDA",
    mid = "white",
    high = "#D95F02",
    midpoint = 0,
    name = "Observed residualized\nDNAm age (years)"
  ) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.35,
    linetype = "dashed",
    color = "gray65"
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.35,
    linetype = "dashed",
    color = "gray65"
  ) +
  labs(
    title = "SOLAR latent score map (Horvath)",
    x = "Latent factor 1 score",
    y = "Latent factor 2 score"
  ) +
  theme_main

ggsave(
  filename = file.path(out_dir, "Plot_M1_latent_score_map_colored_by_residual.png"),
  plot = p_score_continuous,
  width = 7.2,
  height = 6.0,
  dpi = 450,
  bg = "white"
)

# ============================================================
# Plot 3: Latent score map, Factors 1 and 2
# Colored by developmental age group
# ============================================================

p_score_age <- ggplot(
  score_tbl,
  aes(
    x = H_factor_01,
    y = H_factor_02,
    color = Age_group
  )
) +
  geom_point(
    size = 2.1,
    alpha = 0.78
  ) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.35,
    linetype = "dashed",
    color = "gray65"
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.35,
    linetype = "dashed",
    color = "gray65"
  ) +
  labs(
    title = "Latent score map by developmental\nage group (Horvath)",
    x = "Latent factor 1 score",
    y = "Latent factor 2 score",
    color = "Age group"
  ) +
  theme_main

ggsave(
  filename = file.path(out_dir, "Plot_M1_latent_score_map_colored_by_age.png"),
  plot = p_score_age,
  width = 7.2,
  height = 6.0,
  dpi = 450,
  bg = "white"
)

# ============================================================
# Plot 4: Developmental trajectories of latent factors
# Mean +/- standard error by age group
# ============================================================

score_long <- score_tbl %>%
  dplyr::select(sample_index, Age_months, Age_group, all_of(factor_cols)) %>%
  pivot_longer(
    cols = all_of(factor_cols),
    names_to = "factor",
    values_to = "score"
  ) %>%
  dplyr::mutate(
    factor_num = as.integer(str_extract(factor, "\\d+")),
    factor_label = paste0("Factor ", factor_num)
  )

traj_tbl <- score_long %>%
  group_by(Age_months, Age_group, factor_num, factor_label) %>%
  summarise(
    mean_score = mean(score, na.rm = TRUE),
    se_score = sd(score, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

p_traj <- ggplot(
  traj_tbl,
  aes(
    x = Age_months,
    y = mean_score,
    color = factor_label,
    group = factor_label
  )
) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.35,
    linetype = "dashed",
    color = "gray70"
  ) +
  geom_errorbar(
    aes(
      ymin = mean_score - 1.96 * se_score,
      ymax = mean_score + 1.96 * se_score
    ),
    width = 2.0,
    linewidth = 0.45,
    alpha = 0.75
  ) +
  geom_line(
    linewidth = 0.95,
    alpha = 0.90
  ) +
  geom_point(
    size = 2.8,
    alpha = 0.95
  ) +
  scale_x_continuous(
    breaks = sort(unique(score_tbl$Age_months)),
    labels = paste0(sort(unique(score_tbl$Age_months)), "")
  ) +
  guides(
    color = guide_legend(
      nrow = 2,
      byrow = TRUE,
      override.aes = list(linewidth = 0.9, size = 2.8)
    )
  ) +
  labs(
    title = "Developmental latent-factor\ntrajectories (Horvath)",
    x = "Chronological age group (months)",
    y = "Mean latent factor score",
    color = NULL
  ) +
  theme_main +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 15.5,
      hjust = 0.5
    ),
    legend.position = "top",
    legend.box = "horizontal",
    legend.background = element_rect(
      color = "gray75",
      fill = "white",
      linewidth = 0.30
    ),
    legend.margin = margin(3, 5, 3, 5),
    legend.key.width = unit(1.1, "cm"),
    legend.key.height = unit(0.35, "cm"),
    legend.text = element_text(size = 9.5)
  )

p_traj

ggsave(
  filename = file.path(out_dir, "Plot_M1_latent_factor_developmental_profiles.png"),
  plot = p_traj,
  width = 5.1,
  height = 4.2,
  dpi = 450,
  bg = "white"
)

# ============================================================
# Plot 5: Developmental distributions of first few latent factors
# Optional alternative to line trajectories
# ============================================================

top_factor_show <- min(3, q_hat)

score_long_top <- score_long %>%
  filter(factor_num <= top_factor_show)

p_factor_box <- ggplot(
  score_long_top,
  aes(
    x = Age_group,
    y = score,
    fill = Age_group
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.45,
    color = "gray55"
  ) +
  geom_violin(
    width = 0.85,
    alpha = 0.45,
    color = NA,
    trim = FALSE
  ) +
  geom_boxplot(
    width = 0.16,
    outlier.shape = NA,
    fill = "white",
    color = "black",
    linewidth = 0.40
  ) +
  facet_wrap(
    ~ factor_label,
    nrow = 1,
    scales = "free_y"
  ) +
  labs(
    title = "Latent-factor distributions across age groups (Horvath)",
    x = "Chronological age group",
    y = "Latent factor score"
  ) +
  theme_no_legend +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 13.5)
  )

ggsave(
  filename = file.path(out_dir, "Plot_M1_latent_factor_distributions_by_age.png"),
  plot = p_factor_box,
  width = 15.4,
  height = 4.5,
  dpi = 450,
  bg = "white"
)

# ============================================================
# Plot 6: Observed vs fitted Horvath residual acceleration
# ============================================================

label_metrics <- paste0(
  "corr. = ", sprintf("%.2f", corr_val),
  "\nRMSE = ", sprintf("%.2f", rmse_val)
)

p_obs_fit <- ggplot(
  score_tbl,
  aes(
    x = y_fitted,
    y = y_observed
  )
) +
  geom_point(
    size = 1.9,
    alpha = 0.45,
    color = col_dark
  ) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    linewidth = 0.85,
    color = "black",
    fill = "gray80"
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    linewidth = 0.60,
    color = "gray45"
  ) +
  annotate(
    "label",
    x = Inf,
    y = -Inf,
    label = label_metrics,
    hjust = 1.05,
    vjust = -0.25,
    size = 3.6,
    fill = "white",
    color = "black",
    linewidth = 0.35,
    label.padding = unit(0.18, "lines"),
    label.r = unit(0.12, "lines")
  ) +
  labs(
    title = "Observed versus fitted\nresidualized Horvath DNAm age",
    x = "Fitted residualized DNAm age (years)",
    y = "Observed residualized DNAm age (years)"
  ) +
  theme_main +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 15.5,
      hjust = 0.5
    ),
    axis.text.x = element_text(size = 10.5),
    axis.text.y = element_text(size = 10.5),
    axis.title.x = element_text(size = 12.0),
    axis.title.y = element_text(size = 12.0),
    panel.border = element_blank()
  )

p_obs_fit

ggsave(
  filename = file.path(out_dir, "Plot_M1_observed_vs_fitted.png"),
  plot = p_obs_fit,
  width = 5.1,
  height = 4.2,
  dpi = 450,
  bg = "white"
)
# ============================================================
# Plot 7: Residual diagnostics by fitted value
# Optional supporting plot
# ============================================================

p_resid_fit <- ggplot(
  score_tbl,
  aes(
    x = y_fitted,
    y = residual
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.55,
    color = "gray50"
  ) +
  geom_point(
    size = 1.8,
    alpha = 0.42,
    color = col_dark
  ) +
  geom_smooth(
    method = "loess",
    formula = y ~ x,
    se = TRUE,
    linewidth = 0.80,
    color = "black",
    fill = "gray82"
  ) +
  labs(
    title = "Residuals versus fitted values (Horvath)",
    x = "Fitted residualized DNAm age (years)",
    y = "Model residual (years)"
  ) +
  theme_main

ggsave(
  filename = file.path(out_dir, "Plot_M1_residuals_vs_fitted.png"),
  plot = p_resid_fit,
  width = 6.4,
  height = 5.4,
  dpi = 450,
  bg = "white"
)

# ============================================================
# Plot 8: Top CpG importance scores
# Optional if you decide to preview Figure 3 style
# ============================================================

top_n_plot <- 15

top_plot_tbl <- top_cpg_tbl %>%
  slice_head(n = top_n_plot) %>%
  dplyr::mutate(
    CpG_label = ifelse(is.na(CpG), paste0("CpG_", CpG_index), CpG),
    CpG_label = fct_reorder(CpG_label, relative_importance)
  )

p_top_cpg <- ggplot(
  top_plot_tbl,
  aes(
    x = CpG_label,
    y = relative_importance
  )
) +
  geom_col(
    width = 0.72,
    fill = col_main,
    alpha = 0.88
  ) +
  coord_flip() +
  labs(
    title = "Top CpGs by importance score (Horvath)",
    x = NULL,
    y = "Relative importance score"
  ) +
  theme_no_legend +
  theme(
    plot.title = element_text(
      size = 15,
      face = "bold"
    ),
    axis.title.x = element_text(
      size = 15
    ),
    axis.title.y = element_text(
      size = 15
    ),
    axis.text.x = element_text(
      size = 12
    ),
    axis.text.y = element_text(
      size = 12
    )
  )

print(p_top_cpg)

ggsave(
  filename = file.path(out_dir, "Plot_M1_top_CpG_importance_preview.png"),
  plot = p_top_cpg,
  width = 7.0,
  height = 4.8,
  dpi = 450,
  bg = "white"
)

# ============================================================
# Plot 9: Loading heatmap for top CpGs
# Optional preview for Figure 3
# ============================================================

loading_cols <- grep("^loading_V_factor_", names(top_cpg_tbl), value = TRUE)

top_heat_tbl <- top_cpg_tbl %>%
  slice_head(n = min(15, nrow(top_cpg_tbl))) %>%
  dplyr::mutate(
    CpG_label = ifelse(is.na(CpG), paste0("CpG_", CpG_index), CpG),
    CpG_label = factor(CpG_label, levels = rev(CpG_label))
  ) %>%
  dplyr::select(CpG_label, all_of(loading_cols)) %>%
  pivot_longer(
    cols = all_of(loading_cols),
    names_to = "factor",
    values_to = "loading"
  ) %>%
  dplyr::mutate(
    factor_num = as.integer(str_extract(factor, "\\d+")),
    factor_label = paste0("Factor ", factor_num),
    factor_label = factor(
      factor_label,
      levels = paste0("Factor ", sort(unique(factor_num)))
    )
  )

p_loading_heat <- ggplot(
  top_heat_tbl,
  aes(
    x = factor_label,
    y = CpG_label,
    fill = loading
  )
) +
  geom_tile(color = "gray92",
            linewidth = 0.15) +
  scale_fill_gradient2(
    low = "#2B7BDA",
    mid = "white",
    high = "#D95F02",
    midpoint = 0,
    name = "Factor\nloading"
  ) +
  labs(
    title = "CpG loading patterns across selected\nlatent factors (Horvath)",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    axis.text.x = element_text(size = 13.1, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    
    legend.background = element_rect(
      color = "gray70",
      fill = "white",
      linewidth = 0.45
    ),
    
    legend.box.background = element_rect(
      color = "gray70",
      fill = "white",
      linewidth = 0.45
    ),
    
    legend.margin = margin(5, 6, 9, 6),
    
    panel.grid = element_blank(),
    plot.margin = margin(10, 12, 10, 12)
  )

print(p_loading_heat)

ggsave(
  filename = file.path(out_dir, "Plot_M1_top_CpG_loading_heatmap_preview.png"),
  plot = p_loading_heat,
  width = 7.0,
  height = 4.8,
  dpi = 450,
  bg = "white"
)

# ============================================================
# Plot 10: Combined Figure 2 panel from saved PNG files
# Row 1: three separate panels
# Row 2: latent-factor distribution panel spanning full width
# ============================================================

pkg_needed_img <- c("magick", "grid", "png")

for (pkg in pkg_needed_img) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(magick)

# ------------------------------------------------------------
# Read saved PNG files
# ------------------------------------------------------------

img_factor <- image_read(
  file.path(out_dir, "Plot_M1_factor_strength.png")
)

img_traj <- image_read(
  file.path(out_dir, "Plot_M1_latent_factor_developmental_profiles.png")
)

img_obsfit <- image_read(
  file.path(out_dir, "Plot_M1_observed_vs_fitted.png")
)

img_dist <- image_read(
  file.path(out_dir, "Plot_M1_latent_factor_distributions_by_age.png")
)

# ------------------------------------------------------------
# Standardize dimensions
# ------------------------------------------------------------

single_w <- 2300
single_h <- 1890

row2_w <- single_w * 3
row2_h <- 2025

img_factor <- image_resize(img_factor, paste0(single_w, "x", single_h, "!"))
img_traj   <- image_resize(img_traj,   paste0(single_w, "x", single_h, "!"))
img_obsfit <- image_resize(img_obsfit, paste0(single_w, "x", single_h, "!"))

img_dist <- image_resize(img_dist, paste0(row2_w, "x", row2_h, "!"))

# ------------------------------------------------------------
# Create row 1 and row 2
# ------------------------------------------------------------

row1 <- image_append(
  c(img_factor, img_traj, img_obsfit),
  stack = FALSE
)

# Add a small white separator between rows
separator <- image_blank(
  width = row2_w,
  height = 150,
  color = "white"
)

fig2_combined <- image_append(
  c(row1, separator, img_dist, separator),
  stack = TRUE
)

# ------------------------------------------------------------
# Save combined panel
# ------------------------------------------------------------

image_write(
  fig2_combined,
  path = file.path(out_dir, "Plot_M1_combined_Figure2_panel.png"),
  format = "png",
  density = "450x450"
)

# Display in RStudio viewer / plot window
print(fig2_combined)

cat("\nSaved Figure 2 candidate plots in folder:", out_dir, "\n")
cat("Selected rank:", q_hat, "\n")
cat("RMSE:", round(rmse_val, 4), "\n")
cat("Correlation:", round(corr_val, 4), "\n")
cat("R2:", round(r2_val, 4), "\n")
cat("Runtime (sec):", round(runtime_val, 2), "\n")
