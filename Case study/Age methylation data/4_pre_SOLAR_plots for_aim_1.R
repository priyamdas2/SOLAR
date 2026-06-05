# ============================================================
# Figure: Pre-SOLAR residual heterogeneity
# Panels:
#   (a) Residualized Horvath DNAm age by age group
#   (b) Residualized Levine DNAm age by age group
#   (c) Residualized Horvath vs Levine DNAm age
# ============================================================

setwd("U:/SOLAR/Case study/Age methylation data")

# ------------------------------------------------------------
# Packages
# ------------------------------------------------------------

pkg_needed <- c(
  "ggplot2",
  "dplyr",
  "patchwork"
)

for (pkg in pkg_needed) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

if (!requireNamespace("gghalves", quietly = TRUE)) {
  
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }
  
  remotes::install_github("erocoar/gghalves")
}

library(ggplot2)
library(dplyr)
library(patchwork)
library(gghalves)

# ------------------------------------------------------------
# Read phenotype file
# ------------------------------------------------------------

pheno <- read.csv(
  "Pheno_for_MATLAB_with_epi_ages.csv",
  check.names = FALSE
)

plot_dir <- "EpiAge_diagnostic_plots"

if (!dir.exists(plot_dir)) {
  dir.create(plot_dir)
}

# ------------------------------------------------------------
# Prepare data
# ------------------------------------------------------------

set.seed(123)

dat <- pheno %>%
  mutate(
    age_group = factor(
      Age_months,
      levels = sort(unique(Age_months)),
      labels = paste0(sort(unique(Age_months)), " mo")
    )
  ) %>%
  filter(
    !is.na(AgeAccelResidual_Horvath),
    !is.na(AgeAccelResidual_Levine),
    !is.na(age_group)
  ) %>%
  mutate(
    age_num = as.numeric(age_group),
    age_num_jitter = age_num - 0.20 + runif(n(), -0.065, 0.065)
  )

median_dat <- dat %>%
  group_by(age_group, age_num) %>%
  summarise(
    med_horvath = median(AgeAccelResidual_Horvath, na.rm = TRUE),
    med_levine  = median(AgeAccelResidual_Levine, na.rm = TRUE),
    .groups = "drop"
  )

r_hl <- cor(
  dat$AgeAccelResidual_Horvath,
  dat$AgeAccelResidual_Levine,
  use = "complete.obs"
)

# ------------------------------------------------------------
# Common theme
# ------------------------------------------------------------

theme_case <- theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 17, hjust = 0.5),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12, color = "black"),
    axis.line = element_line(linewidth = 0.45),
    axis.ticks = element_line(linewidth = 0.45),
    plot.margin = margin(8, 8, 8, 8),
    legend.position = "none"
  )

# ------------------------------------------------------------
# Colors
# ------------------------------------------------------------

col_horvath <- "#2CA25F"
col_levine  <- "#2B7BDA"
col_point   <- "#303030"

# ------------------------------------------------------------
# Panel (a): Horvath raincloud
# ------------------------------------------------------------

p1 <- ggplot(
  dat,
  aes(
    x = age_num,
    y = AgeAccelResidual_Horvath,
    group = age_group
  )
) +
  gghalves::geom_half_violin(
    aes(fill = age_group),
    side = "r",
    alpha = 0.45,
    color = NA,
    trim = FALSE,
    width = 0.85
  ) +
  geom_point(
    aes(
      x = age_num_jitter,
      y = AgeAccelResidual_Horvath,
      color = age_group
    ),
    alpha = 0.35,
    size = 1.15
  ) +
  geom_segment(
    data = median_dat,
    aes(
      x = age_num - 0.075,
      xend = age_num + 0.075,
      y = med_horvath,
      yend = med_horvath
    ),
    inherit.aes = FALSE,
    linewidth = 0.8,
    color = "black"
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.55,
    color = "gray35"
  ) +
  scale_fill_manual(
    values = rep(col_horvath, length(levels(dat$age_group)))
  ) +
  scale_color_manual(
    values = rep(col_horvath, length(levels(dat$age_group)))
  ) +
  scale_x_continuous(
    breaks = seq_along(levels(dat$age_group)),
    labels = levels(dat$age_group)
  ) +
  labs(
    title = "Residualized Horvath DNAm age",
    x = "Chronological age group",
    y = "Residualized DNAm age (years)"
  ) +
  theme_case

# ------------------------------------------------------------
# Panel (b): Levine raincloud
# ------------------------------------------------------------

p2 <- ggplot(
  dat,
  aes(
    x = age_num,
    y = AgeAccelResidual_Levine,
    group = age_group
  )
) +
  gghalves::geom_half_violin(
    aes(fill = age_group),
    side = "r",
    alpha = 0.45,
    color = NA,
    trim = FALSE,
    width = 0.85
  ) +
  geom_point(
    aes(
      x = age_num_jitter,
      y = AgeAccelResidual_Levine,
      color = age_group
    ),
    alpha = 0.35,
    size = 1.15
  ) +
  geom_segment(
    data = median_dat,
    aes(
      x = age_num - 0.075,
      xend = age_num + 0.075,
      y = med_levine,
      yend = med_levine
    ),
    inherit.aes = FALSE,
    linewidth = 0.8,
    color = "black"
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.55,
    color = "gray35"
  ) +
  scale_fill_manual(
    values = rep(col_levine, length(levels(dat$age_group)))
  ) +
  scale_color_manual(
    values = rep(col_levine, length(levels(dat$age_group)))
  ) +
  scale_x_continuous(
    breaks = seq_along(levels(dat$age_group)),
    labels = levels(dat$age_group)
  ) +
  labs(
    title = "Residualized Levine/DNAm PhenoAge",
    x = "Chronological age group",
    y = "Residualized DNAm age (years)"
  ) +
  theme_case

# ------------------------------------------------------------
# Panel (c): Horvath vs Levine cross-clock residual agreement
# ------------------------------------------------------------

p3 <- ggplot(
  dat,
  aes(
    x = AgeAccelResidual_Levine,
    y = AgeAccelResidual_Horvath
  )
) +
  geom_point(
    color = col_point,
    alpha = 0.35,
    size = 1.6
  ) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    linewidth = 0.8,
    color = "black",
    fill = "gray80"
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.45,
    color = "gray55"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.45,
    color = "gray55"
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = paste0("r = ", sprintf("%.2f", r_hl)),
    hjust = 1.15,
    vjust = 1.5,
    size = 4.3
  ) +
  labs(
    title = "Cross-clock residualized\nDNAm age concordance",
    x = "Residualized Levine/DNAm PhenoAge (years)",
    y = "Residualized Horvath DNAm age (years)"
  ) +
  theme_case

# ------------------------------------------------------------
# Combine 1 x 3 panel
# ------------------------------------------------------------

fig_presolar <- p1 + p2 + p3 +
  plot_layout(nrow = 1, widths = c(1, 1, 1.05))

# Display plot in R plotting window
print(fig_presolar)

# ------------------------------------------------------------
# Save combined figure
# ------------------------------------------------------------

ggsave(
  filename = file.path(plot_dir, "Plot_PreSOLAR_residual_heterogeneity_1x3.png"),
  plot = fig_presolar,
  width = 15.5,
  height = 4.8,
  dpi = 450,
  bg = "white"
)

ggsave(
  filename = file.path(plot_dir, "Plot_PreSOLAR_residual_heterogeneity_1x3.pdf"),
  plot = fig_presolar,
  width = 15.5,
  height = 4.8,
  device = cairo_pdf,
  bg = "white"
)

# ------------------------------------------------------------
# Save separate panels too
# ------------------------------------------------------------

ggsave(
  filename = file.path(plot_dir, "Plot_Horvath_residual_raincloud.png"),
  plot = p1,
  width = 5.2,
  height = 4.8,
  dpi = 450,
  bg = "white"
)

ggsave(
  filename = file.path(plot_dir, "Plot_Levine_residual_raincloud.png"),
  plot = p2,
  width = 5.2,
  height = 4.8,
  dpi = 450,
  bg = "white"
)

ggsave(
  filename = file.path(plot_dir, "Plot_Horvath_Levine_residual_scatter.png"),
  plot = p3,
  width = 5.4,
  height = 4.8,
  dpi = 450,
  bg = "white"
)

cat("Saved plots in folder:", plot_dir, "\n")