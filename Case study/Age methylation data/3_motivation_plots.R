# ============================================================
# SOLAR Section 2 / pre-SOLAR real-data plots
# Creates:
# 1. Combined clock-age boxplot, width:height = 3:1
# 2. Clock calibration correlation bar plot, width:height = 5:4
# 3. Age-group sample-size bar plot, width:height = 5:4
# ============================================================

setwd("U:/SOLAR/Case study/Age methylation data")

cat("Reading phenotype file...\n")

pheno <- read.csv(
  "Pheno_for_MATLAB_with_epi_ages.csv",
  check.names = FALSE
)

chron_age_years <- pheno$Age_months / 12

plot_dir <- "Motivation plots"
if (!dir.exists(plot_dir)) dir.create(plot_dir)

# ------------------------------------------------------------
# Color theme matched to methylation-clock overview figure
# ------------------------------------------------------------

col_horvath <- "#2ca25f"   # green
col_levine  <- "#2b7bde"   # blue
col_hannum  <- "#8e44ad"   # purple
col_orange  <- "#e85d1f"   # age acceleration orange
col_gray    <- "#4d4d4d"

age_group <- factor(
  pheno$Age_months,
  levels = c(3, 9, 48, 72),
  labels = c("3 mo", "9 mo", "48 mo", "72 mo")
)

# ============================================================
# 1. Combined boxplot: DNAm clock ages by developmental window
#    width:height = 3:1
# ============================================================

png(
  file.path(plot_dir, "RealData_clock_age_boxplots.png"),
  width = 3400,
  height = 1000,
  res = 300
)

par(
  mfrow = c(1, 3),
  mar = c(4.5, 4.8, 3.0, 1.2),
  oma = c(0, 0, 0, 0),
  cex.axis = 1.35,
  cex.lab = 1.55,
  cex.main = 1.85
)

boxplot(
  pheno$DNAmAge_Horvath ~ age_group,
  col = adjustcolor(col_horvath, alpha.f = 0.22),
  border = col_horvath,
  boxwex = 0.55,
  lwd = 1.6,
  medlwd = 2.2,
  outline = FALSE,
  xlab = "Age group",
  ylab = "DNAm age (years)",
  main = "Horvath"
)
stripchart(
  pheno$DNAmAge_Horvath ~ age_group,
  vertical = TRUE,
  method = "jitter",
  pch = 16,
  cex = 0.45,
  col = adjustcolor(col_horvath, alpha.f = 0.30),
  add = TRUE
)

boxplot(
  pheno$DNAmPhenoAge_Levine ~ age_group,
  col = adjustcolor(col_levine, alpha.f = 0.22),
  border = col_levine,
  boxwex = 0.55,
  lwd = 1.6,
  medlwd = 2.2,
  outline = FALSE,
  xlab = "Age group",
  ylab = "DNAm age (years)",
  main = "Levine"
)
stripchart(
  pheno$DNAmPhenoAge_Levine ~ age_group,
  vertical = TRUE,
  method = "jitter",
  pch = 16,
  cex = 0.45,
  col = adjustcolor(col_levine, alpha.f = 0.30),
  add = TRUE
)

boxplot(
  pheno$DNAmAge_Hannum ~ age_group,
  col = adjustcolor(col_hannum, alpha.f = 0.22),
  border = col_hannum,
  boxwex = 0.55,
  lwd = 1.6,
  medlwd = 2.2,
  outline = FALSE,
  xlab = "Age group",
  ylab = "DNAm age (years)",
  main = "Hannum"
)
stripchart(
  pheno$DNAmAge_Hannum ~ age_group,
  vertical = TRUE,
  method = "jitter",
  pch = 16,
  cex = 0.45,
  col = adjustcolor(col_hannum, alpha.f = 0.30),
  add = TRUE
)

dev.off()

# ============================================================
# 2. Clock calibration bar plot
#    width:height = 5:4
# ============================================================

clock_corr <- c(
  "Horvath age" = cor(
    pheno$DNAmAge_Horvath,
    chron_age_years,
    use = "complete.obs"
  ),
  
  "Levine age" = cor(
    pheno$DNAmPhenoAge_Levine,
    chron_age_years,
    use = "complete.obs"
  ),
  
  "Hannum age" = cor(
    pheno$DNAmAge_Hannum,
    chron_age_years,
    use = "complete.obs"
  )
)

png(
  file.path(plot_dir, "RealData_clock_calibration_bar_5to4.png"),
  width = 2000,
  height = 1600,
  res = 300
)

par(
  mar = c(5.2, 5.8, 1.0, 1.0),
  mgp = c(3.8, 1.2, 0),
  cex.axis = 1.55,
  cex.lab = 1.75
)

bp <- barplot(
  clock_corr,
  ylim = c(0, 1.05),
  col = c(col_horvath, col_levine, col_hannum),
  border = NA,
  ylab = "Correlation with chron. age",
  xlab = "",
  las = 1
)

text(
  x = bp,
  y = clock_corr + 0.04,
  labels = sprintf("%.2f", clock_corr),
  cex = 1.85,
  font = 2
)

abline(h = seq(0, 1, by = 0.2), col = adjustcolor("gray70", 0.45), lty = 3)

dev.off()

# ============================================================
# 3. Age-group sample-size bar plot
#    width:height = 5:4
# ============================================================

age_counts <- table(age_group)

png(
  file.path(plot_dir, "RealData_age_group_counts_5to4.png"),
  width = 2000,
  height = 1600,
  res = 300
)

par(
  mar = c(5.2, 6.6, 1.0, 1.0),
  mgp = c(3.8, 1.2, 0),
  cex.axis = 1.55,
  cex.lab = 1.75
)

bp2 <- barplot(
  age_counts,
  ylim = c(0, max(age_counts) * 1.18),
  col = col_orange,
  border = NA,
  ylab = "Number of samples",
  xlab = "Age group",
  las = 1
)

text(
  x = bp2,
  y = as.numeric(age_counts) + max(age_counts) * 0.04,
  labels = as.numeric(age_counts),
  cex = 1.85,
  font = 2
)

abline(h = pretty(c(0, max(age_counts))), col = adjustcolor("gray70", 0.45), lty = 3)

dev.off()

cat("Saved plots in:", plot_dir, "\n")