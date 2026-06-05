# ============================================================
# Diagnostics for epigenetic age calculations
# Standalone script
# ============================================================

setwd("U:/SOLAR/Case study/Age methylation data")

# ------------------------------------------------------------
# Read saved objects
# ------------------------------------------------------------

cat("Reading saved RDS and phenotype files...\n")

X_clock <- readRDS("X_clock_CpGs_by_samples_4dp.rds")

pheno <- read.csv(
  "Pheno_for_MATLAB_with_epi_ages.csv",
  check.names = FALSE
)

chron_age_years <- pheno$Age_months / 12

# ============================================================
# Diagnostics
# ============================================================

cat("\n============================================================\n")
cat("Diagnostics for epigenetic age calculations\n")
cat("============================================================\n")

# ------------------------------------------------------------
# 1. Basic matrix diagnostics
# ------------------------------------------------------------

cat("\n[1] X_clock dimensions:\n")
print(dim(X_clock))

cat("\n[2] First few CpG row names:\n")
print(rownames(X_clock)[1:10])

cat("\n[3] First few sample column names:\n")
print(colnames(X_clock)[1:10])

cat("\n[4] Beta-value range:\n")
print(range(X_clock, na.rm = TRUE))

cat("\n[5] Number of missing beta values:\n")
print(sum(is.na(X_clock)))

cat("\n[6] Proportion of missing beta values:\n")
print(mean(is.na(X_clock)))

# ------------------------------------------------------------
# 2. Phenotype and chronological-age diagnostics
# ------------------------------------------------------------

cat("\n[7] Phenotype dimensions:\n")
print(dim(pheno))

cat("\n[8] Age_months summary:\n")
print(summary(pheno$Age_months))

cat("\n[9] Chronological age in years summary:\n")
print(summary(chron_age_years))

cat("\n[10] Check whether existing pheno$age equals Age_months / 12:\n")

if ("age" %in% names(pheno)) {
  print(summary(pheno$age - chron_age_years))
} else {
  cat("No variable named 'age' found in phenotype file.\n")
}

# ------------------------------------------------------------
# 3. Epigenetic-age summaries
# ------------------------------------------------------------

cat("\n[11] Epigenetic age summaries:\n")

epi_vars <- c(
  "DNAmAge_Horvath",
  "DNAmAge_Hannum",
  "DNAmPhenoAge_Levine"
)

for (v in epi_vars) {
  if (v %in% names(pheno)) {
    cat("\n", v, ":\n", sep = "")
    print(summary(pheno[[v]]))
  }
}

cat("\n[12] Age-acceleration difference summaries:\n")

diff_vars <- c(
  "AgeAccelDiff_Horvath",
  "AgeAccelDiff_Hannum",
  "AgeAccelDiff_Levine"
)

for (v in diff_vars) {
  if (v %in% names(pheno)) {
    cat("\n", v, ":\n", sep = "")
    print(summary(pheno[[v]]))
  }
}

cat("\n[13] Residualized age-acceleration summaries:\n")

resid_vars <- c(
  "AgeAccelResidual_Horvath",
  "AgeAccelResidual_Hannum",
  "AgeAccelResidual_Levine"
)

for (v in resid_vars) {
  if (v %in% names(pheno)) {
    cat("\n", v, ":\n", sep = "")
    print(summary(pheno[[v]]))
  }
}

# ------------------------------------------------------------
# 4. Check agreement with methylclock-generated columns
# ------------------------------------------------------------

cat("\n[14] Agreement between manually computed and methylclock-generated columns:\n")

safe_cor <- function(a, b) {
  if (all(is.na(a)) || all(is.na(b))) return(NA_real_)
  cor(a, b, use = "complete.obs")
}

comparison_pairs <- list(
  c("DNAmAge_Horvath", "Horvath"),
  c("DNAmAge_Hannum", "Hannum"),
  c("DNAmPhenoAge_Levine", "Levine"),
  c("AgeAccelDiff_Horvath", "ageAcc.Horvath"),
  c("AgeAccelDiff_Hannum", "ageAcc.Hannum"),
  c("AgeAccelDiff_Levine", "ageAcc.Levine"),
  c("AgeAccelResidual_Horvath", "ageAcc2.Horvath"),
  c("AgeAccelResidual_Hannum", "ageAcc2.Hannum"),
  c("AgeAccelResidual_Levine", "ageAcc2.Levine")
)

for (pair in comparison_pairs) {
  
  v1 <- pair[1]
  v2 <- pair[2]
  
  if (all(c(v1, v2) %in% names(pheno))) {
    
    cat("\nCorrelation:", v1, "vs", v2, "\n")
    
    print(
      safe_cor(
        pheno[[v1]],
        pheno[[v2]]
      )
    )
    
  } else {
    
    cat(
      "\nSkipping:",
      v1,
      "vs",
      v2,
      "because one or both variables are missing.\n"
    )
  }
}

# ------------------------------------------------------------
# 5. Correlation with chronological age
# ------------------------------------------------------------

cat("\n[15] Correlations with chronological age:\n")

age_assoc_vars <- c(
  "DNAmAge_Horvath",
  "DNAmAge_Hannum",
  "DNAmPhenoAge_Levine",
  "AgeAccelDiff_Horvath",
  "AgeAccelDiff_Hannum",
  "AgeAccelDiff_Levine",
  "AgeAccelResidual_Horvath",
  "AgeAccelResidual_Hannum",
  "AgeAccelResidual_Levine"
)

for (v in age_assoc_vars) {
  
  if (v %in% names(pheno)) {
    
    cat(
      "\nCorrelation: ",
      v,
      " vs chronological age\n",
      sep = ""
    )
    
    print(
      safe_cor(
        pheno[[v]],
        chron_age_years
      )
    )
  }
}

# ------------------------------------------------------------
# 6. Save diagnostic summaries to text file
# ------------------------------------------------------------

cat("\n[16] Writing diagnostic report to text file...\n")

sink("epi_age_diagnostics_report.txt")

cat("Diagnostics for epigenetic age calculations\n")
cat("Generated on:", as.character(Sys.time()), "\n\n")

cat("X_clock dimensions:\n")
print(dim(X_clock))

cat("\nBeta-value range:\n")
print(range(X_clock, na.rm = TRUE))

cat("\nNumber of missing beta values:\n")
print(sum(is.na(X_clock)))

cat("\nProportion of missing beta values:\n")
print(mean(is.na(X_clock)))

cat("\nPhenotype dimensions:\n")
print(dim(pheno))

cat("\nAge_months summary:\n")
print(summary(pheno$Age_months))

cat("\nChronological age in years summary:\n")
print(summary(chron_age_years))

if ("age" %in% names(pheno)) {
  
  cat("\nCheck pheno$age - Age_months/12:\n")
  
  print(summary(pheno$age - chron_age_years))
}

cat("\nEpigenetic age summaries:\n")

for (v in epi_vars) {
  
  if (v %in% names(pheno)) {
    
    cat("\n", v, ":\n", sep = "")
    print(summary(pheno[[v]]))
  }
}

cat("\nAge-acceleration difference summaries:\n")

for (v in diff_vars) {
  
  if (v %in% names(pheno)) {
    
    cat("\n", v, ":\n", sep = "")
    print(summary(pheno[[v]]))
  }
}

cat("\nResidualized age-acceleration summaries:\n")

for (v in resid_vars) {
  
  if (v %in% names(pheno)) {
    
    cat("\n", v, ":\n", sep = "")
    print(summary(pheno[[v]]))
  }
}

cat("\nAgreement between manually computed and methylclock-generated columns:\n")

for (pair in comparison_pairs) {
  
  v1 <- pair[1]
  v2 <- pair[2]
  
  if (all(c(v1, v2) %in% names(pheno))) {
    
    cat("\nCorrelation:", v1, "vs", v2, "\n")
    
    print(
      safe_cor(
        pheno[[v1]],
        pheno[[v2]]
      )
    )
  }
}

cat("\nCorrelations with chronological age:\n")

for (v in age_assoc_vars) {
  
  if (v %in% names(pheno)) {
    
    cat(
      "\nCorrelation: ",
      v,
      " vs chronological age\n",
      sep = ""
    )
    
    print(
      safe_cor(
        pheno[[v]],
        chron_age_years
      )
    )
  }
}

sink()

cat("Saved: epi_age_diagnostics_report.txt\n")

# ------------------------------------------------------------
# 7. Basic diagnostic plots
# ------------------------------------------------------------

cat("\n[17] Creating diagnostic plots...\n")

plot_dir <- "EpiAge_diagnostic_plots"

if (!dir.exists(plot_dir)) {
  dir.create(plot_dir)
}

# ------------------------------------------------------------
# Chronological age vs DNAm age
# ------------------------------------------------------------

png(
  file.path(plot_dir, "chron_age_vs_DNAmAge_Horvath.png"),
  width = 1800,
  height = 1400,
  res = 220
)

plot(
  chron_age_years,
  pheno$DNAmAge_Horvath,
  pch = 16,
  col = rgb(0, 0, 0, 0.35),
  xlab = "Chronological age (years)",
  ylab = "Horvath DNAm age (years)",
  main = "Horvath DNAm age vs chronological age"
)

abline(0, 1, lty = 2, lwd = 2)

dev.off()

# ------------------------------------------------------------

png(
  file.path(plot_dir, "chron_age_vs_DNAmAge_Hannum.png"),
  width = 1800,
  height = 1400,
  res = 220
)

plot(
  chron_age_years,
  pheno$DNAmAge_Hannum,
  pch = 16,
  col = rgb(0, 0, 0, 0.35),
  xlab = "Chronological age (years)",
  ylab = "Hannum DNAm age (years)",
  main = "Hannum DNAm age vs chronological age"
)

abline(0, 1, lty = 2, lwd = 2)

dev.off()

# ------------------------------------------------------------

png(
  file.path(plot_dir, "chron_age_vs_DNAmPhenoAge_Levine.png"),
  width = 1800,
  height = 1400,
  res = 220
)

plot(
  chron_age_years,
  pheno$DNAmPhenoAge_Levine,
  pch = 16,
  col = rgb(0, 0, 0, 0.35),
  xlab = "Chronological age (years)",
  ylab = "Levine / DNAm PhenoAge (years)",
  main = "Levine / DNAm PhenoAge vs chronological age"
)

abline(0, 1, lty = 2, lwd = 2)

dev.off()

# ------------------------------------------------------------
# Residual acceleration by age group
# ------------------------------------------------------------

age_group <- factor(
  pheno$Age_months,
  levels = sort(unique(pheno$Age_months)),
  labels = paste0(sort(unique(pheno$Age_months)), " months")
)

# ------------------------------------------------------------

png(
  file.path(plot_dir, "AgeAccelResidual_Horvath_by_age_group.png"),
  width = 1800,
  height = 1400,
  res = 220
)

boxplot(
  pheno$AgeAccelResidual_Horvath ~ age_group,
  xlab = "Chronological age group",
  ylab = "Horvath residual age acceleration",
  main = "Horvath residual acceleration by age group"
)

abline(h = 0, lty = 2, lwd = 2)

dev.off()

# ------------------------------------------------------------

png(
  file.path(plot_dir, "AgeAccelResidual_Hannum_by_age_group.png"),
  width = 1800,
  height = 1400,
  res = 220
)

boxplot(
  pheno$AgeAccelResidual_Hannum ~ age_group,
  xlab = "Chronological age group",
  ylab = "Hannum residual age acceleration",
  main = "Hannum residual acceleration by age group"
)

abline(h = 0, lty = 2, lwd = 2)

dev.off()

# ------------------------------------------------------------

png(
  file.path(plot_dir, "AgeAccelResidual_Levine_by_age_group.png"),
  width = 1800,
  height = 1400,
  res = 220
)

boxplot(
  pheno$AgeAccelResidual_Levine ~ age_group,
  xlab = "Chronological age group",
  ylab = "Levine residual age acceleration",
  main = "Levine residual acceleration by age group"
)

abline(h = 0, lty = 2, lwd = 2)

dev.off()

# ------------------------------------------------------------
# Cross-clock residual acceleration scatterplots
# ------------------------------------------------------------

png(
  file.path(plot_dir, "cross_clock_residual_Levine_vs_Horvath.png"),
  width = 1800,
  height = 1400,
  res = 220
)

plot(
  pheno$AgeAccelResidual_Levine,
  pheno$AgeAccelResidual_Horvath,
  pch = 16,
  col = rgb(0, 0, 0, 0.35),
  xlab = "Levine residual acceleration",
  ylab = "Horvath residual acceleration",
  main = "Cross-clock residual acceleration: Levine vs Horvath"
)

abline(
  lm(
    pheno$AgeAccelResidual_Horvath ~
      pheno$AgeAccelResidual_Levine
  ),
  lwd = 2
)

dev.off()

# ------------------------------------------------------------

png(
  file.path(plot_dir, "cross_clock_residual_Levine_vs_Hannum.png"),
  width = 1800,
  height = 1400,
  res = 220
)

plot(
  pheno$AgeAccelResidual_Levine,
  pheno$AgeAccelResidual_Hannum,
  pch = 16,
  col = rgb(0, 0, 0, 0.35),
  xlab = "Levine residual acceleration",
  ylab = "Hannum residual acceleration",
  main = "Cross-clock residual acceleration: Levine vs Hannum"
)

abline(
  lm(
    pheno$AgeAccelResidual_Hannum ~
      pheno$AgeAccelResidual_Levine
  ),
  lwd = 2
)

dev.off()

# ------------------------------------------------------------

png(
  file.path(plot_dir, "cross_clock_residual_Horvath_vs_Hannum.png"),
  width = 1800,
  height = 1400,
  res = 220
)

plot(
  pheno$AgeAccelResidual_Horvath,
  pheno$AgeAccelResidual_Hannum,
  pch = 16,
  col = rgb(0, 0, 0, 0.35),
  xlab = "Horvath residual acceleration",
  ylab = "Hannum residual acceleration",
  main = "Cross-clock residual acceleration: Horvath vs Hannum"
)

abline(
  lm(
    pheno$AgeAccelResidual_Hannum ~
      pheno$AgeAccelResidual_Horvath
  ),
  lwd = 2
)

dev.off()

cat("Saved diagnostic plots in folder:", plot_dir, "\n")

cat("\nDiagnostics complete.\n")