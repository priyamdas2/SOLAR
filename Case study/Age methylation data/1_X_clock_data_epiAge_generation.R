# Runtime ~ 7-8 minutes
setwd("U:/SOLAR/Case study/Age methylation data")

library(hdf5r)

# ------------------------------------------------------------
# Read CpG names
# ------------------------------------------------------------

X_cpg_names <- readRDS("X_cpg_names.rds")

head(X_cpg_names)
length(X_cpg_names)

# ------------------------------------------------------------
# Read MATLAB v7.3 matrix using HDF5
# ------------------------------------------------------------

cat("Opening MATLAB v7.3 .mat file...\n")

t0 <- Sys.time()

h5 <- H5File$new("X_matrix_final.mat", mode = "r")

# Show object names inside the .mat file
print(h5$ls())

# ------------------------------------------------------------
# Extract matrix
# ------------------------------------------------------------
# Usually there is only one main matrix object.
# This automatically takes the first object name.

obj_name <- h5$ls()$name[1]

cat("Reading object:", obj_name, "\n")

ds <- h5[[obj_name]]
X_values <- ds[1:ds$dims[1], 1:ds$dims[2]]

h5$close_all()

cat("Time to read .mat file:\n")
print(Sys.time() - t0)

dim(X_values)

# ------------------------------------------------------------
# Fix orientation if needed
# Desired format here:
# rows = samples
# cols = CpGs
# expected: 1051 x 865859
# ------------------------------------------------------------

if (nrow(X_values) == length(X_cpg_names)) {
  cat("Detected CpGs as rows; transposing to samples x CpGs...\n")
  X_values <- t(X_values)
}

dim(X_values)

# ------------------------------------------------------------
# Check dimensions
# ------------------------------------------------------------

if (ncol(X_values) != length(X_cpg_names)) {
  stop(
    paste0(
      "Mismatch after orientation check: X has ", ncol(X_values),
      " CpG columns, but CpG-name vector has ",
      length(X_cpg_names), " entries."
    )
  )
}

# ------------------------------------------------------------
# Reduce precision to 4 decimals
# ------------------------------------------------------------

cat("Rounding X values to 4 decimals...\n")
X_values <- round(X_values, 4)

# ------------------------------------------------------------
# Attach CpG names
# Current format:
# rows = samples
# cols = CpGs
# ------------------------------------------------------------

colnames(X_values) <- X_cpg_names

# ------------------------------------------------------------
# Convert to clock-ready format
# Required later by epigenetic age packages:
# rows = CpGs
# cols = samples
# ------------------------------------------------------------

cat("Transposing to clock-ready format: CpGs x samples...\n")

X_clock <- t(X_values)

colnames(X_clock) <- paste0("Sample_", seq_len(ncol(X_clock)))

# Quick checks
dim(X_clock)
rownames(X_clock)[1:5]
X_clock[1:3, 1:3]


# ------------------------------------------------------------
# Calculate epigenetic ages and save augmented phenotype file
# ------------------------------------------------------------

cat("Reading phenotype file...\n")

pheno <- read.csv(
  "Pheno_for_MATLAB.csv",
  check.names = FALSE
)

if (nrow(pheno) != ncol(X_clock)) {
  stop(
    paste0(
      "Mismatch: phenotype has ", nrow(pheno),
      " rows, but X_clock has ", ncol(X_clock),
      " samples."
    )
  )
}

# Chronological age
# Your MATLAB code used Age_months, so convert to years for methylclock.
chron_age_years <- pheno$Age_months / 12

# Install methylclock if needed:
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("methylclock", quietly = TRUE)) BiocManager::install("methylclock", force = TRUE)

library(methylclock)

cat("Calculating DNAm ages using methylclock...\n")

t_epi <- Sys.time()

epi_age <- DNAmAge(
  x = X_clock,
  clocks = c("Horvath", "Hannum", "Levine"),
  age = chron_age_years,
  toBetas = FALSE,
  normalize = FALSE,
  fastImp = TRUE,
  min.perc = 0.8
)

cat("Time to calculate epigenetic ages:\n")
print(Sys.time() - t_epi)

epi_age <- as.data.frame(epi_age)

# Check output columns
print(names(epi_age))

# ------------------------------------------------------------
# Add selected epigenetic age columns to phenotype file
# ------------------------------------------------------------

# Horvath
if ("Horvath" %in% names(epi_age)) {
  pheno$DNAmAge_Horvath <- epi_age$Horvath
  pheno$AgeAccelDiff_Horvath <- epi_age$Horvath - chron_age_years
  pheno$AgeAccelResidual_Horvath <- resid(lm(epi_age$Horvath ~ chron_age_years))
}

# Hannum
if ("Hannum" %in% names(epi_age)) {
  pheno$DNAmAge_Hannum <- epi_age$Hannum
  pheno$AgeAccelDiff_Hannum <- epi_age$Hannum - chron_age_years
  pheno$AgeAccelResidual_Hannum <- resid(lm(epi_age$Hannum ~ chron_age_years))
}

# PhenoAge / Levine
if ("Levine" %in% names(epi_age)) {
  pheno$DNAmPhenoAge_Levine <- epi_age$Levine
  pheno$AgeAccelDiff_Levine <- epi_age$Levine - chron_age_years
  pheno$AgeAccelResidual_Levine <- resid(lm(epi_age$Levine ~ chron_age_years))
}

# Also keep any methylclock-provided acceleration columns if available
extra_cols <- setdiff(names(epi_age), names(pheno))
pheno <- cbind(pheno, epi_age[, extra_cols, drop = FALSE])

# ------------------------------------------------------------
# Save augmented phenotype file
# ------------------------------------------------------------

write.csv(
  pheno,
  file = "Pheno_for_MATLAB_with_epi_ages.csv",
  row.names = FALSE
)

cat("Saved phenotype file with epigenetic ages:\n")
cat("Pheno_for_MATLAB_with_epi_ages.csv\n")



# ------------------------------------------------------------
# Save compact R object  ~ 5 minutes
# ------------------------------------------------------------

cat("Saving clock-ready RDS object...\n")

saveRDS(
  X_clock,
  file = "X_clock_CpGs_by_samples_4dp.rds",
  compress = FALSE
)

cat("Saved:\n")
cat("X_clock_CpGs_by_samples_4dp.rds\n")