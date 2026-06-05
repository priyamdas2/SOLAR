# ============================================================
# Create supplementary table of leading SOLAR CpGs
# Columns:
# CpG ID
# Relative importance
# Dominant factor
# Dominant loading
# Annotated gene(s)
# ============================================================

setwd("U:/SOLAR/Case study")

out_dir  <- "Output"
data_dir <- "Age methylation data"

# ------------------------------------------------------------
# Packages
# ------------------------------------------------------------

cran_pkgs <- c(
  "dplyr",
  "readr",
  "stringr",
  "tibble",
  "tidyr"
)

for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_pkgs <- c(
  "IlluminaHumanMethylationEPICanno.ilm10b4.hg19"
)

for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

library(dplyr)
library(readr)
library(stringr)
library(tibble)
library(tidyr)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# ------------------------------------------------------------
# Input files
# ------------------------------------------------------------

importance_file <- file.path(
  out_dir,
  "Output_from_M1_top_CpG_importance.csv"
)

cpg_names_file <- file.path(
  data_dir,
  "X_cpg_names.rds"
)

if (!file.exists(importance_file)) {
  stop("Cannot find: ", importance_file)
}

if (!file.exists(cpg_names_file)) {
  stop("Cannot find: ", cpg_names_file)
}

# ------------------------------------------------------------
# Read data
# ------------------------------------------------------------

top_tbl <- readr::read_csv(
  importance_file,
  show_col_types = FALSE
)

x_cpg_names <- readRDS(cpg_names_file)

# ------------------------------------------------------------
# Attach CpG IDs if missing
# ------------------------------------------------------------

if (!"CpG" %in% names(top_tbl)) {
  
  if (!"CpG_index" %in% names(top_tbl)) {
    stop("Need CpG_index column.")
  }
  
  top_tbl <- top_tbl %>%
    dplyr::slice_head(n = 50)
  
  top_tbl <- top_tbl %>%
    dplyr::mutate(
      CpG = x_cpg_names[CpG_index]
    )
  
}

# ------------------------------------------------------------
# Identify loading columns
# ------------------------------------------------------------

loading_cols <- grep(
  "^loading_V_factor_",
  names(top_tbl),
  value = TRUE
)

if (length(loading_cols) == 0) {
  stop("No loading_V_factor columns found.")
}

# ------------------------------------------------------------
# Compute dominant factor and dominant loading
# ------------------------------------------------------------

loading_mat <- as.matrix(top_tbl[, loading_cols])

dominant_factor_idx <- apply(
  abs(loading_mat),
  1,
  which.max
)

dominant_loading <- sapply(
  seq_len(nrow(loading_mat)),
  function(i) {
    loading_mat[i, dominant_factor_idx[i]]
  }
)

dominant_factor <- paste0(
  "Factor ",
  dominant_factor_idx
)

# ------------------------------------------------------------
# Load Illumina EPIC annotation
# ------------------------------------------------------------

anno <- getAnnotation(
  IlluminaHumanMethylationEPICanno.ilm10b4.hg19
)

anno_tbl <- as.data.frame(anno) %>%
  tibble::rownames_to_column("CpG")

if (!"UCSC_RefGene_Name" %in% names(anno_tbl)) {
  stop("UCSC_RefGene_Name column not found.")
}

# ------------------------------------------------------------
# Create CpG-to-gene mapping
# ------------------------------------------------------------

gene_tbl <- anno_tbl %>%
  dplyr::transmute(
    CpG = as.character(CpG),
    gene_raw = as.character(UCSC_RefGene_Name)
  ) %>%
  dplyr::mutate(
    gene_raw = ifelse(
      is.na(gene_raw) | gene_raw == "",
      NA,
      gene_raw
    )
  ) %>%
  tidyr::separate_rows(gene_raw, sep = ";") %>%
  dplyr::mutate(
    gene_symbol = stringr::str_trim(gene_raw)
  ) %>%
  dplyr::filter(
    !is.na(gene_symbol),
    gene_symbol != ""
  ) %>%
  dplyr::group_by(CpG) %>%
  dplyr::summarise(
    annotated_genes = paste(
      unique(gene_symbol),
      collapse = "; "
    ),
    .groups = "drop"
  )

# ------------------------------------------------------------
# Create final supplementary table
# ------------------------------------------------------------



supp_tbl <- top_tbl %>%
  dplyr::mutate(
    Rank = dplyr::row_number(),
    relative_importance = importance / max(importance, na.rm = TRUE),
    Dominant_factor = dominant_factor,
    Dominant_loading = dominant_loading
  ) %>%
  dplyr::left_join(
    gene_tbl,
    by = "CpG"
  ) %>%
  dplyr::transmute(
    Rank,
    CpG_ID = CpG,
    Relative_importance = round(relative_importance, 3),
    Dominant_factor,
    Dominant_loading = round(Dominant_loading, 5),
    Annotated_genes = ifelse(
      is.na(annotated_genes),
      "--",
      annotated_genes
    )
  )



# ------------------------------------------------------------
# Save output
# ------------------------------------------------------------

output_file <- file.path(
  out_dir,
  "Output_from_M1p3_Rcode_top_CpG_summary_table.csv"
)

readr::write_csv(
  supp_tbl,
  output_file
)

message("Saved supplementary CpG summary table:")
message(output_file)