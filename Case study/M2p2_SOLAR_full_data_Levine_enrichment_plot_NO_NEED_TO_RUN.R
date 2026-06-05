# ============================================================
# CpG annotation / enrichment plot for leading SOLAR CpGs
# Main-paper enrichment plot using top 1000 CpGs
# extracted from the saved top-2000 SOLAR importance file
# ============================================================

setwd("U:/SOLAR/Case study")
set.seed(123)

out_dir  <- "Output"
data_dir <- "Age methylation data"

if (!dir.exists(out_dir)) dir.create(out_dir)

# ------------------------------------------------------------
# Packages
# ------------------------------------------------------------

cran_pkgs <- c(
  "dplyr", "tidyr", "ggplot2", "readr",
  "forcats", "stringr", "scales", "tibble"
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
  "clusterProfiler",
  "org.Hs.eg.db",
  "IlluminaHumanMethylationEPICanno.ilm10b4.hg19"
)

for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(forcats)
library(stringr)
library(scales)
library(tibble)
library(clusterProfiler)
library(org.Hs.eg.db)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# ------------------------------------------------------------
# Inputs
# ------------------------------------------------------------

importance_file <- file.path(
  out_dir,
  "Output_from_M2_top_CpG_importance_for_enrichment_2000.csv"
)

cpg_names_file <- file.path(data_dir, "X_cpg_names.rds")

if (!file.exists(importance_file)) {
  stop("Cannot find: ", importance_file)
}

if (!file.exists(cpg_names_file)) {
  stop("Cannot find: ", cpg_names_file)
}

top_cpg_tbl <- readr::read_csv(
  importance_file,
  show_col_types = FALSE
)

x_cpg_names <- readRDS(cpg_names_file)

# ------------------------------------------------------------
# Attach CpG names if needed
# ------------------------------------------------------------

if (!"CpG" %in% names(top_cpg_tbl)) {
  
  if (!"CpG_index" %in% names(top_cpg_tbl)) {
    stop("Need either CpG or CpG_index column in importance file.")
  }
  
  top_cpg_tbl <- top_cpg_tbl %>%
    dplyr::mutate(CpG = x_cpg_names[CpG_index])
}

# ------------------------------------------------------------
# Select leading CpGs for main enrichment analysis
# ------------------------------------------------------------

top_n_enrich <- 1000

leading_cpgs <- top_cpg_tbl %>%
  dplyr::arrange(rank) %>%
  dplyr::filter(!is.na(CpG)) %>%
  dplyr::slice_head(n = top_n_enrich) %>%
  dplyr::pull(CpG) %>%
  unique()

# ------------------------------------------------------------
# Load Illumina EPIC CpG annotation
# ------------------------------------------------------------

anno <- getAnnotation(
  IlluminaHumanMethylationEPICanno.ilm10b4.hg19
)

anno_tbl <- as.data.frame(anno) %>%
  tibble::rownames_to_column("CpG")

if (!"UCSC_RefGene_Name" %in% names(anno_tbl)) {
  stop("UCSC_RefGene_Name column not found in EPIC annotation table.")
}

message("Annotation source: IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
message("CpG column used: CpG")
message("Gene column used: UCSC_RefGene_Name")

# ------------------------------------------------------------
# Map leading CpGs to annotated genes
# One CpG can map to multiple genes.
# ------------------------------------------------------------

cpg_gene_tbl <- anno_tbl %>%
  dplyr::transmute(
    CpG = as.character(CpG),
    gene_raw = as.character(UCSC_RefGene_Name)
  ) %>%
  dplyr::filter(CpG %in% leading_cpgs) %>%
  dplyr::filter(!is.na(gene_raw), gene_raw != "") %>%
  tidyr::separate_rows(gene_raw, sep = ";") %>%
  dplyr::mutate(
    gene_symbol = stringr::str_trim(gene_raw)
  ) %>%
  dplyr::filter(gene_symbol != "") %>%
  dplyr::distinct(CpG, gene_symbol)

gene_list <- unique(cpg_gene_tbl$gene_symbol)

message("Leading CpGs selected: ", length(leading_cpgs))
message("Leading CpGs mapped to at least one gene: ", length(unique(cpg_gene_tbl$CpG)))
message("Unique annotated genes used for enrichment: ", length(gene_list))

if (length(gene_list) < 10) {
  stop("Too few annotated genes for enrichment analysis.")
}

# ------------------------------------------------------------
# GO biological-process enrichment
# ------------------------------------------------------------

ego <- clusterProfiler::enrichGO(
  gene          = gene_list,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.20,
  readable      = TRUE
)

ego_tbl <- as.data.frame(ego)

if (nrow(ego_tbl) == 0) {
  
  message("No GO terms passed the adjusted enrichment cutoffs for Levine.")
  message("Re-running enrichment without strict cutoffs to save exploratory results.")
  
  ego <- clusterProfiler::enrichGO(
    gene          = gene_list,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 1,
    qvalueCutoff  = 1,
    readable      = TRUE
  )
  
  ego_tbl <- as.data.frame(ego)
  
  readr::write_csv(
    ego_tbl,
    file.path(out_dir, "Output_from_M2_CpG_GO_enrichment_top1000_exploratory_all_terms.csv")
  )
  
  stop("No adjusted significant GO BP terms found for Levine; exploratory all-term table saved.")
}

# ------------------------------------------------------------
# Prepare plot table
# ------------------------------------------------------------

plot_n_terms <- 10

plot_tbl <- ego_tbl %>%
  dplyr::mutate(
    neglog10_padj = -log10(p.adjust),
    GeneRatio_num = sapply(GeneRatio, function(x) {
      xx <- as.numeric(strsplit(x, "/")[[1]])
      xx[1] / xx[2]
    }),
    term_label = stringr::str_wrap(Description, width = 28)
  ) %>%
  dplyr::arrange(p.adjust) %>%
  dplyr::slice_head(n = plot_n_terms) %>%
  dplyr::mutate(
    term_label = forcats::fct_reorder(term_label, neglog10_padj)
  )

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------

x_min <- max(0, floor(min(plot_tbl$neglog10_padj) * 10) / 10 - 0.05)
x_max <- ceiling(max(plot_tbl$neglog10_padj) * 10) / 10 + 0.05

p_enrich <- ggplot(
  plot_tbl,
  aes(
    x = neglog10_padj,
    y = term_label,
    size = Count,
    color = GeneRatio_num
  )
) +
  geom_point(alpha = 0.92) +
  geom_vline(
    xintercept = -log10(0.05),
    linetype = "dashed",
    linewidth = 0.45,
    color = "gray55"
  ) +
  scale_x_continuous(
    limits = c(x_min, x_max),
    breaks = pretty_breaks(n = 4),
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  scale_color_viridis_c(
    name = "Gene ratio",
    option = "viridis"
  ) +
  scale_size_continuous(
    name = "Mapped genes",
    range = c(2.2, 8.5)
  ) +
  labs(
    title = "GO enrichment of leading SOLAR CpGs (Levine)",
    x = expression("Enrichment strength, " * -log[10] * "(adjusted " * italic(p) * ")"),
    y = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11, lineheight = 0.92),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    legend.position = "right",
    legend.box = "vertical",
    legend.box.just = "center",
    legend.spacing.y = unit(0.35, "cm"),
    legend.background = element_rect(color = "gray80", fill = "white", linewidth = 0.35),
    legend.margin = margin(6, 6, 6, 6),
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.35)
  ) +
  guides(
    size = guide_legend(
      order = 1,
      title.position = "top",
      title.hjust = 0.5,
      label.position = "right"
    ),
    color = guide_colorbar(
      order = 2,
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(0.45, "cm"),
      barheight = unit(2.2, "cm")
    )
  )

p_enrich

ggsave(
  filename = file.path(out_dir, "Plot_M2_CpG_enrichment_dotplot_top1000.png"),
  plot = p_enrich,
  width = 8.5,
  height = 4.8,
  dpi = 450,
  bg = "white"
)

# ------------------------------------------------------------
# Save enrichment table and CpG-gene mapping
# ------------------------------------------------------------

readr::write_csv(
  ego_tbl,
  file.path(out_dir, "Output_from_M2_CpG_GO_enrichment_top1000.csv")
)

readr::write_csv(
  cpg_gene_tbl,
  file.path(out_dir, "Output_from_M2_CpG_gene_mapping_top1000.csv")
)

message("Saved enrichment plot and tables for top 1000 CpGs.")