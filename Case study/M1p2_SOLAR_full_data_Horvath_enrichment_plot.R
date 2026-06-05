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
  "Output_from_M1_top_CpG_importance_for_enrichment_2000.csv"
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
  stop("No enriched GO terms found under current cutoffs.")
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
    title = "GO enrichment of leading SOLAR CpGs (Horvath)",
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
  filename = file.path(out_dir, "Plot_M1_CpG_enrichment_dotplot_top1000.png"),
  plot = p_enrich,
  width = 7.5, # 7
  height = 4.8,
  dpi = 450,
  bg = "white"
)

# ------------------------------------------------------------
# Save enrichment table and CpG-gene mapping
# ------------------------------------------------------------

readr::write_csv(
  ego_tbl,
  file.path(out_dir, "Output_from_M1_CpG_GO_enrichment_top1000.csv")
)

readr::write_csv(
  cpg_gene_tbl,
  file.path(out_dir, "Output_from_M1_CpG_gene_mapping_top1000.csv")
)

message("Saved enrichment plot and tables for top 1000 CpGs.")
#############################################################################

# ------------------------------------------------------------
# Additional supplementary GO table for transparency
# Safe to run after the existing code above.
# Does not change enrichment analysis or plot.
# ------------------------------------------------------------

parse_ratio <- function(x) {
  sapply(x, function(z) {
    zz <- as.numeric(strsplit(as.character(z), "/")[[1]])
    zz[1] / zz[2]
  })
}

ego_tbl_supp <- ego_tbl %>%
  dplyr::mutate(
    GeneRatio_num = parse_ratio(GeneRatio),
    BgRatio_num = parse_ratio(BgRatio),
    EnrichmentRatio = GeneRatio_num / BgRatio_num
  )

# Add CpG-level mapping for genes appearing in each GO term
go_supp_single_csv <- ego_tbl_supp %>%
  dplyr::select(
    GO_ID = ID,
    GO_term = Description,
    GeneRatio,
    BgRatio,
    EnrichmentRatio,
    pvalue,
    p.adjust,
    qvalue,
    Count,
    geneID
  ) %>%
  tidyr::separate_rows(geneID, sep = "/") %>%
  dplyr::rename(gene_symbol = geneID) %>%
  dplyr::left_join(cpg_gene_tbl, by = "gene_symbol") %>%
  dplyr::left_join(
    top_cpg_tbl %>%
      dplyr::mutate(CpG = as.character(CpG)) %>%
      dplyr::select(any_of(c("CpG", "CpG_index", "rank", "RIS", "Imp", "Importance"))),
    by = "CpG"
  ) %>%
  dplyr::group_by(
    GO_ID,
    GO_term,
    GeneRatio,
    BgRatio,
    EnrichmentRatio,
    pvalue,
    p.adjust,
    qvalue,
    Count
  ) %>%
  dplyr::summarise(
    unique_genes = dplyr::n_distinct(gene_symbol, na.rm = TRUE),
    unique_CpGs = dplyr::n_distinct(CpG, na.rm = TRUE),
    gene_symbols = paste(sort(unique(gene_symbol[!is.na(gene_symbol)])), collapse = ";"),
    CpGs = paste(sort(unique(CpG[!is.na(CpG)])), collapse = ";"),
    CpG_gene_pairs = paste(
      sort(unique(paste0(CpG[!is.na(CpG)], ":", gene_symbol[!is.na(gene_symbol)]))),
      collapse = ";"
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    CpG_to_gene_strategy = "CpGs were mapped to genes using UCSC_RefGene_Name from IlluminaHumanMethylationEPICanno.ilm10b4.hg19; one CpG may map to multiple genes.",
    GO_background = "Default clusterProfiler::enrichGO background used in the current analysis.",
    CpG_set = paste0("Top ", top_n_enrich, " CpGs from the SOLAR supervised importance ranking."),
    Factor_note = "SOLAR supervised importance scores aggregate contributions across the inferred latent factors."
  ) %>%
  dplyr::arrange(p.adjust, GO_ID)

readr::write_csv(
  go_supp_single_csv,
  file.path(out_dir, "Supplement_GO_enrichment_CpG_gene_mapping_top1000.csv")
)

message("Saved single supplementary GO table: Supplement_GO_enrichment_CpG_gene_mapping_top1000.csv")

# ------------------------------------------------------------
# Supplementary GO table
# One row per GO term with GO ID, enrichment ratio,
# gene names, and corresponding CpGs.
# Safe to run after the existing code.
# ------------------------------------------------------------

parse_ratio <- function(x) {
  sapply(x, function(z) {
    zz <- as.numeric(strsplit(as.character(z), "/")[[1]])
    zz[1] / zz[2]
  })
}

# Use the same top GO terms shown in the main plot
# Change to Inf if you want all enriched GO terms in the supplement.
supp_n_terms <- plot_n_terms

ego_supp_base <- ego_tbl %>%
  dplyr::mutate(
    GeneRatio_num = parse_ratio(GeneRatio),
    BgRatio_num   = parse_ratio(BgRatio),
    EnrichmentRatio = GeneRatio_num / BgRatio_num
  ) %>%
  dplyr::arrange(p.adjust) %>%
  dplyr::slice_head(n = supp_n_terms)

# Create gene-to-CpG mapping strings within each GO term
go_gene_cpg_pairs <- ego_supp_base %>%
  dplyr::select(
    GO_ID = ID,
    GO_term = Description,
    geneID
  ) %>%
  tidyr::separate_rows(geneID, sep = "/") %>%
  dplyr::rename(gene_symbol = geneID) %>%
  dplyr::left_join(cpg_gene_tbl, by = "gene_symbol") %>%
  dplyr::filter(!is.na(gene_symbol), !is.na(CpG)) %>%
  dplyr::group_by(GO_ID, GO_term, gene_symbol) %>%
  dplyr::summarise(
    CpGs_for_gene = paste(sort(unique(CpG)), collapse = ", "),
    n_CpGs_for_gene = dplyr::n_distinct(CpG),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    Gene_CpG_mapping = paste0(gene_symbol, " [", CpGs_for_gene, "]")
  )

go_gene_cpg_summary <- go_gene_cpg_pairs %>%
  dplyr::group_by(GO_ID, GO_term) %>%
  dplyr::summarise(
    n_unique_genes = dplyr::n_distinct(gene_symbol),
    n_unique_CpGs = sum(n_CpGs_for_gene),
    gene_names = paste(sort(unique(gene_symbol)), collapse = "; "),
    genes_with_corresponding_CpGs = paste(Gene_CpG_mapping, collapse = "; "),
    .groups = "drop"
  )

# Final paper-ready supplementary table
go_supp_paper_ready <- ego_supp_base %>%
  dplyr::select(
    GO_ID = ID,
    GO_term = Description,
    GeneRatio,
    BgRatio,
    EnrichmentRatio,
    Count,
    pvalue,
    p_adjust = p.adjust,
    qvalue
  ) %>%
  dplyr::left_join(
    go_gene_cpg_summary,
    by = c("GO_ID", "GO_term")
  ) %>%
  dplyr::mutate(
    EnrichmentRatio = round(EnrichmentRatio, 3),
    pvalue = signif(pvalue, 3),
    p_adjust = signif(p_adjust, 3),
    qvalue = signif(qvalue, 3)
  ) %>%
  dplyr::arrange(p_adjust, GO_ID)

readr::write_csv(
  go_supp_paper_ready,
  file.path(out_dir, "Supplement_GO_enrichment_table_top1000_CpG_gene_mapping.csv")
)

message("Saved paper-ready supplementary GO table: Supplement_GO_enrichment_table_top1000_CpG_gene_mapping.csv")