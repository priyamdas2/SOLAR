# ============================================================
# Bootstrap selected-rank distribution plot: Horvath
# ============================================================

setwd("U:/SOLAR/Case study")

out_dir <- "Output"

pkg_needed <- c("ggplot2", "dplyr", "readr", "scales")

for (pkg in pkg_needed) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(ggplot2)
library(dplyr)
library(readr)
library(scales)

# ------------------------------------------------------------
# Read bootstrap rank output
# ------------------------------------------------------------

rank_file <- file.path(
  out_dir,
  "Output_from_M1boot_onlyRanks_NumBoot_500.csv"
)

if (!file.exists(rank_file)) {
  stop("Cannot find: ", rank_file)
}

rank_tbl <- readr::read_csv(
  rank_file,
  show_col_types = FALSE
)

# ------------------------------------------------------------
# Prepare percentage distribution
# ------------------------------------------------------------

rank_dist <- rank_tbl %>%
  dplyr::filter(!is.na(selected_rank)) %>%
  dplyr::count(selected_rank, name = "n") %>%
  dplyr::mutate(
    percent = n / sum(n),
    selected_rank = factor(
      selected_rank,
      levels = sort(unique(selected_rank))
    )
  )

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------

col_horvath <- "#2CA25F"

p_rank <- ggplot(
  rank_dist,
  aes(
    x = selected_rank,
    y = percent
  )
) +
  geom_col(
    width = 0.68,
    fill = col_horvath,
    alpha = 0.85
  ) +
  geom_text(
    aes(label = scales::percent(percent, accuracy = 1)),
    vjust = -0.35,
    size = 4.1,
    color = "black"
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, max(rank_dist$percent) * 1.18),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title = "Bootstrap rank stability (Horvath)",
    x = "Selected latent rank",
    y = "Bootstrap frequency"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 16,
      hjust = 0.5
    ),
    axis.title = element_text(size = 13.5),
    axis.text = element_text(size = 12, color = "black"),
    axis.line = element_line(linewidth = 0.45),
    axis.ticks = element_line(linewidth = 0.45),
    panel.grid.major.y = element_line(
      color = "gray90",
      linewidth = 0.35
    ),
    plot.margin = margin(10, 12, 10, 12)
  )

print(p_rank)

ggsave(
  filename = file.path(out_dir, "Plot_M1boot_rank_distribution.png"),
  plot = p_rank,
  width = 6.5,
  height = 4.8,
  dpi = 450,
  bg = "white"
)

# # Optional: save summarized distribution
# readr::write_csv(
#   rank_dist,
#   file.path(out_dir, "Output_from_M1boot_rank_distribution_summary.csv")
# )