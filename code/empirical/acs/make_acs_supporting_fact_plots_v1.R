#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args[grep(file_arg, args)])
if (length(script_path) == 0) {
  stop("Could not determine script path from commandArgs().")
}
script_dir <- dirname(normalizePath(script_path))
root <- normalizePath(file.path(script_dir, "..", ".."))
cross_dest_path <- file.path(root, "output/acs_income_crossmetro_v1/acs_cross_metro_dest_gap_v1.csv")
cross_stock_path <- file.path(root, "output/acs_income_crossmetro_v1/acs_cross_metro_center_gap_v1.csv")
elastic_path <- file.path(root, "output/acs_renter_elasticity_v1/acs_renter_elasticity_summary_v1.csv")
outdir <- file.path(root, "output/acs_supporting_facts_v1")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cross_dest <- read.csv(cross_dest_path)
cross_stock <- read.csv(cross_stock_path)
elastic <- read.csv(elastic_path)

navy <- "#1b4965"
red <- "#ae2012"
quartile_cols <- c("#d8e2eb", "#9bb8cb", "#5f8aa6", "#1b4965")

weighted_mean <- function(x, w) {
  sum(x * w, na.rm = TRUE) / sum(w, na.rm = TRUE)
}

rank_id <- rank(cross_stock$log_rent_gap, ties.method = "first", na.last = "keep")
cross_stock$gap_quartile <- ceiling(4 * rank_id / sum(!is.na(rank_id)))
cross_stock$gap_quartile <- factor(
  cross_stock$gap_quartile,
  levels = 1:4,
  labels = c("Q1", "Q2", "Q3", "Q4")
)

cross_dest <- merge(
  cross_dest,
  cross_stock[, c("met2013", "gap_quartile")],
  by = "met2013",
  all.x = TRUE,
  sort = FALSE
)

quartile_levels <- levels(cross_stock$gap_quartile)

stock_plot <- do.call(
  rbind,
  lapply(quartile_levels, function(lbl) {
    df <- cross_stock[cross_stock$gap_quartile == lbl, ]
    data.frame(
      gap_quartile = lbl,
      series = "Stock center share",
      gap_pp = weighted_mean(df$gap_pp, df$metro_obs)
    )
  })
)

dest_plot <- do.call(
  rbind,
  lapply(quartile_levels, function(lbl) {
    df <- cross_dest[cross_dest$gap_quartile == lbl, ]
    data.frame(
      gap_quartile = lbl,
      series = "Mover center destination",
      gap_pp = weighted_mean(df$gap_pp, df$mover_obs)
    )
  })
)

quartile_plot <- rbind(stock_plot, dest_plot)
quartile_plot$gap_quartile <- factor(quartile_plot$gap_quartile, levels = quartile_levels)
quartile_plot$series <- factor(
  quartile_plot$series,
  levels = c("Stock center share", "Mover center destination")
)

p1 <- ggplot(quartile_plot, aes(x = gap_quartile, y = gap_pp, color = series, group = series)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 3) +
  geom_text(aes(label = sprintf("%.1f", gap_pp)), vjust = -0.7, size = 3.5, show.legend = FALSE) +
  scale_color_manual(values = c("Stock center share" = navy, "Mover center destination" = red)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
  labs(
    x = "Metro rent-gap quartile (Q1 low, Q4 high)",
    y = "Non-parent \u2212 new-parent center share (pp)",
    color = ""
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.title.x = element_text(size = 11),
    axis.title.y = element_text(size = 11),
    legend.position = "top"
  )

ggsave(
  filename = file.path(outdir, "acs_cross_metro_dest_gap_quartiles_v1.png"),
  plot = p1,
  width = 8,
  height = 5,
  dpi = 150
)

rooms_row <- elastic[elastic$outcome == "ln_rooms", ]
elastic_plot <- data.frame(
  group = factor(c("Childless", "Parents"), levels = c("Childless", "Parents")),
  elasticity = c(rooms_row$childless_slope, rooms_row$parent_slope),
  se = c(rooms_row$childless_slope_se, rooms_row$parent_slope_se),
  fill = c(navy, red)
)

p2 <- ggplot(elastic_plot, aes(x = group, y = elasticity, fill = group)) +
  geom_col(width = 0.62, show.legend = FALSE) +
  geom_errorbar(
    aes(ymin = elasticity - 1.96 * se, ymax = elasticity + 1.96 * se),
    width = 0.12,
    linewidth = 0.7
  ) +
  geom_text(aes(label = sprintf("%.3f", elasticity)), vjust = -0.35, size = 4) +
  scale_fill_manual(values = c("Childless" = navy, "Parents" = red)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    x = "",
    y = "Income elasticity of rooms"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.title.y = element_text(size = 11)
  )

ggsave(
  filename = file.path(outdir, "acs_renter_rooms_elasticity_bars_v1.png"),
  plot = p2,
  width = 8,
  height = 5,
  dpi = 150
)

cat(file.path(outdir, "acs_cross_metro_dest_gap_quartiles_v1.png"), "\n")
cat(file.path(outdir, "acs_renter_rooms_elasticity_bars_v1.png"), "\n")
