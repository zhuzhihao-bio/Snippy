#!/usr/bin/env Rscript
# =============================================================================
# plot_heatmap.R
# Publication-style core-SNP distance heatmap for Enterococcus faecalis
# host 22M21 (Snippy core.snpdists.tsv).
#
# Author : Zhu Zhihao (朱志豪)
# Date   : 2026-10-13
# 中文说明: 读取 snp-dists 矩阵，将 Reference 重命名为 X019，按克隆群排序，
#           用 sqrt 色阶绘制核 SNP 距离热图，并在格子内标注原始 SNP 数。
#
# Inputs : core.snpdists.tsv (result/ or cwd or ../result/)
#          metadata.txt or metadata (optional; clone groups are hardcoded)
# Outputs: figures/heatmap.pdf, figures/heatmap.png
#          (or <outdir>/ if given as first CLI arg)
# Panel size: compact for main-figure Example strip (~3.1 x 3.6 in; labels only for SNPs <= 100)
#
# Packages: ggplot2, scales
#   # install.packages(c("ggplot2", "scales"))
#   # install.packages(c("pheatmap"))
#   # if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#   # BiocManager::install("ComplexHeatmap")
# =============================================================================

message("=== plot_heatmap.R: Core SNP distance heatmap (host 22M21) ===")

# --- Packages ----------------------------------------------------------------
need <- c("ggplot2", "scales")
missing <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0L) {
  stop(
    "Missing packages: ", paste(missing, collapse = ", "),
    "\nInstall with: install.packages(c(\"ggplot2\", \"scales\"))",
    call. = FALSE
  )
}
library(ggplot2)
library(scales)

# --- Helpers -----------------------------------------------------------------
find_file <- function(candidates) {
  for (p in candidates) {
    if (file.exists(p) && !dir.exists(p)) {
      message("  Found: ", p)
      return(normalizePath(p, winslash = "/", mustWork = TRUE))
    }
  }
  NULL
}

# Clone groups (fixed order so within-/between-clone blocks are visible)
clone_levels <- c("near-ref", "Cluster A", "Cluster B")
clone_order <- c(
  "X019", "Y1620", "Y1621", "Y1624", "Y1626",  # near-ref / clone with X019
  "Y093", "Y097", "Y101", "Y104",              # Cluster A
  "Y1464", "Y1615"                             # Cluster B
)
clone_map <- c(
  X019 = "near-ref", Y1620 = "near-ref", Y1621 = "near-ref",
  Y1624 = "near-ref", Y1626 = "near-ref",
  Y093 = "Cluster A", Y097 = "Cluster A", Y101 = "Cluster A", Y104 = "Cluster A",
  Y1464 = "Cluster B", Y1615 = "Cluster B"
)
# Okabe-Ito colorblind-friendly palette
clone_colors <- c(
  "near-ref"  = "#0072B2",
  "Cluster A" = "#E69F00",
  "Cluster B" = "#009E73"
)

# --- CLI / paths -------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1L) args[[1]] else "figures"

message("Searching for snp-dists matrix ...")
mat_path <- find_file(c(
  "result/core.snpdists.tsv",
  "core.snpdists.tsv",
  "../result/core.snpdists.tsv"
))
if (is.null(mat_path)) {
  stop(
    "Cannot find core.snpdists.tsv. Tried:\n",
    "  result/core.snpdists.tsv\n",
    "  core.snpdists.tsv\n",
    "  ../result/core.snpdists.tsv\n",
    "Set working directory to the github folder or the example repo root.",
    call. = FALSE
  )
}

message("Searching for metadata (optional) ...")
meta_path <- find_file(c(
  "metadata.txt",
  "metadata",
  "result/../metadata.txt",
  "../metadata.txt",
  "../metadata"
))
if (!is.null(meta_path)) {
  message("  Metadata available at: ", meta_path)
} else {
  message("  Metadata not found; using hardcoded clone groups.")
}

# --- Read matrix -------------------------------------------------------------
message("Reading SNP distance matrix ...")
raw <- utils::read.delim(
  mat_path,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)

# First column = rownames; header row starts with "snp-dists ..." then sample names
rn <- as.character(raw[[1]])
mat <- as.matrix(raw[, -1, drop = FALSE])
storage.mode(mat) <- "numeric"
rownames(mat) <- rn
colnames(mat) <- colnames(raw)[-1]

# Rename Reference -> X019 on rows and columns (Snippy reference isolate)
rename_ref <- function(x) {
  x[x == "Reference"] <- "X019"
  x
}
rownames(mat) <- rename_ref(rownames(mat))
colnames(mat) <- rename_ref(colnames(mat))
message("  Renamed 'Reference' -> 'X019' (display name of Snippy reference).")

# Order by clone groups (NOT alphabetical)
keep <- clone_order[clone_order %in% rownames(mat) & clone_order %in% colnames(mat)]
extra <- unique(c(
  setdiff(rownames(mat), clone_order),
  setdiff(colnames(mat), clone_order)
))
if (length(extra)) {
  message("  Samples not in predefined order appended: ", paste(extra, collapse = ", "))
  keep <- c(keep, extra)
}
mat <- mat[keep, keep, drop = FALSE]
n <- nrow(mat)
message("  Matrix size: ", n, " x ", n)

# --- Long format -------------------------------------------------------------
ylabels <- rownames(mat)
xlabels <- colnames(mat)

df <- expand.grid(
  SampleY = factor(ylabels, levels = rev(ylabels)),
  SampleX = factor(xlabels, levels = xlabels),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
df$SNPs <- as.vector(mat)

map_clone <- function(ids) {
  out <- unname(clone_map[as.character(ids)])
  factor(out, levels = clone_levels)
}
df$CloneY <- map_clone(df$SampleY)

# Cell label colour: white on darker fill, black on lighter (sqrt scale)
sqrt_max <- sqrt(max(df$SNPs, na.rm = TRUE))
df$use_white <- sqrt(df$SNPs) > sqrt_max * 0.45

# Left annotation strip
ann <- data.frame(
  SampleY = factor(ylabels, levels = rev(ylabels)),
  Clone   = map_clone(ylabels),
  stringsAsFactors = FALSE
)

# Axis text colours (Okabe-Ito), aligned to factor level order
x_cols <- unname(clone_colors[as.character(map_clone(xlabels))])
y_cols <- unname(clone_colors[as.character(map_clone(rev(ylabels)))])
x_cols[is.na(x_cols)] <- "grey40"
y_cols[is.na(y_cols)] <- "grey40"

message("Building compact ggplot heatmap for main-figure panel ...")

# Compact panel: label ONLY small SNP counts (within/near clone); hide large between-clone numbers
# Colour bar + Clone legend carry the scale; keep both readable
label_max <- 150L  # show integers only when SNPs <= this
df_lab <- subset(df, SNPs <= label_max)
message("  Cell labels: only SNPs <= ", label_max, " (", nrow(df_lab), " cells); large distances hidden.")

legend_df <- data.frame(
  SampleX = factor(xlabels[1], levels = xlabels),
  SampleY = factor(ylabels[1], levels = rev(ylabels)),
  Clone   = factor(clone_levels, levels = clone_levels)
)

p_heat <- ggplot(df, aes(x = SampleX, y = SampleY)) +
  geom_tile(aes(fill = SNPs), color = "white", linewidth = 0.2) +
  geom_text(
    data = df_lab,
    aes(label = SNPs),
    colour = "black", size = 2.0, fontface = "bold", show.legend = FALSE
  ) +
  scale_fill_gradientn(
    name = "SNPs",
    colours = c("#F7F7F7", "#C6DBEF", "#6BAED6", "#2171B5"),
    trans = "sqrt",
    breaks = c(0, 1000, 5000, 20000),
    labels = c("0", "1,000", "5,000", "20,000"),
    guide = guide_colorbar(
      barwidth = unit(2.6, "cm"),
      barheight = unit(0.28, "cm"),
      title.position = "top",
      title.hjust = 0.5,
      title.theme = element_text(size = 8, face = "bold"),
      label.theme = element_text(size = 6.5),
      ticks = TRUE,
      order = 1
    )
  ) +
  geom_point(
    data = legend_df,
    aes(x = SampleX, y = SampleY, colour = Clone),
    inherit.aes = FALSE,
    size = 0, stroke = 0, alpha = 0, show.legend = TRUE
  ) +
  scale_colour_manual(
    name = "Clone",
    values = clone_colors,
    breaks = clone_levels,
    drop = FALSE,
    guide = guide_legend(
      override.aes = list(size = 3.2, shape = 15, alpha = 1, stroke = 0),
      direction = "horizontal",
      nrow = 1,
      title.position = "top",
      title.hjust = 0.5,
      title.theme = element_text(size = 8, face = "bold"),
      label.theme = element_text(size = 7),
      keywidth = unit(0.45, "cm"),
      keyheight = unit(0.45, "cm"),
      order = 2
    )
  ) +
  coord_fixed() +
  labs(
    title = "Core SNP distances",
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 8) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5, size = 9),
    axis.text.x   = element_text(angle = 45, hjust = 1, vjust = 1, size = 6, colour = x_cols),
    axis.text.y   = element_text(size = 6, colour = y_cols),
    axis.ticks    = element_line(linewidth = 0.25),
    panel.grid    = element_blank(),
    legend.position = "bottom",
    legend.box    = "vertical",
    legend.box.just = "center",
    legend.margin = margin(2, 4, 0, 4),
    legend.spacing.y = unit(3, "pt"),
    plot.margin   = margin(4, 4, 2, 2)
  )

final_plot <- p_heat
message("  Labels only for small SNP counts; SNPs + Clone legends at bottom.")

# --- Save --------------------------------------------------------------------
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
  message("Created output directory: ", outdir)
}

pdf_path <- file.path(outdir, "heatmap.pdf")
png_path <- file.path(outdir, "heatmap.png")

# Main-figure panel footprint (~matches Example strip slot)
panel_w <- 3.1
panel_h <- 3.6

message("Saving ", pdf_path, " (", panel_w, " x ", panel_h, " in) ...")
ggsave(pdf_path, plot = final_plot, width = panel_w, height = panel_h, units = "in", device = "pdf")

message("Saving ", png_path, " (300 dpi) ...")
ggsave(png_path, plot = final_plot, width = panel_w, height = panel_h, units = "in", dpi = 300, device = "png")

message("Done. Outputs:")
message("  ", normalizePath(pdf_path, winslash = "/", mustWork = FALSE))
message("  ", normalizePath(png_path, winslash = "/", mustWork = FALSE))
