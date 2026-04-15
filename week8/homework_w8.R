# Version info: R 4.2.2, Biobase 2.58.0, GEOquery 2.66.0, limma 3.54.0
################################################################
#   Differential expression analysis with limma
library(GEOquery)
library(limma)

# load series and platform data from GEO
gset <- getGEO("GSE34920", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL10558", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- "00112233"
sml <- strsplit(gsms, split="")[[1]]

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
          (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
  exprs(gset) <- log2(ex) }

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("CON","NANOG","SOX2","OCT4"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

gset <- gset[complete.cases(exprs(gset)), ] # skip missing values
feature_gene_ids <- fData(gset)$Gene.symbol
valid_feature_gene_id <- !is.na(feature_gene_ids) & trimws(feature_gene_ids) != ""
gset <- gset[valid_feature_gene_id, ] # keep only probes with annotated gene symbols

fit <- lmFit(gset, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- c("OCT4-CON", "SOX2-CON", "NANOG-CON")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and tables of top significant genes
fit2 <- eBayes(fit2, 0.01)
rank_by <- "B"  # 可改成 "adj.P.Val", "P.Value", "logFC", "AveExpr", "t", "B"
top_n <- 500
output_dir <- "week8/outputs"
venn_circle_cols <- c(OCT4 = "#E6AB02", SOX2 = "#E7298A", NANOG = "#7570B3")
dir.create(output_dir, showWarnings = FALSE)
gene_sets <- list()

for (contrast_name in cts) {
  tT <- topTable(fit2, coef=contrast_name, adjust="fdr", sort.by="B", number=Inf)
  tf_name <- sub("-CON", "", contrast_name, fixed = TRUE)

  if (!rank_by %in% colnames(tT)) {
    stop(
      "rank_by must be one of: ",
      paste(colnames(tT), collapse = ", ")
    )
  }

  tT <- tT[order(
    tT[[rank_by]],
    decreasing = !rank_by %in% c("adj.P.Val", "P.Value"),
    na.last = NA
  ), ]
  tT$Gene.symbol <- as.character(tT$Gene.symbol)
  tT <- tT[!duplicated(tT$Gene.symbol), , drop = FALSE]
  tT <- tT[seq_len(min(top_n, nrow(tT))), , drop = FALSE]

  table_cols <- intersect(
    c("ID", "adj.P.Val", "P.Value", "logFC", "AveExpr", "t", "B", "GI", "Gene.symbol", "Gene.title"),
    colnames(tT)
  )
  contrast_label <- sub("-", "_vs_", contrast_name)

  write.table(
    tT[, table_cols, drop = FALSE],
    file=file.path(output_dir, paste0("GSE34920_", contrast_label, "_top", top_n, ".tsv")),
    row.names=TRUE,
    sep="\t",
    quote=FALSE
  )
  gene_sets[[tf_name]] <- tT$Gene.symbol
}

all_gene_symbols <- sort(unique(unlist(gene_sets, use.names = FALSE)))
venn_input <- cbind(
  OCT4 = as.integer(all_gene_symbols %in% gene_sets$OCT4),
  SOX2 = as.integer(all_gene_symbols %in% gene_sets$SOX2),
  NANOG = as.integer(all_gene_symbols %in% gene_sets$NANOG)
)
rownames(venn_input) <- all_gene_symbols
venn_unts <- vennCounts(venn_input)

png(
  filename = file.path(output_dir, "GSE34920_top500_gene_symbol_venn.png"),
  width = 1600,
  height = 1600,
  res = 200
)
vennDiagram(
  venn_counts,
  circle.col = venn_circle_cols[colnames(venn_input)],
  names = colnames(venn_input),
  main = "GSE34920 Top 500 DE Gene Overlap"
)
dev.off()

draw_sankey_ribbon <- function(x_left, x_right, y_left_bottom, y_left_top,
                               y_right_bottom, y_right_top, fill_col) {
  interp <- seq(0, 1, length.out = 60)
  x_coords <- seq(x_left, x_right, length.out = 60)
  weight <- (1 - cos(pi * interp)) / 2
  lower <- y_left_bottom + (y_right_bottom - y_left_bottom) * weight
  upper <- y_left_top + (y_right_top - y_left_top) * weight
  polygon(
    x = c(x_coords, rev(x_coords)),
    y = c(lower, rev(upper)),
    col = fill_col,
    border = NA
  )
}

pattern_df <- data.frame(
  pattern = apply(venn_input, 1, paste0, collapse = ""),
  stringsAsFactors = FALSE
)
pattern_df$count <- 1L
pattern_df <- aggregate(count ~ pattern, data = pattern_df, FUN = sum)
pattern_df <- pattern_df[pattern_df$count > 0, , drop = FALSE]
pattern_df$OCT4 <- substr(pattern_df$pattern, 1, 1)
pattern_df$SOX2 <- substr(pattern_df$pattern, 2, 2)
pattern_df$NANOG <- substr(pattern_df$pattern, 3, 3)
pattern_df <- pattern_df[order(pattern_df$OCT4, pattern_df$SOX2, pattern_df$NANOG), ]
rownames(pattern_df) <- NULL

stage_names <- c("OCT4", "SOX2", "NANOG")
stage_positions <- c(0.12, 0.5, 0.88)
bar_width <- 0.05
union_size <- length(all_gene_symbols)
stage_gap <- max(20, ceiling(union_size * 0.04))
stage_height <- union_size + stage_gap
flow_cols <- adjustcolor(hcl.colors(nrow(pattern_df), "Dark 3"), alpha.f = 0.5)
stage_segments <- vector("list", length(stage_names))
stage_totals <- vector("list", length(stage_names))

for (stage_idx in seq_along(stage_names)) {
  stage_name <- stage_names[stage_idx]
  stage_segments[[stage_idx]] <- data.frame(
    bottom = numeric(nrow(pattern_df)),
    top = numeric(nrow(pattern_df))
  )
  stage_totals[[stage_idx]] <- list()
  current_top <- stage_height

  for (state in c("1", "0")) {
    state_rows <- which(pattern_df[[stage_name]] == state)
    state_rows <- state_rows[order(pattern_df$pattern[state_rows])]
    state_total <- sum(pattern_df$count[state_rows])
    state_top <- current_top
    state_bottom <- current_top - state_total

    stage_totals[[stage_idx]][[state]] <- c(bottom = state_bottom, top = state_top)

    running_top <- state_top
    for (row_idx in state_rows) {
      segment_top <- running_top
      segment_bottom <- running_top - pattern_df$count[row_idx]
      stage_segments[[stage_idx]]$bottom[row_idx] <- segment_bottom
      stage_segments[[stage_idx]]$top[row_idx] <- segment_top
      running_top <- segment_bottom
    }

    current_top <- state_bottom - stage_gap
  }
}

png(
  filename = file.path(output_dir, "GSE34920_top500_gene_symbol_sankey.png"),
  width = 2200,
  height = 1400,
  res = 200
)
par(mar = c(3, 3, 5, 3))
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, stage_height + stage_gap * 0.5), xaxs = "i", yaxs = "i")

for (stage_idx in seq_len(length(stage_names) - 1)) {
  x_left <- stage_positions[stage_idx] + bar_width / 2
  x_right <- stage_positions[stage_idx + 1] - bar_width / 2

  for (row_idx in seq_len(nrow(pattern_df))) {
    draw_sankey_ribbon(
      x_left = x_left,
      x_right = x_right,
      y_left_bottom = stage_segments[[stage_idx]]$bottom[row_idx],
      y_left_top = stage_segments[[stage_idx]]$top[row_idx],
      y_right_bottom = stage_segments[[stage_idx + 1]]$bottom[row_idx],
      y_right_top = stage_segments[[stage_idx + 1]]$top[row_idx],
      fill_col = flow_cols[row_idx]
    )
  }
}

for (stage_idx in seq_along(stage_names)) {
  x_center <- stage_positions[stage_idx]
  stage_name <- stage_names[stage_idx]

  rect(
    xleft = x_center - bar_width / 2,
    ybottom = stage_totals[[stage_idx]][["1"]]["bottom"],
    xright = x_center + bar_width / 2,
    ytop = stage_totals[[stage_idx]][["1"]]["top"],
    col = venn_circle_cols[stage_name],
    border = NA
  )
  rect(
    xleft = x_center - bar_width / 2,
    ybottom = stage_totals[[stage_idx]][["0"]]["bottom"],
    xright = x_center + bar_width / 2,
    ytop = stage_totals[[stage_idx]][["0"]]["top"],
    col = "#D9D9D9",
    border = NA
  )

  text(
    x = x_center,
    y = stage_height + stage_gap * 0.25,
    labels = stage_name,
    font = 2,
    cex = 1.2
  )
  text(
    x = x_center,
    y = mean(stage_totals[[stage_idx]][["1"]]),
    labels = paste0("In top 500\nn=", sum(pattern_df$count[pattern_df[[stage_name]] == "1"])),
    cex = 0.8
  )
  text(
    x = x_center,
    y = mean(stage_totals[[stage_idx]][["0"]]),
    labels = paste0("Not in top 500\nn=", sum(pattern_df$count[pattern_df[[stage_name]] == "0"])),
    cex = 0.8
  )
}

title(
  main = "Sankey View of Top 500 DE Gene Membership",
  sub = "Union of unique gene symbols across OCT4, SOX2 and NANOG contrasts"
)
dev.off()

co_regulated_genes <- sort(intersect(
  intersect(gene_sets$OCT4, gene_sets$SOX2),
  gene_sets$NANOG
))

write.table(
  data.frame(Gene.symbol = co_regulated_genes),
  file = file.path(output_dir, "GSE34920_co_regulated_genes.tsv"),
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)
# BNC2, CCDC34, CCL26, CHEK2, CHST9, GJA1 ,GLTP, HESX1,
# KCNK12, KLF6, KLKB1, LFNG, LYPD1, MAOA, MLLT6, PQLC1,
# PROM1, RGS17, SAT1, SHTN1, SSBP2, STC1


overlap_summary <- data.frame(
  comparison = c(
    "OCT4",
    "SOX2",
    "NANOG",
    "OCT4_AND_SOX2",
    "OCT4_AND_NANOG",
    "SOX2_AND_NANOG",
    "OCT4_AND_SOX2_AND_NANOG"
  ),
  count = c(
    length(gene_sets$OCT4),
    length(gene_sets$SOX2),
    length(gene_sets$NANOG),
    length(intersect(gene_sets$OCT4, gene_sets$SOX2)),
    length(intersect(gene_sets$OCT4, gene_sets$NANOG)),
    length(intersect(gene_sets$SOX2, gene_sets$NANOG)),
    length(co_regulated_genes)
  )
)

write.table(
  overlap_summary,
  file = file.path(output_dir, "GSE34920_overlap_summary.tsv"),
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

cat("Co-regulated genes shared by OCT4, SOX2 and NANOG:", length(co_regulated_genes), "\n")
