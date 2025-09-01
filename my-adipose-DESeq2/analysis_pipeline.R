# -------------------------  RNA-seq Pipeline (GitHub-ready) -------------------------
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(DESeq2)
  library(readr)
  library(dplyr)
  library(ComplexHeatmap)
  library(circlize)
  library(biomaRt)
  library(RColorBrewer)
  library(viridis)
  library(limma)
  library(MetaCycle)
  library(ggplot2)
  library(gtools)
  library(matrixStats)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(enrichplot)
  library(sva)
})

set.seed(42)

# -------------------------  0. Paths -------------------------
dir_counts   <- "data/counts"
path_meta    <- "data/metadata.csv"
dir_results  <- "results"
dir_plots    <- file.path(dir_results, "plots")
dir_deseq2   <- file.path(dir_results, "deseq2")
dir_go       <- file.path(dir_results, "go")
dir.create(dir_plots, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_deseq2, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_go, recursive = TRUE, showWarnings = FALSE)

# -------------------------  1. Load Metadata -------------------------
if (!file.exists(path_meta)) stop("Missing data/metadata.csv. See data/metadata_template.csv")
metadata <- read.csv(path_meta, header = TRUE)
# Ensure Filename has no .txt extension
metadata$Filename <- gsub("\\.txt$", "", metadata$Filename)

# QC for required columns
req_cols <- c("Filename", "Group", "Timepoint")
if (!all(req_cols %in% names(metadata))) {
  stop("metadata.csv must contain columns: ", paste(req_cols, collapse = ", "))
}

# Factors
metadata$Group     <- factor(metadata$Group, levels = c("Control","CMVS","CSDS"))
metadata$Timepoint <- factor(metadata$Timepoint,
                             levels = c("zt_2","zt_6","zt_10","zt_14","zt_18","zt_22"))

# -------------------------  2. Load Counts -------------------------
count_files <- list.files(dir_counts, pattern = "\\.txt$", full.names = TRUE)
if (length(count_files) == 0) stop("No .txt count files found in data/counts")

load_counts <- function(file_path) {
  # Expect two columns: GeneID, Count (no header)
  readr::read_tsv(file_path,
                  col_names = FALSE,
                  col_types = readr::cols(readr::col_character(), readr::col_double())) %>%
    setNames(c("GeneID", "Count"))
}

count_data_list <- lapply(count_files, load_counts)

# Gene IDs (from the first file)
gene_ids <- count_data_list[[1]]$GeneID

# Combine counts
combined_counts <- do.call(cbind, lapply(count_data_list, function(df) df$Count))
combined_counts <- as.matrix(combined_counts)
rownames(combined_counts) <- gene_ids

# Map columns to metadata order using file basenames
file_keys <- tools::file_path_sans_ext(basename(count_files))
colnames(combined_counts) <- file_keys

# Order counts to match metadata$Filename
common <- intersect(metadata$Filename, colnames(combined_counts))
if (length(common) == 0) stop("No sample names overlap between metadata$Filename and counts")
combined_counts <- combined_counts[, match(metadata$Filename, colnames(combined_counts)), drop = FALSE]

# Filter low-support genes
min_samples_nonzero <- 8
keep <- rowSums(combined_counts > 0) >= min_samples_nonzero
combined_counts <- combined_counts[keep, ]

# Final ordered metadata (aligned to columns of combined_counts)
metadata_ordered <- metadata[match(colnames(combined_counts), metadata$Filename), ]
stopifnot(all(metadata_ordered$Filename == colnames(combined_counts)))

# -------------------------  3. DESeq2 base object -------------------------
dds <- DESeqDataSetFromMatrix(
  countData = combined_counts,
  colData   = metadata_ordered,
  design    = ~ Group
)
dds <- DESeq(dds)

# Save a default results table (no contrast) for bookkeeping
res_default <- results(dds)
write.csv(as.data.frame(res_default),
          file.path(dir_deseq2, "DESeq2_results_default.csv"), row.names = TRUE)

# -------------------------  4. SVA (surrogate variables) -------------------------
# Use normalized counts for svaseq
norm_counts <- counts(dds, normalized = TRUE)
norm_counts <- norm_counts[rowSums(norm_counts > 0) > 1, ]
norm_counts[is.na(norm_counts)] <- 0

mod  <- model.matrix(~ Group + Timepoint, data = metadata_ordered)
mod0 <- model.matrix(~ 1, data = metadata_ordered)

svseq <- svaseq(norm_counts, mod, mod0)
message("Surrogate variables detected: ", svseq$n.sv)

# Add SVs into metadata
for (i in seq_len(svseq$n.sv)) {
  metadata_ordered[[paste0("SV", i)]] <- svseq$sv[, i]
}

sva_info <- list(n_sv = svseq$n.sv, sv = svseq$sv)

# -------------------------  5. Helper: DESeq2 pipeline with optional SVs -------------------------
run_deseq2_pipeline <- function(counts, metadata,
                                groups_keep = NULL,
                                times_keep  = NULL,
                                contrast    = NULL,    # e.g., c("Group","CSDS","CMVS")
                                out_prefix  = "DESeq2",
                                sv_terms    = NULL) {
  # 1) subset metadata
  meta_sub <- metadata
  if (!is.null(groups_keep)) meta_sub <- meta_sub[meta_sub$Group %in% groups_keep, , drop = FALSE]
  if (!is.null(times_keep))  meta_sub <- meta_sub[meta_sub$Timepoint %in% times_keep, , drop = FALSE]
  if (nrow(meta_sub) == 0) stop("Subset produced zero samples; check groups_keep/times_keep")
  
  # 2) match counts
  counts_sub <- counts[, meta_sub$Filename, drop = FALSE]
  keep <- rowSums(counts_sub > 0) >= 8
  counts_sub <- counts_sub[keep, , drop = FALSE]
  
  # 3) factors and reference
  meta_sub$Group <- factor(meta_sub$Group, levels = levels(metadata$Group))
  if (!is.null(contrast)) {
    ref_group <- contrast[3]
    meta_sub$Group <- relevel(meta_sub$Group, ref = ref_group)
  }
  
  # 4) design formula
  if (!is.null(sv_terms) && length(sv_terms) > 0) {
    sv_terms <- sv_terms[sv_terms %in% colnames(meta_sub)]
    design_formula <- as.formula(paste("~", paste(c(sv_terms, "Group"), collapse = " + ")))
  } else {
    design_formula <- ~ Group
  }
  
  dds_sub <- DESeqDataSetFromMatrix(countData = counts_sub, colData = meta_sub, design = design_formula)
  dds_sub <- DESeq(dds_sub)
  
  # 5) results
  if (!is.null(contrast)) {
    res <- results(dds_sub, contrast = contrast)
  } else {
    res <- results(dds_sub)
  }
  res <- res[order(res$padj), ]
  
  # 6) annotate gene symbols via biomaRt
  genes <- rownames(res)
  genes_clean <- sub("\\..*$", "", genes)
  mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  gene_info <- getBM(
    attributes = c("ensembl_gene_id", "external_gene_name"),
    filters = "ensembl_gene_id",
    values = genes_clean,
    mart = mart
  )
  res$GeneSymbol <- gene_info$external_gene_name[match(genes_clean, gene_info$ensembl_gene_id)]
  
  # 7) VST matrix (regress out SVs for visualization only)
  vst_mat <- assay(vst(dds_sub, blind = TRUE))
  design_vis <- model.matrix(~ Group + Timepoint, data = meta_sub)
  if (!is.null(sv_terms) && length(sv_terms) > 0) {
    sv_mat <- as.matrix(meta_sub[, sv_terms, drop = FALSE])
    vst_mat <- limma::removeBatchEffect(vst_mat, covariates = sv_mat, design = design_vis)
  }
  
  # 8) save results
  out_csv <- file.path(dir_deseq2, paste0(out_prefix, "_results.csv"))
  write.csv(as.data.frame(res), out_csv, row.names = TRUE)
  
  list(dds = dds_sub, res = res, vst_mat = vst_mat, metadata = meta_sub)
}

# -------------------------  6. Plots & helpers -------------------------

plot_pca_subset <- function(vst_mat, metadata, ntop = 500, pc_x = 1, pc_y = 2,
                            color_by = "Group", shape_by = "Timepoint", add_labels = FALSE) {
  metadata[[color_by]] <- factor(metadata[[color_by]])
  if (!is.null(shape_by)) metadata[[shape_by]] <- factor(metadata[[shape_by]])
  vst_mat <- vst_mat[, metadata$Filename, drop = FALSE]
  rv <- matrixStats::rowVars(vst_mat)
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  vst_top <- vst_mat[select, ]
  pca <- prcomp(t(vst_top), center = TRUE, scale. = FALSE)
  var_expl <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 1)
  pca_df <- data.frame(PCX = pca$x[, pc_x], PCY = pca$x[, pc_y], metadata)
  p <- ggplot(pca_df, aes(x = PCX, y = PCY, color = .data[[color_by]])) + geom_point(size = 3)
  if (!is.null(shape_by)) p <- p + aes(shape = .data[[shape_by]])
  if (add_labels) p <- p + geom_text(aes(label = Filename), vjust = -0.5, size = 2.5, check_overlap = TRUE)
  p + theme_minimal() +
    labs(title = paste0("PCA: PC", pc_x, " vs PC", pc_y, " (Top ", ntop, ")"),
         x = paste0("PC", pc_x, " (", var_expl[pc_x], "%)"),
         y = paste0("PC", pc_y, " (", var_expl[pc_y], "%)"))
}

count_degs <- function(res, alpha = 0.05) {
  sig <- res[!is.na(res$padj) & res$padj < alpha, ]
  n_up   <- sum(sig$log2FoldChange > 0, na.rm = TRUE)
  n_down <- sum(sig$log2FoldChange < 0, na.rm = TRUE)
  list(up = n_up, down = n_down, total = n_up + n_down)
}

plot_volcano <- function(res, alpha = 0.05, title = "Volcano Plot") {
  df <- as.data.frame(res)
  df$Significant <- "Not Sig"
  df$Significant[!is.na(df$padj) & df$padj < alpha & df$log2FoldChange > 0]  <- "Up"
  df$Significant[!is.na(df$padj) & df$padj < alpha & df$log2FoldChange < 0]  <- "Down"
  ggplot(df, aes(x = log2FoldChange, y = -log10(padj), color = Significant)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Up" = "firebrick", "Down" = "steelblue", "Not Sig" = "grey70")) +
    geom_hline(yintercept = -log10(alpha), linetype = "dashed", color = "black") +
    labs(title = title, x = "log2 Fold Change", y = "-log10(padj)") +
    theme_minimal()
}

summarise_deg_counts <- function(counts, metadata, group1, group2, alpha = 0.05) {
  out_list <- list()
  run_global <- run_deseq2_pipeline(
    counts, metadata,
    groups_keep = c(group1, group2),
    contrast = c("Group", group1, group2),
    out_prefix = paste0(group1, "_vs_", group2, "_AllTimes")
  )
  out_list$Global <- count_degs(run_global$res, alpha)
  day_points   <- c("zt_2","zt_6","zt_10")
  night_points <- c("zt_14","zt_18","zt_22")
  run_day <- run_deseq2_pipeline(
    counts, metadata,
    groups_keep = c(group1, group2),
    times_keep = day_points,
    contrast = c("Group", group1, group2),
    out_prefix = paste0(group1, "_vs_", group2, "_Day")
  )
  out_list$Day <- count_degs(run_day$res, alpha)
  run_night <- run_deseq2_pipeline(
    counts, metadata,
    groups_keep = c(group1, group2),
    times_keep = night_points,
    contrast = c("Group", group1, group2),
    out_prefix = paste0(group1, "_vs_", group2, "_Night")
  )
  out_list$Night <- count_degs(run_night$res, alpha)
  tp_counts <- list()
  for (tp in levels(metadata$Timepoint)) {
    if (tp %in% c(day_points, night_points)) {
      run_tp <- run_deseq2_pipeline(
        counts, metadata,
        groups_keep = c(group1, group2),
        times_keep = tp,
        contrast = c("Group", group1, group2),
        out_prefix = paste0(group1, "_vs_", group2, "_", tp)
      )
      tp_counts[[tp]] <- count_degs(run_tp$res, alpha)
    }
  }
  out_list$Timepoints <- tp_counts
  out_list
}

plot_deg_summary <- function(deg_summary, group1, group2) {
  tp_labels <- names(deg_summary$Timepoints)
  tp_labels_pretty <- gsub("zt_(\\d+)", "ZT \\1", tp_labels)
  df <- data.frame(
    Category = c("Global", "Day", "Night", tp_labels_pretty),
    Up   = c(deg_summary$Global$up,   deg_summary$Day$up,   deg_summary$Night$up,
             sapply(deg_summary$Timepoints, function(x) x$up)),
    Down = c(deg_summary$Global$down, deg_summary$Day$down, deg_summary$Night$down,
             sapply(deg_summary$Timepoints, function(x) x$down))
  )
  df$Total <- df$Up + df$Down
  category_order <- c("Global","Day","Night", paste0("ZT ", c(2,6,10,14,18,22)))
  df$Category <- factor(df$Category, levels = category_order)
  df_melt <- tidyr::pivot_longer(df, cols = c("Up","Down"),
                                 names_to = "Direction", values_to = "Count")
  ggplot(df_melt, aes(x = Category, y = Count, fill = Direction)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = c("Up" = "firebrick", "Down" = "steelblue")) +
    labs(title = paste0("DEG Summary: ", group1, " vs ", group2),
         y = "Number of DEGs", x = "Comparison") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

plot_deg_heatmap <- function(vst_mat, res, metadata, top_n = 50,
                             color_by = "Group", shape_by = "Timepoint",
                             scale_rows = TRUE, title = "Top DEGs Heatmap") {
  if (!"SampleID" %in% colnames(metadata)) {
    metadata$SampleID <- paste(metadata$Group, metadata$Timepoint, sep = "_")
  }
  degs <- res[!is.na(res$padj), ]
  top_degs <- head(degs[order(degs$padj), ], n = top_n)
  top_genes <- rownames(top_degs)
  mat_sub <- vst_mat[top_genes, metadata$Filename, drop = FALSE]
  if ("GeneSymbol" %in% colnames(res)) {
    gene_map <- data.frame(Ensembl = rownames(res), Symbol = res$GeneSymbol)
    symbols <- gene_map$Symbol[match(top_genes, gene_map$Ensembl)]
    symbols[is.na(symbols) | symbols == ""] <- top_genes[is.na(symbols) | symbols == ""]
    rownames(mat_sub) <- symbols
  }
  colnames(mat_sub) <- metadata$SampleID[match(colnames(mat_sub), metadata$Filename)]
  if (scale_rows) mat_sub <- t(scale(t(mat_sub)))
  metadata[[color_by]] <- factor(metadata[[color_by]])
  metadata[[shape_by]] <- factor(metadata[[shape_by]])
  annotation_colors <- list()
  annotation_colors[[color_by]] <- setNames(
    colorRampPalette(brewer.pal(9, "Set1"))(nlevels(metadata[[color_by]])),
    levels(metadata[[color_by]])
  )
  annotation_colors[[shape_by]] <- setNames(
    viridis::viridis(nlevels(metadata[[shape_by]])),
    levels(metadata[[shape_by]])
  )
  ha <- HeatmapAnnotation(df = metadata[, c(color_by, shape_by), drop = FALSE], col = annotation_colors)
  Heatmap(mat_sub, name = ifelse(scale_rows, "z-score", "VST"),
          top_annotation = ha, show_row_names = TRUE, show_column_names = FALSE,
          cluster_rows = TRUE, cluster_columns = TRUE,
          column_title = title, row_names_gp = gpar(fontsize = 8),
          column_names_gp = gpar(fontsize = 8))
}

run_go_analysis <- function(res, alpha = 0.05, lfc_threshold = 1,
                            ontology = "BP", top_n = 20, title = "GO Enrichment") {
  sig <- as.data.frame(res)
  sig <- sig[!is.na(sig$padj) & sig$padj < alpha & abs(sig$log2FoldChange) > lfc_threshold, ]
  genes <- unique(sig$GeneSymbol[!is.na(sig$GeneSymbol) & sig$GeneSymbol != ""])
  if (length(genes) < 5) {
    message("Not enough DEGs for GO analysis")
    return(NULL)
  }
  ego <- enrichGO(
    gene          = genes,
    OrgDb         = org.Mm.eg.db,
    keyType       = "SYMBOL",
    ont           = ontology,
    pAdjustMethod = "BH",
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    p <- dotplot(ego, showCategory = top_n, title = title) +
      theme(axis.text.y = element_text(size = 8))
    print(p)
  } else {
    message("No significant GO terms found")
  }
  return(ego)
}

# -------------------------  7. USAGE (examples you can edit) -------------------------

# How many SVs to use in visualization/design
n_sv_keep <- min(3, svaseq(norm_counts, mod, mod0)$n.sv)  # recompute quickly just to be safe
sv_terms  <- paste0("SV", 1:n_sv_keep)

# Example: CSDS vs CMVS, all timepoints, include SVs
Deseq2object_sva <- run_deseq2_pipeline(
  combined_counts, metadata_ordered,
  groups_keep = c("CSDS", "CMVS"),
  contrast   = c("Group", "CSDS", "CMVS"),
  out_prefix = "CSDS_vs_CMVS_AllTimes_SVA",
  sv_terms   = sv_terms
)

# Save example plots
pdf(file.path(dir_plots, "PCA_CSDS_vs_CMVS.pdf"))
print(plot_pca_subset(Deseq2object_sva$vst_mat, Deseq2object_sva$metadata, ntop = 1000, pc_x = 1, pc_y = 2))
dev.off()

pdf(file.path(dir_plots, "Volcano_CSDS_vs_CMVS.pdf"))
print(plot_volcano(Deseq2object_sva$res, title = "CSDS vs CMVS (All Times)"))
dev.off()

pdf(file.path(dir_plots, "Heatmap_top50_CSDS_vs_CMVS.pdf"), width = 8, height = 10)
ht <- plot_deg_heatmap(Deseq2object_sva$vst_mat, Deseq2object_sva$res, Deseq2object_sva$metadata,
                       top_n = 50,
                       title = "CSDS vs CMVS: Top 50 DEGs")
ComplexHeatmap::draw(ht)
dev.off()

# Example GO analysis
pdf(file.path(dir_plots, "GO_CSDS_vs_CMVS.pdf"))
go_global <- run_go_analysis(Deseq2object_sva$res,
                             alpha = 0.05, lfc_threshold = 1,
                             ontology = "BP", top_n = 20,
                             title = "GO Enrichment: CSDS vs CMVS (All Times)")
dev.off()

# Summaries by time (optional)
# deg_summary <- summarise_deg_counts(combined_counts, metadata_ordered, "CSDS", "CMVS")
# pdf(file.path(dir_plots, "DEG_Summary_CSDS_vs_CMVS.pdf"))
# print(plot_deg_summary(deg_summary, "CSDS", "CMVS"))
# dev.off()
