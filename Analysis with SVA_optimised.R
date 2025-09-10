# -------------------------  1. Load Libraries -------------------------
library(DESeq2)
library(readr)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(biomaRt)
library(RColorBrewer)
library(viridis)
library(readxl)
library(limma)
library(MetaCycle)
library(ggplot2)
library(gtools)
library(matrixStats)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(sva)

# -------------------------  2. Set Working Directory -------------------
setwd("~/Documents/adiposebulkstress/BAT")

# -------------------------  3. Connect to Ensembl ----------------------
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# -------------------------  4. Load Count Data -------------------------
count_files <- list.files(pattern = "\\.txt$")

load_counts <- function(file_name) {
  full_path <- file.path(getwd(), file_name)
  read_tsv(full_path,
           col_names = FALSE,
           col_types = cols(col_character(), col_double())) %>%
    setNames(c("GeneID", "Count"))
}

count_data_list <- lapply(count_files, load_counts)

# Extract Gene IDs (should all be the same across files)
gene_ids <- count_data_list[[1]]$GeneID

# Build combined count matrix
combined_counts <- do.call(cbind, lapply(count_data_list, function(df) df$Count))
combined_counts <- as.matrix(combined_counts)  # ensure it's a matrix
rownames(combined_counts) <- gene_ids
colnames(combined_counts) <- sub("\\.txt$", "", count_files)  # sample names

# Now filter
keep <- rowSums(combined_counts > 0) >= 8
combined_counts <- combined_counts[keep, ]

# -------------------------  5. Load Metadata ---------------------------
metadata_path <- "Analysis metadata_gn_BAT.csv"
metadata <- read.csv(metadata_path, header = TRUE)
metadata$Filename <- gsub("\\.txt$", "", metadata$Filename)

metadata_ordered <- metadata[match(colnames(combined_counts), metadata$Filename), ]
metadata_ordered$Batch <- factor(metadata_ordered$Batch)
metadata_ordered$Group <- factor(metadata_ordered$Group)
metadata_ordered$Timepoint <- factor(metadata_ordered$Timepoint)

if (any(is.na(metadata_ordered$Filename))) {
  stop("Mismatch between metadata and count data files")
}

# -------------------------  6. DESeq2 & VST ---------------------------

# Create composite condition: e.g., "ZT2_Control", "ZT6_E", etc.
metadata_ordered$condition <- factor(
  paste(metadata_ordered$Timepoint, metadata_ordered$Group, sep = "_")
)

dds <- DESeqDataSetFromMatrix(
  countData = combined_counts,
  colData   = metadata_ordered,
  design    = ~ condition
)

dds <- DESeq(dds)
dds$sample_id <- rownames(colData(dds))

# Optional — size factor & dispersion estimates are already part of DESeq(), so this is redundant:
# dds <- estimateSizeFactors(dds)
# dds <- estimateDispersions(dds)

# Get results
res <- results(dds)

genes <- rownames(res)
gene_info <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'), 
                   filters = 'ensembl_gene_id', 
                   values = genes, 
                   mart = mart)
res$GeneSymbol <- gene_info[match(genes, gene_info$ensembl_gene_id), 'external_gene_name']
write.csv(as.data.frame(res), "DESeq2_Results_with_Symbols_2.csv")

vst_transformed <- vst(dds, blind = FALSE)
vst_matrix <- assay(vst_transformed)

#batch <- metadata_ordered$Batch
#vst_matrix <- removeBatchEffect(vst_matrix, batch = batch)

vst_gene_symbols <- gene_info$external_gene_name[match(rownames(vst_matrix), gene_info$ensembl_gene_id)]

################ SVA

# Use normalized counts or vst matrix for SVA
norm_counts <- counts(dds, normalized = TRUE)

# Make sure Group and Timepoint are factors
metadata_ordered$Group <- factor(metadata_ordered$Group)
metadata_ordered$Timepoint <- factor(metadata_ordered$Timepoint)

# Full model: includes your variables of interest
mod <- model.matrix(~ Group + Timepoint, data = metadata_ordered)

# Null model: intercept only
mod0 <- model.matrix(~ 1, data = metadata_ordered)

# Make sure we use a clean count matrix
norm_counts <- counts(dds, normalized = TRUE)

# Remove genes with all zeros (or nearly all zeros)
keep <- rowSums(norm_counts > 0) > 1
norm_counts <- norm_counts[keep, ]

# Replace any NA with 0 (shouldn’t normally happen, but safe)
norm_counts[is.na(norm_counts)] <- 0

# Re-run SVA
mod <- model.matrix(~ Group + Timepoint, data = metadata_ordered)
mod0 <- model.matrix(~ 1, data = metadata_ordered)

svseq <- svaseq(norm_counts, mod, mod0)

# Run svaseq (on counts, not vst, because it models mean-variance)
svseq <- svaseq(norm_counts, mod, mod0)

# Number of surrogate variables detected
cat("Number of surrogate variables:", svseq$n.sv, "\n")

# Add surrogate variables back to metadata
for (i in 1:svseq$n.sv) {
  metadata_ordered[[paste0("SV", i)]] <- svseq$sv[, i]
}

# Quick look at SV correlations
sv_cor <- cor(metadata_ordered[, paste0("SV", 1:svseq$n.sv)], 
              model.matrix(~ Group + Timepoint, data = metadata_ordered))
print(sv_cor)


# -------------------------  DESEQ2 pipeline -------------------------


run_deseq2_pipeline <- function(counts, metadata,
                                groups_keep = NULL,
                                times_keep  = NULL,
                                contrast    = NULL,
                                out_prefix  = "DESeq2",
                                sv_terms    = NULL) {
  
  # --- 1. Subset metadata ---
  meta_sub <- metadata
  if (!is.null(groups_keep)) {
    meta_sub <- meta_sub[meta_sub$Group %in% groups_keep, ]
  }
  if (!is.null(times_keep)) {
    meta_sub <- meta_sub[meta_sub$Timepoint %in% times_keep, ]
  }
  
  # drop unused factor levels and ensure Filename order
  meta_sub <- droplevels(meta_sub)
  rownames(meta_sub) <- meta_sub$Filename
  
  # --- 2. Match counts ---
  counts_sub <- counts[, rownames(meta_sub), drop = FALSE]
  if (ncol(counts_sub) == 0) stop("No samples left after subsetting counts/metadata.")
  
  # Filter out weak genes (require ≥8 samples with counts OR adapt threshold)
  keep <- rowSums(counts_sub > 0) >= 8
  counts_sub <- counts_sub[keep, ]
  if (nrow(counts_sub) == 0) stop("No genes left after filtering by expression threshold.")
  
  # --- 3. Ensure Group factor and set reference ---
  # Important: enforce factors and drop unused levels
  meta_sub$Group <- factor(as.character(meta_sub$Group))
  if (!is.null(contrast)) {
    ref_group <- contrast[3]
    if (! ref_group %in% levels(meta_sub$Group)) {
      stop(sprintf("Reference group '%s' not present in subset meta_sub$Group levels: %s",
                   ref_group, paste(levels(meta_sub$Group), collapse = ",")))
    }
    meta_sub$Group <- relevel(meta_sub$Group, ref = ref_group)
    contrast <- c(contrast[1], contrast[2], ref_group)
  }
  # ensure Timepoint is factor (but we will drop it from colData unless used)
  if ("Timepoint" %in% colnames(meta_sub)) {
    meta_sub$Timepoint <- factor(as.character(meta_sub$Timepoint))
  }
  
  # --- 4. Build dynamic design formula & prepare colData for DESeq ---
  # Only include columns in colData that are used in the design (Group and sv_terms)
  design_terms <- c()
  if (!is.null(sv_terms) && length(sv_terms) > 0) {
    # only keep sv_terms that exist in meta_sub
    sv_terms <- sv_terms[sv_terms %in% colnames(meta_sub)]
    if (length(sv_terms) > 0) design_terms <- c(design_terms, sv_terms)
  }
  design_terms <- c(design_terms, "Group")
  design_formula <- as.formula(paste("~", paste(design_terms, collapse = " + ")))
  
  # prepare colData with only needed columns (and ensure types)
  colData <- meta_sub[, design_terms, drop = FALSE]
  # coerce any SVs to numeric (they may be character if messed up)
  for (st in sv_terms) {
    colData[[st]] <- as.numeric(as.character(colData[[st]]))
  }
  
  # DEBUG prints (remove or comment out if noisy)
  message("DESeq run: samples = ", nrow(colData),
          "; design = ", deparse(design_formula))
  message("Group levels: ", paste(levels(colData$Group), collapse = ", "))
  if (!is.null(sv_terms) && length(sv_terms) > 0) {
    message("SV terms used: ", paste(sv_terms, collapse = ", "))
  }
  
  # --- 5. Create DESeq dataset and run ---
  dds <- DESeqDataSetFromMatrix(
    countData = counts_sub,
    colData   = colData,
    design    = design_formula
  )
  
  # sanity check: ensure design variables have >1 level or are numeric
  if (any(sapply(design_terms, function(x) {
    is.factor(colData[[x]]) && length(levels(colData[[x]])) < 2
  }))) {
    stop("One or more design factors have <2 levels in this subset; aborting.")
  }
  
  dds <- DESeq(dds)
  
  # --- 6. Extract results ---
  if (!is.null(contrast)) {
    res <- results(dds, contrast = contrast)
  } else {
    res <- results(dds)
  }
  res <- res[order(res$padj), ]
  
  # --- 7. Annotate gene symbols (same as before) ---
  genes <- rownames(res)
  genes_clean <- sub("\\..*$", "", genes)  # remove version numbers
  
  mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  gene_info <- getBM(
    attributes = c("ensembl_gene_id", "external_gene_name"),
    filters = "ensembl_gene_id",
    values = genes_clean,
    mart = mart
  )
  
  res$GeneSymbol <- gene_info$external_gene_name[
    match(genes_clean, gene_info$ensembl_gene_id)
  ]
  
  # --- 8. VST and optional SV regression on VST ---
  vst_mat <- assay(vst(dds, blind = TRUE))
  if (!is.null(sv_terms) && all(sv_terms %in% colnames(meta_sub))) {
    sv_mat <- as.matrix(meta_sub[, sv_terms, drop = FALSE])
    design_for_remove <- model.matrix(~ Group, data = meta_sub)  # keep Group structure
    vst_mat <- limma::removeBatchEffect(vst_mat, covariates = sv_mat, design = design_for_remove)
  }
  
  return(list(
    dds      = dds,
    res      = res,
    vst_mat  = vst_mat,
    metadata = meta_sub
  ))
}
# ------------------------- PCA Function -------------------------
plot_pca_subset <- function(vst_mat, metadata, ntop = 500, pc_x = 1, pc_y = 2,
                            color_by = "Group", shape_by = "Timepoint", add_labels = FALSE) {
  
  # Ensure metadata columns are factors
  metadata[[color_by]] <- factor(metadata[[color_by]])
  if (!is.null(shape_by)) metadata[[shape_by]] <- factor(metadata[[shape_by]])
  
  # Subset vst_mat columns to match metadata exactly
  vst_mat <- vst_mat[, metadata$Filename, drop = FALSE]
  
  # Remove zero-variance genes
  rv <- rowVars(vst_mat)
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  vst_top <- vst_mat[select, ]
  
  # Run PCA
  pca <- prcomp(t(vst_top), center = TRUE, scale. = FALSE)
  var_expl <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 1)
  
  # Prepare dataframe for plotting (carry through *all* metadata cols)
  pca_df <- data.frame(
    PCX = pca$x[, pc_x],
    PCY = pca$x[, pc_y],
    metadata
  )
  
  # Plot
  p <- ggplot(pca_df, aes(x = PCX, y = PCY, color = .data[[color_by]])) +
    geom_point(size = 3)
  
  if (!is.null(shape_by)) p <- p + aes(shape = .data[[shape_by]])
  if (add_labels) p <- p + geom_text(aes(label = Filename), vjust = -0.5, size = 2.5, check_overlap = TRUE)
  
  p <- p +
    theme_minimal() +
    labs(
      title = paste0("PCA: PC", pc_x, " vs PC", pc_y, " (Top ", ntop, " variable genes)"),
      x = paste0("PC", pc_x, " (", var_expl[pc_x], "%)"),
      y = paste0("PC", pc_y, " (", var_expl[pc_y], "%)")
    )
  
  return(p)
}


# -------------------------  DEG Counting Helper -------------------------

# Function: Count DEGs based on padj threshold
count_degs <- function(res, alpha = 0.05) {
  sig <- res[!is.na(res$padj) & res$padj < alpha, ]
  n_up   <- sum(sig$log2FoldChange > 0, na.rm = TRUE)
  n_down <- sum(sig$log2FoldChange < 0, na.rm = TRUE)
  return(list(up = n_up, down = n_down, total = n_up + n_down))
}

# Volcano plot using SV-regressed VST fold changes
plot_volcano_sv <- function(deseq2_obj, alpha = 0.05, title = "Volcano Plot") {
  res <- deseq2_obj$res
  vst_mat <- deseq2_obj$vst_mat
  
  # Compute log2 fold change from SV-adjusted VST
  groups <- unique(deseq2_obj$metadata$Group)
  if(length(groups) != 2) stop("Volcano SV plot only supports pairwise comparison")
  
  grp1 <- groups[1]
  grp2 <- groups[2]
  
  samples_grp1 <- deseq2_obj$metadata$Filename[deseq2_obj$metadata$Group == grp1]
  samples_grp2 <- deseq2_obj$metadata$Filename[deseq2_obj$metadata$Group == grp2]
  
  log2FC_sv <- rowMeans(vst_mat[, samples_grp2, drop = FALSE]) -
    rowMeans(vst_mat[, samples_grp1, drop = FALSE])
  
  df <- as.data.frame(res)
  df$log2FoldChange <- log2FC_sv  # replace DESeq2 fold change with SV-adjusted
  df$Significant <- "Not Sig"
  df$Significant[!is.na(df$padj) & df$padj < alpha & df$log2FoldChange > 0] <- "Up"
  df$Significant[!is.na(df$padj) & df$padj < alpha & df$log2FoldChange < 0] <- "Down"
  
  ggplot(df, aes(x = log2FoldChange, y = -log10(padj), color = Significant)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Up" = "firebrick", "Down" = "steelblue", "Not Sig" = "grey70")) +
    geom_hline(yintercept = -log10(alpha), linetype = "dashed", color = "black") +
    labs(title = title, x = "log2 Fold Change (SV-adjusted)", y = "-log10(padj)") +
    theme_minimal()
}

# MA plot using SV-regressed VST fold changes
plot_MA_sv <- function(deseq2_obj, alpha = 0.05, title = "MA Plot") {
  res <- deseq2_obj$res
  vst_mat <- deseq2_obj$vst_mat
  
  # Compute log2 fold change from SV-adjusted VST
  groups <- unique(deseq2_obj$metadata$Group)
  if(length(groups) != 2) stop("MA SV plot only supports pairwise comparison")
  
  grp1 <- groups[1]
  grp2 <- groups[2]
  
  samples_grp1 <- deseq2_obj$metadata$Filename[deseq2_obj$metadata$Group == grp1]
  samples_grp2 <- deseq2_obj$metadata$Filename[deseq2_obj$metadata$Group == grp2]
  
  log2FC_sv <- rowMeans(vst_mat[, samples_grp2, drop = FALSE]) -
    rowMeans(vst_mat[, samples_grp1, drop = FALSE])
  
  baseMean_sv <- rowMeans(vst_mat[, c(samples_grp1, samples_grp2), drop = FALSE])
  df <- as.data.frame(res)
  df$log2FoldChange <- log2FC_sv
  df$logBaseMean <- log2(baseMean_sv + 1)
  df$Significant <- !is.na(df$padj) & df$padj < alpha
  
  ggplot(df, aes(x = logBaseMean, y = log2FoldChange, color = Significant)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("TRUE" = "firebrick", "FALSE" = "grey70")) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    labs(title = title, x = "log2(SV-regressed mean + 1)", y = "log2 Fold Change (SV-adjusted)") +
    theme_minimal()
}
# -------------------------  DEG Summary by Time -------------------------

summarise_deg_counts <- function(counts, metadata, group1, group2, alpha = 0.05, sv_terms = NULL) {
  out_list <- list()
  
  # Global comparison
  run_global <- run_deseq2_pipeline(
    counts, metadata,
    groups_keep = c(group1, group2),
    contrast = c("Group", group1, group2),
    out_prefix = paste0(group1, "_vs_", group2, "_AllTimes"),
    sv_terms = sv_terms
  )
  out_list$Global <- list(counts = count_degs(run_global$res, alpha),
                          res = run_global$res)
  
  # Day vs Night
  day_points   <- c("zt_2","zt_6","zt_10")
  night_points <- c("zt_14","zt_18","zt_22")
  
  run_day <- run_deseq2_pipeline(
    counts, metadata,
    groups_keep = c(group1, group2),
    times_keep = day_points,
    contrast = c("Group", group1, group2),
    out_prefix = paste0(group1, "_vs_", group2, "_Day"),
    sv_terms = sv_terms
  )
  out_list$Day <- list(counts = count_degs(run_day$res, alpha),
                       res = run_day$res)
  
  run_night <- run_deseq2_pipeline(
    counts, metadata,
    groups_keep = c(group1, group2),
    times_keep = night_points,
    contrast = c("Group", group1, group2),
    out_prefix = paste0(group1, "_vs_", group2, "_Night"),
    sv_terms = sv_terms
  )
  out_list$Night <- list(counts = count_degs(run_night$res, alpha),
                         res = run_night$res)
  
  # Each timepoint
  tp_list <- list()
  for (tp in levels(metadata$Timepoint)) {
    if (tp %in% c(day_points, night_points)) {
      meta_tp <- droplevels(metadata[metadata$Timepoint == tp & metadata$Group %in% c(group1, group2), ])
      counts_tp <- counts[, meta_tp$Filename, drop = FALSE]
      
      run_tp <- run_deseq2_pipeline(
        counts_tp, meta_tp,
        groups_keep = c(group1, group2),
        contrast = c("Group", group1, group2),
        out_prefix = paste0(group1, "_vs_", group2, "_", tp),
        sv_terms = sv_terms
      )
      tp_list[[tp]] <- list(counts = count_degs(run_tp$res, alpha),
                            res = run_tp$res)
    }
  }
  out_list$Timepoints <- tp_list
  
  return(out_list)
}


# -------------------------  Plot DEG Summary -------------------------

plot_deg_summary <- function(deg_summary, group1, group2) {
  # Convert timepoint names into pretty labels
  tp_labels <- names(deg_summary$Timepoints)
  tp_labels_pretty <- gsub("zt_(\\d+)", "ZT \\1", tp_labels)  # "zt_2" -> "ZT 2"
  
  # Extract counts correctly
  up_counts   <- vapply(deg_summary$Timepoints, function(x) x$counts$up, numeric(1))
  down_counts <- vapply(deg_summary$Timepoints, function(x) x$counts$down, numeric(1))
  
  # Flatten results into dataframe
  df <- data.frame(
    Category = c("Global", "Day", "Night", tp_labels_pretty),
    Up   = c(deg_summary$Global$counts$up, deg_summary$Day$counts$up, deg_summary$Night$counts$up, up_counts),
    Down = c(deg_summary$Global$counts$down, deg_summary$Day$counts$down, deg_summary$Night$counts$down, down_counts)
  )
  df$Total <- df$Up + df$Down
  
  # Desired order
  category_order <- c("Global", "Day", "Night", 
                      paste0("ZT ", c(2, 6, 10, 14, 18, 22)))
  
  # Apply factor ordering
  df$Category <- factor(df$Category, levels = category_order)
  
  # Melt for plotting
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

#--------------Heatmap---------------


plot_deg_heatmap <- function(vst_mat, res, metadata, top_n = 50,
                             color_by = "Group", shape_by = "Timepoint",
                             scale_rows = TRUE, title = "Top DEGs Heatmap") {
  
  # --- 0. Add SampleID (Group + Timepoint) if not present ---
  if (!"SampleID" %in% colnames(metadata)) {
    metadata$SampleID <- paste(metadata$Group, metadata$Timepoint, sep = "_")
  }
  
  # --- 1. Select top N DEGs by padj ---
  degs <- res[!is.na(res$padj), ]
  top_degs <- head(degs[order(degs$padj), ], n = top_n)
  top_genes <- rownames(top_degs)  # Ensembl IDs
  
  # --- 2. Subset vst_mat ---
  mat_sub <- vst_mat[top_genes, metadata$Filename, drop = FALSE]
  
  # --- 3. Replace Ensembl IDs with GeneSymbols if available ---
  if ("GeneSymbol" %in% colnames(res)) {
    gene_map <- data.frame(
      Ensembl = rownames(res),
      Symbol  = res$GeneSymbol
    )
    
    # Ensure mapping covers all top genes
    symbols <- gene_map$Symbol[match(top_genes, gene_map$Ensembl)]
    
    # Fallback: keep Ensembl if Symbol is missing
    symbols[is.na(symbols) | symbols == ""] <- top_genes[is.na(symbols) | symbols == ""]
    
    rownames(mat_sub) <- symbols
  }
  
  # --- 4. Replace sample filenames with SampleIDs ---
  colnames(mat_sub) <- metadata$SampleID[match(colnames(mat_sub), metadata$Filename)]
  
  # --- 5. Scale by row (z-score) if requested ---
  if (scale_rows) {
    mat_sub <- t(scale(t(mat_sub)))
  }
  
  # --- 6. Ensure factors ---
  metadata[[color_by]] <- factor(metadata[[color_by]])
  metadata[[shape_by]] <- factor(metadata[[shape_by]])
  
  # --- 7. Prepare named colors dynamically ---
  annotation_colors <- list()
  annotation_colors[[color_by]] <- setNames(
    colorRampPalette(brewer.pal(9, "Set1"))(nlevels(metadata[[color_by]])),
    levels(metadata[[color_by]])
  )
  annotation_colors[[shape_by]] <- setNames(
    viridis::viridis(nlevels(metadata[[shape_by]])),
    levels(metadata[[shape_by]])
  )
  
  # --- 8. Prepare annotations ---
  ha <- HeatmapAnnotation(
    df = metadata[, c(color_by, shape_by), drop = FALSE],
    col = annotation_colors
  )
  
  # --- 9. Plot heatmap ---
  Heatmap(
    mat_sub,
    name = ifelse(scale_rows, "z-score", "VST"),
    top_annotation = ha,
    show_row_names = TRUE,
    show_column_names = FALSE,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    column_title = title,
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 8)
  )
}

# ------------------------- Run GO enrichment on DEGs with barplot -------------------------
run_go_analysis <- function(res, 
                            alpha = 0.05, 
                            lfc_threshold = 1,
                            ontology = "BP", 
                            top_n = 20, 
                            title = "GO Enrichment",
                            sv_terms = sv_terms) {
  
  # --- 1. Filter DEGs ---
  sig <- as.data.frame(res)
  sig <- sig[!is.na(sig$padj) & sig$padj < alpha & abs(sig$log2FoldChange) > lfc_threshold, ]
  
  # --- 2. Extract gene symbols ---
  genes <- unique(sig$GeneSymbol[!is.na(sig$GeneSymbol) & sig$GeneSymbol != ""])
  
  if (length(genes) < 5) {
    message("⚠️ Not enough DEGs for GO analysis")
    return(NULL)
  }
  
  # --- 3. Define background (all tested genes) ---
  background_symbols <- unique(res$GeneSymbol[!is.na(res$GeneSymbol) & res$GeneSymbol != ""])
  
  # --- 4. Map SYMBOL → ENTREZ IDs ---
  entrez_ids <- mapIds(org.Mm.eg.db,
                       keys = genes,
                       column = "ENTREZID",
                       keytype = "SYMBOL",
                       multiVals = "first") %>% na.omit()
  
  bg_entrez <- mapIds(org.Mm.eg.db,
                      keys = background_symbols,
                      column = "ENTREZID",
                      keytype = "SYMBOL",
                      multiVals = "first") %>% na.omit()
  
  # --- 5. Run enrichment ---
  ego <- enrichGO(gene          = entrez_ids,
                  universe      = bg_entrez,
                  OrgDb         = org.Mm.eg.db,
                  keyType       = "ENTREZID",
                  ont           = ontology,
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.2)
  
  # --- 6. Plot results ---
  if (nrow(as.data.frame(ego)) > 0) {
    p <- barplot(ego, showCategory = top_n, 
                 title = paste(title, "(", ontology, ")"))
    print(p)
  } else {
    message("⚠️ No significant GO terms found")
  }
  
  return(ego)
}


#######-----------USAGE-------

# 1. Overall analysis (all groups, all times)
#run1 <- run_deseq2_pipeline(combined_counts, metadata_ordered, 
#out_prefix = "AllGroups_AllTimes")

# 2. Pairwise group comparison across all times (Control vs C)
#run2 <- run_deseq2_pipeline(combined_counts, metadata_ordered,
#groups_keep = c("Control", "CMVS"),
#contrast = c("Group", "Control", "CMVS"),
#out_prefix = "Control_vs_C_AllTimes")

# 3. Time-specific analysis (ZT6 only, Control vs E)
#run3 <- run_deseq2_pipeline(combined_counts, metadata_ordered,
#groups_keep = c("Control", "CSDS"),
#times_keep = c("zt_6"),
#contrast = c("Group", "Control", "CSDS"),
#out_prefix = "Control_vs_E_ZT6")


# pick how many surrogate variables to keep
n_sv_keep <- min(3, svseq$n.sv)
sv_terms <- paste0("SV", 1:n_sv_keep)

# ---- Define contrast and group labels once ----
contrast <- c("Group", "CSDS", "Control")
group1 <- contrast[2]
group2 <- contrast[3]
comparison_label <- paste(group1, "vs", group2)

# ---- Main DESeq2 pipeline ----
Deseq2object_sva <- run_deseq2_pipeline(
  combined_counts, metadata_ordered,
  groups_keep = c(group1, group2),
  contrast = contrast,
  out_prefix = paste0(comparison_label, "_AllTimes_SVA"),
  sv_terms = sv_terms
)

# PCA
plot_pca_subset(Deseq2object_sva$vst_mat, Deseq2object_sva$metadata,
                ntop = 1000, pc_x = 1, pc_y = 2)

# Summarise DEGs 
deg_summary <- summarise_deg_counts(
  combined_counts, metadata_ordered,
  group1, group2,
  sv_terms = sv_terms
)

# DEG summary plot
plot_deg_summary(deg_summary, group1, group2)

# Volcano
plot_volcano_sv(Deseq2object_sva, title = paste(comparison_label, "(SV-adjusted)"))

# MA plot
p_ma_sv <- plot_MA_sv(Deseq2object_sva, title = paste(comparison_label, "(SV-adjusted)"))
print(p_ma_sv)

# Heatmap
plot_deg_heatmap(Deseq2object_sva$vst_mat, Deseq2object_sva$res, Deseq2object_sva$metadata,
                 top_n = 50,
                 color_by = "Group",
                 shape_by = "Timepoint",
                 title = paste(comparison_label, ": Top 50 DEGs"))

# GO analyses
go_global <- run_go_analysis(deg_summary$Global$res, 
                             title = paste("GO Enrichment:", comparison_label, "(All Times)"))
go_day    <- run_go_analysis(deg_summary$Day$res, 
                             title = paste("GO Enrichment:", comparison_label, "(Day)"))
go_night  <- run_go_analysis(deg_summary$Night$res, 
                             title = paste("GO Enrichment:", comparison_label, "(Night)"))

# Timepoints
for (tp in names(deg_summary$Timepoints)) {
  run_go_analysis(deg_summary$Timepoints[[tp]]$res,
                  title = paste("GO Enrichment:", comparison_label, toupper(tp)))
}

# ------------------------- Use to diagnose outliers -------------------------

library(ggplot2)
library(matrixStats)

# Use rlog matrix from your DESeq2 run
vst_matrix <- Deseq2object_sva$vst_mat
metadata <- Deseq2object_sva$metadata

# Compute variance for top genes
rv <- rowVars(vst_matrix)
top_genes <- order(rv, decreasing = TRUE)[1:min(500, length(rv))] # top 500 variable genes
rlog_top <- vst_matrix[top_genes, ]

# PCA
pca <- prcomp(t(rlog_top), center = TRUE, scale. = FALSE)
var_expl <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 1)

# Prepare data frame for plotting
pca_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  Filename = metadata$Filename,
  Group = metadata$Group,
  Timepoint = metadata$Timepoint
)

# Plot PCA with sample labels
ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, shape = Timepoint)) +
  geom_point(size = 3) +
  geom_text(aes(label = Filename), vjust = -0.5, size = 2.5, check_overlap = TRUE) +
  labs(
    title = "Quick PCA with rlog (Top 500 variable genes)",
    x = paste0("PC1 (", var_expl[1], "%)"),
    y = paste0("PC2 (", var_expl[2], "%)")
  ) +
  theme_minimal()
