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
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)

# -------------------------  2. Set Working Directory -------------------
# CHANGE FILE PATH - search CHANGE to see all the places u need to change the path in
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
# CHANGE FILE PATH
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

# Relevel first
metadata_ordered$Group <- relevel(metadata_ordered$Group, ref = "Control")

# Create and run DESeq
dds <- DESeqDataSetFromMatrix(countData = combined_counts, colData = metadata_ordered, design = ~ Group)
dds <- DESeq(dds)

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

batch <- metadata_ordered$Batch
adjusted_vst <- removeBatchEffect(vst_matrix, batch = batch)

vst_gene_symbols <- gene_info$external_gene_name[match(rownames(adjusted_vst), gene_info$ensembl_gene_id)]

#plotPCA(dds_transformed)

# -------------------------  7. Reorder Expression Matrix ---------------
metadata_ordered$Timepoint <- gsub("zt_", "", tolower(trimws(as.character(metadata_ordered$Timepoint))))
metadata_ordered$Timepoint <- as.numeric(metadata_ordered$Timepoint)

col_order <- order(metadata_ordered$Timepoint, metadata_ordered$Group)
expr_mat <- adjusted_vst[, col_order]
metadata_ordered <- metadata_ordered[col_order, ]

metadata_ordered$SampleLabel <- paste(metadata_ordered$Group,
                                      metadata_ordered$Timepoint,
                                      ave(seq_along(metadata_ordered$Group),
                                          interaction(metadata_ordered$Group, metadata_ordered$Timepoint),
                                          FUN = seq_along),
                                      sep = "_")
colnames(expr_mat) <- metadata_ordered$SampleLabel

stopifnot(ncol(expr_mat) == nrow(metadata_ordered))


# -------------------------  3-method Rhythmicity Analysis -------------------------
run_3method_rhythm <- function(expr_matrix, metadata, group_name = NULL, period = 24, amp_cutoff = 0.3) {
  # --- choose sample column ---
  if ("SampleLabel" %in% colnames(metadata)) {
    sample_col <- "SampleLabel"
  } else if ("Filename" %in% colnames(metadata)) {
    sample_col <- "Filename"
  } else if ("SampleID" %in% colnames(metadata)) {
    sample_col <- "SampleID"
  } else {
    stop("metadata must contain one of: SampleLabel, Filename, or SampleID")
  }
  sample_labels <- as.character(metadata[[sample_col]])
  times_raw <- metadata$Timepoint
  times <- as.numeric(as.character(times_raw))
  
  # sanity check
  if (!all(sample_labels %in% colnames(expr_matrix))) {
    missing_samps <- sample_labels[!sample_labels %in% colnames(expr_matrix)]
    stop("Missing samples in expr_matrix: ", paste(head(missing_samps, 5), collapse = ", "))
  }
  mat <- expr_matrix[, sample_labels, drop = FALSE]
  stopifnot(all(colnames(mat) == sample_labels))
  
  # temp files for JTK & ARSER
  tmp_dir <- tempdir()
  infile <- file.path(tmp_dir, paste0("Meta_infile_", group_name, ".txt"))
  infile_avg <- file.path(tmp_dir, paste0("Meta_infile_avg_", group_name, ".txt"))
  outdir <- file.path(tmp_dir, paste0("MetaCycle_out_", group_name))
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  # write replicate matrix for JTK
  write.table(data.frame(CycID = rownames(mat), mat, check.names = FALSE),
              file = infile, sep = "\t", row.names = FALSE, quote = FALSE
  )
  
  # collapse replicates for ARSER & Cosinor
  unique_times_sorted <- sort(unique(times))
  avg_expr <- sapply(unique_times_sorted, function(tp) {
    cols <- which(times == tp)
    rowMeans(mat[, cols, drop = FALSE], na.rm = TRUE)
  })
  colnames(avg_expr) <- unique_times_sorted
  avg_expr <- as.matrix(avg_expr)
  
  write.table(data.frame(CycID = rownames(avg_expr), avg_expr, check.names = FALSE),
              file = infile_avg, sep = "\t", row.names = FALSE, quote = FALSE
  )
  
  combined <- data.frame(CycID = rownames(mat), stringsAsFactors = FALSE)
  
  # -------------------------
  # Run JTK + ARSER
  # -------------------------
  meta2d(
    infile = infile, filestyle = "txt", timepoints = times,
    outdir = outdir, minper = period, maxper = period,
    cycMethod = "JTK", outIntegration = "noIntegration"
  )
  
  meta2d(
    infile = infile_avg, filestyle = "txt", timepoints = unique_times_sorted,
    outdir = outdir, minper = period, maxper = period,
    cycMethod = "ARS", outIntegration = "noIntegration", ARSdefaultPer = period
  )
  
  files <- list.files(outdir, full.names = TRUE)
  jtk_file <- files[grepl("JTKresult", files, ignore.case = TRUE)]
  arser_file <- files[grepl("ARSresult", files, ignore.case = TRUE)]
  
  # JTK
  if (length(jtk_file) > 0) {
    jtk_res <- read.delim(jtk_file[1], stringsAsFactors = FALSE, check.names = FALSE)
    phase_name_jtk <- names(jtk_res)[grepl("phase|lag|peak", names(jtk_res), ignore.case = TRUE)]
    adj_name_jtk <- names(jtk_res)[tolower(names(jtk_res)) %in% c("bh.q", "adjp", "qvalue", "adj.p")]
    
    JTK_FDR <- if (length(adj_name_jtk) > 0) as.numeric(jtk_res[[adj_name_jtk[1]]]) else NA
    JTK_phase <- if (length(phase_name_jtk) > 0) as.numeric(jtk_res[[phase_name_jtk[1]]]) else NA
    
    jtk_map <- data.frame(CycID = as.character(jtk_res[[1]]),
                          JTK_FDR, JTK_phase,
                          stringsAsFactors = FALSE
    )
    combined <- left_join(combined, jtk_map, by = "CycID")
  } else {
    combined$JTK_FDR <- NA
    combined$JTK_phase <- NA
  }
  
  # ARSER
  if (length(arser_file) > 0) {
    arser_res <- read.delim(arser_file[1], stringsAsFactors = FALSE, check.names = FALSE)
    phase_name_ars <- names(arser_res)[grepl("phase|lag|peak", names(arser_res), ignore.case = TRUE)]
    adj_name_ars <- names(arser_res)[tolower(names(arser_res)) %in% c("bh.q", "adjp", "qvalue", "adj.p", "pvalue")]
    
    ARSER_q <- if (length(adj_name_ars) > 0) as.numeric(arser_res[[adj_name_ars[1]]]) else NA
    ARSER_phase <- if (length(phase_name_ars) > 0) as.numeric(arser_res[[phase_name_ars[1]]]) else NA
    
    arser_map <- data.frame(CycID = as.character(arser_res[[1]]), ARSER_q, ARSER_phase, stringsAsFactors = FALSE)
    combined <- left_join(combined, arser_map, by = "CycID")
  } else {
    combined$ARSER_q <- NA
    combined$ARSER_phase <- NA
  }
  
  # -------------------------
  # Robust Cosinor
  # -------------------------
  cosinor_results <- lapply(rownames(avg_expr), function(gene_id) {
    gene_data <- data.frame(
      Time = as.numeric(colnames(avg_expr)),
      Expression = as.numeric(avg_expr[gene_id, ])
    )
    if (length(unique(gene_data$Time)) < 3 || sd(gene_data$Expression) == 0) {
      return(data.frame(CycID = gene_id, COS_p = NA, COS_phase = NA))
    }
    tryCatch(
      {
        fit <- lm(Expression ~ cos(2 * pi * Time / period) + sin(2 * pi * Time / period), data = gene_data)
        fit0 <- lm(Expression ~ 1, data = gene_data)
        an <- tryCatch(anova(fit0, fit), error = function(e) NULL)
        pval <- if (!is.null(an) && nrow(an) >= 2) an$`Pr(>F)`[2] else NA
        coefs <- coef(fit)
        a <- coefs[2]
        b <- coefs[3]
        phase_hours <- ((atan2(-b, a) %% (2 * pi)) * period) / (2 * pi)
        data.frame(CycID = gene_id, COS_p = pval, COS_phase = phase_hours)
      },
      error = function(e) {
        data.frame(CycID = gene_id, COS_p = NA, COS_phase = NA)
      }
    )
  })
  cosinor_df <- do.call(rbind, cosinor_results)
  combined <- left_join(combined, cosinor_df, by = "CycID")
  
  # New Amplitude and Log2FC calculation
  # This returns the simple difference between max and min VST expression.
  amps <- apply(avg_expr, 1, function(x) {
    mx <- max(x, na.rm = TRUE)
    mn <- min(x, na.rm = TRUE)
    mx - mn
  })
  
  # We'll use this new amplitude value for our filter.
  # The `Daily_log2FC` variable can now be removed.
  Daily_log2FC <- amps
  
  # --- Rhythmicity Summary & Filters ---
  combined <- combined %>%
    left_join(data.frame(CycID = names(amps), VST_Amp = amps), by = "CycID") %>%
    dplyr::as_tibble() %>%
    mutate(
      Phase = ifelse(!is.na(JTK_phase), JTK_phase, ARSER_phase),
      JTK_hit = !is.na(JTK_FDR) & JTK_FDR < 0.2,
      ARSER_hit = !is.na(ARSER_q) & ARSER_q < 0.1,
      COS_hit = !is.na(COS_p) & COS_p < 0.05,
      # The Rhythmic filter now uses the simpler VST amplitude cutoff.
      Rhythmic = rowSums(cbind(JTK_hit, ARSER_hit, COS_hit), na.rm = TRUE) >= 2 & VST_Amp >= amp_cutoff
    ) %>%
    dplyr::select(-JTK_hit, -ARSER_hit, -COS_hit) # Remove temporary columns
  
  # add group suffix
  if (!is.null(group_name)) {
    suffix <- paste0("_", group_name)
    names(combined)[-1] <- paste0(names(combined)[-1], suffix)
  }
  
  return(combined)
}

# The loop and merging code remains the same, as the problem was in the function definition.
all_rhythm <- list()
for (grp in unique(metadata_ordered$Group)) {
  cat("Running for group:", grp, "\n")
  grp_meta <- metadata_ordered %>% filter(Group == grp)
  grp_meta <- grp_meta[order(as.numeric(as.character(grp_meta$Timepoint))), ]
  
  res_grp <- run_3method_rhythm(expr_mat, grp_meta, group_name = as.character(grp), period = 24, amp_cutoff = 0.3)
  all_rhythm[[as.character(grp)]] <- res_grp
}

# Merge by CycID
rhythm_compare <- Reduce(function(x, y) full_join(x, y, by = "CycID"), all_rhythm)
group_list <- unique(metadata_ordered$Group)

# ------------------------- Merge results across groups -------------------------
rhythm_compare <- Reduce(function(x, y) full_join(x, y, by = "CycID"), all_rhythm)
names(rhythm_compare) <- gsub("\\.x$", "_Control", names(rhythm_compare))
names(rhythm_compare) <- gsub("\\.y$", "_CMVS", names(rhythm_compare))

# ------------------------- Add Gene Symbols -------------------------
ens2symbol <- gene_info %>%
  filter(!duplicated(ensembl_gene_id)) %>%
  dplyr::select(CycID = ensembl_gene_id, GeneSymbol = external_gene_name)

rhythm_compare <- rhythm_compare %>%
  left_join(ens2symbol, by = "CycID") %>%
  mutate(GeneSymbol = ifelse(is.na(GeneSymbol), CycID, GeneSymbol))

# -------------------------  Total number of genes barplot ----


# Count rhythmic genes per group
rhythmic_counts <- rhythm_compare %>%
  summarise(
    Control = sum(Rhythmic_Control, na.rm = TRUE),
    CMVS    = sum(Rhythmic_CMVS, na.rm = TRUE),
    CSDS    = sum(Rhythmic_CSDS, na.rm = TRUE)
  ) %>%
  pivot_longer(cols = everything(), names_to = "Group", values_to = "Rhythmic_Genes")

# Reorder the Group column
rhythmic_counts$Group <- factor(rhythmic_counts$Group, levels = c("Control", "CMVS", "CSDS"))

# Create the plot object
p <- ggplot(rhythmic_counts, aes(x = Group, y = Rhythmic_Genes, fill = Group)) +
  geom_bar(stat = "identity") +
  labs(title = "Total Number of Rhythmic Genes per Group",
       x = "Group", y = "Number of Rhythmic Genes") +
  theme_minimal(base_size = 14) +
  scale_fill_brewer(palette = "Set2") +
  geom_text(aes(label = Rhythmic_Genes), vjust = -0.5)

# Print the plot to display it
print(p)

# -------------------------  9. Heatmaps of Top 50 Rhythmic Genes ----

plot_top50_heatmap <- function(group_name, rhythm_compare, expr_mat, metadata, gene_info) {
  message("Generating heatmap for group: ", group_name)
  
  rhythmic_col <- paste0("Rhythmic_", group_name)
  if (!rhythmic_col %in% colnames(rhythm_compare)) {
    stop("Column ", rhythmic_col, " not found in rhythm_compare. ",
         "Check run_3method_rhythm() output for group_name = '", group_name, "'.")
  }
  
  # --- filter rhythmic genes for this group
  rhythmic_genes <- rhythm_compare %>%
    filter(.data[[rhythmic_col]] == TRUE)
  
  if (nrow(rhythmic_genes) == 0) {
    warning("No rhythmic genes found for ", group_name, " — skipping heatmap")
    return(NULL)
  }
  
  # --- choose p-value columns
  jtk_col   <- paste0("JTK_FDR_", group_name)
  arser_col <- paste0("ARSER_q_", group_name)
  cos_col  <- paste0("COS_p_", group_name)
  
  # --- compute combined p (take min across methods)
  rhythmic_genes <- rhythmic_genes %>%
    mutate(
      combined_p = pmin(
        coalesce(.data[[jtk_col]],   1),
        coalesce(.data[[arser_col]], 1),
        coalesce(.data[[cos_col]], 1),  # use COS instead
        na.rm = TRUE
      )
    )
  
  # --- top 50 lowest combined p
  top_genes <- rhythmic_genes %>%
    arrange(combined_p) %>%
    slice_head(n = 50)
  
  gene_ids <- top_genes$CycID
  
  # --- subset expression
  sample_col <- if ("SampleLabel" %in% colnames(metadata)) "SampleLabel" else "Filename"
  meta_sub <- metadata %>% arrange(as.numeric(as.character(Timepoint)))
  sub_expr <- expr_mat[gene_ids, meta_sub[[sample_col]], drop = FALSE]
  
  # --- z-score per row
  sub_expr_z <- t(scale(t(sub_expr)))  # row-scale
  sub_expr_z[is.na(sub_expr_z)] <- 0   # safety
  
  # --- gene labels
  gene_labels <- top_genes %>%
    left_join(gene_info, by = c("CycID" = "ensembl_gene_id")) %>%
    mutate(label = ifelse(is.na(external_gene_name), CycID, external_gene_name)) %>%
    pull(label)
  
  # --- heatmap
  Heatmap(sub_expr_z,
          name = paste0("Zscore_", group_name),
          cluster_rows = TRUE, cluster_columns = FALSE,
          show_row_names = TRUE, show_column_names = FALSE,
          row_names_gp = gpar(fontsize = 6),
          column_split = factor(meta_sub$Timepoint, levels = sort(unique(meta_sub$Timepoint))),
          row_labels = gene_labels,
          col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")))
}



# Example usage
plot_top50_heatmap("Control", rhythm_compare, expr_mat, metadata_ordered, gene_info)



#------rhytmicity changes vs control

library(tidyverse)
phase_diff <- function(phase1, phase2) {
  diff <- abs(phase1 - phase2)
  ifelse(diff > 12, 24 - diff, diff)
}

rhythm_compare <- rhythm_compare %>%
  mutate(
    Status_CMVS = case_when(
      Rhythmic_Control == TRUE  & Rhythmic_CMVS == FALSE ~ "Lost Rhythm",
      Rhythmic_Control == FALSE & Rhythmic_CMVS == TRUE  ~ "Gained Rhythm",
      Rhythmic_Control == TRUE  & Rhythmic_CMVS == TRUE &
        (abs(log2(VST_Amp_CMVS / VST_Amp_Control)) > log2(1.5) |
           phase_diff(Phase_CMVS, Phase_Control) > 4) ~ "Altered Rhythm",
      TRUE ~ "No Change"
    ),
    Status_CSDS = case_when(
      Rhythmic_Control == TRUE  & Rhythmic_CSDS == FALSE ~ "Lost Rhythm",
      Rhythmic_Control == FALSE & Rhythmic_CSDS == TRUE  ~ "Gained Rhythm",
      Rhythmic_Control == TRUE  & Rhythmic_CSDS == TRUE &
        (abs(log2(VST_Amp_CSDS / VST_Amp_Control)) > log2(1.5) |
           phase_diff(Phase_CSDS, Phase_Control) > 4) ~ "Altered Rhythm",
      TRUE ~ "No Change"
    )
  )
#barplot 
status_long <- rhythm_compare %>%
  dplyr::select(CycID, Status_CMVS, Status_CSDS) %>%
  tidyr::pivot_longer(cols = starts_with("Status"),
               names_to = "Group", values_to = "Status") %>%
  mutate(Group = ifelse(Group == "Status_CMVS", "CMVS", "CSDS")) %>%
  filter(Status != "No Change")

ggplot(status_long, aes(x = Status, fill = Group)) +
  geom_bar(position = "dodge") +
  labs(title = "Rhythmicity Changes vs Control",
       x = "Change Type", y = "Number of Genes") +
  theme_minimal(base_size = 14) +
  scale_fill_brewer(palette = "Set1")

# ------------------------- GO Enrichment -------------------------
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)

run_go_plot <- function(gene_symbols, background_symbols, group_name, ont = "BP", showCategory = 20) {
  entrez_ids <- mapIds(org.Mm.eg.db,
                       keys = gene_symbols,
                       column = "ENTREZID",
                       keytype = "SYMBOL",
                       multiVals = "first") %>% na.omit()
  
  bg_entrez <- mapIds(org.Mm.eg.db,
                      keys = background_symbols,
                      column = "ENTREZID",
                      keytype = "SYMBOL",
                      multiVals = "first") %>% na.omit()
  
  ego <- enrichGO(gene = entrez_ids, universe = bg_entrez,
                  OrgDb = org.Mm.eg.db, keyType = "ENTREZID",
                  ont = ont, pAdjustMethod = "BH",
                  pvalueCutoff = 0.05, qvalueCutoff = 0.2)
  
  if (nrow(as.data.frame(ego)) == 0) {
    message("⚠️ No significant GO terms for ", group_name)
    return(NULL)
  }
  barplot(ego, showCategory = showCategory, title = paste(group_name, "GO Enrichment (", ont, ")"))
}

# Example: lost/gained/altered sets
genes_lost_CMVS    <- rhythm_compare %>% filter(Status_CMVS == "Lost Rhythm") %>% pull(GeneSymbol)
genes_gained_CMVS  <- rhythm_compare %>% filter(Status_CMVS == "Gained Rhythm") %>% pull(GeneSymbol)
genes_altered_CMVS <- rhythm_compare %>% filter(Status_CMVS == "Altered Rhythm") %>% pull(GeneSymbol)

background_genes <- rhythm_compare$GeneSymbol

run_go_plot(genes_lost_CMVS,    background_genes, "CMVS Lost Rhythm")
run_go_plot(genes_gained_CMVS,  background_genes, "CMVS Gained Rhythm")
run_go_plot(genes_altered_CMVS, background_genes, "CMVS Altered Rhythm")

genes_lost_CSDS    <- rhythm_compare %>% filter(Status_CSDS == "Lost Rhythm") %>% pull(GeneSymbol)
genes_gained_CSDS  <- rhythm_compare %>% filter(Status_CSDS == "Gained Rhythm") %>% pull(GeneSymbol)
genes_altered_CSDS <- rhythm_compare %>% filter(Status_CSDS == "Altered Rhythm") %>% pull(GeneSymbol)

run_go_plot(genes_lost_CSDS,    background_genes, "CSDS Lost Rhythm")
run_go_plot(genes_gained_CSDS,  background_genes, "CSDS Gained Rhythm")
run_go_plot(genes_altered_CSDS, background_genes, "CSDS Altered Rhythm")

