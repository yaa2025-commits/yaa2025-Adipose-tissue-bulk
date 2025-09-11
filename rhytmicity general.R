# -------------------------  1. Load Libraries -------------------------
library(DESeq2)
library(readr)
library(dplyr)
library(biomaRt)
library(RColorBrewer)
library(ggplot2)
library(limma)

# -------------------------  Set Working Directory -------------------
setwd("~/Documents/adiposebulkstress/vW")

# -------------------------  2. Connect to Ensembl ----------------------
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# -------------------------  3. Load Count Data -------------------------
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
combined_counts <- as.matrix(combined_counts)
rownames(combined_counts) <- gene_ids
colnames(combined_counts) <- sub("\\.txt$", "", count_files)

# Filter low counts
keep <- rowSums(combined_counts > 0) >= 8
combined_counts <- combined_counts[keep, ]

# -------------------------  4. Load Metadata ---------------------------
metadata <- read.csv("Analysis metadata_gn_vW.csv", header = TRUE)
metadata$Filename <- gsub("\\.txt$", "", metadata$Filename)
metadata_ordered <- metadata[match(colnames(combined_counts), metadata$Filename), ]
metadata_ordered$Batch <- factor(metadata_ordered$Batch)
metadata_ordered$Group <- factor(metadata_ordered$Group)
metadata_ordered$Timepoint <- factor(metadata_ordered$Timepoint)

# -------------------------  5. DESeq2 & VST ----------------------------
metadata_ordered$Group <- relevel(metadata_ordered$Group, ref = "Control")

dds <- DESeqDataSetFromMatrix(countData = combined_counts,
                              colData = metadata_ordered,
                              design = ~ Group)
dds <- DESeq(dds)

vst_transformed <- vst(dds, blind = FALSE)
vst_matrix <- assay(vst_transformed)

# Batch correction
batch <- metadata_ordered$Batch
adjusted_vst <- removeBatchEffect(vst_matrix, batch = batch)

# Gene mapping
genes <- rownames(adjusted_vst)
gene_info <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                   filters = 'ensembl_gene_id',
                   values = genes,
                   mart = mart)

# -------------------------  6. Reorder Expression Matrix ---------------
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

# -------------------------  Plot gene rhythm function - single gene ---------------

plot_gene_rhythm <- function(gene_symbol, 
                             expr_matrix = expr_mat, 
                             metadata = metadata_ordered,
                             gene_map = gene_info) {
  ens_ids <- gene_map$ensembl_gene_id[gene_map$external_gene_name == gene_symbol]
  if (length(ens_ids) == 0) stop(paste0("Gene symbol '", gene_symbol, "' not found."))
  ens_id <- ens_ids[1]
  if (!(ens_id %in% rownames(expr_matrix))) stop(paste0("Ensembl ID '", ens_id, "' not in matrix."))
  
  meta <- metadata
  meta$Expression <- as.numeric(expr_matrix[ens_id, meta$SampleLabel])
  
  summary_data <- meta %>%
    group_by(Group, Timepoint) %>%
    summarise(Mean = mean(Expression, na.rm = TRUE),
              SD = sd(Expression, na.rm = TRUE), .groups = "drop")
  
  ggplot(summary_data, aes(x = Timepoint, y = Mean, color = Group, group = Group)) +
    geom_line(size = 1.2) +
    geom_point(size = 2) +
    geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD, fill = Group),
                alpha = 0.2, color = NA) +
    scale_x_continuous(breaks = sort(unique(summary_data$Timepoint))) +
    labs(title = gene_symbol,
         x = "Zeitgeber Time (ZT)",
         y = "VST Expression") +
    theme_minimal(base_size = 14) +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2")
}

# -------------------------  Plot gene rhythm function - core clock genes ---------------


plot_core_clock_rhythms <- function(genes = c("Clock", "Bmal1", "Nr1d1", "Nr1d2", "Dbp",
                                              "Per1", "Per2", "Per3", "Cry1", "Cry2", "Rora", "Lep"),
                                    expr_matrix = expr_mat,
                                    metadata = metadata_ordered,
                                    gene_map = gene_info) {
  all_data <- list()
  
  for (gene_symbol in genes) {
    ens_ids <- gene_map$ensembl_gene_id[gene_map$external_gene_name == gene_symbol]
    if (length(ens_ids) == 0) next
    ens_id <- ens_ids[1]
    if (!(ens_id %in% rownames(expr_matrix))) next
    
    meta <- metadata
    meta$Expression <- as.numeric(expr_matrix[ens_id, meta$SampleLabel])
    meta$Gene <- gene_symbol
    
    all_data[[gene_symbol]] <- meta
  }
  
  plot_data <- bind_rows(all_data)
  
  summary_data <- plot_data %>%
    group_by(Gene, Group, Timepoint) %>%
    summarise(Mean = mean(Expression, na.rm = TRUE),
              SD = sd(Expression, na.rm = TRUE), .groups = "drop")
  
  ggplot(summary_data, aes(x = Timepoint, y = Mean, color = Group, group = Group)) +
    geom_line(size = 1) +
    geom_point(size = 1.5) +
    geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD, fill = Group),
                alpha = 0.15, color = NA) +
    facet_wrap(~ Gene, scales = "free_y", ncol = 3) +
    scale_x_continuous(breaks = sort(unique(summary_data$Timepoint))) +
    labs(title = "Core Clock Gene Rhythms",
         x = "Zeitgeber Time (ZT)",
         y = "VST Expression") +
    theme_minimal(base_size = 14) +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2")
}

# Single gene
plot_gene_rhythm("Ucp1")

# Core clock panel
plot_core_clock_rhythms()


