# Example runner for common comparisons
source("analysis_pipeline.R")

# Example 1: CSDS vs CMVS (all times) with SVA covariates
n_sv_keep <- min(3, sva_info$n_sv)
sv_terms  <- paste0("SV", 1:n_sv_keep)

ex1 <- run_deseq2_pipeline(
  combined_counts, metadata_ordered,
  groups_keep = c("CSDS", "CMVS"),
  contrast = c("Group", "CSDS", "CMVS"),
  out_prefix = "CSDS_vs_CMVS_AllTimes_SVA",
  sv_terms = sv_terms
)

pdf("results/plots/pca_CSDS_vs_CMVS.pdf")
print(plot_pca_subset(ex1$vst_mat, ex1$metadata, ntop = 1000, pc_x = 1, pc_y = 2))
dev.off()

pdf("results/plots/volcano_CSDS_vs_CMVS.pdf")
print(plot_volcano(ex1$res, title = "CSDS vs CMVS (All Times)"))
dev.off()

pdf("results/plots/heatmap_top50_CSDS_vs_CMVS.pdf", width = 8, height = 10)
ht <- plot_deg_heatmap(ex1$vst_mat, ex1$res, ex1$metadata,
                       top_n = 50, title = "CSDS vs CMVS: Top 50 DEGs")
ComplexHeatmap::draw(ht)
dev.off()

# GO analysis
go1 <- run_go_analysis(ex1$res, title = "GO: CSDS vs CMVS")
if (!is.null(go1)) {
  openxlsx::write.xlsx(as.data.frame(go1), "results/go/GO_CSDS_vs_CMVS.xlsx")
}
