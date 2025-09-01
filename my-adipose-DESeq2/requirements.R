# Install CRAN + Bioconductor packages in one go
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

cran_pkgs <- c(
  "readr","dplyr","RColorBrewer","viridis","readxl",
  "ggplot2","gtools","matrixStats","circlize","ComplexHeatmap",
  "MetaCycle"
)

bio_pkgs <- c(
  "DESeq2","biomaRt","limma","sva",
  "clusterProfiler","org.Mm.eg.db","enrichplot"
)

for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, dependencies = TRUE)
  }
}

for (p in bio_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    BiocManager::install(p, ask = FALSE, update = FALSE)
  }
}

message("\nAll set. Session info:\n")
print(sessionInfo())
