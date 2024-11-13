library(Seurat)
library(ggplot2)

args <- commandArgs(trailingOnly=TRUE)
input_file <- args[1]
output_dir <- args[2]

# Load data
seurat_obj <- readRDS(input_file)

# Ensure output directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# PCA Plot
DimPlot(seurat_obj, reduction = "pca") + ggtitle("PCA Plot")
ggsave(file.path(output_dir, "pca_plot.pdf"))

# Feature Plot
# FeaturePlot(seurat_obj, features = c("gene1", "gene2"))
# ggsave(file.path(output_dir, "feature_plot.pdf"))

# Export Cluster Assignment and PCA plot coordinates
meta_extract <- data.frame(
  Cell = colnames(seurat_obj),
  Cluster = Idents(seurat_obj),
  PC1 = seurat_obj@reductions$pca@cell.embeddings[, 1],  # First dimension of PCA
  PC2 = seurat_obj@reductions$pca@cell.embeddings[, 2]   # Second dimension of PCA
)
write.csv(meta_extract, file.path(output_dir, "meta_extract.csv"), row.names = FALSE)
