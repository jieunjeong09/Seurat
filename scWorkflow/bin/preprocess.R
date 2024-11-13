library(Seurat)

args <- commandArgs(trailingOnly=TRUE)
input_dir <- args[1]
output_file <- args[2]

# Load and preprocess data
data <- Read10X(input.dir = input_dir)
seurat_obj <- CreateSeuratObject(counts = data)

# add filtering steps as needed

# Perform SCTransform normalization
seurat_obj <- SCTransform(seurat_obj)
saveRDS(seurat_obj, file = output_file)
