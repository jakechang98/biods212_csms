# plotting UMAPs of initial and new subtypes
library(data.table)
library(dplyr)
library(ggplot2)
library(Seurat)

# Load scRNA-seq data and prepare for input
data <- fread("Desktop/Stanford/Winter2024/STEMREM/Project/scRNAseq_2/GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt", data.table = FALSE)
rownames(data) <- data$Index
data$Index <- NULL

# Load annotation data and prepare for input
annotation <- read.csv("/Users/savagyan/Desktop/Stanford/Winter2024/STEMREM/Project/scRNAseq_2/GSE132465_GEO_processed_CRC_10X_cell_annotation.txt", sep = "\t", row.names = 1)
new_clusters <- read.csv("/Users/savagyan/Desktop/Stanford/Winter2024/STEMREM/Project/scRNAseq_2/GSE132465_GEO_processed_CRC_10X_ct2_1d_type.txt", sep = "\t", row.names = 1)
annotation_tumor <- annotation[annotation$Class == "Tumor",]

# filtering the scRNA seq to only tumor patient samples
tumor_patients <- c(rownames(annotation[annotation$Class == "Tumor",]))  
data_tumor <- data[, tumor_patients]

new_clusters <- new_clusters[rownames(annotation_tumor),]
annotation_tumor$Cell_cluster <- new_clusters

# creating the Seurat object for visualization
seurat <- CreateSeuratObject(counts = as.matrix(data_tumor), project = "nn")
seurat <- NormalizeData(seurat, normalization.method = "LogNormalize", scale.factor = 10000)
seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)
seurat <- ScaleData(seurat, features = rownames(seurat))
seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat))
seurat <- RunUMAP(seurat, dims = 1:pc_dims, seed.use = seed)

# add annotations to metadata
seurat <- AddMetaData(object = seurat, metadata = annotation_tumor)


# plot UMAP
# initial cell types (major)
DimPlot(seurat, reduction = "umap", group.by = "Cell_type", label = FALSE) +
  xlab("UMAP1") + ylab("UMAP2") + ggtitle("Major Cell Types") +
  theme(legend.text = element_text(size = 8), legend.title = element_text(size = 12),
        axis.text = element_text(size = 10), axis.title = element_text(size = 10),            
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5, margin = margin(b = 20)),
        aspect.ratio = 1)

DimPlot(seurat, reduction = "umap", group.by = "Cell_cluster", label = FALSE) +
  xlab("UMAP1") + ylab("UMAP2") + ggtitle("Potency-informed Cell Types") +
  theme(legend.text = element_text(size = 8), legend.title = element_text(size = 12),
        axis.text = element_text(size = 10), axis.title = element_text(size = 10),            
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5, margin = margin(b = 20)),
        aspect.ratio = 1)

