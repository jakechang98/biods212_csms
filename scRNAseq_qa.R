# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #
# Code to process scRNA-seq for downstream analysis
#
#
# Author: Jake Chang
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #

library(BiocManager)
BiocManager::install()

library(SingleCellExperiment)
library(data.table)
library(ggplot2)
library(scater)
library(scran)
library(scuttle)

plotting_dir <- file.path(getwd(), "plots")

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Read files and create SingleCellExperiment object
cell_annotation_fp <- "data/GSE132465_GEO_processed_CRC_10X_cell_annotation.txt"
counts_matrix_fp <- "data/GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt"

cell_annotation <- read.delim(cell_annotation_fp)
rownames(cell_annotation) <- cell_annotation$Index

counts_matrix <- fread(counts_matrix_fp, showProgress = FALSE)
counts_matrix <- as.data.frame(counts_matrix)
rownames(counts_matrix) <- counts_matrix$Index


sce <- SingleCellExperiment(
  list(counts=counts_matrix[, -1]),
  colData=cell_annotation,
  metadata=list(study="GSE132465")
)
rownames(sce) <- counts_matrix$Index
colnames(sce) <- colnames(counts_matrix)[-1]

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Exploratory plotting
p1 <- ggplot(cell_annotation, aes(x = Patient, fill = Cell_subtype)) +
  geom_bar(position = "stack") +
  labs(title = "Cell subtypes by donor",
       x = "Donors",
       y = "Count",
       fill = "Cell Subtype") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p1 + facet_wrap(~ Class)
ggsave(file.path(plotting_dir, "cell_subtypes_per_donor.png"),
       plot = last_plot(),
       width = 15,
       height = 8)


p2 <- ggplot(cell_annotation, aes(x = Patient, fill = Cell_type)) +
  geom_bar(position = "stack") +
  labs(title = "Cell types by donor",
       x = "Donors",
       y = "Count",
       fill = "Cell Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p2 + facet_wrap(~ Class)
ggsave(file.path(plotting_dir, "cell_types_per_donor.png"),
       plot = last_plot(),
       width = 15,
       height = 8)

rm(cell_annotation, counts_matrix)

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# lib.sf.sce <- librarySizeFactors(sce)
# hist(log10(lib.sf.sce), xlab="Log10[Size factor]", col='grey80')

sce <- logNormCounts(sce)
saveRDS(sce, file = "sce.rds")

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


