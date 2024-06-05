# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #
# Differential expression analysis performed on scRNA-seq cohort
# Normal vs. tumor
#
#
# Author: Jake Chang
# Code adapted from this tutorial:
# https://hbctraining.github.io/scRNA-seq/lessons/pseudobulk_DESeq2_scrnaseq.html
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #

# Load libraries
library(scater)
library(Seurat)
library(tidyverse)
library(cowplot)
library(Matrix.utils)
library(edgeR)
library(dplyr)
library(magrittr)

# install.packages("grr")
# install.packages("https://cran.r-project.org/src/contrib/Archive/Matrix.utils/Matrix.utils_0.9.8.tar.gz", type = "source", repos = NULL)
library(Matrix)

library(purrr)
library(reshape2)
library(S4Vectors)
library(tibble)
library(SingleCellExperiment)
library(pheatmap)

# BiocManager::install("apeglm")
library(apeglm)
library(png)

# BiocManager::install("DESeq2")
library(DESeq2)
library(RColorBrewer)
library(dplyr)

# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #
# Load SingleCellExperiment object
sce <- readRDS("sce.rds")

# Named vector of cluster names
kids <- purrr::set_names(levels(factor((sce$Cell_type))))
kids

# Total number of clusters
nk <- length(kids)
nk

# Named vector of sample names
sids <- purrr::set_names(levels(factor(sce$Sample)))

# Total number of samples
ns <- length(sids)
ns

# Generate sample level metadata

## Determine the number of cells per sample
table(sce$Sample)

## Turn named vector into a numeric vector of number of cells per sample
n_cells <- as.numeric(table(sce$Sample))

## Determine how to reoder the samples (rows) of the metadata to match the order of sample names in sids vector
m <- match(sids, sce$Sample)

## Create the sample level metadata by combining the reordered metadata with the number of cells corresponding to each sample.
ei <- data.frame(colData(sce)[m, ],
                 n_cells, row.names = NULL) |>
  # select(-"cell_subtype")
  select(-"Cell_type")
ei

# Aggregate the counts per sample_id and cluster_id

# Subset metadata to only include the cluster and sample IDs to aggregate across
groups <- colData(sce)[, c("Cell_type", "Sample")]
# groups <- colData(sce)[, c("cell_subtype", "Sample")]


# Aggregate across cluster-sample groups
# install.packages("grr")
# install.packages("https://cran.r-project.org/src/contrib/Archive/Matrix.utils/Matrix.utils_0.9.8.tar.gz", type = "source", repos = NULL)
# pb <- aggregate.Matrix(t(counts(sce)),
#                        groupings = groups, fun = "sum")
pb <- aggregate.Matrix(t(sce@assays@data@listData$counts),
                       groupings = groups, fun = "sum")

class(pb)
dim(pb)

counts(sce)[1:5, 1:5]
(pb[1:6, 1:6])

# Not every cluster is present in all samples; create a vector that represents how to split samples
splitf <- sapply(stringr::str_split(rownames(pb),
                                    pattern = "_",
                                    n = 2),
                 `[`, 1)

# Turn into a list and split the list into components for each cluster and transform, so rows are genes and columns are samples and make rownames as the sample IDs
# foo <- pb

pb <- split.data.frame(pb,
                       factor(splitf)) %>%
  lapply(function(u)
    set_colnames(t(u),
                 # stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+"))
                 sub(".*_", "", rownames(u))
    ))

class(pb)

# Explore the different components of list
str(pb)

# Print out the table of cells in each cluster-sample group
options(width = 100)
table(sce$Cell_type, sce$Sample)
# table(sce$cell_subtype, sce$Sample)


# Get sample names for each of the cell type clusters

# prep. data.frame for plotting
get_sample_ids <- function(x){
  pb[[x]] |>
    colnames()
}

de_samples <- purrr::map(1:length(kids), get_sample_ids) %>%
  unlist()

# Get cluster IDs for each of the samples

samples_list <- map(1:length(kids), get_sample_ids)

get_cluster_ids <- function(x){
  rep(names(pb)[x],
      each = length(samples_list[[x]]))
}

de_cluster_ids <- map(1:length(kids), get_cluster_ids) %>%
  unlist()

# Create a data frame with the sample IDs, cluster IDs and condition

gg_df <- data.frame(cluster_id = de_cluster_ids,
                    sample_id = de_samples)

gg_df <- left_join(gg_df, ei[, c("Sample", "Class")], by = c("sample_id" = "Sample"))


metadata <- gg_df |>
  dplyr::select(cluster_id, sample_id, Class)

metadata
#::::::::::::::::::::::::::::
# Generate vector of cluster IDs
clusters <- levels(factor(metadata$cluster_id))
clusters

# Subset the metadata to only the B cells
cluster_metadata <- metadata[which(metadata$cluster_id == clusters[1]), ]
head(cluster_metadata)

# Assign the rownames of the metadata to be the sample IDs
rownames(cluster_metadata) <- cluster_metadata$sample_id
rownames(cluster_metadata) <- gsub("-", ".", rownames(cluster_metadata))
cluster_metadata$sample_id <- gsub("-", ".", cluster_metadata$sample_id)
head(cluster_metadata)

#############################################
ei <- data.frame(colData(sce)[m, ],
                 n_cells, row.names = NULL) %>%
  select(-"Cell_type")
# ei
for(cluster in clusters){
  print(cluster)
  cluster_metadata <- metadata[which(metadata$cluster_id == cluster), ]
  head(cluster_metadata)

  # Assign the rownames of the metadata to be the sample IDs
  rownames(cluster_metadata) <- cluster_metadata$sample_id
  rownames(cluster_metadata) <- gsub("-", ".", rownames(cluster_metadata))
  cluster_metadata$sample_id <- gsub("-", ".", cluster_metadata$sample_id)
  head(cluster_metadata)
  counts <- pb[[cluster]]
  colnames(counts) <- gsub("-", ".", colnames(counts))
  cluster_counts <- data.frame(counts[, which(colnames(counts) %in% rownames(cluster_metadata))])

  all(rownames(cluster_metadata) == colnames(cluster_counts))
  dds <- DESeqDataSetFromMatrix(cluster_counts,
                                colData = cluster_metadata,
                                design = ~ Class)
  dds <- DESeq(dds)
  contrast <- c("Class", levels(factor(cluster_metadata$Class))[2], levels(factor(cluster_metadata$Class))[1])

  # resultsNames(dds)
  res <- results(dds,
                 contrast = contrast,
                 alpha = 0.05)

  res <- lfcShrink(dds,
                   type="normal",
                   contrast =  contrast,
                   res=res)

  # Turn the results object into a tibble for use with tidyverse functions
  res_tbl <- res %>%
    data.frame() %>%
    rownames_to_column(var="gene") %>%
    as_tibble()

  write.csv(res_tbl,
            paste0("results/", cluster, "_", levels(factor(cluster_metadata$Class))[2], "_vs_", levels(factor(cluster_metadata$Class))[1], "_all_genes.csv"),
            quote = FALSE,
            row.names = FALSE)

  # Set thresholds
  padj_cutoff <- 0.05

  # Subset the significant results
  sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>%
    dplyr::arrange(padj)

  # Check significant genes output
  sig_res

  # Write significant results to file
  write.csv(sig_res,
            paste0("results/", cluster, "_", levels(cluster_metadata$Class)[2], "_vs_", levels(cluster_metadata$Class)[1], "_sig_genes.csv"),
            quote = FALSE,
            row.names = FALSE)

  ## ggplot of top genes
  normalized_counts <- counts(dds,
                              normalized = TRUE)

  ## Order results by padj values
  top20_sig_genes <- sig_res %>%
    dplyr::arrange(padj) %>%
    dplyr::pull(gene) %>%
    head(n=20)


  top20_sig_norm <- data.frame(normalized_counts) %>%
    rownames_to_column(var = "gene") %>%
    dplyr::filter(gene %in% top20_sig_genes)

  gathered_top20_sig <- top20_sig_norm %>%
    gather(colnames(top20_sig_norm)[2:length(colnames(top20_sig_norm))], key = "samplename", value = "normalized_counts")

  eii <- ei |> rename("sample_id" = "Sample", "group_id" = "Class")
  eii$sample_id =  gsub("-", ".", eii$sample_id)
  gathered_top20_sig <- inner_join(eii[, c("sample_id", "group_id" )], gathered_top20_sig, by = c("sample_id" = "samplename"))

  ## plot using ggplot2
  ggplot(gathered_top20_sig) +
    geom_point(aes(x = gene,
                   y = normalized_counts,
                   color = group_id),
               position=position_jitter(w=0.1,h=0)) +
    scale_y_log10() +
    xlab("Genes") +
    ylab("log10 Normalized Counts") +
    ggtitle("Top 20 Significant DE Genes") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(plot.title = element_text(hjust = 0.5))

  ggsave(paste0("results/",
                cluster,
                "_",
                levels(cluster_metadata$Class)[2],
                "_vs_",
                levels(cluster_metadata$Class)[1],
                "_top20sig_genes.png"),
         last_plot())

  # Extract normalized counts for only the significant genes
  sig_norm <- data.frame(normalized_counts) %>%
    rownames_to_column(var = "gene") %>%
    dplyr::filter(gene %in% sig_res$gene)

  # Set a color palette
  heat_colors <- brewer.pal(6, "YlOrRd")

  # Run pheatmap using the metadata data frame for the annotation
  pheatmap(sig_norm[ , 2:length(colnames(sig_norm))],
           color = heat_colors,
           cluster_rows = T,
           show_rownames = F,
           annotation = cluster_metadata[, c("Class", "cluster_id")],
           border_color = NA,
           fontsize = 10,
           scale = "row",
           fontsize_row = 10,
           height = 20)

  ggsave(paste0("results/",
                cluster,
                "_",
                levels(cluster_metadata$Class)[2],
                "_vs_",
                levels(cluster_metadata$Class)[1],
                "_heatmap.png"),
         last_plot())

  ## Obtain logical vector where TRUE values denote padj values < 0.05 and fold change > 1.5 in either direction
  res_table_thres <- res_tbl %>%
    mutate(threshold = padj < 0.05 & abs(log2FoldChange) >= 0.58)

  ## Volcano plot
  ggplot(res_table_thres) +
    geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold)) +
    ggtitle("Volcano plot of stimulated B cells relative to control") +
    xlab("log2 fold change") +
    ylab("-log10 adjusted p-value") +
    scale_y_continuous(limits = c(0,50)) +
    theme(legend.position = "none",
          plot.title = element_text(size = rel(1.5), hjust = 0.5),
          axis.title = element_text(size = rel(1.25)))

  ggsave(paste0("results/",
                cluster,
                "_",
                levels(cluster_metadata$Class)[2],
                "_vs_",
                levels(cluster_metadata$Class)[1],
                "_volcano.png"),
         last_plot())
}
