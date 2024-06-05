# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #
# Code for performing CSIDE's deconvolution and spatially variable genes
#
#
# Author: Jake Chang
# Code adapted from package vignette:
# https://raw.githack.com/dmcable/spacexr/master/vignettes/visium_full_regions.html
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #

# To install C-SIDE (package is actually called spacexr)
# options(timeout = 600000000) ### set this to avoid timeout error
# devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)

library(spacexr)
library(dplyr)
library(SingleCellExperiment)
library(beepr)

cluster_type <- "one_dimensional"
if(cluster_type == "multi_dimensional"){
  cluster_path <- "data/multi_dimensional_cluster_data.csv"
}
if(cluster_type == "one_dimensional"){
  cluster_path <- "data/cluster_data.csv"
}

sce <- readRDS("sce.rds")
non_normal_cells <- colData(sce)$Class != "Normal"
sce <- sce[, non_normal_cells]
mast_cells <- colData(sce)$Cell_type != "Mast cells"
sce <- sce[, mast_cells]

samples <- list.files("data/spatial_cyto_inputs")
priority_samples <- c(1, 2, 3, 4, 5, 13, 14)
samples <- samples[priority_samples]

for(samp in samples){
  datadir <- file.path("data/spatial_cyto_inputs", samp)
  counts <- read.csv(file.path(datadir,"ST_data.csv"))
  rownames(counts) <- counts$GENES
  counts <- counts[, -1]
  coords <- read.csv(file.path(datadir,"Coordinates.csv"))
  rownames(coords) <- gsub("-", ".", coords$Spot.ID)
  coords <- coords[, -1]
  coords <- coords |> select(row, col) |> rename(x = "row", y = "col")
  nUMI <- colSums(counts)
  puck <- SpatialRNA(coords, counts, nUMI)
  barcodes <- colnames(puck@counts)
  # plot_puck_continuous(puck, barcodes, puck@nUMI, ylimit = c(0,round(quantile(puck@nUMI,0.9))),
  #                      title ='plot of nUMI')

  reference <- read.csv(cluster_path)
  reference$cell_subtype <- paste(reference$phenotype, reference$cluster, sep = "_")
  reference <- select(reference, barcode, cell_subtype)

  sce_copy <- sce
  colData(sce_copy) <- DataFrame(left_join(data.frame(colData(sce_copy)), reference, by = c("Index" = "barcode")))
  counts <- sce_copy@assays@data@listData$counts
  meta_data <- data.frame(colData(sce_copy))
  cell_types <- meta_data$cell_subtype
  names(cell_types) <- meta_data$Index
  cell_types <- as.factor(cell_types)
  nUMI <- colSums(counts)
  reference <- Reference(counts, cell_types, nUMI)

  start_time <- Sys.time()
  myRCTD <- create.RCTD(puck, reference, max_cores = parallel::detectCores() - 2)
  myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')
  end_time <- Sys.time()
  beep()
  elapsed_time <- end_time - start_time
  print(elapsed_time)

  saveRDS(myRCTD, paste0("results/", samp, "/", samp, ".rds"))
  saveRDS(myRCTD, file.path("results", samp, paste0(samp, ".rds")))

  dir_path <- file.path("results", cluster_type, samp)
  file_path <- file.path(dir_path, paste0(samp, ".rds"))

  # Create the directory if it does not exist
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }

  # Save the object
  saveRDS(myRCTD, file_path)


  barcodes <- colnames(myRCTD@spatialRNA@counts)
  weights <- myRCTD@results$weights
  norm_weights <- normalize_weights(weights)
  cell_types <- unique(colData(sce_copy)$cell_subtype)

  for(type in cell_types){
    print(type)
    plot_puck_continuous(myRCTD@spatialRNA, barcodes, norm_weights[,type], ylimit = c(0,0.5),
                         title =paste('plot of', type, 'weights'))
    ggplot2::ggsave(paste0("plots/visium/", cluster_type, "/", samp, "/", type, "_deconvolution.png"))
  }

  # Split region into quadrants
  region_1_1 <- barcodes[which((myRCTD@spatialRNA@coords$x < quantile(myRCTD@spatialRNA@coords$x, 1/2)) & (myRCTD@spatialRNA@coords$y < quantile(myRCTD@spatialRNA@coords$y, 1/2)))]
  region_1_2 <- barcodes[which((myRCTD@spatialRNA@coords$x < quantile(myRCTD@spatialRNA@coords$x, 1/2)) & (myRCTD@spatialRNA@coords$y > quantile(myRCTD@spatialRNA@coords$y, 1/2)))]
  region_2_1 <- barcodes[which((myRCTD@spatialRNA@coords$x > quantile(myRCTD@spatialRNA@coords$x, 1/2)) & (myRCTD@spatialRNA@coords$y < quantile(myRCTD@spatialRNA@coords$y, 1/2)))]
  region_2_2 <- barcodes[which((myRCTD@spatialRNA@coords$x > quantile(myRCTD@spatialRNA@coords$x, 1/2)) & (myRCTD@spatialRNA@coords$y > quantile(myRCTD@spatialRNA@coords$y, 1/2)))]

  region_list <- list(region_1_1, region_1_2, region_2_1, region_2_2)

  myRCTD@config$max_cores <- parallel::detectCores() - 2
  start_time <- Sys.time()
  myRCTD <- run.CSIDE.regions(myRCTD, region_list,
                              cell_types = c("Epithelial cells_1", "Epithelial cells_2",
                                             "Stromal cells_1", "Stromal cells_2", "Stromal cells_3",
                                             "Myeloids_1", "Myeloids_2"),
                              cell_type_threshold = 1, doublet_mode = F, weight_threshold = 0.75)
  end_time <- Sys.time()
  beep()
  elapsed_time <- end_time - start_time
  print(elapsed_time)

  sig_gene_list <- myRCTD@de_results$sig_gene_list
  saveRDS(myRCTD, file_path)
}

