# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #
# Code for making Visium gene expression plots for final pres and report
#
#
# Author: Jake Chang
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #

# Visualize the expression patterns for some of the interesting
# spatially differentially expressed genes
library(Seurat)
library(SeuratData)

relevant_samples <- c("SN048_A121573_Rep1", "SN048_A121573_Rep2", "SN048_A416371_Rep1",
                      "SN048_A416371_Rep2", "SN123_A551763_Rep1", "SN84_A120838_Rep1",
                      "SN84_A120838_Rep2")
samp <- relevant_samples[1]


visium_data <- Load10X_Spatial(data.dir = file.path(path, samp),
                               assay = "Spatial")
visium_data <- SCTransform(visium_data, assay = "Spatial", verbose = FALSE)
SpatialDimPlot(visium_data)

SpatialFeaturePlot(visium_data,"nCount_Spatial")
SpatialFeaturePlot(visium_data,"nFeature_Spatial")

# Myogenesis genes
SpatialFeaturePlot(visium_data, features = c("COL1A2", "COL6A2", "MMP1", "MMP2", "TAGLN"), image.alpha = 0.25)


# KRAS genes
SpatialFeaturePlot(visium_data, features = c("DDIT4", "GDF15", "MUC4", "MUC5AC"), image.alpha = 0.25)

# EMT genes
SpatialFeaturePlot(visium_data, features = c("MMP1", "CD55", "CEACAM5", "CEACAM6",
                                             "DDIT4", "ERO1A", "KRT19", "LGALS3",
                                             "MET", "S100A6", "SPARC", "TAGLN"), image.alpha = 0.25)

# Plot for slide deck
SpatialFeaturePlot(visium_data, features = c("CEACAM5", "CD55", "SPARC",
                                             "COL1A2", "MMP1", "TAGLN"), image.alpha = 0.25)

SpatialFeaturePlot(visium_data, features = c("MMP1", "CD55", "SPARC"), image.alpha = 0.25)
SpatialFeaturePlot(visium_data, features = c("CEACAM5", "CD55", "SPARC"), image.alpha = 0.25)



SpatialFeaturePlot(visium_data, features = "CD55", image.alpha = 0.25)

SpatialFeaturePlot(visium_data, features = "MUC13", image.alpha = 0.25)
SpatialFeaturePlot(visium_data, features = "REG4", image.alpha = 0.25)
SpatialFeaturePlot(visium_data, features = c("AEBP1", "IGFBP5", "MMP1"), image.alpha = 0.25)
SpatialFeaturePlot(visium_data, features = c("COL1A2", "COL6A2", "DCN", "MMP2", "SPARC", "TAGLN"), image.alpha = 0.25)

SpatialFeaturePlot(visium_data, features = "AEBP1", image.alpha = 0.25)
SpatialFeaturePlot(visium_data, features = "AEBP1", image.alpha = 0.25)
SpatialFeaturePlot(visium_data, features = "AEBP1", image.alpha = 0.25)
SpatialFeaturePlot(visium_data, features = "AEBP1", image.alpha = 0.25)


default_assay <- DefaultAssay(visium_data)

# Extract the gene names from the default assay
gene_list <- rownames(GetAssayData(visium_data, assay = default_assay, slot = "counts"))
