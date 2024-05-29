library(dplyr)
library(ggplot2)

single_cell_data <- readRDS("data/sce.rds")
developmental_potential_scores <- read.csv("data/cytotrace_results.csv")

cell_types <- data.frame(single_cell_data$Index, single_cell_data$Cell_type)

collated_info <- merge(x = developmental_potential_scores, y = cell_types,
                       by.x = "X", by.y = "single_cell_data.Index")

grouped_collation <- collated_info %>% group_by(collated_info$single_cell_data.Cell_type)
grouped_collation <- grouped_collation %>% group_split()

# TODO's to try (notes for myself)

# fit a gaussian mixture model (but how different will this be from the 1D case), use BIC / VIC metrics to determine validity
# try adding it as another dimension to the gene counts then cluster everything
# or just look at BIC for different k-means models to find minimum number of clusters that may be needed to explain all the data
# maximize inter-cluster variability, minimize intra-cluster variability

# return spreadsheet with barcoded cell, original cell type, and cluster label

set.seed(212)

column_names <- c("barcode", "phenotype", "value", "cluster")
cluster_data <- data.frame(matrix(ncol = length(column_names), nrow = 0))
colnames(cluster_data) <- column_names

for (grouped_data in grouped_collation) {
    if (unique(grouped_data$single_cell_data.Cell_type) == "B cells" | 
        unique(grouped_data$single_cell_data.Cell_type) == "Stromal cells") {
        k_means_result <- kmeans(grouped_data$CytoTRACE2_Score, centers = 3)
        
    } else if (unique(grouped_data$single_cell_data.Cell_type) == "Epithelial cells" |
               unique(grouped_data$single_cell_data.Cell_type) == "Myeloids" | 
               unique(grouped_data$single_cell_data.Cell_type) == "T cells") {
        k_means_result <- kmeans(grouped_data$CytoTRACE2_Score, centers = 2)
        
    } else if (unique(grouped_data$single_cell_data.Cell_type) == "Mast cells") {
        k_means_result <- kmeans(grouped_data$CytoTRACE2_Score, centers = 1)
    }
    
    x <- data.frame(barcode = grouped_data$X, phenotype = grouped_data$single_cell_data.Cell_type,
                    value = grouped_data$CytoTRACE2_Score, cluster = factor(k_means_result$cluster))
    
    cluster_data <- rbind(cluster_data, x)

    ggplot(x, aes(x = value, y = 0, color = cluster)) +
        geom_point(size = 3)  +
        annotate("segment", x = 0, xend = 1, y = 0, yend = 0, size = 2) +
        annotate("segment", x = 0, xend = 0, y = -0.1, yend = 0.1, size = 2) +
        annotate("segment", x = 1, xend = 1, y = -0.1, yend = 0.1, size = 2) +
        scale_x_continuous(limits = c(0,1)) +
        scale_y_continuous(limits = c(-1,1)) +
        theme(panel.background = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.title = element_blank())

    ggsave(sprintf("k_means_clustering_%s.pdf", unique(grouped_data$single_cell_data.Cell_type)))
}

write.csv(cluster_data, file = "cluster_data.csv", row.names = FALSE)