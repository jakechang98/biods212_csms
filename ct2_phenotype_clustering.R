library(dplyr)
library(ggplot2)
library(scater)
library(tibble)
library(umap)

multi_dimensional_k_means <- function() {
    single_cell_data <- readRDS("data/sce.rds")
    
    tumor_cells <- colData(single_cell_data)$Class != "Normal"
    single_cell_data <- single_cell_data[, tumor_cells]
    
    single_cell_data <- runPCA(single_cell_data, ncomponents = 10, exprs_values = "logcounts")
    
    cell_types <- data.frame(single_cell_data$Index, single_cell_data$Cell_type)
    
    developmental_potential_scores <- read.csv("data/cytotrace_results.csv")
    
    reduced_dimensions <- data.frame(reducedDims(single_cell_data)$PCA)
    reduced_dimensions <- rownames_to_column(reduced_dimensions, "X")
    
    reduced_dimensions <- merge(x = developmental_potential_scores[, c("X", "CytoTRACE2_Score")], y = reduced_dimensions, by = "X")
    reduced_dimensions <- merge(x = reduced_dimensions, y = cell_types, by.x = "X", by.y = "single_cell_data.Index")
    
    grouped_collation <- reduced_dimensions %>% group_by(reduced_dimensions$single_cell_data.Cell_type)
    grouped_collation <- grouped_collation %>% group_split()
    
    set.seed(212)
    
    column_names <- c("barcode", "phenotype", "value", "cluster")
    cluster_data <- data.frame(matrix(ncol = length(column_names), nrow = 0))
    colnames(cluster_data) <- column_names
    
    k_values <- 1:30
    
    for (grouped_data in grouped_collation) {
        if (unique(grouped_data$single_cell_data.Cell_type) != "Mast cells") {
            wss <- sapply(k_values, function(k) {
                kmeans(grouped_data[, -which(names(grouped_data) %in% c("X", "single_cell_data.Cell_type", 
                                                                        "reduced_dimensions$single_cell_data.Cell_type"))], centers = k)$tot.withinss
            })
            
            elbow_data <- data.frame(k = k_values, wss = wss)
            
            ggplot(elbow_data, aes(x = k, y = wss)) +
                geom_line() +
                geom_point() +
                labs(title = sprintf("elbow_plot_%s.pdf", unique(grouped_data$single_cell_data.Cell_type)),
                     x = "number of clusters (k)",
                     y = "within-cluster sum of squares (wss)") +
                theme_minimal()
            
            ggsave(sprintf("figures/k_means_clustering_elbow_%s.pdf", unique(grouped_data$single_cell_data.Cell_type)))
        }
        
        if (unique(grouped_data$single_cell_data.Cell_type) == "B cells") {
            k_means_result <- kmeans(grouped_data$CytoTRACE2_Score, centers = 2)
            
        } else if (unique(grouped_data$single_cell_data.Cell_type) == "Epithelial cells" | 
                   unique(grouped_data$single_cell_data.Cell_type) == "Stromal cells" |
                   unique(grouped_data$single_cell_data.Cell_type) == "T cells") {
            
            k_means_result <- kmeans(grouped_data$CytoTRACE2_Score, centers = 5)
            
        } else if (unique(grouped_data$single_cell_data.Cell_type) == "Myeloids") {
            k_means_result <- kmeans(grouped_data$CytoTRACE2_Score, centers = 3)
            
        } else if (unique(grouped_data$single_cell_data.Cell_type) == "Mast cells") {
            k_means_result <- kmeans(grouped_data$CytoTRACE2_Score, centers = 1)
        }
        
        phenotype_clusters <- data.frame(barcode = grouped_data$X, phenotype = grouped_data$single_cell_data.Cell_type,
                                         value = grouped_data$CytoTRACE2_Score, cluster = factor(k_means_result$cluster))
        
        cluster_data <- rbind(cluster_data, phenotype_clusters)
        
        if (unique(grouped_data$single_cell_data.Cell_type) != "Mast cells") {
            umap_result <- umap(grouped_data %>% select("CytoTRACE2_Score", starts_with("PC")))
            umap_result <- as.data.frame(umap_result$layout)
            
            colnames(umap_result) <- c("UMAP1", "UMAP2")
            umap_result$cluster <- phenotype_clusters$cluster
            umap_result$score <- phenotype_clusters$value
            
            ggplot(umap_result, aes(x = UMAP1, y = UMAP2, color = cluster)) +
                geom_point(size = 2) +
                theme_minimal() +
                labs(title = sprintf("%s", unique(grouped_data$single_cell_data.Cell_type)),
                     x = "UMAP1",
                     y = "UMAP2",
                     color = "cluster") +
                theme(plot.title = element_text(hjust = 0.5))
            
            ggsave(sprintf("figures/umap_cluster_%s.pdf", unique(grouped_data$single_cell_data.Cell_type)))
            
            ggplot(umap_result, aes(x = UMAP1, y = UMAP2, color = score)) +
                geom_point(size = 2) +
                scale_color_gradient(low = "blue", high = "red") +
                theme_minimal() +
                labs(title = sprintf("%s", unique(grouped_data$single_cell_data.Cell_type)),
                     x = "UMAP1",
                     y = "UMAP2",
                     color = "score") +
                theme(plot.title = element_text(hjust = 0.5))
            
            ggsave(sprintf("figures/umap_score_%s.pdf", unique(grouped_data$single_cell_data.Cell_type)))
        }
    }
    
    write.csv(cluster_data, file = "cluster_data.csv", row.names = FALSE)
}

one_dimensional_k_means <- function() {
    single_cell_data <- readRDS("data/sce.rds")
    developmental_potential_scores <- read.csv("data/cytotrace_results.csv")
    
    cell_types <- data.frame(single_cell_data$Index, single_cell_data$Cell_type)

    collated_info <- merge(x = developmental_potential_scores, y = cell_types,
                             by.x = "X", by.y = "single_cell_data.Index")

    grouped_collation <- collated_info %>% group_by(collated_info$single_cell_data.Cell_type)
    grouped_collation <- grouped_collation %>% group_split()

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

        ggsave(sprintf("figures/k_means_clustering_%s.pdf", unique(grouped_data$single_cell_data.Cell_type)))
    }

    write.csv(cluster_data, file = "cluster_data.csv", row.names = FALSE)
}

multi_dimensional_k_means()