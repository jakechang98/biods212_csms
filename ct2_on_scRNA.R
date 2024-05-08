# Running CytoTRACE 2 on 23 primary CRC patient data

# Load scRNA-seq data and prepare for input
# CytoTRACE 2 requires a single-cell RNA-sequencing gene expression object as input, where genes are rows and cells are columns 
data <- fread("/Users/savagyan/Desktop/Stanford/Winter2024/STEMREM/Project/scRNAseq copy/GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt", data.table = FALSE)
rownames(data) <- data$Index
data$Index <- NULL

# Load annotation data and prepare for input
# CytoTRACE 2 accepts cell phenotype annotations as an optional input. It should be passed as a loaded dataframe containing 
# cell IDs as rownames and phenotype labels (string formatted without special characters) in the first column.
annotation <- fread("/Users/savagyan/Desktop/Stanford/Winter2024/STEMREM/Project/scRNAseq copy/GSE132465_GEO_processed_CRC_10X_cell_annotation.txt", data.table = FALSE)
annotation <- as.data.frame(annotation)
rownames(annotation)  <- annotation$Index
annotation_tumor <- annotation[annotation$Class == "Tumor",]
annotation_tumor$Index <- annotation_tumor$Cell_type

# filtering the scRNA seq to only tumor patient samples
tumor_patients <-c(annotation[annotation$Class == "Tumor",]$Index)  
data_tumor <- data[, tumor_patients]

# running CytoTRACE 2
ct2_results <- cytotrace2(data_tumor_df, species = "human") 
write.csv(ct2_results, "/Users/savagyan/Desktop/Spring 2024/BMI212/biods212_csms/ct2_results.csv") #saving results

# generating plots
plots <- plotData(ct2_results, annotation = annotation_tumor, expression_data = data_tumor_df) 


# the rds overflows, saved the plots manually from Rstudio consolde
# final_Rds <- list(cytotrace2_results = ct2_results, plots = plots)
# saveRDS(object = final_Rds, file = "cytotrace2_results.rds")
