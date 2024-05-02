#!/bin/bash

# Path to the base directory for spatial cyto inputs
base_dir="spatial_cyto_inputs"

# Find and sort directories within the base directory
directories=$(find "$base_dir" -mindepth 1 -maxdepth 1 -type d | sort)

# Loop through each directory
for dir in $directories; do
    # Extract the directory name
    dir_name=$(basename "$dir")

    # Construct the STP and CP arguments
    stp_path="$dir/ST_data.csv"
    cp_path="$dir/Coordinates.csv"
    output_dir="cytospace_results/${dir_name}"

    mkdir -p "$output_dir"
    # Run the cytospace command with the dynamically constructed paths
    cytospace \
    -sp scRNAseq/GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt \
    -ctp scRNAseq/GSE132465_GEO_processed_CRC_10X_cell_cell_type.txt \
    -stp "$stp_path" \
    -cp "$cp_path" \
    -o  "$output_dir" 

    echo "Processed $dir_name"
done

echo "All directories processed."

