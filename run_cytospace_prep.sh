#!/bin/bash

# Directory where spatial data directories are located
input_dir="spatial"
output_dir="spatial_cyto_inputs"

# Check if the input directory exists
if [ ! -d "$input_dir" ]; then
  echo "Input directory does not exist: $input_dir"
  exit 1
fi

# Iterate over each directory in the spatial directory
for dir in "$input_dir"/*; do
  # Extract the base name of the directory
  base_name=$(basename "$dir")
  
  # Construct the input and output paths
  input_path="$input_dir/$base_name"
  output_path="$output_dir/$base_name"
  
  # Check if the directory exists in the output directory, if not, create it
  if [ ! -d "$output_path" ]; then
    mkdir -p "$output_path"
  fi

  # Execute the R script with the constructed paths
  Rscript generate_cytospace_input_from_spaceranger_output.R "$input_path" "$output_path"
done

