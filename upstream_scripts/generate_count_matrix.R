# ------------------------------------------------------------------------------
# Script: generate_count_matrix.R
# Description: Merge Kallisto est_counts from multiple samples into one CSV matrix(Count Matrix)
# Author: Shiti Ghosh
# ------------------------------------------------------------------------------

# Load required libraries
library(readr)
library(dplyr)

# Set the main Kallisto output directory
kallisto_dir <- "/data/output/kallisto_output"

# Get a list of all sample subdirectories
sample_dirs <- list.dirs(path = kallisto_dir, full.names = TRUE, recursive = FALSE)

# Construct full paths to each .tsv file (renamed from abundance.tsv earlier)
tsv_files <- sapply(sample_dirs, function(dir) {
  file.path(dir, paste0(basename(dir), ".tsv"))
})

# Function to read a single Kallisto .tsv file and extract est_counts
read_kallisto_counts <- function(file) {
  data <- read_tsv(file, show_col_types = FALSE)
  counts <- data$est_counts
  names(counts) <- data$target_id
  return(counts)
}

# Read and merge data from all samples
count_list <- lapply(tsv_files, read_kallisto_counts)

# Use transcript IDs from the first file as rownames
transcript_ids <- names(count_list[[1]])

# Combine into matrix
count_matrix <- do.call(cbind, count_list)
rownames(count_matrix) <- transcript_ids
colnames(count_matrix) <- basename(sample_dirs)

# Save the final count matrix as CSV
write.csv(count_matrix, file = file.path(kallisto_dir, "count_matrix.csv"), row.names = TRUE)

cat("count_matrix.csv has been successfully created.\n")
