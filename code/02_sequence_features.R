# ==========================================================================================================
# Sequence Feature Engineering
# Adam Wait
# Project: RNA Stability Prediction
# Description: This script loads transcript sequences from RefSeq FASTA files,
#              matches them to the dataset, and computes sequence-based features
#              such as GC content.
# ==========================================================================================================

# Load libraries
library(Biostrings)
library(dplyr)

# ----------------------------------------------------------------------------------------------------------
# Set directories
# ----------------------------------------------------------------------------------------------------------

input_data_dir = "Output"
fasta_dir <- "Input/fasta"
output_dir <- "Output"

# ----------------------------------------------------------------------------------------------------------
# Load cleaned dataset
# ----------------------------------------------------------------------------------------------------------

combined_data <- read.csv(
  file.path(input_data_dir, "stability_combined.csv"),
  stringsAsFactors = FALSE
)

head(combined_data)
nrow(combined_data)

# ----------------------------------------------------------------------------------------------------------
# Find all FASTA files
# ----------------------------------------------------------------------------------------------------------

fasta_files <- list.files(
  path = fasta_dir,
  pattern = "\\.fna(\\.gz)?$",
  full.names = TRUE
)

fasta_files


# ----------------------------------------------------------------------------------------------------------
# Read and combine all transcript FASTA files
# ----------------------------------------------------------------------------------------------------------

all_sequences <- do.call(c, lapply(fasta_files, readDNAStringSet))

all_sequences
length(all_sequences)

# ----------------------------------------------------------------------------------------------------------
# Clean FASTA header names
# ----------------------------------------------------------------------------------------------------------

# Save original FASTA names
seq_ids <- names(all_sequences)

# Keep only the first token before the first space
seq_ids_clean <- sub(" .*", "", seq_ids)

# Remove version numbers (e.g., NM_000014.6 -> NM_000014)
seq_ids_clean <- sub("\\..*", "", seq_ids_clean)

# Assign cleaned names back to the sequence object
names(all_sequences) <- seq_ids_clean

# Check first few cleaned IDs
head(names(all_sequences))

# ----------------------------------------------------------------------------------------------------------
# Match FASTA sequences to transcript IDs in combined dataset
# ----------------------------------------------------------------------------------------------------------

# Match sequences to transcript IDs in the dataset
matched_sequences <- all_sequences[names(all_sequences) %in% combined_data$transcript_id]

# Check how many sequences matched
length(matched_sequences)

# Preview some matched IDs
head(names(matched_sequences))

# Remove duplicate transcript IDs so each transcript appears only once
combined_data <- combined_data[!duplicated(combined_data$transcript_id), ]

# Rebuild dataset so it only contains rows with matching sequences
combined_data_seq <- combined_data[combined_data$transcript_id %in% names(matched_sequences), ]

# Check dimensions again
length(matched_sequences)
nrow(combined_data_seq)
head(combined_data_seq)

# ----------------------------------------------------------------------------------------------------------
# Compute GC content
# ----------------------------------------------------------------------------------------------------------

# Calculate proportion of G and C nucleotides in each matched sequence
gc_freq <- letterFrequency(
  matched_sequences,
  letters = c("G", "C"),
  as.prob = TRUE
)

# Build a small dataframe with transcript IDs and GC content
gc_df <- data.frame(
  transcript_id = names(matched_sequences),
  gc_content = rowSums(gc_freq),
  stringsAsFactors = FALSE
)

# Check GC content table
head(gc_df)
summary(gc_df$gc_content)

# ----------------------------------------------------------------------------------------------------------
# Merge GC content into RNA stability dataset
# ----------------------------------------------------------------------------------------------------------

combined_data_seq <- merge(
  combined_data_seq,
  gc_df,
  by = "transcript_id"
)

# Check merged result
nrow(combined_data_seq)
head(combined_data_seq)
summary(combined_data_seq$gc_content)

write.csv(
  combined_data_seq,
  file.path(output_dir, "stability_with_gc.csv"),
  row.names = FALSE
)

combined_data_seq$seq_length <- width(matched_sequences[combined_data_seq$transcript_id])

summary(combined_data_seq$length_diff)