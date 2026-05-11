# ==========================================================================================================
# RNA Structure Feature Engineering
# Adam Wait
# Project: RNA Stability Prediction
# Description: This script exports matched transcript sequences, runs RNAfold externally,
#              reads the RNAfold output, extracts minimum free energy (MFE), and merges
#              structure-based features into the dataset.
# ==========================================================================================================

library(Biostrings)

# ----------------------------------------------------------------------------------------------------------
# Set directories
# ----------------------------------------------------------------------------------------------------------

input_data_dir <- "Output"
fasta_dir <- "Input/fasta"
structure_dir <- "Input/structure"
output_dir <- "Output"

# Load dataset with sequence features
combined_data_seq <- read.csv(
  file.path(output_dir, "stability_with_gc.csv"),
  stringsAsFactors = FALSE
)

head(combined_data_seq)
nrow(combined_data_seq)

# Load all transcript FASTA files
fasta_files <- list.files(
  path = fasta_dir,
  pattern = "\\.fna$",
  full.names = TRUE
)

all_sequences <- do.call(c, lapply(fasta_files, readDNAStringSet))

# Clean FASTA names
seq_ids <- names(all_sequences)
seq_ids_clean <- sub(" .*", "", seq_ids)
seq_ids_clean <- sub("\\..*", "", seq_ids_clean)
names(all_sequences) <- seq_ids_clean

# Match sequences to transcript IDs in the dataset
matched_sequences <- all_sequences[names(all_sequences) %in% combined_data_seq$transcript_id]

length(matched_sequences)
head(names(matched_sequences))

# Write matched sequences to a FASTA file for RNAfold
writeXStringSet(
  matched_sequences,
  filepath = file.path(structure_dir, "matched_transcripts.fa")
)