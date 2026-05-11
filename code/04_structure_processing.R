# ==========================================================================================================
# Project B528
# Adam Wait
# RNA Structure Feature Extraction (RNAfold Output Processing)
# ==========================================================================================================

# Set input and output directories
input_dir <- "Input"
output_dir <- "Output"

# ----------------------------------------------------------------------------------------------------------
# Read RNAfold output file
# This file contains repeating blocks of:
#   1) transcript header
#   2) sequence
#   3) structure + MFE
# ----------------------------------------------------------------------------------------------------------
rnafold_lines <- readLines("Input/structure/rnafold_output.txt")

# Total number of lines (sanity check)
length(rnafold_lines)

# ----------------------------------------------------------------------------------------------------------
# Extract header lines and structure/MFE lines
# Headers occur every 3rd line starting at line 1
# Structure + MFE occur every 3rd line starting at line 3
# ----------------------------------------------------------------------------------------------------------
header_lines <- rnafold_lines[seq(1, length(rnafold_lines), by = 3)]
energy_lines <- rnafold_lines[seq(3, length(rnafold_lines), by = 3)]

# Check extraction worked correctly
length(header_lines)
length(energy_lines)
head(header_lines)
head(energy_lines)

# ----------------------------------------------------------------------------------------------------------
# Extract transcript IDs from header lines
# Remove ">" and any extra description text
# ----------------------------------------------------------------------------------------------------------
transcript_id <- sub("^>", "", header_lines)
transcript_id <- sub(" .*", "", transcript_id)

head(transcript_id)

# ----------------------------------------------------------------------------------------------------------
# Extract Minimum Free Energy (MFE) values
# Values are located at the end of structure lines in parentheses
# Example: "...((((....)))) (-461.60)"
# ----------------------------------------------------------------------------------------------------------
mfe <- sub(".*\\(([[:space:]]*[-0-9\\.]+)\\)$", "\\1", energy_lines)
mfe <- as.numeric(trimws(mfe))

# Inspect extracted MFE values
head(mfe)
summary(mfe)

# ----------------------------------------------------------------------------------------------------------
# Create dataframe containing transcript IDs and MFE values
# ----------------------------------------------------------------------------------------------------------
mfe_df <- data.frame(
  transcript_id = transcript_id,
  mfe = mfe,
  stringsAsFactors = FALSE
)

# Check dataframe
head(mfe_df)
nrow(mfe_df)

# ----------------------------------------------------------------------------------------------------------
# Load sequence-based dataset (from previous step)
# This contains GC content, expression, half-life, etc.
# ----------------------------------------------------------------------------------------------------------
combined_data_seq <- read.csv(
  "Output/stability_with_gc.csv",
  stringsAsFactors = FALSE
)

# Inspect dataset
head(combined_data_seq)
nrow(combined_data_seq)
colnames(combined_data_seq)

# ----------------------------------------------------------------------------------------------------------
# Extract sequence lines from RNAfold output (every 3rd line starting at line 2)
# These are used to compute transcript length
# ----------------------------------------------------------------------------------------------------------
sequence_lines <- rnafold_lines[seq(2, length(rnafold_lines), by = 3)]

# Compute sequence length for each transcript
seq_length <- nchar(sequence_lines)

# Inspect sequence lengths
head(seq_length)
summary(seq_length)
length(seq_length)

# ----------------------------------------------------------------------------------------------------------
# Add sequence length and normalized MFE to dataframe
# MFE per nucleotide allows comparison across different transcript lengths
# ----------------------------------------------------------------------------------------------------------
mfe_df$seq_length <- seq_length
mfe_df$mfe_per_nt <- mfe_df$mfe / mfe_df$seq_length

# Inspect updated dataframe
head(mfe_df)
summary(mfe_df$mfe_per_nt)

# ----------------------------------------------------------------------------------------------------------
# Merge sequence features with structure features
# Only transcripts present in both datasets will be retained
# ----------------------------------------------------------------------------------------------------------
combined_data_struct <- merge(
  combined_data_seq,
  mfe_df,
  by = "transcript_id"
)

# Inspect merged dataset
head(combined_data_struct)
nrow(combined_data_struct)
summary(combined_data_struct$mfe)
summary(combined_data_struct$mfe_per_nt)

# ----------------------------------------------------------------------------------------------------------
# Save final dataset for modeling
# ----------------------------------------------------------------------------------------------------------
write.csv(
  combined_data_struct,
  file.path(output_dir, "stability_with_structure.csv"),
  row.names = FALSE
)