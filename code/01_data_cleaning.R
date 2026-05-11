# ==========================================================================================================
# RNA Stability Data Cleaning
# Adam Wait
# Project: RNA Stability Prediction
# Description: This script reads and cleans RNA stability data from supplemental Excel tables (S1 and S2),
#              formats the data, labels coding vs noncoding transcripts, and prepares a combined dataset
#              for downstream feature engineering and modeling.
# ==========================================================================================================

# Load required libraries
library(readxl)
library(stringr)

# ----------------------------------------------------------------------------------------------------------
# Set Directories
# ----------------------------------------------------------------------------------------------------------

input_dir <- "Input/paper1_raw"
output_dir <- "Output"

# ----------------------------------------------------------------------------------------------------------
# Read Excel Sheets
# ----------------------------------------------------------------------------------------------------------

# Table S1: Coding transcripts (RefSeq NM)
sheet1 <- read_excel(
  file.path(input_dir, "Tani_Supp_Tables_revised2.xls"),
  sheet = 1,
  skip = 2
)

# Table S2: Noncoding transcripts (RefSeq NR)
sheet2 <- read_excel(
  file.path(input_dir, "Tani_Supp_Tables_revised2.xls"),
  sheet = 2,
  skip = 2
)

# ----------------------------------------------------------------------------------------------------------
# Rename Columns
# ----------------------------------------------------------------------------------------------------------

colnames(sheet1) <- c("transcript_id", "length", "expression", "half_life")
colnames(sheet2) <- c("transcript_id", "length", "expression", "half_life")

# ----------------------------------------------------------------------------------------------------------
# Clean Transcript IDs
# ----------------------------------------------------------------------------------------------------------

# Remove trailing commas (e.g., "NM_12345," → "NM_12345")
sheet1$transcript_id <- sub(",.*", "", sheet1$transcript_id)
sheet2$transcript_id <- sub(",.*", "", sheet2$transcript_id)
# ----------------------------------------------------------------------------------------------------------
# Convert Columns to Numeric
# ----------------------------------------------------------------------------------------------------------

sheet1$length <- as.numeric(sheet1$length)
sheet1$expression <- as.numeric(sheet1$expression)
sheet1$half_life <- as.numeric(sheet1$half_life)

sheet2$length <- as.numeric(sheet2$length)
sheet2$expression <- as.numeric(sheet2$expression)
sheet2$half_life <- as.numeric(sheet2$half_life)

# ----------------------------------------------------------------------------------------------------------
# Remove Missing / Invalid Rows
# ----------------------------------------------------------------------------------------------------------

sheet1 <- subset(sheet1,
                 !is.na(transcript_id) &
                 !is.na(length) &
                 !is.na(expression) &
                 !is.na(half_life))

sheet2 <- subset(sheet2,
                 !is.na(transcript_id) &
                 !is.na(length) &
                 !is.na(expression) &
                 !is.na(half_life))

# ----------------------------------------------------------------------------------------------------------
# Add Transcript Type Labels
# ----------------------------------------------------------------------------------------------------------

# Label each dataset so we can distinguish coding vs noncoding later
sheet1$type <- "coding"
sheet2$type <- "noncoding"

# ----------------------------------------------------------------------------------------------------------
# Ensure Column Order Matches Before Merge
# ----------------------------------------------------------------------------------------------------------

sheet1 <- sheet1[, c("transcript_id", "length", "expression", "half_life", "type")]
sheet2 <- sheet2[, c("transcript_id", "length", "expression", "half_life", "type")]

# ----------------------------------------------------------------------------------------------------------
# Combine Datasets
# ----------------------------------------------------------------------------------------------------------

# Merge coding and noncoding transcripts into one dataset
combined_data <- rbind(sheet1, sheet2)

# ----------------------------------------------------------------------------------------------------------
# Self Checks
# ----------------------------------------------------------------------------------------------------------

# Check total number of rows
nrow(combined_data)

# Check distribution of transcript types
table(combined_data$type)

# Summary statistics
summary(combined_data)

# Preview combined dataset
head(combined_data)

# ----------------------------------------------------------------------------------------------------------
# Save Cleaned Dataset
# ----------------------------------------------------------------------------------------------------------

write.csv(
  combined_data,
  file.path(output_dir, "stability_combined.csv"),
  row.names = FALSE
)

# ==========================================================================================================
# END OF DATA CLEANING PHASE
# ==========================================================================================================