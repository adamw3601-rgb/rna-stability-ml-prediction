rna-stability-ml-prediction


**Author:** Adam Wait  
**Date:** 04/29/2026  

---

## Languages, Tools, and Packages

**Language**

- R (v4.5.2)

**Packages**

- randomForest  
- ggplot2  
- dplyr  

**Tools**

- Positron (Posit IDE)  
- ViennaRNA RNAfold  

---

## Program Overview

This project evaluates whether integrating sequence, structure, and expression features improves RNA stability prediction. RNA stability is represented by transcript half-life, which varies widely across genes and plays an important role in gene regulation.

The program processes transcript-level data and generates features from multiple sources. Sequence-based features include transcript length and GC content, while structural features are derived from RNA secondary structure predictions using RNAfold. Expression values are also incorporated into the dataset.

The workflow consists of data cleaning, feature generation, structure processing, and model development. Regression models are used to predict continuous RNA half-life values, while classification models group transcripts into stability categories (Short, Moderate, Stable, Long).

Model performance is evaluated using R² for regression and accuracy and confusion matrices for classification.

---

## Inputs

**Files**:

- `stability_combined.csv`  
- `stability_with_gc.csv`  
- `stability_with_structure.csv`  
- FASTA files (located in `Input/fasta/`)  
- RNAfold output file (`rnafold_output.txt`)  

**Description**:

- RNA half-life values (continuous)
- Expression values
- Transcript sequences
- RNAfold structure predictions (minimum free energy)

**Parameters**:

- Expression values are log-transformed using log(1 + x)
- RNA half-life values are log-transformed for regression
- RNA stability classes are defined as:
  - Short-lived: < 2 hours  
  - Moderate: 2–4 hours  
  - Stable: 4–8 hours  
  - Long-lived: > 8 hours  

---

## Outputs

The program generates both data files and figures in the `Output` folder.

**Data Outputs**:

- `regression_results.csv`
- `classification_results.csv`
- `classification_confusion_matrix.csv`
- `classification_comparison_summary.csv`

**Intermediate Data**:

- `stability_with_gc.csv`
- `stability_with_structure.csv`
- `stability_combined.csv`

**Figures** (stored in `Output/Figures/`):

- `predicted_vs_actual_half_life.png`
- `confusion_matrix_heatmap.png`
- `actual_vs_predicted_classes.png`

Each output summarizes model performance or visualizes prediction results.

---

## Scripts

- `01_data_cleaning.R`  
  Cleans and prepares the RNA stability dataset for analysis.

- `02_sequence_features.R`  
  Calculates sequence-based features such as transcript length and GC content.

- `03_structure_features.R`  
  Prepares transcript sequences for RNAfold structure prediction.

- `04_structure_processing.R`  
  Processes RNAfold output and extracts structural features (MFE, MFE per nucleotide).

- `05_modeling.R`  
  Performs regression and classification modeling, evaluates performance, and generates output files and figures.

---
