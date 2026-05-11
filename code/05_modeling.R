# ==========================================================================================================
# Project B528
# Adam Wait
# Modeling RNA Stability
# ==========================================================================================================

# ----------------------------------------------------------------------------------------------------------
# Libraries
# ----------------------------------------------------------------------------------------------------------
library(randomForest)
library(ggplot2)
library(dplyr)

# ----------------------------------------------------------------------------------------------------------
# Directories
# ----------------------------------------------------------------------------------------------------------
input_dir <- "Input"
output_dir <- "Output"
fig_dir <- file.path(output_dir, "Figures")

# Create Figures folder if it does not already exist
if (!dir.exists(fig_dir)) {
  dir.create(fig_dir)
}

# ==========================================================================================================
# 1. Load Dataset
# ==========================================================================================================

df <- read.csv(
  file.path(output_dir, "stability_with_structure.csv"),
  stringsAsFactors = FALSE
)

# View basic dataset information
head(df)
nrow(df)
colnames(df)

# ==========================================================================================================
# 2. Prepare Variables
# ==========================================================================================================

# Convert transcript type to factor for modeling
df$type <- as.factor(df$type)

# Log-transform skewed variables
# log1p(x) = log(1 + x), which avoids issues with zero values
df$log_half_life <- log1p(df$half_life)
df$log_expression <- log1p(df$expression)

# Create RNA stability classes for classification
# Short: < 2 hours
# Moderate: 2–4 hours
# Stable: 4–8 hours
# Long: > 8 hours
df$stability_class <- cut(
  df$half_life,
  breaks = c(-Inf, 2, 4, 8, Inf),
  labels = c("Short", "Moderate", "Stable", "Long")
)

df$stability_class <- droplevels(as.factor(df$stability_class))

# Check stability class distribution
print(table(df$stability_class))

# ==========================================================================================================
# 3. Evaluation Functions
# ==========================================================================================================

# R-squared for regression models
r2 <- function(y, yhat) {
  1 - sum((y - yhat)^2) / sum((y - mean(y))^2)
}

# Accuracy for classification models
accuracy <- function(actual, predicted) {
  mean(actual == predicted)
}

# ==========================================================================================================
# 4. Train/Test Split
# ==========================================================================================================

set.seed(528)

idx <- sample(seq_len(nrow(df)), size = 0.8 * nrow(df))

train <- df[idx, ]
test  <- df[-idx, ]

# Remove unused factor levels after splitting
train$stability_class <- droplevels(train$stability_class)
test$stability_class  <- droplevels(test$stability_class)

# Check class distributions in train and test sets
print(table(train$stability_class))
print(table(test$stability_class))

# ==========================================================================================================
# 5. Regression Models
# Goal: Predict continuous log-transformed RNA half-life
# ==========================================================================================================

# Baseline model: simple features only
m_base <- lm(
  log_half_life ~ log_expression + length,
  data = train
)

# Sequence model: baseline + sequence-related features
m_seq <- lm(
  log_half_life ~ log_expression + length + gc_content + type,
  data = train
)

# Structure model: structure-related features only
m_struct <- lm(
  log_half_life ~ mfe_per_nt + seq_length,
  data = train
)

# Combined linear model: all available feature types
m_full <- lm(
  log_half_life ~ log_expression + length + gc_content + type + mfe_per_nt + seq_length,
  data = train
)

# Random forest regression model: nonlinear model using all features
rf_reg <- randomForest(
  log_half_life ~ log_expression + length + gc_content + type + mfe_per_nt + seq_length,
  data = train,
  ntree = 100
)

# ----------------------------------------------------------------------------------------------------------
# Regression Predictions
# ----------------------------------------------------------------------------------------------------------
pred_base   <- predict(m_base,   newdata = test)
pred_seq    <- predict(m_seq,    newdata = test)
pred_struct <- predict(m_struct, newdata = test)
pred_full   <- predict(m_full,   newdata = test)
pred_rf     <- predict(rf_reg,   newdata = test)

# ----------------------------------------------------------------------------------------------------------
# Regression Results
# ----------------------------------------------------------------------------------------------------------
regression_results <- data.frame(
  Model = c("Baseline", "Sequence", "Structure", "Combined", "Random Forest"),
  R2 = c(
    r2(test$log_half_life, pred_base),
    r2(test$log_half_life, pred_seq),
    r2(test$log_half_life, pred_struct),
    r2(test$log_half_life, pred_full),
    r2(test$log_half_life, pred_rf)
  )
)

print(regression_results)

write.csv(
  regression_results,
  file.path(output_dir, "regression_results.csv"),
  row.names = FALSE
)

# ==========================================================================================================
# 6. Regression Figure: Predicted vs Actual Half-Life
# ==========================================================================================================

png(
  filename = file.path(fig_dir, "predicted_vs_actual_half_life.png"),
  width = 800,
  height = 600
)

plot(
  test$log_half_life,
  pred_full,
  xlab = "Actual Log Half-life",
  ylab = "Predicted Log Half-life",
  main = "Predicted vs Actual RNA Half-life",
  pch = 16,
  col = rgb(0, 0, 0, 0.5),
  cex = 0.7
)

# Red = perfect prediction
abline(0, 1, col = "red", lwd = 2)

# Blue = model trend
lines(
  lowess(test$log_half_life, pred_full),
  col = "blue",
  lwd = 2
)

legend(
  "topleft",
  legend = c("Ideal Fit (y = x)", "Model Trend"),
  col = c("red", "blue"),
  lwd = 2,
  bty = "n"
)
dev.off()

# ==========================================================================================================
# 7. Classification Models
# Goal: Predict RNA stability category instead of exact half-life
# ==========================================================================================================

# Baseline random forest classifier
rf_class_base <- randomForest(
  stability_class ~ log_expression + length,
  data = train,
  ntree = 100
)

# Sequence random forest classifier
rf_class_seq <- randomForest(
  stability_class ~ log_expression + length + gc_content + type,
  data = train,
  ntree = 100
)

# Structure random forest classifier
rf_class_struct <- randomForest(
  stability_class ~ mfe_per_nt + seq_length,
  data = train,
  ntree = 100
)

# Combined random forest classifier
rf_class_full <- randomForest(
  stability_class ~ log_expression + length + gc_content + type + mfe_per_nt + seq_length,
  data = train,
  ntree = 100
)

# ----------------------------------------------------------------------------------------------------------
# Classification Predictions
# ----------------------------------------------------------------------------------------------------------
pred_class_base   <- predict(rf_class_base,   newdata = test)
pred_class_seq    <- predict(rf_class_seq,    newdata = test)
pred_class_struct <- predict(rf_class_struct, newdata = test)
pred_class_full   <- predict(rf_class_full,   newdata = test)

# ----------------------------------------------------------------------------------------------------------
# Classification Results
# ----------------------------------------------------------------------------------------------------------
classification_results <- data.frame(
  Model = c("RF - Baseline", "RF - Sequence", "RF - Structure", "RF - Combined"),
  Accuracy = c(
    accuracy(test$stability_class, pred_class_base),
    accuracy(test$stability_class, pred_class_seq),
    accuracy(test$stability_class, pred_class_struct),
    accuracy(test$stability_class, pred_class_full)
  )
)

print(classification_results)

write.csv(
  classification_results,
  file.path(output_dir, "classification_results.csv"),
  row.names = FALSE
)

# ==========================================================================================================
# 8. Confusion Matrix
# ==========================================================================================================

# Confusion matrix compares actual stability classes to predicted stability classes
confusion_matrix <- table(
  Actual = test$stability_class,
  Predicted = pred_class_full
)

print(confusion_matrix)

write.csv(
  confusion_matrix,
  file.path(output_dir, "classification_confusion_matrix.csv")
)

# Create a long-format summary table of actual vs predicted combinations
comparison_table <- data.frame(
  Actual = test$stability_class,
  Predicted = pred_class_full
)

comparison_summary <- comparison_table %>%
  group_by(Actual, Predicted) %>%
  summarise(Count = n(), .groups = "drop") %>%
  arrange(desc(Count))

print(comparison_summary)

write.csv(
  comparison_summary,
  file.path(output_dir, "classification_comparison_summary.csv"),
  row.names = FALSE
)

# ==========================================================================================================
# 9. Classification Figure: Confusion Matrix Heatmap
# ==========================================================================================================

conf_df <- as.data.frame(confusion_matrix)

# Convert factors to character so ggplot does not show numeric factor codes
conf_df$Actual <- as.character(conf_df$Actual)
conf_df$Predicted <- as.character(conf_df$Predicted)

conf_heatmap <- ggplot(conf_df, aes(x = Predicted, y = Actual, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Freq), size = 5) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(
    title = "Confusion Matrix Heatmap",
    x = "Predicted Class",
    y = "Actual Class"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(fig_dir, "confusion_matrix_heatmap.png"),
  plot = conf_heatmap,
  width = 8,
  height = 6,
  dpi = 300
)

# ==========================================================================================================
# 10. Classification Figure: Actual vs Predicted Class Distribution
# ==========================================================================================================

dist_plot_df <- data.frame(
  Class = c(as.character(test$stability_class), as.character(pred_class_full)),
  Group = c(
    rep("Actual", length(test$stability_class)),
    rep("Predicted", length(pred_class_full))
  )
)

dist_plot_df$Class <- factor(
  dist_plot_df$Class,
  levels = c("Short", "Moderate", "Stable", "Long")
)

class_dist_plot <- ggplot(dist_plot_df, aes(x = Class, fill = Group)) +
  geom_bar(position = "dodge") +
  labs(
    title = "Actual vs Predicted RNA Stability Classes",
    x = "RNA Stability Class",
    y = "Number of Transcripts"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(fig_dir, "actual_vs_predicted_classes.png"),
  plot = class_dist_plot,
  width = 8,
  height = 6,
  dpi = 300
)

# ==========================================================================================================
# End of Script
# ==========================================================================================================