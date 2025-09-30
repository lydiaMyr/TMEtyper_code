# Load required libraries
library(caret)
library(gbm)      # For GBM
library(xgboost)  # For XGBoost
library(pls)      # For PLS
library(pamr)     # For PAM

# Set random seed for reproducibility
set.seed(123)

# Data preparation function
prepare_data <- function(train_data) {
  # Ensure the last column is the label and named 'Label'
  if (colnames(train_data)[ncol(train_data)] != "Label") {
    colnames(train_data)[ncol(train_data)] <- "Label"
  }
  
  # Convert label to factor with 7 levels
  train_data$Label <- as.factor(train_data$Label)
  
  # Remove rows with missing values
  train_data <- na.omit(train_data)
  
  return(train_data)
}

# GBM model with feature importance
train_gbm <- function(train_data) {
  cat("Training GBM model...\n")
  
  # GBM tuning parameters
  gbm_grid <- expand.grid(
    n.trees = 150,
    interaction.depth = 3,
    shrinkage = 0.1,
    n.minobsinnode = 10
  )
  
  # Train control settings
  ctrl <- trainControl(
    method = "cv",
    number = 5,
    classProbs = TRUE,
    summaryFunction = multiClassSummary,
    verboseIter = FALSE
  )
  
  # Train GBM model
  gbm_model <- train(
    Label ~ .,
    data = train_data,
    method = "gbm",
    tuneGrid = gbm_grid,
    trControl = ctrl,
    verbose = FALSE
  )
  
  # Extract variable importance
  gbm_importance <- varImp(gbm_model, scale = FALSE)$importance
  gbm_importance$Feature <- rownames(gbm_importance)
  gbm_importance$Model <- "GBM"
  
  return(list(model = gbm_model, importance = gbm_importance))
}

# XGBoost model with feature importance
train_xgboost <- function(train_data) {
  cat("Training XGBoost model...\n")
  
  # XGBoost tuning parameters
  xgb_grid <- expand.grid(
    nrounds = 150,
    max_depth = 6,
    eta = 0.3,
    gamma = 0,
    colsample_bytree = 0.8,
    subsample = 0.8,
    min_child_weight = 1
  )
  
  # Train control settings
  ctrl <- trainControl(
    method = "cv",
    number = 5,
    classProbs = TRUE,
    summaryFunction = multiClassSummary,
    verboseIter = FALSE
  )
  
  # Train XGBoost model
  xgb_model <- train(
    Label ~ .,
    data = train_data,
    method = "xgbTree",
    tuneGrid = xgb_grid,
    trControl = ctrl,
    verbose = FALSE
  )
  
  # Extract variable importance
  xgb_importance <- varImp(xgb_model, scale = FALSE)$importance
  xgb_importance$Feature <- rownames(xgb_importance)
  xgb_importance$Model <- "XGBoost"
  
  return(list(model = xgb_model, importance = xgb_importance))
}

# Random Forest model with feature importance
train_rf <- function(train_data) {
  cat("Training Random Forest model...\n")
  
  # Calculate mtry (sqrt of number of features)
  n_features <- ncol(train_data) - 1  # Exclude label column
  mtry_value <- round(sqrt(n_features))
  
  # Random Forest tuning parameters
  rf_grid <- expand.grid(
    mtry = mtry_value
  )
  
  # Train control settings
  ctrl <- trainControl(
    method = "cv",
    number = 5,
    classProbs = TRUE,
    summaryFunction = multiClassSummary,
    verboseIter = FALSE
  )
  
  # Train Random Forest model
  rf_model <- train(
    Label ~ .,
    data = train_data,
    method = "rf",
    tuneGrid = rf_grid,
    trControl = ctrl,
    ntree = 500,
    importance = TRUE  # Ensure importance calculation
  )
  
  # Extract variable importance
  rf_importance <- varImp(rf_model, scale = FALSE)$importance
  rf_importance$Feature <- rownames(rf_importance)
  rf_importance$Model <- "RandomForest"
  
  return(list(model = rf_model, importance = rf_importance))
}

# Logit Boost model with feature importance
train_logitboost <- function(train_data) {
  cat("Training Logit Boost model...\n")
  
  # Logit Boost tuning parameters
  lb_grid <- expand.grid(
    nIter = 150
  )
  
  # Train control settings
  ctrl <- trainControl(
    method = "cv",
    number = 5,
    classProbs = TRUE,
    summaryFunction = multiClassSummary,
    verboseIter = FALSE
  )
  
  # Train Logit Boost model
  lb_model <- train(
    Label ~ .,
    data = train_data,
    method = "LogitBoost",
    tuneGrid = lb_grid,
    trControl = ctrl
  )
  
  # Extract variable importance
  lb_importance <- varImp(lb_model, scale = FALSE)$importance
  lb_importance$Feature <- rownames(lb_importance)
  lb_importance$Model <- "LogitBoost"
  
  return(list(model = lb_model, importance = lb_importance))
}

# PLS model with feature importance
train_pls <- function(train_data) {
  cat("Training PLS model...\n")
  
  # PLS tuning parameters - ncomp will be optimized
  pls_grid <- expand.grid(
    ncomp = 1:min(20, ncol(train_data) - 1)
  )
  
  # Train control settings
  ctrl <- trainControl(
    method = "cv",
    number = 5,
    classProbs = TRUE,
    summaryFunction = multiClassSummary,
    verboseIter = FALSE
  )
  
  # Train PLS model
  pls_model <- train(
    Label ~ .,
    data = train_data,
    method = "pls",
    tuneGrid = pls_grid,
    trControl = ctrl,
    probMethod = "softmax"
  )
  
  # Extract variable importance
  pls_importance <- varImp(pls_model, scale = FALSE)$importance
  pls_importance$Feature <- rownames(pls_importance)
  pls_importance$Model <- "PLS"
  
  return(list(model = pls_model, importance = pls_importance))
}

# PAM model with feature importance
train_pam <- function(train_data) {
  cat("Training PAM model...\n")
  
  # Create threshold sequence from 0 to 2
  threshold_seq <- seq(0, 2, by = 0.1)
  
  # PAM tuning parameters
  pam_grid <- expand.grid(
    threshold = threshold_seq
  )
  
  # Train control settings
  ctrl <- trainControl(
    method = "cv",
    number = 10,  # 10-fold CV as specified
    classProbs = TRUE,
    summaryFunction = multiClassSummary,
    verboseIter = FALSE
  )
  
  # Train PAM model
  pam_model <- train(
    Label ~ .,
    data = train_data,
    method = "pam",
    tuneGrid = pam_grid,
    trControl = ctrl
  )
  
  # Extract variable importance
  pam_importance <- varImp(pam_model, scale = FALSE)$importance
  pam_importance$Feature <- rownames(pam_importance)
  pam_importance$Model <- "PAM"
  
  return(list(model = pam_model, importance = pam_importance))
}

# Main function to run all models and extract feature importance
run_feature_selection <- function(train_data, output_file = "feature_importance_results.csv") {
  
  # Prepare data
  cat("=== Preparing Data ===\n")
  train_data <- prepare_data(train_data)
  cat("Data dimensions:", nrow(train_data), "samples,", ncol(train_data), "features (including label)\n")
  cat("Label distribution:\n")
  print(table(train_data$Label))
  
  # List to store all results
  all_importance <- list()
  all_models <- list()
  
  # Train all models
  tryCatch({
    # GBM
    gbm_result <- train_gbm(train_data)
    all_importance$GBM <- gbm_result$importance
    all_models$GBM <- gbm_result$model
    
    # XGBoost
    xgb_result <- train_xgboost(train_data)
    all_importance$XGBoost <- xgb_result$importance
    all_models$XGBoost <- xgb_result$model
    
    # Random Forest
    rf_result <- train_rf(train_data)
    all_importance$RandomForest <- rf_result$importance
    all_models$RandomForest <- rf_result$model
    
    # Logit Boost
    lb_result <- train_logitboost(train_data)
    all_importance$LogitBoost <- lb_result$importance
    all_models$LogitBoost <- lb_result$model
    
    # PLS
    pls_result <- train_pls(train_data)
    all_importance$PLS <- pls_result$importance
    all_models$PLS <- pls_result$model
    
    # PAM
    pam_result <- train_pam(train_data)
    all_importance$PAM <- pam_result$importance
    all_models$PAM <- pam_result$model
    
  }, error = function(e) {
    cat("Error in model training:", e$message, "\n")
  })
  
  # Combine all importance results
  combined_importance <- do.call(rbind, all_importance)
  
  # Reshape to wide format for easier comparison
  importance_wide <- combined_importance %>%
    group_by(Feature) %>%
    summarise(
      GBM = ifelse("GBM" %in% Model, Overall[Model == "GBM"], NA),
      XGBoost = ifelse("XGBoost" %in% Model, Overall[Model == "XGBoost"], NA),
      RandomForest = ifelse("RandomForest" %in% Model, Overall[Model == "RandomForest"], NA),
      LogitBoost = ifelse("LogitBoost" %in% Model, Overall[Model == "LogitBoost"], NA),
      PLS = ifelse("PLS" %in% Model, Overall[Model == "PLS"], NA),
      PAM = ifelse("PAM" %in% Model, Overall[Model == "PAM"], NA),
      .groups = 'drop'
    )
  
  # Calculate average importance across all models
  importance_wide$Average_Importance <- rowMeans(importance_wide[, -1], na.rm = TRUE)
  
  # Sort by average importance
  importance_wide <- importance_wide[order(-importance_wide$Average_Importance), ]
  
  # Save results to CSV
  write.csv(importance_wide, output_file, row.names = FALSE)
  cat("Feature importance results saved to:", output_file, "\n")
  
  # Print top 20 most important features
  cat("\n=== Top 20 Most Important Features ===\n")
  print(head(importance_wide, 20))
  
  # Return results
  return(list(
    importance_table = importance_wide,
    models = all_models,
    combined_importance = combined_importance
  ))
}

# Function to select top features based on average importance
select_top_features <- function(importance_results, top_n = 50) {
  importance_table <- importance_results$importance_table
  
  # Select top N features
  top_features <- head(importance_table$Feature, top_n)
  
  cat("Selected top", top_n, "features:\n")
  print(top_features)
  
  return(top_features)
}

# Function to create feature importance plot
plot_feature_importance <- function(importance_results, top_n = 30) {
  library(ggplot2)
  library(tidyr)
  
  importance_table <- importance_results$importance_table
  
  # Select top N features for plotting
  top_features <- head(importance_table, top_n)
  
  # Convert to long format for plotting
  plot_data <- top_features %>%
    select(Feature, GBM, XGBoost, RandomForest, LogitBoost, PLS, PAM) %>%
    pivot_longer(cols = -Feature, names_to = "Model", values_to = "Importance")
  
  # Create plot
  p <- ggplot(plot_data, aes(x = reorder(Feature, Importance), y = Importance, fill = Model)) +
    geom_bar(stat = "identity", position = "dodge") +
    coord_flip() +
    labs(title = paste("Top", top_n, "Feature Importance Across Models"),
         x = "Features", y = "Importance Score") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  print(p)
  
  return(p)
}

# Usage example:
# results <- run_feature_selection(train_data, "my_feature_importance.csv")
# top_50_features <- select_top_features(results, 50)
# feature_plot <- plot_feature_importance(results, 30)

# Save the top features for downstream analysis
# write.csv(top_50_features, "top_50_features.csv", row.names = FALSE)