################################### ROC curve ####################################
data <- read.xlsx('DATASET_PRE_TREATMENT_334.xlsx')

#convert to binary
data$Response_bin <- ifelse(data$Response == "R", 1,
                            ifelse(data$Response == "NR", 0, NA))
table(data$Response, useNA = "ifany")
table(data$Response_bin, useNA = "ifany")
unique(data$Response_bin)

library(dplyr)
library(pROC)

data <- data %>%
  mutate(
    LA_TAM = ifelse(LA_TAM == 0, 1, LA_TAM)
  )

roc_data <- data %>%
  select(Response_bin, LA_TAM) %>%
  filter(!is.na(Response_bin), !is.na(LA_TAM))

table(roc_data$Response_bin)

roc_obj <- roc(response = roc_data$Response_bin,
               predictor = roc_data$LA_TAM)

auc(roc_obj)
plot(roc_obj, col = "blue", lwd = 2,
     main = paste0("ROC Curve for LA_TAM (AUC = ", round(auc(roc_obj), 3), ")"))

#################################### XgBoost ##################################
library(dplyr)
library(tidyr)
library(xgboost)
library(pROC)

features <- c("LA_TAM")

data2 <- data %>%
  mutate(
    Response_bin = factor(Response, levels = c("NR", "R"))
  ) %>%
  drop_na(all_of(features), Response_bin)

cohorts <- unique(data2$Cohort)

results <- list()
all_preds <- data.frame()
best_params_all <- list()

for (cohort in cohorts) {
  
  cat("Testing on cohort:", cohort, "\n")
  
  train_data <- data2 %>% filter(Cohort != cohort)
  test_data  <- data2 %>% filter(Cohort == cohort)
  
  x_train <- as.matrix(train_data[, features])
  y_train <- train_data$Response_bin
  
  x_test <- as.matrix(test_data[, features])
  y_test <- test_data$Response_bin
  
  # numeric labels for cv
  y_train_num <- ifelse(y_train == "R", 1, 0)
  
  dtrain <- xgb.DMatrix(data = x_train, label = y_train_num)
  
  ##############################################################################
  # Hyperparameter tuning on training data only
  ##############################################################################
  
  grid <- expand.grid(
    nrounds = c(50, 100, 200, 300),
    max_depth = c(1, 2, 3),
    eta = c(0.01, 0.05, 0.1),
    subsample = c(0.7, 0.8, 1.0),
    colsample_bytree = c(0.7, 0.8, 1.0),
    min_child_weight = c(1, 3, 5),
    gamma = c(0, 0.5, 1)
  )
  
  tuning_results <- data.frame()
  
  for (i in seq_len(nrow(grid))) {
    
    params <- list(
      objective = "binary:logistic",
      eval_metric = "auc",
      max_depth = grid$max_depth[i],
      eta = grid$eta[i],
      subsample = grid$subsample[i],
      colsample_bytree = grid$colsample_bytree[i],
      min_child_weight = grid$min_child_weight[i],
      gamma = grid$gamma[i]
    )
    
    cv_fit <- xgb.cv(
      params = params,
      data = dtrain,
      nrounds = grid$nrounds[i],
      nfold = 5,
      stratified = TRUE,
      early_stopping_rounds = 20,
      maximize = TRUE,
      verbose = 0
    )
    
    best_iter <- cv_fit$best_iteration
    best_auc <- cv_fit$evaluation_log$test_auc_mean[best_iter]
    
    tuning_results <- rbind(
      tuning_results,
      data.frame(
        nrounds = grid$nrounds[i],
        max_depth = grid$max_depth[i],
        eta = grid$eta[i],
        subsample = grid$subsample[i],
        colsample_bytree = grid$colsample_bytree[i],
        min_child_weight = grid$min_child_weight[i],
        gamma = grid$gamma[i],
        best_iteration = best_iter,
        cv_auc = best_auc
      )
    )
  }
  
  best_row <- tuning_results %>%
    arrange(desc(cv_auc)) %>%
    slice(1)
  
  print(best_row)
  
  best_params_all[[cohort]] <- best_row
  
  ##############################################################################
  # Final model using best params on full training set
  ##############################################################################
  
  model <- xgb.train(
    params = list(
      objective = "binary:logistic",
      eval_metric = "auc",
      max_depth = best_row$max_depth,
      eta = best_row$eta,
      subsample = best_row$subsample,
      colsample_bytree = best_row$colsample_bytree,
      min_child_weight = best_row$min_child_weight,
      gamma = best_row$gamma
    ),
    data = dtrain,
    nrounds = best_row$best_iteration,
    verbose = 0
  )
  
  preds <- predict(model, newdata = x_test)
  
  y_test_num <- ifelse(y_test == "R", 1, 0)
  
  roc_obj <- roc(response = y_test_num, predictor = preds, quiet = TRUE)
  auc_val <- as.numeric(auc(roc_obj))
  
  cat("AUC:", round(auc_val, 3), "\n")
  
  results[[cohort]] <- auc_val
  
  all_preds <- rbind(
    all_preds,
    data.frame(
      Cohort = cohort,
      Response = y_test_num,
      Prediction = preds
    )
  )
}
############################### Extract importance #################################
  imp <- xgb.importance(feature_names = features, model = model)
  imp$Cohort <- cohort
  
  importance_list[[cohort]] <- imp
}

#Combine all cohorts
importance_df <- bind_rows(importance_list)

importance_summary <- importance_df %>%
  group_by(Feature) %>%
  summarise(
    MeanGain = mean(Gain),
    SDGain = sd(Gain)
  ) %>%
  arrange(desc(MeanGain))

importance_summary
library(ggplot2)

ggplot(importance_summary,
       aes(x = reorder(Feature, MeanGain), y = MeanGain)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  theme_minimal() +
  ylab("Mean Importance (Gain)") +
  xlab("") +
  ggtitle("Feature Importance Across Cohorts")
