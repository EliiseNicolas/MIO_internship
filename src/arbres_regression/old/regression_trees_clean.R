
# ============================================================
# 01 - Parameters
# ============================================================

library(rpart)
library(rpart.plot)
library(ggplot2)
library(randomForest)
library(mgcv)
library(xgboost)
library(dplyr)

rm(list = ls())

freq <- 18
diurnal_period <- 3                 # 3: day, 1: night
dp <- "day"
split <- "LOYO"                     # "LOYO" or "RS (80/20)"
data_leakage_spatio_temp <- TRUE
dist_thr <- 0.5
pigment_type <- "chla_ratio"        # "chla_ratio", "total_ratio", "conc"
n_trees <- 300
log <- TRUE
FOD_0 <- FALSE
model_type <- "RF"                  # "CART", "RF", "XGBOOST"

n_folds <- 5
n_repeats <- 10
set.seed(123)

thresholds <- c(0.25, 0.5, 0.75, 1, 1.25, 1.5)

path_datas <- paste0("F:/data_elise/ds_NASC_pig_ftle_fod/ds_NASC_mean_pig_grid_all/ds_NASC_mean_pig_grid_pig_ftle_fod_2018_2021_2023_transect_", freq, "kHz_mask9.rds")

datas <- readRDS(path_datas)
str(datas)


# ============================================================
# 02 - Data preparation
# ============================================================

# Filter day/night
datas <- datas[datas$day == diurnal_period, ]
str(datas)


# Filter NASC extreme values
print(range(datas$nasc))

q <- quantile(datas$nasc, probs = c(0.05, 0.95), na.rm = TRUE)
datas <- datas |> filter(nasc >= q[1], nasc <= q[2])

print(range(datas$nasc))


# Log NASC
if (log) datas$nasc <- log10(datas$nasc)

print(range(datas$nasc))


# NASC distribution
if (log) {
  nasc_scale <- "log10(NASC)"
  x_label <- "log10(NASC)"
  p <- ggplot(datas, aes(x = nasc, y = after_stat(count / sum(count) * 100))) + geom_histogram(bins = 200)
} else {
  nasc_scale <- "NASC"
  x_label <- "NASC"
  p <- ggplot(datas, aes(x = nasc, y = after_stat(count / sum(count) * 100))) + geom_histogram(bins = 200)
}

p +
  theme_bw() +
  labs(x = x_label, y = "Percentage (%)", title = paste(nasc_scale, "distribution - 2018, 2021, 2023 -", dp, "-", freq, "kHz"), subtitle = paste(nasc_scale, "without extreme values (0.05 - 0.95)"))


# Build regression dataset
df <- data.frame(NASC = datas$nasc, year = format(datas$time_nasc, "%Y"), time = datas$time_nasc, fod = datas$fod, ftle = datas$ftle)

if (pigment_type == "chla_ratio") {
  print("Chla")
  vars_num <- c("NASC", "per_ratio_chla", "but_ratio_chla", "fuco_ratio_chla", "hex_ratio_chla", "allo_ratio_chla", "zea_ratio_chla", "chlb_ratio_chla", "total_chla", "ftle")
  df$total_chla <- datas$Chla
  df$per_ratio_chla <- datas$Per_Chla
  df$but_ratio_chla <- datas$But_Chla
  df$fuco_ratio_chla <- datas$Fuco_Chla
  df$hex_ratio_chla <- datas$Hex_Chla
  df$allo_ratio_chla <- datas$Allo_Chla
  df$zea_ratio_chla <- datas$Zea_Chla
  df$chlb_ratio_chla <- datas$Chlb_Chla
}

if (pigment_type == "total_ratio") {
  print("total")
  vars_num <- c("NASC", "chla_ratio_total", "per_ratio_total", "but_ratio_total", "fuco_ratio_total", "hex_ratio_total", "allo_ratio_total", "zea_ratio_total", "chlb_ratio_total", "total_pig", "ftle")
  df$chla_ratio_total <- datas$Chla_total
  df$per_ratio_total <- datas$Per_total
  df$but_ratio_total <- datas$But_total
  df$fuco_ratio_total <- datas$Fuco_total
  df$hex_ratio_total <- datas$Hex_total
  df$allo_ratio_total <- datas$Allo_total
  df$zea_ratio_total <- datas$Zea_total
  df$chlb_ratio_total <- datas$Chlb_total
  df$total_pig <- datas$total_pig
}

if (pigment_type == "conc") {
  print("conc")
  vars_num <- c("NASC", "chla", "per", "but", "fuco", "hex", "allo", "zea", "chlb", "total_pig", "ftle")
  df$chla <- datas$Chla
  df$per <- datas$Per
  df$but <- datas$But
  df$fuco <- datas$Fuco
  df$hex <- datas$Hex
  df$allo <- datas$Allo
  df$zea <- datas$Zea
  df$chlb <- datas$Chlb
  df$total_pig <- datas$total_pig
}

str(df)


# ============================================================
# 03 - FOD & missing data
# ============================================================

if (FOD_0) {
  fod_0 <- as.character(datas$fod)
  fod_0[as.numeric(fod_0) > 6] <- "0"
  datas$fod <- fod_0
}

df$fod <- as.factor(df$fod)
df$fod[df$fod == "NA"] <- NA
df$fod <- droplevels(df$fod)

table(df$fod, useNA = "ifany")

cat("Nombre total de données avant filtrage :", nrow(df), "\n")
print(colSums(is.na(df)))

df <- df |> filter(if_all(all_of(setdiff(vars_num, "NASC")), ~ !is.na(.)), !is.na(fod), fod != "NA")

cat("Nombre total de données après filtrage :", nrow(df), "\n")
print(colSums(is.na(df)))


# ============================================================
# 04 - Spatio-temporal leakage
# ============================================================

remove_close_pairs <- function(df_input, vars_num, dist_thr, time_thr = 3600) {
  df_work <- df_input
  n_initial <- nrow(df_work)
  n_pairs_initial <- NA
  n_iterations <- 0
  removed_rows <- integer(0)
  pairs_history <- data.frame(iteration = integer(0), n_pairs = integer(0), removed = integer(0))
  
  repeat {
    if (nrow(df_work) < 2) break
    
    df_scaled <- scale(df_work[, vars_num])
    dist_mat <- as.matrix(dist(df_scaled))
    time_diff <- abs(outer(as.numeric(df_work$time), as.numeric(df_work$time), FUN = "-"))
    
    pairs <- which(time_diff <= time_thr & dist_mat < dist_thr, arr.ind = TRUE)
    pairs <- pairs[pairs[, 1] < pairs[, 2], , drop = FALSE]
    n_pairs <- nrow(pairs)
    
    if (is.na(n_pairs_initial)) n_pairs_initial <- n_pairs
    if (n_pairs == 0) break
    
    pair <- pairs[1, ]
    remove_row <- pair[2]
    original_row <- as.integer(rownames(df_work)[remove_row])
    
    removed_rows <- c(removed_rows, original_row)
    n_iterations <- n_iterations + 1
    pairs_history <- rbind(pairs_history, data.frame(iteration = n_iterations, n_pairs = n_pairs, removed = original_row))
    
    df_work <- df_work[-remove_row, , drop = FALSE]
  }
  
  list(df = df_work, n_initial = n_initial, n_remaining = nrow(df_work), n_removed = n_initial - nrow(df_work), n_pairs_initial = ifelse(is.na(n_pairs_initial), 0, n_pairs_initial), n_pairs_final = 0, n_iterations = n_iterations, removed_rows = removed_rows, pairs_history = pairs_history)
}


if (data_leakage_spatio_temp) {
  
  cat("Lignes dupliquées :", sum(duplicated(df)), "\n")
  
  df_scaled <- scale(df[, vars_num])
  dist_mat <- as.matrix(dist(df_scaled))
  time_diff <- abs(outer(as.numeric(df$time), as.numeric(df$time), FUN = "-"))
  
  dist_mat_plot <- dist_mat
  dist_mat_plot[time_diff > 3600] <- NA
  dist_mat_plot[lower.tri(dist_mat_plot, diag = TRUE)] <- NA
  dist_close <- dist_mat_plot[time_diff <= 3600 & !is.na(dist_mat_plot)]
  
  print(quantile(dist_close, probs = c(0.01, 0.05, 0.10, 0.25, 0.50, 0.75), na.rm = TRUE))
  
  ggplot(data.frame(distance = dist_close), aes(x = distance)) +
    geom_histogram(bins = 100) +
    geom_vline(xintercept = dist_thr, linetype = "dashed") +
    theme_bw() +
    labs(x = "Distance", y = "Count", title = "Distribution of distances for observations < 1 h apart", subtitle = paste("Selected threshold:", dist_thr))
}


# ============================================================
# 05 - Leakage sensitivity
# ============================================================

sensitivity_results <- data.frame()
sensitivity_data <- list()

if (data_leakage_spatio_temp) {
  
  for (thr in thresholds) {
    
    cat("\n========================================\n")
    cat("Threshold :", thr, "\n")
    cat("========================================\n")
    
    leakage_thr <- remove_close_pairs(df, vars_num, thr)
    df_thr <- leakage_thr$df
    
    cat("Paires initiales :", leakage_thr$n_pairs_initial, "\n")
    cat("Observations supprimées :", leakage_thr$n_removed, "\n")
    cat("Observations restantes :", leakage_thr$n_remaining, "\n")
    cat("Nombre d'itérations :", leakage_thr$n_iterations, "\n")
    
    years_all_thr <- sort(unique(df_thr$year))
    loyo_thr <- data.frame()
    
    for (yr in years_all_thr) {
      
      train_yr <- df_thr[df_thr$year != yr, ]
      test_yr <- df_thr[df_thr$year == yr, ]
      
      train_yr$year <- NULL
      test_yr$year <- NULL
      train_yr$time <- NULL
      test_yr$time <- NULL
      
      train_yr_scaled <- scale(train_yr[, vars_num])
      center_yr <- attr(train_yr_scaled, "scaled:center")
      scale_yr <- attr(train_yr_scaled, "scaled:scale")
      test_yr_scaled <- scale(test_yr[, vars_num], center = center_yr, scale = scale_yr)
      
      ds_train_yr <- train_yr
      ds_train_yr[, vars_num] <- train_yr_scaled
      
      ds_test_yr <- test_yr
      ds_test_yr[, vars_num] <- test_yr_scaled
      
      if (model_type == "CART") model_yr <- rpart(NASC ~ ., data = ds_train_yr, method = "anova")
      
      if (model_type == "RF") model_yr <- randomForest(NASC ~ ., data = ds_train_yr, ntree = n_trees, importance = TRUE)
      
      if (model_type == "XGBOOST") {
        X_train <- as.matrix(ds_train_yr[, setdiff(vars_num, "NASC")])
        y_train <- ds_train_yr$NASC
        X_test <- as.matrix(ds_test_yr[, setdiff(vars_num, "NASC")])
        model_yr <- xgboost(data = X_train, label = y_train, nrounds = 100, max_depth = 2, eta = 0.05, min_child_weight = 3, subsample = 0.8, colsample_bytree = 0.8, objective = "reg:squarederror", verbose = 0)
        pred_yr <- predict(model_yr, X_test)
      } else {
        pred_yr <- predict(model_yr, newdata = ds_test_yr)
      }
      
      obs_yr <- ds_test_yr$NASC
      RMSE_yr <- sqrt(mean((pred_yr - obs_yr)^2))
      R2_yr <- 1 - sum((obs_yr - pred_yr)^2) / sum((obs_yr - mean(obs_yr))^2)
      MAE_yr <- mean(abs(pred_yr - obs_yr))
      
      loyo_thr <- rbind(loyo_thr, data.frame(year = yr, RMSE = RMSE_yr, R2 = R2_yr, MAE = MAE_yr, n_train = nrow(train_yr), n_test = nrow(test_yr)))
    }
    
    sensitivity_results <- rbind(sensitivity_results, data.frame(threshold = thr, n_pairs_initial = leakage_thr$n_pairs_initial, n_removed = leakage_thr$n_removed, percent_removed = 100 * leakage_thr$n_removed / leakage_thr$n_initial, n_remaining = leakage_thr$n_remaining, n_iterations = leakage_thr$n_iterations, RMSE = mean(loyo_thr$RMSE), RMSE_sd = sd(loyo_thr$RMSE), R2 = mean(loyo_thr$R2), R2_sd = sd(loyo_thr$R2), MAE = mean(loyo_thr$MAE), MAE_sd = sd(loyo_thr$MAE)))
    
    sensitivity_data[[as.character(thr)]] <- leakage_thr
  }
  
  print(sensitivity_results)
  
  ggplot(sensitivity_results, aes(x = threshold, y = RMSE)) +
    geom_line() +
    geom_point(size = 3) +
    geom_vline(xintercept = dist_thr, linetype = "dashed") +
    theme_bw() +
    labs(x = "Distance threshold", y = "RMSE", title = paste("Sensitivity analysis -", model_type), subtitle = paste("LOYO -", dp, "-", freq, "kHz"))
  
  ggplot(sensitivity_results, aes(x = threshold, y = R2)) +
    geom_line() +
    geom_point(size = 3) +
    geom_vline(xintercept = dist_thr, linetype = "dashed") +
    theme_bw() +
    labs(x = "Distance threshold", y = "R²", title = paste("Sensitivity analysis -", model_type), subtitle = paste("LOYO -", dp, "-", freq, "kHz"))
  
  ggplot(sensitivity_results, aes(x = threshold, y = percent_removed)) +
    geom_line() +
    geom_point(size = 3) +
    geom_vline(xintercept = dist_thr, linetype = "dashed") +
    theme_bw() +
    labs(x = "Distance threshold", y = "Observations removed (%)", title = "Sensitivity of data removal to distance threshold", subtitle = paste("Observations < 1 hour apart -", dp, "-", freq, "kHz"))
  
  leakage_final <- remove_close_pairs(df, vars_num, dist_thr)
  df <- leakage_final$df
  
  cat("\n========================================\n")
  cat("NETTOYAGE FINAL\n")
  cat("========================================\n")
  cat("Threshold :", dist_thr, "\n")
  cat("Paires initiales :", leakage_final$n_pairs_initial, "\n")
  cat("Observations initiales :", leakage_final$n_initial, "\n")
  cat("Observations supprimées :", leakage_final$n_removed, "\n")
  cat("Pourcentage supprimé :", round(100 * leakage_final$n_removed / leakage_final$n_initial, 2), "%\n")
  cat("Observations restantes :", leakage_final$n_remaining, "\n")
  cat("Nombre d'itérations :", leakage_final$n_iterations, "\n")
  
  df_scaled <- scale(df[, vars_num])
  dist_mat <- as.matrix(dist(df_scaled))
  time_diff <- abs(outer(as.numeric(df$time), as.numeric(df$time), FUN = "-"))
  
  final_pairs <- which(time_diff <= 3600 & dist_mat < dist_thr, arr.ind = TRUE)
  final_pairs <- final_pairs[final_pairs[, 1] < final_pairs[, 2], , drop = FALSE]
  
  cat("Paires restantes après nettoyage :", nrow(final_pairs), "\n")
  
  if (nrow(final_pairs) == 0) {
    cat("OK : aucune paire répondant aux critères de leakage ne reste.\n")
  } else {
    warning("Des paires de leakage sont encore présentes.")
  }
}


# ============================================================
# 06 - LOYO cross-validation
# ============================================================

years_all <- sort(unique(df$year))
loyo_results <- data.frame()
loyo_predictions <- data.frame()

for (yr in years_all) {
  
  train_yr <- df[df$year != yr, ]
  test_yr <- df[df$year == yr, ]
  
  train_yr$year <- NULL
  test_yr$year <- NULL
  train_yr$time <- NULL
  test_yr$time <- NULL
  
  train_yr_scaled <- scale(train_yr[, vars_num])
  center_yr <- attr(train_yr_scaled, "scaled:center")
  scale_yr <- attr(train_yr_scaled, "scaled:scale")
  test_yr_scaled <- scale(test_yr[, vars_num], center = center_yr, scale = scale_yr)
  
  ds_train_yr <- train_yr
  ds_train_yr[, vars_num] <- train_yr_scaled
  
  ds_test_yr <- test_yr
  ds_test_yr[, vars_num] <- test_yr_scaled
  
  if (model_type == "CART") model_yr <- rpart(NASC ~ ., data = ds_train_yr, method = "anova")
  
  if (model_type == "RF") model_yr <- randomForest(NASC ~ ., data = ds_train_yr, ntree = n_trees, importance = TRUE)
  
  if (model_type == "XGBOOST") {
    X_train <- as.matrix(ds_train_yr[, setdiff(vars_num, "NASC")])
    y_train <- ds_train_yr$NASC
    X_test <- as.matrix(ds_test_yr[, setdiff(vars_num, "NASC")])
    model_yr <- xgboost(data = X_train, label = y_train, nrounds = 100, max_depth = 2, eta = 0.05, min_child_weight = 3, subsample = 0.8, colsample_bytree = 0.8, objective = "reg:squarederror", verbose = 0)
    pred_yr <- predict(model_yr, X_test)
  } else {
    pred_yr <- predict(model_yr, newdata = ds_test_yr)
  }
  
  obs_yr <- ds_test_yr$NASC
  
  RMSE_yr <- sqrt(mean((pred_yr - obs_yr)^2))
  R2_yr <- 1 - sum((obs_yr - pred_yr)^2) / sum((obs_yr - mean(obs_yr))^2)
  MAE_yr <- mean(abs(pred_yr - obs_yr))
  
  loyo_results <- rbind(loyo_results, data.frame(year = yr, n_train = nrow(train_yr), n_test = nrow(test_yr), RMSE = RMSE_yr, R2 = R2_yr, MAE = MAE_yr))
  
  loyo_predictions <- rbind(loyo_predictions, data.frame(year = yr, observed = obs_yr, predicted = pred_yr))
}

print(loyo_results)

cat("\n--- LOYO : moyenne sur toutes les années ---\n")
cat("RMSE moyen :", mean(loyo_results$RMSE), "(sd =", sd(loyo_results$RMSE), ")\n")
cat("R² moyen :", mean(loyo_results$R2), "(sd =", sd(loyo_results$R2), ")\n")
cat("MAE moyen :", mean(loyo_results$MAE), "(sd =", sd(loyo_results$MAE), ")\n")

ggplot(loyo_results, aes(x = year, y = RMSE)) +
  geom_col() +
  theme_bw() +
  labs(x = "Année test", y = "RMSE", title = paste("RMSE par année - LOYO -", model_type), subtitle = paste("Transect 2018, 2021, 2023,", dp, ",", freq, "kHz, pigments:", pigment_type))

ggplot(loyo_predictions, aes(x = observed, y = predicted, color = year)) +
  geom_point(alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  theme_bw() +
  labs(x = "True NASC", y = "Predicted NASC", title = paste("Real NASC vs. Predicted NASC - LOYO -", model_type), color = "Année test")


# ============================================================
# 07 - Random k-fold cross-validation
# ============================================================

n_folds <- 5
set.seed(123)

fold_assign <- sample(rep(seq_len(n_folds), length.out = nrow(df)))

rs_results <- data.frame()
rs_predictions <- data.frame()

for (k in seq_len(n_folds)) {
  
  train_k <- df[fold_assign != k, ]
  test_k <- df[fold_assign == k, ]
  
  train_k$year <- NULL
  test_k$year <- NULL
  train_k$time <- NULL
  test_k$time <- NULL
  
  train_k_scaled <- scale(train_k[, vars_num])
  center_k <- attr(train_k_scaled, "scaled:center")
  scale_k <- attr(train_k_scaled, "scaled:scale")
  test_k_scaled <- scale(test_k[, vars_num], center = center_k, scale = scale_k)
  
  ds_train_k <- train_k
  ds_train_k[, vars_num] <- train_k_scaled
  
  ds_test_k <- test_k
  ds_test_k[, vars_num] <- test_k_scaled
  
  if (model_type == "CART") model_k <- rpart(NASC ~ ., data = ds_train_k, method = "anova")
  
  if (model_type == "RF") model_k <- randomForest(NASC ~ ., data = ds_train_k, ntree = n_trees, importance = TRUE)
  
  if (model_type == "XGBOOST") {
    X_train <- as.matrix(ds_train_k[, setdiff(vars_num, "NASC")])
    y_train <- ds_train_k$NASC
    X_test <- as.matrix(ds_test_k[, setdiff(vars_num, "NASC")])
    model_k <- xgboost(data = X_train, label = y_train, nrounds = 100, max_depth = 2, eta = 0.05, min_child_weight = 3, subsample = 0.8, colsample_bytree = 0.8, objective = "reg:squarederror", verbose = 0)
    pred_k <- predict(model_k, X_test)
  } else {
    pred_k <- predict(model_k, newdata = ds_test_k)
  }
  
  obs_k <- ds_test_k$NASC
  
  RMSE_k <- sqrt(mean((pred_k - obs_k)^2))
  R2_k <- 1 - sum((obs_k - pred_k)^2) / sum((obs_k - mean(obs_k))^2)
  MAE_k <- mean(abs(pred_k - obs_k))
  
  rs_results <- rbind(rs_results, data.frame(fold = k, n_train = nrow(train_k), n_test = nrow(test_k), RMSE = RMSE_k, R2 = R2_k, MAE = MAE_k))
  
  rs_predictions <- rbind(rs_predictions, data.frame(fold = k, observed = obs_k, predicted = pred_k))
}

print(rs_results)

cat("\n--- Random Split : moyenne sur les", n_folds, "folds ---\n")
cat("RMSE moyen :", mean(rs_results$RMSE), "(sd =", sd(rs_results$RMSE), ")\n")
cat("R² moyen :", mean(rs_results$R2), "(sd =", sd(rs_results$R2), ")\n")
cat("MAE moyen :", mean(rs_results$MAE), "(sd =", sd(rs_results$MAE), ")\n")

ggplot(rs_results, aes(x = factor(fold), y = RMSE)) +
  geom_col() +
  theme_bw() +
  labs(x = "Fold", y = "RMSE", title = paste("RMSE par fold - Random Split (k =", n_folds, ") -", model_type), subtitle = paste("Transect 2018, 2021, 2023,", dp, ",", freq, "kHz, pigments:", pigment_type))

ggplot(rs_predictions, aes(x = observed, y = predicted, color = factor(fold))) +
  geom_point(alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  theme_bw() +
  labs(x = "True NASC", y = "Predicted NASC", title = paste("Real NASC vs. Predicted NASC - Random Split -", model_type), color = "Fold")


# ============================================================
# 08 - Final model
# ============================================================

if (split == "RS (80/20)") {
  
  train <- df[train_index, ]
  test <- df[test_index, ]
  
} else {
  
  test_year <- names(which.min(table(df$year)))
  train <- df[df$year != test_year, ]
  test <- df[df$year == test_year, ]
}

train$year <- NULL
test$year <- NULL
train$time <- NULL
test$time <- NULL

train_scaled <- scale(train[, vars_num])
train_center <- attr(train_scaled, "scaled:center")
train_scale <- attr(train_scaled, "scaled:scale")
test_scaled <- scale(test[, vars_num], center = train_center, scale = train_scale)

ds_train_scaled <- train
ds_train_scaled[, vars_num] <- train_scaled

ds_test_scaled <- test
ds_test_scaled[, vars_num] <- test_scaled

cat("\n--- Dimensions finales ---\n")
cat("Train :", nrow(ds_train_scaled), "observations\n")
cat("Test :", nrow(ds_test_scaled), "observations\n")


if (model_type == "CART") {
  
  model <- rpart(NASC ~ ., data = ds_train_scaled, method = "anova")
  
  rpart.plot(model, extra = 101, main = paste("Regression tree - CART -", split, "- scaled data"))
  
}


if (model_type == "RF") {
  
  model <- randomForest(NASC ~ ., data = ds_train_scaled, ntree = n_trees, importance = TRUE)
  
  print(importance(model))
  
}


if (model_type == "XGBOOST") {
  
  X_train <- as.matrix(ds_train_scaled[, setdiff(vars_num, "NASC")])
  y_train <- ds_train_scaled$NASC
  
  X_test <- as.matrix(ds_test_scaled[, setdiff(vars_num, "NASC")])
  
  model <- xgboost(data = X_train, label = y_train, nrounds = 100, max_depth = 2, eta = 0.05, min_child_weight = 3, subsample = 0.8, colsample_bytree = 0.8, objective = "reg:squarederror", verbose = 0)
  
}


# ============================================================
# 09 - Predictions & metrics
# ============================================================

observed <- ds_test_scaled$NASC

if (model_type == "XGBOOST") {
  prediction <- predict(model, X_test)
} else {
  prediction <- predict(model, newdata = ds_test_scaled)
}

resultats <- data.frame(NASC_reel = observed, NASC_predit = prediction)

RMSE <- sqrt(mean((prediction - observed)^2))
R2 <- 1 - sum((observed - prediction)^2) / sum((observed - mean(observed))^2)
MAE <- mean(abs(prediction - observed))

cat("\n--- Final model performance ---\n")
cat("RMSE :", RMSE, "\n")
cat("R² :", R2, "\n")
cat("MAE :", MAE, "\n")


# ============================================================
# 10 - Variable importance
# ============================================================

if (model_type == "CART") {
  
  importance <- model$variable.importance
  importance_pct <- importance / sum(importance) * 100
  
  importance_df <- data.frame(variable = names(importance_pct), importance = as.numeric(importance_pct))
  importance_df <- importance_df[order(importance_df$importance), ]
  
}

if (model_type == "RF") {
  
  importance_df <- data.frame(variable = rownames(importance(model)), importance = importance(model)[, "%IncMSE"])
  importance_df <- importance_df[order(importance_df$importance), ]
  
}

if (model_type == "XGBOOST") {
  
  importance_df <- xgb.importance(model = model)
  importance_df <- data.frame(variable = importance_df$Feature, importance = importance_df$Gain * 100)
  importance_df <- importance_df[order(importance_df$importance), ]
  
}

if (model_type != "XGBOOST") {
  
  ggplot(importance_df, aes(x = reorder(variable, importance), y = importance)) +
    geom_col() +
    coord_flip() +
    theme_bw() +
    labs(x = "Variable", y = "Importance (%)", title = paste("Variable importance for NASC prediction -", model_type), subtitle = paste("Transect 2018, 2021, 2023,", dp, ",", freq, "kHz, split:", split, "- scaled data\nPigment type:", pigment_type))
  
} else {
  
  ggplot(importance_df, aes(x = reorder(variable, importance), y = importance)) +
    geom_col() +
    coord_flip() +
    theme_bw() +
    labs(x = "Variable", y = "Importance (%)", title = paste("Variable importance for NASC prediction -", model_type), subtitle = paste("Transect 2018, 2021, 2023,", dp, ",", freq, "kHz\nPigment type:", pigment_type))
  
}


# ============================================================
# 11 - Prediction plot
# ============================================================

lim <- range(c(observed, prediction), na.rm = TRUE)

plot(observed, prediction, xlim = lim, ylim = lim, xlab = "True NASC", ylab = "Predicted NASC", main = paste("Real NASC vs. Predicted NASC -", model_type), pch = 4, cex = 0.5)

mtext(paste("Transect 2018, 2021, 2023,", dp, ",", freq, "kHz, split:", split, "- scaled data\nPigment type:", pigment_type), side = 3, line = 0.5, cex = 0.9)

abline(0, 1, col = "red")

legend("bottomright", legend = c(paste0("Predictions (RMSE = ", round(RMSE, 2), ", R² = ", round(R2, 2), ")"), "Identity line (y = x)"), col = c("black", "red"), pch = c(4, NA), lty = c(NA, 1))


# ============================================================
# 12 - Learning curve
# ============================================================

train_fractions <- seq(0.1, 1, by = 0.1)
learning_curve <- data.frame()

set.seed(123)

for (frac in train_fractions) {
  
  n_sub <- floor(frac * nrow(ds_train_scaled))
  sub_index <- sample(seq_len(nrow(ds_train_scaled)), size = n_sub)
  train_sub <- ds_train_scaled[sub_index, ]
  
  if (model_type == "CART") model_sub <- rpart(NASC ~ ., data = train_sub, method = "anova")
  
  if (model_type == "RF") model_sub <- randomForest(NASC ~ ., data = train_sub, ntree = n_trees, importance = TRUE)
  
  if (model_type == "XGBOOST") {
    X_train <- as.matrix(train_sub[, setdiff(vars_num, "NASC")])
    y_train <- train_sub$NASC
    X_test <- as.matrix(ds_test_scaled[, setdiff(vars_num, "NASC")])
    model_sub <- xgboost(data = X_train, label = y_train, nrounds = 100, max_depth = 2, eta = 0.05, min_child_weight = 3, subsample = 0.8, colsample_bytree = 0.8, objective = "reg:squarederror", verbose = 0)
    pred_train_sub <- predict(model_sub, X_train)
    pred_test_sub <- predict(model_sub, X_test)
  } else {
    pred_train_sub <- predict(model_sub, newdata = train_sub)
    pred_test_sub <- predict(model_sub, newdata = ds_test_scaled)
  }
  
  RMSE_train_sub <- sqrt(mean((pred_train_sub - train_sub$NASC)^2))
  RMSE_test_sub <- sqrt(mean((pred_test_sub - ds_test_scaled$NASC)^2))
  
  learning_curve <- rbind(learning_curve, data.frame(n_train = n_sub, RMSE_train = RMSE_train_sub, RMSE_test = RMSE_test_sub))
}

print(learning_curve)

learning_curve_long <- rbind(
  data.frame(n_train = learning_curve$n_train, RMSE = learning_curve$RMSE_train, set = "Train"),
  data.frame(n_train = learning_curve$n_train, RMSE = learning_curve$RMSE_test, set = "Test")
)

ggplot(learning_curve_long, aes(x = n_train, y = RMSE, color = set)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  labs(x = "Taille de l'échantillon d'entraînement", y = "RMSE", title = paste("Courbe d'apprentissage -", model_type), subtitle = paste("Transect 2018, 2021, 2023,", dp, ",", freq, "kHz, split:", split), color = "Dataset")


# ============================================================
# 13 - Validation curve
# ============================================================

if (model_type == "RF") {
  
  hp_values <- c(10, 50, 100, 200, 300, 500, 800)
  hp_label <- "Nombre d'arbres (ntree)"
  validation_curve <- data.frame()
  
  for (hp in hp_values) {
    
    model_val <- randomForest(NASC ~ ., data = ds_train_scaled, ntree = hp, importance = TRUE)
    
    pred_train_val <- predict(model_val, newdata = ds_train_scaled)
    pred_test_val <- predict(model_val, newdata = ds_test_scaled)
    
    RMSE_train_val <- sqrt(mean((pred_train_val - ds_train_scaled$NASC)^2))
    RMSE_test_val <- sqrt(mean((pred_test_val - ds_test_scaled$NASC)^2))
    
    validation_curve <- rbind(validation_curve, data.frame(hyperparameter = hp, RMSE_train = RMSE_train_val, RMSE_test = RMSE_test_val))
  }
}


if (model_type == "CART") {
  
  hp_values <- c(0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2)
  hp_label <- "Complexité (cp)"
  validation_curve <- data.frame()
  
  for (hp in hp_values) {
    
    model_val <- rpart(NASC ~ ., data = ds_train_scaled, method = "anova", control = rpart.control(cp = hp))
    
    pred_train_val <- predict(model_val, newdata = ds_train_scaled)
    pred_test_val <- predict(model_val, newdata = ds_test_scaled)
    
    RMSE_train_val <- sqrt(mean((pred_train_val - ds_train_scaled$NASC)^2))
    RMSE_test_val <- sqrt(mean((pred_test_val - ds_test_scaled$NASC)^2))
    
    validation_curve <- rbind(validation_curve, data.frame(hyperparameter = hp, RMSE_train = RMSE_train_val, RMSE_test = RMSE_test_val))
  }
}


if (model_type == "XGBOOST") {
  
  hp_values <- c(10, 50, 100, 200, 300, 500)
  hp_label <- "Nombre d'itérations (nrounds)"
  validation_curve <- data.frame()
  
  X_train <- as.matrix(ds_train_scaled[, setdiff(vars_num, "NASC")])
  y_train <- ds_train_scaled$NASC
  X_test <- as.matrix(ds_test_scaled[, setdiff(vars_num, "NASC")])
  
  for (hp in hp_values) {
    
    model_val <- xgboost(data = X_train, label = y_train, nrounds = hp, max_depth = 2, eta = 0.05, min_child_weight = 3, subsample = 0.8, colsample_bytree = 0.8, objective = "reg:squarederror", verbose = 0)
    
    pred_train_val <- predict(model_val, X_train)
    pred_test_val <- predict(model_val, X_test)
    
    RMSE_train_val <- sqrt(mean((pred_train_val - ds_train_scaled$NASC)^2))
    RMSE_test_val <- sqrt(mean((pred_test_val - ds_test_scaled$NASC)^2))
    
    validation_curve <- rbind(validation_curve, data.frame(hyperparameter = hp, RMSE_train = RMSE_train_val, RMSE_test = RMSE_test_val))
  }
}

validation_curve_long <- rbind(
  data.frame(hyperparameter = validation_curve$hyperparameter, RMSE = validation_curve$RMSE_train, set = "Train"),
  data.frame(hyperparameter = validation_curve$hyperparameter, RMSE = validation_curve$RMSE_test, set = "Test")
)

ggplot(validation_curve_long, aes(x = hyperparameter, y = RMSE, color = set)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  labs(x = hp_label, y = "RMSE", title = paste("Courbe de validation -", model_type), subtitle = paste("Transect 2018, 2021, 2023,", dp, ",", freq, "kHz, split:", split), color = "Dataset")
