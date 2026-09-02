# ============================================================
# RANDOM FOREST vs XGBOOST
# + SPATIO-TEMPORAL LEAKAGE REMOVAL
# + LOYO
# + RANDOM 80/20
# + FOD_0 TRUE/FALSE
# + LOG TRUE/FALSE
# + FOD AS PREDICTOR
# ============================================================

library(randomForest)
library(xgboost)
library(dplyr)
library(ggplot2)

set.seed(123)


# ============================================================
# 01 - PARAMETRES
# ============================================================

freq <- 120
diurnal_period <- 3
dp <- "day"

pigment_type <- "chla_ratio"

n_repeats <- 10
train_fraction <- 0.8

# Leakage
time_thr <- 3600

distance_thresholds <- c(0.25, 0.5, 0.75, 1)

# FOD / log
FOD_values <- c(FALSE, TRUE)
log_values <- c(FALSE, TRUE)


# ============================================================
# RANDOM FOREST - HYPERPARAMETRES
# ============================================================

ntree_values <- c(100, 300, 500)

nodesize_values <- c(1, 3, 5, 10)

mtry_fraction <- c(0.5, 1, 1.5)


# ============================================================
# XGBOOST - HYPERPARAMETRES
# ============================================================

max_depth_values <- c(3, 5, 7)

eta_values <- c(0.03, 0.1, 0.3)

min_child_weight_values <- c(1, 3, 5)

subsample_values <- c(0.8, 1)

colsample_bytree_values <- c(0.8, 1)

nrounds <- 500


# ============================================================
# DONNEES
# ============================================================

path_datas <- paste0("F:/data_elise/ds_NASC_pig_ftle_fod/ds_NASC_mean_pig_grid_all/ds_NASC_mean_pig_grid_pig_ftle_fod_2018_2021_2023_transect_", freq, "kHz_mask9.rds")

datas <- readRDS(path_datas)


# ============================================================
# 02 - PREPARATION DES DONNEES
# ============================================================

prepare_data <- function(datas, pigment_type = "chla_ratio", FOD_0 = FALSE, log = TRUE, diurnal_period = 3) {
  
  dat <- datas
  
  # ----------------------------------------------------------
  # Day / night
  # ----------------------------------------------------------
  
  dat <- dat[dat$day == diurnal_period, ]
  
  # ----------------------------------------------------------
  # NASC extreme values
  # ----------------------------------------------------------
  
  q <- quantile(dat$nasc, probs = c(0.05, 0.95), na.rm = TRUE)
  
  dat <- dat |>
    filter(nasc >= q[1], nasc <= q[2])
  
  # ----------------------------------------------------------
  # Log NASC
  # ----------------------------------------------------------
  
  if (log) {
    
    # Vérification des valeurs avant log10
    if (any(!is.finite(dat$nasc))) {
      stop("ERREUR : NASC contient des valeurs non finies avant log10().")
    }
    
    if (any(dat$nasc <= 0, na.rm = TRUE)) {
      stop("ERREUR : NASC contient des valeurs <= 0. Impossible d'appliquer log10().")
    }
    
    cat("  Vérification LOG : toutes les valeurs NASC sont > 0 -> OK\n")
    
    dat$nasc <- log10(dat$nasc)
    
    # Vérification après transformation
    if (any(!is.finite(dat$nasc))) {
      stop("ERREUR : log10(NASC) produit des valeurs non finies.")
    }
    
    cat("  Vérification LOG : log10(NASC) contient uniquement des valeurs finies -> OK\n")
  }
  
  # ----------------------------------------------------------
  # FOD_0
  # ----------------------------------------------------------
  
  if (FOD_0) {
    
    fod_0 <- as.character(dat$fod)
    
    fod_num <- suppressWarnings(as.numeric(fod_0))
    
    fod_0[!is.na(fod_num) & fod_num > 6] <- "0"
    
    dat$fod <- fod_0
  }
  
  # ----------------------------------------------------------
  # Base dataframe
  # ----------------------------------------------------------
  
  df <- data.frame(
    NASC = dat$nasc,
    year = format(dat$time_nasc, "%Y"),
    time = dat$time_nasc,
    fod = dat$fod,
    ftle = dat$ftle
  )
  
  # ----------------------------------------------------------
  # Pigments
  # ----------------------------------------------------------
  
  if (pigment_type == "chla_ratio") {
    
    vars_num <- c("NASC", "per_ratio_chla", "but_ratio_chla", "fuco_ratio_chla", "hex_ratio_chla", "allo_ratio_chla", "zea_ratio_chla", "chlb_ratio_chla", "total_chla", "ftle")
    
    df$total_chla <- dat$Chla
    df$per_ratio_chla <- dat$Per_Chla
    df$but_ratio_chla <- dat$But_Chla
    df$fuco_ratio_chla <- dat$Fuco_Chla
    df$hex_ratio_chla <- dat$Hex_Chla
    df$allo_ratio_chla <- dat$Allo_Chla
    df$zea_ratio_chla <- dat$Zea_Chla
    df$chlb_ratio_chla <- dat$Chlb_Chla
  }
  
  if (pigment_type == "total_ratio") {
    
    vars_num <- c("NASC", "chla_ratio_total", "per_ratio_total", "but_ratio_total", "fuco_ratio_total", "hex_ratio_total", "allo_ratio_total", "zea_ratio_total", "chlb_ratio_total", "total_pig", "ftle")
    
    df$chla_ratio_total <- dat$Chla_total
    df$per_ratio_total <- dat$Per_total
    df$but_ratio_total <- dat$But_total
    df$fuco_ratio_total <- dat$Fuco_total
    df$hex_ratio_total <- dat$Hex_total
    df$allo_ratio_total <- dat$Allo_total
    df$zea_ratio_total <- dat$Zea_total
    df$chlb_ratio_total <- dat$Chlb_total
    df$total_pig <- dat$total_pig
  }
  
  if (pigment_type == "conc") {
    
    vars_num <- c("NASC", "chla", "per", "but", "fuco", "hex", "allo", "zea", "chlb", "total_pig", "ftle")
    
    df$chla <- dat$Chla
    df$per <- dat$Per
    df$but <- dat$But
    df$fuco <- dat$Fuco
    df$hex <- dat$Hex
    df$allo <- dat$Allo
    df$zea <- dat$Zea
    df$chlb <- dat$Chlb
    df$total_pig <- dat$total_pig
  }
  
  # ----------------------------------------------------------
  # FOD factor
  # ----------------------------------------------------------
  
  df$fod <- as.factor(df$fod)
  
  df$fod[df$fod == "NA"] <- NA
  
  df$fod <- droplevels(df$fod)
  
  # ----------------------------------------------------------
  # Missing values
  # ----------------------------------------------------------
  
  df <- df |>
    filter(if_all(all_of(setdiff(vars_num, "NASC")), ~ !is.na(.)), !is.na(fod))
  
  list(df = df, vars_num = vars_num)
}


# ============================================================
# 03 - SUPPRESSION DU DATA LEAKAGE
# ============================================================

remove_close_pairs <- function(df_input, predictor_vars, dist_thr = 0.5, time_thr = 3600) {
  
  df_work <- df_input
  
  n_initial <- nrow(df_work)
  
  removed_rows <- integer(0)
  
  n_iterations <- 0
  
  repeat {
    
    if (nrow(df_work) < 2) {
      break
    }
    
    # --------------------------------------------------------
    # Distance uniquement sur les prédicteurs numériques
    # --------------------------------------------------------
    
    df_scaled <- scale(df_work[, predictor_vars, drop = FALSE])
    
    dist_mat <- as.matrix(dist(df_scaled))
    
    # --------------------------------------------------------
    # Différence temporelle
    # --------------------------------------------------------
    
    time_diff <- abs(outer(as.numeric(df_work$time), as.numeric(df_work$time), FUN = "-"))
    
    # --------------------------------------------------------
    # Paires proches
    # --------------------------------------------------------
    
    pairs <- which(time_diff <= time_thr & dist_mat < dist_thr, arr.ind = TRUE)
    
    pairs <- pairs[pairs[, 1] < pairs[, 2], , drop = FALSE]
    
    n_pairs <- nrow(pairs)
    
    if (n_pairs == 0) {
      break
    }
    
    # --------------------------------------------------------
    # Suppression du deuxième élément
    # --------------------------------------------------------
    
    pair <- pairs[1, ]
    
    remove_row <- pair[2]
    
    removed_rows <- c(removed_rows, remove_row)
    
    n_iterations <- n_iterations + 1
    
    df_work <- df_work[-remove_row, , drop = FALSE]
  }
  
  list(
    df = df_work,
    n_initial = n_initial,
    n_remaining = nrow(df_work),
    n_removed = n_initial - nrow(df_work),
    percent_removed = 100 * (n_initial - nrow(df_work)) / n_initial,
    n_iterations = n_iterations,
    removed_rows = removed_rows
  )
}


# ============================================================
# 04 - METRIQUES
# ============================================================

calculate_metrics <- function(obs, pred) {
  
  RMSE <- sqrt(mean((pred - obs)^2))
  
  R2 <- 1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2)
  
  MAE <- mean(abs(pred - obs))
  
  data.frame(RMSE = RMSE, R2 = R2, MAE = MAE)
}


# ============================================================
# 05 - RANDOM FOREST
# ============================================================

fit_rf <- function(train, test, vars_num, mtry, nodesize, ntree) {
  
  # ----------------------------------------------------------
  # Scaling des variables numériques
  # ----------------------------------------------------------
  
  train_scaled <- scale(train[, vars_num])
  
  center <- attr(train_scaled, "scaled:center")
  
  scale_values <- attr(train_scaled, "scaled:scale")
  
  test_scaled <- scale(test[, vars_num], center = center, scale = scale_values)
  
  ds_train <- train
  
  ds_test <- test
  
  ds_train[, vars_num] <- train_scaled
  
  ds_test[, vars_num] <- test_scaled
  
  # ----------------------------------------------------------
  # Variables utilisées par le modèle
  # ----------------------------------------------------------
  
  predictor_vars <- c(setdiff(vars_num, "NASC"), "fod")
  
  x_train <- ds_train[, predictor_vars, drop = FALSE]
  
  x_test <- ds_test[, predictor_vars, drop = FALSE]
  
  y_train <- ds_train$NASC
  
  # ----------------------------------------------------------
  # Random Forest
  # ----------------------------------------------------------
  
  model <- randomForest(
    x = x_train,
    y = y_train,
    ntree = ntree,
    mtry = mtry,
    nodesize = nodesize,
    importance = FALSE
  )
  
  # ----------------------------------------------------------
  # Prediction
  # ----------------------------------------------------------
  
  prediction <- predict(model, newdata = x_test)
  
  observed <- ds_test$NASC
  
  calculate_metrics(observed, prediction)
}


# ============================================================
# 06 - XGBOOST
# ============================================================

fit_xgb <- function(train, test, vars_num, max_depth, eta, min_child_weight, subsample, colsample_bytree, nrounds) {
  
  # ----------------------------------------------------------
  # Scaling des variables numériques
  # ----------------------------------------------------------
  
  train_scaled <- scale(train[, vars_num])
  
  center <- attr(train_scaled, "scaled:center")
  
  scale_values <- attr(train_scaled, "scaled:scale")
  
  test_scaled <- scale(test[, vars_num], center = center, scale = scale_values)
  
  ds_train <- train
  
  ds_test <- test
  
  ds_train[, vars_num] <- train_scaled
  
  ds_test[, vars_num] <- test_scaled
  
  # ----------------------------------------------------------
  # Variables numériques
  # ----------------------------------------------------------
  
  predictor_vars <- setdiff(vars_num, "NASC")
  
  X_train_num <- as.matrix(ds_train[, predictor_vars, drop = FALSE])
  
  X_test_num <- as.matrix(ds_test[, predictor_vars, drop = FALSE])
  
  # ----------------------------------------------------------
  # FOD en variables dummy
  # ----------------------------------------------------------
  
  fod_levels <- union(levels(ds_train$fod), levels(ds_test$fod))
  
  ds_train$fod <- factor(ds_train$fod, levels = fod_levels)
  
  ds_test$fod <- factor(ds_test$fod, levels = fod_levels)
  
  fod_train <- model.matrix(~ fod - 1, data = ds_train)
  
  fod_test <- model.matrix(~ fod - 1, data = ds_test)
  
  common_cols <- union(colnames(fod_train), colnames(fod_test))
  
  fod_train_complete <- matrix(0, nrow = nrow(fod_train), ncol = length(common_cols))
  
  fod_test_complete <- matrix(0, nrow = nrow(fod_test), ncol = length(common_cols))
  
  colnames(fod_train_complete) <- common_cols
  
  colnames(fod_test_complete) <- common_cols
  
  fod_train_complete[, colnames(fod_train)] <- fod_train
  
  fod_test_complete[, colnames(fod_test)] <- fod_test
  
  # ----------------------------------------------------------
  # Matrices finales
  # ----------------------------------------------------------
  
  X_train <- cbind(X_train_num, fod_train_complete)
  
  X_test <- cbind(X_test_num, fod_test_complete)
  
  y_train <- ds_train$NASC
  
  # ----------------------------------------------------------
  # XGBoost
  # ----------------------------------------------------------
  
  dtrain <- xgb.DMatrix(data = X_train, label = y_train)
  
  dtest <- xgb.DMatrix(data = X_test)
  
  model <- xgb.train(
    params = list(
      objective = "reg:squarederror",
      eval_metric = "rmse",
      max_depth = max_depth,
      eta = eta,
      min_child_weight = min_child_weight,
      subsample = subsample,
      colsample_bytree = colsample_bytree
    ),
    data = dtrain,
    nrounds = nrounds,
    verbose = 0
  )
  
  # ----------------------------------------------------------
  # Prediction
  # ----------------------------------------------------------
  
  prediction <- predict(model, dtest)
  
  observed <- ds_test$NASC
  
  calculate_metrics(observed, prediction)
}


# ============================================================
# 07 - RESULTATS
# ============================================================

results <- data.frame()

combination_id <- 0


# ============================================================
# 08 - BOUCLE PRINCIPALE
# ============================================================

for (FOD_0 in FOD_values) {
  
  for (log in log_values) {
    
    cat("\n\n")
    cat("####################################################\n")
    cat("FOD_0 =", FOD_0, " | LOG =", log, "\n")
    cat("####################################################\n")
    
    # --------------------------------------------------------
    # Préparation
    # --------------------------------------------------------
    
    prepared <- prepare_data(
      datas = datas,
      pigment_type = pigment_type,
      FOD_0 = FOD_0,
      log = log,
      diurnal_period = diurnal_period
    )
    
    df_raw <- prepared$df
    
    vars_num <- prepared$vars_num
    
    # Variables numériques utilisées pour le leakage
    predictor_vars <- setdiff(vars_num, "NASC")
    
    cat("Observations avant leakage:", nrow(df_raw), "\n")
    
    # ========================================================
    # DIFFERENTS SEUILS DE LEAKAGE
    # ========================================================
    
    for (dist_thr in distance_thresholds) {
      
      cat("\n")
      cat("-----------------------------------------------\n")
      cat("Distance threshold =", dist_thr, "\n")
      cat("-----------------------------------------------\n")
      
      # ------------------------------------------------------
      # Suppression leakage
      # ------------------------------------------------------
      
      leakage <- remove_close_pairs(
        df_input = df_raw,
        predictor_vars = predictor_vars,
        dist_thr = dist_thr,
        time_thr = time_thr
      )
      
      df <- leakage$df
      
      cat("Observations supprimées:", leakage$n_removed, "\n")
      cat("Observations restantes:", leakage$n_remaining, "\n")
      cat("Pourcentage supprimé:", round(leakage$percent_removed, 1), "%\n")
      
      
      # ======================================================
      # GRID RANDOM FOREST
      # ======================================================
      
      p <- length(c(predictor_vars, "fod"))
      
      mtry_values <- unique(
        pmin(
          p,
          pmax(
            1,
            round(p * mtry_fraction)
          )
        )
      )
      
      grid_rf <- expand.grid(
        mtry = mtry_values,
        nodesize = nodesize_values,
        ntree = ntree_values
      )
      
      cat("Combinaisons RF:", nrow(grid_rf), "\n")
      
      
      # ======================================================
      # LOYO - RANDOM FOREST
      # ======================================================
      
      for (i in seq_len(nrow(grid_rf))) {
        
        hp <- grid_rf[i, ]
        
        cat(
          "\n[PROGRESSION] RF LOYO : combinaison",
          i,
          "/",
          nrow(grid_rf),
          "| mtry =", hp$mtry,
          "| nodesize =", hp$nodesize,
          "| ntree =", hp$ntree,
          "\n"
        )
        
        years <- sort(unique(df$year))
        
        loyo_metrics <- data.frame()
        
        for (yr in years) {
          
          cat("   -> Année test :", yr, "\n")
          
          train <- df[df$year != yr, ]
          
          test <- df[df$year == yr, ]
          
          if (nrow(train) < 10 || nrow(test) < 2) {
            next
          }
          
          metrics <- fit_rf(
            train = train,
            test = test,
            vars_num = vars_num,
            mtry = hp$mtry,
            nodesize = hp$nodesize,
            ntree = hp$ntree
          )
          
          metrics$year <- yr
          
          loyo_metrics <- rbind(loyo_metrics, metrics)
        }
        
        if (nrow(loyo_metrics) > 0) {
          
          combination_id <- combination_id + 1
          
          results <- rbind(
            results,
            data.frame(
              combination = combination_id,
              model = "RF",
              validation = "LOYO",
              FOD_0 = FOD_0,
              log = log,
              dist_thr = dist_thr,
              n_initial = leakage$n_initial,
              n_removed = leakage$n_removed,
              percent_removed = leakage$percent_removed,
              n_remaining = leakage$n_remaining,
              mtry = hp$mtry,
              nodesize = hp$nodesize,
              ntree = hp$ntree,
              max_depth = NA,
              eta = NA,
              min_child_weight = NA,
              subsample = NA,
              colsample_bytree = NA,
              nrounds = NA,
              RMSE = mean(loyo_metrics$RMSE),
              RMSE_sd = sd(loyo_metrics$RMSE),
              R2 = mean(loyo_metrics$R2),
              R2_sd = sd(loyo_metrics$R2),
              MAE = mean(loyo_metrics$MAE),
              MAE_sd = sd(loyo_metrics$MAE)
            )
          )
        }
      }
      
      
      # ======================================================
      # RANDOM SPLIT - RANDOM FOREST
      # ======================================================
      
      for (i in seq_len(nrow(grid_rf))) {
        
        hp <- grid_rf[i, ]
        
        cat(
          "\n[PROGRESSION] RF RANDOM SPLIT : combinaison",
          i,
          "/",
          nrow(grid_rf),
          "| mtry =", hp$mtry,
          "| nodesize =", hp$nodesize,
          "| ntree =", hp$ntree,
          "\n"
        )
        
        rs_metrics <- data.frame()
        
        for (rep in seq_len(n_repeats)) {
          
          cat(
            "   -> Répétition",
            rep,
            "/",
            n_repeats,
            "\n"
          )
          
          set.seed(1000 + rep)
          
          n_train <- floor(train_fraction * nrow(df))
          
          train_index <- sample(seq_len(nrow(df)), size = n_train)
          
          test_index <- setdiff(seq_len(nrow(df)), train_index)
          
          train <- df[train_index, ]
          
          test <- df[test_index, ]
          
          metrics <- fit_rf(
            train = train,
            test = test,
            vars_num = vars_num,
            mtry = hp$mtry,
            nodesize = hp$nodesize,
            ntree = hp$ntree
          )
          
          rs_metrics <- rbind(rs_metrics, metrics)
        }
        
        combination_id <- combination_id + 1
        
        results <- rbind(
          results,
          data.frame(
            combination = combination_id,
            model = "RF",
            validation = "RS",
            FOD_0 = FOD_0,
            log = log,
            dist_thr = dist_thr,
            n_initial = leakage$n_initial,
            n_removed = leakage$n_removed,
            percent_removed = leakage$percent_removed,
            n_remaining = leakage$n_remaining,
            mtry = hp$mtry,
            nodesize = hp$nodesize,
            ntree = hp$ntree,
            max_depth = NA,
            eta = NA,
            min_child_weight = NA,
            subsample = NA,
            colsample_bytree = NA,
            nrounds = NA,
            RMSE = mean(rs_metrics$RMSE),
            RMSE_sd = sd(rs_metrics$RMSE),
            R2 = mean(rs_metrics$R2),
            R2_sd = sd(rs_metrics$R2),
            MAE = mean(rs_metrics$MAE),
            MAE_sd = sd(rs_metrics$MAE)
          )
        )
      }
      
      
      # ======================================================
      # GRID XGBOOST
      # ======================================================
      
      grid_xgb <- expand.grid(
        max_depth = max_depth_values,
        eta = eta_values,
        min_child_weight = min_child_weight_values,
        subsample = subsample_values,
        colsample_bytree = colsample_bytree_values
      )
      
      cat("Combinaisons XGBoost:", nrow(grid_xgb), "\n")
      
      
      # ======================================================
      # LOYO - XGBOOST
      # ======================================================
      
      for (i in seq_len(nrow(grid_xgb))) {
        
        hp <- grid_xgb[i, ]
        
        cat(
          "\n[PROGRESSION] XGBOOST LOYO : combinaison",
          i,
          "/",
          nrow(grid_xgb),
          "| max_depth =", hp$max_depth,
          "| eta =", hp$eta,
          "| min_child_weight =", hp$min_child_weight,
          "| subsample =", hp$subsample,
          "| colsample =", hp$colsample_bytree,
          "\n"
        )
        
        years <- sort(unique(df$year))
        
        loyo_metrics <- data.frame()
        
        for (yr in years) {
          
          cat("   -> Année test :", yr, "\n")
          
          train <- df[df$year != yr, ]
          
          test <- df[df$year == yr, ]
          
          if (nrow(train) < 10 || nrow(test) < 2) {
            next
          }
          
          metrics <- fit_xgb(
            train = train,
            test = test,
            vars_num = vars_num,
            max_depth = hp$max_depth,
            eta = hp$eta,
            min_child_weight = hp$min_child_weight,
            subsample = hp$subsample,
            colsample_bytree = hp$colsample_bytree,
            nrounds = nrounds
          )
          
          metrics$year <- yr
          
          loyo_metrics <- rbind(loyo_metrics, metrics)
        }
        
        if (nrow(loyo_metrics) > 0) {
          
          combination_id <- combination_id + 1
          
          results <- rbind(
            results,
            data.frame(
              combination = combination_id,
              model = "XGBoost",
              validation = "LOYO",
              FOD_0 = FOD_0,
              log = log,
              dist_thr = dist_thr,
              n_initial = leakage$n_initial,
              n_removed = leakage$n_removed,
              percent_removed = leakage$percent_removed,
              n_remaining = leakage$n_remaining,
              mtry = NA,
              nodesize = NA,
              ntree = NA,
              max_depth = hp$max_depth,
              eta = hp$eta,
              min_child_weight = hp$min_child_weight,
              subsample = hp$subsample,
              colsample_bytree = hp$colsample_bytree,
              nrounds = nrounds,
              RMSE = mean(loyo_metrics$RMSE),
              RMSE_sd = sd(loyo_metrics$RMSE),
              R2 = mean(loyo_metrics$R2),
              R2_sd = sd(loyo_metrics$R2),
              MAE = mean(loyo_metrics$MAE),
              MAE_sd = sd(loyo_metrics$MAE)
            )
          )
        }
      }
      
      
      # ======================================================
      # RANDOM SPLIT - XGBOOST
      # ======================================================
      
      for (i in seq_len(nrow(grid_xgb))) {
        
        hp <- grid_xgb[i, ]
        
        cat(
          "\n[PROGRESSION] XGBOOST RANDOM SPLIT : combinaison",
          i,
          "/",
          nrow(grid_xgb),
          "| max_depth =", hp$max_depth,
          "| eta =", hp$eta,
          "| min_child_weight =", hp$min_child_weight,
          "| subsample =", hp$subsample,
          "| colsample =", hp$colsample_bytree,
          "\n"
        )
        
        rs_metrics <- data.frame()
        
        for (rep in seq_len(n_repeats)) {
          
          cat(
            "   -> Répétition",
            rep,
            "/",
            n_repeats,
            "\n"
          )
          
          set.seed(1000 + rep)
          
          n_train <- floor(train_fraction * nrow(df))
          
          train_index <- sample(seq_len(nrow(df)), size = n_train)
          
          test_index <- setdiff(seq_len(nrow(df)), train_index)
          
          train <- df[train_index, ]
          
          test <- df[test_index, ]
          
          metrics <- fit_xgb(
            train = train,
            test = test,
            vars_num = vars_num,
            max_depth = hp$max_depth,
            eta = hp$eta,
            min_child_weight = hp$min_child_weight,
            subsample = hp$subsample,
            colsample_bytree = hp$colsample_bytree,
            nrounds = nrounds
          )
          
          rs_metrics <- rbind(rs_metrics, metrics)
        }
        
        combination_id <- combination_id + 1
        
        results <- rbind(
          results,
          data.frame(
            combination = combination_id,
            model = "XGBoost",
            validation = "RS",
            FOD_0 = FOD_0,
            log = log,
            dist_thr = dist_thr,
            n_initial = leakage$n_initial,
            n_removed = leakage$n_removed,
            percent_removed = leakage$percent_removed,
            n_remaining = leakage$n_remaining,
            mtry = NA,
            nodesize = NA,
            ntree = NA,
            max_depth = hp$max_depth,
            eta = hp$eta,
            min_child_weight = hp$min_child_weight,
            subsample = hp$subsample,
            colsample_bytree = hp$colsample_bytree,
            nrounds = nrounds,
            RMSE = mean(rs_metrics$RMSE),
            RMSE_sd = sd(rs_metrics$RMSE),
            R2 = mean(rs_metrics$R2),
            R2_sd = sd(rs_metrics$R2),
            MAE = mean(rs_metrics$MAE),
            MAE_sd = sd(rs_metrics$MAE)
          )
        )
      }
    }
  }
}


# ============================================================
# 09 - MEILLEUR MODELE GLOBAL LOYO
# ============================================================

best_LOYO <- results |>
  filter(validation == "LOYO") |>
  arrange(desc(R2)) |>
  slice(1)

cat("\n\n")
cat("====================================================\n")
cat("MEILLEUR MODELE GLOBAL - LOYO\n")
cat("====================================================\n")

print(best_LOYO)


# ============================================================
# 10 - MEILLEUR MODELE GLOBAL RANDOM SPLIT
# ============================================================

best_RS <- results |>
  filter(validation == "RS") |>
  arrange(desc(R2)) |>
  slice(1)

cat("\n\n")
cat("====================================================\n")
cat("MEILLEUR MODELE GLOBAL - RANDOM SPLIT\n")
cat("====================================================\n")

print(best_RS)


# ============================================================
# 11 - MEILLEUR RF LOYO
# ============================================================

best_RF_LOYO <- results |>
  filter(model == "RF", validation == "LOYO", dist_thr == 1) |>
  arrange(desc(R2)) |>
  slice(1)

cat("\n\n")
cat("====================================================\n")
cat("MEILLEUR RANDOM FOREST - LOYO\n")
cat("====================================================\n")

print(best_RF_LOYO)


# ============================================================
# 12 - MEILLEUR XGBOOST LOYO
# ============================================================

best_XGB_LOYO <- results |>
  filter(model == "XGBoost", validation == "LOYO", dist_thr == 1) |>
  arrange(desc(R2)) |>
  slice(1)

cat("\n\n")
cat("====================================================\n")
cat("MEILLEUR XGBOOST - LOYO\n")
cat("====================================================\n")

print(best_XGB_LOYO)


# ============================================================
# 13 - COMPARAISON RF vs XGBOOST
# ============================================================

model_comparison <- results |>
  group_by(model, validation) |>
  summarise(
    best_RMSE = min(RMSE),
    best_R2 = R2[which.min(RMSE)],
    best_MAE = MAE[which.min(RMSE)],
    .groups = "drop"
  ) |>
  arrange(validation, best_RMSE)

cat("\n\n")
cat("====================================================\n")
cat("COMPARAISON RF vs XGBOOST\n")
cat("====================================================\n")

print(model_comparison)


# ============================================================
# 14 - MEILLEUR MODELE PAR CONDITION
# ============================================================

best_by_condition <- results |>
  group_by(validation, FOD_0, log, dist_thr) |>
  arrange(RMSE) |>
  slice(1) |>
  ungroup() |>
  arrange(validation, RMSE)

cat("\n\n")
cat("====================================================\n")
cat("MEILLEUR MODELE PAR CONDITION\n")
cat("====================================================\n")

print(best_by_condition)


# ============================================================
# 15 - TOP 10 LOYO
# ============================================================

cat("\n\n")
cat("====================================================\n")
cat("TOP 10 LOYO - RF + XGBOOST\n")
cat("====================================================\n")

results |>
  filter(validation == "LOYO") |>
  arrange(desc(R2)) |>
  select(
    model,
    FOD_0,
    log,
    dist_thr,
    n_removed,
    n_remaining,
    mtry,
    nodesize,
    ntree,
    max_depth,
    eta,
    min_child_weight,
    subsample,
    colsample_bytree,
    RMSE,
    RMSE_sd,
    R2,
    MAE
  ) |>
  slice_head(n = 10) |>
  print()


# ============================================================
# 16 - TOP 10 RANDOM SPLIT
# ============================================================

cat("\n\n")
cat("====================================================\n")
cat("TOP 10 RANDOM SPLIT - RF + XGBOOST\n")
cat("====================================================\n")

results |>
  filter(validation == "RS") |>
  arrange(desc(R2)) |>
  select(
    model,
    FOD_0,
    log,
    dist_thr,
    n_removed,
    n_remaining,
    mtry,
    nodesize,
    ntree,
    max_depth,
    eta,
    min_child_weight,
    subsample,
    colsample_bytree,
    RMSE,
    RMSE_sd,
    R2,
    MAE
  ) |>
  slice_head(n = 10) |>
  print()


# ============================================================
# 17 - RESUME DES CONDITIONS
# ============================================================

summary_conditions <- results |>
  group_by(model, validation, FOD_0, log, dist_thr) |>
  summarise(
    RMSE_mean = mean(RMSE),
    RMSE_sd = mean(RMSE_sd),
    R2_mean = mean(R2),
    MAE_mean = mean(MAE),
    n_removed = first(n_removed),
    n_remaining = first(n_remaining),
    .groups = "drop"
  ) |>
  arrange(validation, RMSE_mean)

cat("\n\n")
cat("====================================================\n")
cat("RESUME DES CONDITIONS\n")
cat("====================================================\n")

print(summary_conditions, n=500)


# ============================================================
# 18 - COMPARAISON RF / XGBOOST PAR DISTANCE
# ============================================================

ggplot(
  summary_conditions,
  aes(x = dist_thr, y = RMSE_mean, color = model)
) +
  geom_line() +
  geom_point(size = 3) +
  facet_grid(
    FOD_0 ~ log + validation,
    labeller = labeller(
      FOD_0 = function(x) paste("FOD_0 =", x),
      log = function(x) paste("log =", x),
      validation = function(x) paste("Validation =", x)
    )
  ) +
  theme_bw() +
  labs(
    x = "Distance threshold",
    y = "RMSE moyen",
    color = "Modèle",
    title = "Comparaison Random Forest vs XGBoost"
  )


# ============================================================
# 19 - COMPARAISON DES MEILLEURS MODELES
# ============================================================

best_models_plot <- results |>
  group_by(model, validation) |>
  slice_min(RMSE, n = 1, with_ties = FALSE) |>
  ungroup()

ggplot(
  best_models_plot,
  aes(x = model, y = RMSE, fill = model)
) +
  geom_col() +
  facet_wrap(~ validation) +
  theme_bw() +
  labs(
    x = "Modèle",
    y = "Meilleur RMSE",
    title = "Meilleur RMSE obtenu par modèle"
  )


# ============================================================
# 20 - SAUVEGARDE
# ============================================================

outpath <- "F:/data_elise/tune_trees/"

dir.create(
  outpath,
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  results,
  paste0(outpath, "RF_vs_XGBoost_tuning_results.csv"),
  row.names = FALSE
)

write.csv(
  best_by_condition,
  paste0(outpath, "RF_vs_XGBoost_best_by_condition.csv"),
  row.names = FALSE
)

write.csv(
  model_comparison,
  paste0(outpath, "RF_vs_XGBoost_model_comparison.csv"),
  row.names = FALSE
)

write.csv(
  summary_conditions,
  paste0(outpath, "RF_vs_XGBoost_summary_conditions.csv"),
  row.names = FALSE
)

cat("\n\n")
cat("====================================================\n")
cat("SAUVEGARDE TERMINEE\n")
cat("====================================================\n")
cat("Dossier :", outpath, "\n")
cat("Fichier 1 : RF_vs_XGBoost_tuning_results.csv\n")
cat("Fichier 2 : RF_vs_XGBoost_best_by_condition.csv\n")
cat("Fichier 3 : RF_vs_XGBoost_model_comparison.csv\n")
cat("Fichier 4 : RF_vs_XGBoost_summary_conditions.csv\n")