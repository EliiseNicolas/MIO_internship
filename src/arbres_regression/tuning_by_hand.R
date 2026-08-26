library(xgboost)
library(dplyr)
library(ggplot2)
library(sf)
library(geosphere)
library(tidyr)
library(patchwork)
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
log <- TRUE
# CART
max_depth <- 7
min_leafs <- 
# Random Forest
ntree_values <- c(200, 500)
nodesize_values <- c(3, 5, 10)
mtry_fraction <- c(0.3, 0.5, 0.8, 1.0)

# XGBoost
max_depth_values <- c(2, 3, 4, 5)
eta_values <- c(0.03, 0.05, 0.1)
min_child_weight_values <- c(3, 5, 10)
subsample_values <- c(0.8, 1)
colsample_bytree_values <- c(0.8, 1)

nrounds <- 1000


path_datas <- paste0("F:/data_elise/ds_NASC_pig_ftle_fod/ds_NASC_mean_pig_grid_all/ds_NASC_mean_pig_grid_pig_ftle_fod_2018_2021_2023_transect_", freq, "kHz_mask9.rds")

datas <- readRDS(path_datas)
str(datas)

# ============================================================
# 02 - PREPARATION DES DONNEES
# ============================================================

prepare_data <- function(datas, pigment_type = "chla_ratio", FOD_0 = FALSE, log = TRUE, diurnal_period = 3) {
  
  dat <- datas
  
  # Day / night
  dat <- dat[dat$day == diurnal_period, ]
  
  # NASC extreme values
  q <- quantile(dat$nasc, probs = c(0.05, 0.95), na.rm = TRUE)
  
  dat <- dat |>
    dplyr::filter(
      .data$nasc >= q[1],
      .data$nasc <= q[2]
    )
  
  # Log NASC

  if (log) {
    if (any(dat$nasc <= 0, na.rm = TRUE)) {
      stop("ERREUR : NASC contient des valeurs <= 0. Impossible d'appliquer log10().")
    }
    dat$nasc <- log10(dat$nasc)
  }

  # FOD_0
  if (FOD_0) {
    fod_0 <- as.character(dat$fod)
    fod_num <- suppressWarnings(as.numeric(fod_0))
    fod_0[!is.na(fod_num) & fod_num > 6] <- "0"
    dat$fod <- fod_0
  }

  # Base dataframe

  df <- data.frame(
    NASC = dat$nasc,
    year = format(dat$time_nasc, "%Y"),
    time = dat$time_nasc,
    lat <- dat$lat_nasc,
    lon <- dat$lon_nasc,
    fod = dat$fod,
    ftle = dat$ftle
  )

  # Pigments

  if (pigment_type == "chla_ratio") {
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

  # FOD factor

  df$fod <- as.factor(df$fod)
  df$fod[df$fod == "NA"] <- NA
  df$fod <- droplevels(df$fod)

  # Missing values

  df <- df[complete.cases(df), ] # retirer toutes les lignes ou il y a un NA dans la colonne
  print(colSums(is.na(df)))
  
  return(df)
}
df <- prepare_data(datas)
str(df)

# ============================================================
# 03 - DISTRIBUTION DES NOUVELLES DONNEES
# ============================================================
distributions_datas <- function(ds){
  # distrib pigs
  pig_vars <- c(
    "total_chla",
    "per_ratio_chla",
    "but_ratio_chla",
    "fuco_ratio_chla",
    "hex_ratio_chla",
    "allo_ratio_chla",
    "zea_ratio_chla",
    "chlb_ratio_chla"
  )
  
  pig_long <- ds %>%
    select(all_of(pig_vars)) %>%
    pivot_longer(
      cols = everything(),
      names_to = "pigment",
      values_to = "value"
    ) %>%
    filter(
      is.finite(value)
    )
  
  p1 <- ggplot(
    pig_long,
    aes(x = value)
  ) +
    geom_histogram(
      bins = 50
    ) +
    facet_wrap(
      ~ pigment,
      scales = "free"
    ) +
    theme_bw() +
    scale_x_continuous(
      n.breaks = 3
    ) +
    labs(
      x = "Concentration",
      y = "Count"
    )

  n_total <- nrow(ds)
  n_complete <- sum(complete.cases(ds[pig_vars]))
  n_missing <- n_total - n_complete
  pct_missing <- n_missing / n_total * 100
  pct_complete <- n_complete / n_total * 100
  

  p2 <- ggplot() +
    annotate(
      "text",
      x = 0,
      y = 0,
      label = paste0(
        "Total rows :\n",
        n_total,
        
        "\n\nRows with missing values :\n",
        n_missing,
        " (",
        round(pct_missing, 2),
        "%)",
        
        "\n\nRows with no missing values :\n",
        n_complete,
        " (",
        round(pct_complete, 2),
        "%)"
      ),
      size = 4
    ) +
    theme_void()
  
  p <- p1 + p2 +
    plot_layout(
      widths = c(3, 1)
    ) +
    plot_annotation(
      title = paste("Distribution of pigment variables - 2018, 2021, 2023 - Transect -", dp, "-", freq, "kHz"),
      subtitle = "Missing values are excluded from histograms"
    )
  print(p) 
  
  # FOD
  fod_distribution <- df %>%
    count(fod) %>%
    mutate(
      percentage = n / sum(n) * 100
    )
  
  print(fod_distribution)
  
  
  ggplot(
    fod_distribution,
    aes(
      x = fod,
      y = percentage
    )
  ) +
    geom_col() +
    labs(
      x = "FOD cluster",
      y = "Percentage (%)",
      title = paste("FOD distribution - 2018, 2021, 2023 - Transect -", dp, "-", freq, "kHz")
    ) +
    theme_bw()
  
  # NASC
  if (log) {
    nasc_scale <- "log10(NASC)"
    x_label <- "log10(NASC)"
    p <- ggplot(df, aes(x = NASC,
                        y = after_stat(count / sum(count) * 100))) +
      geom_histogram(bins = 200)
    
  } else {
    nasc_scale <- "NASC"
    x_label <- "NASC"
    
    p <- ggplot(df, aes(x = NASC,
                        y = after_stat(count / sum(count) * 100))) +
      geom_histogram(bins = 200)
  }
  
  p <- p +
    theme_bw() +
    labs(
      x = x_label,
      y = "Percentage (%)",
      title = paste(
        nasc_scale,
        "distribution - 2018, 2021, 2023 -",
        dp, "-", freq, "kHz"
      ),
      subtitle = paste(
        nasc_scale,
        "without extreme values (0.05 - 0.95)"
      )
    )
  
  print(p)
  
  # ftle
  ggplot(
    df,
    aes(
      x = ftle,
      y = after_stat(count / sum(count) * 100)
    )
  ) +
    geom_histogram(
      bins = 200
    ) +
    theme_bw() +
    labs(
      x = "FTLE",
      y = "Percentage (%)",
      title = paste("FTLE distribution - 2018, 2021, 2023 - Transect -", dp, "-", freq, "kHz")
    )
  
  # distrib spatiale et temporelle des données
  df2 <- df %>%
    mutate(
      year = as.integer(year),
      lat = lat....dat.lat_nasc,
      lon = lon....dat.lon_nasc
    ) %>%
    arrange(year, time)
  
  df_dist <- df2 %>%
    group_by(year) %>%
    arrange(time, .by_group = TRUE) %>%
    mutate(
      time_diff_h = as.numeric(difftime(time, lag(time), units = "hours")),
      
      distance_km = geosphere::distHaversine(
        cbind(lag(lon), lag(lat)),
        cbind(lon, lat)
      ) / 1000
    ) %>%
    ungroup()
  
  head(df_dist[, c(
    "year", "time", "lat", "lon",
    "time_diff_h", "distance_km"
  )])
  
  distance_moyenne <- df_dist %>%
    group_by(year) %>%
    summarise(
      n = n(),
      distance_moyenne_km = mean(distance_km, na.rm = TRUE),
      distance_mediane_km = median(distance_km, na.rm = TRUE),
      temps_moyen_h = mean(time_diff_h, na.rm = TRUE),
      temps_median_h = median(time_diff_h, na.rm = TRUE)
    )
  
  distance_moyenne
  
  library(viridis)
  
  df2 <- df2 %>%
    mutate(
      day = as.Date(time)
    )
  
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(viridis)
  
  
  # Fonction pour créer un plot par année
  make_map <- function(data, year_value, palette_option) {
    
    d <- data %>%
      filter(year == year_value) %>%
      arrange(day) %>%
      mutate(
        day = factor(day, levels = unique(day))
      )
    
    ggplot(d, aes(x = lon, y = lat)) +
      geom_point(
        aes(color = day),
        size = 3,
        alpha = 0.8
      ) +
      coord_equal() +
      scale_color_viridis_d(option = palette_option) +
      theme_bw() +
      labs(
        x = "Longitude",
        y = "Latitude",
        color = "Jour",
        title = as.character(year_value)
      ) +
      theme(
        legend.position = "right"
      )
  }
  
  p2018 <- make_map(df2, 2018, "D")
  p2021 <- make_map(df2, 2021, "C")
  p2023 <- make_map(df2, 2023, "A")
  
  p_map <- p2018 + p2021 + p2023 +
    plot_layout(ncol = 1)
  
  print(p_map)

    
  df_dist2 <- df_dist %>%
    filter(!is.na(distance_km)) %>%
    group_by(year) %>%
    mutate(
      n_total = n()
    ) %>%
    ungroup()
  
  p <- ggplot(
    df_dist2 %>% filter(distance_km <= 100),
    aes(
      x = distance_km,
      y = after_stat(count / unique(n_total) * 100)
    )
  ) +
    geom_histogram(
      bins = 5,
      boundary = 0
    ) +
    facet_wrap(~year) +
    scale_x_continuous(
      limits = c(0, 100)
    ) +
    theme_bw() +
    labs(
      x = "Distance entre observations successives (km)",
      y = "Pourcentage de toutes les observations (%)",
      title = "Distribution des distances spatiales",
      subtitle = "Distances ≤ 100 km"
    )
  
  print(p)
  
  # Intervalle temporel entre observations successives
  df_dist2 <- df_dist %>%
    filter(!is.na(time_diff_h)) %>%
    group_by(year) %>%
    mutate(
      time_diff_min = time_diff_h * 60,
      n_total = n()
    ) %>%
    ungroup()
  
  p <- ggplot(
    df_dist2 %>% 
      filter(time_diff_min <= 120),
    aes(
      x = time_diff_min,
      y = after_stat(count / unique(n_total) * 100)
    )
  ) +
    geom_histogram(
      bins = 12,
      boundary = 0
    ) +
    facet_wrap(~year) +
    scale_x_continuous(
      limits = c(0, 120),
      breaks = seq(0, 120, 20)
    ) +
    theme_bw() +
    labs(
      x = "Intervalle entre observations successives (minutes)",
      y = "Pourcentage de toutes les observations (%)",
      title = "Distribution des intervalles temporels",
      subtitle = "Intervalles ≤ 120 minutes"
    )
  
  print(p)
  
}
distributions_datas(df)

################ je me suis arrêtée là
# ============================================================
# 04 - SUPPRESSION DU DATA LEAKAGE
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
    
    # Distance calculée UNIQUEMENT sur les prédicteurs
    
    df_scaled <- scale(df_work[, predictor_vars, drop = FALSE])
    dist_mat <- as.matrix(dist(df_scaled))
    
    # Différence temporelle
    
    time_diff <- abs(outer(as.numeric(df_work$time), as.numeric(df_work$time), FUN = "-"))
    
    # Paires proches
    
    pairs <- which(time_diff <= time_thr & dist_mat < dist_thr, arr.ind = TRUE)
    
    pairs <- pairs[pairs[, 1] < pairs[, 2], , drop = FALSE]
    
    n_pairs <- nrow(pairs)
    
    if (n_pairs == 0) {
      break
    }
    
    # Suppression du deuxième élément de la première paire
    
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
str(df)
predictor_vars <- setdiff(names(df),
  c("time", "lat", "lon", "year", "NASC", "fod")
)
print(predictor_vars)
df_cleaned <- remove_close_pairs(df, predictor_vars = predictor_vars)
distributions_datas(df_cleaned$df)
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
# 05 - XGBOOST
# ============================================================

fit_xgb <- function(train, test, vars_num, max_depth, eta, min_child_weight, subsample, colsample_bytree, nrounds) {
  
  
  # ----------------------------------------------------------
  # Matrices X / y
  # ----------------------------------------------------------
  
  predictor_vars <- setdiff(vars_num, "NASC")
  
  X_train <- as.matrix(ds_train[, predictor_vars, drop = FALSE])
  X_test <- as.matrix(ds_test[, predictor_vars, drop = FALSE])
  
  y_train <- ds_train$NASC
  
  # ----------------------------------------------------------
  # FOD
  # ----------------------------------------------------------
  
  if ("fod" %in% names(ds_train)) {
    
    fod_train <- model.matrix(~ fod - 1, data = ds_train)
    fod_test <- model.matrix(~ fod - 1, data = ds_test)
    
    common_cols <- intersect(colnames(fod_train), colnames(fod_test))
    
    fod_train <- fod_train[, common_cols, drop = FALSE]
    fod_test <- fod_test[, common_cols, drop = FALSE]
    
    X_train <- cbind(X_train, fod_train)
    X_test <- cbind(X_test, fod_test)
  }
  
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
# 06 - RESULTATS
# ============================================================

results <- data.frame()
combination_id <- 0


# ============================================================
# 07 - BOUCLE FOD / LOG / LEAKAGE
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
    
    prepared <- prepare_data(datas = datas, pigment_type = pigment_type, FOD_0 = FOD_0, log = log, diurnal_period = diurnal_period)
    
    df_raw <- prepared$df
    vars_num <- prepared$vars_num
    
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
      
      leakage <- remove_close_pairs(df_input = df_raw, predictor_vars = predictor_vars, dist_thr = dist_thr, time_thr = time_thr)
      
      df <- leakage$df
      
      cat("Observations supprimées:", leakage$n_removed, "\n")
      cat("Observations restantes:", leakage$n_remaining, "\n")
      cat("Pourcentage supprimé:", round(leakage$percent_removed, 1), "%\n")
      
      # ------------------------------------------------------
      # Grid XGBoost
      # ------------------------------------------------------
      
      grid <- expand.grid(
        max_depth = max_depth_values,
        eta = eta_values,
        min_child_weight = min_child_weight_values,
        subsample = subsample_values,
        colsample_bytree = colsample_bytree_values
      )
      
      cat("Nombre de combinaisons:", nrow(grid), "\n")
      
      # ======================================================
      # LOYO
      # ======================================================
      
      for (i in seq_len(nrow(grid))) {
        
        hp <- grid[i, ]
        
        years <- sort(unique(df$year))
        loyo_metrics <- data.frame()
        
        for (yr in years) {
          
          train <- df[df$year != yr, ]
          test <- df[df$year == yr, ]
          
          if (nrow(train) < 10 || nrow(test) < 2) {
            next
          }
          
          metrics <- fit_xgb(train = train, test = test, vars_num = vars_num, max_depth = hp$max_depth, eta = hp$eta, min_child_weight = hp$min_child_weight, subsample = hp$subsample, colsample_bytree = hp$colsample_bytree, nrounds = nrounds)
          
          metrics$year <- yr
          
          loyo_metrics <- rbind(loyo_metrics, metrics)
        }
        
        if (nrow(loyo_metrics) > 0) {
          
          combination_id <- combination_id + 1
          
          results <- rbind(
            results,
            data.frame(
              combination = combination_id,
              validation = "LOYO",
              FOD_0 = FOD_0,
              log = log,
              dist_thr = dist_thr,
              n_initial = leakage$n_initial,
              n_removed = leakage$n_removed,
              percent_removed = leakage$percent_removed,
              n_remaining = leakage$n_remaining,
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
      # RANDOM SPLIT 80/20
      # ======================================================
      
      for (i in seq_len(nrow(grid))) {
        
        hp <- grid[i, ]
        rs_metrics <- data.frame()
        
        for (rep in seq_len(n_repeats)) {
          
          set.seed(1000 + rep)
          
          n_train <- floor(train_fraction * nrow(df))
          
          train_index <- sample(seq_len(nrow(df)), size = n_train)
          test_index <- setdiff(seq_len(nrow(df)), train_index)
          
          train <- df[train_index, ]
          test <- df[test_index, ]
          
          metrics <- fit_xgb(train = train, test = test, vars_num = vars_num, max_depth = hp$max_depth, eta = hp$eta, min_child_weight = hp$min_child_weight, subsample = hp$subsample, colsample_bytree = hp$colsample_bytree, nrounds = nrounds)
          
          rs_metrics <- rbind(rs_metrics, metrics)
        }
        
        combination_id <- combination_id + 1
        
        results <- rbind(
          results,
          data.frame(
            combination = combination_id,
            validation = "RS",
            FOD_0 = FOD_0,
            log = log,
            dist_thr = dist_thr,
            n_initial = leakage$n_initial,
            n_removed = leakage$n_removed,
            percent_removed = leakage$percent_removed,
            n_remaining = leakage$n_remaining,
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
# 08 - MEILLEUR MODELE LOYO
# ============================================================

best_LOYO <- results |>
  filter(validation == "LOYO") |>
  arrange(RMSE) |>
  slice(1)

cat("\n\n")
cat("====================================================\n")
cat("MEILLEUR MODELE LOYO\n")
cat("====================================================\n")

print(best_LOYO)


# ============================================================
# 09 - MEILLEUR MODELE RANDOM SPLIT
# ============================================================

best_RS <- results |>
  filter(validation == "RS") |>
  arrange(RMSE) |>
  slice(1)

cat("\n\n")
cat("====================================================\n")
cat("MEILLEUR MODELE RANDOM SPLIT\n")
cat("====================================================\n")

print(best_RS)


# ============================================================
# 10 - TOP 10 LOYO
# ============================================================

cat("\n\n")
cat("====================================================\n")
cat("TOP 10 LOYO\n")
cat("====================================================\n")

results |>
  filter(validation == "LOYO") |>
  arrange(RMSE) |>
  select(validation, FOD_0, log, dist_thr, n_removed, n_remaining, max_depth, eta, min_child_weight, subsample, colsample_bytree, nrounds, RMSE, RMSE_sd, R2, MAE) |>
  slice_head(n = 10) |>
  print()


# ============================================================
# 11 - TOP 10 RANDOM SPLIT
# ============================================================

cat("\n\n")
cat("====================================================\n")
cat("TOP 10 RANDOM SPLIT\n")
cat("====================================================\n")

results |>
  filter(validation == "RS") |>
  arrange(RMSE) |>
  select(validation, FOD_0, log, dist_thr, n_removed, n_remaining, max_depth, eta, min_child_weight, subsample, colsample_bytree, nrounds, RMSE, RMSE_sd, R2, MAE) |>
  slice_head(n = 10) |>
  print()


# ============================================================
# 12 - COMPARAISON LOG / FOD / LEAKAGE
# ============================================================

summary_conditions <- results |>
  group_by(validation, FOD_0, log, dist_thr) |>
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

print(summary_conditions)

best_by_condition <- results |>
  group_by(validation, FOD_0, log, dist_thr) |>
  arrange(RMSE) |>
  slice(1) |>
  ungroup() |>
  arrange(validation, RMSE)

print(best_by_condition)


# ============================================================
# 13 - PLOT RMSE EN FONCTION DU LEAKAGE
# ============================================================

ggplot(summary_conditions, aes(x = dist_thr, y = RMSE_mean, color = validation)) +
  geom_line() +
  geom_point(size = 3) +
  facet_grid(
    FOD_0 ~ log,
    labeller = labeller(
      FOD_0 = function(x) paste("FOD_0 =", x),
      log = function(x) paste("log =", x)
    )
  ) +
  theme_bw() +
  labs(
    x = "Distance threshold",
    y = "RMSE moyen",
    color = "Validation",
    title = "Sensibilité du tuning XGBoost au data leakage"
  )


# ============================================================
# 14 - SAUVEGARDE
# ============================================================

write.csv(results, "XGBoost_tuning_LOYO_RS_FOD_LOG_leakage.csv", row.names = FALSE)

write.csv(summary_conditions, "XGBoost_tuning_summary_conditions.csv", row.names = FALSE)