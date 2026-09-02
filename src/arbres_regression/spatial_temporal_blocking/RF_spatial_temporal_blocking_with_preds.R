# =====================================================================
# Random Forest : CV bloquée -- DEUX SCHÉMAS INDÉPENDANTS
# =====================================================================
# 1) Blocage TEMPOREL : plusieurs résolutions testées (1 jour, 3 jours),
#    chacune avec son propre buffer.
# 2) Blocage SPATIAL : plusieurs résolutions testées, CHACUNE avec son
#    propre buffer (buffer proportionnel à la taille de la cellule --
#    un buffer de 100km n'a pas le même sens à 20km qu'à 1500km de
#    résolution).
#
# install.packages(c("dplyr","tidyr","purrr","ggplot2","ranger","FNN","tibble","caret"))

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(ranger)
library(FNN)
library(tibble)
library(caret)

set.seed(42)

# ---------------------------------------------------------------------
# 0. DONNÉES (pipeline de nettoyage inchangé)
# ---------------------------------------------------------------------
freq <- 38
diurnal_period <- 3 # 3 : day, 1: night
dp <- "day"
path_ds <- paste0("F:/data_elise/ds_NASC_pig_ftle_fod/ds_NASC_per_esu_all/ds_NASC_per_esu_pig_ftle_fod_2018_2021_2022_2023_transect_", freq, "kHz_mask9.rds")
datas <- readRDS(path_ds)

datas <- datas[datas$day == diurnal_period, ]

q <- quantile(datas$nasc, probs = c(0.05, 0.95), na.rm = TRUE)
datas <- datas |> dplyr::filter(nasc >= q[1], nasc <= q[2])
datas$nasc <- log10(datas$nasc)

datas$fod <- as.factor(datas$fod)
datas$fod[datas$fod == "NA"] <- NA
datas$fod <- droplevels(datas$fod)

df <- data.frame(
  NASC            = datas$nasc,
  year            = format(datas$time_nasc, "%Y"),
  time            = datas$time_nasc,
  lat             = datas$lat_nasc,
  lon             = datas$lon_nasc,
  fod             = datas$fod,
  ftle            = datas$ftle,
  total_chla      = datas$Chla,
  per_ratio_chla  = datas$Per_Chla,
  but_ratio_chla  = datas$But_Chla,
  fuco_ratio_chla = datas$Fuco_Chla,
  hex_ratio_chla  = datas$Hex_Chla,
  allo_ratio_chla = datas$Allo_Chla,
  zea_ratio_chla  = datas$Zea_Chla,
  chlb_ratio_chla = datas$Chlb_Chla
)

vars_num <- c("NASC", "per_ratio_chla", "but_ratio_chla", "fuco_ratio_chla",
              "hex_ratio_chla", "allo_ratio_chla", "zea_ratio_chla",
              "chlb_ratio_chla", "total_chla", "ftle")

fod_levels <- levels(df$fod)  # sauvegardé explicitement, réutilisé en section 11

df <- df |>
  dplyr::filter(
    if_all(all_of(setdiff(vars_num, "NASC")), ~ !is.na(.)),
    !is.na(fod), fod != "NA"
  )

covariates_num <- setdiff(vars_num, "NASC")
covariates_all <- c(covariates_num, "fod")
response_var   <- "NASC"

cat("Nombre d'observations après filtrage :", nrow(df), "\n")

# ---------------------------------------------------------------------
# 1. CONVERSION LON/LAT -> KM (nécessaire pour le blocage spatial en km)
# ---------------------------------------------------------------------
lon0 <- mean(df$lon, na.rm = TRUE)
lat0 <- mean(df$lat, na.rm = TRUE)
km_per_deg_lat <- 110.574
km_per_deg_lon <- 111.320 * cos(lat0 * pi / 180)
df$x_km <- (df$lon - lon0) * km_per_deg_lon
df$y_km <- (df$lat - lat0) * km_per_deg_lat

# ---------------------------------------------------------------------
# 2. FONCTIONS UTILITAIRES COMMUNES
# ---------------------------------------------------------------------
rmse <- function(obs, pred) sqrt(mean((obs - pred)^2, na.rm = TRUE))

r2_vs_train_mean <- function(obs, pred, train_mean) {
  1 - sum((obs - pred)^2, na.rm = TRUE) / sum((obs - train_mean)^2, na.rm = TRUE)
}

fit_and_eval <- function(train_idx, test_idx, params, data, response = response_var,
                         covs = covariates_all) {
  train_df <- data[train_idx, c(response, covs)]
  test_df  <- data[test_idx,  c(response, covs)]
  train_df <- train_df[stats::complete.cases(train_df), ]
  test_df  <- test_df[stats::complete.cases(test_df), ]
  form  <- stats::as.formula(paste(response, "~ ."))
  model <- ranger(form, data = train_df, mtry = params$mtry,
                  min.node.size = params$min.node.size, num.trees = params$num.trees,
                  num.threads = max(1, parallel::detectCores() - 1))
  
  train_mean <- mean(train_df[[response]], na.rm = TRUE)
  pred_test  <- predict(model, test_df)$predictions
  pred_train <- predict(model, train_df)$predictions
  
  tibble(
    rmse_test  = rmse(test_df[[response]],  pred_test),
    rmse_train = rmse(train_df[[response]], pred_train),
    r2_test    = r2_vs_train_mean(test_df[[response]],  pred_test,  train_mean),
    r2_train   = r2_vs_train_mean(train_df[[response]], pred_train, train_mean)
  )
}

fast_min_gap_days <- function(train_dates, test_dates) {
  train_num   <- as.numeric(as.Date(train_dates))
  test_sorted <- sort(unique(as.numeric(as.Date(test_dates))))
  idx       <- findInterval(train_num, test_sorted)
  left_val  <- test_sorted[pmax(idx, 1)]
  right_val <- test_sorted[pmin(idx + 1, length(test_sorted))]
  pmin(abs(train_num - left_val), abs(train_num - right_val))
}

# ---------------------------------------------------------------------
# 3A. BLOCAGE TEMPOREL : blocs de N jours + tampon
#     -- inchangé dans sa logique, block_days et buffer_days sont
#        maintenant passés par résolution (voir section 4).
# ---------------------------------------------------------------------
assign_temporal_block <- function(data, block_days) {
  t0 <- min(as.numeric(as.Date(data$time)))
  data %>% mutate(temporal_block = floor((as.numeric(as.Date(time)) - t0) / block_days))
}

make_temporal_fold <- function(data, block_id, buffer_days) {
  test_idx            <- which(data$temporal_block == block_id)
  candidate_train_idx <- which(data$temporal_block != block_id)
  gaps      <- fast_min_gap_days(data$time[candidate_train_idx], data$time[test_idx])
  train_idx <- candidate_train_idx[gaps > buffer_days]
  list(train = train_idx, test = test_idx)
}

build_temporal_folds <- function(data, block_days, buffer_days,
                                 min_block_n = 50,
                                 max_folds_fraction = 0.3,
                                 min_folds = 5,
                                 max_folds_abs = 30,
                                 seed = 42) {
  data_blocked <- assign_temporal_block(data, block_days)
  blocks_count <- data_blocked %>% count(temporal_block)
  blocks_used  <- blocks_count %>% filter(n >= min_block_n) %>% pull(temporal_block)
  
  # Nombre de folds cible = fraction des blocs disponibles, borné
  n_target <- round(length(blocks_used) * max_folds_fraction)
  n_target <- max(min_folds, min(n_target, max_folds_abs, length(blocks_used)))
  
  set.seed(seed)
  if (length(blocks_used) > n_target) blocks_used <- sample(blocks_used, n_target)
  
  folds <- map(blocks_used, ~ make_temporal_fold(data_blocked, .x, buffer_days))
  names(folds) <- paste0("t_", blocks_used)
  
  valid <- map_lgl(folds, ~ length(.x$train) > 0 && length(.x$test) > 0)
  if (any(!valid)) {
    cat("  [!] ", sum(!valid), "bloc(s) temporel(s) ignoré(s) : train ou test vide après buffer\n")
  }
  folds <- folds[valid]
  
  list(data = data_blocked, folds = folds,
       n_blocks_total = nrow(blocks_count), n_blocks_used = length(folds),
       n_folds_target = n_target)
}

# ---------------------------------------------------------------------
# 3B. BLOCAGE SPATIAL : grille en km (résolutions testées), CHACUNE avec
#     son propre buffer + nombre de folds adapté au nombre de blocs
#     disponibles.
# ---------------------------------------------------------------------
assign_spatial_block <- function(data, cellsize_x, cellsize_y) {
  data %>%
    mutate(
      block_x       = floor(x_km / cellsize_x),
      block_y       = floor(y_km / cellsize_y),
      spatial_block = paste(block_x, block_y, sep = "_")
    )
}

make_spatial_fold_buffered <- function(data, block_id, buffer_km) {
  test_idx            <- which(data$spatial_block == block_id)
  candidate_train_idx <- which(data$spatial_block != block_id)
  test_mat      <- as.matrix(data[test_idx, c("x_km", "y_km")])
  candidate_mat <- as.matrix(data[candidate_train_idx, c("x_km", "y_km")])
  dist_to_test <- FNN::knnx.dist(data = test_mat, query = candidate_mat, k = 1)[, 1]
  train_idx <- candidate_train_idx[dist_to_test > buffer_km]
  list(train = train_idx, test = test_idx)
}

build_spatial_folds <- function(data, cellsize_x, cellsize_y, buffer_km,
                                min_block_n = 50,
                                max_folds_fraction = 0.3,
                                min_folds = 5,
                                max_folds_abs = 30,
                                seed = 42) {
  data_blocked <- assign_spatial_block(data, cellsize_x, cellsize_y)
  blocks_count <- data_blocked %>% count(spatial_block)
  blocks_used  <- blocks_count %>% filter(n >= min_block_n) %>% pull(spatial_block)
  
  n_target <- round(length(blocks_used) * max_folds_fraction)
  n_target <- max(min_folds, min(n_target, max_folds_abs, length(blocks_used)))
  
  set.seed(seed)
  if (length(blocks_used) > n_target) blocks_used <- sample(blocks_used, n_target)
  
  folds <- map(blocks_used, ~ make_spatial_fold_buffered(data_blocked, .x, buffer_km))
  names(folds) <- paste0("s_", blocks_used)
  
  valid <- map_lgl(folds, ~ length(.x$train) > 0 && length(.x$test) > 0)
  if (any(!valid)) {
    cat("  [!] ", sum(!valid), "bloc(s) spatial/spatiaux ignoré(s) : train ou test vide après buffer\n")
  }
  folds <- folds[valid]
  
  list(data = data_blocked, folds = folds,
       n_blocks_total = nrow(blocks_count), n_blocks_used = length(folds),
       n_folds_target = n_target, buffer_km = buffer_km)
}

# ---------------------------------------------------------------------
# 4. CONSTRUCTION DES DEUX SCHÉMAS + COMPARAISON BASELINE
# ---------------------------------------------------------------------

# Paramètres d'adaptation du nombre de folds à la résolution
max_folds_fraction <- 0.3
min_folds          <- 5
max_folds_abs       <- 30

# ---- Résolutions SPATIALES testées, chacune avec son propre buffer ----
# Le buffer est choisi proportionnellement à la taille de cellule : trop
# petit par rapport à la cellule, il ne protège pas vraiment contre la
# fuite spatiale ; trop grand par rapport à une cellule fine, il laisse
# quasiment plus de train disponible.
resolutions_to_test <- list(
  list(lon_km = 20,   lat_km = 20,   buffer_km = 5),
  list(lon_km = 60,   lat_km = 60,   buffer_km = 5),
  list(lon_km = 500,  lat_km = 500,  buffer_km = 100),
  list(lon_km = 1500, lat_km = 1000, buffer_km = 100)
)

# ---- Résolutions TEMPORELLES testées, chacune avec son propre buffer ----
# 1 jour de bloc -> buffer plus court (ex. 1 jour) ; 3 jours de bloc ->
# buffer un peu plus large (ex. 2 jours), à ajuster selon la portée
# temporelle réelle de l'autocorrélation résiduelle (cf. variogrammes
# calculés précédemment).
temporal_resolutions_to_test <- list(
  list(block_days = 1, buffer_days = 1),
  list(block_days = 3, buffer_days = 2)
)

built_spatial_list <- map(resolutions_to_test, ~
                            build_spatial_folds(df, .x$lon_km, .x$lat_km, .x$buffer_km,
                                                max_folds_fraction = max_folds_fraction,
                                                min_folds = min_folds,
                                                max_folds_abs = max_folds_abs)
)
names(built_spatial_list) <- map_chr(resolutions_to_test, ~ paste0(.x$lon_km, "x", .x$lat_km, "km_buffer", .x$buffer_km, "km"))

built_temporal_list <- map(temporal_resolutions_to_test, ~
                             build_temporal_folds(df, .x$block_days, .x$buffer_days,
                                                  max_folds_fraction = max_folds_fraction,
                                                  min_folds = min_folds,
                                                  max_folds_abs = max_folds_abs)
)
names(built_temporal_list) <- map_chr(temporal_resolutions_to_test, ~ paste0(.x$block_days, "j_buffer", .x$buffer_days, "j"))

# Diagnostic : combien de folds ont été retenus par résolution
cat("\n=== Nombre de folds retenus par résolution SPATIALE (adaptatif) ===\n")
walk2(built_spatial_list, names(built_spatial_list), function(built, label) {
  cat(sprintf("  %-25s : %d blocs disponibles -> %d folds retenus (cible : %d)\n",
              label, built$n_blocks_total, built$n_blocks_used, built$n_folds_target))
})

cat("\n=== Nombre de folds retenus par résolution TEMPORELLE (adaptatif) ===\n")
walk2(built_temporal_list, names(built_temporal_list), function(built, label) {
  cat(sprintf("  %-25s : %d blocs disponibles -> %d folds retenus (cible : %d)\n",
              label, built$n_blocks_total, built$n_blocks_used, built$n_folds_target))
})

all_built <- c(built_temporal_list, built_spatial_list)

baseline_params <- list(mtry = 3, min.node.size = 5, num.trees = 300)

comparison_summary <- imap_dfr(all_built, function(built, config_name) {
  if (length(built$folds) == 0) {
    cat("  [!] Schéma '", config_name, "' ignoré dans la comparaison : 0 fold valide\n")
    return(tibble())
  }
  fold_perf <- map_dfr(names(built$folds), function(fid) {
    f <- built$folds[[fid]]
    fit_and_eval(f$train, f$test, baseline_params, built$data) %>%
      mutate(fold_id = fid, n_train = length(f$train), n_test = length(f$test))
  })
  fold_perf %>%
    summarise(
      n_blocks_total  = built$n_blocks_total,
      n_folds_used    = length(built$folds),
      mean_n_train    = mean(n_train),
      mean_n_test     = mean(n_test),
      mean_rmse_test  = mean(rmse_test),
      sd_rmse_test    = sd(rmse_test),
      mean_rmse_train = mean(rmse_train),
      mean_r2_test    = mean(r2_test),
      sd_r2_test      = sd(r2_test)
    ) %>%
    mutate(config = config_name)
})

cat("\n=== Comparaison de tous les schémas de blocage (temporel x2, spatial x4) ===\n")
print(comparison_summary)

ggplot(comparison_summary, aes(x = config, y = mean_rmse_test)) +
  geom_col(fill = "steelblue") +
  geom_errorbar(aes(ymin = mean_rmse_test - sd_rmse_test, ymax = mean_rmse_test + sd_rmse_test), width = 0.2) +
  labs(title = "RMSE test par schéma de blocage",
       subtitle = "Barres d'erreur = écart-type inter-fold -- buffer et nb de folds adaptés à chaque résolution",
       x = NULL, y = "RMSE test (moyenne inter-fold)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggplot(comparison_summary, aes(x = config, y = mean_r2_test)) +
  geom_col(fill = "darkgreen") +
  geom_errorbar(aes(ymin = mean_r2_test - sd_r2_test, ymax = mean_r2_test + sd_r2_test), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "R² test (vs moyenne du train) par schéma de blocage",
       subtitle = "Sous la ligne rouge = pas mieux qu'une prédiction naïve (moyenne du train)",
       x = NULL, y = "R² test") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# ---------------------------------------------------------------------
# 5. PIPELINE COMPLET (tuning, courbe d'apprentissage, extrapolation,
#    importance des variables) -- réutilisable pour n'importe quel schéma
# ---------------------------------------------------------------------
tuning_grid <- expand.grid(mtry = c(2, 3, 4), min.node.size = c(5, 10), num.trees = 300)

run_pipeline <- function(built, label, tuning_grid, fractions = seq(0.2, 1, by = 0.2)) {
  folds <- built$folds
  data  <- built$data
  
  tuning_results <- expand_grid(grid_id = seq_len(nrow(tuning_grid)), fold_id = names(folds)) %>%
    mutate(params = map(grid_id, ~ tuning_grid[.x, ])) %>%
    mutate(res = map2(fold_id, params, ~ fit_and_eval(folds[[.x]]$train, folds[[.x]]$test, .y, data))) %>%
    unnest(res)
  
  tuning_summary <- tuning_results %>%
    group_by(grid_id) %>%
    summarise(mean_rmse_test = mean(rmse_test), sd_rmse_test = sd(rmse_test),
              mean_rmse_train = mean(rmse_train),
              mean_r2_test = mean(r2_test), .groups = "drop") %>%
    left_join(tuning_grid %>% mutate(grid_id = row_number()), by = "grid_id") %>%
    arrange(mean_rmse_test)
  
  best_params <- tuning_grid[tuning_summary$grid_id[1], ]
  
  learning_curves <- map_dfr(folds, function(f) {
    map_dfr(fractions, function(frac) {
      n_sub   <- max(20, round(frac * length(f$train)))
      sub_idx <- sample(f$train, n_sub)
      fit_and_eval(sub_idx, f$test, best_params, data) %>% mutate(fraction = frac, n_train = n_sub)
    })
  }, .id = "fold_id")
  
  learning_curve_summary <- learning_curves %>%
    group_by(fraction) %>%
    summarise(mean_rmse_train = mean(rmse_train), mean_rmse_test = mean(rmse_test), .groups = "drop")
  
  extrapolation_results <- map_dfr(folds, function(f) {
    train_raw <- as.matrix(data[f$train, covariates_num])
    test_raw  <- as.matrix(data[f$test,  covariates_num])
    center <- colMeans(train_raw, na.rm = TRUE)
    scale_ <- apply(train_raw, 2, sd, na.rm = TRUE)
    train_mat <- scale(train_raw, center = center, scale = scale_)
    test_mat  <- scale(test_raw,  center = center, scale = scale_)
    train_nn <- knn.dist(train_mat, k = 1)[, 1]
    test_nn  <- knnx.dist(train_mat, test_mat, k = 1)[, 1]
    tibble(extrapolation_index = mean(test_nn) / mean(train_nn))
  }, .id = "fold_id")
  
  ctrl <- trainControl(method = "cv", index = map(folds, "train"), indexOut = map(folds, "test"),
                       savePredictions = "final")
  
  model_caret <- train(
    x = as.data.frame(data[, covariates_all]),
    y = data[[response_var]],
    method = "ranger", trControl = ctrl,
    tuneGrid = data.frame(mtry = best_params$mtry, splitrule = "variance",
                          min.node.size = best_params$min.node.size),
    importance = "permutation", num.trees = best_params$num.trees
  )
  
  importance_df <- as.data.frame(caret::varImp(model_caret)$importance) %>%
    tibble::rownames_to_column("variable") %>%
    mutate(config = label)
  
  data_final <- data[, c(response_var, covariates_all)]
  data_final <- data_final[stats::complete.cases(data_final), ]
  form_final <- stats::as.formula(paste(response_var, "~ ."))
  model_final <- ranger(
    form_final, data = data_final,
    mtry = best_params$mtry, min.node.size = best_params$min.node.size,
    num.trees = best_params$num.trees,
    importance = "permutation",
    num.threads = max(1, parallel::detectCores() - 1)
  )
  
  list(label = label, tuning_summary = tuning_summary, best_params = best_params,
       learning_curve_summary = learning_curve_summary,
       extrapolation_results = extrapolation_results,
       model_caret = model_caret, model_final = model_final,
       importance_df = importance_df)
}

# ---------------------------------------------------------------------
# 6. APPLICATION : UNE résolution temporelle vs UNE résolution spatiale
# ---------------------------------------------------------------------
SELECTED_TEMPORAL_LABEL <- "3j_buffer2j"
SELECTED_SPATIAL_LABEL  <- "60x60km_buffer5km"

results_temporal <- run_pipeline(built_temporal_list[[SELECTED_TEMPORAL_LABEL]],
                                 paste0("temporel_", SELECTED_TEMPORAL_LABEL),
                                 tuning_grid)
results_spatial  <- run_pipeline(built_spatial_list[[SELECTED_SPATIAL_LABEL]],
                                 paste0("spatial_", SELECTED_SPATIAL_LABEL),
                                 tuning_grid)

cat("\nMeilleurs hyperparamètres -- bloc temporel :\n"); print(results_temporal$best_params)
cat("Meilleurs hyperparamètres -- bloc spatial :\n");   print(results_spatial$best_params)

# ---- Sauvegarde des modèles finaux (RF, entraînés sur toutes les données) ----
model_output_dir <- "outputs_cv_diagnostics/models"
dir.create(model_output_dir, showWarnings = FALSE, recursive = TRUE)

saveRDS(results_temporal$model_final, file.path(model_output_dir, paste0("rf_model_temporel_", SELECTED_TEMPORAL_LABEL, ".rds")))
saveRDS(results_spatial$model_final,  file.path(model_output_dir, paste0("rf_model_spatial_", SELECTED_SPATIAL_LABEL, ".rds")))
saveRDS(fod_levels, file.path(model_output_dir, "fod_levels.rds"))

cat("\nModèles sauvegardés dans :", normalizePath(model_output_dir), "\n")

# ---------------------------------------------------------------------
# 7. COMPARAISON DE LA STABILITÉ DE L'IMPORTANCE DES VARIABLES
# ---------------------------------------------------------------------
importance_comparison <- bind_rows(results_temporal$importance_df, results_spatial$importance_df)

ggplot(importance_comparison, aes(x = variable, y = Overall, fill = config)) +
  geom_col(position = "dodge") +
  coord_flip() +
  labs(title = "Stabilité de l'importance des variables : bloc temporel vs bloc spatial",
       x = NULL, y = "Importance (permutation)", fill = NULL) +
  theme_minimal()

# =====================================================================
# 8. DIAGNOSTICS COMPLETS PAR FOLD
# =====================================================================

built           <- built_spatial_list[[SELECTED_SPATIAL_LABEL]]   # ou built_temporal_list[[...]]
results_selected <- results_spatial                                # ou results_temporal

folds  <- built$folds
data   <- built$data
params <- results_selected$best_params

diagnostic_label <- paste0("blocage spatial ", SELECTED_SPATIAL_LABEL,
                           " (", length(folds), " folds)")
# Si tu utilises built_temporal_list à la place, remplace par :
# diagnostic_label <- paste0("blocage temporel ", SELECTED_TEMPORAL_LABEL,
#                             " (", length(folds), " folds)")

compute_fold_diagnostics <- function(fold, data, params,
                                     response = response_var,
                                     covs_model = covariates_all,
                                     covs_num   = covariates_num) {
  
  needed_cols <- unique(c(response, covs_model, "x_km", "y_km"))
  train_df <- data[fold$train, needed_cols]
  test_df  <- data[fold$test,  needed_cols]
  train_df <- train_df[stats::complete.cases(train_df[, c(response, covs_model)]), ]
  test_df  <- test_df[stats::complete.cases(test_df[, c(response, covs_model)]), ]
  
  form  <- stats::as.formula(paste(response, "~", paste(covs_model, collapse = " + ")))
  model <- ranger(form, data = train_df[, c(response, covs_model)],
                  mtry = params$mtry, min.node.size = params$min.node.size,
                  num.trees = params$num.trees,
                  num.threads = max(1, parallel::detectCores() - 1))
  
  pred_test  <- predict(model, test_df[, covs_model])$predictions
  train_mean <- mean(train_df[[response]])
  
  obs_pred <- tibble(obs = test_df[[response]], pred = pred_test)
  
  metrics <- tibble(
    n_train              = nrow(train_df),
    n_test               = nrow(test_df),
    rmse_test            = sqrt(mean((obs_pred$obs - obs_pred$pred)^2)),
    r2_test              = 1 - sum((obs_pred$obs - obs_pred$pred)^2) / sum((obs_pred$obs - train_mean)^2),
    var_intra_fold_test  = var(obs_pred$obs),
    sd_intra_fold_test   = sd(obs_pred$obs)
  )
  
  geo_dist <- FNN::knnx.dist(
    data  = as.matrix(train_df[, c("x_km", "y_km")]),
    query = as.matrix(test_df[,  c("x_km", "y_km")]),
    k = 1
  )[, 1]
  metrics$mean_geo_dist_km <- mean(geo_dist)
  
  train_num <- as.matrix(train_df[, covs_num])
  test_num  <- as.matrix(test_df[,  covs_num])
  center <- colMeans(train_num, na.rm = TRUE)
  scale_ <- apply(train_num, 2, sd, na.rm = TRUE)
  train_num_s <- scale(train_num, center = center, scale = scale_)
  test_num_s  <- scale(test_num,  center = center, scale = scale_)
  train_nn <- FNN::knn.dist(train_num_s, k = 1)[, 1]
  test_nn  <- FNN::knnx.dist(train_num_s, test_num_s, k = 1)[, 1]
  metrics$mean_covariate_dist <- mean(test_nn)
  metrics$extrapolation_index <- mean(test_nn) / mean(train_nn)
  
  covariate_stats <- bind_rows(
    tibble(set = "test",  variable = covs_num,
           mean = colMeans(test_num,  na.rm = TRUE), sd = apply(test_num,  2, sd, na.rm = TRUE)),
    tibble(set = "train", variable = covs_num,
           mean = colMeans(train_num, na.rm = TRUE), sd = apply(train_num, 2, sd, na.rm = TRUE))
  )
  
  list(obs_pred = obs_pred, metrics = metrics, covariate_stats = covariate_stats)
}

compute_fod_distribution <- function(fold, data) {
  bind_rows(
    tibble(fod = data$fod[fold$train], set = "train"),
    tibble(fod = data$fod[fold$test],  set = "test")
  ) %>%
    filter(!is.na(fod)) %>%
    count(set, fod, .drop = FALSE) %>%
    group_by(set) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup()
}

compute_ftle_values <- function(fold, data) {
  bind_rows(
    tibble(ftle = data$ftle[fold$train], set = "train"),
    tibble(ftle = data$ftle[fold$test],  set = "test")
  ) %>%
    filter(is.finite(ftle))
}

all_results <- imap(folds, ~ compute_fold_diagnostics(.x, data, params))
fod_dist_all  <- imap_dfr(folds, ~ compute_fod_distribution(.x, data) %>% mutate(fold_id = .y))
ftle_vals_all <- imap_dfr(folds, ~ compute_ftle_values(.x, data)      %>% mutate(fold_id = .y))

obs_pred_all        <- imap_dfr(all_results, ~ mutate(.x$obs_pred, fold_id = .y))
metrics_all         <- imap_dfr(all_results, ~ mutate(.x$metrics, fold_id = .y))
covariate_stats_all <- imap_dfr(all_results, ~ mutate(.x$covariate_stats, fold_id = .y))

cat("\n=== Métriques par fold --", diagnostic_label, "===\n")
print(metrics_all)

global_rmse <- sqrt(mean((obs_pred_all$obs - obs_pred_all$pred)^2))
global_r2   <- 1 - sum((obs_pred_all$obs - obs_pred_all$pred)^2) /
  sum((obs_pred_all$obs - mean(obs_pred_all$obs))^2)

cat("\nRMSE globale (poolée, tous folds) :", round(global_rmse, 3), "\n")
cat("R² global (poolé, tous folds)     :", round(global_r2, 3), "\n")

ggplot(obs_pred_all, aes(x = obs, y = pred, color = fold_id)) +
  geom_point(alpha = 0.4, size = 1.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  annotate("label",
           x = min(obs_pred_all$obs, na.rm = TRUE),
           y = max(obs_pred_all$pred, na.rm = TRUE),
           hjust = 0, vjust = 1,
           label = sprintf("RMSE = %.3f\nR² = %.3f", global_rmse, global_r2),
           fill = "white", alpha = 0.8, size = 3.5) +
  labs(title = "NASC observé vs prédit (CV bloquée)",
       subtitle = diagnostic_label,
       x = "NASC observé (log10)", y = "NASC prédit (log10)", color = "Fold") +
  theme_minimal()

ggplot(results_selected$learning_curve_summary, aes(x = fraction)) +
  geom_line(aes(y = mean_rmse_train, color = "Train")) +
  geom_point(aes(y = mean_rmse_train, color = "Train")) +
  geom_line(aes(y = mean_rmse_test, color = "Test")) +
  geom_point(aes(y = mean_rmse_test, color = "Test")) +
  labs(title = "Courbe d'apprentissage", subtitle = diagnostic_label,
       x = "Fraction du train utilisée",
       y = "RMSE", color = NULL) +
  theme_minimal()

ggplot(metrics_all, aes(x = fold_id, y = r2_test)) +
  geom_col(fill = "darkgreen") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "R² par fold (référencé sur la moyenne du train)",
       subtitle = diagnostic_label,
       x = "Fold", y = "R² test") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(metrics_all, aes(x = fold_id, y = var_intra_fold_test)) +
  geom_col(fill = "steelblue") +
  labs(title = "Variance du NASC observé, intra-fold (test)",
       subtitle = diagnostic_label,
       x = "Fold", y = "Variance (log10 NASC)") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

dist_long <- metrics_all %>%
  select(fold_id, mean_geo_dist_km, mean_covariate_dist) %>%
  pivot_longer(-fold_id, names_to = "metric", values_to = "value")

ggplot(dist_long, aes(x = fold_id, y = value, fill = metric)) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~metric, scales = "free_y", labeller = as_labeller(c(
    mean_geo_dist_km    = "Distance géographique moyenne (km)",
    mean_covariate_dist = "Distance covariables test->train (standardisée)"
  ))) +
  labs(title = "Distance test -> train, par fold", subtitle = diagnostic_label,
       x = "Fold", y = NULL) +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(covariate_stats_all, aes(x = fold_id, y = mean, color = set)) +
  geom_point(position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2,
                position = position_dodge(width = 0.4)) +
  facet_wrap(~variable, scales = "free_y") +
  labs(title = "Moyenne ± écart-type de chaque covariable, par fold (train vs test)",
       subtitle = diagnostic_label,
       x = "Fold", y = NULL, color = NULL) +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(fod_dist_all, aes(x = fod, y = prop, fill = set)) +
  geom_col(position = "dodge") +
  facet_wrap(~fold_id) +
  labs(title = "Distribution des clusters FOD, train vs test, par fold",
       subtitle = diagnostic_label,
       x = "Cluster FOD", y = "Proportion", fill = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(ftle_vals_all, aes(x = ftle, fill = set, color = set)) +
  geom_density(alpha = 0.3) +
  facet_wrap(~fold_id, scales = "free_y") +
  labs(title = "Distribution du FTLE, train vs test, par fold",
       subtitle = diagnostic_label,
       x = "FTLE", y = "Densité", fill = NULL, color = NULL) +
  theme_minimal()

# =====================================================================
# 9. CARTE DES POINTS NASC COLORÉS PAR FOLD
# =====================================================================

built_map_s  <- built_spatial_list[[SELECTED_SPATIAL_LABEL]]
data_map_s   <- built_map_s$data
data_map_s$row_id <- seq_len(nrow(data_map_s))

fold_lookup_s <- imap_dfr(built_map_s$folds, ~ tibble(row_id = .x$test, fold_id = .y))
data_map_s <- data_map_s %>%
  left_join(fold_lookup_s, by = "row_id") %>%
  mutate(fold_id = ifelse(is.na(fold_id), "Non utilisé (bloc trop petit)", fold_id))

ggplot(data_map_s, aes(x = lon, y = lat, color = fold_id)) +
  geom_point(size = 0.8, alpha = 0.7) +
  coord_quickmap() +
  labs(title = paste0("Points NASC colorés par fold -- blocage spatial (", SELECTED_SPATIAL_LABEL, ")"),
       x = "Longitude", y = "Latitude", color = "Fold") +
  theme_minimal()

built_map_t  <- built_temporal_list[[SELECTED_TEMPORAL_LABEL]]
data_map_t   <- built_map_t$data
data_map_t$row_id <- seq_len(nrow(data_map_t))

fold_lookup_t <- imap_dfr(built_map_t$folds, ~ tibble(row_id = .x$test, fold_id = .y))
data_map_t <- data_map_t %>%
  left_join(fold_lookup_t, by = "row_id") %>%
  mutate(fold_id = ifelse(is.na(fold_id), "Non utilisé (bloc trop petit)", fold_id))

ggplot(data_map_t, aes(x = lon, y = lat, color = fold_id)) +
  geom_point(size = 0.8, alpha = 0.7) +
  coord_quickmap() +
  labs(title = paste0("Points NASC colorés par fold -- blocage temporel (", SELECTED_TEMPORAL_LABEL, ")"),
       x = "Longitude", y = "Latitude", color = "Fold") +
  theme_minimal()

status_by_fold <- imap_dfr(built_map_s$folds, function(fold, fid) {
  df_status <- data_map_s %>% select(lon, lat)
  df_status$status  <- "Exclu par le buffer"
  df_status$status[fold$train] <- "Train conservé"
  df_status$status[fold$test]  <- "Test (bloc isolé)"
  df_status$fold_id <- fid
  df_status
})

ggplot(status_by_fold, aes(x = lon, y = lat, color = status)) +
  geom_point(size = 0.5, alpha = 0.6) +
  scale_color_manual(values = c("Test (bloc isolé)"   = "red",
                                "Train conservé"      = "steelblue",
                                "Exclu par le buffer" = "grey80")) +
  coord_quickmap() +
  facet_wrap(~fold_id) +
  labs(title = "Effet du buffer spatial -- tous les folds",
       subtitle = paste0("Résolution ", SELECTED_SPATIAL_LABEL),
       x = "Longitude", y = "Latitude", color = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom")

# =====================================================================
# 10. GÉNÉRATION + SAUVEGARDE DE TOUS LES PLOTS + MODÈLES,
#     POUR TOUTES LES RÉSOLUTIONS (spatiales ET temporelles)
# =====================================================================

output_root <- "outputs_cv_diagnostics"
dir.create(output_root, showWarnings = FALSE)

generate_and_save_diagnostics <- function(built, label, tuning_grid, output_root,
                                          buffer_description = "") {
  
  res_dir <- file.path(output_root, label)
  dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)
  diag_label <- paste0(label, if (nchar(buffer_description)) paste0(" (", buffer_description, ")") else "")
  
  cat("\n==== Résolution :", label, "-- dossier :", res_dir, "====\n")
  
  results_res <- run_pipeline(built, label, tuning_grid)
  
  saveRDS(results_res$model_final, file.path(res_dir, "rf_model_final.rds"))
  cat("  -> Modèle RF sauvegardé :", file.path(res_dir, "rf_model_final.rds"), "\n")
  
  all_fold_results <- imap(built$folds, ~ compute_fold_diagnostics(.x, built$data, results_res$best_params))
  obs_pred_res        <- imap_dfr(all_fold_results, ~ mutate(.x$obs_pred, fold_id = .y))
  metrics_res         <- imap_dfr(all_fold_results, ~ mutate(.x$metrics, fold_id = .y))
  covariate_stats_res <- imap_dfr(all_fold_results, ~ mutate(.x$covariate_stats, fold_id = .y))
  
  global_rmse_res <- sqrt(mean((obs_pred_res$obs - obs_pred_res$pred)^2))
  global_r2_res   <- 1 - sum((obs_pred_res$obs - obs_pred_res$pred)^2) /
    sum((obs_pred_res$obs - mean(obs_pred_res$obs))^2)
  
  p1 <- ggplot(obs_pred_res, aes(x = obs, y = pred, color = fold_id)) +
    geom_point(alpha = 0.4, size = 1.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
    annotate("label", x = min(obs_pred_res$obs, na.rm = TRUE), y = max(obs_pred_res$pred, na.rm = TRUE),
             hjust = 0, vjust = 1, label = sprintf("RMSE = %.3f\nR\u00b2 = %.3f", global_rmse_res, global_r2_res),
             fill = "white", alpha = 0.8, size = 3.5) +
    labs(title = "NASC observé vs prédit", subtitle = diag_label,
         x = "NASC observé (log10)", y = "NASC prédit (log10)", color = "Fold") +
    theme_minimal()
  ggsave(file.path(res_dir, "01_obs_vs_pred.png"), p1, width = 7, height = 5, dpi = 150)
  
  p2 <- ggplot(results_res$learning_curve_summary, aes(x = fraction)) +
    geom_line(aes(y = mean_rmse_train, color = "Train")) +
    geom_point(aes(y = mean_rmse_train, color = "Train")) +
    geom_line(aes(y = mean_rmse_test, color = "Test")) +
    geom_point(aes(y = mean_rmse_test, color = "Test")) +
    labs(title = "Courbe d'apprentissage", subtitle = diag_label,
         x = "Fraction du train utilisée", y = "RMSE", color = NULL) +
    theme_minimal()
  ggsave(file.path(res_dir, "02_learning_curve.png"), p2, width = 7, height = 5, dpi = 150)
  
  p3 <- ggplot(metrics_res, aes(x = fold_id, y = r2_test)) +
    geom_col(fill = "darkgreen") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    labs(title = "R² par fold (référencé sur la moyenne du train)", subtitle = diag_label,
         x = "Fold", y = "R² test") +
    theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(res_dir, "03_r2_par_fold.png"), p3, width = 7, height = 5, dpi = 150)
  
  p4 <- ggplot(metrics_res, aes(x = fold_id, y = var_intra_fold_test)) +
    geom_col(fill = "steelblue") +
    labs(title = "Variance du NASC observé, intra-fold (test)", subtitle = diag_label,
         x = "Fold", y = "Variance (log10 NASC)") +
    theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(res_dir, "04_variance_intra_fold.png"), p4, width = 7, height = 5, dpi = 150)
  
  dist_long_res <- metrics_res %>%
    select(fold_id, mean_geo_dist_km, mean_covariate_dist) %>%
    pivot_longer(-fold_id, names_to = "metric", values_to = "value")
  p5 <- ggplot(dist_long_res, aes(x = fold_id, y = value, fill = metric)) +
    geom_col(show.legend = FALSE) +
    facet_wrap(~metric, scales = "free_y", labeller = as_labeller(c(
      mean_geo_dist_km    = "Distance géographique moyenne (km)",
      mean_covariate_dist = "Distance covariables test->train (standardisée)"
    ))) +
    labs(title = "Distance test -> train, par fold", subtitle = diag_label, x = "Fold", y = NULL) +
    theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(res_dir, "05_distance_test_train.png"), p5, width = 9, height = 5, dpi = 150)
  
  p6 <- ggplot(covariate_stats_res, aes(x = fold_id, y = mean, color = set)) +
    geom_point(position = position_dodge(width = 0.4)) +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2, position = position_dodge(width = 0.4)) +
    facet_wrap(~variable, scales = "free_y") +
    labs(title = "Moyenne ± écart-type de chaque covariable, par fold (train vs test)", subtitle = diag_label,
         x = "Fold", y = NULL, color = NULL) +
    theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(res_dir, "06_covariable_stats.png"), p6, width = 10, height = 8, dpi = 150)
  
  fod_dist_res  <- imap_dfr(built$folds, ~ compute_fod_distribution(.x, built$data) %>% mutate(fold_id = .y))
  p6b <- ggplot(fod_dist_res, aes(x = fod, y = prop, fill = set)) +
    geom_col(position = "dodge") +
    facet_wrap(~fold_id) +
    labs(title = "Distribution des clusters FOD, train vs test, par fold", subtitle = diag_label,
         x = "Cluster FOD", y = "Proportion", fill = NULL) +
    theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(res_dir, "06b_fod_distribution.png"), p6b, width = 9, height = 7, dpi = 150)
  
  ftle_vals_res <- imap_dfr(built$folds, ~ compute_ftle_values(.x, built$data) %>% mutate(fold_id = .y))
  p6c <- ggplot(ftle_vals_res, aes(x = ftle, fill = set, color = set)) +
    geom_density(alpha = 0.3) +
    facet_wrap(~fold_id, scales = "free_y") +
    labs(title = "Distribution du FTLE, train vs test, par fold", subtitle = diag_label,
         x = "FTLE", y = "Densité", fill = NULL, color = NULL) +
    theme_minimal()
  ggsave(file.path(res_dir, "06c_ftle_distribution.png"), p6c, width = 9, height = 7, dpi = 150)
  
  data_map_res <- built$data
  data_map_res$row_id <- seq_len(nrow(data_map_res))
  fold_lookup_res <- imap_dfr(built$folds, ~ tibble(row_id = .x$test, fold_id = .y))
  data_map_res <- data_map_res %>%
    left_join(fold_lookup_res, by = "row_id") %>%
    mutate(fold_id = ifelse(is.na(fold_id), "Non utilisé (bloc trop petit)", fold_id))
  p7 <- ggplot(data_map_res, aes(x = lon, y = lat, color = fold_id)) +
    geom_point(size = 0.8, alpha = 0.7) +
    coord_quickmap() +
    labs(title = "Points NASC colorés par fold", subtitle = diag_label,
         x = "Longitude", y = "Latitude", color = "Fold") +
    theme_minimal()
  ggsave(file.path(res_dir, "07_carte_points_par_fold.png"), p7, width = 8, height = 6, dpi = 150)
  
  status_by_fold_res <- imap_dfr(built$folds, function(fold, fid) {
    df_status <- data_map_res %>% select(lon, lat)
    df_status$status  <- "Exclu par le buffer"
    df_status$status[fold$train] <- "Train conservé"
    df_status$status[fold$test]  <- "Test (bloc isolé)"
    df_status$fold_id <- fid
    df_status
  })
  p8 <- ggplot(status_by_fold_res, aes(x = lon, y = lat, color = status)) +
    geom_point(size = 0.5, alpha = 0.6) +
    scale_color_manual(values = c("Test (bloc isolé)"   = "red",
                                  "Train conservé"      = "steelblue",
                                  "Exclu par le buffer" = "grey80")) +
    coord_quickmap() +
    facet_wrap(~fold_id) +
    labs(title = "Effet du buffer spatial -- tous les folds", subtitle = diag_label,
         x = "Longitude", y = "Latitude", color = NULL) +
    theme_minimal() + theme(legend.position = "bottom")
  ggsave(file.path(res_dir, "08_effet_buffer_tous_folds.png"), p8, width = 10, height = 8, dpi = 150)
  
  p9 <- ggplot(results_res$importance_df, aes(x = reorder(variable, Overall), y = Overall)) +
    geom_col(fill = "purple") +
    coord_flip() +
    labs(title = "Importance des variables (permutation)", subtitle = diag_label,
         x = NULL, y = "Importance") +
    theme_minimal()
  ggsave(file.path(res_dir, "09_importance_variables.png"), p9, width = 7, height = 5, dpi = 150)
  
  write.csv(metrics_res, file.path(res_dir, "metrics_par_fold.csv"), row.names = FALSE)
  write.csv(results_res$tuning_summary, file.path(res_dir, "tuning_summary.csv"), row.names = FALSE)
  writeLines(sprintf("RMSE globale = %.4f\nR2 global = %.4f\nNombre de folds = %d",
                     global_rmse_res, global_r2_res, length(built$folds)),
             file.path(res_dir, "resume_global.txt"))
  
  cat("  ->", length(list.files(res_dir, pattern = "\\.png$")), "plots enregistrés dans", res_dir, "\n")
  
  invisible(list(results = results_res, metrics = metrics_res,
                 obs_pred = obs_pred_res, covariate_stats = covariate_stats_res))
}

# --- Boucle sur les résolutions spatiales ---
all_diagnostics_spatial <- imap(built_spatial_list, function(built, label) {
  generate_and_save_diagnostics(built, label, tuning_grid, output_root,
                                buffer_description = paste0("buffer ", built$buffer_km, " km, ", length(built$folds), " folds"))
})

# --- Boucle sur les résolutions temporelles ---
all_diagnostics_temporal <- imap(built_temporal_list, function(built, label) {
  generate_and_save_diagnostics(built, label, tuning_grid, output_root,
                                buffer_description = paste0(length(built$folds), " folds"))
})

cat("\nTous les diagnostics et modèles sont enregistrés sous :", normalizePath(output_root), "\n")

# =====================================================================
# 11. PRÉDICTION DU NASC SUR TOUTE LA GRILLE, POUR UNE DATE DONNÉE
#     -- pour TOUTES les résolutions entraînées (spatiales + temporelles)
# =====================================================================

library(ggplot2)

day_ds <- readRDS("F:/data_elise/ds_day_ftle_pig_fod/ds_ftle_pig_fod_20230126.rds")

grid_points_base <- expand.grid(lon = day_ds$lon, lat = day_ds$lat)

grid_points_base$ftle <- as.vector(day_ds$ftle)

for (p in names(day_ds$pig)) {
  grid_points_base[[p]] <- as.vector(day_ds$pig[[p]])
}

grid_points_base$total_chla      <- grid_points_base$Chla
grid_points_base$per_ratio_chla  <- grid_points_base$Per  / grid_points_base$Chla
grid_points_base$but_ratio_chla  <- grid_points_base$But  / grid_points_base$Chla
grid_points_base$fuco_ratio_chla <- grid_points_base$Fuco / grid_points_base$Chla
grid_points_base$hex_ratio_chla  <- grid_points_base$Hex  / grid_points_base$Chla
grid_points_base$allo_ratio_chla <- grid_points_base$Allo / grid_points_base$Chla
grid_points_base$zea_ratio_chla  <- grid_points_base$Zea  / grid_points_base$Chla
grid_points_base$chlb_ratio_chla <- grid_points_base$Chlb / grid_points_base$Chla

for (v in covariates_num) {
  bad <- !is.finite(grid_points_base[[v]])
  grid_points_base[[v]][bad] <- NA
}

grid_points_base$fod <- formatC(as.vector(day_ds$fod), width = 2)
grid_points_base$fod <- factor(grid_points_base$fod, levels = fod_levels)

grid_points_clean <- grid_points_base[stats::complete.cases(grid_points_base[, covariates_all]), ]

cat("Points valides :", nrow(grid_points_clean), "/", nrow(grid_points_base), "\n")

model_configs <- c(
  setNames(
    file.path(output_root, names(built_spatial_list), "rf_model_final.rds"),
    names(built_spatial_list)
  ),
  setNames(
    file.path(output_root, names(built_temporal_list), "rf_model_final.rds"),
    names(built_temporal_list)
  )
)

print(model_configs)

prediction_dir <- file.path(output_root, "predictions")
dir.create(prediction_dir, showWarnings = FALSE, recursive = TRUE)

all_predictions <- list()

for (config_name in names(model_configs)) {
  
  model_path <- model_configs[[config_name]]
  
  if (!file.exists(model_path)) {
    cat("  [!] Modèle introuvable, ignoré :", model_path, "\n")
    next
  }
  
  cat("\n==== Prédiction avec le modèle :", config_name, "====\n")
  
  model_final <- readRDS(model_path)
  
  pred_df <- grid_points_clean
  pred_df$NASC_pred <- predict(model_final, data = pred_df)$predictions
  pred_df$config <- config_name
  
  all_predictions[[config_name]] <- pred_df
  
  p_pred <- ggplot(pred_df, aes(x = lon, y = lat, fill = NASC_pred)) +
    geom_raster() +
    scale_fill_viridis_c() +
    coord_quickmap() +
    theme_bw() +
    labs(
      title = paste("NASC prédit -", format(day_ds$date, "%Y-%m-%d")),
      subtitle = paste0("Modèle : ", config_name, " -- ",
                        nrow(pred_df), " / ", nrow(grid_points_base), " pixels prédits"),
      x = "Longitude", y = "Latitude",
      fill = "log10(NASC)"
    )
  
  print(p_pred)
  
  ggsave(
    file.path(prediction_dir, paste0("NASC_pred_", format(day_ds$date, "%Y%m%d"), "_", config_name, ".png")),
    p_pred, width = 8, height = 6, dpi = 150
  )
  
  cat("  -> Carte sauvegardée :",
      file.path(prediction_dir, paste0("NASC_pred_", format(day_ds$date, "%Y%m%d"), "_", config_name, ".png")), "\n")
}

predictions_all <- dplyr::bind_rows(all_predictions)

shared_limits <- range(predictions_all$NASC_pred, na.rm = TRUE)

p_compare <- ggplot(predictions_all, aes(x = lon, y = lat, fill = NASC_pred)) +
  geom_raster() +
  scale_fill_viridis_c(limits = shared_limits) +
  coord_quickmap() +
  facet_wrap(~config) +
  theme_bw() +
  labs(
    title = paste("NASC prédit -", format(day_ds$date, "%Y-%m-%d"), "-- comparaison des modèles"),
    x = "Longitude", y = "Latitude",
    fill = "log10(NASC)"
  )

print(p_compare)

ggsave(
  file.path(prediction_dir, paste0("NASC_pred_", format(day_ds$date, "%Y%m%d"), "_comparaison_toutes_resolutions.png")),
  p_compare, width = 14, height = 10, dpi = 150
)

cat("\nCarte comparative sauvegardée dans :", prediction_dir, "\n")