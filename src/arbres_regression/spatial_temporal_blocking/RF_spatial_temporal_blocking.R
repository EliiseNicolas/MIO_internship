# =====================================================================
# Random Forest : CV bloquée -- DEUX SCHÉMAS INDÉPENDANTS
# =====================================================================
# 1) Blocage TEMPOREL : blocs de 20 jours (chaque bloc = un fold), avec un
#    tampon de 5 jours autour du bloc de test (aucune donnée d'entraînement
#    à moins de 5 jours de n'importe quelle date du fold de test).
# 2) Blocage SPATIAL : grille en km (résolutions lon/lat testées), avec un
#    tampon de 200 km autour du bloc de test (aucune donnée d'entraînement
#    à moins de 200 km de n'importe quel point du fold de test).
#
# Ces deux schémas sont VOLONTAIREMENT séparés (et non combinés comme dans
# une version précédente) : l'objectif est de pouvoir comparer les deux
# (performance, extrapolation, stabilité de l'importance des variables),
# pas de construire un blocage spatio-temporel joint.
#
# install.packages(c("dplyr","tidyr","purrr","ggplot2","ranger","FNN","tibble","caret","CAST"))

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(ranger)
library(FNN)
library(tibble)
library(caret)
library(CAST)

set.seed(42)

# ---------------------------------------------------------------------
# 0. DONNÉES (pipeline de nettoyage inchangé)
# ---------------------------------------------------------------------
freq <- 18
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

# R² "prédictif" (= Nash-Sutcliffe Efficiency), référencé sur la moyenne du
# TRAIN (pas du test) : mesure si le modèle fait mieux que prédire la
# moyenne apprise sur le train, ce qui est le bon repère quand on évalue de
# l'extrapolation (utiliser la moyenne du test serait optimiste, car cette
# info ne serait pas disponible en usage réel).
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

# Gap temporel minimal (jours), version vectorisée (findInterval) -- utilisé
# pour le tampon du blocage TEMPOREL.
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
                                 min_block_n = 50, max_folds = 10, seed = 42) {
  data_blocked <- assign_temporal_block(data, block_days)
  blocks_count <- data_blocked %>% count(temporal_block)
  blocks_used  <- blocks_count %>% filter(n >= min_block_n) %>% pull(temporal_block)
  set.seed(seed)
  if (length(blocks_used) > max_folds) blocks_used <- sample(blocks_used, max_folds)
  folds <- map(blocks_used, ~ make_temporal_fold(data_blocked, .x, buffer_days))
  names(folds) <- paste0("t_", blocks_used)
  
  # Garde-fou : un fold avec train ou test vide (ex. domaine trop petit pour
  # le nombre de blocs demandé) ferait planter ranger() plus loin -- on
  # l'écarte ici avec un message plutôt que de laisser planter le script.
  valid <- map_lgl(folds, ~ length(.x$train) > 0 && length(.x$test) > 0)
  if (any(!valid)) {
    cat("  [!] ", sum(!valid), "bloc(s) temporel(s) ignoré(s) : train ou test vide après buffer\n")
  }
  folds <- folds[valid]
  
  list(data = data_blocked, folds = folds,
       n_blocks_total = nrow(blocks_count), n_blocks_used = length(folds))
}

# ---------------------------------------------------------------------
# 3B. BLOCAGE SPATIAL : grille en km (résolutions testées) + tampon en km
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
  # distance de chaque point candidat au point de test le plus proche
  # (FNN::knnx.dist, vectorisé -- indispensable vu la taille des données)
  dist_to_test <- FNN::knnx.dist(data = test_mat, query = candidate_mat, k = 1)[, 1]
  train_idx <- candidate_train_idx[dist_to_test > buffer_km]
  list(train = train_idx, test = test_idx)
}

build_spatial_folds <- function(data, cellsize_x, cellsize_y, buffer_km,
                                min_block_n = 50, max_folds = 10, seed = 42) {
  data_blocked <- assign_spatial_block(data, cellsize_x, cellsize_y)
  blocks_count <- data_blocked %>% count(spatial_block)
  blocks_used  <- blocks_count %>% filter(n >= min_block_n) %>% pull(spatial_block)
  set.seed(seed)
  if (length(blocks_used) > max_folds) blocks_used <- sample(blocks_used, max_folds)
  folds <- map(blocks_used, ~ make_spatial_fold_buffered(data_blocked, .x, buffer_km))
  names(folds) <- paste0("s_", blocks_used)
  
  # Garde-fou : un fold avec train ou test vide (ex. le domaine tient tout
  # entier dans une seule cellule à cette résolution -> aucun autre bloc
  # disponible comme train) ferait planter ranger() plus loin -- on l'écarte
  # ici avec un message plutôt que de laisser planter le script.
  valid <- map_lgl(folds, ~ length(.x$train) > 0 && length(.x$test) > 0)
  if (any(!valid)) {
    cat("  [!] ", sum(!valid), "bloc(s) spatial/spatiaux ignoré(s) : train ou test vide après buffer\n")
  }
  folds <- folds[valid]
  
  list(data = data_blocked, folds = folds,
       n_blocks_total = nrow(blocks_count), n_blocks_used = length(folds))
}

# ---------------------------------------------------------------------
# 4. CONSTRUCTION DES DEUX SCHÉMAS + COMPARAISON BASELINE
# ---------------------------------------------------------------------
temporal_block_days  <- 10
temporal_buffer_days <- 5
spatial_buffer_km    <- 5

resolutions_to_test <- list(
  c(lon_km = 60, lat_km = 60)
  # c(lon_km = 1000, lat_km = 700),
  # c(lon_km = 500,  lat_km = 500),
  # c(lon_km = 200,  lat_km = 200)
)

built_temporal <- build_temporal_folds(df, temporal_block_days, temporal_buffer_days)

built_spatial_list <- map(resolutions_to_test, ~
                            build_spatial_folds(df, .x["lon_km"], .x["lat_km"], spatial_buffer_km))
names(built_spatial_list) <- map_chr(resolutions_to_test, ~ paste0(.x["lon_km"], "x", .x["lat_km"], "km"))

all_built <- c(list(temporal_20j_buffer5j = built_temporal), built_spatial_list)

baseline_params <- list(mtry = 3, min.node.size = 5, num.trees = 300)

# NB : si le buffer est grand par rapport à la taille du bloc (ex. 200km de
# buffer avec des cellules de 200x200km), il peut ne rester quasiment plus
# de train dans certains folds -- vérifie mean_n_train ci-dessous.
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

cat("\n=== Comparaison blocage temporel vs blocage spatial (toutes résolutions) ===\n")
print(comparison_summary)

ggplot(comparison_summary, aes(x = config, y = mean_rmse_test)) +
  geom_col(fill = "steelblue") +
  geom_errorbar(aes(ymin = mean_rmse_test - sd_rmse_test, ymax = mean_rmse_test + sd_rmse_test), width = 0.2) +
  labs(title = "RMSE test par schéma de blocage",
       subtitle = "Barres d'erreur = écart-type inter-fold",
       x = NULL, y = "RMSE test (moyenne inter-fold)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# R² test (référencé sur la moyenne du train) : normalisé par la variance,
# donc comparable entre schémas/folds même si leur variance de NASC diffère.
# Une valeur proche de 0 (ou négative) = le modèle ne fait pas mieux (ou
# moins bien) que prédire la moyenne du train -> extrapolation qui échoue.
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
# 5. PIPELINE COMPLET (tuning, courbe d'apprentissage, extrapolation, AOA,
#    importance des variables) -- réutilisable pour n'importe quel schéma
# ---------------------------------------------------------------------
tuning_grid <- expand.grid(mtry = c(2, 3, 4), min.node.size = c(5, 10), num.trees = 300)

run_pipeline <- function(built, label, tuning_grid, fractions = seq(0.2, 1, by = 0.2)) {
  folds <- built$folds
  data  <- built$data
  
  # --- tuning ---
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
  
  # --- courbe d'apprentissage ---
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
  
  # --- extrapolation DIY (covariables numériques uniquement) ---
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
  
  # --- modèle caret + AOA (CAST) + importance des variables ---
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
  
  AOA_result <- aoa(newdata = as.data.frame(data[, covariates_all]), model = model_caret)
  
  importance_df <- as.data.frame(caret::varImp(model_caret)$importance) %>%
    tibble::rownames_to_column("variable") %>%
    mutate(config = label)
  
  list(label = label, tuning_summary = tuning_summary, best_params = best_params,
       learning_curve_summary = learning_curve_summary,
       extrapolation_results = extrapolation_results,
       model_caret = model_caret, AOA_result = AOA_result, importance_df = importance_df)
}

# ---------------------------------------------------------------------
# 6. APPLICATION : bloc temporel vs UNE résolution spatiale choisie
# ---------------------------------------------------------------------
# Ajuste ce label selon ce que montre la section 4 (ex. "500x500km")
SELECTED_SPATIAL_LABEL <- "60x60km"

results_temporal <- run_pipeline(built_temporal, "temporel_20j_buffer5j", tuning_grid)
results_spatial  <- run_pipeline(built_spatial_list[[SELECTED_SPATIAL_LABEL]],
                                 paste0("spatial_", SELECTED_SPATIAL_LABEL, "_buffer200km"),
                                 tuning_grid)

cat("\nMeilleurs hyperparamètres -- bloc temporel :\n"); print(results_temporal$best_params)
cat("Meilleurs hyperparamètres -- bloc spatial :\n");   print(results_spatial$best_params)

# ---------------------------------------------------------------------
# 7. COMPARAISON DE LA STABILITÉ DE L'IMPORTANCE DES VARIABLES
# ---------------------------------------------------------------------
# Une variable importante dans LES DEUX schémas = lien structurel probable.
# Une variable qui s'effondre dans un des deux = probable proxy (spatial ou
# temporel) plutôt qu'une vraie relation mécanistique (point 7 du protocole).
importance_comparison <- bind_rows(results_temporal$importance_df, results_spatial$importance_df)

ggplot(importance_comparison, aes(x = variable, y = Overall, fill = config)) +
  geom_col(position = "dodge") +
  coord_flip() +
  labs(title = "Stabilité de l'importance des variables : bloc temporel vs bloc spatial",
       x = NULL, y = "Importance (permutation)", fill = NULL) +
  theme_minimal()

# ---------------------------------------------------------------------
# 8. AOA -- résumé pour les deux schémas
# ---------------------------------------------------------------------
cat("\nProportion hors AOA -- bloc temporel :",
    mean(results_temporal$AOA_result$AOA == 0), "\n")
cat("Proportion hors AOA -- bloc spatial   :",
    mean(results_spatial$AOA_result$AOA == 0), "\n")

# =====================================================================
# 9. DIAGNOSTICS COMPLETS PAR FOLD (obs/pred, R², variance intra-fold,
#    distances test->train, stats des covariables par fold)
# =====================================================================
# Réutilise les objets déjà construits ci-dessus -- aucune librairie
# supplémentaire à charger, rien à recharger.

# --- 9.0. Choix du schéma à diagnostiquer ---
# Change les DEUX lignes suivantes ensemble (built + results_selected) --
# ne jamais changer l'une sans l'autre, sinon les graphiques mélangeront
# les folds d'un schéma avec les résultats (best_params, learning curve)
# d'un autre schéma.
built           <- built_spatial_list[[SELECTED_SPATIAL_LABEL]]   # ou built_temporal
results_selected <- results_spatial                                # ou results_temporal

folds  <- built$folds
data   <- built$data
params <- results_selected$best_params

# Label dynamique injecté dans tous les titres ci-dessous, pour que chaque
# graphique dise explicitement de quel schéma de blocage il provient.
diagnostic_label <- paste0("blocage spatial ", SELECTED_SPATIAL_LABEL,
                           " (buffer ", spatial_buffer_km, " km)")
# Si tu utilises built_temporal à la place, remplace par :
# diagnostic_label <- paste0("blocage temporel ", temporal_block_days, "j",
#                             " (buffer ", temporal_buffer_days, "j)")

# --- 9.1. Fonction de diagnostic par fold (un seul entraînement par fold) ---
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
  
  # distance GÉOGRAPHIQUE : chaque point de test -> son train le plus proche
  geo_dist <- FNN::knnx.dist(
    data  = as.matrix(train_df[, c("x_km", "y_km")]),
    query = as.matrix(test_df[,  c("x_km", "y_km")]),
    k = 1
  )[, 1]
  metrics$mean_geo_dist_km <- mean(geo_dist)
  
  # distance dans l'ESPACE DES COVARIABLES numériques (indice d'extrapolation)
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
  
  # moyenne / écart-type de chaque covariable numérique, train ET test
  covariate_stats <- bind_rows(
    tibble(set = "test",  variable = covs_num,
           mean = colMeans(test_num,  na.rm = TRUE), sd = apply(test_num,  2, sd, na.rm = TRUE)),
    tibble(set = "train", variable = covs_num,
           mean = colMeans(train_num, na.rm = TRUE), sd = apply(train_num, 2, sd, na.rm = TRUE))
  )
  
  list(obs_pred = obs_pred, metrics = metrics, covariate_stats = covariate_stats)
}
# --- Distribution FOD (catégorielle) train vs test, par fold ---
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

# --- Valeurs brutes FTLE train vs test, par fold (pour distribution complète) ---
compute_ftle_values <- function(fold, data) {
  bind_rows(
    tibble(ftle = data$ftle[fold$train], set = "train"),
    tibble(ftle = data$ftle[fold$test],  set = "test")
  ) %>%
    filter(is.finite(ftle))
}

# --- 9.2. Application à tous les folds ---
all_results <- imap(folds, ~ compute_fold_diagnostics(.x, data, params))
fod_dist_all  <- imap_dfr(folds, ~ compute_fod_distribution(.x, data) %>% mutate(fold_id = .y))
ftle_vals_all <- imap_dfr(folds, ~ compute_ftle_values(.x, data)      %>% mutate(fold_id = .y))

obs_pred_all        <- imap_dfr(all_results, ~ mutate(.x$obs_pred, fold_id = .y))
metrics_all         <- imap_dfr(all_results, ~ mutate(.x$metrics, fold_id = .y))
covariate_stats_all <- imap_dfr(all_results, ~ mutate(.x$covariate_stats, fold_id = .y))

cat("\n=== Métriques par fold --", diagnostic_label, "===\n")
print(metrics_all)

# --- 9.3. NASC observé vs prédit ---
# RMSE et R² globaux = calculés sur l'ensemble des prédictions de TOUS les
# folds combinés (pas la moyenne des métriques par fold -- une valeur
# poolée). Le R² global utilise ici la moyenne du test poolé comme
# référence (convention R² classique) ; c'est différent du r2_test par
# fold (section 9.5) qui utilise la moyenne du train de CE fold pour rester
# honnête vis-à-vis de l'extrapolation. Les deux sont légitimes mais
# répondent à des questions différentes -- garde cette distinction en tête
# si les valeurs ne concordent pas exactement.
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

# --- 9.4. Courbe d'apprentissage (déjà calculée en section 5/6, tracée ici) ---
ggplot(results_selected$learning_curve_summary, aes(x = fraction)) +
  geom_line(aes(y = mean_rmse_train, color = "Train")) +
  geom_point(aes(y = mean_rmse_train, color = "Train")) +
  geom_line(aes(y = mean_rmse_test, color = "Test")) +
  geom_point(aes(y = mean_rmse_test, color = "Test")) +
  labs(title = "Courbe d'apprentissage", subtitle = diagnostic_label,
       x = "Fraction du train utilisée",
       y = "RMSE", color = NULL) +
  theme_minimal()

# --- 9.5. R² par fold ---
ggplot(metrics_all, aes(x = fold_id, y = r2_test)) +
  geom_col(fill = "darkgreen") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "R² par fold (référencé sur la moyenne du train)",
       subtitle = diagnostic_label,
       x = "Fold", y = "R² test") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# --- 9.6. Variance intra-fold (dispersion du NASC observé dans chaque fold test) ---
ggplot(metrics_all, aes(x = fold_id, y = var_intra_fold_test)) +
  geom_col(fill = "steelblue") +
  labs(title = "Variance du NASC observé, intra-fold (test)",
       subtitle = diagnostic_label,
       x = "Fold", y = "Variance (log10 NASC)") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# --- 9.7. Distance test -> train par fold (géographique + covariables) ---
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

# --- 9.8. Moyenne / écart-type de chaque covariable, par fold (train vs test) ---
ggplot(covariate_stats_all, aes(x = fold_id, y = mean, color = set)) +
  geom_point(position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2,
                position = position_dodge(width = 0.4)) +
  facet_wrap(~variable, scales = "free_y") +
  labs(title = "Moyenne ± écart-type de chaque covariable, par fold (train vs test)",
       subtitle = diagnostic_label,
       x = "Fold", y = NULL, color = NULL) +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# =====================================================================
# 10. CARTE DES POINTS NASC COLORÉS PAR FOLD
# =====================================================================
# Attribue à chaque point le fold dont il est le TEST ; NA si son bloc
# n'a pas été retenu comme fold (ex. bloc trop petit, cf. min_block_n).

# --- 10.1. Blocage SPATIAL ---
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

# --- 10.2. Blocage TEMPOREL ---
# Utile pour comparer : contrairement au spatial, les folds temporels ne
# forment PAS des zones géographiques compactes (un même bloc de 20 jours
# peut regrouper des points dispersés dans tout le domaine).
data_map_t <- built_temporal$data
data_map_t$row_id <- seq_len(nrow(data_map_t))

fold_lookup_t <- imap_dfr(built_temporal$folds, ~ tibble(row_id = .x$test, fold_id = .y))
data_map_t <- data_map_t %>%
  left_join(fold_lookup_t, by = "row_id") %>%
  mutate(fold_id = ifelse(is.na(fold_id), "Non utilisé (bloc trop petit)", fold_id))

ggplot(data_map_t, aes(x = lon, y = lat, color = fold_id)) +
  geom_point(size = 0.8, alpha = 0.7) +
  coord_quickmap() +
  labs(title = "Points NASC colorés par fold -- blocage temporel (20j)",
       x = "Longitude", y = "Latitude", color = "Fold") +
  theme_minimal()

# --- 10.3. Effet du buffer sur TOUS les folds spatiaux (facet_wrap) ---
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
       subtitle = paste0("Résolution ", SELECTED_SPATIAL_LABEL, " -- buffer ", spatial_buffer_km, " km"),
       x = "Longitude", y = "Latitude", color = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom")

# Si tu préfères un grand graphique par fold plutôt que des vignettes côte
# à côte, décommente cette boucle (elle affichera chaque fold séparément) :
#
# iwalk(built_map_s$folds, function(fold, fid) {
#   df_status <- data_map_s %>% select(lon, lat)
#   df_status$status <- "Exclu par le buffer"
#   df_status$status[fold$train] <- "Train conservé"
#   df_status$status[fold$test]  <- "Test (bloc isolé)"
#   p <- ggplot(df_status, aes(x = lon, y = lat, color = status)) +
#     geom_point(size = 0.8, alpha = 0.7) +
#     scale_color_manual(values = c("Test (bloc isolé)" = "red",
#                                    "Train conservé"    = "steelblue",
#                                    "Exclu par le buffer" = "grey80")) +
#     coord_quickmap() +
#     labs(title = paste0("Effet du buffer spatial -- fold ", fid),
#          subtitle = paste0("Résolution ", SELECTED_SPATIAL_LABEL, " -- buffer ", spatial_buffer_km, " km"),
#          x = "Longitude", y = "Latitude", color = NULL) +
#     theme_minimal()
#   print(p)
# })

# =====================================================================
# 11. GÉNÉRATION + SAUVEGARDE DE TOUS LES PLOTS, POUR TOUTES LES RÉSOLUTIONS
# =====================================================================
# Un dossier par résolution (ex. outputs_cv_diagnostics/500x500km/), avec
# chaque plot en .png numéroté + les tableaux de métriques en .csv.
# ATTENTION : ceci relance run_pipeline() (tuning + courbe d'apprentissage +
# AOA) pour CHAQUE résolution -- c'est le calcul le plus coûteux du script,
# multiplié par 4 ici. Prévois du temps si le jeu de données est gros.

output_root <- "outputs_cv_diagnostics"
dir.create(output_root, showWarnings = FALSE)

generate_and_save_diagnostics <- function(built, label, tuning_grid, output_root,
                                          buffer_description = "") {
  
  res_dir <- file.path(output_root, label)
  dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)
  diag_label <- paste0(label, if (nchar(buffer_description)) paste0(" (", buffer_description, ")") else "")
  
  cat("\n==== Résolution :", label, "-- dossier :", res_dir, "====\n")
  
  # --- pipeline complet (tuning, learning curve, extrapolation, AOA, importance) ---
  results_res <- run_pipeline(built, label, tuning_grid)
  
  # --- diagnostics par fold ---
  all_fold_results <- imap(built$folds, ~ compute_fold_diagnostics(.x, built$data, results_res$best_params))
  obs_pred_res        <- imap_dfr(all_fold_results, ~ mutate(.x$obs_pred, fold_id = .y))
  metrics_res         <- imap_dfr(all_fold_results, ~ mutate(.x$metrics, fold_id = .y))
  covariate_stats_res <- imap_dfr(all_fold_results, ~ mutate(.x$covariate_stats, fold_id = .y))
  
  global_rmse_res <- sqrt(mean((obs_pred_res$obs - obs_pred_res$pred)^2))
  global_r2_res   <- 1 - sum((obs_pred_res$obs - obs_pred_res$pred)^2) /
    sum((obs_pred_res$obs - mean(obs_pred_res$obs))^2)
  
  # --- 1. obs vs pred ---
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
  
  # --- 2. courbe d'apprentissage ---
  p2 <- ggplot(results_res$learning_curve_summary, aes(x = fraction)) +
    geom_line(aes(y = mean_rmse_train, color = "Train")) +
    geom_point(aes(y = mean_rmse_train, color = "Train")) +
    geom_line(aes(y = mean_rmse_test, color = "Test")) +
    geom_point(aes(y = mean_rmse_test, color = "Test")) +
    labs(title = "Courbe d'apprentissage", subtitle = diag_label,
         x = "Fraction du train utilisée", y = "RMSE", color = NULL) +
    theme_minimal()
  ggsave(file.path(res_dir, "02_learning_curve.png"), p2, width = 7, height = 5, dpi = 150)
  
  # --- 3. R² par fold ---
  p3 <- ggplot(metrics_res, aes(x = fold_id, y = r2_test)) +
    geom_col(fill = "darkgreen") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    labs(title = "R² par fold (référencé sur la moyenne du train)", subtitle = diag_label,
         x = "Fold", y = "R² test") +
    theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(res_dir, "03_r2_par_fold.png"), p3, width = 7, height = 5, dpi = 150)
  
  # --- 4. variance intra-fold ---
  p4 <- ggplot(metrics_res, aes(x = fold_id, y = var_intra_fold_test)) +
    geom_col(fill = "steelblue") +
    labs(title = "Variance du NASC observé, intra-fold (test)", subtitle = diag_label,
         x = "Fold", y = "Variance (log10 NASC)") +
    theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(res_dir, "04_variance_intra_fold.png"), p4, width = 7, height = 5, dpi = 150)
  
  # --- 5. distances test -> train ---
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
  
  # --- 6. moyenne/écart-type des covariables par fold ---
  p6 <- ggplot(covariate_stats_res, aes(x = fold_id, y = mean, color = set)) +
    geom_point(position = position_dodge(width = 0.4)) +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2, position = position_dodge(width = 0.4)) +
    facet_wrap(~variable, scales = "free_y") +
    labs(title = "Moyenne ± écart-type de chaque covariable, par fold (train vs test)", subtitle = diag_label,
         x = "Fold", y = NULL, color = NULL) +
    theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(res_dir, "06_covariable_stats.png"), p6, width = 10, height = 8, dpi = 150)
  
  # --- 7. carte des points colorés par fold ---
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
  
  # --- 8. effet du buffer, tous les folds ---
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
  
  # --- 9. importance des variables (issue du modèle caret de run_pipeline) ---
  p9 <- ggplot(results_res$importance_df, aes(x = reorder(variable, Overall), y = Overall)) +
    geom_col(fill = "purple") +
    coord_flip() +
    labs(title = "Importance des variables (permutation)", subtitle = diag_label,
         x = NULL, y = "Importance") +
    theme_minimal()
  ggsave(file.path(res_dir, "09_importance_variables.png"), p9, width = 7, height = 5, dpi = 150)
  
  # --- tableaux associés ---
  write.csv(metrics_res, file.path(res_dir, "metrics_par_fold.csv"), row.names = FALSE)
  write.csv(results_res$tuning_summary, file.path(res_dir, "tuning_summary.csv"), row.names = FALSE)
  writeLines(sprintf("RMSE globale = %.4f\nR2 global = %.4f\nProportion hors AOA = %.4f",
                     global_rmse_res, global_r2_res, mean(results_res$AOA_result$AOA == 0)),
             file.path(res_dir, "resume_global.txt"))
  
  cat("  ->", length(list.files(res_dir, pattern = "\\.png$")), "plots enregistrés dans", res_dir, "\n")
  
  invisible(list(results = results_res, metrics = metrics_res,
                 obs_pred = obs_pred_res, covariate_stats = covariate_stats_res))
}

# --- Boucle sur les 4 résolutions spatiales ---
all_diagnostics_spatial <- imap(built_spatial_list, function(built, label) {
  generate_and_save_diagnostics(built, label, tuning_grid, output_root,
                                buffer_description = paste0("buffer ", spatial_buffer_km, " km"))
})

# --- Le schéma temporel, dans son propre dossier ---
diagnostics_temporal <- generate_and_save_diagnostics(
  built_temporal, "temporel_20j", tuning_grid, output_root,
  buffer_description = paste0("buffer ", temporal_buffer_days, " j")
)

cat("\nTous les diagnostics sont enregistrés sous :", normalizePath(output_root), "\n")

  
