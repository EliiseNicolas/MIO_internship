library(rpart)
library(rpart.plot)
library(ggplot2)
library(randomForest)
library(mgcv)
library(dplyr)
library(tidyr)
library(purrr)
library(FNN)

freq <- 120
diurnal_period <- 3 # 3 : day, 1: night
dp <- "day"
n_CV <- 10
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

fod_levels <- levels(df$fod)  # sauvegardé pour la prédiction future

cat("Nombre d'observations après filtrage :", nrow(df), "\n") # 47000

# NB : lat/lon/time CONSERVÉS (contrairement à la version précédente) --
# nécessaires pour les distances géographiques et la carte des résidus.
str(df)

# =====================================================================
# CV à 10 folds -- diagnostics complets par fold
# =====================================================================

set.seed(123)
n <- nrow(df)
fold_id <- sample(rep(seq_len(n_CV), length.out = n))

rmse_fn <- function(obs, pred) sqrt(mean((obs - pred)^2, na.rm = TRUE))

compute_fold_diagnostics <- function(k, data, fold_id, covs_num, covs_all, response) {
  
  train_idx <- which(fold_id != k)
  test_idx  <- which(fold_id == k)
  
  train_df <- data[train_idx, ]
  test_df  <- data[test_idx, ]
  
  model_k <- randomForest(
    stats::as.formula(paste(response, "~", paste(covs_all, collapse = " + "))),
    data = train_df[, c(response, covs_all)],
    ntree = 300,
    importance = TRUE
  )
  
  pred_test  <- predict(model_k, newdata = test_df)
  pred_train <- predict(model_k, newdata = train_df)
  
  obs_test  <- test_df[[response]]
  obs_train <- train_df[[response]]
  
  train_mean <- mean(obs_train)
  
  residuals <- obs_test - pred_test
  
  metrics <- tibble(
    fold                = k,
    n_train             = nrow(train_df),
    n_test              = nrow(test_df),
    rmse_test           = rmse_fn(obs_test, pred_test),
    rmse_train          = rmse_fn(obs_train, pred_train),
    r2_test             = 1 - sum((obs_test - pred_test)^2) / sum((obs_test - train_mean)^2),
    r2_train            = 1 - sum((obs_train - pred_train)^2) / sum((obs_train - train_mean)^2),
    var_intra_fold_test = var(obs_test),
    sd_intra_fold_test  = sd(obs_test),
    mean_residual       = mean(residuals),
    sd_residual         = sd(residuals)
  )
  
  # distance géographique : chaque point de test -> son train le plus proche
  geo_dist <- FNN::knnx.dist(
    data  = as.matrix(train_df[, c("lon", "lat")]),
    query = as.matrix(test_df[,  c("lon", "lat")]),
    k = 1
  )[, 1]
  metrics$mean_geo_dist_km <- mean(geo_dist) * 111  # approx degrés -> km
  
  # distance euclidienne dans l'espace des covariables numériques
  # (standardisées sur le train)
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
  
  # stats des covariables numériques, train vs test
  covariate_stats <- bind_rows(
    tibble(set = "test",  variable = covs_num,
           mean = colMeans(test_num,  na.rm = TRUE), sd = apply(test_num,  2, sd, na.rm = TRUE)),
    tibble(set = "train", variable = covs_num,
           mean = colMeans(train_num, na.rm = TRUE), sd = apply(train_num, 2, sd, na.rm = TRUE))
  ) %>% mutate(fold = k)
  
  obs_pred <- tibble(
    obs = obs_test, pred = pred_test, residual = residuals,
    lon = test_df$lon, lat = test_df$lat, fold = k
  )
  
  list(metrics = metrics, covariate_stats = covariate_stats,
       obs_pred = obs_pred, model = model_k)
}

cat("\n=== Entraînement des", n_CV, "folds ===\n")
all_fold_results <- map(seq_len(n_CV), ~ compute_fold_diagnostics(
  .x, df, fold_id, covariates_num, covariates_all, response_var
))

metrics_all         <- map_dfr(all_fold_results, "metrics")
covariate_stats_all <- map_dfr(all_fold_results, "covariate_stats")
obs_pred_all        <- map_dfr(all_fold_results, "obs_pred")

cat("\n=== Métriques par fold ===\n")
print(metrics_all)

cat("\n--- Moyenne / écart-type inter-fold ---\n")
cat("RMSE test :", round(mean(metrics_all$rmse_test), 4), "(sd =", round(sd(metrics_all$rmse_test), 4), ")\n")
cat("R2 test   :", round(mean(metrics_all$r2_test), 4), "(sd =", round(sd(metrics_all$r2_test), 4), ")\n")

# ---------------------------------------------------------------------
# 1. NASC observé vs prédit (toutes prédictions out-of-fold poolées)
# ---------------------------------------------------------------------
global_rmse <- rmse_fn(obs_pred_all$obs, obs_pred_all$pred)
global_r2   <- 1 - sum((obs_pred_all$obs - obs_pred_all$pred)^2) /
  sum((obs_pred_all$obs - mean(obs_pred_all$obs))^2)

ggplot(obs_pred_all, aes(x = obs, y = pred, color = factor(fold))) +
  geom_point(alpha = 0.4, size = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  annotate("label",
           x = min(obs_pred_all$obs), y = max(obs_pred_all$pred),
           hjust = 0, vjust = 1,
           label = sprintf("RMSE = %.3f\nR² = %.3f", global_rmse, global_r2),
           fill = "white", alpha = 0.85, size = 3.5) +
  labs(
    title = "NASC observé vs prédit -- 10-fold CV (out-of-fold)",
    x = "NASC observé (log10)", y = "NASC prédit (log10)", color = "Fold"
  ) +
  theme_bw()

# ---------------------------------------------------------------------
# 2. RMSE et R² par fold
# ---------------------------------------------------------------------
ggplot(metrics_all, aes(x = factor(fold), y = rmse_test)) +
  geom_col(fill = "steelblue") +
  labs(title = "RMSE par fold", x = "Fold", y = "RMSE") +
  theme_bw()

ggplot(metrics_all, aes(x = factor(fold), y = r2_test)) +
  geom_col(fill = "darkgreen") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "R² par fold", x = "Fold", y = "R²") +
  theme_bw()

# ---------------------------------------------------------------------
# 3. Variabilité intra-fold (variance du NASC observé dans chaque fold test)
# ---------------------------------------------------------------------
ggplot(metrics_all, aes(x = factor(fold), y = var_intra_fold_test)) +
  geom_col(fill = "orange") +
  labs(title = "Variance du NASC observé, intra-fold (test)",
       x = "Fold", y = "Variance (log10 NASC)") +
  theme_bw()

# ---------------------------------------------------------------------
# 4. Distance géographique et distance dans l'espace des covariables,
#    par fold
# ---------------------------------------------------------------------
dist_long <- metrics_all %>%
  select(fold, mean_geo_dist_km, mean_covariate_dist) %>%
  pivot_longer(-fold, names_to = "metric", values_to = "value")

ggplot(dist_long, aes(x = factor(fold), y = value, fill = metric)) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~metric, scales = "free_y", labeller = as_labeller(c(
    mean_geo_dist_km    = "Distance géographique moyenne (km)",
    mean_covariate_dist = "Distance covariables test->train (standardisée)"
  ))) +
  labs(title = "Distance test -> train, par fold", x = "Fold", y = NULL) +
  theme_bw()

# ---------------------------------------------------------------------
# 5. Moyenne / écart-type de chaque covariable, train vs test, par fold
# ---------------------------------------------------------------------
ggplot(covariate_stats_all, aes(x = factor(fold), y = mean, color = set)) +
  geom_point(position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2,
                position = position_dodge(width = 0.4)) +
  facet_wrap(~variable, scales = "free_y") +
  labs(title = "Moyenne ± écart-type de chaque covariable, par fold (train vs test)",
       x = "Fold", y = NULL, color = NULL) +
  theme_bw()

# ---------------------------------------------------------------------
# 6. CARTE DES RÉSIDUS SPATIALISÉE
#    -- permet de voir si le modèle est plus performant dans certaines
#       zones géographiques que d'autres
# ---------------------------------------------------------------------
ggplot(obs_pred_all, aes(x = lon, y = lat, color = residual)) +
  geom_point(size = 1.2, alpha = 0.7) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  coord_quickmap() +
  labs(
    title = "Résidus spatialisés (observé - prédit), out-of-fold",
    subtitle = "Bleu = sous-estimation, Rouge = surestimation",
    x = "Longitude", y = "Latitude", color = "Résidu\n(log10 NASC)"
  ) +
  theme_bw()

# Version alternative : valeur absolue des résidus (juste l'ampleur de
# l'erreur, sans distinguer sur/sous-estimation) -- utile pour repérer
# directement les "zones difficiles" pour le modèle
ggplot(obs_pred_all, aes(x = lon, y = lat, color = abs(residual))) +
  geom_point(size = 1.2, alpha = 0.7) +
  scale_color_viridis_c(option = "magma", direction = -1) +
  coord_quickmap() +
  labs(
    title = "Amplitude des erreurs de prédiction, spatialisée (out-of-fold)",
    x = "Longitude", y = "Latitude", color = "|Résidu|\n(log10 NASC)"
  ) +
  theme_bw()

# ---------------------------------------------------------------------
# 7. Modèle final : entraîné sur TOUTES les données, pour la prédiction
#    sur grille (section suivante)
# ---------------------------------------------------------------------
model <- randomForest(
  NASC ~ .,
  data = df[, c(response_var, covariates_all)],
  ntree = 300,
  importance = TRUE
)


# ---- Sauvegarde du modèle et des niveaux fod (nécessaires en prédiction) ----
model_output_dir <- "outputs_rf_simple"
dir.create(model_output_dir, showWarnings = FALSE, recursive = TRUE)

saveRDS(model, file.path(model_output_dir, paste0("rf_model_final_", freq, "kHz.rds")))
saveRDS(fod_levels, file.path(model_output_dir, paste0("fod_levels_", freq, "kHz.rds")))

cat("Modèle sauvegardé :", file.path(model_output_dir, paste0("rf_model_final_", freq, "kHz.rds")), "\n")


importance_df <- data.frame(
  variable = rownames(importance(model)),
  importance = importance(model)[, "%IncMSE"]
)

ggplot(importance_df, aes(x = reorder(variable, importance), y = importance)) +
  geom_col() +
  coord_flip() +
  labs(
    x = "Variable",
    y = "Importance (%IncMSE)",
    title = "Importance des variables - Random Forest (modèle final, toutes données)"
  ) +
  theme_bw()

# =====================================================================
# Prédiction du NASC sur toute la grille, pour une date donnée
# =====================================================================

day_ds <- readRDS("F:/data_elise/ds_day_ftle_pig_fod/ds_ftle_pig_fod_20230126.rds")

grid_points <- expand.grid(lon = day_ds$lon, lat = day_ds$lat)

grid_points$ftle <- as.vector(day_ds$ftle)

for (p in names(day_ds$pig)) {
  grid_points[[p]] <- as.vector(day_ds$pig[[p]])
}

grid_points$total_chla      <- grid_points$Chla
grid_points$per_ratio_chla  <- grid_points$Per  / grid_points$Chla
grid_points$but_ratio_chla  <- grid_points$But  / grid_points$Chla
grid_points$fuco_ratio_chla <- grid_points$Fuco / grid_points$Chla
grid_points$hex_ratio_chla  <- grid_points$Hex  / grid_points$Chla
grid_points$allo_ratio_chla <- grid_points$Allo / grid_points$Chla
grid_points$zea_ratio_chla  <- grid_points$Zea  / grid_points$Chla
grid_points$chlb_ratio_chla <- grid_points$Chlb / grid_points$Chla

grid_points$fod <- formatC(as.vector(day_ds$fod), width = 2)
grid_points$fod <- factor(grid_points$fod, levels = fod_levels)

grid_points_clean <- grid_points[stats::complete.cases(grid_points[, covariates_all]), ]

cat("Points valides :", nrow(grid_points_clean), "/", nrow(grid_points), "\n")

grid_points_clean$NASC_pred <- predict(model, newdata = grid_points_clean)

ggplot(grid_points_clean, aes(x = lon, y = lat, fill = NASC_pred)) +
  geom_raster() +
  scale_fill_viridis_c() +
  coord_quickmap() +
  theme_bw() +
  labs(
    title = paste("NASC prédit -", format(day_ds$date, "%Y-%m-%d")),
    subtitle = paste0(nrow(grid_points_clean), " / ", nrow(grid_points), " pixels prédits"),
    x = "Longitude", y = "Latitude",
    fill = "log10(NASC)"
  )

