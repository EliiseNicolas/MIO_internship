library(ggplot2)
library(xgboost)
library(dplyr)
library(purrr)
library(tibble)

freq <- 120
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

covariates_num <- setdiff(vars_num, "NASC")
covariates_all <- c(covariates_num, "fod")
response_var   <- "NASC"

# ---------------------------------------------------------------------
# IMPORTANT : on NE FILTRE PAS les lignes avec NA sur les pigments --
# XGBoost apprend une direction par défaut pour les valeurs manquantes
# à chaque split. On garde uniquement le filtre sur fod (le facteur qu'on
# veut prédire correctement) et sur NASC lui-même.
# ---------------------------------------------------------------------
df <- df |>
  dplyr::filter(!is.na(fod), fod != "NA", is.finite(NASC))

cat("Nombre d'observations après filtrage :", nrow(df), "\n")

fod_levels <- levels(df$fod)   # à réutiliser tel quel pour la prédiction

# ---------------------------------------------------------------------
# Fonction d'encodage : covariables numériques (NA préservés tels quels)
# + dummies manuels pour fod (NA se propage automatiquement dans chaque
# colonne dummy si fod est NA sur cette ligne -- pas de model.matrix())
# ---------------------------------------------------------------------
build_design_matrix <- function(data, covs_num, fod_levels) {
  mat_num <- as.matrix(data[, covs_num])
  
  fod_dummies <- sapply(fod_levels, function(lvl) as.numeric(data$fod == lvl))
  colnames(fod_dummies) <- paste0("fod_", trimws(fod_levels))
  
  cbind(mat_num, fod_dummies)
}

# ---- Split train/test 80/20 ----
set.seed(123)
n <- nrow(df)
train_index <- sample(seq_len(n), size = 0.8 * n)

train <- df[train_index, ]
test  <- df[-train_index, ]

X_train <- build_design_matrix(train, covariates_num, fod_levels)
y_train <- train$NASC

X_test <- build_design_matrix(test, covariates_num, fod_levels)
y_test <- test$NASC

# =====================================================================
# TUNING SIMPLE DE XGBOOST (max_depth x eta x min_child_weight), via CV
# =====================================================================

dtrain <- xgb.DMatrix(data = X_train, label = y_train, missing = NA)

# Grille de tuning : petite mais couvre les valeurs usuelles
tuning_grid <- expand.grid(
  max_depth        = c(3, 4, 6),
  eta              = c(0.01, 0.05, 0.1),
  min_child_weight = c(1, 3, 5)
)

print(tuning_grid)

# ---------------------------------------------------------------------
# Fonction d'évaluation d'un jeu d'hyperparamètres via xgb.cv
# (CV interne à XGBoost, gère les NA nativement comme en prédiction)
# ---------------------------------------------------------------------
evaluate_params_xgb <- function(max_depth_val, eta_val, min_child_weight_val,
                                dtrain, nfold = 10, nrounds = 500,
                                early_stopping_rounds = 20) {
  
  params <- list(
    max_depth = max_depth_val,
    eta = eta_val,
    min_child_weight = min_child_weight_val,
    subsample = 0.8,
    colsample_bytree = 0.8,
    objective = "reg:squarederror",
    eval_metric = "rmse"
  )
  
  cv_result <- xgb.cv(
    params = params,
    data = dtrain,
    nrounds = nrounds,
    nfold = nfold,
    early_stopping_rounds = early_stopping_rounds,
    verbose = 0
  )
  
  best_iter <- cv_result$best_iteration
  best_row  <- cv_result$evaluation_log[best_iter, ]
  
  tibble(
    max_depth        = max_depth_val,
    eta              = eta_val,
    min_child_weight = min_child_weight_val,
    best_nrounds     = best_iter,
    mean_rmse_test   = best_row$test_rmse_mean,
    sd_rmse_test     = best_row$test_rmse_std,
    mean_rmse_train  = best_row$train_rmse_mean
  )
}

# ---------------------------------------------------------------------
# Boucle sur toute la grille
# ---------------------------------------------------------------------
cat("\n=== Tuning XGBoost en cours (", nrow(tuning_grid), "combinaisons, 10-fold CV interne) ===\n")

tuning_results <- pmap_dfr(tuning_grid, function(max_depth, eta, min_child_weight) {
  cat("  max_depth =", max_depth, "| eta =", eta, "| min_child_weight =", min_child_weight, "\n")
  evaluate_params_xgb(max_depth, eta, min_child_weight, dtrain)
})

tuning_results <- tuning_results %>% arrange(mean_rmse_test)

cat("\n=== Résultats du tuning (triés par RMSE croissant) ===\n")
print(tuning_results)

best_params <- tuning_results[1, ]
cat("\nMeilleurs hyperparamètres :\n")
cat("  max_depth        =", best_params$max_depth, "\n")
cat("  eta              =", best_params$eta, "\n")
cat("  min_child_weight =", best_params$min_child_weight, "\n")
cat("  nrounds optimal  =", best_params$best_nrounds, "\n")
cat("  RMSE moyen (CV)  =", round(best_params$mean_rmse_test, 4), "\n")

# ---------------------------------------------------------------------
# Visualisation : RMSE selon max_depth et eta, facetté par min_child_weight
# ---------------------------------------------------------------------
ggplot(tuning_results, aes(x = factor(eta), y = mean_rmse_test, fill = factor(max_depth))) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(
    aes(ymin = mean_rmse_test - sd_rmse_test, ymax = mean_rmse_test + sd_rmse_test),
    position = position_dodge(width = 0.8), width = 0.2
  ) +
  facet_wrap(~min_child_weight, labeller = label_both) +
  labs(
    title = "Tuning XGBoost (CV interne, nrounds optimal par early stopping)",
    x = "eta", y = "RMSE moyen (± sd CV)", fill = "max_depth"
  ) +
  theme_bw()

# ---------------------------------------------------------------------
# Modèle final avec les meilleurs hyperparamètres, entraîné sur TOUT
# le train (X_train/y_train), avec le nombre de rounds optimal trouvé
# ---------------------------------------------------------------------
final_params <- list(
  max_depth = best_params$max_depth,
  eta = best_params$eta,
  min_child_weight = best_params$min_child_weight,
  subsample = 0.8,
  colsample_bytree = 0.8,
  objective = "reg:squarederror"
)

model <- xgb.train(
  params = final_params,
  data = dtrain,
  nrounds = best_params$best_nrounds,
  verbose = 0
)

cat("\nModèle final entraîné avec les meilleurs hyperparamètres et",
    best_params$best_nrounds, "rounds\n")

# ---- Réévaluation sur le test set avec le modèle tuné ----
dtest <- xgb.DMatrix(data = X_test, label = y_test, missing = NA)
pred <- predict(model, dtest)
obs  <- y_test

RMSE <- sqrt(mean((pred - obs)^2))
R2   <- 1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2)

cat("RMSE (test, après tuning) :", RMSE, "\n")
cat("R2   (test, après tuning) :", R2, "\n")

# ---------------------------------------------------------------------
# CARTE DES RÉSIDUS SPATIALISÉE (test set)
#    -- permet de voir si le modèle est plus performant dans certaines
#       zones géographiques que d'autres
# ---------------------------------------------------------------------

results_spatial <- data.frame(
  lon      = test$lon,
  lat      = test$lat,
  obs      = obs,
  pred     = pred,
  residual = obs - pred
)

ggplot(results_spatial, aes(x = lon, y = lat, color = residual)) +
  geom_point(size = 1.2, alpha = 0.7) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  coord_quickmap() +
  labs(
    title = "Résidus spatialisés (observé - prédit) - XGBoost tuné",
    subtitle = "Bleu = sous-estimation, Rouge = surestimation",
    x = "Longitude", y = "Latitude", color = "Résidu\n(log10 NASC)"
  ) +
  theme_bw()

# Version alternative : amplitude de l'erreur (sans distinguer sur/sous-
# estimation) -- utile pour repérer directement les "zones difficiles"
ggplot(results_spatial, aes(x = lon, y = lat, color = abs(residual))) +
  geom_point(size = 1.2, alpha = 0.7) +
  scale_color_viridis_c(option = "magma", direction = -1) +
  coord_quickmap() +
  labs(
    title = "Amplitude des erreurs de prédiction, spatialisée - XGBoost tuné",
    x = "Longitude", y = "Latitude", color = "|Résidu|\n(log10 NASC)"
  ) +
  theme_bw()

# ---- Plot NASC observé vs prédit ----
results <- data.frame(obs = obs, pred = pred)

ggplot(results, aes(x = obs, y = pred)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  annotate(
    "label",
    x = min(obs), y = max(pred),
    hjust = 0, vjust = 1,
    label = sprintf("RMSE = %.3f\nR² = %.3f", RMSE, R2)
  ) +
  labs(
    x = "NASC observé (log10)",
    y = "NASC prédit (log10)",
    title = "NASC observé vs prédit - XGBoost tuné (80/20)"
  ) +
  theme_bw()

# ---- Importance des variables ----
importance_matrix <- xgb.importance(model = model)

ggplot(importance_matrix, aes(x = reorder(Feature, Gain), y = Gain)) +
  geom_col() +
  coord_flip() +
  labs(
    x = "Variable",
    y = "Importance (Gain)",
    title = "Importance des variables - XGBoost (tuné)"
  ) +
  theme_bw()

# =====================================================================
# Prédiction du NASC sur toute la grille, pour une date donnée
# (avec NA autorisés sur les pigments -- XGBoost les gère nativement)
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

# ratios : Inf/NaN (division par 0 ou 0/0) remplacés par NA pour que
# xgboost les traite comme manquants plutôt que comme des valeurs
# numériques aberrantes
for (v in covariates_num) {
  bad <- !is.finite(grid_points[[v]])
  grid_points[[v]][bad] <- NA
}

grid_points$fod <- formatC(as.vector(day_ds$fod), width = 2)
grid_points$fod <- factor(grid_points$fod, levels = fod_levels)

# ------------------------------------------------------------
# Prédiction AVEC tous les pixels (NA autorisés, gérés par XGBoost)
# ------------------------------------------------------------

X_grid <- build_design_matrix(grid_points, covariates_num, fod_levels)
dgrid  <- xgb.DMatrix(data = X_grid, missing = NA)
grid_points$NASC_pred <- predict(model, dgrid)

cat("Pixels prédits (tous, y compris avec NA partiels) :", nrow(grid_points), "\n")
cat("Dont pixels avec au moins un pigment manquant :",
    sum(!stats::complete.cases(grid_points[, covariates_num])), "\n")

# ------------------------------------------------------------
# Prédiction SANS les pixels incomplets (comparaison, comme avec RF)
# ------------------------------------------------------------

grid_points_clean <- grid_points[stats::complete.cases(grid_points[, covariates_all]), ]

X_grid_clean <- build_design_matrix(grid_points_clean, covariates_num, fod_levels)
dgrid_clean  <- xgb.DMatrix(data = X_grid_clean, missing = NA)
grid_points_clean$NASC_pred <- predict(model, dgrid_clean)

cat("Pixels valides (complete.cases) :", nrow(grid_points_clean), "/", nrow(grid_points), "\n")

# ------------------------------------------------------------
# Cartes comparées
# ------------------------------------------------------------

# Même échelle de couleur pour les deux cartes, pour une comparaison honnête
shared_limits <- range(c(grid_points$NASC_pred, grid_points_clean$NASC_pred), na.rm = TRUE)

p_all <- ggplot(grid_points, aes(x = lon, y = lat, fill = NASC_pred)) +
  geom_raster() +
  scale_fill_viridis_c(limits = shared_limits) +
  coord_quickmap() +
  theme_bw() +
  labs(
    title = "Avec NA (XGBoost gère le manquant)",
    subtitle = paste0(nrow(grid_points), " pixels prédits"),
    x = "Longitude", y = "Latitude",
    fill = "log10(NASC)"
  )

p_clean <- ggplot(grid_points_clean, aes(x = lon, y = lat, fill = NASC_pred)) +
  geom_raster() +
  scale_fill_viridis_c(limits = shared_limits) +
  coord_quickmap() +
  theme_bw() +
  labs(
    title = "Sans NA (pixels complets uniquement)",
    subtitle = paste0(nrow(grid_points_clean), " pixels prédits"),
    x = "Longitude", y = "Latitude",
    fill = "log10(NASC)"
  )

library(patchwork)

(p_all + p_clean) +
  plot_annotation(
    title = paste("NASC prédit -", format(day_ds$date, "%Y-%m-%d"), "-- modèle XGBoost tuné"),
    subtitle = "Comparaison : gestion native des NA (XGBoost) vs filtrage strict"
  )