library(ggplot2)
library(xgboost)
library(patchwork)

# ---- Dossiers de sortie ----
out_dir   <- "naive_XGB_CV_10_folds_mean_prediction"   # figures
model_dir <- "naive_XGB_CV_10_folds_models"            # modèles .rds
dir.create(out_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(model_dir, showWarnings = FALSE, recursive = TRUE)

n_CV <- 10

for (freq in c(38, 70, 120, 200)){
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
  # IMPORTANT : contrairement à la version randomForest, on NE FILTRE PLUS
  # les lignes avec NA sur les pigments -- XGBoost apprend une direction
  # par défaut pour les valeurs manquantes à chaque split. On garde
  # uniquement le filtre sur fod (le facteur qu'on veut prédire correctement)
  # et sur NASC lui-même, qui doit toujours être connu pour l'entraînement.
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
  
  feature_names <- c(covariates_num, paste0("fod_", trimws(fod_levels)))
  
  n <- nrow(df)
  
  params <- list(
    max_depth = 4,
    eta = 0.1,
    min_child_weight = 1,
    objective = "reg:squarederror"
  )
  
  # ============================================================
  # ---- XGBoost répété n_CV fois (Monte-Carlo CV) ----
  # A chaque itération : nouveau tirage aléatoire d'un split
  # train/test 80/20, entraînement d'un modèle sur ce train.
  # - Les n_CV modèles servent à moyenner la prédiction sur la
  #   grille finale.
  # - Pour l'évaluation, on empile les couples (obs, pred) des
  #   n_CV tests (chacun sur son propre tirage) et on calcule
  #   aussi un RMSE/R2 moyen sur les n_CV itérations.
  # ============================================================
  
  models_list     <- vector("list", n_CV)
  importance_list <- vector("list", n_CV)
  results_list    <- vector("list", n_CV)
  RMSE_vec        <- numeric(n_CV)
  R2_vec          <- numeric(n_CV)
  
  for (k in 1:n_CV) {
    train_index <- sample(seq_len(n), size = 0.8 * n)
    
    train <- df[train_index, ]
    test  <- df[-train_index, ]
    
    X_train <- build_design_matrix(train, covariates_num, fod_levels)
    y_train <- train$NASC
    
    X_test <- build_design_matrix(test, covariates_num, fod_levels)
    y_test <- test$NASC
    
    dtrain <- xgb.DMatrix(data = X_train, label = y_train, missing = NA)
    dtest  <- xgb.DMatrix(data = X_test,  label = y_test,  missing = NA)
    
    model_k <- xgb.train(
      params  = params,
      data    = dtrain,
      nrounds = 500,
      verbose = 0
    )
    
    pred_k <- predict(model_k, dtest)
    obs_k  <- y_test
    
    models_list[[k]]     <- model_k
    importance_list[[k]] <- xgb.importance(model = model_k)
    results_list[[k]]    <- data.frame(obs = obs_k, pred = pred_k, fold = k)
    
    RMSE_vec[k] <- sqrt(mean((pred_k - obs_k)^2))
    R2_vec[k]   <- 1 - sum((obs_k - pred_k)^2) / sum((obs_k - mean(obs_k))^2)
  }
  
  RMSE <- mean(RMSE_vec)
  R2   <- mean(R2_vec)
  
  cat("RMSE moyen sur", n_CV, "tirages :", RMSE, "(sd =", sd(RMSE_vec), ")\n")
  cat("R2   moyen sur", n_CV, "tirages :", R2,   "(sd =", sd(R2_vec),   ")\n")
  
  # ---- Sauvegarde du modèle (ensemble des n_CV boosters) ----
  # On sauvegarde aussi covariates_num et fod_levels, nécessaires
  # pour reconstruire la matrice de design lors de prédictions futures.
  xgb_saved <- list(
    models_list    = models_list,
    covariates_num = covariates_num,
    covariates_all = covariates_all,
    fod_levels     = fod_levels,
    feature_names  = feature_names,
    freq           = freq,
    RMSE_mean      = RMSE,
    R2_mean        = R2
  )
  
  saveRDS(
    xgb_saved,
    file = file.path(model_dir, paste0("XGB_model_", freq, "kHz.rds"))
  )
  
  # ---- Plot NASC observé vs prédit (points des n_CV tirages empilés) ----
  results <- do.call(rbind, results_list)
  
  p <- ggplot(results, aes(x = obs, y = pred)) +
    geom_point(alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    annotate(
      "label",
      x = min(results$obs), y = max(results$pred),
      hjust = 0, vjust = 1,
      label = sprintf("RMSE = %.3f\nR² = %.3f", RMSE, R2)
    ) +
    labs(
      x = "NASC observé (log10)",
      y = "NASC prédit (log10)",
      title = paste0("NASC observé vs prédit - ", n_CV, " tirages XGBoost - ", freq, " kHz")
    ) +
    theme_bw()
  print(p)
  ggsave(
    filename = file.path(out_dir, paste0("obs_vs_pred_", freq, "kHz.png")),
    plot = p, width = 7, height = 6, dpi = 300
  )
  
  # ---- Importance des variables moyennée sur les n_CV modèles ----
  # (une feature absente d'un modèle -> Gain = 0 pour ce modèle)
  importance_mat <- sapply(importance_list, function(imp) {
    full <- setNames(rep(0, length(feature_names)), feature_names)
    full[imp$Feature] <- imp$Gain
    full
  })
  
  importance_df <- data.frame(
    variable   = feature_names,
    importance = rowMeans(importance_mat)
  )
  
  p <- ggplot(importance_df, aes(x = reorder(variable, importance), y = importance)) +
    geom_col() +
    coord_flip() +
    labs(
      x = "Variable",
      y = "Importance moyenne (Gain)",
      title = paste0("Importance des variables - moyenne sur ", n_CV, " modèles XGBoost - ", freq, " kHz")
    ) +
    theme_bw()
  print(p)
  ggsave(
    filename = file.path(out_dir, paste0("importance_variables_", freq, "kHz.png")),
    plot = p, width = 7, height = 6, dpi = 300
  )
  
  # ============================================================
  # Prediction du NASC sur toute la grille, pour une date donnée
  # (avec NA autorisés sur les pigments -- XGBoost les gère nativement)
  # ============================================================
  
  day_ds <- readRDS("F:/data_elise/ds_day_ftle_pig_fod/ds_ftle_pig_fod_20230126.rds")
  
  grid_points <- expand.grid(lon = day_ds$lon, lat = day_ds$lat)
  
  grid_points$ftle <- as.vector(day_ds$ftle)
  
  # NB: renommé p -> p_name pour ne pas écraser la variable ggplot "p" utilisée plus haut
  for (p_name in names(day_ds$pig)) {
    grid_points[[p_name]] <- as.vector(day_ds$pig[[p_name]])
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
  # Construction de la matrice de design (NA conservés, plus besoin
  # de filtrer les pixels incomplets)
  # ------------------------------------------------------------
  
  X_grid <- build_design_matrix(grid_points, covariates_num, fod_levels)
  dgrid  <- xgb.DMatrix(data = X_grid, missing = NA)
  
  # ---- Prédiction = moyenne des n_CV modèles ----
  grid_pred_matrix   <- sapply(models_list, function(m) predict(m, dgrid))
  grid_points$NASC_pred <- rowMeans(grid_pred_matrix)
  
  cat("Pixels prédits (tous, y compris avec NA partiels) :", nrow(grid_points), "\n")
  cat("Dont pixels avec au moins un pigment manquant :",
      sum(!stats::complete.cases(grid_points[, covariates_num])), "\n")
  
  # ------------------------------------------------------------
  # Prédiction SANS les pixels incomplets (comparaison, comme avec RF)
  # ------------------------------------------------------------
  
  grid_points_clean <- grid_points[stats::complete.cases(grid_points[, covariates_all]), ]
  
  X_grid_clean <- build_design_matrix(grid_points_clean, covariates_num, fod_levels)
  dgrid_clean  <- xgb.DMatrix(data = X_grid_clean, missing = NA)
  
  grid_pred_clean_matrix    <- sapply(models_list, function(m) predict(m, dgrid_clean))
  grid_points_clean$NASC_pred <- rowMeans(grid_pred_clean_matrix)
  
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
  
  p <- (p_all + p_clean) +
    plot_annotation(
      title = paste0("NASC prédit - moyenne de ", n_CV, " modèles XGBoost - ",
                     format(day_ds$date, "%Y-%m-%d"), " - ", freq, " kHz"),
      subtitle = "Comparaison : gestion native des NA (XGBoost) vs filtrage strict"
    )
  print(p)
  ggsave(
    filename = file.path(out_dir, paste0("carte_NASC_pred_", freq, "kHz.png")),
    plot = p, width = 12, height = 6, dpi = 300
  )
}