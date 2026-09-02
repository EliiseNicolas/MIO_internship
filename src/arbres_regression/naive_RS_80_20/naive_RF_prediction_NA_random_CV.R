library(rpart)
library(rpart.plot)
library(ggplot2)
library(randomForestSRC)
library(mgcv)

# ---- Dossier de sortie pour les plots ----
out_dir <- "naive_RF_CV_10_folds_mean_prediction"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Dossier de sortie pour les modèles entraînés ----
model_dir <- "naive_RF_CV_10_folds_models"
dir.create(model_dir, showWarnings = FALSE, recursive = TRUE)

for (freq in c(18, 38, 70, 120, 200)){
  
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
  
  cat("Nombre d'observations après filtrage :", nrow(df), "\n") # 47000
  df$year <- NULL; df$time<- NULL; df$lat<-NULL; df$lon <- NULL
  str(df)
  
  # Niveaux de fod figés une bonne fois pour toutes : servent à
  # reformater grid_points$fod de la même façon quel que soit le
  # tirage train/test (train/test sont des sous-ensembles de df,
  # donc partagent les mêmes niveaux de facteur).
  fod_levels <- levels(df$fod)
  
  n <- nrow(df)
  
  # ============================================================
  # ---- Random Forest (randomForestSRC) répété n_CV fois ----
  # (Monte-Carlo CV) : nouveau tirage aléatoire d'un split
  # train/test 80/20 à chaque itération.
  # na.action = "na.impute" + nimpute = 2 permet au RF de gérer
  # les NA à l'entraînement, ET on repasse na.action = "na.impute"
  # à predict() pour pouvoir prédire directement sur la grille
  # day_ds, qui contient des valeurs manquantes, sans avoir à la
  # filtrer au préalable.
  # ============================================================
  
  models_list     <- vector("list", n_CV)
  importance_list <- vector("list", n_CV)
  results_list    <- vector("list", n_CV)
  RMSE_vec        <- numeric(n_CV)
  R2_vec          <- numeric(n_CV)
  
  for (k in 1:n_CV) {
    train_index <- sample(seq_len(n), size = 0.8 * n)
    
    train <- df[train_index, c(response_var, covariates_all)]
    test  <- df[-train_index, c(response_var, covariates_all)]
    
    model_k <- rfsrc(
      NASC ~ .,
      data       = train,
      ntree      = 1000,
      mtry       = 4,
      nodesize   = 10,
      na.action  = "na.impute",
      nimpute    = 2,
      importance = TRUE
    )
    
    pred_k <- predict(model_k, newdata = test, na.action = "na.impute")$predicted
    obs_k  <- test$NASC
    
    models_list[[k]]     <- model_k
    importance_list[[k]] <- model_k$importance
    results_list[[k]]    <- data.frame(obs = obs_k, pred = pred_k, fold = k)
    
    RMSE_vec[k] <- sqrt(mean((pred_k - obs_k)^2))
    R2_vec[k]   <- 1 - sum((obs_k - pred_k)^2) / sum((obs_k - mean(obs_k))^2)
  }
  
  RMSE <- mean(RMSE_vec)
  R2   <- mean(R2_vec)
  
  cat("RMSE moyen sur", n_CV, "tirages :", RMSE, "(sd =", sd(RMSE_vec), ")\n")
  cat("R2   moyen sur", n_CV, "tirages :", R2,   "(sd =", sd(R2_vec),   ")\n")
  
  # ---- Sauvegarde du modèle (ensemble des n_CV forêts) ----
  # On sauvegarde aussi fod_levels et covariates_all, nécessaires
  # pour reconstruire correctement les données à prédire plus tard.
  rf_saved <- list(
    models_list    = models_list,
    covariates_all = covariates_all,
    fod_levels     = fod_levels,
    freq           = freq,
    RMSE_mean      = RMSE,
    R2_mean        = R2
  )
  
  saveRDS(
    rf_saved,
    file = file.path(model_dir, paste0("RF_model_", freq, "kHz.rds"))
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
      title = paste0("NASC observé vs prédit - ", n_CV, " tirages RF (randomForestSRC) - ", freq, " kHz")
    ) +
    theme_bw()
  print(p)
  ggsave(
    filename = file.path(out_dir, paste0("obs_vs_pred_", freq, "kHz.png")),
    plot = p, width = 7, height = 6, dpi = 300
  )
  
  # ---- Importance des variables moyennée sur les n_CV forêts ----
  importance_mat <- do.call(cbind, importance_list)
  importance_df <- data.frame(
    variable   = rownames(importance_mat),
    importance = rowMeans(importance_mat)
  )
  
  p <- ggplot(importance_df, aes(x = reorder(variable, importance), y = importance)) +
    geom_col() +
    coord_flip() +
    labs(
      x = "Variable",
      y = "Importance moyenne (VIMP)",
      title = paste0("Importance des variables - moyenne sur ", n_CV, " forêts - ", freq, " kHz")
    ) +
    theme_bw()
  print(p)
  ggsave(
    filename = file.path(out_dir, paste0("importance_variables_", freq, "kHz.png")),
    plot = p, width = 7, height = 6, dpi = 300
  )
  
  day_ds <- readRDS("F:/data_elise/ds_day_ftle_pig_fod/ds_ftle_pig_fod_20230126.rds")
  
  # ============================================================
  # Prédiction du NASC sur toute la grille, pour une date donnée
  # (la grille peut contenir des NA : gérés par imputation)
  # ============================================================
  
  # day_ds = objet créé précédemment (déjà tout sur une grille commune
  # lon/lat, résolution pigments 1080x720)
  
  # ------------------------------------------------------------
  # 1. Assemblage du data.frame de prédiction
  # ------------------------------------------------------------
  # expand.grid fait varier lon en premier (le plus vite), exactement
  # comme as.vector() sur une matrice [lon, lat] -> l'ordre est cohérent
  
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
  
  grid_points$fod <- formatC(as.vector(day_ds$fod), width = 2)
  grid_points$fod <- factor(grid_points$fod, levels = fod_levels)
  
  # ------------------------------------------------------------
  # 2. Prédiction = moyenne des n_CV forêts
  # ------------------------------------------------------------
  # na.action = "na.impute" côté predict() : plus besoin de filtrer
  # les points incomplets de la grille au préalable, le RF impute
  # les covariables manquantes avant de prédire.
  
  grid_pred_matrix <- sapply(
    models_list,
    function(m) predict(m, newdata = grid_points, na.action = "na.impute")$predicted
  )
  grid_points$NASC_pred <- rowMeans(grid_pred_matrix)
  
  # ------------------------------------------------------------
  # 3. Carte
  # ------------------------------------------------------------
  
  p <- ggplot(grid_points, aes(x = lon, y = lat, fill = NASC_pred)) +
    geom_raster() +
    scale_fill_viridis_c() +
    coord_quickmap() +
    theme_bw() +
    labs(
      title = paste0("NASC prédit - moyenne de ", n_CV, " forêts (randomForestSRC) - ", freq, " kHz - ", format(day_ds$date, "%Y-%m-%d")),
      x = "Longitude", y = "Latitude",
      fill = "log10(NASC)"
    )
  print(p)
  ggsave(
    filename = file.path(out_dir, paste0("carte_NASC_pred_", freq, "kHz.png")),
    plot = p, width = 8, height = 6, dpi = 300
  )
}