library(rpart)
library(rpart.plot)
library(ggplot2)
library(randomForest)
library(mgcv)

# ---- Dossier de sortie pour les plots ----
out_dir <- "C:/Users/mmolinet/elisou_ta_stagiaire_pref/MIO_internship_III/figures/trees/naive_RF_XGB_RS_80_20/naive_RF_CV_10_folds_mean_prediction"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

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
  
  n <- nrow(df)
  
  # ============================================================
  # ---- Random Forest répété n_CV fois (Monte-Carlo CV) ----
  # A chaque itération : nouveau tirage aléatoire d'un split
  # train/test 80/20, entraînement d'une forêt sur ce train.
  # - Les n_CV forêts servent à moyenner la prédiction sur la
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
    
    train <- df[train_index, c(response_var, covariates_all)]
    test  <- df[-train_index, c(response_var, covariates_all)]
    
    model_k <- randomForest(
      NASC ~ .,
      data = train,
      ntree = 300,
      mtry = 4,
      nodesize = 10,
      importance = TRUE
    )
    
    pred_k <- predict(model_k, newdata = test)
    obs_k  <- test$NASC
    
    models_list[[k]]     <- model_k
    importance_list[[k]] <- importance(model_k)[, "%IncMSE"]
    results_list[[k]]    <- data.frame(obs = obs_k, pred = pred_k, fold = k)
    
    RMSE_vec[k] <- sqrt(mean((pred_k - obs_k)^2))
    R2_vec[k]   <- 1 - sum((obs_k - pred_k)^2) / sum((obs_k - mean(obs_k))^2)
  }
  
  RMSE <- mean(RMSE_vec)
  R2   <- mean(R2_vec)
  
  cat("RMSE moyen sur", n_CV, "tirages :", RMSE, "(sd =", sd(RMSE_vec), ")\n")
  cat("R2   moyen sur", n_CV, "tirages :", R2,   "(sd =", sd(R2_vec),   ")\n")
  
  # ---- Plot NASC observé vs prédit (points des n_CV tirages empilés) ----
  results <- do.call(rbind, results_list)
  
  p <- ggplot(results, aes(x = obs, y = pred)) +
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
      title = paste0("NASC observé vs prédit - ", n_CV, " tirages RF empilés - ", freq, " kHz")
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
      y = "Importance moyenne (%IncMSE)",
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
  grid_points$fod <- factor(grid_points$fod, levels = models_list[[1]]$forest$xlevels$fod)
  
  # ------------------------------------------------------------
  # 2. Filtrer les points incomplets
  # ------------------------------------------------------------
  # grid_points_clean <- grid_points[stats::complete.cases(grid_points[, covariates_all]), ]
  # 
  # cat("Points valides :", nrow(grid_points_clean), "/", nrow(grid_points), "\n")
  
  # ------------------------------------------------------------
  # 3. Prédiction = moyenne des n_CV forêts
  # ------------------------------------------------------------
  
  grid_pred_matrix <- sapply(models_list, function(m) predict(m, newdata = grid_points))
  grid_points$NASC_pred <- rowMeans(grid_pred_matrix)
  
  # ------------------------------------------------------------
  # 4. Carte
  # ------------------------------------------------------------
  
  p <- ggplot(grid_points, aes(x = lon, y = lat, fill = NASC_pred)) +
    geom_raster() +
    scale_fill_viridis_c() +
    coord_quickmap() +
    theme_bw() +
    labs(
      title = paste0("NASC prédit - moyenne de ", n_CV, " forêts - ", freq, " kHz - ", format(day_ds$date, "%Y-%m-%d")),
      x = "Longitude", y = "Latitude",
      fill = "log10(NASC)"
    )
  print(p)
  ggsave(
    filename = file.path(out_dir, paste0("carte_NASC_pred_", freq, "kHz.png")),
    plot = p, width = 8, height = 6, dpi = 300
  )
}