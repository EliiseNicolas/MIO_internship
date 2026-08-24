# LIbraries
library(rpart)
library(rpart.plot)
library(ggplot2)
library(randomForest)
library(mgcv)

rm(list = ls())

# -------------- Global variables
freq <- 120

path_datas <- paste0("F:/data_elise/ds_NASC_pig_ftle_fod/ds_NASC_mean_pig_grid_all/ds_NASC_mean_pig_grid_pig_ftle_fod_2018_2021_2023_transect_", freq, "kHz_mask9.rds")
datas <- readRDS(path_datas)
str(datas)

# -------------- Hyperparameters

diurnal_period <- 3 # 3 : day, 1: night
dp <- "day"
split <- "RS (80/20)" # or "LOYO" #
data_leakage_spatio_temp <- TRUE
dist_thr <- 0.25
pigment_type <-"total_ratio" # # or  or "total_ratio" "conc"  "chla_ratio" # 
n_trees <- 100

####### 1 - filter datas day/night

datas <- datas[datas$day == diurnal_period,]
str(datas)

####### 2 - filtrer NASC - remove extreme values

q <- quantile(datas$nasc, probs = c(0.05, 0.95), na.rm = TRUE)
datas <- datas |>
  dplyr::filter(nasc >= q[1], nasc <= q[2])

####### 3 - plot distribution du NASC

ggplot(
  datas,
  aes(
    x = nasc,
    y = after_stat(count / sum(count) * 100)
  )
) +
  geom_histogram(
    bins = 200
  ) +
  scale_x_log10(
    labels = scales::label_number()
  ) +
  theme_bw() +
  labs(
    x = "NASC (log10 scale)",
    y = "Percentage (%)",
    title = paste("NASC distribution - 2021, 2022, 2023 - ", dp, "-", freq, "kHz"),
    subtitle="NASC without extreme values (0.05 - 0.95) "
  )

####### 4 - Build regression dataset
df <- data.frame(
  NASC = datas$nasc,
  year = format(datas$time_nasc, "%Y"),
  time = datas$time_nasc,
  fod = datas$fod,
  total_pig = datas$total_pig,
  ftle = datas$ftle
)

if (pigment_type == "chla_ratio"){
  print("Chla")
  vars_num <- c("per", "but", "fuco", "hex", "allo", "zea", "chlb", "total_pig", "ftle")
  # chla = datas$Chla,
  df$per = datas$Per_Chla
  df$but = datas$But_Chla
  df$fuco = datas$Fuco_Chla
  df$hex = datas$Hex_Chla
  df$allo = datas$Allo_Chla
  df$zea = datas$Zea_Chla
  df$chlb = datas$Chlb_Chla
}

if (pigment_type == "total_ratio"){
  vars_num <- c("chla", "per", "but", "fuco", "hex", "allo", "zea", "chlb", "total_pig", "ftle")
  print("total")
  df$chla = datas$Chla_total
  df$per = datas$Per_total
  df$but = datas$But_total
  df$fuco = datas$Fuco_total
  df$hex = datas$Hex_total
  df$allo = datas$Allo_total
  df$zea = datas$Zea_total
  df$chlb = datas$Chlb_total
}

if (pigment_type == "conc"){
  print("conc")
  vars_num <- c("chla", "per", "but", "fuco", "hex", "allo", "zea", "chlb", "total_pig", "ftle")
  df$chla = datas$Chla
  df$per = datas$Per
  df$but = datas$But
  df$fuco = datas$Fuco
  df$hex = datas$Hex
  df$allo = datas$Allo
  df$zea = datas$Zea
  df$chlb = datas$Chlb
}

str(df)

####### 5 - Diagnostique du nombre de données
print(c("Nombre total de données avant filtrage : ", nrow(df)))
colSums(is.na(df))
# df <- df[!is.na(df$total_pig) & !is.na(df$fuco) & !is.na(df$zea) & !is.na(df$chla) & !is.na(df$per) & !is.na(df$allo)& !df$fod=="NA" &!is.na(df$ftle),]#  & df$fod=="NA"
df <- df |>
  dplyr::filter(
    if_all(all_of(setdiff(vars_num, "NASC")), ~ !is.na(.)),
    !is.na(fod),
    fod != "NA"
  )
print(c("Nombre total de données après filtrage : ", nrow(df)))
#  day - 18kHz : 309 données, 120kHz 312 données

colSums(is.na(df))

####### 6 - Dataleakage spatio-temporel : retirer les lignes trop proches
# Deux observations sont considérées comme proches si elles sont à moins d'une heure et qu'elles ont une distance > dist_thr
if(data_leakage_spatio_temp){
  # verif lignes dupliquées
  print(c("Nombre de lignes dupliquées : ", sum(duplicated(df)))) 
  
  # Standardisation
  df_scaled <- scale(df[, vars_num])

  # Matrice de distances
  dist_mat <- as.matrix(dist(df_scaled))

  # Matrice des différences temporelles en secondes
  time_diff <- abs(outer(as.numeric(df$time), as.numeric(df$time),FUN = "-"))

  # Garder les paires suffisemment proches et à moins d'1h
  dist_mat[time_diff > 3600] <- NA
  dist_mat[lower.tri(dist_mat, diag = TRUE)] <- NA
  pairs <- which(!is.na(dist_mat), arr.ind = TRUE)
  
  if (nrow(pairs) == 0) {
    message("Aucune paire d'observations à moins d'une heure.")
  } else {
    res <- data.frame(ligne1 = pairs[, 1], ligne2 = pairs[, 2], distance = dist_mat[pairs])
    # Trier par similarité
    res <- res[order(res$distance), ]
  
    # Les 50 paires les plus proches
    res <- head(res, 50)
    
    # infos sur les paires 
    res$time1 <- df$time[res$ligne1]
    res$time2 <- df$time[res$ligne2]
  
    res$time_diff_min <- as.numeric(
      difftime(res$time2, res$time1, units = "mins")
    )
    
    # créer un label pour chaque paire
    res$obs1_label <- paste0(res$ligne1, " | ", format(res$time1, "%Y-%m-%d %H:%M"))
    res$obs2_label <- paste0(res$ligne2, " | ", format(res$time2, "%Y-%m-%d %H:%M"))
    res$pair <- paste0(res$obs1_label," ↔ ",res$obs2_label)
    
    # Mettre les paires dans l'ordre de distance
    res$pair <- factor(res$pair, levels = rev(res$pair))
    
    # Plot des 50 lignes les plus proches
    p <- ggplot(res, aes(x = distance, y = pair)) +
    geom_point(size = 3) +
    theme_minimal() +
    labs(
      x = "Distance",
      y = "Paire d'observations",
      title = "Observations les plus proches",
    )
    print(p)
    
    # retirer les lignes trop proches
    print(c("nombre de lignes trop proches : ", sum(dist_mat < dist_thr, na.rm = TRUE), "thr = ", dist_thr))
    pairs <- which(dist_mat < dist_thr, arr.ind = TRUE)
  
    # Supprimer uniquement la 2e observation de chaque paire
    rows_to_remove <- unique(pairs[, 2])
    cat("Nombre de lignes supprimées :", length(rows_to_remove), "\n")
    cat("Nombre de lignes restantes :", nrow(df), "\n")
    df <- df[-rows_to_remove, ]
    
    # Heatmap de la matrice de distance après filtrage
    # Standardisation
    df_scaled <- scale(df[, vars_num])
    
    # Matrice de distances
    dist_mat <- as.matrix(dist(df_scaled))
    
    # Matrice des différences temporelles en secondes
    time_diff <- abs(outer(as.numeric(df$time), as.numeric(df$time),FUN = "-"))
    
    # Garder les paires suffisemment proches et à moins d'1h
    dist_mat[time_diff > 3600] <- NA
    dist_mat[lower.tri(dist_mat, diag = TRUE)] <- NA
    
    # Transformer en format long
    heatmap_data <- as.data.frame(as.table(dist_mat))
    print(dim(dist_mat))
    print(class(dist_mat))
    print(length(vars_num))
    
    heatmap_data <- as.data.frame(as.table(dist_mat))
    
    print(names(heatmap_data))
    print(dim(heatmap_data))
    
    colnames(heatmap_data) <- c("obs1", "obs2", "distance")
    heatmap_data <- heatmap_data[!is.na(heatmap_data$distance), ]
  
    # Plot de la heatmap
    p <- ggplot(
      heatmap_data,
      aes(
        x = obs1,
        y = obs2,
        fill = distance
      )
    ) +
      geom_tile() +
      scale_fill_viridis_c(
        option = "viridis",
        na.value = "white"
      ) +
      theme_minimal() +
      labs(
        x = "Observation",
        y = "Observation",
        fill = "Distance",
        title = "Matrice de distances",
        subtitle="Données à moins d'une heure d'intervalle et dont la distance est >0.25"
      ) +
      theme(
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank()
      ) +
      coord_fixed()
    
    print(p)
  }
}


####### 6 - Split train/test

set.seed(123)
n <- nrow(df)

if (split == "RS (80/20)"){
  print("RS")
  train_index <- sample(seq_len(n), size = 0.8*n)
  
  train <- df[train_index, ]
  test  <- df[-train_index, ]
}

if (split == "LOYO"){
  print("LOYO")
  # Nombre d'observations par année
  year_counts <- table(df$year)
  print(year_counts)

  # Année avec le moins d'observations
  test_year <- names(which.min(year_counts))
  cat("Année utilisée pour le test :", test_year, "\n")
  
  # split train/test
  train_index <- df$year != test_year
  test_index <- df$year == test_year
  train <- df[train_index, ]
  test <- df[test_index, ]

  cat("Train :", nrow(train), "observations\n", nrow(train) / nrow(df) * 100, "%\n")
  cat("Test :", nrow(test), "observations\n", nrow(test) / nrow(df) * 100, "%\n")
}

# retirer le time et l'année du df
train$year <- NULL
test$year <- NULL
train$time <- NULL
test$time <- NULL

# scale
train_scaled <- scale(train[, vars_num])
train_center <- attr(train_scaled, "scaled:center")
train_scale <- attr(train_scaled, "scaled:scale")

test_scaled <- scale(test[, vars_num], center = train_center, scale = train_scale)

# remise en forme de dataframe
ds_train_scaled <- train
ds_train_scaled[, vars_num] <- train_scaled
ds_test_scaled <- test
ds_test_scaled[, vars_num] <- test_scaled

# dernières vérifications
colSums(is.na(ds_train_scaled))
cat("\n--- Dimensions finales ---\n")
cat("Train :", nrow(ds_train_scaled), "observations\n")
cat("Test :", nrow(ds_test_scaled), "observations\n")

# verifier le dataleakage
if(data_leakage_spatio_temp){
  is_train <- seq_len(nrow(df)) %in% train_index
  is_test  <- !is_train
  
  # Garder uniquement les paires TRAIN <-> TEST
  train_test_pairs <- outer(is_train, is_test, FUN = "&")
  
  # Nombre de paires train-test avec distance < dist_thr
  sum(dist_mat < dist_thr & train_test_pairs, na.rm = TRUE)
  pairs <- which(dist_mat < dist_thr & train_test_pairs, arr.ind = TRUE)
  
  # dataset de leakage
  leakage <- data.frame(
      train = pairs[, 1],
      test = pairs[, 2],
      distance = dist_mat[pairs]
    )
  leakage$time_train <- df$time[leakage$train]
  leakage$time_test  <- df$time[leakage$test]
  leakage$time_diff_min <- abs(as.numeric(difftime(leakage$time_train,leakage$time_test,units = "mins")))
  print(c("Nombre de donnée leaked :", leakage, "nombre de données totales du test set: ", nrow(ds_test_scaled)))
}


####### 7 - Regression tree
model_type <- "RF"

if (model_type == "CART"){
  
  # train model
  model <- rpart(
    NASC ~ ., # entries = all datas but not NASC
    data = ds_train_scaled,
    method = "anova"
  )
  
  # plot tree
  rpart.plot(
    model,
    extra = 101,
    main = paste("Regression tree - CART algorithm -", split, "- data scaled")
  )
  
  # variables importance 
  # Importance des variables
  importance <- model$variable.importance
  importance_pct <- importance / sum(importance) * 100
  
  importance_df <- data.frame(
    variable = names(importance_pct),
    importance = as.numeric(importance_pct)
  )
  importance_df <- importance_df[order(importance_df$importance), ] # sort par ordre d'importance décroissant
}

if (model_type == "RF"){
  # train
  model <- randomForest(
    NASC ~ .,
    data = ds_train_scaled,
    ntree = n_trees,
    importance = TRUE
  )
  
  # variables importance
  importance(model)
  
  # %IncMSE = importance par permutation
  importance_df <- data.frame(
    variable = rownames(importance(model)),
    importance = importance(model)[, "%IncMSE"]
  )
  
  importance_df <- importance_df[
    order(importance_df$importance),
  ]
}

####### 8 - Compute predictions
observed <- ds_test_scaled$NASC

prediction <- predict(model, newdata = ds_test_scaled)

resultats <- data.frame(
  NASC_reel = observed,
  NASC_predit = prediction
)

####### 9 - Compute statistics

RMSE <- sqrt(mean((prediction - observed)^2))

R2 <- 1 - sum((observed-prediction)^2) / sum((observed-mean(observed))^2)

MAE <- mean(abs(prediction - observed))

####### 10 - Plot prediction results

# Plot variable importance

ggplot(importance_df, aes( # Plot
  x = reorder(variable, importance),
  y = importance
)) +
  geom_col() +
  coord_flip() +
  theme_bw() +
  labs(
    x = "Variable",
    y = "Importance (%)",
    title = paste("Variable importance for NASC prediction -", model_type, "algorithm"),
    subtitle=paste("Transect 2018, 2021, 2023, ", dp, ",", freq, "kHz, split : ", split, "on scaled data","\n", "Pigments types : ", pigment_type)
  )

# Plot predictions
par(mar = c(5, 4, 6, 2) + 0.1)
lim <- range(
  c(observed, prediction),
  na.rm = TRUE
)
print(lim)
plot(
  observed,
  prediction,
  xlim = lim,
  ylim =lim,
  xlab = "True NASC",
  ylab = "Predicted NASC",
  main = paste("Real NASC vs. Predicted NASC -", model_type, "algorithm"),
  pch = 4,
  cex = 0.5
)
mtext(
  paste("Transect 2018, 2021, 2023, ", dp, ",", freq, "kHz, split : ", split, "on scaled datas", "\n", "pigments types : ", pigment_type), # "Split: leave-one-year-out (90/10)",
  side = 3,
  line = 0.5,
  cex = 0.9
)

abline(0, 1, col = "red")

legend(
  "bottomright",
  legend = c(paste0("Predictions (RMSE = ", round(RMSE, 2), ", R² =", round(R2, 2), ")"), "Identity line (y = x)"),
  col = c("black", "red"),
  pch = c(4, NA),
  lty = c(NA, 1)
)

####### 6 - TRUE LOYO (Leave-One-Year-Out)

set.seed(123)

years <- sort(unique(df$year))
print(years)

# Vérification
if (length(years) < 2) {
  stop("Il faut au moins 2 années différentes pour faire un LOYO.")
}

# Liste pour stocker les résultats
results_loyo <- list()
metrics_loyo <- list()
models_loyo <- list()
importance_loyo <- list()


####### 7 - LOYO : une année à la fois en test

for (test_year in years) {
  
  cat("\n========================================\n")
  cat("LOYO - Année test :", test_year, "\n")
  cat("========================================\n")
  
  # -----------------------------
  # Train / test
  # -----------------------------
  
  train <- df[df$year != test_year, ]
  test  <- df[df$year == test_year, ]
  
  cat("Train :", nrow(train), "observations\n")
  cat("Test  :", nrow(test), "observations\n")
  
  cat(
    "Années train :",
    paste(sort(unique(train$year)), collapse = ", "),
    "\n"
  )
  
  
  # -----------------------------
  # Retirer year et time
  # -----------------------------
  
  train_model <- train
  test_model  <- test
  
  train_model$year <- NULL
  test_model$year  <- NULL
  
  train_model$time <- NULL
  test_model$time  <- NULL
  
  
  # -----------------------------
  # Scaling
  # IMPORTANT :
  # calculé uniquement sur le TRAIN
  # -----------------------------
  
  train_scaled <- scale(train_model[, vars_num])
  
  train_center <- attr(train_scaled, "scaled:center")
  train_scale  <- attr(train_scaled, "scaled:scale")
  
  test_scaled <- scale(
    test_model[, vars_num],
    center = train_center,
    scale = train_scale
  )
  
  
  # -----------------------------
  # Reconstruction des datasets
  # -----------------------------
  
  ds_train_scaled <- train_model
  ds_test_scaled  <- test_model
  
  ds_train_scaled[, vars_num] <- train_scaled
  ds_test_scaled[, vars_num]  <- test_scaled
  
  
  # -----------------------------
  # Vérifications
  # -----------------------------
  
  if (any(is.na(ds_train_scaled[, vars_num]))) {
    stop(
      paste(
        "NA dans le train pour l'année test",
        test_year
      )
    )
  }
  
  if (any(is.na(ds_test_scaled[, vars_num]))) {
    stop(
      paste(
        "NA dans le test pour l'année test",
        test_year
      )
    )
  }
  
  
  # -----------------------------
  # Modèle
  # -----------------------------
  
  if (model_type == "CART") {
    
    model <- rpart(
      NASC ~ .,
      data = ds_train_scaled,
      method = "anova"
    )
    
  }
  
  
  if (model_type == "RF") {
    
    model <- randomForest(
      NASC ~ .,
      data = ds_train_scaled,
      ntree = n_trees,
      importance = TRUE
    )
    
  }
  
  
  # -----------------------------
  # Prédictions
  # -----------------------------
  
  observed <- ds_test_scaled$NASC
  
  prediction <- predict(
    model,
    newdata = ds_test_scaled
  )
  
  
  # -----------------------------
  # Statistiques
  # -----------------------------
  
  RMSE <- sqrt(
    mean(
      (prediction - observed)^2
    )
  )
  
  MAE <- mean(
    abs(prediction - observed)
  )
  
  R2 <- 1 -
    sum((observed - prediction)^2) /
    sum((observed - mean(observed))^2)
  
  
  # -----------------------------
  # Modèle nul
  # -----------------------------
  
  # prédiction = moyenne du TRAIN
  prediction_null <- rep(
    mean(ds_train_scaled$NASC),
    length(observed)
  )
  
  RMSE_null <- sqrt(
    mean(
      (prediction_null - observed)^2
    )
  )
  
  R2_null <- 1 -
    sum((observed - prediction_null)^2) /
    sum((observed - mean(observed))^2)
  
  
  # -----------------------------
  # Sauvegarde des métriques
  # -----------------------------
  
  metrics_loyo[[test_year]] <- data.frame(
    test_year = test_year,
    n_train = nrow(ds_train_scaled),
    n_test = nrow(ds_test_scaled),
    RMSE = RMSE,
    MAE = MAE,
    R2 = R2,
    RMSE_null = RMSE_null,
    R2_null = R2_null
  )
  
  
  # -----------------------------
  # Sauvegarde des prédictions
  # -----------------------------
  
  results_loyo[[test_year]] <- data.frame(
    test_year = test_year,
    NASC_reel = observed,
    NASC_predit = prediction
  )
  
  
  # -----------------------------
  # Sauvegarde du modèle
  # -----------------------------
  
  models_loyo[[test_year]] <- model
  
  
  # -----------------------------
  # Importance des variables
  # -----------------------------
  
  if (model_type == "RF") {
    
    imp <- importance(model)
    
    importance_loyo[[test_year]] <- data.frame(
      test_year = test_year,
      variable = rownames(imp),
      IncMSE = imp[, "%IncMSE"]
    )
    
  }
  
  
  # -----------------------------
  # Affichage
  # -----------------------------
  
  cat("\nRésultats année test", test_year, "\n")
  cat("RMSE      :", round(RMSE, 3), "\n")
  cat("MAE       :", round(MAE, 3), "\n")
  cat("R²        :", round(R2, 3), "\n")
  cat("RMSE null :", round(RMSE_null, 3), "\n")
}

####### 8 - Résultats LOYO

metrics_loyo_df <- do.call(
  rbind,
  metrics_loyo
)

print(metrics_loyo_df)