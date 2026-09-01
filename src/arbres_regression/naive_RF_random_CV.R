library(rpart)
library(rpart.plot)
library(ggplot2)
library(randomForest)
library(mgcv)

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

cat("Nombre d'observations après filtrage :", nrow(df), "\n") # 47000
df$year <- NULL; df$time<- NULL; df$lat<-NULL; df$lon <- NULL
str(df)
# ---- Split train/test 80/20 ----
n <- nrow(df)
train_index <- sample(seq_len(n), size = 0.8 * n)

train <- df[train_index, c(response_var, covariates_all)]
test  <- df[-train_index, c(response_var, covariates_all)]

# ---- Random Forest ----
model <- randomForest(
  NASC ~ .,
  data = train,
  ntree = 300,
  importance = TRUE
)

# ---- Prédictions ----
pred <- predict(model, newdata = test)
obs  <- test$NASC

RMSE <- sqrt(mean((pred - obs)^2))
R2   <- 1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2)

cat("RMSE :", RMSE, "\n")
cat("R2   :", R2, "\n")

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
    title = "NASC observé vs prédit - Random Forest (80/20)"
  ) +
  theme_bw()

# ---- Importance des variables ----
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
    title = "Importance des variables - Random Forest"
  ) +
  theme_bw()

day_ds <- readRDS("F:/data_elise/ds_day_ftle_pig_fod/ds_ftle_pig_fod_20230126.rds")

# ============================================================
# Prediction du NASC sur toute la grille, pour une date donnée
# ============================================================

library(ggplot2)

# day_ds = objet créé précédemment (déjà tout sur une grille commune
# lon/lat, résolution pigments 1080x720)

# ------------------------------------------------------------
# 1. Assemblage du data.frame de prédiction
# ------------------------------------------------------------
# expand.grid fait varier lon en premier (le plus vite), exactement
# comme as.vector() sur une matrice [lon, lat] -> l'ordre est cohérent

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
grid_points$fod <- factor(grid_points$fod, levels = model$forest$xlevels$fod)

# ------------------------------------------------------------
# 2. Filtrer les points incomplets
# ------------------------------------------------------------
# grid_points_clean <- grid_points[stats::complete.cases(grid_points[, covariates_all]), ]
# 
# cat("Points valides :", nrow(grid_points_clean), "/", nrow(grid_points), "\n")

# ------------------------------------------------------------
# 3. Prédiction
# ------------------------------------------------------------

grid_points$NASC_pred <- predict(model, newdata = grid_points)

# ------------------------------------------------------------
# 4. Carte
# ------------------------------------------------------------

ggplot(grid_points, aes(x = lon, y = lat, fill = NASC_pred)) +
  geom_raster() +
  scale_fill_viridis_c() +
  coord_quickmap() +
  theme_bw() +
  labs(
    title = paste("NASC prédit -", format(day_ds$date, "%Y-%m-%d")),
    x = "Longitude", y = "Latitude",
    fill = "log10(NASC)"
  )