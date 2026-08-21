# LIbraries
library(rpart)
library(rpart.plot)
library(ggplot2)
library(randomForest)
library(mgcv)
rm(list = ls())

# -------------- Paths
path_datas <- "F:/data_elise/ds_NASC_pig_ftle_fod/ds_NASC_pig_ftle_fod_2018_2021_2023_transect_120kHz_day_mask9.rds"
datas <- readRDS(path_datas)

# ------------- Keep only relevant variables
df <- data.frame(
  NASC = datas$nasc,
  year = format(datas$time_nasc, "%Y"),
  time = datas$time_nasc,
  fod = as.factor(datas$fod),
  chla = datas$Chla,
  per = datas$Per,
  but = datas$But,
  fuco = datas$Fuco,
  hex = datas$Hex,
  allo = datas$Allo,
  zea = datas$Zea,
  chlb = datas$Chlb,
  ftle = datas$ftle
)
# on retire DvChla des predicteurs car sa moyenne est 0 et sa variance est 0 aussi

print(sum(!is.na(df$chla)))
print(sum(!is.na(df$fod)))
print(c("Nombre total de données", nrow(df)))

# garder uniquement les données sans na
df <- df[!is.na(df$chla) & !is.na(df$fod) & !is.na(df$ftle) & !is.na(df$per) & !is.na(df$allo),]
print(c("Nombre total de données", nrow(df))) # 215
colSums(is.na(df))


# ----------------------------------------------- vérifier le dataleakage 

# -----------------------------------------------
# Vérification du data leakage
# Deux observations sont considérées comme proches
# uniquement si elles sont à moins d'une heure
# -----------------------------------------------

# Vérifier s'il existe des lignes dupliquées
any(duplicated(df)) # Aucune ligne identique

# Variables utilisées pour calculer la similarité
vars_num <- c(
  "NASC", "chla", "per", "but", "fuco", "hex",
  "allo", "zea", "chlb", "ftle"
)

# Standardisation
df_scaled <- scale(df[, vars_num])

# Matrice de distances
dist_mat <- as.matrix(dist(df_scaled))

# Matrice des différences temporelles en secondes
time_diff <- abs(
  outer(
    as.numeric(df$time),
    as.numeric(df$time),
    FUN = "-"
  )
)

# -----------------------------------------------
# Ne garder que les paires :
# 1. différentes
# 2. à moins d'une heure
# -----------------------------------------------

dist_mat[time_diff > 3600] <- NA
dist_mat[lower.tri(dist_mat, diag = TRUE)] <- NA

# -----------------------------------------------
# Récupérer les paires restantes
# -----------------------------------------------

pairs <- which(!is.na(dist_mat), arr.ind = TRUE)

res <- data.frame(
  ligne1 = pairs[, 1],
  ligne2 = pairs[, 2],
  distance = dist_mat[pairs]
)

# Trier par similarité
res <- res[order(res$distance), ]

# Les 50 paires les plus proches
res <- head(res, 50)

# -----------------------------------------------
# Informations supplémentaires
# -----------------------------------------------

res$time1 <- df$time[res$ligne1]
res$time2 <- df$time[res$ligne2]

res$time_diff_min <- as.numeric(
  difftime(res$time2, res$time1, units = "mins")
)

res$time_diff_min <- abs(res$time_diff_min)

res$fod1 <- df$fod[res$ligne1]
res$fod2 <- df$fod[res$ligne2]

res$year1 <- df$year[res$ligne1]
res$year2 <- df$year[res$ligne2]

res

# -----------------------------------------------
# Créer un label pour chaque paire
# -----------------------------------------------

res$obs1_label <- paste0(
  res$ligne1,
  " | ",
  format(res$time1, "%Y-%m-%d %H:%M"),
  " | FOD=", res$fod1
)

res$obs2_label <- paste0(
  res$ligne2,
  " | ",
  format(res$time2, "%Y-%m-%d %H:%M"),
  " | FOD=", res$fod2
)

res$pair <- paste0(
  res$obs1_label,
  " ↔ ",
  res$obs2_label
)

# Mettre les paires dans l'ordre de distance
res$pair <- factor(
  res$pair,
  levels = rev(res$pair)
)

ggplot(res, aes(x = distance, y = pair)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(
    x = "Distance",
    y = "Paire d'observations",
    title = "Observations les plus proches",
  )

sum(dist_mat < 0.5, na.rm = TRUE)

# supprimer les doublons 
pairs <- which(
  dist_mat < 0.25,
  arr.ind = TRUE
)

# Supprimer uniquement la 2e observation de chaque paire
rows_to_remove <- unique(pairs[, 2])

# Nettoyage
df <- df[-rows_to_remove, ]

cat("Nombre de lignes supprimées :", length(rows_to_remove), "\n")
cat("Nombre de lignes restantes :", nrow(df), "\n")

df_scaled <- scale(df[, vars_num])

dist_mat <- as.matrix(dist(df_scaled))

time_diff <- abs(
  outer(
    as.numeric(df$time),
    as.numeric(df$time),
    FUN = "-"
  )
)

dist_mat[time_diff > 3600] <- NA
dist_mat[lower.tri(dist_mat, diag = TRUE)] <- NA

# ------------------------------------------------
# HEATMAP DE LA MATRICE DE DISTANCES APRÈS FILTRAGE
# ------------------------------------------------

library(ggplot2)

# Recalcul de la matrice de distances
df_scaled <- scale(df[, vars_num])

dist_mat <- as.matrix(dist(df_scaled))

# Matrice des différences temporelles
time_diff <- abs(
  outer(
    as.numeric(df$time),
    as.numeric(df$time),
    FUN = "-"
  )
)

# Garder uniquement les paires à moins d'une heure
dist_mat[time_diff > 3600] <- NA
# 
# # Supprimer diagonale + partie inférieure
# dist_mat[lower.tri(dist_mat, diag = TRUE)] <- NA

# Transformer en format long
heatmap_data <- as.data.frame(as.table(dist_mat))

colnames(heatmap_data) <- c("obs1", "obs2", "distance")

# Retirer les NA
heatmap_data <- heatmap_data[!is.na(heatmap_data$distance), ]

# Heatmap
ggplot(
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

# ----------------- Split train/test 1 : random
set.seed(123)
n <- nrow(df)
train_index <- sample(seq_len(n), size = 0.8*n)

train <- df[train_index, ]
test  <- df[-train_index, ]

train$year <- NULL
test$year <- NULL
train$time <- NULL
test$time <- NULL

# scale train and no data leakage
vars_num <- c("chla", "per", "but", "fuco", "hex",
              "allo", "zea", "chlb", "ftle")
train_scaled <- scale(train[, vars_num])
train_center <- attr(train_scaled, "scaled:center")
train_scale  <- attr(train_scaled, "scaled:scale")
test_scaled <- scale(
  test[, vars_num],
  center = train_center,
  scale = train_scale
)

# remise en forme 
ds_train_scaled <- train
ds_train_scaled[, vars_num] <- train_scaled

ds_test_scaled <- test
ds_test_scaled[, vars_num] <- test_scaled

colSums(is.na(ds_train_scaled)) # OK, pas de NA

# verifier le dataleakage
# Indiquer à quelle partie appartient chaque observation
is_train <- seq_len(nrow(df)) %in% train_index
is_test  <- !is_train

# Garder uniquement les paires TRAIN <-> TEST
train_test_pairs <- outer(
  is_train,
  is_test,
  FUN = "&"
)

# Nombre de paires train-test avec distance < 0.25
sum(
  dist_mat < 0.25 & train_test_pairs,
  na.rm = TRUE
)

pairs <- which(
  dist_mat < 0.25 & train_test_pairs,
  arr.ind = TRUE
)

leakage <- data.frame(
  train = pairs[, 1],
  test = pairs[, 2],
  distance = dist_mat[pairs]
)

leakage$time_train <- df$time[leakage$train]
leakage$time_test  <- df$time[leakage$test]

leakage$time_diff_min <- abs(
  as.numeric(
    difftime(
      leakage$time_train,
      leakage$time_test,
      units = "mins"
    )
  )
)

leakage
print(nrow(ds_test_scaled)) 

# leakage de 19/43 données 


# --------------- Split train/test 2 : leave-one-year-out
# cat("Nombre total de données :", nrow(df), "\n")
# 
# # Nombre d'observations par année
# year_counts <- table(df$year)
# print(year_counts)
# 
# # Année avec le moins d'observations
# test_year <- names(which.min(year_counts))
# cat("Année utilisée pour le test :", test_year, "\n")
# 
# # Split train / test
# train <- subset(df, year != test_year)
# test  <- subset(df, year == test_year)
# 
# train$year <- NULL
# test$year <- NULL
# cat("Train :", nrow(train), "observations\n", nrow(train)/nrow(df)*100, "%")
# cat("Test  :", nrow(test), "observations\n", nrow(test)/nrow(df)*100, "%")

# ----------------- Multivariate Regression Tree
# Train model on train dataset only
model <- rpart(
  NASC ~ ., # entries = all datas but not NASC
  data = ds_train_scaled,
  method = "anova"
)


rpart.plot(
  model,
  extra = 101,
  main = "Regression tree - CART algorithm - RS (80/20) - data scaled" # LOYO (90/10)
)

# cp_table <- printcp(model)
# best_row <- which.min(model$cptable[, "xerror"])

# Importance des variables
importance <- model$variable.importance
importance_pct <- importance / sum(importance) * 100

importance_df <- data.frame( # transfo en df
  variable = names(importance_pct),
  importance = as.numeric(importance_pct)
)

importance_df <- importance_df[order(importance_df$importance), ] # sort par ordre d'importance

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
    title = "Variable importance for NASC prediction - CART algorithm",
    subtitle="Random split (80/20) on scaled data" #  #
  )
# --------------- Test Model

prediction <- predict(
  model,
  newdata = ds_test_scaled
)

observed <- ds_test_scaled$NASC

RMSE <- sqrt(mean((prediction - observed)^2))

R2 <- 1 - sum((observed - prediction)^2) /
  sum((observed - mean(observed))^2)

MAE <- mean(abs(prediction - observed))

print(RMSE)
print(R2)
print(MAE)

# --------------- Results analysis
resultats <- data.frame(
  NASC_reel = ds_test_scaled$NASC,
  NASC_predit = prediction
)

head(resultats)


plot(
  ds_test_scaled$NASC,
  prediction,
  xlab = "True NASC",
  ylab = "Predicted NASC",
  main = "Real NASC vs. Predicted NASC - CART algorithm",
  pch = 4,
  cex = 0.5
)
mtext(
  "Random split (80/20) on scaled datas", # "Split: leave-one-year-out (90/10)",
  side = 3,
  line = 0.5,
  cex = 0.9
)

abline(0, 1, col = "red")

legend(
  "topleft",
  legend = c(paste0("Predictions (RMSE = ", round(RMSE, 2), ", R² =", round(R2, 2), ")"), "Identity line (y = x)"),
  col = c("black", "red"),
  pch = c(4, NA),
  lty = c(NA, 1)
)


# ------------------------------------------------------------ Random forest
model_rf <- randomForest(
  NASC ~ .,
  data = ds_train_scaled,
  ntree = 500,
  importance = TRUE
)

# ---------------------------------------------------------
# Importance des variables
# ---------------------------------------------------------

importance(model_rf)

# %IncMSE = importance par permutation
importance_df <- data.frame(
  variable = rownames(importance(model_rf)),
  importance = importance(model_rf)[, "%IncMSE"]
)

importance_df <- importance_df[
  order(importance_df$importance),
]

# Plot importance
ggplot(
  importance_df,
  aes(
    x = reorder(variable, importance),
    y = importance
  )
) +
  geom_col() +
  coord_flip() +
  theme_bw() +
  labs(
    x = "Variable",
    y = "% Increase in MSE",
    title = "Variable importance for NASC prediction - Random forest algorithm",
    subtitle= "Random split (80/20) on scaled datas" # "split leave-one-year_out (90/10)"
  )


# ---------------------------------------------------------
# Test Model
# ---------------------------------------------------------

prediction <- predict(
  model_rf,
  newdata = ds_test_scaled
)


# ---------------------------------------------------------
# Results analysis
# ---------------------------------------------------------

resultats <- data.frame(
  NASC_reel = ds_test_scaled$NASC,
  NASC_predit = prediction
)

head(resultats)


# ---------------------------------------------------------
# RMSE
# ---------------------------------------------------------

RMSE <- sqrt(
  mean(
    (prediction - ds_test_scaled$NASC)^2
  )
)

print(RMSE)


# ---------------------------------------------------------
# MAE
# ---------------------------------------------------------

MAE <- mean(
  abs(prediction - ds_test_scaled$NASC)
)

print(MAE)


# ---------------------------------------------------------
# R²
# ---------------------------------------------------------

R2 <- 1 -
  sum(
    (ds_test_scaled$NASC - prediction)^2
  ) /
  sum(
    (ds_test_scaled$NASC - mean(ds_test_scaled$NASC))^2
  )

print(R2)


# ---------------------------------------------------------
# Real NASC vs Predicted NASC
# ---------------------------------------------------------

plot(
  ds_test_scaled$NASC,
  prediction,
  xlab = "True NASC",
  ylab = "Predicted NASC",
  main = "Real NASC vs. Predicted NASC - Random forest algorithm",
  pch = 4,
  cex = 0.5
)

abline(
  0,
  1,
  col = "red"
)
mtext(
  "Random split (80/20) on scaled datas", # "Split: leave-one-year-out (90/10)",
  side = 3,
  line = 0.5,
  cex = 0.9
)
legend(
  "topleft",
  legend = c(
    paste0("Predictions (RMSE = ", round(RMSE, 2), ")"),
    paste0("R² = ", round(R2, 2)),
    "Identity line (y = x)"
  ),
  col = c(
    "black",
    "black",
    "red"
  ),
  pch = c(
    4,
    NA,
    NA
  ),
  lty = c(
    NA,
    NA,
    1
  )
)


