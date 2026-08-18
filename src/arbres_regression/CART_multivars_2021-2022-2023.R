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
  fod = as.factor(datas$fod),
  chla = datas$Chla,
  per = datas$Per,
  but = datas$But,
  fuco = datas$Fuco,
  hex = datas$Hex,
  allo = datas$Allo,
  zea = datas$Zea,
  chlb = datas$Chlb,
  dvchla = datas$DvChla,
  ftle = datas$ftle
)
print(sum(!is.na(df$chla)))
print(sum(!is.na(df$fod)))
print(c("Nombre total de données", nrow(df)))

# garder uniquement les données sans na
df <- df[!is.na(df$chla) & !is.na(df$fod) & !is.na(df$ftle) & !is.na(df$per) & !is.na(df$allo),]
print(c("Nombre total de données", nrow(df))) # 215
colSums(is.na(df))

# ----------------- Split train/test 1 : random
set.seed(123)
n <- nrow(df)
train_index <- sample(seq_len(n), size = 0.8*n)

train <- df[train_index, ]
test  <- df[-train_index, ]
train$year <- NULL
test$year <- NULL

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
  data = train,
  method = "anova"
)


rpart.plot(
  model,
  extra = 101,
  main = "Regression tree - CART algorithm - LOYO(90/10)" # RS (80/20)
)

cp_table <- printcp(model)
best_row <- which.min(model$cptable[, "xerror"])
R2_approx <- 1 - model$cptable[best_row, "xerror"]
print(R2_approx)

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
    subtitle="split leave-one-year_out (90/10)" # "Random split (80/20)" #
  )
# --------------- Test Model

prediction <- predict(
  model,
  newdata = test
)

# --------------- Results analysis
resultats <- data.frame(
  NASC_reel = test$NASC,
  NASC_predit = prediction
)

head(resultats)

RMSE <- sqrt(mean((prediction - test$NASC)^2))
print(RMSE)

plot(
  test$NASC,
  prediction,
  xlab = "True NASC",
  ylab = "Predicted NASC",
  main = "Real NASC vs. Predicted NASC - CART algorithm",
  pch = 4,
  cex = 0.5
)
mtext(
  "split leave-one-year_out (90/10)", # "Random split (80/20)", # "Split: leave-one-year-out (90/10)",
  side = 3,
  line = 0.5,
  cex = 0.9
)

abline(0, 1, col = "red")

legend(
  "topleft",
  legend = c(paste0("Predictions (RMSE = ", round(RMSE, 2), ")"), "Identity line (y = x)"),
  col = c("black", "red"),
  pch = c(4, NA),
  lty = c(NA, 1)
)
print(model)

# ------------------------------------------------------------ Random forest
model_rf <- randomForest(
  NASC ~ .,
  data = train,
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
    subtitle= "split leave-one-year_out (90/10)"# "Random split (80/20)" # "split leave-one-year_out (90/10)"
  )


# ---------------------------------------------------------
# Test Model
# ---------------------------------------------------------

prediction <- predict(
  model_rf,
  newdata = test
)


# ---------------------------------------------------------
# Results analysis
# ---------------------------------------------------------

resultats <- data.frame(
  NASC_reel = test$NASC,
  NASC_predit = prediction
)

head(resultats)


# ---------------------------------------------------------
# RMSE
# ---------------------------------------------------------

RMSE <- sqrt(
  mean(
    (prediction - test$NASC)^2
  )
)

print(RMSE)


# ---------------------------------------------------------
# MAE
# ---------------------------------------------------------

MAE <- mean(
  abs(prediction - test$NASC)
)

print(MAE)


# ---------------------------------------------------------
# R²
# ---------------------------------------------------------

R2 <- 1 -
  sum(
    (test$NASC - prediction)^2
  ) /
  sum(
    (test$NASC - mean(test$NASC))^2
  )

print(R2)


# ---------------------------------------------------------
# Real NASC vs Predicted NASC
# ---------------------------------------------------------

plot(
  test$NASC,
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
  "split leave-one-year_out (90/10)", # "Random split (80/20)", # "Split: leave-one-year-out (90/10)",
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


