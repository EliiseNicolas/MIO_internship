# LIbraries
library(rpart)
library(rpart.plot)
library(ggplot2)
rm(list = ls())

# -------------- Paths
path_datas <- "/run/media/mmolinet/KER22/données elise/elisou_ta_stagiaire_pref/prepross/regression_ds_2021_2022_2023_m100.rds"
datas <- readRDS(path_datas)
print(length(unique(as.Date(datas$time_nasc))))
print(unique(as.Date(datas$time_nasc)))

# filtrer par diurnal period
diurnal_periods <- c("night", "sunrise", "day", "sunset")
period <- "day"
period_idx <- which(diurnal_periods == period)
datas_filtered <- datas[datas$day == period_idx,]
str(datas_filtered) # OK
print(length(unique(as.Date(datas_filtered$time_nasc))))
# moyenner par jour
# datas_filtered$date <- as.Date(datas_filtered$time_nasc)
# 
# daily <- aggregate(
#   cbind(
#     nasc,
#     Chla, Per, But, Fuco, Hex, Allo, Zea, Chlb, DvChla,
#     ftle
#   ) ~ date,
#   data = datas_filtered,
#   FUN = function(x) mean(x, na.rm = TRUE)
# )
# 
# str(daily)

# library(dplyr)

daily <- datas_filtered %>%
  mutate(date = as.Date(time_nasc)) %>%
  group_by(date) %>%
  summarise(
    across(
      c(nasc, fod, Chla, Per, But, Fuco, Hex, Allo, Zea, Chlb, DvChla, ftle),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

# ------------- Keep only relevant variables
df <- data.frame(
  NASC = daily$nasc,
  year = format(daily$date, "%Y"),
  fod = daily$fod,
  chla = daily$Chla,
  per = daily$Per,
  but = daily$But,
  fuco = daily$Fuco,
  hex = daily$Hex,
  allo = daily$Allo,
  zea = daily$Zea,
  chlb = daily$Chlb,
  dvchla = daily$DvChla,
  ftle = daily$ftle
)
print(sum(!is.na(df$chla)))
print(df)
print(c("Nombre total de données", nrow(df)))

# ----------------- Split train/test 1 : random
set.seed(123)
n <- nrow(df)
train_index <- sample(seq_len(n), size = 0.8*n)

train <- df[train_index, ]
test  <- df[-train_index, ]
train$year <- NULL
test$year <- NULL

# --------------- Split train/test 2 : leave-one-year-out
cat("Nombre total de données :", nrow(df), "\n")

# Nombre d'observations par année
year_counts <- table(df$year)
print(year_counts)

# Année avec le moins d'observations
test_year <- names(which.min(year_counts))
cat("Année utilisée pour le test :", test_year, "\n")

# Split train / test
train <- subset(df, year != test_year)
test  <- subset(df, year == test_year)

train$year <- NULL
test$year <- NULL
cat("Train :", nrow(train), "observations\n", nrow(train)/nrow(df)*100, "%")
cat("Test  :", nrow(test), "observations\n", nrow(test)/nrow(df)*100, "%")

# ----------------- Multivariate Regression Tree
# Train model on train dataset only
model <- rpart(
  NASC ~ ., # entries = all datas but not NASC
  data = train,
  method = "anova"
)

rpart.plot(model, extra=101)

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
    title = "Variable importance for NASC prediction"
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
  main = "Real NASC vs. Predicted NASC",
  pch = 4,
  cex = 0.5
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
