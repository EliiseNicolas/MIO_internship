# =====================================================================
# TUNING SIMPLE DU RANDOM FOREST (mtry x nodesize), via CV à 10 folds
# =====================================================================

library(dplyr)
library(purrr)
library(ggplot2)

set.seed(123)
n <- nrow(df)
fold_id <- sample(rep(seq_len(n_CV), length.out = n))

# Grille de tuning : petite mais couvre les valeurs usuelles
tuning_grid <- expand.grid(
  mtry     = c(2, 3, 4, floor(sqrt(length(covariates_all)))) |> unique(),
  nodesize = c(1, 5, 10),
  ntree    = 300
)

print(tuning_grid)

# ---------------------------------------------------------------------
# Fonction d'évaluation d'un jeu d'hyperparamètres, sur les 10 folds
# ---------------------------------------------------------------------
evaluate_params <- function(mtry_val, nodesize_val, ntree_val, data, fold_id, n_CV,
                            response, covs) {
  
  fold_perf <- map_dfr(seq_len(n_CV), function(k) {
    
    train_idx <- which(fold_id != k)
    test_idx  <- which(fold_id == k)
    
    train_df <- data[train_idx, c(response, covs)]
    test_df  <- data[test_idx,  c(response, covs)]
    
    model_k <- randomForest(
      stats::as.formula(paste(response, "~ .")),
      data = train_df,
      ntree = ntree_val,
      mtry = mtry_val,
      nodesize = nodesize_val
    )
    
    pred_test <- predict(model_k, newdata = test_df)
    obs_test  <- test_df[[response]]
    
    tibble(
      fold      = k,
      rmse_test = sqrt(mean((pred_test - obs_test)^2)),
      r2_test   = 1 - sum((obs_test - pred_test)^2) / sum((obs_test - mean(train_df[[response]]))^2)
    )
  })
  
  tibble(
    mtry           = mtry_val,
    nodesize       = nodesize_val,
    ntree          = ntree_val,
    mean_rmse_test = mean(fold_perf$rmse_test),
    sd_rmse_test   = sd(fold_perf$rmse_test),
    mean_r2_test   = mean(fold_perf$r2_test),
    sd_r2_test     = sd(fold_perf$r2_test)
  )
}

# ---------------------------------------------------------------------
# Boucle sur toute la grille
# ---------------------------------------------------------------------
cat("\n=== Tuning en cours (", nrow(tuning_grid), "combinaisons x", n_CV, "folds) ===\n")

tuning_results <- pmap_dfr(tuning_grid, function(mtry, nodesize, ntree) {
  cat("  mtry =", mtry, "| nodesize =", nodesize, "| ntree =", ntree, "\n")
  evaluate_params(mtry, nodesize, ntree, df, fold_id, n_CV, response_var, covariates_all)
})

tuning_results <- tuning_results %>% arrange(mean_rmse_test)

cat("\n=== Résultats du tuning (triés par RMSE croissant) ===\n")
print(tuning_results)

best_params <- tuning_results[1, ]
cat("\nMeilleurs hyperparamètres : mtry =", best_params$mtry,
    "| nodesize =", best_params$nodesize,
    "| ntree =", best_params$ntree, "\n")
cat("RMSE moyen :", round(best_params$mean_rmse_test, 4), "\n")
cat("R² moyen   :", round(best_params$mean_r2_test, 4), "\n")

# ---------------------------------------------------------------------
# Visualisation : RMSE selon mtry et nodesize
# ---------------------------------------------------------------------
ggplot(tuning_results, aes(x = factor(mtry), y = mean_rmse_test, fill = factor(nodesize))) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(
    aes(ymin = mean_rmse_test - sd_rmse_test, ymax = mean_rmse_test + sd_rmse_test),
    position = position_dodge(width = 0.8), width = 0.2
  ) +
  labs(
    title = "Tuning du Random Forest (10-fold CV)",
    x = "mtry", y = "RMSE moyen (± sd inter-fold)", fill = "nodesize"
  ) +
  theme_bw()

ggplot(tuning_results, aes(x = factor(mtry), y = mean_r2_test, fill = factor(nodesize))) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(
    aes(ymin = mean_r2_test - sd_r2_test, ymax = mean_r2_test + sd_r2_test),
    position = position_dodge(width = 0.8), width = 0.2
  ) +
  labs(
    title = "Tuning du Random Forest (10-fold CV)",
    x = "mtry", y = "R² moyen (± sd inter-fold)", fill = "nodesize"
  ) +
  theme_bw()

# ---------------------------------------------------------------------
# Modèle final avec les meilleurs hyperparamètres, entraîné sur toutes
# les données
# ---------------------------------------------------------------------
model <- randomForest(
  NASC ~ .,
  data = df[, c(response_var, covariates_all)],
  ntree = best_params$ntree,
  mtry = best_params$mtry,
  nodesize = best_params$nodesize,
  importance = TRUE
)

cat("\nModèle final entraîné avec mtry =", best_params$mtry,
    ", nodesize =", best_params$nodesize, ", ntree =", best_params$ntree, "\n")