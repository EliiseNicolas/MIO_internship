# =====================================================================
# 06_tuning.R -- tuning générique, réutilisé pour CART / RF / XGB
# =====================================================================
# PRINCIPE IMPORTANT (répond à la question "faut-il un tuning différent
# par type de split ?") :
#
#   OUI. Les hyperparamètres doivent être choisis avec le MÊME schéma de
#   validation que celui utilisé ensuite pour évaluer/déployer le modèle.
#   Un jeu d'hyperparamètres optimisé en CV naive (split aléatoire) est
#   presque toujours trop "souple" (trop profond / trop peu régularisé)
#   pour un contexte de prédiction spatiale ou temporelle réelle, parce
#   que la CV naive laisse fuiter de l'information locale (points très
#   proches en espace/temps dans le train ET le test). A l'inverse, des
#   hyperparamètres tunés en blocage spatial/temporel sont probablement
#   plus conservateurs (plus régularisés) -- ce qui est justement ce
#   qu'on veut pour une vraie extrapolation spatiale.
#
#   Idéalement on retune même SEPAREMENT pour chaque résolution de
#   blocage (20x20km n'a pas le même volume de train disponible par
#   fold que 1500x1000km -- le nombre optimal d'arbres / la profondeur
#   raisonnable ne sont pas les mêmes). C'est ce que fait ce pipeline :
#   `tune_model()` est appelée une fois par (fréquence x modèle x schéma).
#
# Ce fichier ne fait QUE le tuning ; l'entraînement final avec les
# meilleurs paramètres se fait dans 07_train_final.R, sur les MEMES
# folds (même objet `scheme`) pour éviter toute incohérence.

tune_model <- function(scheme, model_type, grid, backend,
                        response = RESPONSE_VAR, covs_model = COVARIATES_ALL) {

  folds <- scheme$folds
  data  <- scheme$data

  cat(sprintf("Tuning %s sur schema '%s' (%d folds, %d combinaisons)\n",
              model_type, scheme$scheme %||% "?", length(folds), nrow(grid)))

  tuning_results <- purrr::pmap_dfr(grid, function(...) {
    params <- list(...)
    fold_perf <- imap_dfr(folds, function(f, fid) {
      train_df <- prep_fold_data(data, f$train, backend, response, covs_model)
      test_df  <- prep_fold_data(data, f$test,  backend, response, covs_model)

      model <- backend$fit(train_df, params)
      pred_test  <- backend$predict(model, test_df)
      pred_train <- backend$predict(model, train_df)

      tibble(
        fold_id    = fid,
        rmse_test  = rmse_fn(test_df[[response]],  pred_test),
        rmse_train = rmse_fn(train_df[[response]], pred_train)
      )
    })
    fold_perf %>%
      summarise(
        mean_rmse_test  = mean(rmse_test),  sd_rmse_test  = sd(rmse_test),
        mean_rmse_train = mean(rmse_train), sd_rmse_train = sd(rmse_train)
      ) %>%
      bind_cols(as_tibble(params))
  })

  tuning_results <- tuning_results %>% arrange(mean_rmse_test)
  best_params <- as.list(tuning_results[1, names(grid)])

  list(tuning_results = tuning_results, best_params = best_params, model_type = model_type)
}

# ---------------------------------------------------------------------
# Grilles de tuning par défaut -- à affiner selon les premiers résultats
# ---------------------------------------------------------------------
default_tuning_grid <- function(model_type) {
  switch(model_type,
    "cart" = expand.grid(
      cp       = c(0.0005, 0.001, 0.005, 0.01),
      minsplit = c(10, 20, 40),
      maxdepth = c(5, 10, 15)
    ),
    "rf" = expand.grid(
      mtry          = c(2, 3, 4),
      min.node.size = c(5, 10, 20),
      num.trees     = c(300, 500)
    ),
    "xgb" = expand.grid(
      max_depth        = c(3, 4, 6),
      eta              = c(0.01, 0.05, 0.1),
      min_child_weight = c(1, 3, 5),
      nrounds          = c(300, 500)
    ),
    "rfsrc" = expand.grid(
      mtry     = c(2, 3, 4),
      nodesize = c(5, 15, 30),
      ntree    = c(300, 500)
    ),
    stop("model_type inconnu : ", model_type)
  )
}

# ---------------------------------------------------------------------
# Variante XGBoost plus efficace : utilise xgb.cv avec les folds
# personnalisés du schéma (spatial/temporel/naive) et l'early stopping
# natif pour trouver `nrounds` optimal, au lieu de le mettre dans la
# grille (plus rapide, et plus proche de la pratique standard XGBoost).
# ---------------------------------------------------------------------
tune_xgb_early_stopping <- function(scheme, grid, fod_levels,
                                     covs_num = COVARIATES_NUM, response = RESPONSE_VAR,
                                     nrounds_max = 1000, early_stopping_rounds = 20) {
  data <- scheme$data
  X <- build_design_matrix(data, covs_num, fod_levels)
  y <- data[[response]]
  dall <- xgboost::xgb.DMatrix(data = X, label = y, missing = NA)

  # xgb.cv attend une liste de vecteurs d'indices de TEST (positions
  # dans `dall`), un par fold
  cv_folds <- map(scheme$folds, "test")

  cat(sprintf("Tuning xgb (early stopping) sur schema '%s' (%d folds, %d combinaisons)\n",
              scheme$scheme %||% "?", length(cv_folds), nrow(grid)))

  tuning_results <- purrr::pmap_dfr(grid, function(max_depth, eta, min_child_weight) {
    params <- list(
      max_depth = max_depth, eta = eta, min_child_weight = min_child_weight,
      subsample = 0.8, colsample_bytree = 0.8,
      objective = "reg:squarederror", eval_metric = "rmse"
    )
    cv_result <- xgboost::xgb.cv(
      params = params, data = dall, nrounds = nrounds_max,
      folds = cv_folds, early_stopping_rounds = early_stopping_rounds, verbose = 0
    )

    # NE PAS utiliser cv_result$best_iteration : selon la version de
    # xgboost, ce champ peut être NULL au premier niveau (il a été
    # déplacé dans une sous-liste $early_stop dans certaines versions
    # recentes, ex. 3.2.x). On calcule donc `best_iter` nous-mêmes,
    # directement a partir de evaluation_log (RMSE test minimal) --
    # robuste quelle que soit la version du package installee.
    log_names <- names(cv_result$evaluation_log)
    find_col <- function(pattern) {
      m <- grep(pattern, log_names, value = TRUE)
      if (length(m) == 0) {
        stop(
          "Colonne introuvable dans evaluation_log pour le motif '", pattern,
          "'. Colonnes disponibles : ", paste(log_names, collapse = ", ")
        )
      }
      m[1]
    }
    col_test_mean  <- find_col("test.*rmse.*mean")
    col_test_std   <- find_col("test.*rmse.*std")
    col_train_mean <- find_col("train.*rmse.*mean")

    best_iter <- which.min(cv_result$evaluation_log[[col_test_mean]])
    best_row  <- cv_result$evaluation_log[best_iter, ]

    tibble(
      max_depth = max_depth, eta = eta, min_child_weight = min_child_weight,
      nrounds = best_iter,
      mean_rmse_test  = best_row[[col_test_mean]],
      sd_rmse_test    = best_row[[col_test_std]],
      mean_rmse_train = best_row[[col_train_mean]]
    )
  })

  tuning_results <- tuning_results %>% arrange(mean_rmse_test)
  best_params <- as.list(tuning_results[1, c("max_depth", "eta", "min_child_weight", "nrounds")])

  list(tuning_results = tuning_results, best_params = best_params, model_type = "xgb")
}
