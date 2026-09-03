# =====================================================================
# 04_diagnostics.R -- diagnostics par fold, communs à tous les modèles
# =====================================================================
# Le même code produit TOUTES les métriques et TOUTES les tables dont on
# a besoin pour les plots (peu importe CART/RF/XGB, naive/blocked) :
#  - métriques de performance (RMSE, R², par fold)
#  - distance géographique test -> train le plus proche
#  - distance euclidienne (covariables standardisées) test -> train
#  - stats des covariables numériques (moyenne/sd/variance) train vs test
#  - distribution de fod (train vs test)
#  - distribution complète de chaque covariable numérique (train vs test)
#  - table obs/pred/résidus spatialisée (pour les cartes)
#  - importance des variables

prep_fold_data <- function(data, idx, backend, response = RESPONSE_VAR, covs = COVARIATES_ALL) {
  d <- data[idx, ]
  if (isTRUE(backend$needs_complete)) {
    d <- d[stats::complete.cases(d[, c(response, covs)]), ]
  }
  d
}

compute_fold_diagnostics <- function(fold, data, params, backend,
                                      response   = RESPONSE_VAR,
                                      covs_model = COVARIATES_ALL,
                                      covs_num   = COVARIATES_NUM) {

  train_df <- prep_fold_data(data, fold$train, backend, response, covs_model)
  test_df  <- prep_fold_data(data, fold$test,  backend, response, covs_model)

  model <- backend$fit(train_df, params)

  pred_test  <- backend$predict(model, test_df)
  pred_train <- backend$predict(model, train_df)

  obs_test   <- test_df[[response]]
  obs_train  <- train_df[[response]]
  train_mean <- mean(obs_train, na.rm = TRUE)
  residuals  <- obs_test - pred_test

  metrics <- tibble(
    n_train    = nrow(train_df), n_test = nrow(test_df),
    rmse_test  = rmse_fn(obs_test, pred_test),
    rmse_train = rmse_fn(obs_train, pred_train),
    r2_test    = 1 - sum((obs_test  - pred_test)^2,  na.rm = TRUE) / sum((obs_test  - train_mean)^2, na.rm = TRUE),
    r2_train   = 1 - sum((obs_train - pred_train)^2, na.rm = TRUE) / sum((obs_train - train_mean)^2, na.rm = TRUE),
    var_intra_fold_test = var(obs_test, na.rm = TRUE),
    sd_intra_fold_test  = sd(obs_test, na.rm = TRUE),
    mean_residual = mean(residuals, na.rm = TRUE),
    sd_residual   = sd(residuals, na.rm = TRUE)
  )

  # ---- distance géographique : chaque test -> train le plus proche (km) ----
  # Moyenne ET écart-type/variance de cette distance A L'INTERIEUR du fold
  # (dispersion entre les points de test, pas seulement une moyenne unique).
  geo_dist <- FNN::knnx.dist(
    data  = as.matrix(train_df[, c("x_km", "y_km")]),
    query = as.matrix(test_df[,  c("x_km", "y_km")]), k = 1
  )[, 1]
  metrics$mean_geo_dist_km <- mean(geo_dist)
  metrics$sd_geo_dist_km   <- sd(geo_dist)
  metrics$var_geo_dist_km  <- var(geo_dist)

  # ---- distance euclidienne dans l'espace des covariables numériques ----
  # (standardisées sur le train) -- même logique : moyenne + dispersion.
  train_num <- as.matrix(train_df[, covs_num])
  test_num  <- as.matrix(test_df[,  covs_num])
  center <- colMeans(train_num, na.rm = TRUE)
  scale_ <- apply(train_num, 2, sd, na.rm = TRUE)
  train_num_s <- scale(train_num, center = center, scale = scale_)
  test_num_s  <- scale(test_num,  center = center, scale = scale_)
  train_nn <- FNN::knn.dist(train_num_s, k = 1)[, 1]
  test_nn  <- FNN::knnx.dist(train_num_s, test_num_s, k = 1)[, 1]
  metrics$mean_covariate_dist <- mean(test_nn)
  metrics$sd_covariate_dist   <- sd(test_nn)
  metrics$var_covariate_dist  <- var(test_nn)
  metrics$extrapolation_index <- mean(test_nn) / mean(train_nn)

  # ---- variance / moyenne / sd des covariables numériques, train vs test ----
  covariate_stats <- bind_rows(
    tibble(set = "test",  variable = covs_num,
           mean = colMeans(test_num,  na.rm = TRUE), sd = apply(test_num,  2, sd,  na.rm = TRUE),
           var  = apply(test_num,  2, var, na.rm = TRUE)),
    tibble(set = "train", variable = covs_num,
           mean = colMeans(train_num, na.rm = TRUE), sd = apply(train_num, 2, sd,  na.rm = TRUE),
           var  = apply(train_num, 2, var, na.rm = TRUE))
  )

  # ---- distribution de fod, train vs test ----
  fod_dist <- bind_rows(
    tibble(fod = train_df$fod, set = "train"),
    tibble(fod = test_df$fod,  set = "test")
  ) %>%
    filter(!is.na(fod)) %>%
    count(set, fod, .drop = FALSE) %>%
    group_by(set) %>% mutate(prop = n / sum(n)) %>% ungroup()

  # ---- distribution complète de chaque covariable numérique, train vs test --
  numeric_dist <- bind_rows(
    pivot_longer(train_df[, covs_num, drop = FALSE], everything(),
                 names_to = "variable", values_to = "value") %>% mutate(set = "train"),
    pivot_longer(test_df[,  covs_num, drop = FALSE], everything(),
                 names_to = "variable", values_to = "value") %>% mutate(set = "test")
  )

  # ---- table obs/pred spatialisée (pour cartes de résidus) ----
  obs_pred <- tibble(
    obs = obs_test, pred = pred_test, residual = residuals,
    lon = test_df$lon, lat = test_df$lat
  )

  importance_df <- backend$importance(model)

  list(
    metrics = metrics, covariate_stats = covariate_stats, fod_dist = fod_dist,
    numeric_dist = numeric_dist, obs_pred = obs_pred, importance = importance_df,
    model = model, train_idx = fold$train, test_idx = fold$test
  )
}

# ---------------------------------------------------------------------
# Boucle sur tous les folds d'un schéma, avec un jeu d'hyperparamètres
# fixé -- utilisée à la fois pendant le tuning (pour évaluer une
# combinaison d'hyperparamètres) et pendant l'entraînement final (avec
# les meilleurs hyperparamètres trouvés).
# ---------------------------------------------------------------------
run_cv_scheme <- function(scheme, params, backend, label = "") {
  folds <- scheme$folds
  data  <- scheme$data
  cat(sprintf("  -> %d folds (%s)\n", length(folds), label))

  results <- imap(folds, ~ compute_fold_diagnostics(.x, data, params, backend))

  list(
    label           = label,
    scheme          = scheme,
    results         = results,
    metrics         = imap_dfr(results, ~ mutate(.x$metrics, fold_id = .y)),
    obs_pred        = imap_dfr(results, ~ mutate(.x$obs_pred, fold_id = .y)),
    covariate_stats = imap_dfr(results, ~ mutate(.x$covariate_stats, fold_id = .y)),
    fod_dist        = imap_dfr(results, ~ mutate(.x$fod_dist, fold_id = .y)),
    numeric_dist    = imap_dfr(results, ~ mutate(.x$numeric_dist, fold_id = .y)),
    importance      = imap_dfr(results, ~ mutate(.x$importance, fold_id = .y)),
    models          = map(results, "model")
  )
}

# ---------------------------------------------------------------------
# Courbe d'apprentissage : sous-échantillonnage progressif du train de
# chaque fold, réentraînement, RMSE train/test en fonction de la
# fraction de train utilisée.
# ---------------------------------------------------------------------
compute_learning_curve <- function(scheme, params, backend,
                                    fractions = seq(0.2, 1, by = 0.2),
                                    response = RESPONSE_VAR, covs_model = COVARIATES_ALL,
                                    seed = 1) {
  folds <- scheme$folds
  data  <- scheme$data
  set.seed(seed)

  out <- imap_dfr(folds, function(f, fid) {
    map_dfr(fractions, function(frac) {
      n_sub   <- max(20, round(frac * length(f$train)))
      sub_idx <- sample(f$train, min(n_sub, length(f$train)))

      train_df <- prep_fold_data(data, sub_idx, backend, response, covs_model)
      test_df  <- prep_fold_data(data, f$test,  backend, response, covs_model)

      model <- backend$fit(train_df, params)
      pred_test  <- backend$predict(model, test_df)
      pred_train <- backend$predict(model, train_df)

      tibble(
        fraction   = frac,
        n_train    = nrow(train_df),
        rmse_train = rmse_fn(train_df[[response]], pred_train),
        rmse_test  = rmse_fn(test_df[[response]],  pred_test)
      )
    }) %>% mutate(fold_id = fid)
  })

  summary <- out %>%
    group_by(fraction) %>%
    summarise(
      mean_rmse_train = mean(rmse_train), sd_rmse_train = sd(rmse_train),
      mean_rmse_test  = mean(rmse_test),  sd_rmse_test  = sd(rmse_test),
      .groups = "drop"
    )

  list(detail = out, summary = summary)
}
