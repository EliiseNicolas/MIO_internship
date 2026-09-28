# =====================================================================
# 08_robustness.R -- robustesse au bruit gaussien sur les covariables
# =====================================================================
# Equivalent tabulaire du test "robustesse a la perturbation de texture"
# en vision : au lieu de perturber des pixels, on ajoute un bruit gaussien
# croissant aux covariables NUMERIQUES du jeu de TEST (jamais au train),
# et on mesure a quelle vitesse chaque modele se degrade.
#
# Principe : le modele est entraine UNE SEULE FOIS par fold (sur le train
# non perturbe), puis reevalue sur des copies de plus en plus bruitees du
# meme jeu de test -- pas de reentrainement par niveau de bruit, donc le
# cout de calcul reste raisonnable.
#
# L'amplitude du bruit est exprimee en FRACTION de l'ecart-type de
# chaque covariable (calcule sur le TRAIN, jamais sur le test, pour ne
# pas fuiter d'information) : niveau 0.5 = bruit gaussien d'ecart-type
# egal a 50% de l'ecart-type naturel de la covariable. Ca permet de
# comparer des covariables a des echelles tres differentes (ftle vs un
# ratio *_totpig) sur un pied d'egalite.

# ---------------------------------------------------------------------
# Ajoute un bruit gaussien independant a chaque covariable numerique
# d'un data.frame. `noise_sd_vec` est un vecteur nomme (une valeur
# d'ecart-type de bruit par covariable) -- typiquement
# noise_level * sd(train_df[[v]]), calcule une fois par fold.
# ---------------------------------------------------------------------
add_gaussian_noise <- function(df, covs_num, noise_sd_vec, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  for (v in covs_num) {
    sdv <- noise_sd_vec[[v]]
    if (!is.na(sdv) && sdv > 0) {
      df[[v]] <- df[[v]] + stats::rnorm(nrow(df), mean = 0, sd = sdv)
    }
  }
  df
}

# ---------------------------------------------------------------------
# Pour UN fold : entraine une fois, puis reevalue a plusieurs niveaux de
# bruit croissants sur le test. Retourne une ligne par niveau de bruit.
# ---------------------------------------------------------------------
compute_fold_noise_robustness <- function(fold, data, params, backend,
                                           noise_levels = c(0, 0.1, 0.25, 0.5, 1, 2),
                                           response = RESPONSE_VAR,
                                           covs_model = COVARIATES_ALL,
                                           covs_num = COVARIATES_NUM,
                                           seed = 1) {

  train_df <- prep_fold_data(data, fold$train, backend, response, covs_model)
  test_df  <- prep_fold_data(data, fold$test,  backend, response, covs_model)

  model <- backend$fit(train_df, params)

  # Ecart-type "naturel" de chaque covariable, calcule sur le TRAIN
  # uniquement (jamais sur le test -- eviterait de fuiter de l'info)
  base_sd <- vapply(covs_num, function(v) stats::sd(train_df[[v]], na.rm = TRUE), numeric(1))

  purrr::map_dfr(noise_levels, function(nl) {
    noise_sd_vec <- base_sd * nl
    test_noisy <- add_gaussian_noise(test_df, covs_num, noise_sd_vec, seed = seed)
    pred_noisy <- backend$predict(model, test_noisy)

    tibble(
      noise_level = nl,
      rmse_test   = rmse_fn(test_df[[response]], pred_noisy),
      r2_test     = 1 - sum((test_df[[response]] - pred_noisy)^2, na.rm = TRUE) /
        sum((test_df[[response]] - mean(train_df[[response]], na.rm = TRUE))^2, na.rm = TRUE)
    )
  })
}

# ---------------------------------------------------------------------
# Boucle sur tous les folds d'un schema -- meme convention que
# run_cv_scheme() dans 04_diagnostics.R.
# ---------------------------------------------------------------------
run_noise_robustness_scheme <- function(scheme, params, backend,
                                         noise_levels = c(0, 0.1, 0.25, 0.5, 1, 2),
                                         label = "") {
  folds <- scheme$folds
  data  <- scheme$data
  cat(sprintf("  -> robustesse au bruit : %d folds x %d niveaux (%s)\n",
              length(folds), length(noise_levels), label))

  imap_dfr(folds, function(f, fid) {
    compute_fold_noise_robustness(f, data, params, backend, noise_levels = noise_levels) %>%
      mutate(fold_id = fid)
  })
}
