# =====================================================================
# 20_run_regularization_stress_test.R -- SCRIPT 10 : stress-test régul.
# =====================================================================
# Repond a la question "mes hyperparametres sont-ils trop severement
# regularises ?" en comparant, pour chaque freq x modele x schema :
#   - "actuel"  : les hyperparametres tunes par 10_run_tuning.R (deja
#                 sur disque, aucun re-tuning ici)
#   - "extreme" : des hyperparametres delibrement poussés au maximum de
#                 flexibilite (quasi pas de regularisation) --
#                 cp=0/minsplit=2/maxdepth=30 (CART), mtry=max/
#                 min.node.size=1 (RF), max_depth=10/min_child_weight=1/
#                 eta=0.3 (XGB)
#
# Un seul entrainement par fold pour chaque regime (pas un vrai tuning,
# juste UN point de comparaison) -- resultat : R2 train (le modele
# arrive-t-il a fitter ses propres donnees sans contrainte ?) et RMSE
# test (est-ce que cette flexibilite supplementaire AIDE ou NUIT ?).
#
# Grille de lecture (verdict automatique, colonne `verdict`) :
#   - delta R2 train faible, meme en extreme -> PLAFOND D'INFORMATION :
#     ce n'est pas un probleme d'hyperparametres, le signal disponible
#     dans les covariables (ou le bruit de NASC lui-meme) est le
#     facteur limitant, pas la regularisation.
#   - delta R2 train fort ET RMSE test s'ameliore en extreme -> SOUS-
#     REGULARISE : le tuning actuel n'explore pas assez loin cote
#     flexible, elargir la grille (06_tuning.R::default_tuning_grid).
#   - delta R2 train fort MAIS RMSE test se degrade en extreme ->
#     REGULARISATION ACTUELLE JUSTIFIEE (overfit confirme des qu'on
#     relache la contrainte).
#
# Sorties, sous outputs_pipeline/regularization_stress_test/ :
#   - stress_test_results.csv
#   - stress_test_r2_train.png, stress_test_rmse_test.png (par frequence)

source("R/00_config.R")
source("R/01_data_prep.R")
source("R/02_folds.R")
source("R/03_models.R")
source("R/04_diagnostics.R")
source("R/05_plots.R")

tuning_dir <- path_out("tuning")
out_root   <- path_out("regularization_stress_test")
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

# Schemas testes : naive (reference) + tous les schemas bloques par
# defaut -- reduis cette liste si tu veux juste un aller-retour rapide
# (ex. STRESS_TEST_SCHEMES <- c("naive_RS_80_20", "blocked_temporal_1j")).
STRESS_TEST_SCHEMES <- c(
  "naive_RS_80_20",
  paste0("blocked_spatial_", map_chr(SPATIAL_RESOLUTIONS, "label")),
  paste0("blocked_temporal_", map_chr(TEMPORAL_RESOLUTIONS, "label"))
)
MODELS <- c("cart", "rf", "xgb")

# Seuils du verdict automatique -- ajustables si les resultats sont
# systematiquement dans la zone "ambigu".
DELTA_R2_TRAIN_THRESHOLD  <- 0.10   # a partir de quand un gain de R2 train est "significatif"
DELTA_RMSE_TEST_THRESHOLD <- 0.01   # a partir de quand un changement de RMSE test est "significatif"

extreme_params <- function(model) {
  switch(model,
    "cart" = list(cp = 0, minsplit = 2, maxdepth = 30),
    "rf"   = list(mtry = length(COVARIATES_ALL), min.node.size = 1, num.trees = 500),
    "xgb"  = list(max_depth = 10, eta = 0.3, min_child_weight = 1, nrounds = 500,
                  subsample = 0.8, colsample_bytree = 0.8),
    stop("model inconnu : ", model)
  )
}

compute_verdict <- function(delta_r2_train, delta_rmse_test) {
  if (delta_r2_train < DELTA_R2_TRAIN_THRESHOLD) {
    "Plafond d'information (pas un probleme d'hyperparametres)"
  } else if (delta_rmse_test < -DELTA_RMSE_TEST_THRESHOLD) {
    "SOUS-REGULARISE : elargir la grille de tuning"
  } else if (delta_rmse_test > DELTA_RMSE_TEST_THRESHOLD) {
    "Regularisation actuelle justifiee (overfit confirme en mode extreme)"
  } else {
    "Ambigu : train ameliore mais test quasi inchange"
  }
}

results <- list()

for (freq in FREQS) {

  cat("\n============================================================\n")
  cat("STRESS-TEST REGULARISATION -- FREQUENCE :", freq, "kHz\n")
  cat("============================================================\n")

  prep_complete <- load_and_clean(freq, drop_na_numeric = TRUE)   # pour CART/RF
  prep_xgb      <- load_and_clean(freq, drop_na_numeric = FALSE)  # pour XGB
  fod_levels    <- prep_xgb$fod_levels

  schemes_rf  <- build_all_schemes(prep_complete$df)
  schemes_xgb <- build_all_schemes(prep_xgb$df)

  for (model in MODELS) {
    backend <- if (model == "xgb") make_backend("xgb", fod_levels = fod_levels) else make_backend(model)
    schemes <- if (model == "xgb") schemes_xgb else schemes_rf

    for (scheme_name in STRESS_TEST_SCHEMES) {

      tuning_path <- file.path(tuning_dir, sprintf("%s_%dkHz_%s.rds", model, freq, scheme_name))
      if (!file.exists(tuning_path)) {
        cat("  [!] tuning introuvable pour", model, scheme_name, "-- saute.\n")
        next
      }
      params_actuel  <- readRDS(tuning_path)$best_params
      params_extreme <- extreme_params(model)

      scheme <- schemes[[scheme_name]]
      cat(" ", model, "-", scheme_name, "\n")

      cv_actuel  <- run_cv_scheme(scheme, params_actuel,  backend, label = "actuel")
      cv_extreme <- run_cv_scheme(scheme, params_extreme, backend, label = "extreme")

      r2_train_actuel   <- mean(cv_actuel$metrics$r2_train,  na.rm = TRUE)
      r2_train_extreme  <- mean(cv_extreme$metrics$r2_train, na.rm = TRUE)
      rmse_test_actuel  <- mean(cv_actuel$metrics$rmse_test,  na.rm = TRUE)
      rmse_test_extreme <- mean(cv_extreme$metrics$rmse_test, na.rm = TRUE)

      delta_r2_train  <- r2_train_extreme  - r2_train_actuel
      delta_rmse_test <- rmse_test_extreme - rmse_test_actuel

      verdict <- compute_verdict(delta_r2_train, delta_rmse_test)

      results[[length(results) + 1]] <- tibble(
        freq = freq, model = model, scheme = scheme_name,
        r2_train_actuel = r2_train_actuel, r2_train_extreme = r2_train_extreme,
        delta_r2_train = delta_r2_train,
        rmse_test_actuel = rmse_test_actuel, rmse_test_extreme = rmse_test_extreme,
        delta_rmse_test = delta_rmse_test,
        verdict = verdict
      )

      cat(sprintf("    R2 train : %.3f -> %.3f (%+.3f) | RMSE test : %.3f -> %.3f (%+.3f)\n",
                  r2_train_actuel, r2_train_extreme, delta_r2_train,
                  rmse_test_actuel, rmse_test_extreme, delta_rmse_test))
      cat("    Verdict :", verdict, "\n")
    }
  }
}

stress_df <- bind_rows(results)
write.csv(stress_df, file.path(out_root, "stress_test_results.csv"), row.names = FALSE)

for (freq in FREQS) {
  sub <- stress_df %>% filter(freq == !!freq)
  if (nrow(sub) == 0) next

  p_r2 <- plot_stress_test(sub, metric = "r2_train", subtitle = sprintf("%d kHz", freq))
  ggsave(file.path(out_root, sprintf("stress_test_r2_train_%dkHz.png", freq)), p_r2, width = 11, height = 6, dpi = 150)

  p_rmse <- plot_stress_test(sub, metric = "rmse_test", subtitle = sprintf("%d kHz", freq))
  ggsave(file.path(out_root, sprintf("stress_test_rmse_test_%dkHz.png", freq)), p_rmse, width = 11, height = 6, dpi = 150)
}

cat("\n============================================================\n")
cat("RESUME DES VERDICTS\n")
cat("============================================================\n")
print(stress_df %>% select(freq, model, scheme, delta_r2_train, delta_rmse_test, verdict))

cat("\nStress-test termine. Resultats dans :", normalizePath(out_root), "\n")
