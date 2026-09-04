# =====================================================================
# test_fix_xgb_tuning.R -- verification rapide du fix xgb.cv
# =====================================================================
# Teste tune_xgb_early_stopping() sur UNE SEULE combinaison (38 kHz,
# naive) avec une grille reduite a 2 combinaisons, pour verifier que le
# fix (colonnes evaluation_log retrouvees par motif) fonctionne, sans
# attendre tout le tuning complet (30 combinaisons x 5 schemas x 2 freq).
# Ne touche a AUCUN fichier de outputs_pipeline/tuning/ -- juste un test
# en memoire, affiche a l'ecran.

source("R/00_config.R")
source("R/01_data_prep.R")
source("R/02_folds.R")
source("R/03_models.R")
source("R/04_diagnostics.R")
source("R/06_tuning.R")

cat("Chargement des donnees 38 kHz...\n")
prep_xgb <- load_and_clean(38, drop_na_numeric = FALSE)
scheme   <- build_naive_folds(prep_xgb$df, n_folds = 3)  # 3 folds seulement, pour aller vite

grille_test <- expand.grid(max_depth = c(3, 6), eta = 0.1, min_child_weight = 1)  # 2 combinaisons

cat("\nLancement du test (2 combinaisons, 3 folds) ...\n")
resultat <- tune_xgb_early_stopping(scheme, grille_test, fod_levels = prep_xgb$fod_levels)

cat("\n=== SUCCES : voici le resultat ===\n")
print(resultat$tuning_results)
cat("\nMeilleurs parametres trouves :\n")
print(resultat$best_params)

cat("\nSi tu vois un tableau avec les colonnes mean_rmse_test/sd_rmse_test/",
    "mean_rmse_train remplies (pas d'erreur), le fix fonctionne -- tu peux",
    "relancer 10_run_tuning.R (ou run_all_pipeline.R) en confiance.\n")

