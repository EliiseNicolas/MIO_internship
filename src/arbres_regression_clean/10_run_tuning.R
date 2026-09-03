# =====================================================================
# 10_run_tuning.R -- SCRIPT 1 : tuning propre CART / RF / XGB
# =====================================================================
# Pour CHAQUE fréquence (38, 120 kHz) et CHAQUE schéma de split (naive
# RS 80/20, blocage spatial x3 résolutions, blocage temporel x1), on
# tune séparément CART, RF, XGB. Voir la note dans 06_tuning.R pour la
# justification de ce choix (tuning différent par schéma).
#
# Sorties : un fichier RDS par (freq, modele, schema) contenant
# `tuning_results` (toute la grille) et `best_params`, + un plot de la
# grille de tuning. Tout est stocké sous outputs_pipeline/tuning/.

source("R/00_config.R")
source("R/01_data_prep.R")
source("R/02_folds.R")
source("R/03_models.R")
source("R/04_diagnostics.R")
source("R/05_plots.R")
source("R/06_tuning.R")

tuning_dir <- path_out("tuning")
dir.create(tuning_dir, showWarnings = FALSE, recursive = TRUE)

for (freq in FREQS) {

  cat("\n============================================================\n")
  cat("FREQUENCE :", freq, "kHz\n")
  cat("============================================================\n")

  # deux versions du jeu de données : sans NA (CART/RF) et avec NA (XGB)
  prep_complete <- load_and_clean(freq, drop_na_numeric = TRUE)
  prep_xgb      <- load_and_clean(freq, drop_na_numeric = FALSE)

  schemes_complete <- build_all_schemes(prep_complete$df)
  schemes_xgb      <- build_all_schemes(prep_xgb$df)

  fod_levels <- prep_xgb$fod_levels

  for (scheme_name in names(schemes_complete)) {

    cat("\n---- Schema :", scheme_name, "----\n")

    # ---- CART ----
    backend_cart <- make_backend("cart")
    tune_cart <- tune_model(schemes_complete[[scheme_name]], "cart",
                             default_tuning_grid("cart"), backend_cart)
    saveRDS(tune_cart, file.path(tuning_dir, sprintf("cart_%dkHz_%s.rds", freq, scheme_name)))

    p_cart <- plot_validation_curve(tune_cart$tuning_results, "cp", facet_var = "maxdepth",
                                     subtitle = sprintf("CART - %d kHz - %s", freq, scheme_name))
    ggsave(file.path(tuning_dir, sprintf("cart_%dkHz_%s.png", freq, scheme_name)),
           p_cart, width = 8, height = 6, dpi = 150)

    # ---- RF ----
    backend_rf <- make_backend("rf")
    tune_rf <- tune_model(schemes_complete[[scheme_name]], "rf",
                           default_tuning_grid("rf"), backend_rf)
    saveRDS(tune_rf, file.path(tuning_dir, sprintf("rf_%dkHz_%s.rds", freq, scheme_name)))

    p_rf <- plot_validation_curve(tune_rf$tuning_results, "mtry", facet_var = "min.node.size",
                                   subtitle = sprintf("RF - %d kHz - %s", freq, scheme_name))
    ggsave(file.path(tuning_dir, sprintf("rf_%dkHz_%s.png", freq, scheme_name)),
           p_rf, width = 8, height = 6, dpi = 150)

    # ---- XGB (early stopping sur nrounds, folds du schema "avec NA") ----
    tune_xgb <- tune_xgb_early_stopping(
      schemes_xgb[[scheme_name]],
      expand.grid(max_depth = c(3, 4, 6), eta = c(0.01, 0.05, 0.1), min_child_weight = c(1, 3, 5)),
      fod_levels = fod_levels
    )
    saveRDS(tune_xgb, file.path(tuning_dir, sprintf("xgb_%dkHz_%s.rds", freq, scheme_name)))

    p_xgb <- plot_validation_curve(tune_xgb$tuning_results, "eta", facet_var = "max_depth",
                                    subtitle = sprintf("XGB - %d kHz - %s", freq, scheme_name))
    ggsave(file.path(tuning_dir, sprintf("xgb_%dkHz_%s.png", freq, scheme_name)),
           p_xgb, width = 8, height = 6, dpi = 150)

    cat(sprintf("  CART best: cp=%.4f minsplit=%d maxdepth=%d | RMSE=%.4f\n",
                tune_cart$best_params$cp, tune_cart$best_params$minsplit,
                tune_cart$best_params$maxdepth, min(tune_cart$tuning_results$mean_rmse_test)))
    cat(sprintf("  RF   best: mtry=%d min.node.size=%d num.trees=%d | RMSE=%.4f\n",
                tune_rf$best_params$mtry, tune_rf$best_params$min.node.size,
                tune_rf$best_params$num.trees, min(tune_rf$tuning_results$mean_rmse_test)))
    cat(sprintf("  XGB  best: max_depth=%d eta=%.3f min_child_weight=%d nrounds=%d | RMSE=%.4f\n",
                tune_xgb$best_params$max_depth, tune_xgb$best_params$eta,
                tune_xgb$best_params$min_child_weight, tune_xgb$best_params$nrounds,
                min(tune_xgb$tuning_results$mean_rmse_test)))
  }
}

cat("\nTuning termine. Resultats dans :", normalizePath(tuning_dir), "\n")
