# =====================================================================
# 11_run_training.R -- SCRIPT 2 : entraînement 10 folds RF & XGB
# =====================================================================
# Pour chaque fréquence (38, 120 kHz), pour RF et XGB, pour :
#   - naive RS 80/20 (10 tirages)
#   - blocage spatial 1500x1000km, 200x200km, 20x20km (10 folds cible)
#   - blocage temporel 1j (10 folds cible)
# on charge les meilleurs hyperparamètres trouvés en 10_run_tuning.R et
# on entraîne 10 modèles (un par fold), en sauvegardant TOUS les
# diagnostics et plots demandés.
#
# Sorties, sous outputs_pipeline/training/<freq>kHz/<model>/<schema>/ :
#   - models.rds              (liste des N_CV modèles, un par fold)
#   - metrics_par_fold.csv
#   - obs_pred_all.csv
#   - importance_all.csv
#   - plots (voir generate_and_save_all_plots ci-dessous)

source("R/00_config.R")
source("R/01_data_prep.R")
source("R/02_folds.R")
source("R/03_models.R")
source("R/04_diagnostics.R")
source("R/05_plots.R")

training_dir <- path_out("training")
tuning_dir   <- path_out("tuning")
dir.create(training_dir, showWarnings = FALSE, recursive = TRUE)

generate_and_save_all_plots <- function(cv_res, lc, scheme, out_dir, label) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  saveRDS(cv_res$models, file.path(out_dir, "models.rds"))
  write.csv(cv_res$metrics,  file.path(out_dir, "metrics_par_fold.csv"), row.names = FALSE)
  write.csv(cv_res$obs_pred, file.path(out_dir, "obs_pred_all.csv"), row.names = FALSE)
  write.csv(cv_res$importance, file.path(out_dir, "importance_all.csv"), row.names = FALSE)

  save_plot <- function(p, name, w = 8, h = 6) {
    ggsave(file.path(out_dir, name), p, width = w, height = h, dpi = 150)
  }

  # -- performance / entraînement --
  save_plot(plot_learning_curve(lc$summary, subtitle = label), "01_learning_curve.png")
  save_plot(plot_obs_vs_pred_scatter(cv_res$obs_pred, subtitle = label), "02_obs_vs_pred_scatter.png")
  save_plot(plot_obs_vs_pred_hist(cv_res$obs_pred, subtitle = label), "03_obs_vs_pred_hist.png")
  save_plot(plot_obs_vs_pred_by_fold(cv_res$obs_pred, subtitle = label), "04_obs_vs_pred_by_fold.png", 10, 8)
  save_plot(plot_metric_bar(cv_res$metrics, "rmse_test", subtitle = label), "05_rmse_par_fold.png")
  save_plot(plot_metric_bar(cv_res$metrics, "r2_test", fill = "darkgreen", hline0 = TRUE, subtitle = label),
            "06_r2_par_fold.png")
  save_plot(plot_importance_mean_sd(cv_res$importance, subtitle = label), "07_importance_variables.png")

  # -- cartes de résidus --
  save_plot(plot_residual_map(cv_res$obs_pred, subtitle = label), "08_carte_residus.png")
  save_plot(plot_abs_residual_map(cv_res$obs_pred, subtitle = label), "09_carte_residus_abs.png")

  # -- métadonnées : variance/distances/distributions --
  save_plot(plot_metric_bar(cv_res$metrics, "var_intra_fold_test", fill = "orange", subtitle = label),
            "10_variance_intra_fold.png")
  save_plot(plot_distances_per_fold(cv_res$metrics, subtitle = label), "11_distances_test_train.png", 9, 5)
  save_plot(plot_covariate_stats(cv_res$covariate_stats, subtitle = label), "12_covariable_stats.png", 10, 8)
  save_plot(plot_covariate_variance(cv_res$covariate_stats, subtitle = label), "13_covariable_variance.png", 10, 8)
  save_plot(plot_fod_distribution(cv_res$fod_dist, subtitle = label), "14_fod_distribution.png", 9, 7)
  save_plot(plot_distributions_per_fold(cv_res$numeric_dist, subtitle = label),
            "15_distributions_par_fold.png", 12, max(6, 1.2 * length(scheme$folds)))

  # -- localisation spatiale train/test (surtout utile pour le blocage) --
  save_plot(plot_points_colored_by_fold(scheme, subtitle = label), "16_carte_points_par_fold.png")
  if (scheme$scheme != "naive_RS_80_20") {
    save_plot(plot_spatial_train_test(scheme, subtitle = label), "17_carte_train_test_buffer.png", 10, 8)
  }

  global_rmse <- rmse_fn(cv_res$obs_pred$obs, cv_res$obs_pred$pred)
  global_r2   <- 1 - sum((cv_res$obs_pred$obs - cv_res$obs_pred$pred)^2) /
    sum((cv_res$obs_pred$obs - mean(cv_res$obs_pred$obs))^2)
  writeLines(
    sprintf("RMSE globale (out-of-fold) = %.4f\nR2 global (out-of-fold) = %.4f\nNombre de folds = %d",
            global_rmse, global_r2, length(scheme$folds)),
    file.path(out_dir, "resume_global.txt")
  )

  cat(sprintf("  -> %s : RMSE=%.4f R2=%.4f (%d folds) enregistre dans %s\n",
              label, global_rmse, global_r2, length(scheme$folds), out_dir))
}

MODEL_SCHEMES <- c("naive_RS_80_20",
                   paste0("blocked_spatial_", map_chr(SPATIAL_RESOLUTIONS, "label")),
                   paste0("blocked_temporal_", map_chr(TEMPORAL_RESOLUTIONS, "label")))

for (freq in FREQS) {

  cat("\n============================================================\n")
  cat("ENTRAINEMENT -- FREQUENCE :", freq, "kHz\n")
  cat("============================================================\n")

  prep_complete <- load_and_clean(freq, drop_na_numeric = TRUE)   # pour RF
  prep_xgb      <- load_and_clean(freq, drop_na_numeric = FALSE)  # pour XGB
  fod_levels    <- prep_xgb$fod_levels

  schemes_rf  <- build_all_schemes(prep_complete$df)
  schemes_xgb <- build_all_schemes(prep_xgb$df)

  for (scheme_name in MODEL_SCHEMES) {

    # ---- RF ----
    cat("\n-- RF --", scheme_name, "\n")
    tune_rf <- readRDS(file.path(tuning_dir, sprintf("rf_%dkHz_%s.rds", freq, scheme_name)))
    backend_rf <- make_backend("rf")
    scheme_rf  <- schemes_rf[[scheme_name]]

    cv_rf <- run_cv_scheme(scheme_rf, tune_rf$best_params, backend_rf, label = scheme_name)
    lc_rf <- compute_learning_curve(scheme_rf, tune_rf$best_params, backend_rf)

    out_dir_rf <- file.path(training_dir, paste0(freq, "kHz"), "rf", scheme_name)
    generate_and_save_all_plots(cv_rf, lc_rf, scheme_rf, out_dir_rf,
                                 label = sprintf("RF - %d kHz - %s", freq, scheme_name))

    # ---- XGB ----
    cat("\n-- XGB --", scheme_name, "\n")
    tune_xgb <- readRDS(file.path(tuning_dir, sprintf("xgb_%dkHz_%s.rds", freq, scheme_name)))
    backend_xgb <- make_backend("xgb", fod_levels = fod_levels)
    scheme_xgb  <- schemes_xgb[[scheme_name]]

    cv_xgb <- run_cv_scheme(scheme_xgb, tune_xgb$best_params, backend_xgb, label = scheme_name)
    lc_xgb <- compute_learning_curve(scheme_xgb, tune_xgb$best_params, backend_xgb)

    out_dir_xgb <- file.path(training_dir, paste0(freq, "kHz"), "xgb", scheme_name)
    generate_and_save_all_plots(cv_xgb, lc_xgb, scheme_xgb, out_dir_xgb,
                                 label = sprintf("XGB - %d kHz - %s", freq, scheme_name))
  }
}

cat("\nEntrainement termine. Resultats dans :", normalizePath(training_dir), "\n")
