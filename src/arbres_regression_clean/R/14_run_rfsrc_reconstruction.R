# =====================================================================
# 14_run_rfsrc_reconstruction.R -- SCRIPT 4 : test de reconstruction RF
# =====================================================================
# Objectif : contrairement a RF (ranger), qui doit jeter tout pixel ou
# au moins une covariable manque (-> "trous" NA sur la carte),
# randomForestSRC gere le manquant NATIVEMENT via imputation
# (`na.action = "na.impute"`), a l'entrainement ET a la prediction. On
# teste ici s'il peut "reconstruire" une carte ENTIERE (2023-01-26),
# sans aucun trou, a partir des covariables NORMALISEES (COVARIATES_ALL
# -- ftle, Chla_total, *_totpig, fod), les memes que le pipeline
# principal CART/RF/XGB, extraites du meme fichier grille unique que
# les scripts 12/13.
#
# C'est un script de TEST isole (pas de tuning complet, pas de CV
# multi-fold) : le but est de comparer la "reconstruction" (nombre de
# pixels predits) entre RF classique (complete-case) et randomForestSRC
# (imputation), pas de comparer des schemas de validation croisee.
#
# Sorties, sous outputs_pipeline/rfsrc_reconstruction/<freq>kHz/ :
#   - model.rds
#   - metrics_test.txt (RMSE/R2 sur un split naive 80/20, pour info)
#   - importance_variables.png
#   - reconstruction_map.png (carte complete, aucun pixel manquant)
#   - reconstruction_stats.txt (% de pixels "sauves" par l'imputation,
#     cad qui auraient ete NA/jetes dans le pipeline RF ranger classique)

source("R/00_config.R")
source("R/01_data_prep.R")
source("R/02_folds.R")
source("R/03_models.R")
source("R/04_diagnostics.R")
source("R/05_plots.R")

out_root <- path_out("rfsrc_reconstruction")
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

# Hyperparametres par defaut (pas de tuning complet pour ce test
# exploratoire -- utilise 06_tuning.R::tune_model() avec
# make_backend("rfsrc") et default_tuning_grid("rfsrc") si tu veux
# un vrai tuning avant de conclure quoi que ce soit).
RFSRC_PARAMS <- list(mtry = 3, nodesize = 15, ntree = 500)

day_ds <- load_grid_all_dates()
date_idx <- get_date_index(day_ds, TARGET_DATE_SINGLE)

# Echelle de couleur partagee avec 12_/13_ (meme logique, cf.
# SHARED_SCALE_SCOPE dans 00_config.R) : plage de NASC observee a
# l'entrainement, pour que la carte de reconstruction soit directement
# comparable aux cartes RF/XGB produites par 12_run_grid_prediction.R.
training_nasc_range <- setNames(
  lapply(FREQS, function(f) range(load_and_clean(f, drop_na_numeric = FALSE)$df$NASC, na.rm = TRUE)),
  as.character(FREQS)
)
shared_limits <- if (SHARED_SCALE_SCOPE == "per_freq") {
  training_nasc_range
} else {
  setNames(rep(list(range(unlist(training_nasc_range))), length(FREQS)), as.character(FREQS))
}

for (freq in FREQS) {

  cat("\n============================================================\n")
  cat("RECONSTRUCTION RFSRC -- FREQUENCE :", freq, "kHz\n")
  cat("============================================================\n")

  out_dir <- file.path(out_root, paste0(freq, "kHz"))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  lims <- shared_limits[[as.character(freq)]]

  # drop_na_numeric = FALSE : randomForestSRC gere le manquant nativement,
  # on garde donc les lignes avec NA sur les covariables numeriques
  # (seuls fod/NASC restent filtres, cf. 01_data_prep.R::load_and_clean).
  prep <- load_and_clean(freq, drop_na_numeric = FALSE)
  fod_levels <- prep$fod_levels
  backend <- make_backend("rfsrc")   # covariables normalisees (COVARIATES_ALL) par defaut

  # ---- Petit split naive 80/20 juste pour un chiffre de performance ----
  # (le but du script est la reconstruction de carte, pas la comparaison
  # de schemas de CV -- pour ca, utiliser le pipeline principal)
  scheme <- build_naive_folds(prep$df, n_folds = 1, frac_train = 0.8)
  diag   <- compute_fold_diagnostics(scheme$folds[[1]], scheme$data, RFSRC_PARAMS, backend)

  writeLines(
    sprintf("RMSE test = %.4f\nR2 test = %.4f\n(split naive 80/20, ntree=%d, mtry=%d, nodesize=%d)",
            diag$metrics$rmse_test, diag$metrics$r2_test,
            RFSRC_PARAMS$ntree, RFSRC_PARAMS$mtry, RFSRC_PARAMS$nodesize),
    file.path(out_dir, "metrics_test.txt")
  )
  cat(sprintf("  RMSE test = %.4f | R2 test = %.4f\n", diag$metrics$rmse_test, diag$metrics$r2_test))

  p_imp <- plot_importance_mean_sd(
    diag$importance %>% mutate(fold_id = 1),
    subtitle = sprintf("randomForestSRC (covariables normalisees) - %d kHz", freq)
  )
  ggsave(file.path(out_dir, "importance_variables.png"), p_imp, width = 7, height = 5, dpi = 150)

  # ---- Modele final : entraine sur TOUTES les donnees disponibles ----
  cat("  Entrainement du modele final (toutes les donnees)...\n")
  model_final <- backend$fit(prep$df, RFSRC_PARAMS)
  saveRDS(model_final, file.path(out_dir, "model.rds"))

  # ---- Reconstruction de la carte complete (2023-01-26) ----
  extracted <- extract_grid_for_date(day_ds, date_idx, fod_levels)
  grid_all  <- extracted$grid

  n_total      <- nrow(grid_all)
  n_incomplete <- sum(!stats::complete.cases(grid_all[, COVARIATES_ALL]))

  cat(sprintf("  Grille %s : %d pixels total, %d incomplets (%.1f%%) -- perdus avec RF classique\n",
              format(extracted$date), n_total, n_incomplete, 100 * n_incomplete / n_total))

  # predict_rfsrc() impute les covariables manquantes de `newdata` avant
  # de predire -> une prediction pour TOUS les pixels, y compris ceux
  # avec une ou plusieurs covariables manquantes.
  grid_all$NASC_pred <- backend$predict(model_final, grid_all)

  n_predicted <- sum(!is.na(grid_all$NASC_pred))
  cat(sprintf("  -> %d / %d pixels reconstruits (vs %d avec RF ranger classique, +%.1f%%)\n",
              n_predicted, n_total, n_total - n_incomplete,
              100 * (n_predicted - (n_total - n_incomplete)) / max(1, n_total - n_incomplete)))

  p_map <- plot_prediction_map(
    grid_all,
    title = "NASC reconstruit - randomForestSRC (imputation native des NA)",
    subtitle = sprintf("%d kHz - %s - %d/%d pixels reconstruits (%.1f%%), dont %d etaient incomplets",
                        freq, format(extracted$date, "%Y-%m-%d"),
                        n_predicted, n_total, 100 * n_predicted / n_total, n_incomplete),
    limits = lims
  )
  ggsave(file.path(out_dir, "reconstruction_map.png"), p_map, width = 8, height = 6, dpi = 150)

  # Sauvegarde numerique -- necessaire pour la comparaison quantitative
  # des cartes dans 17_run_map_comparison.R
  grid_to_save <- grid_all %>%
    transmute(lon, lat, NASC_pred, freq = freq, scheme = "naive_RS_80_20", model = "rfsrc",
              layer_id = paste(freq, "rfsrc_reconstruction", sep = "_"))
  saveRDS(grid_to_save, file.path(out_dir, "reconstruction_grid.rds"))

  writeLines(
    sprintf(paste0(
      "Date : %s\n",
      "Covariables : normalisees (COVARIATES_ALL -- ftle, Chla_total, *_totpig, fod)\n",
      "Pixels total : %d\n",
      "Pixels incomplets (au moins 1 covariable manquante) : %d (%.1f%%)\n",
      "Pixels reconstruits par randomForestSRC (imputation) : %d (%.1f%%)\n",
      "Pixels qui auraient ete predits par RF ranger (complete-case uniquement) : %d (%.1f%%)\n"
    ),
    format(extracted$date), n_total, n_incomplete, 100 * n_incomplete / n_total,
    n_predicted, 100 * n_predicted / n_total,
    n_total - n_incomplete, 100 * (n_total - n_incomplete) / n_total),
    file.path(out_dir, "reconstruction_stats.txt")
  )

  cat("  -> carte de reconstruction enregistree dans", out_dir, "\n")
}

cat("\nTest de reconstruction termine. Resultats dans :", normalizePath(out_root), "\n")
