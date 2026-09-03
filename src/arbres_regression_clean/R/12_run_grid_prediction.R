# =====================================================================
# 12_run_grid_prediction.R -- SCRIPT 3 : prédiction sur grille (mono-date)
# =====================================================================
# Utilise le fichier grille unique (`PATH_GRID_ALL_DATES`, structure
# corrigée -- voir data_generation/generate_ds_ftle_pig_fod_all_dates.R),
# dont on extrait UNE SEULE date (`TARGET_DATE_SINGLE`, 2023-01-26 par
# défaut) avec `extract_grid_for_date()` : les ratios *_totpig/Chla_total
# sont déjà calculés dans le fichier, aucune formule à deviner. Pour les
# 133 dates d'un coup, voir `13_run_grid_prediction_multidate.R`.
# Pour chaque (fréquence x modèle x schéma), on recharge les N_CV
# modèles entraînés en 11_run_training.R, on prédit sur la grille
# spatiale du jour choisi avec CHACUN des 10 modèles, et on moyenne.
#
#   - RF   : uniquement les pixels complets (pas de NA) -> une seule
#            carte, avec des "trous" (NA) là où au moins un pigment
#            manque.
#   - XGB  : deux cartes : avec NA (XGBoost gère nativement le
#            manquant) et sans NA (pixels complets uniquement), pour
#            comparaison directe avec RF.
#
# Les diagnostics "réel vs prédit" (nuage de points, histogramme, par
# fold) et "importance des variables" sont déjà produits en
# 11_run_training.R à partir des prédictions out-of-fold ; ils ne sont
# PAS recalculés ici car la grille n'a pas de NASC observé (c'est une
# vraie zone de prédiction, pas un jeu de test).

source("R/00_config.R")
source("R/01_data_prep.R")
source("R/02_folds.R")
source("R/03_models.R")
source("R/05_plots.R")

training_dir  <- path_out("training")
prediction_dir <- path_out("predictions")
dir.create(prediction_dir, showWarnings = FALSE, recursive = TRUE)

MODEL_SCHEMES <- c("naive_RS_80_20",
                   paste0("blocked_spatial_", map_chr(SPATIAL_RESOLUTIONS, "label")),
                   paste0("blocked_temporal_", map_chr(TEMPORAL_RESOLUTIONS, "label")))

predict_mean_rf <- function(models, grid_clean) {
  preds <- sapply(models, function(m) predict_rf(m, grid_clean))
  rowMeans(preds)
}
predict_mean_xgb <- function(models, grid_df, fod_levels) {
  preds <- sapply(models, function(m) predict_xgb(m, grid_df, fod_levels = fod_levels))
  rowMeans(preds)
}

all_predictions <- list()
day_ds <- load_grid_all_dates()
date_idx_single <- get_date_index(day_ds, TARGET_DATE_SINGLE)

for (freq in FREQS) {

  cat("\n============================================================\n")
  cat("PREDICTION GRILLE -- FREQUENCE :", freq, "kHz\n")
  cat("============================================================\n")

  prep_xgb <- load_and_clean(freq, drop_na_numeric = FALSE)  # juste pour fod_levels
  fod_levels <- prep_xgb$fod_levels

  grid_info  <- extract_grid_for_date(day_ds, date_idx_single, fod_levels)
  grid_all   <- grid_info$grid
  grid_clean <- grid_all[stats::complete.cases(grid_all[, COVARIATES_ALL]), ]
  date_label <- format(grid_info$date, "%Y-%m-%d")

  cat(sprintf("Grille : %d pixels total, %d complets (%.1f%%)\n",
              nrow(grid_all), nrow(grid_clean), 100 * nrow(grid_clean) / nrow(grid_all)))

  for (scheme_name in MODEL_SCHEMES) {

    out_dir <- file.path(prediction_dir, paste0(freq, "kHz"), scheme_name)
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    label_base <- sprintf("%d kHz - %s - %s", freq, scheme_name, date_label)

    # ---- RF : pixels complets uniquement ----
    rf_models_path <- file.path(training_dir, paste0(freq, "kHz"), "rf", scheme_name, "models.rds")
    if (file.exists(rf_models_path)) {
      rf_models <- readRDS(rf_models_path)
      grid_rf <- grid_clean
      grid_rf$NASC_pred <- predict_mean_rf(rf_models, grid_rf)

      p_rf <- plot_prediction_map(
        grid_rf, title = "NASC predit - RF (moyenne 10 folds)",
        subtitle = paste0(label_base, " -- ", nrow(grid_rf), "/", nrow(grid_all), " pixels (complets uniquement)")
      )
      ggsave(file.path(out_dir, "rf_prediction_map_no_NA.png"), p_rf, width = 8, height = 6, dpi = 150)

      all_predictions[[paste(freq, scheme_name, "rf", sep = "_")]] <-
        grid_rf %>% transmute(lon, lat, NASC_pred, freq = freq, scheme = scheme_name, model = "RF (sans NA)")

      cat("  RF  :", file.path(out_dir, "rf_prediction_map_no_NA.png"), "\n")
    } else {
      cat("  [!] modele RF introuvable pour", scheme_name, "-- avez-vous lance 11_run_training.R ?\n")
    }

    # ---- XGB : avec ET sans NA ----
    xgb_models_path <- file.path(training_dir, paste0(freq, "kHz"), "xgb", scheme_name, "models.rds")
    if (file.exists(xgb_models_path)) {
      xgb_models <- readRDS(xgb_models_path)

      grid_xgb_all   <- grid_all
      grid_xgb_clean <- grid_clean
      grid_xgb_all$NASC_pred   <- predict_mean_xgb(xgb_models, grid_xgb_all,   fod_levels)
      grid_xgb_clean$NASC_pred <- predict_mean_xgb(xgb_models, grid_xgb_clean, fod_levels)

      shared_limits <- range(c(grid_xgb_all$NASC_pred, grid_xgb_clean$NASC_pred), na.rm = TRUE)

      p_all <- plot_prediction_map(grid_xgb_all, title = "Avec NA (XGBoost gere le manquant)",
                                    subtitle = paste0(nrow(grid_xgb_all), " pixels predits"),
                                    limits = shared_limits)
      p_clean <- plot_prediction_map(grid_xgb_clean, title = "Sans NA (pixels complets uniquement)",
                                      subtitle = paste0(nrow(grid_xgb_clean), " pixels predits"),
                                      limits = shared_limits)

      p_combo <- (p_all + p_clean) +
        plot_annotation(title = paste("NASC predit - XGB (moyenne 10 folds) -", label_base),
                         subtitle = "Comparaison : gestion native des NA (XGBoost) vs filtrage strict")
      ggsave(file.path(out_dir, "xgb_prediction_map_with_and_without_NA.png"),
             p_combo, width = 12, height = 6, dpi = 150)

      all_predictions[[paste(freq, scheme_name, "xgb_all", sep = "_")]] <-
        grid_xgb_all   %>% transmute(lon, lat, NASC_pred, freq = freq, scheme = scheme_name, model = "XGB (avec NA)")
      all_predictions[[paste(freq, scheme_name, "xgb_clean", sep = "_")]] <-
        grid_xgb_clean %>% transmute(lon, lat, NASC_pred, freq = freq, scheme = scheme_name, model = "XGB (sans NA)")

      cat("  XGB :", file.path(out_dir, "xgb_prediction_map_with_and_without_NA.png"), "\n")
    } else {
      cat("  [!] modele XGB introuvable pour", scheme_name, "-- avez-vous lance 11_run_training.R ?\n")
    }
  }
}

# ---- Carte comparative globale : tous modeles / schemas / frequences ----
predictions_all <- bind_rows(all_predictions)
if (nrow(predictions_all) > 0) {
  for (f in FREQS) {
    sub <- predictions_all %>% filter(freq == f)
    if (nrow(sub) == 0) next
    p_compare <- ggplot(sub, aes(x = lon, y = lat, fill = NASC_pred)) +
      geom_raster() +
      scale_fill_viridis_c() +
      coord_quickmap() +
      facet_grid(model ~ scheme) +
      theme_bw() +
      labs(title = paste0("NASC predit -- comparaison tous schemas/modeles -- ", f, " kHz"),
           x = "Longitude", y = "Latitude", fill = "log10(NASC)")
    ggsave(file.path(prediction_dir, paste0("comparaison_globale_", f, "kHz.png")),
           p_compare, width = 16, height = 10, dpi = 150)
  }
}

cat("\nPredictions terminees. Resultats dans :", normalizePath(prediction_dir), "\n")
