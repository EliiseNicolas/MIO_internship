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
# ECHELLE DE COULEUR PARTAGEE (cf. SHARED_SCALE_SCOPE, 00_config.R) :
# toutes les cartes (RF, XGB avec/sans NA, toutes fréquences/schémas)
# utilisent la MEME échelle viridis, calculée comme l'union de la plage
# de NASC observée à l'entraînement et de la plage effectivement prédite
# -- pour qu'une même couleur représente la même valeur d'une carte à
# l'autre. Structure en deux phases : Phase 1 calcule toutes les
# prédictions (léger : seulement lon/lat/NASC_pred sont gardés en
# mémoire, pas les covariables complètes) ; Phase 2 calcule l'échelle
# partagée puis génère tous les PNG.
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

# ---------------------------------------------------------------------
# PHASE 1 : calcule toutes les predictions, stocke seulement
# lon/lat/NASC_pred (leger) + une liste de "plot_jobs" a generer en
# phase 2 une fois l'echelle commune connue.
# ---------------------------------------------------------------------
all_predictions <- list()
plot_jobs <- list()
training_nasc_range <- list()   # plage NASC observee a l'entrainement, par freq

day_ds <- load_grid_all_dates()
date_idx_single <- get_date_index(day_ds, TARGET_DATE_SINGLE)

for (freq in FREQS) {

  cat("\n============================================================\n")
  cat("PREDICTION GRILLE (PHASE 1) -- FREQUENCE :", freq, "kHz\n")
  cat("============================================================\n")

  prep_xgb <- load_and_clean(freq, drop_na_numeric = FALSE)
  fod_levels <- prep_xgb$fod_levels
  training_nasc_range[[as.character(freq)]] <- range(prep_xgb$df$NASC, na.rm = TRUE)

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

      layer_rf <- paste(freq, scheme_name, "RF", sep = "_")
      all_predictions[[layer_rf]] <- grid_rf %>%
        transmute(lon, lat, NASC_pred, freq = freq, scheme = scheme_name, model = "RF (sans NA)", layer_id = layer_rf)

      plot_jobs[[length(plot_jobs) + 1]] <- list(
        type = "single", layer_id = layer_rf, out_dir = out_dir, filename = "rf_prediction_map_no_NA.png",
        title = "NASC predit - RF (moyenne 10 folds)",
        subtitle = paste0(label_base, " -- ", nrow(grid_rf), "/", nrow(grid_all), " pixels (complets uniquement)")
      )
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

      layer_all   <- paste(freq, scheme_name, "XGB_avecNA", sep = "_")
      layer_clean <- paste(freq, scheme_name, "XGB_sansNA", sep = "_")
      all_predictions[[layer_all]] <- grid_xgb_all %>%
        transmute(lon, lat, NASC_pred, freq = freq, scheme = scheme_name, model = "XGB (avec NA)", layer_id = layer_all)
      all_predictions[[layer_clean]] <- grid_xgb_clean %>%
        transmute(lon, lat, NASC_pred, freq = freq, scheme = scheme_name, model = "XGB (sans NA)", layer_id = layer_clean)

      plot_jobs[[length(plot_jobs) + 1]] <- list(
        type = "xgb_combo", layer_id_all = layer_all, layer_id_clean = layer_clean,
        out_dir = out_dir, filename = "xgb_prediction_map_with_and_without_NA.png",
        n_all = nrow(grid_xgb_all), n_clean = nrow(grid_xgb_clean),
        title = paste("NASC predit - XGB (moyenne 10 folds) -", label_base)
      )
    } else {
      cat("  [!] modele XGB introuvable pour", scheme_name, "-- avez-vous lance 11_run_training.R ?\n")
    }
  }
}

cat(sprintf("\nPhase 1 terminee : %d couches de prediction calculees.\n", length(all_predictions)))

predictions_all <- bind_rows(all_predictions)
saveRDS(predictions_all, file.path(prediction_dir, "predictions_all_combined.rds"))

# ---------------------------------------------------------------------
# PHASE 2 : echelle de couleur partagee, puis generation de tous les PNG
# ---------------------------------------------------------------------
cat("\n============================================================\n")
cat("GENERATION DES CARTES (PHASE 2) -- echelle de couleur partagee\n")
cat("============================================================\n")

scope_key_freq <- function(f) if (SHARED_SCALE_SCOPE == "per_freq") as.character(f) else "ALL"

compute_shared_limits <- function(freqs_in_scope) {
  pred_vals <- predictions_all %>% filter(freq %in% freqs_in_scope) %>% pull(NASC_pred)
  train_vals <- unlist(training_nasc_range[as.character(freqs_in_scope)])
  range(c(pred_vals, train_vals), na.rm = TRUE)
}

if (SHARED_SCALE_SCOPE == "per_freq") {
  limits_by_scope <- setNames(lapply(FREQS, function(f) compute_shared_limits(f)), as.character(FREQS))
} else {
  limits_by_scope <- list(ALL = compute_shared_limits(FREQS))
}

get_limits_for_freq <- function(f) limits_by_scope[[scope_key_freq(f)]]

for (job in plot_jobs) {

  if (job$type == "single") {
    grid_df <- all_predictions[[job$layer_id]]
    lims <- get_limits_for_freq(grid_df$freq[1])
    p <- plot_prediction_map(grid_df, title = job$title, subtitle = job$subtitle, limits = lims)
    ggsave(file.path(job$out_dir, job$filename), p, width = 8, height = 6, dpi = 150)
    cat("  ->", file.path(job$out_dir, job$filename), "\n")

  } else if (job$type == "xgb_combo") {
    grid_all_j   <- all_predictions[[job$layer_id_all]]
    grid_clean_j <- all_predictions[[job$layer_id_clean]]
    lims <- get_limits_for_freq(grid_all_j$freq[1])

    p_all <- plot_prediction_map(grid_all_j, title = "Avec NA (XGBoost gere le manquant)",
                                  subtitle = paste0(job$n_all, " pixels predits"), limits = lims)
    p_clean <- plot_prediction_map(grid_clean_j, title = "Sans NA (pixels complets uniquement)",
                                    subtitle = paste0(job$n_clean, " pixels predits"), limits = lims)
    p_combo <- (p_all + p_clean) +
      plot_annotation(title = job$title,
                       subtitle = "Comparaison : gestion native des NA (XGBoost) vs filtrage strict")
    ggsave(file.path(job$out_dir, job$filename), p_combo, width = 12, height = 6, dpi = 150)
    cat("  ->", file.path(job$out_dir, job$filename), "\n")
  }
}

# ---- Carte comparative globale : tous modeles / schemas / frequences ----
if (nrow(predictions_all) > 0) {
  for (f in FREQS) {
    sub <- predictions_all %>% filter(freq == f)
    if (nrow(sub) == 0) next
    lims <- get_limits_for_freq(f)
    p_compare <- ggplot(sub, aes(x = lon, y = lat, fill = NASC_pred)) +
      geom_raster() +
      scale_fill_viridis_c(limits = lims) +
      coord_quickmap() +
      facet_grid(model ~ scheme) +
      theme_bw() +
      labs(title = paste0("NASC predit -- comparaison tous schemas/modeles -- ", f, " kHz"),
           x = "Longitude", y = "Latitude", fill = "log10(NASC)")
    ggsave(file.path(prediction_dir, paste0("comparaison_globale_", f, "kHz.png")),
           p_compare, width = 16, height = 10, dpi = 150)
  }
  cat("\nGrilles de prediction (numeriques) sauvegardees dans :",
      file.path(prediction_dir, "predictions_all_combined.rds"), "\n")
}

cat("\nPredictions terminees. Resultats dans :", normalizePath(prediction_dir), "\n")
