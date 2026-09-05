# =====================================================================
# 13_run_grid_prediction_multidate.R -- SCRIPT 3bis : prédiction 133 dates
# =====================================================================
# Utilise le même fichier grille unique que 12_run_grid_prediction.R
# (`PATH_GRID_ALL_DATES`, structure corrigée), en bouclant sur TOUTES
# les dates disponibles au lieu d'une seule. Les ratios *_totpig/Chla_total
# sont déjà calculés dans le fichier -- pas de formule à deviner.
#
# Pour chaque fréquence x schéma (par défaut seulement
# `MULTIDATE_PREDICTION_SCHEMES`, cf. 00_config.R -- étends cette liste
# si tu veux aussi les 133 cartes pour les schémas bloqués) x date :
#   - RF  : carte (pixels complets uniquement)
#   - XGB : carte avec NA gérés nativement
#
# ATTENTION AU VOLUME : 133 dates x 2 modèles x N schémas x N fréquences
# peut représenter plusieurs centaines de fichiers PNG. Réduis
# MULTIDATE_PREDICTION_SCHEMES/FREQS si tu veux limiter la sortie.
#
# ECHELLE DE COULEUR PARTAGEE (cf. SHARED_SCALE_SCOPE, 00_config.R) :
# chaque carte (RF, XGB, toutes dates/schémas/fréquences selon la
# portée choisie) utilise la plage de NASC OBSERVEE A L'ENTRAINEMENT
# comme échelle fixe -- calculée une seule fois par fréquence (pas
# besoin de calculer d'abord les 133x2 prédictions pour connaître leur
# étendue, ce qui éviterait un passage supplémentaire coûteux sur toute
# la grille). Une prédiction ponctuelle qui dépasserait légèrement cette
# plage sera simplement écrêtée visuellement (couleur extrême de
# l'échelle) -- acceptable pour une série de 133 cartes.
#
# Sorties, sous outputs_pipeline/predictions_multidate/<freq>kHz/<schema>/ :
#   - rf_YYYYMMDD.png, xgb_YYYYMMDD.png (une carte par date)
#   - timeseries_mean_prediction.csv (moyenne spatiale prédite par date)
#   - timeseries_mean_prediction.png (évolution temporelle, RF vs XGB)

source("R/00_config.R")
source("R/01_data_prep.R")
source("R/02_folds.R")
source("R/03_models.R")
source("R/05_plots.R")

training_dir <- path_out("training")
out_root     <- path_out("predictions_multidate")
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

# Chargée UNE SEULE fois (fichier potentiellement volumineux, identique
# pour toutes les fréquences)
day_ds  <- load_grid_all_dates()
n_dates <- length(day_ds$date)

# Echelle de couleur partagee : plage de NASC observee a l'entrainement,
# calculee une fois par frequence (peu couteux, pas besoin de charger la
# grille pour ca).
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
  cat("PREDICTION MULTI-DATE -- FREQUENCE :", freq, "kHz\n")
  cat("============================================================\n")

  prep_xgb   <- load_and_clean(freq, drop_na_numeric = FALSE)  # juste pour fod_levels
  fod_levels <- prep_xgb$fod_levels
  lims       <- shared_limits[[as.character(freq)]]

  cat(sprintf("Grille multi-date : %d dates, %d x %d pixels\n",
              n_dates, length(day_ds$lon), length(day_ds$lat)))

  for (scheme_name in MULTIDATE_PREDICTION_SCHEMES) {

    out_dir <- file.path(out_root, paste0(freq, "kHz"), scheme_name)
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

    rf_models_path  <- file.path(training_dir, paste0(freq, "kHz"), "rf",  scheme_name, "models.rds")
    xgb_models_path <- file.path(training_dir, paste0(freq, "kHz"), "xgb", scheme_name, "models.rds")

    rf_models  <- if (file.exists(rf_models_path))  readRDS(rf_models_path)  else NULL
    xgb_models <- if (file.exists(xgb_models_path)) readRDS(xgb_models_path) else NULL

    if (is.null(rf_models) && is.null(xgb_models)) {
      cat("  [!] aucun modele trouve pour", scheme_name, "-- avez-vous lance 11_run_training.R ?\n")
      next
    }

    timeseries <- vector("list", n_dates)

    for (d in seq_len(n_dates)) {

      extracted  <- extract_grid_for_date(day_ds, d, fod_levels)
      grid_all   <- extracted$grid
      grid_clean <- grid_all[stats::complete.cases(grid_all[, COVARIATES_ALL]), ]
      date_str   <- format(extracted$date, "%Y%m%d")

      row <- tibble(date = extracted$date)

      if (!is.null(rf_models) && nrow(grid_clean) > 0) {
        preds_rf <- sapply(rf_models, function(m) predict_rf(m, grid_clean))
        grid_rf  <- grid_clean
        grid_rf$NASC_pred <- rowMeans(preds_rf)

        p_rf <- plot_prediction_map(
          grid_rf, title = paste0("NASC predit - RF - ", format(extracted$date, "%Y-%m-%d")),
          subtitle = sprintf("%d kHz - %s - %d/%d pixels complets", freq, scheme_name,
                              nrow(grid_rf), nrow(grid_all)),
          limits = lims
        )
        ggsave(file.path(out_dir, sprintf("rf_%s.png", date_str)), p_rf, width = 7, height = 5, dpi = 100)

        row$rf_mean_pred <- mean(grid_rf$NASC_pred, na.rm = TRUE)
        row$rf_pct_pixels_predits <- 100 * nrow(grid_rf) / nrow(grid_all)
      }

      if (!is.null(xgb_models)) {
        preds_xgb <- sapply(xgb_models, function(m) predict_xgb(m, grid_all, fod_levels = fod_levels))
        grid_xgb  <- grid_all
        grid_xgb$NASC_pred <- rowMeans(preds_xgb)

        p_xgb <- plot_prediction_map(
          grid_xgb, title = paste0("NASC predit - XGB (avec NA) - ", format(extracted$date, "%Y-%m-%d")),
          subtitle = sprintf("%d kHz - %s", freq, scheme_name),
          limits = lims
        )
        ggsave(file.path(out_dir, sprintf("xgb_%s.png", date_str)), p_xgb, width = 7, height = 5, dpi = 100)

        row$xgb_mean_pred <- mean(grid_xgb$NASC_pred, na.rm = TRUE)
      }

      timeseries[[d]] <- row
      if (d %% 10 == 0) cat(sprintf("  ... %d / %d dates traitees\n", d, n_dates))
    }

    ts_df <- bind_rows(timeseries)
    write.csv(ts_df, file.path(out_dir, "timeseries_mean_prediction.csv"), row.names = FALSE)

    ts_long <- ts_df %>%
      select(date, any_of(c("rf_mean_pred", "xgb_mean_pred"))) %>%
      pivot_longer(-date, names_to = "modele", values_to = "mean_pred") %>%
      filter(!is.na(mean_pred))

    p_ts <- ggplot(ts_long, aes(x = date, y = mean_pred, color = modele)) +
      geom_line() + geom_point(size = 0.8) +
      labs(title = "Evolution temporelle du NASC predit (moyenne spatiale)",
           subtitle = sprintf("%d kHz - %s", freq, scheme_name),
           x = "Date", y = "NASC predit moyen (log10)", color = NULL) +
      theme_bw()
    ggsave(file.path(out_dir, "timeseries_mean_prediction.png"), p_ts, width = 9, height = 5, dpi = 150)

    cat(sprintf("  -> %s : %d cartes/date enregistrees dans %s\n", scheme_name, n_dates, out_dir))
  }
}

cat("\nPrediction multi-date terminee. Resultats dans :", normalizePath(out_root), "\n")
