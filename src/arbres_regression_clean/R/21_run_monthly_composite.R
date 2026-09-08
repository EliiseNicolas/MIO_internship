# =====================================================================
# 21_run_monthly_composite.R -- SCRIPT 11 : composite mensuel
# =====================================================================
# Pour chaque freq x schéma (MONTHLY_COMPOSITE_SCHEMES) x mois calendaire
# présent dans la grille multi-date : prédit avec XGB "avec NA" (le
# modèle qui produit une prédiction sur TOUS les pixels, sans les trous
# de RF complete-case) pour CHAQUE jour du mois, puis agrège pixel par
# pixel sur l'ensemble des jours :
#
#   - carte MOYENNE mensuelle : moyenne des prédictions journalières
#     valides (NA ignorés, pas propagés -- un pixel n'est NA sur la
#     carte mensuelle que si TOUS les jours du mois étaient NA pour ce
#     pixel).
#   - carte de PURETE : fraction de jours du mois où ce pixel avait une
#     prédiction valide (1 = tous les jours du mois ont une valeur non-NA
#     pour ce pixel, 0 = aucun jour). Équivalent du concept de "purity"
#     en compositing satellite -- indique la CONFIANCE à accorder à
#     chaque pixel de la carte moyenne, pas juste sa valeur.
#
# ATTENTION AU VOLUME : un mois de 30 jours = 30 prédictions sur toute
# la grille par freq x schéma -- coûteux si beaucoup de schémas/mois
# sont sélectionnés. Réduis MONTHLY_COMPOSITE_SCHEMES/MONTHLY_COMPOSITE_MONTHS
# si besoin.
#
# Sorties, sous outputs_pipeline/monthly_composite/<freq>kHz/<schema>/ :
#   - monthly_mean_<YYYY-MM>.png, monthly_purity_<YYYY-MM>.png
#   - monthly_composite_<YYYY-MM>.rds (grille numerique : lon, lat, mean, purity, n_valid, n_days)

source("R/00_config.R")
source("R/01_data_prep.R")
source("R/02_folds.R")
source("R/03_models.R")
source("R/05_plots.R")
source("R/10_basemap.R")

training_dir <- path_out("training")
out_root     <- path_out("monthly_composite")
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

basemap <- load_basemap_sf()

# Schémas à traiter (par défaut seulement naive -- même logique que
# MULTIDATE_PREDICTION_SCHEMES, étends si besoin).
MONTHLY_COMPOSITE_SCHEMES <- c("naive_RS_80_20")

# Mois à traiter : NULL = tous les mois présents dans la grille
# multi-date (déduits automatiquement) ; sinon un vecteur de chaînes
# "YYYY-MM", ex. c("2018-01", "2018-02").
MONTHLY_COMPOSITE_MONTHS <- NULL

predict_mean_xgb <- function(models, grid_df, fod_levels) {
  preds <- sapply(models, function(m) predict_xgb(m, grid_df, fod_levels = fod_levels))
  rowMeans(preds)
}

day_ds <- load_grid_all_dates()
all_months <- format(day_ds$date, "%Y-%m")
months_to_process <- if (is.null(MONTHLY_COMPOSITE_MONTHS)) sort(unique(all_months)) else MONTHLY_COMPOSITE_MONTHS

cat(sprintf("Mois disponibles dans la grille : %d (%s a %s)\n",
            length(unique(all_months)), min(all_months), max(all_months)))
cat(sprintf("Mois traites : %s\n", paste(months_to_process, collapse = ", ")))

for (freq in FREQS) {

  cat("\n============================================================\n")
  cat("COMPOSITE MENSUEL -- FREQUENCE :", freq, "kHz\n")
  cat("============================================================\n")

  prep_xgb   <- load_and_clean(freq, drop_na_numeric = FALSE)  # pour fod_levels + echelle couleur
  fod_levels <- prep_xgb$fod_levels
  lims       <- get_prediction_color_limits(freq, range(prep_xgb$df$NASC, na.rm = TRUE))

  for (scheme_name in MONTHLY_COMPOSITE_SCHEMES) {

    xgb_models_path <- file.path(training_dir, paste0(freq, "kHz"), "xgb", scheme_name, "models.rds")
    if (!file.exists(xgb_models_path)) {
      cat("  [!] modele XGB introuvable pour", scheme_name, "-- avez-vous lance 11_run_training.R ?\n")
      next
    }
    xgb_models <- readRDS(xgb_models_path)

    out_dir <- file.path(out_root, paste0(freq, "kHz"), scheme_name)
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

    for (ym in months_to_process) {

      date_idx_month <- which(all_months == ym)
      if (length(date_idx_month) == 0) {
        cat("  [!] aucune date pour", ym, "-- saute.\n")
        next
      }
      n_days <- length(date_idx_month)
      cat(sprintf("  %s - %s : %d jours\n", scheme_name, ym, n_days))

      # Matrice [n_pixels x n_days] des predictions journalieres
      n_pixels <- length(day_ds$lon) * length(day_ds$lat)
      daily_preds <- matrix(NA_real_, nrow = n_pixels, ncol = n_days)
      lon_vec <- lat_vec <- NULL

      for (j in seq_along(date_idx_month)) {
        extracted <- extract_grid_for_date(day_ds, date_idx_month[j], fod_levels)
        grid_all  <- extracted$grid
        if (is.null(lon_vec)) { lon_vec <- grid_all$lon; lat_vec <- grid_all$lat }
        daily_preds[, j] <- predict_mean_xgb(xgb_models, grid_all, fod_levels)
      }

      n_valid    <- rowSums(!is.na(daily_preds))
      mean_pred  <- rowMeans(daily_preds, na.rm = TRUE)
      mean_pred[n_valid == 0] <- NA_real_   # aucun jour valide -> NA (pas de composite possible)
      purity     <- n_valid / n_days

      composite_df <- tibble(
        lon = lon_vec, lat = lat_vec,
        mean_pred = mean_pred, purity = purity, n_valid = n_valid, n_days = n_days
      )
      saveRDS(composite_df, file.path(out_dir, sprintf("monthly_composite_%s.rds", ym)))

      # ---- carte moyenne mensuelle ----
      p_mean <- ggplot(composite_df, aes(x = lon, y = lat, fill = mean_pred)) +
        geom_raster()
      if (!is.null(basemap)) {
        p_mean <- p_mean +
          geom_sf(data = basemap, inherit.aes = FALSE, fill = "grey40", color = "grey20", linewidth = 0.2) +
          coord_sf(xlim = range(lon_vec), ylim = range(lat_vec), expand = FALSE)
      } else {
        p_mean <- p_mean + coord_quickmap()
      }
      p_mean <- p_mean +
        scale_fill_viridis_c(limits = lims) +
        theme_pipeline +
        labs(title = "NASC predit -- composite mensuel (moyenne, XGB avec NA)",
             subtitle = sprintf("%d kHz - %s - %s (%d jours)", freq, scheme_name, ym, n_days),
             x = "Longitude", y = "Latitude", fill = "log10(NASC)\nmoyen")
      ggsave(file.path(out_dir, sprintf("monthly_mean_%s.png", ym)), p_mean, width = 8, height = 6, dpi = 150)

      # ---- carte de purete ----
      p_purity <- ggplot(composite_df, aes(x = lon, y = lat, fill = purity)) +
        geom_raster()
      if (!is.null(basemap)) {
        p_purity <- p_purity +
          geom_sf(data = basemap, inherit.aes = FALSE, fill = "grey40", color = "grey20", linewidth = 0.2) +
          coord_sf(xlim = range(lon_vec), ylim = range(lat_vec), expand = FALSE)
      } else {
        p_purity <- p_purity + coord_quickmap()
      }
      p_purity <- p_purity +
        scale_fill_viridis_c(option = "magma", limits = c(0, 1)) +
        theme_pipeline +
        labs(title = "Purete du composite mensuel (fraction de jours valides par pixel)",
             subtitle = sprintf("%d kHz - %s - %s (%d jours)", freq, scheme_name, ym, n_days),
             x = "Longitude", y = "Latitude", fill = "Purete\n(0-1)")
      ggsave(file.path(out_dir, sprintf("monthly_purity_%s.png", ym)), p_purity, width = 8, height = 6, dpi = 150)

      cat(sprintf("    -> purete moyenne = %.2f (mediane = %.2f)\n",
                  mean(purity, na.rm = TRUE), stats::median(purity, na.rm = TRUE)))
    }
  }
}

cat("\nComposite mensuel termine. Resultats dans :", normalizePath(out_root), "\n")
