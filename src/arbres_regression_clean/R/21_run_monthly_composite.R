# =====================================================================
# 21_run_monthly_composite.R -- SCRIPT 11 : composite mensuel
# =====================================================================
# Pour chaque freq x schéma (MONTHLY_COMPOSITE_SCHEMES) x mois calendaire
# présent dans la grille multi-date : prédit pour CHAQUE jour du mois
# avec le modèle choisi (MONTHLY_COMPOSITE_MODEL, "rf" ou "xgb"), puis
# agrège pixel par pixel sur l'ensemble des jours :
#
#   - carte MOYENNE mensuelle SEULE (monthly_mean_<YYYY-MM>.png) :
#     moyenne des prédictions journalières valides (NA ignorés, pas
#     propagés -- un pixel n'est NA que si TOUS les jours du mois
#     étaient NA pour lui). La carte n'a PAS besoin d'être complète
#     pour être utile.
#   - carte MOYENNE + PURETE EN TRANSPARENCE, même carte
#     (monthly_mean_avec_purete_<YYYY-MM>.png) : canal alpha = purete --
#     un pixel avec peu de jours valides apparaît plus pâle/transparent.
#   - carte de PURETE SEULE (monthly_purity_<YYYY-MM>.png) : nombre de
#     jours valides / nombre de jours total du mois, par pixel (0 à 1) --
#     pour une lecture précise de la valeur exacte, moins ambiguë que
#     la transparence seule.
#
# RF ("rf") vs XGB ("xgb") -- IMPORTANT, ça change le sens de la purete :
#   - RF (ranger) ne gere PAS les NA du tout : chaque jour, seuls les
#     pixels COMPLETS ce jour-la recoivent une prediction ; les autres
#     restent NA pour ce jour precis (variable jour apres jour selon la
#     couverture satellite). La purete resultante varie donc REELLEMENT
#     pixel par pixel -- c'est le mode a utiliser si on veut un vrai
#     indice de confiance/transparence.
#   - XGB gere le NA nativement : CHAQUE pixel a une prediction TOUS les
#     jours (purete quasi toujours = 1, peu informatif comme indice de
#     confiance -- mais carte finale complete, sans aucun trou).
# Aucun réentraînement necessaire dans les deux cas : on recharge juste
# les modeles deja entraines par 11_run_training.R.
#
# ATTENTION AU VOLUME : un mois de 30 jours = 30 prédictions sur toute
# la grille par freq x schéma -- coûteux si beaucoup de schémas/mois
# sont sélectionnés. Réduis MONTHLY_COMPOSITE_SCHEMES/MONTHLY_COMPOSITE_MONTHS
# si besoin.
#
# Sorties, sous outputs_pipeline/monthly_composite/<freq>kHz/<model>/<schema>/ :
#   - monthly_mean_<YYYY-MM>.png, monthly_mean_avec_purete_<YYYY-MM>.png, monthly_purity_<YYYY-MM>.png
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

# Modele utilise pour le composite : "rf" (purete reellement variable,
# carte potentiellement incomplete -- cf. note ci-dessus) ou "xgb"
# (carte complete, purete peu informative).
MONTHLY_COMPOSITE_MODEL <- "rf"

# Schémas à traiter (par défaut seulement naive -- même logique que
# MULTIDATE_PREDICTION_SCHEMES, étends si besoin).
MONTHLY_COMPOSITE_SCHEMES <- c("naive_RS_80_20")

# Mois à traiter : NULL = tous les mois présents dans la grille
# multi-date (déduits automatiquement) ; sinon un vecteur de chaînes
# "YYYY-MM", ex. c("2018-01", "2018-02").
MONTHLY_COMPOSITE_MONTHS <- NULL

# ---- Prediction journaliere, moyenne des K folds -----------------------
# RF : seuls les pixels COMPLETS ce jour-la recoivent une prediction --
# le reste reste NA pour CE jour (pas pour tout le mois). C'est ce qui
# fait varier la purete mensuelle de facon informative.
predict_mean_rf_partial <- function(models, grid_df) {
  preds_full <- rep(NA_real_, nrow(grid_df))
  complete_idx <- stats::complete.cases(grid_df[, COVARIATES_ALL])
  if (any(complete_idx)) {
    grid_complete <- grid_df[complete_idx, , drop = FALSE]
    preds <- sapply(models, function(m) predict_rf(m, grid_complete))
    preds_full[complete_idx] <- rowMeans(preds)
  }
  preds_full
}
# XGB : gere le NA nativement -- une prediction pour TOUS les pixels,
# tous les jours (purete quasi toujours 1, cf. note en tete de fichier).
predict_mean_xgb <- function(models, grid_df, fod_levels) {
  preds <- sapply(models, function(m) predict_xgb(m, grid_df, fod_levels = fod_levels))
  rowMeans(preds)
}

day_ds <- load_grid_all_dates()
all_months <- format(day_ds$date, "%Y-%m")
months_to_process <- if (is.null(MONTHLY_COMPOSITE_MONTHS)) sort(unique(all_months)) else MONTHLY_COMPOSITE_MONTHS

cat(sprintf("Modele utilise pour le composite : %s\n", toupper(MONTHLY_COMPOSITE_MODEL)))
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

    models_path <- file.path(training_dir, paste0(freq, "kHz"), MONTHLY_COMPOSITE_MODEL, scheme_name, "models.rds")
    if (!file.exists(models_path)) {
      cat("  [!] modele", toupper(MONTHLY_COMPOSITE_MODEL), "introuvable pour", scheme_name,
          "-- avez-vous lance 11_run_training.R ?\n")
      next
    }
    models <- readRDS(models_path)

    out_dir <- file.path(out_root, paste0(freq, "kHz"), MONTHLY_COMPOSITE_MODEL, scheme_name)
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

    for (ym in months_to_process) {

      date_idx_month <- which(all_months == ym)
      if (length(date_idx_month) == 0) {
        cat("  [!] aucune date pour", ym, "-- saute.\n")
        next
      }
      n_days <- length(date_idx_month)
      cat(sprintf("  %s - %s - %s : %d jours\n", toupper(MONTHLY_COMPOSITE_MODEL), scheme_name, ym, n_days))

      # Matrice [n_pixels x n_days] des predictions journalieres
      n_pixels <- length(day_ds$lon) * length(day_ds$lat)
      daily_preds <- matrix(NA_real_, nrow = n_pixels, ncol = n_days)
      lon_vec <- lat_vec <- NULL

      for (j in seq_along(date_idx_month)) {
        extracted <- extract_grid_for_date(day_ds, date_idx_month[j], fod_levels)
        grid_all  <- extracted$grid
        if (is.null(lon_vec)) { lon_vec <- grid_all$lon; lat_vec <- grid_all$lat }
        daily_preds[, j] <- if (MONTHLY_COMPOSITE_MODEL == "rf") {
          predict_mean_rf_partial(models, grid_all)
        } else {
          predict_mean_xgb(models, grid_all, fod_levels)
        }
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

      lon_res <- if (length(unique(lon_vec)) > 1) stats::median(diff(sort(unique(lon_vec)))) else NA_real_
      lat_res <- if (length(unique(lat_vec)) > 1) stats::median(diff(sort(unique(lat_vec)))) else NA_real_

      # ---- carte moyenne mensuelle SEULE (sans purete) ----
      p_mean_only <- ggplot(composite_df, aes(x = lon, y = lat, fill = mean_pred)) +
        geom_tile(width = lon_res, height = lat_res)
      if (!is.null(basemap)) {
        p_mean_only <- p_mean_only +
          geom_sf(data = basemap, inherit.aes = FALSE, fill = "grey40", color = "grey20", linewidth = 0.2) +
          coord_sf(xlim = range(lon_vec), ylim = range(lat_vec), expand = FALSE)
      } else {
        p_mean_only <- p_mean_only + coord_quickmap()
      }
      p_mean_only <- p_mean_only +
        scale_fill_viridis_c(limits = lims) +
        theme_pipeline +
        labs(title = sprintf("NASC predit -- composite mensuel (moyenne, %s)", toupper(MONTHLY_COMPOSITE_MODEL)),
             subtitle = sprintf("%d kHz - %s - %s (%d jours)", freq, scheme_name, ym, n_days),
             x = "Longitude", y = "Latitude", fill = "log10(NASC)\nmoyen")
      ggsave(file.path(out_dir, sprintf("monthly_mean_%s.png", ym)), p_mean_only, width = 8, height = 6, dpi = 150)

      # ---- carte moyenne mensuelle, PURETE EN TRANSPARENCE (meme carte) ----
      # fill = valeur moyenne (couleur), alpha = purete (transparence) --
      # un pixel avec peu de jours valides apparait plus transparent/pale.
      p_mean_purity <- ggplot(composite_df, aes(x = lon, y = lat, fill = mean_pred, alpha = purity)) +
        geom_tile(width = lon_res, height = lat_res)
      if (!is.null(basemap)) {
        p_mean_purity <- p_mean_purity +
          geom_sf(data = basemap, inherit.aes = FALSE, fill = "grey40", color = "grey20", linewidth = 0.2) +
          coord_sf(xlim = range(lon_vec), ylim = range(lat_vec), expand = FALSE)
      } else {
        p_mean_purity <- p_mean_purity + coord_quickmap()
      }
      p_mean_purity <- p_mean_purity +
        scale_fill_viridis_c(limits = lims) +
        scale_alpha_continuous(range = c(0, 1), limits = c(0, 1)) +  # range=c(0,1) : purete=0 -> totalement transparent
        theme_pipeline +
        labs(title = sprintf("NASC predit -- composite mensuel (moyenne, %s)", toupper(MONTHLY_COMPOSITE_MODEL)),
             subtitle = sprintf("%d kHz - %s - %s (%d jours) -- transparence = purete", freq, scheme_name, ym, n_days),
             x = "Longitude", y = "Latitude", fill = "log10(NASC)\nmoyen", alpha = "Purete\n(0-1)")
      ggsave(file.path(out_dir, sprintf("monthly_mean_avec_purete_%s.png", ym)), p_mean_purity, width = 8, height = 6, dpi = 150)

      # ---- carte de purete seule (lecture precise, en complement de la
      #      transparence -- moins ambigu pour lire une valeur exacte) ----
      p_purity <- ggplot(composite_df, aes(x = lon, y = lat, fill = purity)) +
        geom_tile(width = lon_res, height = lat_res)
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
        labs(title = "Purete du composite mensuel (nb jours valides / nb jours du mois, par pixel)",
             subtitle = sprintf("%d kHz - %s - %s (%d jours) - modele %s",
                                 freq, scheme_name, ym, n_days, toupper(MONTHLY_COMPOSITE_MODEL)),
             x = "Longitude", y = "Latitude", fill = "Purete\n(0-1)")
      ggsave(file.path(out_dir, sprintf("monthly_purity_%s.png", ym)), p_purity, width = 8, height = 6, dpi = 150)

      cat(sprintf("    -> %d/%d pixels avec au moins 1 jour valide (%.1f%%) | purete moyenne = %.2f (mediane = %.2f)\n",
                  sum(n_valid > 0), n_pixels, 100 * sum(n_valid > 0) / n_pixels,
                  mean(purity, na.rm = TRUE), stats::median(purity, na.rm = TRUE)))
    }
  }
}

cat("\nComposite mensuel termine. Resultats dans :", normalizePath(out_root), "\n")
