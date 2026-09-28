# =====================================================================
# 19_run_variogram_analysis.R -- SCRIPT 9 : variogrammes empiriques
# =====================================================================
# Calcule un variogramme SPATIAL et un variogramme TEMPOREL des residus
# detendances (NASC ~ latitude), pour juger si les buffers actuels
# (SPATIAL_RESOLUTIONS$buffer_km, TEMPORAL_RESOLUTIONS$buffer_days dans
# 00_config.R) sont coherents avec la vraie portee d'autocorrelation, ou
# fixes trop bas/trop haut "a dire d'expert". Les buffers actuels sont
# superposes sur les plots (ligne verticale rouge) pour comparaison
# visuelle directe.
#
# Ne re-entraine aucun modele -- calcul purement descriptif sur les
# donnees d'entrainement, independant de CART/RF/XGB.
#
# Sorties, sous outputs_pipeline/variogram/<freq>kHz/ :
#   - variogram_spatial.csv, variogram_spatial.png
#   - variogram_temporal.csv, variogram_temporal.png

source("R/00_config.R")
source("R/01_data_prep.R")
source("R/02_folds.R")
source("R/05_plots.R")
source("R/09_variogram.R")

out_root <- path_out("variogram")
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

for (freq in FREQS) {

  cat("\n============================================================\n")
  cat("VARIOGRAM ANALYSIS -- FREQUENCE :", freq, "kHz\n")
  cat("============================================================\n")

  out_dir <- file.path(out_root, paste0(freq, "kHz"))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  prep <- load_and_clean(freq, drop_na_numeric = TRUE)
  df   <- detrend_residuals(prep$df, response_col = RESPONSE_VAR, trend_formula = "lat")

  cat(sprintf("  Tendance NASC ~ lat : R2 = %.3f (variance expliquee par le seul gradient latitudinal)\n",
              summary(attr(df, "trend_model"))$r.squared))

  # ---- Variogramme SPATIAL (residus detendances, coord. x_km/y_km) ----
  variogram_spatial <- compute_empirical_variogram(
    df, value_col = "residual_detrended", coord_cols = c("x_km", "y_km"), n_lags = 15
  )
  write.csv(variogram_spatial, file.path(out_dir, "variogram_spatial.csv"), row.names = FALSE)

  # superpose les 3 buffers spatiaux actuels (100/20/5 km par defaut)
  p_spatial <- plot_variogram(variogram_spatial, subtitle = sprintf("%d kHz - NASC detendance sur lat", freq),
                               x_label = "Distance spatiale (km)")
  for (r in SPATIAL_RESOLUTIONS) {
    p_spatial <- p_spatial +
      geom_vline(xintercept = r$buffer_km, linetype = "dotted", color = "grey40") +
      annotate("text", x = r$buffer_km, y = min(variogram_spatial$semivariance), vjust = 0,
               label = paste0(r$label, " buffer=", r$buffer_km, "km"), angle = 90, size = 2.8, color = "grey40")
  }
  ggsave(file.path(out_dir, "variogram_spatial.png"), p_spatial, width = 9, height = 6, dpi = 150)
  cat("  -> variogram_spatial.png (buffers actuels superposes en pointilles gris)\n")

  # ---- Variogramme TEMPOREL (residus detendances, coord. = jour) ----
  df$day_num <- as.numeric(as.Date(df$time))
  variogram_temporal <- compute_empirical_variogram(
    df, value_col = "residual_detrended", coord_cols = "day_num", n_lags = 15
  )
  write.csv(variogram_temporal, file.path(out_dir, "variogram_temporal.csv"), row.names = FALSE)

  p_temporal <- plot_variogram(variogram_temporal, subtitle = sprintf("%d kHz - NASC detendance sur lat", freq),
                                x_label = "Distance temporelle (jours)")
  for (r in TEMPORAL_RESOLUTIONS) {
    p_temporal <- p_temporal +
      geom_vline(xintercept = r$buffer_days, linetype = "dotted", color = "grey40") +
      annotate("text", x = r$buffer_days, y = min(variogram_temporal$semivariance), vjust = 0,
               label = paste0(r$label, " buffer=", r$buffer_days, "j"), angle = 90, size = 2.8, color = "grey40")
  }
  ggsave(file.path(out_dir, "variogram_temporal.png"), p_temporal, width = 9, height = 6, dpi = 150)
  cat("  -> variogram_temporal.png (buffer actuel superpose en pointille gris)\n")
}

cat("\n============================================================\n")
cat("LECTURE DU VARIOGRAMME : la semi-variance doit normalement CROITRE\n")
cat("avec la distance puis PLAFONNER (le 'palier'/sill) a une distance dite\n")
cat("'portee' (range) -- c'est la distance au-dela de laquelle deux points\n")
cat("ne sont plus autocorreles. Si le buffer actuel (ligne pointillee) est\n")
cat("BIEN AVANT que la courbe plafonne, il est probablement TROP COURT\n")
cat("(fuite residuelle au-dela du buffer). S'il est BIEN APRES, il est\n")
cat("probablement TROP LARGE (train inutilement reduit).\n")
cat("============================================================\n")

cat("\nAnalyse variogramme terminee. Resultats dans :", normalizePath(out_root), "\n")
