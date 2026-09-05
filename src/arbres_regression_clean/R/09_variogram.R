# =====================================================================
# 09_variogram.R -- variogramme empirique des résidus détendancés
# =====================================================================
# Objectif : calibrer (ou au moins questionner) les distances de buffer
# spatial (100km/20km/5km) et temporel (1j) fixées "a dire d'expert"
# dans 00_config.R, en regardant la vraie portee d'autocorrelation des
# residus -- une fois le gradient a grande echelle (ex. latitudinal)
# retire, pour ne pas confondre tendance et autocorrelation locale (cf.
# discussion sur la stationnarite de 1er/2nd ordre).
#
# Deux etapes :
#   1) detrend_residuals()          -- retire une tendance simple
#      (regression lineaire sur lat/lon) de la variable d'interet.
#   2) compute_empirical_variogram() -- calcule la semi-variance par
#      classe de distance, sur un sous-echantillon (le calcul est en
#      O(n^2), infaisable sur la totalite des donnees).

# ---------------------------------------------------------------------
# Retire une tendance simple (regression lineaire) de `response_col`,
# et ajoute une colonne "residual_detrended" au data.frame.
# `trend_formula` par defaut = juste la latitude (le gradient identifie
# dans la discussion) ; passer "lat + lon" pour une tendance 2D complete.
# ---------------------------------------------------------------------
detrend_residuals <- function(df, response_col = RESPONSE_VAR, trend_formula = "lat") {
  form <- stats::as.formula(paste(response_col, "~", trend_formula))
  trend_model <- stats::lm(form, data = df)
  df$residual_detrended <- stats::residuals(trend_model)
  attr(df, "trend_model") <- trend_model
  df
}

# ---------------------------------------------------------------------
# Variogramme empirique : semi-variance = 0.5 * moyenne((z_i - z_j)^2),
# par classe ("lag") de distance entre paires de points. `coord_cols`
# peut etre 1 colonne (variogramme temporel, ex. "day_num") ou 2
# colonnes (variogramme spatial, ex. c("x_km","y_km")).
#
# Sous-echantillonne a `n_sample` points AVANT de calculer les distances
# par paires (le nombre de paires croit en O(n^2) -- 3000 points ->
# ~4.5M paires, deja beaucoup mais gerable en memoire ; augmenter
# n_sample degrade vite les performances).
# ---------------------------------------------------------------------
compute_empirical_variogram <- function(df, value_col = "residual_detrended",
                                         coord_cols = c("x_km", "y_km"),
                                         n_lags = 15, max_dist = NULL,
                                         n_sample = 3000, seed = 1) {
  df <- df[stats::complete.cases(df[, c(value_col, coord_cols)]), ]

  set.seed(seed)
  if (nrow(df) > n_sample) {
    df <- df[sample(seq_len(nrow(df)), n_sample), ]
  }

  coords <- as.matrix(df[, coord_cols, drop = FALSE])
  values <- df[[value_col]]

  dist_mat  <- as.matrix(stats::dist(coords))
  value_mat <- outer(values, values, function(a, b) (a - b)^2)

  # ne garder que la moitie superieure (hors diagonale) pour ne pas
  # compter chaque paire deux fois
  upper_idx <- upper.tri(dist_mat)
  dists     <- dist_mat[upper_idx]
  sq_diffs  <- value_mat[upper_idx]

  if (is.null(max_dist)) max_dist <- stats::quantile(dists, 0.5, na.rm = TRUE)  # moitie des paires, portee typique

  lag_breaks <- seq(0, max_dist, length.out = n_lags + 1)
  lag_bin    <- cut(dists, breaks = lag_breaks, include.lowest = TRUE, labels = FALSE)

  tibble(lag_bin = lag_bin, dist = dists, sq_diff = sq_diffs) %>%
    filter(!is.na(lag_bin)) %>%
    group_by(lag_bin) %>%
    summarise(
      lag_mid      = mean(dist),
      semivariance = 0.5 * mean(sq_diff),
      n_pairs      = dplyr::n(),
      .groups = "drop"
    ) %>%
    arrange(lag_mid)
}
