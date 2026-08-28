# ============================================================
# 0. Packages
# ============================================================

library(gstat)
library(ggplot2)
library(dplyr)
library(patchwork)

# ============================================================
# 1. Conversion lon/lat -> km
# ------------------------------------------------------------
# Approximation planaire locale (comme pour l'analyse pigments).
# Valable sur l'etendue du transect ; si le transect couvre une
# tres large gamme de latitudes, une projection plus rigoureuse
# (ex. package sf) serait preferable.
# ============================================================

lonlat_to_km <- function(dat) {
  lat0 <- mean(dat$lat, na.rm = TRUE)
  dat %>% mutate(x_km = lon * 111.32 * cos(lat0 * pi / 180), y_km = lat * 111.32)
}

# ============================================================
# 2. Sv (dB) -> moyenne par ping
# ------------------------------------------------------------
# Le Sv est une grandeur logarithmique : la moyenne doit se faire
# en lineaire (10^(Sv/10)) puis etre reconvertie en dB, sinon le
# resultat est biaise. Traitement par blocs pour limiter le pic
# memoire (la matrice profiles fait ~230 millions de cellules).
# ============================================================

mean_sv_by_ping <- function(profiles_mat, chunk_size = 50000) {
  n <- nrow(profiles_mat)
  out_mean <- numeric(n)
  out_nvalid <- integer(n)
  for (start in seq(1, n, by = chunk_size)) {
    end <- min(start + chunk_size - 1, n)
    block <- profiles_mat[start:end, , drop = FALSE]
    lin <- 10^(block / 10)
    out_mean[start:end] <- 10 * log10(rowMeans(lin, na.rm = TRUE))
    out_nvalid[start:end] <- rowSums(!is.na(block))
  }
  list(value = out_mean, n_valid = out_nvalid)
}

# ============================================================
# 3. Construction des donnees ping-level (equivalent de make_spatial_data)
# ------------------------------------------------------------
# depth_range = c(min, max) pour restreindre la couche analysee
# (ex : couche epipelagique). NULL = toute la colonne d'eau.
# min_depth_obs = nb minimum de bins de profondeur valides par ping
# pour garder le ping (evite les moyennes calculees sur 1-2 valeurs).
# ============================================================

make_ping_data <- function(profiles, depth_range = NULL, min_depth_obs = 5, n_sample = 3000, day_value = NULL) {
  
  # Filtrage jour/nuit (ou tout autre sous-ensemble de jours) AVANT
  # toute autre operation, pour garder mat/lon/lat/time synchronises.
  if (!is.null(day_value)) {
    mask_day <- profiles$day %in% day_value
    mat <- profiles$profiles[mask_day, , drop = FALSE]
    lon <- profiles$lon[mask_day]
    lat <- profiles$lat[mask_day]
    time <- profiles$time[mask_day]
  } else {
    mat <- profiles$profiles
    lon <- profiles$lon
    lat <- profiles$lat
    time <- profiles$time
  }
  
  depth <- profiles$depth
  
  if (!is.null(depth_range)) {
    keep_cols <- depth >= depth_range[1] & depth <= depth_range[2]
    mat <- mat[, keep_cols, drop = FALSE]
  }
  n_depth <- ncol(mat)
  
  sv <- mean_sv_by_ping(mat)
  
  dat <- data.frame(
    lon = lon, lat = lat, time = time,
    value = sv$value, n_valid = sv$n_valid, frac_valid = sv$n_valid / n_depth
  )
  
  pct_na_global <- 1 - sum(dat$n_valid) / (nrow(dat) * n_depth)
  n_pings_total <- nrow(dat)
  
  dat <- dat %>% filter(is.finite(value), n_valid >= min_depth_obs)
  n_pings_kept <- nrow(dat)
  
  if (nrow(dat) > n_sample) {
    set.seed(123)
    dat <- dat %>% slice_sample(n = n_sample)
  }
  
  attr(dat, "coverage") <- data.frame(
    n_depth_bins = n_depth, pct_na_global = pct_na_global,
    n_pings_total = n_pings_total, n_pings_kept = n_pings_kept, n_pings_sampled = nrow(dat)
  )
  dat
}

# ============================================================
# 4. Diagnostic de couverture verticale (mesure de l'effet des NA)
# ============================================================

make_coverage_plot <- function(dat, label) {
  cov <- attr(dat, "coverage")
  subtitle <- sprintf(
    "%.0f%% de NA (profondeur) | %d/%d pings retenus (min_depth_obs) | %d echantillonnes",
    100 * cov$pct_na_global, cov$n_pings_kept, cov$n_pings_total, cov$n_pings_sampled
  )
  ggplot(dat, aes(x = lon, y = lat, color = n_valid)) +
    geom_point(size = 0.8, alpha = 0.7) +
    scale_color_viridis_c(name = "N bins\nvalides") +
    labs(title = paste("Couverture verticale -", label), subtitle = subtitle, x = "Longitude", y = "Latitude") +
    theme_minimal()
}

#  ============================================================
# 5. VARIOGRAMME DIRECTIONNEL
# ------------------------------------------------------------
# Attention : un transect echantillonne surtout la direction de la
# route du navire. Les directions peu suivies par le bateau seront
# mal renseignees (peu de paires) -> a interpreter avec prudence,
# contrairement a une grille reguliere couvrant tout l'espace.
# ============================================================

make_variogram <- function(dat, label, cutoff = NULL, alpha = c(0, 45, 90, 135), directional = TRUE) {
  
  dat <- lonlat_to_km(dat)
  if (nrow(dat) < 30) return(patchwork::wrap_elements(grid::textGrob(paste("Pas assez de données -", label))))
  
  if (is.null(cutoff)) {
    dx <- diff(range(dat$x_km)); dy <- diff(range(dat$y_km))
    cutoff <- 0.5 * sqrt(dx^2 + dy^2)
  }
  
  if (!directional) {
    # Variogramme omnidirectionnel : une seule courbe, toutes les
    # directions sont regroupees (pas d'argument alpha).
    v <- variogram(value ~ 1, locations = ~x_km + y_km, data = dat, cutoff = cutoff, width = cutoff / 15)
    return(
      ggplot(v, aes(x = dist, y = gamma)) +
        geom_point(size = 1.8) +
        geom_line() +
        labs(title = paste("Variogramme omnidirectionnel -", label), x = "Distance (km)", y = "Semi-variance") +
        theme_minimal()
    )
  }
  
  v <- variogram(value ~ 1, locations = ~x_km + y_km, data = dat, cutoff = cutoff, width = cutoff / 15, alpha = alpha)
  
  # Convention gstat : alpha en degres, sens horaire depuis le Nord.
  alpha_labels <- c("0" = "Méridional", "45" = "Diagonale NE-SO", "90" = "Zonal", "135" = "Diagonale NO-SE")
  v$direction <- factor(alpha_labels[as.character(v$dir.hor)], levels = alpha_labels)
  
  ggplot(v, aes(x = dist, y = gamma, color = direction)) +
    geom_point(size = 1.6) +
    geom_line() +
    labs(title = paste("Variogramme directionnel -", label), x = "Distance (km)", y = "Semi-variance", color = "Direction") +
    theme_minimal()
}

# ============================================================
# 6. CORRELOGRAMME SPATIAL DIRECTIONNEL
# ============================================================

make_correlogram_spatial <- function(dat, label, n_pairs = 50000, n_bins = 20, min_pairs_per_bin = 30, directional = TRUE) {
  
  dat <- lonlat_to_km(dat)
  n <- nrow(dat)
  if (n < 30) return(patchwork::wrap_elements(grid::textGrob(paste("Pas assez de données -", label))))
  
  set.seed(123)
  i <- sample(seq_len(n), n_pairs, replace = TRUE)
  j <- sample(seq_len(n), n_pairs, replace = TRUE)
  
  dx <- dat$x_km[i] - dat$x_km[j]
  dy <- dat$y_km[i] - dat$y_km[j]
  distance <- sqrt(dx^2 + dy^2)
  
  z1 <- dat$value[i]
  z2 <- dat$value[j]
  
  keep <- distance > 0 & is.finite(z1) & is.finite(z2)
  pairs <- data.frame(distance = distance[keep], z1 = z1[keep], z2 = z2[keep])
  pairs$bin <- cut(pairs$distance, breaks = n_bins, include.lowest = TRUE)
  
  if (!directional) {
    # Corrélogramme omnidirectionnel : toutes les directions regroupees.
    cor_df <- pairs %>%
      group_by(bin) %>%
      summarise(distance = mean(distance), correlation = cor(z1, z2), n = n(), .groups = "drop") %>%
      filter(n >= min_pairs_per_bin)
    
    return(
      ggplot(cor_df, aes(x = distance, y = correlation)) +
        geom_hline(yintercept = 0, linetype = "dashed") +
        geom_point(size = 1.8) +
        geom_line() +
        labs(title = paste("Correlogramme spatial omnidirectionnel -", label), x = "Distance (km)", y = "Corrélation") +
        ylim(-1, 1) +
        theme_minimal()
    )
  }
  
  angle <- (atan2(dy, dx) * 180 / pi) %% 180
  pairs$angle <- angle[keep]
  pairs$direction <- cut(
    pairs$angle, breaks = c(-Inf, 22.5, 67.5, 112.5, 157.5, Inf),
    labels = c("Zonal (E-O)", "Diagonale", "Méridional (N-S)", "Diagonale", "Zonal (E-O)")
  )
  
  cor_df <- pairs %>%
    group_by(direction, bin) %>%
    summarise(distance = mean(distance), correlation = cor(z1, z2), n = n(), .groups = "drop") %>%
    filter(n >= min_pairs_per_bin)
  
  ggplot(cor_df, aes(x = distance, y = correlation, color = direction)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point(size = 1.6) +
    geom_line() +
    labs(title = paste("Correlogramme spatial -", label), x = "Distance (km)", y = "Corrélation", color = "Direction") +
    ylim(-1, 1) +
    theme_minimal()
}

# ============================================================
# 7. CORRELOGRAMME TEMPOREL
# ------------------------------------------------------------
# Remplace l'ACF classique (qui exige un pas de temps regulier) par
# un correlogramme par paires, comme pour l'espace : compatible avec
# l'echantillonnage irregulier des pings (quelques secondes d'ecart,
# eventuels trous). Bins bases sur les quantiles (et non largeur
# fixe) car la distribution des ecarts de temps entre paires
# aleatoires est tres asymetrique (beaucoup de grands ecarts, peu de
# petits) sur un transect qui dure plusieurs jours.
#
# ATTENTION : le navire se deplacant en continu, un faible ecart de
# temps implique presque toujours un faible ecart spatial. Ce
# correlogramme reflete donc un melange de decorrelation temporelle
# et spatiale, pas un effet du temps pur.
# ============================================================

make_correlogram_temporal <- function(dat, label, n_pairs = 50000, n_bins = 20, min_pairs_per_bin = 30, time_unit = "hours") {
  
  n <- nrow(dat)
  if (n < 30) return(patchwork::wrap_elements(grid::textGrob(paste("Pas assez de données -", label))))
  
  set.seed(123)
  i <- sample(seq_len(n), n_pairs, replace = TRUE)
  j <- sample(seq_len(n), n_pairs, replace = TRUE)
  
  dt <- abs(as.numeric(difftime(dat$time[i], dat$time[j], units = time_unit)))
  z1 <- dat$value[i]
  z2 <- dat$value[j]
  
  keep <- dt > 0 & is.finite(z1) & is.finite(z2)
  pairs <- data.frame(dt = dt[keep], z1 = z1[keep], z2 = z2[keep])
  
  breaks <- unique(quantile(pairs$dt, probs = seq(0, 1, length.out = n_bins + 1)))
  pairs$bin <- cut(pairs$dt, breaks = breaks, include.lowest = TRUE)
  
  cor_df <- pairs %>%
    group_by(bin) %>%
    summarise(dt = mean(dt), correlation = cor(z1, z2), n = n(), .groups = "drop") %>%
    filter(n >= min_pairs_per_bin)
  
  ggplot(cor_df, aes(x = dt, y = correlation)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point(size = 1.8) +
    geom_line() +
    labs(title = paste("Correlogramme temporel -", label), x = paste0("Écart de temps (", time_unit, ")"), y = "Corrélation") +
    ylim(-1, 1) +
    theme_minimal()
}

# ============================================================
# 8. FONCTION QUI PRODUIT LES 4 PLOTS
# ============================================================

make_profile_plots <- function(profiles, label = "Sv", depth_range = NULL, min_depth_obs = 5, n_sample = 3000, time_unit = "hours", day_value = NULL) {
  
  cat("Traitement :", label, "\n")
  if (is.null(day_value)) {
    day_label <- ""
  } else if (day_value == 3) {
    day_label <- "day"
  } else if (day_value == 1) {
    day_label <- "night"
  } else {
    day_label <- ""
  }
  dat <- make_ping_data(profiles, depth_range = depth_range, min_depth_obs = min_depth_obs, n_sample = n_sample, day_value = day_value)
  cov <- attr(dat, "coverage")
  cat(sprintf("  -> %.0f%% NA (profondeur) | %d/%d pings retenus | %d échantillonnés\n",
              100 * cov$pct_na_global, cov$n_pings_kept, cov$n_pings_total, cov$n_pings_sampled))
  
  p0 <- make_coverage_plot(dat, label)
  p1 <- make_variogram(dat, label, directional = FALSE)
  p2 <- make_correlogram_spatial(dat, label, directional = FALSE)
  p3 <- make_correlogram_temporal(dat, label, time_unit = time_unit)
  
  plot_acoustique <- (p0 + p1) / (p2 + p3) +
    patchwork::plot_annotation(title = paste("Dépendances spatio-temporelles -", label, freq, "kHz", day_label))
  
  ggsave(
    filename = paste0("C:/Users/mmolinet/elisou_ta_stagiaire_pref/MIO_internship_III/figures/diag_autocorrelation_spatio_temp/diag_autocorrelation_spatio_temp_Sv_", freq, "kHz_", day_label, ".png"),
    plot = plot_acoustique,
    width = 14, height = 10, dpi = 300, units = "in"
  )
  
  plot_acoustique
}

# ============================================================
# 9. EXECUTION
# ------------------------------------------------------------
# depth_range = NULL -> moyenne integree sur toute la colonne d'eau.
# Pour une couche specifique, ex : depth_range = c(50, 200).
# ============================================================
freqs <- c(18) # , 38, 70, 120, 200
for (freq in freqs){
  for (day_val in c(NULL, 3, 1)){
    path_sv <- paste0("F:/data_elise/sv_cropped/sv_cropped_all_years/Sv_2018_2021_2022_2023_", freq, "kHz.rds")
    profiles <- readRDS(path_sv)
    plot_acoustique <- make_profile_plots(profiles, label = "Sv", depth_range = NULL, min_depth_obs = 5, n_sample = 3000, time_unit = "hours", day_value = day_val)
    plot_acoustique
  }
}

