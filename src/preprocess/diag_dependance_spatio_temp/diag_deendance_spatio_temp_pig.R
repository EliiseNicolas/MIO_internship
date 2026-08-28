# ============================================================
# 0. Packages
# ============================================================

library(gstat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# ============================================================
# 1. Pigments à analyser
# ============================================================

pigments <- c("Chla", "Per", "But", "Fuco", "Hex", "Allo", "Zea", "Chlb", "DvChla")
path_pig <- "F:/data_elise/pigmeann/pigments_2018_2021_2022_2023_crop.rds"
pigs <- readRDS(path_pig)
print(range(pigs$lat))
print(range(pigs$lon))
freq <- 18
path_sv <- paste0("F:/data_elise/sv_cropped/sv_cropped_all_years/Sv_2018_2021_2022_2023_", freq, "kHz.rds")
profiles <- readRDS(path_sv)
str(profiles)

# ============================================================
# 2. Helper : orientation de l'array (time, lon, lat)
# ------------------------------------------------------------
# Le code original supposait que dim(x) = (time, lon, lat) sans le
# vérifier. Si jamais l'array est (time, lat, lon), expand.grid(lon=...,
# lat=...) + as.vector(z) désaligne silencieusement les valeurs avec
# les coordonnées. On détecte et corrige automatiquement.
# ============================================================

fix_lonlat_order <- function(x, pigs) {
  d <- dim(x)
  n_lon <- length(pigs$lon)
  n_lat <- length(pigs$lat)
  if (d[2] == n_lon && d[3] == n_lat) return(x)
  if (d[2] == n_lat && d[3] == n_lon) return(aperm(x, c(1, 3, 2)))
  stop("Dimensions de l'array incompatibles avec pigs$lon / pigs$lat.")
}

# ============================================================
# 3. Construction des données spatiales + suivi de la couverture NA
# ------------------------------------------------------------
# ~80-90% de NA (nuages) : on calcule, pour chaque pixel, le nombre
# d'observations valides utilisées dans la moyenne temporelle
# (n_valid) et on exige un minimum (min_obs) pour garder le pixel,
# afin d'éviter des moyennes instables calculées sur 1 ou 2 valeurs.
# ============================================================

make_spatial_data <- function(pigment, pigs, date_index = NULL, n_sample = 3000, min_obs = 3) {
  
  x <- pigs[[paste0("c_cond_", pigment)]]
  x <- fix_lonlat_order(x, pigs)
  if (!is.null(date_index)) x <- x[date_index, , , drop = FALSE]
  
  n_time <- dim(x)[1]
  n_valid <- apply(x, c(2, 3), function(v) sum(!is.na(v)))
  z <- apply(x, c(2, 3), mean, na.rm = TRUE)
  z[is.nan(z)] <- NA
  
  dat <- expand.grid(lon = pigs$lon, lat = pigs$lat)
  dat$value <- as.vector(z)
  dat$n_valid <- as.vector(n_valid)
  dat$frac_valid <- dat$n_valid / n_time
  
  pct_na_global <- 1 - sum(dat$n_valid) / (n_time * nrow(dat))
  n_pixels_total <- nrow(dat)
  
  dat <- dat %>% filter(is.finite(value), n_valid >= min_obs)
  n_pixels_kept <- nrow(dat)
  
  if (nrow(dat) > n_sample) {
    set.seed(123)
    dat <- dat %>% slice_sample(n = n_sample)
  }
  
  attr(dat, "coverage") <- data.frame(
    pigment = pigment, n_time = n_time, pct_na_global = pct_na_global,
    n_pixels_total = n_pixels_total, n_pixels_kept = n_pixels_kept,
    n_pixels_sampled = nrow(dat)
  )
  dat
}

# ============================================================
# 4. Diagnostic de couverture (mesure de l'effet des NA)
# ============================================================

make_coverage_plot <- function(dat, pigment_name) {
  cov <- attr(dat, "coverage")
  subtitle <- sprintf(
    "%.0f%% de NA (nuages) | %d/%d pixels retenus (min_obs) | %d échantillonnés",
    100 * cov$pct_na_global, cov$n_pixels_kept, cov$n_pixels_total, cov$n_pixels_sampled
  )
  ggplot(dat, aes(x = lon, y = lat, color = n_valid)) +
    geom_point(size = 1.2) +
    scale_color_viridis_c(name = "N obs\nvalides") +
    labs(title = paste("Couverture temporelle -", pigment_name), subtitle = subtitle, x = "Longitude", y = "Latitude") +
    theme_minimal()
}

# ============================================================
# 5. Conversion lon/lat -> kilomètres
# ============================================================

lonlat_to_km <- function(dat) {
  lat0 <- mean(dat$lat, na.rm = TRUE)
  dat %>% mutate(x_km = lon * 111.32 * cos(lat0 * pi / 180), y_km = lat * 111.32)
}

# ============================================================
# 5bis. Palette et libellés partagés pour les directions
# ------------------------------------------------------------
# Utilises a la fois par make_variogram et make_correlogram, pour
# que "Meridional", "Zonal" etc. aient la meme couleur et le meme
# texte dans tous les plots.
# ============================================================

direction_levels <- c("Méridional (N-S)", "Diagonale NE-SO", "Zonal (E-O)", "Diagonale NO-SE")
direction_colors <- c(
  "Méridional (N-S)" = "#1b9e77",
  "Diagonale NE-SO"  = "#7570b3",
  "Zonal (E-O)"      = "#d95f02",
  "Diagonale NO-SE"  = "#e7298a"
)


#============================================================
  # 6. VARIOGRAMME DIRECTIONNEL (anisotropie)
  # ------------------------------------------------------------
# alpha = c(0, 45, 90, 135) : 0/180 = E-W (zonal), 90 = N-S (méridional).
# Une forte anisotropie lon/lat se traduit par des courbes très
# différentes selon la direction.
# ============================================================

make_variogram <- function(dat, pigment_name, cutoff = NULL, alpha = c(0, 45, 90, 135)) {
  
  dat <- lonlat_to_km(dat)
  if (nrow(dat) < 30) return(patchwork::wrap_elements(grid::textGrob(paste("Pas assez de données -", pigment_name))))
  
  if (is.null(cutoff)) {
    dx <- diff(range(dat$x_km)); dy <- diff(range(dat$y_km))
    cutoff <- 0.5 * sqrt(dx^2 + dy^2)
  }
  
  v <- variogram(value ~ 1, locations = ~x_km + y_km, data = dat, cutoff = cutoff, width = cutoff / 15, alpha = alpha)
  
  # gstat mesure alpha en degres, sens horaire depuis le Nord (axe y) :
  # 0 = Nord-Sud (meridional), 90 = Est-Ouest (zonal).
  alpha_labels <- c("0" = "Méridional (N-S)", "45" = "Diagonale NE-SO", "90" = "Zonal (E-O)", "135" = "Diagonale NO-SE")
  v$direction <- factor(alpha_labels[as.character(v$dir.hor)], levels = direction_levels)
  
  ggplot(v, aes(x = dist, y = gamma, color = direction)) +
    geom_point(size = 1.6) +
    geom_line() +
    scale_color_manual(values = direction_colors, drop = FALSE) +
    labs(title = paste("Variogramme directionnel -", pigment_name), x = "Distance (km)", y = "Semi-variance", color = "Direction") +
    theme_minimal()
}

# ============================================================
# 7. CORRELOGRAMME SPATIAL DIRECTIONNEL
# ------------------------------------------------------------
# On classe chaque paire par direction (zonal / méridional / diagonal)
# avant de la classer par distance, pour visualiser l'anisotropie.
# ============================================================

make_correlogram <- function(dat, pigment_name, n_pairs = 50000, n_bins = 20, min_pairs_per_bin = 30) {
  
  dat <- lonlat_to_km(dat)
  n <- nrow(dat)
  if (n < 30) return(patchwork::wrap_elements(grid::textGrob(paste("Pas assez de données -", pigment_name))))
  
  set.seed(123)
  i <- sample(seq_len(n), n_pairs, replace = TRUE)
  j <- sample(seq_len(n), n_pairs, replace = TRUE)
  
  dx <- dat$x_km[i] - dat$x_km[j]
  dy <- dat$y_km[i] - dat$y_km[j]
  distance <- sqrt(dx^2 + dy^2)
  angle <- (atan2(dy, dx) * 180 / pi) %% 180
  
  z1 <- dat$value[i]
  z2 <- dat$value[j]
  
  keep <- distance > 0 & is.finite(z1) & is.finite(z2)
  pairs <- data.frame(distance = distance[keep], angle = angle[keep], z1 = z1[keep], z2 = z2[keep])
  
  # 4 categories identiques au variogramme (au lieu de fusionner les
  # deux diagonales en une seule "Diagonal").
  pairs$direction <- cut(
    pairs$angle, breaks = c(-Inf, 22.5, 67.5, 112.5, 157.5, Inf),
    labels = c("Zonal (E-O)", "Diagonale NE-SO", "Méridional (N-S)", "Diagonale NO-SE", "Zonal (E-O)")
  )
  pairs$direction <- factor(pairs$direction, levels = direction_levels)
  pairs$bin <- cut(pairs$distance, breaks = n_bins, include.lowest = TRUE)
  
  cor_df <- pairs %>%
    group_by(direction, bin) %>%
    summarise(distance = mean(distance), correlation = cor(z1, z2), n = n(), .groups = "drop") %>%
    filter(n >= min_pairs_per_bin)
  
  ggplot(cor_df, aes(x = distance, y = correlation, color = direction)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point(size = 1.6) +
    geom_line() +
    scale_color_manual(values = direction_colors, drop = FALSE) +
    labs(title = paste("Correlogramme spatial directionnel -", pigment_name), x = "Distance (km)", y = "Corrélation", color = "Direction") +
    ylim(-1, 1) +
    theme_minimal()
}

# ============================================================
# 8. AUTOCORRÉLATION TEMPORELLE (ACF)
# ============================================================

make_acf <- function(pigment, pigs, lag_max = 60) {
  
  x <- pigs[[paste0("c_cond_", pigment)]]
  x <- fix_lonlat_order(x, pigs)
  
  temporal_mean <- apply(x, 1, mean, na.rm = TRUE)
  temporal_mean[is.nan(temporal_mean)] <- NA
  
  if (sum(!is.na(temporal_mean)) < lag_max + 5) {
    return(patchwork::wrap_elements(grid::textGrob(paste("Série trop courte -", pigment))))
  }
  
  acf_obj <- acf(temporal_mean, lag.max = lag_max, na.action = na.pass, plot = FALSE)
  acf_df <- data.frame(lag = acf_obj$lag[, 1, 1], acf = acf_obj$acf[, 1, 1])
  
  ggplot(acf_df, aes(x = lag, y = acf)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_segment(aes(xend = lag, yend = 0)) +
    labs(title = paste("Autocorrélation temporelle -", pigment), x = "Lag (jours)", y = "ACF") +
    theme_minimal()
}

# ============================================================
# 9. FONCTION QUI PRODUIT LES 4 PLOTS (couverture + variogramme + correlogramme + ACF)
# ============================================================

make_four_plots <- function(pigment, pigs, n_sample = 3000, min_obs = 3) {
  
  cat("Traitement :", pigment, "\n")
  
  spatial_data <- make_spatial_data(pigment = pigment, pigs = pigs, n_sample = n_sample, min_obs = min_obs)
  cov <- attr(spatial_data, "coverage")
  cat(sprintf("  -> %.0f%% NA global | %d/%d pixels retenus (min_obs=%d) | %d échantillonnés\n",
              100 * cov$pct_na_global, cov$n_pixels_kept, cov$n_pixels_total, min_obs, cov$n_pixels_sampled))
  
  p0 <- make_coverage_plot(spatial_data, pigment)
  p1 <- make_variogram(spatial_data, pigment)
  p2 <- make_correlogram(spatial_data, pigment)
  p3 <- make_acf(pigment, pigs)
  
  # list(plot = (p0 + p1) / (p2 + p3), coverage = cov)
  plot_acoustique <- (p0 + p1) / (p2 + p3) +
    patchwork::plot_annotation(title = paste("Dépendances spatio-temporelles -", pigment))
  
  ggsave(
    filename = paste0("C:/Users/mmolinet/elisou_ta_stagiaire_pref/MIO_internship_III/figures/diag_autocorrelation_spatio_temp/diag_autocorrelation_spatio_temp_", pigment, ".png"),
    plot = plot_acoustique,
    width = 14, height = 10, dpi = 300, units = "in"
  )
  
  plot_acoustique
}

# ============================================================
# 10. BOUCLE SUR TOUS LES PIGMENTS (robuste aux erreurs individuelles)
# ============================================================

plots_pigments <- list()
coverage_summary <- list()

for (pig in pigments) {
  res <- tryCatch(
    make_four_plots(pigment = pig, pigs = pigs, n_sample = 5000, min_obs = 10),
    error = function(e) { cat("  ERREUR pour", pig, ":", conditionMessage(e), "\n"); NULL }
  )
  if (!is.null(res)) {
    plots_pigments[[pig]] <- res$plot
    coverage_summary[[pig]] <- res$coverage
  }
}

# Table récapitulative de l'effet des NA sur tous les pigments
coverage_summary_df <- bind_rows(coverage_summary)
print(coverage_summary_df)

# ============================================================
# 11. AFFICHAGE
# ============================================================

plots_pigments$Chla
plots_pigments$Per
plots_pigments$But
plots_pigments$Fuco
plots_pigments$Hex
plots_pigments$Allo
plots_pigments$Zea
plots_pigments$Chlb
plots_pigments$DvChla

# ============================================================
# 12. VARIOGRAMME + CORRELOGRAMME PAR AXE (LAT vs LON), TOUS PIGMENTS
# ------------------------------------------------------------
# alpha = 0 (gstat) -> Nord-Sud (latitude) ; alpha = 90 -> Est-Ouest (longitude).
# tol.hor = 22.5 restreint aux paires quasi-alignees sur l'axe (meme
# tolerance que la classification "Zonal"/"Meridional" du correlogramme).
# Standardisation (z-score) par pigment pour rendre les courbes comparables.
# ============================================================

get_variogram_axis_df <- function(dat, pigment_name, axis = c("lat", "lon"), cutoff = NULL, tol.hor = 22.5) {
  axis <- match.arg(axis)
  dat <- lonlat_to_km(dat)
  if (nrow(dat) < 30) return(NULL)
  if (is.null(cutoff)) {
    dx <- diff(range(dat$x_km)); dy <- diff(range(dat$y_km))
    cutoff <- 0.5 * sqrt(dx^2 + dy^2)
  }
  alpha_val <- if (axis == "lat") 0 else 90
  v <- variogram(value ~ 1, locations = ~x_km + y_km, data = dat, cutoff = cutoff, width = cutoff / 15, alpha = alpha_val, tol.hor = tol.hor)
  data.frame(pigment = pigment_name, dist = v$dist, gamma = v$gamma)
}

get_correlogram_axis_df <- function(dat, pigment_name, axis = c("lat", "lon"), n_pairs = 50000, n_bins = 20, min_pairs_per_bin = 30, tol.hor = 22.5) {
  axis <- match.arg(axis)
  dat <- lonlat_to_km(dat)
  n <- nrow(dat)
  if (n < 30) return(NULL)
  set.seed(123)
  i <- sample(seq_len(n), n_pairs, replace = TRUE)
  j <- sample(seq_len(n), n_pairs, replace = TRUE)
  dx <- dat$x_km[i] - dat$x_km[j]; dy <- dat$y_km[i] - dat$y_km[j]
  distance <- sqrt(dx^2 + dy^2)
  angle <- (atan2(dy, dx) * 180 / pi) %% 180
  z1 <- dat$value[i]; z2 <- dat$value[j]
  
  in_axis <- if (axis == "lat") {
    angle >= (90 - tol.hor) & angle <= (90 + tol.hor)          # meridional (N-S)
  } else {
    angle <= tol.hor | angle >= (180 - tol.hor)                 # zonal (E-O)
  }
  
  keep <- distance > 0 & is.finite(z1) & is.finite(z2) & in_axis
  pairs <- data.frame(distance = distance[keep], z1 = z1[keep], z2 = z2[keep])
  if (nrow(pairs) < min_pairs_per_bin) return(NULL)
  pairs$bin <- cut(pairs$distance, breaks = n_bins, include.lowest = TRUE)
  
  cor_df <- pairs %>%
    group_by(bin) %>%
    summarise(distance = mean(distance), correlation = cor(z1, z2), n = n(), .groups = "drop") %>%
    filter(n >= min_pairs_per_bin)
  
  data.frame(pigment = pigment_name, distance = cor_df$distance, correlation = cor_df$correlation)
}

# ============================================================
# 13. BOUCLE TOUS PIGMENTS, PAR AXE
# ============================================================

variogram_lat_all <- list(); variogram_lon_all <- list()
correlogram_lat_all <- list(); correlogram_lon_all <- list()

for (pig in pigments) {
  sp <- tryCatch(make_spatial_data(pigment = pig, pigs = pigs, n_sample = 3000, min_obs = 3), error = function(e) NULL)
  if (is.null(sp) || nrow(sp) < 30) next
  sp$value <- as.numeric(scale(sp$value))  # standardisation z-score
  
  variogram_lat_all[[pig]] <- get_variogram_axis_df(sp, pig, axis = "lat")
  variogram_lon_all[[pig]] <- get_variogram_axis_df(sp, pig, axis = "lon")
  correlogram_lat_all[[pig]] <- get_correlogram_axis_df(sp, pig, axis = "lat")
  correlogram_lon_all[[pig]] <- get_correlogram_axis_df(sp, pig, axis = "lon")
}

variogram_lat_df <- bind_rows(variogram_lat_all)
variogram_lon_df <- bind_rows(variogram_lon_all)
correlogram_lat_df <- bind_rows(correlogram_lat_all)
correlogram_lon_df <- bind_rows(correlogram_lon_all)

# ============================================================
# 13bis. Palette de couleurs fixe pour les pigments
# ------------------------------------------------------------
# Necessaire pour garantir la meme couleur par pigment sur les 4
# plots avant de collecter une legende unique (sinon un pigment
# absent d'un sous-ensemble peut decaler les couleurs des autres).
# ============================================================

pigment_colors <- setNames(scales::hue_pal()(length(pigments)), pigments)
# ============================================================
# 14. PLOTS
# ============================================================

p_variogram_lat <- ggplot(variogram_lat_df, aes(x = dist, y = gamma, color = pigment)) +
  geom_line(alpha = 0.6) + geom_point(size = 1, alpha = 0.6) +
  geom_smooth(aes(group = 1), color = "black", se = FALSE, linewidth = 1.2, method = "loess") +
  scale_color_manual(values = pigment_colors, limits = pigments) +
  labs(title = "Variogramme - axe Latitude (Méridional)", x = "Distance (km)", y = "Semi-variance (z-score)", color = "Pigment") +
  theme_minimal()

p_variogram_lon <- ggplot(variogram_lon_df, aes(x = dist, y = gamma, color = pigment)) +
  geom_line(alpha = 0.6) + geom_point(size = 1, alpha = 0.6) +
  geom_smooth(aes(group = 1), color = "black", se = FALSE, linewidth = 1.2, method = "loess") +
  scale_color_manual(values = pigment_colors, limits = pigments) +
  labs(title = "Variogramme - axe Longitude (Zonal)", x = "Distance (km)", y = "Semi-variance (z-score)", color = "Pigment") +
  theme_minimal()

p_correlogram_lat <- ggplot(correlogram_lat_df, aes(x = distance, y = correlation, color = pigment)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(alpha = 0.6) + geom_point(size = 1, alpha = 0.6) +
  geom_smooth(aes(group = 1), color = "black", se = FALSE, linewidth = 1.2, method = "loess") +
  scale_color_manual(values = pigment_colors, limits = pigments) +
  labs(title = "Correlogramme - axe Latitude (Méridional)", x = "Distance (km)", y = "Corrélation", color = "Pigment") +
  ylim(-1, 1) + theme_minimal()

p_correlogram_lon <- ggplot(correlogram_lon_df, aes(x = distance, y = correlation, color = pigment)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(alpha = 0.6) + geom_point(size = 1, alpha = 0.6) +
  geom_smooth(aes(group = 1), color = "black", se = FALSE, linewidth = 1.2, method = "loess") +
  scale_color_manual(values = pigment_colors, limits = pigments) +
  labs(title = "Correlogramme - axe Longitude (Zonal)", x = "Distance (km)", y = "Corrélation", color = "Pigment") +
  ylim(-1, 1) + theme_minimal()

# guides = "collect" fusionne les 4 legendes identiques en une seule,
# affichee une fois sur le cote (position par defaut = "right").
p_combined_axes <- (p_variogram_lat + p_variogram_lon) / (p_correlogram_lat + p_correlogram_lon) +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(title = "Dépendance spatiale par axe - tous pigments") &
  theme(legend.position = "right")

p_combined_axes

ggsave(
  filename = "C:/Users/mmolinet/elisou_ta_stagiaire_pref/MIO_internship_III/figures/diag_autocorrelation_spatio_temp/diag_autocorrelation_lat_lon_tous_pigments.png",
  plot = p_combined_axes, width = 16, height = 10, dpi = 300, units = "in"
)
