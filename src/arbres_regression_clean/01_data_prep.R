# =====================================================================
# 01_data_prep.R -- chargement + nettoyage
# =====================================================================
# Utilisé par TOUS les scripts (tuning, entraînement, prédiction), pour
# TOUS les modèles (CART, RF, XGB) et TOUS les schémas (naive, blocked).
# Garantit que les 3 familles de modèles voient exactement les mêmes
# lignes filtrées de la même façon (à l'exception du traitement des NA
# sur les covariables numériques, cf. `drop_na_numeric`).
#
# - CART et RF (ranger) ne gèrent pas nativement le manquant dans ce
#   pipeline -> drop_na_numeric = TRUE (on filtre les lignes incomplètes).
# - XGBoost apprend une direction par défaut à chaque split pour les NA
#   -> drop_na_numeric = FALSE (on garde les NA, seul fod/NASC sont filtrés).

load_and_clean <- function(freq,
                            drop_na_numeric = TRUE,
                            diurnal_period = DIURNAL_PERIOD) {

  datas <- readRDS(PATH_TEMPLATE(freq))
  datas <- datas[datas$day == diurnal_period, ]

  q <- quantile(datas$nasc, probs = c(0.05, 0.95), na.rm = TRUE)
  datas <- datas |> dplyr::filter(nasc >= q[1], nasc <= q[2])
  datas$nasc <- log10(datas$nasc)

  datas$fod <- as.factor(datas$fod)
  datas$fod[datas$fod == "NA"] <- NA
  datas$fod <- droplevels(datas$fod)

  df <- data.frame(
    NASC            = datas$nasc,
    year            = format(datas$time_nasc, "%Y"),
    time            = datas$time_nasc,
    lat             = datas$lat_nasc,
    lon             = datas$lon_nasc,
    fod             = datas$fod,
    ftle            = datas$ftle,
    total_chla      = datas$Chla,
    per_ratio_chla  = datas$Per_Chla,
    but_ratio_chla  = datas$But_Chla,
    fuco_ratio_chla = datas$Fuco_Chla,
    hex_ratio_chla  = datas$Hex_Chla,
    allo_ratio_chla = datas$Allo_Chla,
    zea_ratio_chla  = datas$Zea_Chla,
    chlb_ratio_chla = datas$Chlb_Chla
  )

  fod_levels <- levels(df$fod)  # sauvegardé AVANT filtrage -> réutilisé en prédiction

  if (drop_na_numeric) {
    df <- df |>
      dplyr::filter(
        if_all(all_of(COVARIATES_NUM), ~ !is.na(.)),
        !is.na(fod), fod != "NA"
      )
  } else {
    df <- df |> dplyr::filter(!is.na(fod), fod != "NA", is.finite(NASC))
  }

  # conversion lon/lat -> km (centrée), nécessaire pour le blocage spatial
  # et pour les distances géographiques dans les diagnostics
  lon0 <- mean(df$lon, na.rm = TRUE)
  lat0 <- mean(df$lat, na.rm = TRUE)
  km_per_deg_lat <- 110.574
  km_per_deg_lon <- 111.320 * cos(lat0 * pi / 180)
  df$x_km <- (df$lon - lon0) * km_per_deg_lon
  df$y_km <- (df$lat - lat0) * km_per_deg_lat

  cat(sprintf(
    "[%d kHz] %d observations apres filtrage (drop_na_numeric=%s)\n",
    freq, nrow(df), drop_na_numeric
  ))

  list(df = df, fod_levels = fod_levels, freq = freq)
}

# ---------------------------------------------------------------------
# Construction du jeu de covariables pour la grille de prédiction
# spatiale (date donnée), commun à tous les modèles.
# ---------------------------------------------------------------------
build_grid_covariates <- function(fod_levels, path_grid = PATH_GRID_DAY) {
  day_ds <- readRDS(path_grid)

  grid_points <- expand.grid(lon = day_ds$lon, lat = day_ds$lat)
  grid_points$ftle <- as.vector(day_ds$ftle)

  for (p in names(day_ds$pig)) {
    grid_points[[p]] <- as.vector(day_ds$pig[[p]])
  }

  grid_points$total_chla      <- grid_points$Chla
  grid_points$per_ratio_chla  <- grid_points$Per  / grid_points$Chla
  grid_points$but_ratio_chla  <- grid_points$But  / grid_points$Chla
  grid_points$fuco_ratio_chla <- grid_points$Fuco / grid_points$Chla
  grid_points$hex_ratio_chla  <- grid_points$Hex  / grid_points$Chla
  grid_points$allo_ratio_chla <- grid_points$Allo / grid_points$Chla
  grid_points$zea_ratio_chla  <- grid_points$Zea  / grid_points$Chla
  grid_points$chlb_ratio_chla <- grid_points$Chlb / grid_points$Chla

  # Inf/NaN (division par 0) -> NA pour que XGBoost les traite comme
  # manquants plutôt que comme des valeurs numériques aberrantes
  for (v in COVARIATES_NUM) {
    bad <- !is.finite(grid_points[[v]])
    grid_points[[v]][bad] <- NA
  }

  grid_points$fod <- formatC(as.vector(day_ds$fod), width = 2)
  grid_points$fod <- factor(grid_points$fod, levels = fod_levels)

  list(grid = grid_points, date = day_ds$date)
}
