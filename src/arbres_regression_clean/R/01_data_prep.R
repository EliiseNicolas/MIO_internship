# =====================================================================
# 01_data_prep.R -- chargement + nettoyage
# =====================================================================
# Utilisé par TOUS les scripts (tuning, entraînement, prédiction), pour
# TOUS les modèles (CART, RF, XGB, randomForestSRC) et TOUS les schémas
# (naive, blocked). Garantit que tous voient exactement les mêmes lignes
# filtrées de la même façon (à l'exception du traitement des NA sur les
# covariables numériques, cf. `drop_na_numeric`).
#
# - CART et RF (ranger) ne gèrent pas nativement le manquant dans ce
#   pipeline -> drop_na_numeric = TRUE (on filtre les lignes incomplètes).
# - XGBoost et randomForestSRC gèrent le manquant nativement (direction
#   par défaut à chaque split pour XGB, imputation pour rfsrc)
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

  # ---------------------------------------------------------------------
  # Structure de données mise à jour (colonnes déjà normalisées, plus de
  # calcul de ratio manuel) : NASC vient de `nasc`, lat/lon/time restent
  # `*_nasc` (les colonnes `*_fod`, `*_pig`, `*_ftle` sont les positions
  # -- potentiellement différentes -- des capteurs sources et ne sont
  # PAS utilisées comme coordonnées du point NASC lui-même).
  # ---------------------------------------------------------------------
  df <- data.frame(
    NASC          = datas$nasc,
    year          = format(datas$time_nasc, "%Y"),
    time          = datas$time_nasc,
    lat           = datas$lat_nasc,
    lon           = datas$lon_nasc,
    fod           = datas$fod,
    ftle          = datas$ftle,
    Chla_total    = datas$Chla_total,
    Chla_totpig   = datas$Chla_totpig,
    Per_totpig    = datas$Per_totpig,
    But_totpig    = datas$But_totpig,
    Fuco_totpig   = datas$Fuco_totpig,
    Hex_totpig    = datas$Hex_totpig,
    Allo_totpig   = datas$Allo_totpig,
    Zea_totpig    = datas$Zea_totpig,
    Chlb_totpig   = datas$Chlb_totpig,
    DvChla_totpig = datas$DvChla_totpig
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

# =====================================================================
# GRILLE DE PREDICTION -- fichier unique toutes-dates (structure
# corrigée, voir data_generation/generate_ds_ftle_pig_fod_all_dates.R) :
#   list(date[n_date], lon[1080], lat[720],
#        ftle[n_date,1080,720],
#        pig = list(Chla,...,DvChla (bruts) + chla_total, chla_totpig,
#                   ..., dvchla_totpig (normalisés, noms en minuscules),
#                   chacun [n_date,1080,720]),
#        fod[n_date,1080,720])
# `ftle`, `pig` ET `fod` ont maintenant TOUS le même ordre de dimensions
# -- plus de cas particulier à gérer pour fod.
# =====================================================================

load_grid_all_dates <- function(path_grid = PATH_GRID_ALL_DATES) {
  readRDS(path_grid)
}

get_date_index <- function(day_ds, target_date) {
  idx <- which(day_ds$date == as.Date(target_date))
  if (length(idx) == 0) {
    stop("Date ", target_date, " introuvable dans le fichier grille (dates disponibles : ",
         format(min(day_ds$date)), " a ", format(max(day_ds$date)), ").")
  }
  idx[1]
}

# ---- Covariables normalisées (COVARIATES_NUM), pour RF ranger / XGB / rfsrc --
# Utilise directement les champs déjà normalisés -- aucune formule à
# deviner (contrairement à l'ancienne grille mono-date incomplète).
extract_grid_for_date <- function(day_ds, date_idx, fod_levels) {
  grid_points <- expand.grid(lon = day_ds$lon, lat = day_ds$lat)

  grid_points$ftle <- as.vector(day_ds$ftle[date_idx, , ])

  for (src_name in names(MULTIDATE_PIG_NAME_MAP)) {
    target_name <- MULTIDATE_PIG_NAME_MAP[[src_name]]
    grid_points[[target_name]] <- as.vector(day_ds$pig[[src_name]][date_idx, , ])
  }

  for (v in COVARIATES_NUM) {
    bad <- !is.finite(grid_points[[v]])
    grid_points[[v]][bad] <- NA
  }

  grid_points$fod <- factor(formatC(as.vector(day_ds$fod[date_idx, , ]), width = 2),
                             levels = fod_levels)

  list(grid = grid_points, date = day_ds$date[date_idx])
}

