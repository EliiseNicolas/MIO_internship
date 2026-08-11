# librairies
library(ncdf4)
library(dplyr)
library(ggplot2)
# clear var
rm(list = ls())

# hyperparameters
path2018 <- "/run/media/mmolinet/KER22/données elise/raw/acoustic/LOCEAN_SOOP-BA_A_20180105T121559Z_MARIONDUFRESNE_FV02_EchointegrationAcoustic-18-38-70-120-200_END-20180201T103636Z_C-20260522T153853Z.nc"
path2021 <- "/run/media/mmolinet/KER22/données elise/raw/acoustic/LOCEAN_SOOP-BA_A_20210122T143044Z_MARIONDUFRESNE_FV02_EchointegrationAcoustic-18-38-70-120-200_END-20210307T041919Z_C-20260609T160140Z.nc"
path2023 <- "/run/media/mmolinet/KER22/données elise/raw/acoustic/LOCEAN_SOOP-BA_A_20230123T103153Z_MARIONDUFRESNE_FV02_EchointegrationAcoustic-18-38-70-120-200_END-20230227T021804Z_C-20260728T105027Z.nc"

freq <- 200 # kHz

# ---------------------------------------------- Diagnostic NA

diagnostic_na <- function(path, freq, year){
  par(mfrow = c(2,2), oma = c(0,0,3,0))
  
  # Ouverture du fichier
  ds <- nc_open(path)
  
  Sv <- ncvar_get(ds, "Sv")
  idx_freq <- which(ncvar_get(ds, "instrument_frequency") == freq)
  Sv <- Sv[,,idx_freq]
  time <- ncvar_get(ds, "time")
  time <- as.POSIXct(
    time * 86400,
    origin = "1950-01-01",
    tz = "UTC"
  )
  depth <- ds$dim$depth$vals
  print(depth)
  nc_close(ds)
  
  # nombre total NA
  n_na <- sum(is.na(Sv))
  pct_na <- 100 * n_na / length(Sv)
  
  cat("Total NA :", n_na, "\n")
  cat("Percentage :", round(pct_na, 2), "%\n")
  
  # nombre de NA par profondeur
  pct_depth <- 100 * rowMeans(is.na(Sv))
  
  plot(
    x = depth,
    pct_depth,
    type = "l",
    xlab = "Depth (m)",
    ylab = "% NA", 
    main = "Nombre de NA par profondeur",
    cex.main = 0.8
  )
  
  mtext(paste0("Diagnostic NA ", year, " transect dataset"), outer = TRUE, cex = 1.5)
  # nombre de NA par profil
  print(depth)
  idx <- which(depth>25 & depth <165)
  Sv_crop <- Sv[idx,]
  na_time <- colSums(is.na(Sv_crop))
  
  plot(
    x = time,
    na_time,
    type = "l",
    xlab = "Time",
    ylab = "Number of NA", 
    main = "Nombre de NA par profil (cropped 25-165m)",
    cex.main = 0.8
  )
  
  # Profondeurs totalement manquantes
  print(c("Profondeurs totalement manquantes :  ", depth[which(rowSums(is.na(Sv)) == ncol(Sv))]))
  
  # nombre de profils totalement NA
  print(c("Nombre de profils totalement NA : ", which(colSums(is.na(Sv)) == nrow(Sv))))
  
  # image des NA
  na_matrix <- is.na(Sv)
  
  image(
    x = time,
    y = depth,
    z = t(na_matrix),
    col = c("black", "white"),
    xlab = "Time",
    ylab = "Depth (m)",
    main = "Location of missing values",
    cex.main = 0.8,
    ylim = rev(range(depth))
  )
  
  par(mfrow = c(1,1))
}

dev.off()
diagnostic_na(path2018, freq, 2018) # cut depth 165m Un profil avec 3 NA, le reste -> 0 NA
diagnostic_na(path2021, freq, 2021) # cut depth 165m Un profil avec 3 NA, le reste -> 0 NA
dev.off()
diagnostic_na(path2022, freq, 2022)# cut depth 165m Un profil avec 1 NA, le reste -> 0 NA
dev.off()
diagnostic_na(path2023, freq, 2023)# cut depth 165m Un profil avec 1 NA, le reste -> 0 NA


# ------------------------------------------------ NASC function
compute_nasc <- function(path, freq){
  ds <- nc_open(path)
  lat <- ncvar_get(ds, "latitude")
  lon <- ncvar_get(ds, "longitude")
  day <- ncvar_get(ds, "day")
  Sv <- ncvar_get(ds, "Sv") #(depth, time ("days since 1950-01-01 00:00:00 UTC"), channel)
  time <- ncvar_get(ds, "time") 
  time <- as.POSIXct(
    time * 86400,
    origin = "1950-01-01",
    tz = "UTC"
  )
  idx_freq <- which(ncvar_get(ds, "instrument_frequency")==freq)
  depth <- ds$dim$depth$vals
  idx <- which(depth>25 & depth <200)
  Sv <- Sv[idx,,idx_freq] # (depth, time)
  
  dim(Sv)
  nc_close(ds)
  
  # Compute NASC
  sv <- 10^(Sv/10) # linear sv
  int <- colSums(sv, na.rm = TRUE)
  sa <- int * 10 # car on a un step de profondeur de 10m 
  NASC <- 4 * pi * 1852**2 * sa # integration sur la profondeur
  print(length(NASC))
  
  print(dim(lon))
  
  return(data.frame(
    time = time,
    lat = lat,
    lon = lon,
    day = day,
    NASC = NASC
  ))
}

# Computation NASC
nasc2021 <- compute_nasc(path2021, freq)
nasc2022 <- compute_nasc(path2022, freq)
nasc2023 <- compute_nasc(path2023, freq)

# ---------------------------------------------- Diagnostic mean and error NASC
diagnostic_boxplot <- function(nasc, year, diurnal_period = "") {
  
  nasc$date <- as.Date(nasc$time)
  
  day_list <- c("night", "sunrise", "day", "sunset")
  
  # Filtrage selon la période
  if (diurnal_period != "") {
    if (!diurnal_period %in% day_list) {
      stop("Période inconnue : choisir night, sunrise, day ou sunset")
    }
    
    nasc <- nasc[nasc$day == which(day_list == diurnal_period), ]
  }
  
  # Boxplot NASC journalier
  ggplot(nasc, aes(x = factor(date), y = NASC)) +
    geom_boxplot(fill = "lightblue", outlier.colour = "red") +
    theme_bw() +
    labs(
      title = paste("Distribution de NASC par jour -", 
                    diurnal_period, "-", year),
      x = "Date",
      y = "NASC"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

diagnostic_boxplot(nasc2021, 2021, "night")
diagnostic_boxplot(nasc2022, 2022, "night")
diagnostic_boxplot(nasc2023, 2023, "night")


diagnostic_diurnal_variance <- function(nasc, year) {
  
  nasc$date <- as.Date(nasc$time)
  
  day_list <- c("all", "night", "sunrise", "day", "sunset")
  
  results <- list()
  
  for (period in day_list) {
    
    # Filtrage selon la période
    if (period == "all") {
      data_period <- nasc
    } else {
      period_index <- which(c("night", "sunrise", "day", "sunset") == period)
      data_period <- nasc[nasc$day == period_index, ]
    }
    
    # Statistiques journalières
    daily_stats <- data_period %>%
      group_by(date) %>%
      summarise(
        NASC_mean_day = mean(NASC, na.rm = TRUE),
        NASC_var_day = var(NASC, na.rm = TRUE),
        n = n(),
        .groups = "drop"
      ) %>%
      mutate(period = period)
    
    results[[period]] <- daily_stats
  }
  
  daily_all <- bind_rows(results)
  
  
  # Résumé annuel des statistiques journalières
  summary <- daily_all %>%
    group_by(period) %>%
    summarise(
      mean_daily_NASC = mean(NASC_mean_day, na.rm = TRUE),
      variance_daily_mean = var(NASC_mean_day, na.rm = TRUE),
      mean_daily_variance = mean(NASC_var_day, na.rm = TRUE),
      mean_n = mean(n),
      .groups = "drop"
    )
  
  print(summary)
  
  
  # Graphique : variabilité des moyennes journalières
  ggplot(daily_all, aes(x = period, y = NASC_var_day)) +
    geom_boxplot(fill = "lightblue") +
    theme_bw() +
    # ylim(0, 50000) +
    labs(
      title = paste("variabilité intra-journalière du NASC. -", year),
      x = "Période",
      y = "Variance journalière NASC"
    )
}
diagnostic_diurnal_variance(nasc2021, 2021)
diagnostic_diurnal_variance(nasc2022, 2022)
diagnostic_diurnal_variance(nasc2023, 2023)

# ---------------------------------------------- Concatenation
str(nasc2021)
ds_nasc <- bind_rows(
  data.frame(
    time = nasc2021$time,
    lat = nasc2021$lat,
    lon = nasc2021$lon,
    NASC = nasc2021$NASC,
    day = nasc2021$day,
    year = 2021
  ),
  data.frame(
    time = nasc2022$time,
    lat = nasc2022$lat,
    lon = nasc2022$lon,
    NASC = nasc2022$NASC,
    day = nasc2022$day,
    year = 2022
  ),
  data.frame(
    time = nasc2023$time,
    lat = nasc2023$lat,
    lon = nasc2023$lon,
    NASC = nasc2023$NASC,
    day = nasc2023$day,
    year = 2023
  )
)
str(nasc2023)

ds_nasc_day <- ds_nasc %>%
  filter(day == 3)

ds_nasc_night <- ds_nasc %>%
  filter(day == 1)

saveRDS(ds_nasc, "/run/media/mmolinet/KER22/données elise/elisou_ta_stagiaire_pref/prepross/NASC/ds_nasc_2021_2022_2023.rds")
saveRDS(ds_nasc_night, "/run/media/mmolinet/KER22/données elise/elisou_ta_stagiaire_pref/prepross/NASC/ds_nasc_2021_2022_2023_night.rds")
saveRDS(ds_nasc_day, "/run/media/mmolinet/KER22/données elise/elisou_ta_stagiaire_pref/prepross/NASC/ds_nasc_2021_2022_2023_day.rds")


saved_ds <- readRDS("/run/media/mmolinet/KER22/données elise/elisou_ta_stagiaire_pref/prepross/NASC/ds_nasc_2021_2022_2023.rds")

# verifs 

# verfi toutes données
n2021 <- length(nasc2021$nasc)
n2022 <- length(nasc2022$nasc)
n2023 <- length(nasc2023$nasc)
n_expected <- n2021 + n2022 + n2023
nrow(saved_ds) == n_expected # ok

# verif que lat non modifiées
lat_global_2021 <- saved_ds$lat[saved_ds$year == 2021]
lat_nasc_2021 <- as.vector(nasc2021$lat)
lat_global_2021 <- as.vector(lat_global_2021)
all.equal(sort(lat_nasc_2021),sort(lat_global_2021)) # OK