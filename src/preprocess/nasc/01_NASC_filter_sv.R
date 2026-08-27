# Description
# From obsaustral ncdf transect data 2018, 2021, 2023 : 
#   1) filter by channel, keep only 200kHz data
#   2) cut data on depth dimension to remove most NA
#   3) filter lon >40°E and lat <-30°N
#   4) filter day/night 
# => save intermediary rds file : 2018_day.rds, 2018_night.rds, 2021_day.rds,...

#   5) Concatenate each year of same period (i.e. 2018_day, 2021_day, 2023_day)
# => save intermediary rds file : 2018_2021_2022_2023_day.rds, 2018_2021_2022_2023_night.rds

#   6) get pigmeann grid (check if pigmeann grid consistent on every year of data)
#   7) mean profile by grid (for day and night data)
# => save intermediary rds file : mean_pig_grid_2018_2021_2023_day.rds, ...

#   8) compute NASC for each mean profile 
# => save intermediary rds file : nasc_mean_pig_grid_2018_2021_2023_day.rds,..

# Packages
library(ncdf4)

rm(list=ls())

# Global variables

path2018 <- "G:/data_elise/raw/acooustic/LOCEAN_SOOP-BA_A_20180105T121559Z_MARIONDUFRESNE_FV02_EchointegrationAcoustic-18-38-70-120-200_END-20180201T103636Z_C-20260522T153853Z.nc"
path2021 <- "G:/data_elise/raw/acooustic/LOCEAN_SOOP-BA_A_20210122T143044Z_MARIONDUFRESNE_FV02_EchointegrationAcoustic-18-38-70-120-200_END-20210307T041919Z_C-20260609T160140Z.nc"
path2022 <- "G:/data_elise/raw/acooustic/LOCEAN_SOOP-BA_A_20220201T075726Z_MARIONDUFRESNE_FV02_EchointegrationAcoustic-18-38-70-120-200_END-20220305T041935Z_C-20260827T064134Z.nc"
path2023 <-"G:/data_elise/raw/acooustic/LOCEAN_SOOP-BA_A_20230123T103153Z_MARIONDUFRESNE_FV02_EchointegrationAcoustic-18-38-70-120-200_END-20230227T021804Z_C-20260728T105027Z.nc"
years <- c("2018", "2021", "2022", "2023")
paths <- c("2018"=path2018, "2022"=path2022, "2021"=path2021, "2023"=path2023)

freqs <- c(18, 38, 70, 120, 200)
outpath <- paste0("F:/data_elise/sv_cropped/sv_cropped_per_year")
dir.create(outpath, recursive = TRUE, showWarnings = FALSE) # creer le dir s'il n'existe pas deja


depth_min <- 25; 
max_depths <- c(700, 700, 500, 300, 150)
lat_min <- -60; lat_max <- -30
lon_min <- 45; lon_max <- 90

# ------------------------------------------part I - find depth max value 
# define depth max value relative to nans
diagnostic_na <- function(path, freq, year){
  
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
  nc_close(ds)
  
  # nombre total NA
  n_na <- sum(is.na(Sv))
  pct_na <- 100 * n_na / (nrow(Sv) * ncol(Sv))
  
  cat("Total NA :", n_na, "\n")
  cat("Percentage :", round(pct_na, 2), "%\n")
  
  # nombre de NA par profondeur
  pct_depth <- 100 * rowMeans(is.na(Sv))
  
  # nombre de NA par profil
  idx <- which(depth > depth_min & depth < depth_max)
  Sv_crop <- Sv[idx,]
  
  n_na_crop <- sum(is.na(Sv_crop))
  pct_na_crop <- 100 * n_na_crop / 
    (nrow(Sv_crop) * ncol(Sv_crop))
  
  cat("Total NA cropped :", n_na_crop, "\n")
  cat("Percentage cropped :", round(pct_na_crop, 2), "%\n")
  
  na_time <- colSums(is.na(Sv_crop))
  
  
  
  # Plot 1 : NA par profondeur
  # Ouverture du fichier PNG
  filename <- paste0(
    "C:/Users/mmolinet/elisou_ta_stagiaire_pref/MIO_internship_III/figures/diag_NA_sv_profiles/",
    freq, "kHz/diagnostic_NA_", year, "_", freq,
    "kHz_transect_all_dataset.png"
  )
  
  print(filename)
  png(
    filename = paste0(
      "C:/Users/mmolinet/elisou_ta_stagiaire_pref/MIO_internship_III/figures/diag_NA_sv_profiles/",
      freq, "kHz/diagnostic_NA_", year, "_", freq,
      "kHz_transect_all_dataset.png"
    ),
    width = 1600,
    height = 900,
    res = 150
  )

  # Configuration de la fenêtre graphique
  par(
    mfrow = c(1, 2),
    oma = c(0, 0, 3, 0)
  )

  plot(
    x = depth,
    y = pct_depth,
    type = "l",
    xlab = "Depth (m)",
    ylab = "% NA",
    main = "Nombre de NA par profondeur",
    cex.main = 0.8
  )

  # Plot 2 : NA par profil
  plot(
    x = time,
    y = na_time,
    type = "l",
    xlab = "Time",
    ylab = "Number of NA",
    main = paste0(
      "Nombre de NA par profil\n",
      "(cropped ", depth_min, "-", depth_max, " m)"
    ),
    cex.main = 0.8
  )

  # Titre général de la fenêtre
  mtext(
    paste0("Diagnostic NA ", year, " transect dataset ", freq, " kHz"),
    outer = TRUE,
    cex = 1.5
  )
  dev.off()
  # Profondeurs totalement manquantes
  print(
    c(
      "Profondeurs totalement manquantes :",
      depth[which(rowSums(is.na(Sv)) == ncol(Sv))]
    )
  )
  
  # profils totalement NA
  print(
    c(
      "Nombre de profils totalement NA :",
      length(which(colSums(is.na(Sv)) == nrow(Sv)))
    )
  )
}
year <- "2022"
for (i in 1:5){
  freq <- freqs[i]
  print(freq)
  depth_max <- max_depths[i]
  diagnostic_na(paths[year], freq, year) 
  # for (year in years){
  #   print(year)
  #   diagnostic_na(paths[year], freq, year) 
  # }
}


# cut depth 1000m for 18kHz, 300m for 120kHz and 165m for 200kHz
# 18kHz - après crop a 700m : 2018 4.03%, 2021 17.79%, 2023 4.14%
# 38kHz - après crop a 700m : 2018 2.94%, 2021 17.66%, 2023 2.07%
# 70kHz - après crop a 500m : 2018 3.28%, 2021 14.08%, 2023 1.51%
# 120kHz - après crop 300m : 2018 3.96% de NA, 2021 9.47% de NA, 2023 1.2% de NA
# 200kHz - après crop 150m : 2018 3.69%, 2021 3.83%, 2023 0.69%



# --------------------------------------------------------------- Part II - filtering
years <- "2022"
for (i in 1:5){
  freq <- freqs[i]
  print(freq)
  depth_max <- max_depths[i]
  for (year in years){
    print(year)
    path <- paths[year]
    ds <- nc_open(path)
    Sv <- ncvar_get(ds, "Sv") # (depth, time, frequency)
    print(dim(Sv))
    time <- ncvar_get(ds, "time") # (time, )
    lat <- ncvar_get(ds, "latitude") # (time,)
    lon <- ncvar_get(ds, "longitude") # (time,)
    depth <- ncvar_get(ds, "depth")
    day <- ncvar_get(ds, "day") # (time,)
    
    # filter by channel
    idx_freq <- which(ncvar_get(ds, "instrument_frequency")==freq)
    print(idx_freq)
    Sv <- Sv[,,idx_freq] # (depth, time)
    
    # crop data
    depth_idx <- which(depth>depth_min & depth<depth_max)
    time_idx <- which(lat>lat_min & lat<lat_max & lon>lon_min & lon<lon_max & (day==1 | day==3))
    
    lat <- lat[time_idx]
    lon <- lon[time_idx]
    depth <- depth[depth_idx]
    day <- day[time_idx]
    print(c(range(depth), range(lat), range(lon), unique(day)))
    
    print(dim(Sv))
    Sv <- Sv[depth_idx, time_idx]
    print(dim(Sv))
    time <- time[time_idx]
    
    # crop depth for 2022 ds
    if (year == "2022"){
      print("yes")
      n_depth <- length(depth)
      depth <- round(depth[2:n_depth])
      Sv <- Sv[2:n_depth,]
      print(c(dim(Sv), depth)) # 337 depth et commence a 27m
    }
    
    # Create new dataset
    time_posix <- as.POSIXct("1950-01-01", tz = "UTC") + time * 86400
    
    df <- list(
      profiles = t(Sv),
      lat = lat,
      lon = lon,
      depth = depth,
      time = time_posix, 
      day = day
    )
    
    nc_close(ds)
    
    # save
    saveRDS(df, paste0(outpath, "/Sv_", year, "_", freq, "kHz" , ".rds"))
  }
}





