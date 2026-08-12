# Description

# 01_NASC_filter_sv
      # From obsaustral ncdf transect data 2018, 2021, 2023 : 
      #   1) filter by channel, keep only 200kHz data
      #   2) cut data on depth dimension to remove most NA
      #   3) filter lon >40°E and lat <-30°N
      #   4) filter day/night 
      # => save intermediary rds file : 2018_day.rds, 2018_night.rds, 2021_day.rds,...

# 02_NASC_concat_sv_diurnal_period
      #   5) Concatenate each year of same period (i.e. 2018_day, 2021_day, 2023_day)
      # => save intermediary rds file : 2018_2021_2022_2023_day.rds, 2018_2021_2022_2023_night.rds

# 03_NASC_check_pigmeann_grid
      #   6) get pigmeann grid (check if pigmeann grid consistent on every year of data)
      #   7) mean profile by grid (for day and night data)
      # => save intermediary rds file : mean_pig_grid_2018_2021_2023_day.rds, ...

#   8) compute NASC for each mean profile 
# => save intermediary rds file : nasc_mean_pig_grid_2018_2021_2023_day.rds,..

# rm(list=ls())
# library
library(ncdf4)

folder_path <- "/mnt/KER22/data_elise/raw/PIGMeANN/daily"

years <- c(2018, 2021, 2022, 2023)

# Tous les fichiers .nc
fichiers <- list.files(
  folder_path,
  pattern = "\\.nc$",
  recursive = TRUE,
  full.names = TRUE
)

# Stocker les grilles
grilles <- list()

for (y in years) {
  
  # Fichiers correspondant à l'année
  files_y <- fichiers[grepl(as.character(y), basename(fichiers))]
  
  print(paste("Année :", y))
  print(files_y)
  
  # Vérifier qu'on a bien un fichier
  if (length(files_y) == 0) {
    warning(paste("Aucun fichier trouvé pour", y))
    next
  }
  
  # Ouvrir le premier fichier
  nc <- nc_open(files_y[1])
  
  lat <- nc$dim$lat$vals
  lon <- nc$dim$lon$vals
  
  nc_close(nc)
  
  # Stocker
  grilles[[as.character(y)]] <- list(
    lat = lat[lat >= -60 & lat <= -30],
    lon = lon[lon > 44.97 & lon <= 90]
  )
}

reference_year <- as.character(years[1])

for (y in years[-1]) {
  
  y <- as.character(y)
  cat(
    "\n", reference_year, "vs", y, "\n",
    "lat :", length(grilles[[reference_year]]$lat), "vs", length(grilles[[y]]$lat), "\n",
    "lon :", length(grilles[[reference_year]]$lon), "vs", length(grilles[[y]]$lon), "\n"
  )
  
  
  
  lat_diff <- grilles[[y]]$lat - grilles[[reference_year]]$lat
  lon_diff <- grilles[[y]]$lon - grilles[[reference_year]]$lon
  
  
  cat(
    "LAT : max =", max(abs(lat_diff), na.rm = TRUE),
    " | mean =", mean(abs(lat_diff), na.rm = TRUE),
    "\n"
  )
  
  cat(
    "LON : max =", max(abs(lon_diff), na.rm = TRUE),
    " | mean =", mean(abs(lon_diff), na.rm = TRUE),
    "\n"
  )
  
  same_lat <- identical(
    round(grilles[[reference_year]]$lat, 4),
    round(grilles[[y]]$lat, 4)
  )
  
  same_lon <- identical(
    round(grilles[[reference_year]]$lon, 4),
    round(grilles[[y]]$lon, 4)
  )
  
  cat(
    reference_year, "vs", y,
    "\n  lat :", same_lat,
    "\n  lon :", same_lon,
    "\n"
  )
  
}

# save reference year grid
pig_grid <- list(
  lat = grilles[["2018"]]$lat,
  lon = grilles[["2018"]]$lon
)
saveRDS(
  pig_grid,
  file = "/home/elise/Documents/stage_MIO/pt_III/data_preprocessed/NASC/transect_2018_2022_2023/pigmeann_grid.rds"
)


# Test OK => same grid if we filter enough to keep 
# lat = lat[lat >= -60 & lat <= -30],
# lon = lon[lon > 44.97 & lon <= 90]