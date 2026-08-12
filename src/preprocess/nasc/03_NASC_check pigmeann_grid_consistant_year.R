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

# library
library(ncdf4)

folder_path <- ""

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
    lat = lat,
    lon = lon
  )
}

reference_year <- as.character(years[1])

for (y in years[-1]) {
  
  y <- as.character(y)
  
  same_lat <- identical(
    grilles[[reference_year]]$lat,
    grilles[[y]]$lat
  )
  
  same_lon <- identical(
    grilles[[reference_year]]$lon,
    grilles[[y]]$lon
  )
  
  cat(
    reference_year, "vs", y,
    "\n  lat :", same_lat,
    "\n  lon :", same_lon,
    "\n"
  )
  
  # save reference year grid
  pig_grid <- list(
    lat = grilles[["2018"]]$lat,
    lon = grilles[["2018"]]$lon
  )
  
  saveRDS(
    pig_grid,
    file = "pigmeann_grid.rds"
  )
  
}