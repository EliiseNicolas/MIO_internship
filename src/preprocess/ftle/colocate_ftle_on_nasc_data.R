# Description 

# from NASC dataset (computed with 05_compute_nasc_mean_profile_pig_grid.R),
# NASC dataset of shape (n_ESU, 4) with 4 columns time, lon, lat, NASC
# 
# So from this dataset, we want to associate FTLE data that matches time, lat, long
# and create a new ds : nasc_ftle_pig_grid_day_200kHz.rds containing (n_ESU, 5) 
# with 5 columns time, lon, lat, NASC, FTLE

# libraries
library(ncdf4)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# Global variables 
nasc_path <- "F:/data_elise/NASC/120kHz/NASC_mean_pig_grid_by_year_2018_2021_2023_night_120kHz.rds"
folder_path <- "F:/data_elise/ftle"

# ---------------------------------------------------------
# NASC
# ---------------------------------------------------------
nasc_ds <- readRDS(nasc_path)
print(head(nasc_ds$time))
dates_all_str <- format(as.Date(nasc_ds$time), "%Y-%m-%d")
dates_unique_str <- (unique(dates_all_str))

print(length(dates_unique_str))

# ---------------------------------------------------------
# FICHIERS FTLE
# ---------------------------------------------------------

files <- list.files(
  path = folder_path,
  pattern = "\\.nc$",
  full.names = TRUE
)
print(files)

# On garde uniquement les fichiers correspondant aux dates
files <- files[
  grepl(
    paste(dates_unique_str, collapse = "|"),
    basename(files)
  )
]

print(length(files)) # 68 sur 76 pour 120kHz
print(basename(files))


# ---------------------------------------------------------
# EXTRACTION FTLE
# ---------------------------------------------------------

ftle <- data.frame(
  time = nasc_ds$time,
  lat_sv = nasc_ds$lat,
  lon_sv = nasc_ds$lon,
  lat_ftle = NA_real_,
  lon_ftle = NA_real_, 
  ftle = NA_real_
)

p <- 3
for(date_i in dates_unique_str){
  
  print(date_i)
  
  # =====================================================
  # INDICES ESU À CETTE DATE
  # =====================================================
  
  ind <- which(dates_all_str == date_i)
  
  if(length(ind) == 0)
    next
  
  
  # =====================================================
  # OUVRIR LE FICHIER DU JOUR
  # =====================================================
  
  pattern <- paste0("map_", date_i, ".*\\.nc$")
  
  filename <- list.files(
    folder_path,
    pattern = pattern,
    full.names = TRUE
  )
  
  if(length(filename) == 0)
    next
  
  ds <- nc_open(filename[1])
  
  
  # =====================================================
  # COORDONNÉES FTLE
  # =====================================================
  
  lat_ftle <- ds$dim$lat$vals
  lon_ftle <- ds$dim$lon$vals
  
  nlat <- length(lat_ftle)
  nlon <- length(lon_ftle)
  
  
  # =====================================================
  # TROUVER LES PIXELS FTLE LES PLUS PROCHES
  # =====================================================
  
  idx_lon <- sapply(
    nasc_ds$lon[ind],
    function(x) which.min(abs(lon_ftle - x))
  )
  
  idx_lat <- sapply(
    nasc_ds$lat[ind],
    function(x) which.min(abs(lat_ftle - x))
  )
  
  
  # =====================================================
  # ENREGISTRER LES COORDONNÉES FTLE ASSOCIÉES
  # =====================================================
  
  ftle$lat_ftle[ind] <- lat_ftle[idx_lat]
  ftle$lon_ftle[ind] <- lon_ftle[idx_lon]
  
  
  # =====================================================
  # FENÊTRE LAT/LON À LIRE
  # =====================================================
  
  r <- (p - 1) / 2
  
  lon_start <- max(1, min(idx_lon) - r)
  lon_end   <- min(nlon, max(idx_lon) + r)
  
  lat_start <- max(1, min(idx_lat) - r)
  lat_end   <- min(nlat, max(idx_lat) + r)
  
  lon_count <- lon_end - lon_start + 1
  lat_count <- lat_end - lat_start + 1
  
  
  # =====================================================
  # LIRE UNIQUEMENT LA FENÊTRE FTLE
  # =====================================================
  
  ftle_data <- ncvar_get(
    ds,
    "FTLE",
    start = c(lon_start, lat_start, 1),
    count = c(lon_count, lat_count, 1)
  )
  
  
  # =====================================================
  # EXTRACTION AU VOISINAGE DES POINTS
  # =====================================================
  
  for(j in seq_along(ind)){
    
    i <- ind[j]
    
    # coordonnées de la station dans la fenêtre lue
    ilon <- idx_lon[j] - lon_start + 1
    ilat <- idx_lat[j] - lat_start + 1
    
    r <- (p - 1) / 2
    
    lon_win <- max(
      1,
      ilon - r
    ):min(
      nrow(ftle_data),
      ilon + r
    )
    
    lat_win <- max(
      1,
      ilat - r
    ):min(
      ncol(ftle_data),
      ilat + r
    )
    
    ftle$ftle[i] <- mean(
      ftle_data[lon_win, lat_win],
      na.rm = TRUE
    )
  }
  
  
  # =====================================================
  # FERMER / LIBÉRER
  # =====================================================
  
  nc_close(ds)
  
  rm(ds, ftle_data)
  gc()
}

str(ftle)

saveRDS(
  ftle,
  file = paste0("F:/ftle_colocated_transect/120kHz/ftle",
                p*p,
                "_1d_",
                "_colocated_with_NASC_transect_2018_2021_2023_120kHz_night.rds"
  )
)
