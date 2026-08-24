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
folder_path <- "F:/data_elise/ftle/ftle_raw"

# ---------------------------------------------------------
# FICHIERS FTLE
# ---------------------------------------------------------

files <- list.files(
  path = folder_path,
  pattern = "\\.nc$",
  full.names = TRUE
)
print(files)


for (freq in c(18, 38, 70, 120, 200)){
  # ---------------------------------------------------------
  # NASC
  # ---------------------------------------------------------
  # nasc_path <- paste0("F:/data_elise/NASC/NASC_pig_mean/mean_Sv_pig_grid_by_date_2018_2021_2023_", freq, "kHz.rds") # mean NASC per pig grid
  nasc_path <- paste0("F:/data_elise/NASC/NASC_all_ESU/NASC_per_ESU_2018_2021_2023_", freq, "kHz.rds") # NASC per ESU
  nasc_ds <- readRDS(nasc_path)
  print(head(nasc_ds$time))
  dates_all_str <- format(as.Date(nasc_ds$time), "%Y-%m-%d")
  dates_unique_str <- (unique(dates_all_str))
  
  print(length(dates_unique_str))
  
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
  
  
  for(date_i in dates_unique_str){
    
    print(date_i)
    
    # =====================================================
    # INDICES ESU À CETTE DATE
    # =====================================================
    
    ind <- which(dates_all_str == date_i)
    
    if(length(ind) == 0) next
    
    
    # =====================================================
    # OUVRIR LE FICHIER DU JOUR
    # =====================================================
    
    pattern <- paste0("map_", date_i, ".*\\.nc$")
    
    filename <- list.files(folder_path, pattern = pattern, full.names = TRUE)
    
    if(length(filename) == 0){
      warning(paste("Aucun fichier FTLE trouvé pour", date_i))
      next
    }
    
    ds <- nc_open(filename[1])
    
    
    # =====================================================
    # COORDONNÉES FTLE
    # =====================================================
    
    lat_ftle <- ds$dim$lat$vals
    lon_ftle <- ds$dim$lon$vals
    
    # =====================================================
    # TROUVER LES PIXELS FTLE LES PLUS PROCHES
    # =====================================================
    
    idx_lon <- sapply(nasc_ds$lon[ind], function(x) which.min(abs(lon_ftle - x)))
    idx_lat <- sapply(nasc_ds$lat[ind], function(x) which.min(abs(lat_ftle - x)))
    
    
    # =====================================================
    # ENREGISTRER LES COORDONNÉES FTLE ASSOCIÉES
    # =====================================================
    
    ftle$lat_ftle[ind] <- lat_ftle[idx_lat]
    ftle$lon_ftle[ind] <- lon_ftle[idx_lon]
    
    
    # =====================================================
    # LIRE LES FTLE
    # =====================================================
    
    ftle_data <- ncvar_get(ds, "FTLE")
    
    
    # =====================================================
    # EXTRACTION AU VOISINAGE DES POINTS
    # =====================================================
    
    for(j in seq_along(ind)){
      
      i <- ind[j]
      
      ftle$ftle[i] <- ftle_data[idx_lon[j], idx_lat[j]]
    }
    
    
    # =====================================================
    # FERMER / LIBÉRER
    # =====================================================
    
    nc_close(ds)
    
    rm(ds, ftle_data, idx_lon, idx_lat, lat_ftle, lon_ftle)
    gc()
  }
  
  dlat_ftle <- nasc_ds$lat - ftle$lat_ftle
  dlon_ftle <- nasc_ds$lon - ftle$lon_ftle
  summary(dlat_ftle)
  summary(dlon_ftle) # OK
  str(ftle)
  
  saveRDS(
    ftle,
    file = paste0("F:/data_elise/ftle/ftle_colocated_transect/ftle_colocated_NASC_per_esu/ftle_colocated_with_NASC_per_esu_2018_2021_2023_", freq, "kHz.rds"
    )
  )
}
