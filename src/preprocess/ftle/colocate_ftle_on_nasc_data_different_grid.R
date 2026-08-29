# Description 
#
# from NASC dataset (computed with 05_compute_nasc_mean_profile_pig_grid.R),
# NASC dataset of shape (n_ESU, 4) with 4 columns time, lon, lat, NASC
#
# So from this dataset, we want to associate FTLE data that matches time, lat, long
# and create a new ds : nasc_ftle_pig_grid_day_200kHz.rds containing (n_ESU, 5)
# with 5 columns time, lon, lat, NASC, FTLE

# libraries
library(dplyr)
library(tidyr)

# global variables
freqs <- c(18, 38, 70, 120, 200)
lat_res <- c(200, 500, 700, 1000)
lon_res <- c(200, 500, 1000, 1500)

for (g in seq_along(lat_res)) {
  lon_res_i <- lon_res[g]
  lat_res_i <- lat_res[g]
  
  ftle_path <- paste0("F:/data_elise/ftle/grids_custom/ftle_pigmeann_grid_lon", lon_res_i, "_lat", lat_res_i, ".rds")
  ftle_grid <- readRDS(ftle_path)   # list: date, lon, lat, ftle[date,lon,lat], n_valid[date,lon,lat]
  str(ftle_grid)
  
  ftle_dates_str <- format(as.Date(ftle_grid$date), "%Y-%m-%d")
  
  grid_label <- paste0("lon", lon_res_i, "_lat", lat_res_i)
  
  for (freq in freqs){
    # ---------------------------------------------------------
    # NASC
    # ---------------------------------------------------------
    nasc_path <- paste0("F:/data_elise/NASC/NASC_pig_mean/NASC_mean_Sv_pig_grid_lon", lon_res_i, "_lat", lat_res_i,"_2018_2022_2021_2023_", freq, "kHz.rds")
    nasc_ds <- readRDS(nasc_path)
    print(head(nasc_ds$time))
    
    dates_all_str <- format(as.Date(nasc_ds$time), "%Y-%m-%d")
    dates_unique_str <- unique(dates_all_str)
    
    print(length(dates_unique_str))
    
    # ---------------------------------------------------------
    # EXTRACTION FTLE
    # ---------------------------------------------------------
    
    ftle_out <- data.frame(
      time     = nasc_ds$time,
      lat_sv   = nasc_ds$lat,
      lon_sv   = nasc_ds$lon,
      lat_ftle = NA_real_,
      lon_ftle = NA_real_,
      ftle     = NA_real_,
      n_valid  = NA_real_
    )
    
    for (date_i in dates_unique_str){
      
      print(date_i)
      
      # =====================================================
      # INDICES ESU À CETTE DATE
      # =====================================================
      
      ind <- which(dates_all_str == date_i)
      
      if (length(ind) == 0) next
      
      # =====================================================
      # INDICE DU JOUR DANS LA GRILLE FTLE
      # =====================================================
      
      idx_date <- which(ftle_dates_str == date_i)
      
      if (length(idx_date) == 0){
        warning(paste("Aucune date FTLE trouvée pour", date_i))
        next
      }
      idx_date <- idx_date[1]  # sécurité si doublon
      
      # =====================================================
      # COORDONNÉES FTLE
      # =====================================================
      
      lat_ftle <- ftle_grid$lat
      lon_ftle <- ftle_grid$lon
      
      # =====================================================
      # TROUVER LES PIXELS FTLE LES PLUS PROCHES
      # =====================================================
      
      idx_lon <- sapply(nasc_ds$lon[ind], function(x) which.min(abs(lon_ftle - x)))
      idx_lat <- sapply(nasc_ds$lat[ind], function(x) which.min(abs(lat_ftle - x)))
      
      # =====================================================
      # ENREGISTRER LES COORDONNÉES FTLE ASSOCIÉES
      # =====================================================
      
      ftle_out$lat_ftle[ind] <- lat_ftle[idx_lat]
      ftle_out$lon_ftle[ind] <- lon_ftle[idx_lon]
      
      # =====================================================
      # EXTRACTION AU VOISINAGE DES POINTS
      # =====================================================
      
      for (j in seq_along(ind)){
        i <- ind[j]
        ftle_out$ftle[i]    <- ftle_grid$ftle[idx_date, idx_lon[j], idx_lat[j]]
        ftle_out$n_valid[i] <- ftle_grid$n_valid[idx_date, idx_lon[j], idx_lat[j]]
      }
      
      rm(idx_lon, idx_lat, lat_ftle, lon_ftle, idx_date)
    }
    
    dlat_ftle <- nasc_ds$lat - ftle_out$lat_ftle
    dlon_ftle <- nasc_ds$lon - ftle_out$lon_ftle
    summary(dlat_ftle)
    summary(dlon_ftle) # OK
    str(ftle_out)
    
    saveRDS(
      ftle_out,
      file = paste0("F:/data_elise/ftle/ftle_colocated_transect/ftle_colocated_NASC_mean_grid_pig/ftle_colocated_with_NASC_mean_grid_2018_2021_2022_2023_", freq, "kHz_lon", lon_res_i, "_lat", lat_res_i, ".rds")
    )
  }
}
print(ftle_dates_str)
ftle_raw <- readRDS("F:/data_elise/ftle/ftle_2018_2021_2022_2023_cropped.rds")
print(ftle_raw$date)
# # diagnostic des NAs
# rm(list=ls())
# freq <- 18
# ftle_out <- readRDS(paste0("F:/data_elise/ftle/ftle_colocated_transect/ftle_colocated_NASC_per_esu/ftle_colocated_with_NASC_per_esu_2018_2021_2022_2023_", freq, "kHz_lon", lon_res_i, "_lat", lat_res_i, ".rds"))
# print(sum(is.na(ftle_out)))
# idx_na <- which(is.na(ftle_out$ftle))
# time_na <- ftle_out$time[idx_na]
# print(time_na)
# print(unique(as.Date(time_na)))
# unique(is.na(ftle_out$time))