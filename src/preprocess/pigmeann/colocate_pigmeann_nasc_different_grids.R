library(dplyr)
library(tidyr)

rm(list = ls())

# Global Variables

list_pigs <- c("Chla", "Per", "But", "Fuco", "Hex", "Allo", "Zea", "Chlb", "DvChla")

freqs <- c(18, 38, 70, 120, 200)
lat_res <- c(200, 500, 700, 1000)
lon_res <- c(200, 500, 1000, 1500)

for (g in seq_along(lat_res)) {
  lon_res_i <- lon_res[g]
  lat_res_i <- lat_res[g]
  
  pig_path <- paste0("F:/data_elise/pigmeann/grids_custom/pigments_pigmeann_grid_lon", lon_res_i, "_lat", lat_res_i, ".rds")
  pigs <- readRDS(pig_path)   # list: date, lon, lat, c_cond_<pig>[date,lon,lat], n_valid_<pig>[date,lon,lat]
  str(pigs)
  
  pigs_dates_str <- format(as.Date(pigs$date), "%Y-%m-%d")
  
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
    # Initialisation dataframe
    # ---------------------------------------------------------
    
    pig_cond <- data.frame(
      time    = nasc_ds$time,
      lat_sv  = nasc_ds$lat,
      lon_sv  = nasc_ds$lon,
      lat_pig = NA_real_,
      lon_pig = NA_real_
    )
    
    # Ajouter les colonnes pigments (concentration + n_valid)
    pig_cond[list_pigs] <- NA_real_
    pig_cond[paste0("n_valid_", list_pigs)] <- NA_real_
    
    # ---------------------------------------------------------
    # EXTRACTION PIGMENTS (point le plus proche)
    # ---------------------------------------------------------
    
    for (date_i in dates_unique_str){
      
      print(date_i)
      
      # =====================================================
      # INDICES ESU À CETTE DATE
      # =====================================================
      
      ind <- which(dates_all_str == date_i)
      
      if (length(ind) == 0) next
      
      # =====================================================
      # INDICE DU JOUR DANS LA GRILLE PIGMENTS
      # =====================================================
      
      idx_date <- which(pigs_dates_str == date_i)
      
      if (length(idx_date) == 0){
        warning(paste("Aucune date PIGMeANN trouvée pour", date_i))
        next
      }
      idx_date <- idx_date[1]  # sécurité si doublon
      
      # =====================================================
      # COORDONNÉES PIGMENTS
      # =====================================================
      
      lat_pig <- pigs$lat
      lon_pig <- pigs$lon
      
      # =====================================================
      # TROUVER LES PIXELS PIGMeANN LES PLUS PROCHES
      # =====================================================
      
      idx_lon <- sapply(nasc_ds$lon[ind], function(x) which.min(abs(lon_pig - x)))
      idx_lat <- sapply(nasc_ds$lat[ind], function(x) which.min(abs(lat_pig - x)))
      
      # =====================================================
      # ENREGISTRER LES COORDONNÉES PIGMeANN ASSOCIÉES
      # =====================================================
      
      pig_cond$lat_pig[ind] <- lat_pig[idx_lat]
      pig_cond$lon_pig[ind] <- lon_pig[idx_lon]
      
      # =====================================================
      # EXTRACTION AU POINT LE PLUS PROCHE, POUR CHAQUE PIGMENT
      # =====================================================
      
      for (pig in list_pigs){
        
        c_cond_arr  <- pigs[[paste0("c_cond_", pig)]]
        n_valid_arr <- pigs[[paste0("n_valid_", pig)]]
        
        for (j in seq_along(ind)){
          i <- ind[j]
          pig_cond[i, pig]                          <- c_cond_arr[idx_date, idx_lon[j], idx_lat[j]]
          pig_cond[i, paste0("n_valid_", pig)]       <- n_valid_arr[idx_date, idx_lon[j], idx_lat[j]]
        }
      }
      
      rm(idx_lon, idx_lat, lat_pig, lon_pig, idx_date)
    }
    
    dlat_pig <- nasc_ds$lat - pig_cond$lat_pig
    dlon_pig <- nasc_ds$lon - pig_cond$lon_pig
    summary(dlat_pig)
    summary(dlon_pig) # OK
    str(pig_cond)
    
    # ---------------------------------------------------------
    # AJOUT DES RATIOS
    # ---------------------------------------------------------
    
    # On part des concentrations originales (point le plus proche)
    pig_cond_ratio <- pig_cond
    
    # =========================================================
    # SOMME DES CONCENTRATIONS DE TOUS LES PIGMENTS
    # =========================================================
    
    total_pig <- rowSums(
      pig_cond[list_pigs],
      na.rm = TRUE
    )
    
    # Si aucun pigment n'est disponible, on met NA
    total_pig[total_pig == 0] <- NA_real_
    
    # Ajouter la concentration totale
    pig_cond_ratio$total_pig <- total_pig
    
    # =========================================================
    # RATIOS PAR RAPPORT A LA ChlA
    # =========================================================
    
    for (pig in list_pigs) {
      pig_cond_ratio[[paste0(pig, "_Chla")]] <-
        pig_cond[[pig]] / pig_cond$Chla
    }
    
    # =========================================================
    # RATIOS PAR RAPPORT AU TOTAL DES PIGMENTS
    # =========================================================
    
    for (pig in list_pigs) {
      pig_cond_ratio[[paste0(pig, "_total")]] <-
        pig_cond[[pig]] / total_pig
    }
    
    # ---------------------------------------------------------
    # SAUVEGARDE
    # ---------------------------------------------------------
    
    saveRDS(
      pig_cond_ratio,
      file = paste0(
        "F:/data_elise/pigmeann/pigs_colocated_NASC_mean_pig_grid/",
        "NASC_mean_pig_grid_pig_conc_ratio_nearest_point_",
        "2018_2021_2022_2023_", freq, "kHz_", grid_label, ".rds"
      )
    )
  }
}
str(pig_cond_ratio)