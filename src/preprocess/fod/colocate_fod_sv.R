# Description

# Create .rds with following columns : time, lat, lon, cluster 1, cluster 2, ..., cluster n 
# The number of rows corresponds to the number of ESU collected during oceanographic 
# campaigns of interest (2021 2022 2023)

rm(list = ls())

# Librairies
library(ncdf4)

# Paths
freqs <- c(18, 38, 70, 120, 200)
lat_res <- c(200, 500, 700, 1000)
lon_res <- c(200, 500, 1000, 1500)

path_clusters_2018_2021_2022_2023 <- "F:/data_elise/fod_elise_2018_2021_2022_2023/mclust_results_2018_2021_2022_2023/transitions_upgraded/cluster_transition_map_renamed.rds"
lon_path <- "F:/data_elise/fod_elise_2018_2021_2022_2023/mclust_results_2018_2021_2022_2023/lon.rds"
lat_path <- "F:/data_elise/fod_elise_2018_2021_2022_2023/mclust_results_2018_2021_2022_2023/lat.rds"
time_path <- "F:/data_elise/fod_elise_2018_2021_2022_2023/mclust_results_2018_2021_2022_2023/time.rds"

# Open RDS
clusters <- readRDS(path_clusters_2018_2021_2022_2023)
lat_fod <- readRDS(lat_path)
lon_fod <- readRDS(lon_path)
time_fod <- readRDS(time_path)
date_fod <- as.Date(time_fod)

for (g in seq_along(lat_res)) {
  lon_res_i <- lon_res[g]
  lat_res_i <- lat_res[g]
  grid_label <- paste0("lon", lon_res_i, "_lat", lat_res_i)
  
  for (freq in freqs){
    print(freq)
    path_nasc <- paste0("F:/data_elise/NASC/NASC_pig_mean/NASC_mean_Sv_pig_grid_lon", lon_res_i, "_lat", lat_res_i,"_2018_2022_2021_2023_", freq, "kHz.rds")
    # path_nasc <- paste0("F:/data_elise/NASC/NASC_all_ESU/NASC_per_ESU_2018_2021_2023_", freq, "kHz.rds") # NASC per ESU
    
    # Open NASC
    nasc_ds <- readRDS(path_nasc)
    time_nasc <- nasc_ds$time
    
    date_nasc_all <- as.Date(time_nasc)
    # print(date_nasc_all)
    date_nasc_unique <- unique(date_nasc_all)
    date_nasc_unique
    lat_nasc <- nasc_ds$lat
    lon_nasc <- nasc_ds$lon
    
    # For every ping in nasc, find the closest fod cluster (in time and space)
    # Since FOD has daily resolution, loop is on days
    fod_all <- rep(NA, length(time_nasc))
    lat_fod_match <- rep(NA, length(time_nasc))
    lon_fod_match <- rep(NA, length(time_nasc))
    
    for (d in date_nasc_unique){
      # print(head(time_fod))
      d <- as.Date(d, origin = "1970-01-01")
      # print(d)
      idx_time_fod <- which(date_fod == d)
      
      if(length(idx_time_fod) == 0){
        print(paste("Jour non trouvé dans FOD :", format(d, "%Y-%m-%d")))
        next
      } # regarder les jours pas présents dans le ds
      
      # keep only ESU at day d
      mask_nasc_day <- date_nasc_all == d
      idx_nasc_day <- which(mask_nasc_day)
      
      lat_nasc_day <- lat_nasc[mask_nasc_day]
      lon_nasc_day <- lon_nasc[mask_nasc_day]
      
      # find closest lat/lon nasc coords in FOD
      for(i in seq_along(lat_nasc_day)){
        idx_lon_day_fod <- which.min(abs(lon_fod - lon_nasc_day[i]))
        idx_lat_day_fod <- which.min(abs(lat_fod - lat_nasc_day[i]))
        fod_cl <- clusters[idx_lon_day_fod, idx_lat_day_fod, idx_time_fod]
        
        fod_all[idx_nasc_day[i]] <- clusters[
          idx_lon_day_fod,
          idx_lat_day_fod,
          idx_time_fod
        ]
        
        # coordonnées FOD associées
        # print(c(
        #   "NASC lat =", lat_nasc_day[i],
        #   "FOD lat =", lat_fod[idx_lat_day_fod],
        #   "diff lat =", lat_fod[idx_lat_day_fod] - lat_nasc_day[i]
        # ))
        # 
        # print(c(
        #   "NASC lon =", lon_nasc_day[i],
        #   "FOD lon =", lon_fod[idx_lon_day_fod],
        #   "diff lon =", lon_fod[idx_lon_day_fod] - lon_nasc_day[i]
        # ))
        lat_fod_match[idx_nasc_day[i]] <- lat_fod[idx_lat_day_fod]
        lon_fod_match[idx_nasc_day[i]] <- lon_fod[idx_lon_day_fod]
      }
    }
  # str(fod_all)
  # print(unique(fod_all))
  length(lat_fod_match)
  length(lat_nasc)
  r <- lat_fod_match - lat_nasc
  # summary(r)
  
  nasc_fod_match <- data.frame(
    time_nasc = time_nasc,
    lat_nasc = lat_nasc,
    lon_nasc = lon_nasc,
    lat_fod = lat_fod_match,
    lon_fod = lon_fod_match,
    fod_cluster = fod_all
  )
  
  # str(nasc_fod_match)
  
  
  # VERIF
  print(unique(nasc_fod_match$fod_cluster)) # pas de cluster 8 c-a-d de transition 1-4
  print(unique(as.vector(clusters)))
  
  # verif du match des lon/lat
  dlat_fod <- nasc_fod_match$lat_fod - nasc_fod_match$lat_nasc
  dlon_fod <- nasc_fod_match$lon_fod - nasc_fod_match$lon_nasc
  print(summary(dlat_fod))
  print(summary(dlon_fod)) # OK
  
  
  saveRDS(
    nasc_fod_match,
    paste0("F:/data_elise/fod_elise_2018_2021_2022_2023/fod_colocated_nasc_2018_2021_2022_2023_transect/fod_colocated_NASC_mean_pig_grid/NASC_mean_pig_FOD_with_transitions_cluster_match_2018_2021_2023_", freq, "kHz_lon", lon_res_i, "_lat", lat_res_i, ".rds")
    )
  }
}




