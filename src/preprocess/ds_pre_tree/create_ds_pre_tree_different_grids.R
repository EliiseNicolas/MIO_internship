# Description

# From pigmeann rds, fod rds and nasc rds, create a new ds (entry to regression tree)
# columns : time_nasc, lat_nasc, lon_nasc, lat_fod, lon_fod, lat_pig, lon_pig, nasc, fod, pigm1, pigm2, ... pigm9, ftle

# Libraries
library(rpart)

# Global variables
rm(list = ls())

freqs <- c(18, 38, 70, 120, 200)
lat_res <- c(200, 500, 700, 1000)
lon_res <- c(200, 500, 1000, 1500)

for (g in seq_along(lat_res)) {
  lon_res_i <- lon_res[g]
  lat_res_i <- lat_res[g]
  grid_label <- paste0("lon", lon_res_i, "_lat", lat_res_i)
  
  for (freq in freqs){
    
    # ds NASC per ESU, par grille
    path_pig  <- paste0("F:/data_elise/pigmeann/pigs_colocated_NASC_mean_pig_grid/NASC_mean_pig_grid_pig_conc_ratio_nearest_point_2018_2021_2022_2023_", freq, "kHz_", grid_label, ".rds")
    path_fod  <- paste0("F:/data_elise/fod_elise_2018_2021_2022_2023/fod_colocated_nasc_2018_2021_2022_2023_transect/fod_colocated_NASC_mean_pig_grid/NASC_mean_pig_FOD_with_transitions_cluster_match_2018_2021_2023_", freq, "kHz_", grid_label, ".rds")
    path_nasc <- paste0("F:/data_elise/NASC/NASC_pig_mean/NASC_mean_Sv_pig_grid_", grid_label, "_2018_2022_2021_2023_", freq, "kHz.rds")
    path_ftle <- paste0("F:/data_elise/ftle/ftle_colocated_transect/ftle_colocated_NASC_mean_grid_pig/ftle_colocated_with_NASC_mean_grid_2018_2021_2022_2023_", freq, "kHz_", grid_label, ".rds")
    
    print(path_nasc)
    
    pig  <- readRDS(path_pig)
    fod  <- readRDS(path_fod)
    nasc <- readRDS(path_nasc)
    ftle <- readRDS(path_ftle)
    
    str(pig)
    str(fod)
    str(nasc)
    str(ftle)
    print(length(unique(as.Date(nasc$time))))
    
    # --------------- Creation du dataset final
    regression_ds <- data.frame(
      time_nasc = nasc$time,
      lat_nasc  = nasc$lat,
      lon_nasc  = nasc$lon,
      day       = nasc$day,
      nasc      = nasc$NASC
    )
    
    print(length(regression_ds$time_nasc))
    print(length(fod$time_nasc))
    
    # ---- Diagnostic taille
    nrow(nasc)
    nrow(fod)
    length(regression_ds$time_nasc)
    length(fod$time_nasc)
    
    # ---- Construire une clé unique par ESU (time + lat + lon)
    key_nasc <- paste(
      regression_ds$time_nasc,
      regression_ds$lat_nasc,
      regression_ds$lon_nasc,
      sep = "_"
    )
    
    key_fod <- paste(
      fod$time_nasc,
      fod$lat_nasc,
      fod$lon_nasc,
      sep = "_"
    )
    
    # ---- Trouver les lignes présentes dans nasc mais absentes de fod
    missing_in_fod <- key_nasc[!(key_nasc %in% key_fod)]
    print(missing_in_fod)
    length(missing_in_fod)
    
    # ---- Trouver les lignes présentes dans fod mais absentes de nasc (si jamais)
    missing_in_nasc <- key_fod[!(key_fod %in% key_nasc)]
    print(missing_in_nasc)
    length(missing_in_nasc)
    
    # ---- Voir les lignes complètes correspondantes (plus lisible)
    regression_ds[!(key_nasc %in% key_fod), c("time_nasc", "lat_nasc", "lon_nasc")]
    #----------------  Vérification alignement
    stopifnot(all( # verif fod
      regression_ds$time_nasc == fod$time_nasc &
        regression_ds$lat_nasc == fod$lat_nasc &
        regression_ds$lon_nasc == fod$lon_nasc
    ))
    
    stopifnot(all( # verif pig
      regression_ds$lat_nasc == pig$lat_sv &
        regression_ds$lon_nasc == pig$lon_sv
    ))
    
    stopifnot(all( # verif ftle
      regression_ds$time_nasc == ftle$time &
        regression_ds$lat_nasc == ftle$lat_sv &
        regression_ds$lon_nasc == ftle$lon_sv
    ))
    
    # ------------------- Ajout des variables ftle, pig et fod dans le dataset final
    # ---- FOD
    print(sum(!is.na(fod$lat_fod)))
    regression_ds$lat_fod <- fod$lat_fod
    regression_ds$lon_fod <- fod$lon_fod
    regression_ds$fod     <- format(fod$fod_cluster)
    
    # ----- Pigments
    vars_pig <- setdiff(
      names(pig),
      c("time", "lat_sv", "lon_sv")
    )
    print(vars_pig)
    regression_ds[vars_pig] <- pig[vars_pig]
    
    # Vérifier que l'ordre NASC a été conservé
    stopifnot(all(
      regression_ds$time_nasc == nasc$time &
        regression_ds$lat_nasc == nasc$lat &
        regression_ds$lon_nasc == nasc$lon
    ))
    
    # ---- FTLE
    regression_ds$lat_ftle <- ftle$lat_ftle
    regression_ds$lon_ftle <- ftle$lon_ftle
    regression_ds$ftle     <- ftle$ftle
    
    # ---------------------- VERIFS
    str(regression_ds)
    print(all(regression_ds$time_nasc == nasc$time))
    print(all(regression_ds$lat_nasc == nasc$lat))
    print(all(regression_ds$lon_nasc == nasc$lon))
    
    print(all(regression_ds$time_nasc == ftle$time))
    print(all(regression_ds$lat_nasc == ftle$lat_sv))
    print(all(regression_ds$lon_nasc == ftle$lon_sv))
    
    # verif années contenues
    print(unique(format(regression_ds$time_nasc, "%Y")))
    
    # verif nb de ESU
    print(nrow(nasc) == nrow(regression_ds)) # OK toutes les ESU conservées
    
    # verif nasc non modifié
    print(all.equal(
      regression_ds$nasc,
      nasc$NASC
    ))
    
    # différence de lat/lon entre nasc, fod et pig
    dlat_fod <- regression_ds$lat_fod - regression_ds$lat_nasc
    dlon_fod <- regression_ds$lon_fod - regression_ds$lon_nasc
    print(summary(dlat_fod))
    print(summary(dlon_fod))
    
    dlat_pig <- regression_ds$lat_pig - regression_ds$lat_nasc
    dlon_pig <- regression_ds$lon_pig - regression_ds$lon_nasc
    print(summary(dlat_pig))
    print(summary(dlon_pig))
    
    dlat_ftle <- regression_ds$lat_ftle - regression_ds$lat_nasc
    dlon_ftle <- regression_ds$lon_ftle - regression_ds$lon_nasc
    print(summary(dlat_ftle))
    print(summary(dlon_ftle))
    
    # ------------ SAVE (un fichier par grille et par fréquence)
    saveRDS(
      regression_ds,
      paste0(
       "F:/data_elise/ds_NASC_pig_ftle_fod/ds_NASC_mean_pig_grid_all/ds_NASC_per_esu_pig_ftle_fod_2018_2021_2022_2023_transect_",
        freq, "kHz_", grid_label, ".rds"
      )
    )
  }
}