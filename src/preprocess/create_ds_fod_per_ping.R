# Description

# Create .rds with following columns : time, lat, lon, cluster 1, cluster 2, ..., cluster n 
# The number of rows corresponds to the number of ESU collected during oceanographic 
# campaigns of interest (2021 2022 2023)

rm(list = ls())
# Librairies
library(ncdf4)

# Paths
path_gmm_2021_2022_2023 <- "/run/media/mmolinet/KER22/données elise/elisou_ta_stagiaire_pref/prepross/fod/FOD_2021_2022_2023/FOD_results_complete_2021_2022_2023.rds"
path_nasc_2021_2022_2023 <- "/run/media/mmolinet/KER22/données elise/elisou_ta_stagiaire_pref/prepross/NASC/ds_nasc_2021_2022_2023.rds"

# Open RDS
gmm_2021_2022_2023 <- readRDS(path_gmm_2021_2022_2023)
cls <- gmm_2021_2022_2023$cluster_map #(lon, lat, time)

lon_fod <- gmm_2021_2022_2023$lon
lat_fod <- gmm_2021_2022_2023$lat
time_fod <- gmm_2021_2022_2023$time
time_fod <- as.Date(time_fod)


# Open NASC
nasc_ds <- readRDS(path_nasc_2021_2022_2023)
nasc_ds <- nasc_ds[nasc_ds$lat < -30 & nasc_ds$lon > 40, ]
time_nasc <- nasc_ds$time

date_nasc_all <- as.Date(time_nasc)
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
  d <- as.Date(d, origin = "1970-01-01") # remettre dans le bon format car boucle for le détruit
  
  idx_time_fod <- which(time_fod == d)
  
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
    fod_cl <- cls[idx_lon_day_fod, idx_lat_day_fod, idx_time_fod]
   
    fod_all[idx_nasc_day[i]] <- cls[
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
print(unique(fod_all))
length(lat_fod_match)
length(lat_nasc)
r <- lat_fod_match - lat_nasc
summary(r)

nasc_fod_match <- data.frame(
  time_nasc = time_nasc,
  lat_nasc = lat_nasc,
  lon_nasc = lon_nasc,
  lat_fod = lat_fod_match,
  lon_fod = lon_fod_match,
  fod_cluster = fod_all
)

saveRDS(
  nasc_fod_match,
  "/run/media/mmolinet/KER22/données elise/elisou_ta_stagiaire_pref/prepross/fod/NASC_FOD_cluster_match_2021_2023.rds"
)

# rm(list = ls())
fod <- readRDS("/run/media/mmolinet/KER22/données elise/elisou_ta_stagiaire_pref/prepross/fod/NASC_FOD_cluster_match_2021_2023.rds")
str(fod)
fod<- fod[fod$lat_nasc < -30, ]
fod<- fod[fod$lon_nasc > 40, ]

dlat_fod <- fod$lat_fod - fod$lat_nasc
dlon_fod <- fod$lon_fod - fod$lon_nasc
summary(dlat_fod)
summary(dlon_fod) # OK
unique(ds$fod_cluster)

