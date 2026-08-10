# Descritption

# From pigmeann rds, fod rds and nasc rds, create a new ds (entry to regression tree)
# columns : time_nasc, lat_nasc, lon_nasc, lat_fod, lon_fod, lat_pig, lon_pig, nasc, fod1, fod2, ..., fod6, pigm1, pigm2, ... pigm9

rm(list = ls())
library(rpart)
# ---------- Path et ouverture des ds
path_pig <- "/run/media/mmolinet/KER22/données elise/elisou_ta_stagiaire_pref/prepross/pigmeann/pig_mask100_1d_2021_2022_2023.rds"
path_fod <- "/run/media/mmolinet/KER22/données elise/elisou_ta_stagiaire_pref/prepross/fod/NASC_FOD_cluster_match_2021_2023.rds"
path_nasc <- "/run/media/mmolinet/KER22/données elise/elisou_ta_stagiaire_pref/prepross/NASC/ds_nasc_2021_2022_2023.rds"
path_ftle <- "/run/media/mmolinet/KER22/données elise/elisou_ta_stagiaire_pref/prepross/FTLE/ftle_mask9_1d_2021_2022_2023.rds"

pig <- readRDS(path_pig)
fod <- readRDS(path_fod)
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
  lat_nasc = nasc$lat,
  lon_nasc = nasc$lon,
  nasc = nasc$NASC,
  day = nasc$day
)

#-------------- Filtrer les lats et lons qui ne nous intéressent pas
# ds final
regression_ds<- regression_ds[regression_ds$lat_nasc < -30, ]
regression_ds<- regression_ds[regression_ds$lon_nasc > 40, ]

# ds nasc
nasc<- nasc[nasc$lat < -30, ]
nasc<- nasc[nasc$lon > 40, ]
print(length(unique(as.Date(nasc$time))))
# ds pigmeann
pig<- pig[pig$lat_sv < -30, ]
pig<- pig[pig$lon_sv > 40, ]

# ds ftle
ftle<- ftle[ftle$lat_sv < -30, ]
ftle<- ftle[ftle$lon_sv > 40, ]

# ds fod
fod<- fod[fod$lat_nasc < -30, ]
fod<- fod[fod$lon_nasc > 40, ]


#----------------  Vérification alignement
all( # verif fod
  regression_ds$time_nasc == fod$time_nasc &
    regression_ds$lat_nasc == fod$lat_nasc &
    regression_ds$lon_nasc == fod$lon_nasc
) # TRUE

all( # verif pig
    regression_ds$lat_nasc == pig$lat_sv &
    regression_ds$lon_nasc == pig$lon_sv
) # TRUE

all( # verif ftle
  regression_ds$time_nasc == ftle$time &
    regression_ds$lat_nasc == ftle$lat_sv &
    regression_ds$lon_nasc == ftle$lon_sv 
) # TRUE

# ------------------- Ajout des variables ftle, pig et fod dans le dataset final
# ---- FOD
# Ajout coordonnées FOD + cluster
regression_ds$lat_fod <- fod$lat_fod
regression_ds$lon_fod <- fod$lon_fod
regression_ds$fod <- fod$fod_cluster

# # Transformer cluster en variables fod1...fod6
# nclust <- 6
# for(i in 0:nclust){
#   regression_ds[[paste0("fod", i)]] <- 
#     as.integer(fod$fod_cluster == i)
# }
# sum(is.na(fod$fod_cluster))
# s <- rowSums(regression_ds[paste0("fod", 0:6)], na.rm = TRUE)
# which(s != 1)
# idx <- which(s != 1)
# regression_ds[idx, paste0("fod", 0:6)] # que des NA donc OK

# ----- Pigments
regression_ds <- merge(
  regression_ds,
  pig[, c(
    "time",
    "lat_pig",
    "lon_pig",
    "Chla",
    "Per",
    "But",
    "Fuco",
    "Hex",
    "Allo",
    "Zea",
    "Chlb",
    "DvChla"
  )],
  by.x = "time_nasc",
  by.y = "time",
  all.x = TRUE
)

# ---- FTLE 
regression_ds$lat_ftle <- ftle$lat_ftle
regression_ds$lon_ftle <- ftle$lon_ftle
regression_ds$ftle <- ftle$ftle

# ---------------------- VERIFS
str(regression_ds)

# verif années contenues
unique(format(regression_ds$time_nasc, "%Y")) # 2021 2022 2023 OK

# verif nb de ESU
nrow(nasc) == nrow(regression_ds) # OK toutes les ESU conservées

# verif nasc non modifié 
all.equal(
  regression_ds$nasc,
  nasc$NASC
) # TRUE

# différence de lat/lon entre nasc fod et pig
str(regression_ds)
dlat_fod <- regression_ds$lat_fod - regression_ds$lat_nasc
dlon_fod <- regression_ds$lon_fod - regression_ds$lon_nasc
summary(dlat_fod)
summary(dlon_fod) # OK

dlat_pig <- regression_ds$lat_pig - regression_ds$lat_nasc
dlon_pig <- regression_ds$lon_pig - regression_ds$lon_nasc
summary(dlat_pig)
summary(dlon_pig) # OK

dlat_ftle <- regression_ds$lat_ftle - regression_ds$lat_nasc
dlon_ftle <- regression_ds$lon_ftle - regression_ds$lon_nasc
summary(dlat_ftle)
summary(dlon_ftle) # OK

# ------------ SAVE
saveRDS(regression_ds, "/run/media/mmolinet/KER22/données elise/elisou_ta_stagiaire_pref/prepross/regression_ds_2021_2022_2023_m100.rds")
