# Descritption

# From pigmeann rds, fod rds and nasc rds, create a new ds (entry to regression tree)
# columns : time_nasc, lat_nasc, lon_nasc, lat_fod, lon_fod, lat_pig, lon_pig, nasc, fod1, fod2, ..., fod6, pigm1, pigm2, ... pigm9

rm(list = ls())
library(rpart)
# ---------- Path et ouverture des ds
path_pig <- "F:/data_elise/pigmeann/pigs_transect_2018_2021_2023_120kHz/pig_mask9_1d_2018_2021_2023_120kHz_day.rds"
path_fod <- "F:/data_elise/fod_elise_2018_2021_2022_2023/fod_colocated_nasc_2018_2021_2023_transect/NASC_FOD_cluster_match_2018_2021_2023_day_120kHz.rds"
path_nasc <- "F:/data_elise/NASC/120kHz/NASC_mean_pig_grid_by_year_2018_2021_2023_day_120kHz.rds"
path_ftle <- "F:/data_elise/ftle/ftle_colocated_transect/120kHz/ftle_9_1d__colocated_with_NASC_transect_2018_2021_2023_120kHz_day.rds"

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
  nasc = nasc$NASC
)

#-------------- Filtrer les lats et lons qui ne nous intéressent pas
# ds final
# regression_ds<- regression_ds[regression_ds$lat_nasc < -30, ]
# regression_ds<- regression_ds[regression_ds$lon_nasc > 40, ]
# 
# # ds nasc
# nasc<- nasc[nasc$lat < -30, ]
# nasc<- nasc[nasc$lon > 40, ]
# print(length(unique(as.Date(nasc$time))))
# # ds pigmeann
# pig<- pig[pig$lat_sv < -30, ]
# pig<- pig[pig$lon_sv > 40, ]
# 
# # ds ftle
# ftle<- ftle[ftle$lat_sv < -30, ]
# ftle<- ftle[ftle$lon_sv > 40, ]
# 
# # ds fod
# fod<- fod[fod$lat_nasc < -30, ]
# fod<- fod[fod$lon_nasc > 40, ]


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
print(sum(!is.na(fod$lat_fod)))
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
regression_ds$lat_pig <- pig$lat_pig
regression_ds$lon_pig <- pig$lon_pig
regression_ds$Chla <- pig$Chla
regression_ds$Per <- pig$Per
regression_ds$But <- pig$But
regression_ds$Fuco <- pig$Fuco
regression_ds$Hex <- pig$Hex
regression_ds$Allo <- pig$Allo
regression_ds$Zea <- pig$Zea
regression_ds$Chlb <- pig$Chlb
regression_ds$DvChla <- pig$DvChla

# Vérifier que l'ordre NASC a été conservé
all(
  regression_ds$time_nasc == nasc$time &
    regression_ds$lat_nasc == nasc$lat &
    regression_ds$lon_nasc == nasc$lon
)

# ---- FTLE 
regression_ds$lat_ftle <- ftle$lat_ftle
regression_ds$lon_ftle <- ftle$lon_ftle
regression_ds$ftle <- ftle$ftle

# ---------------------- VERIFS
str(regression_ds)
all(regression_ds$time_nasc == nasc$time)
all(regression_ds$lat_nasc == nasc$lat)
all(regression_ds$lon_nasc == nasc$lon)

all(regression_ds$time_nasc == ftle$time)
all(regression_ds$lat_nasc == ftle$lat_sv)
all(regression_ds$lon_nasc == ftle$lon_sv)

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
saveRDS(regression_ds, "F:/data_elise/ds_NASC_pig_ftle_fod/ds_NASC_pig_ftle_fod_2018_2021_2023_transect_120kHz_day_mask9.rds")
