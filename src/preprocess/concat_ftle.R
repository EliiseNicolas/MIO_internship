# Description

# Concatenation of ftle rds of 2021 2022 2023, created with ftle.Rmd script.

library(dplyr)

# Paths
path2021 <- "/run/media/elise/KER22/données elise/elisou_ta_stagiaire_pref/prepross/FTLE/ftle2021.rds"
path2022 <- "/run/media/elise/KER22/données elise/elisou_ta_stagiaire_pref/prepross/FTLE/ftle2022.rds"
path2023 <- "/run/media/elise/KER22/données elise/elisou_ta_stagiaire_pref/prepross/FTLE/ftle2023.rds"

# read ds and concat
ftle2021 <- readRDS(path2021)
ftle2022 <- readRDS(path2022)
ftle2023 <- readRDS(path2023)

ftle_all <- bind_rows(
  ftle2021,
  ftle2022,
  ftle2023
)
ftle_all

# filtrer : retirer sv_points qui ont lat > -30°N et lon <40°E
ftle_all<- ftle_all[ftle_all$lat_sv < -30, ]
ftle_all<- ftle_all[ftle_all$lon_sv > 40, ]
print(ftle_all)

# VERIFS
range(ftle_all$time)

nrow(ftle2021)
nrow(ftle2022)
nrow(ftle_all)

# verif ecart de lat/lon entre sv et ftle
dlat <- ftle_all$lat_ftle - ftle_all$lat_sv
dlon <- ftle_all$lon_ftle - ftle_all$lon_sv

summary(dlon) # OK
summary(dlat) # OK

saveRDS(
  ftle_all,
  file = "/run/media/elise/KER22/données elise/elisou_ta_stagiaire_pref/prepross/FTLE/ftle_mask9_1d_2021_2022_2023.rds"
)
