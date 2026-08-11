# Description
# From obsaustral ncdf transect data 2018, 2021, 2023 : 
#   1) filter by channel, keep only 200kHz data
#   2) cut data on depth dimension to remove most NA
#   3) filter lon >40°E and lat <-30°N
#   4) filter day/night 
# => save intermediary rds file : 2018_day.rds, 2018_night.rds, 2021_day.rds,...

#   5) Concatenate each year of same period (i.e. 2018_day, 2021_day, 2023_day)
# => save intermediary rds file : 2018_2021_2022_2023_day.rds, 2018_2021_2022_2023_night.rds

#   6) get pigmeann grid (check if pigmeann grid consistent on every year of data)
#   7) mean profile by grid (for day and night data)
# => save intermediary rds file : mean_pig_grid_2018_2021_2023_day.rds, ...

#   8) compute NASC for each mean profile 
# => save intermediary rds file : nasc_mean_pig_grid_2018_2021_2023_day.rds,..

# Packages
library(ncdf4)
rm (list=ls())
# Paths
path2018 <- "/run/media/mmolinet/KER22/data_elise/raw/acooustic/LOCEAN_SOOP-BA_A_20180105T121559Z_MARIONDUFRESNE_FV02_EchointegrationAcoustic-18-38-70-120-200_END-20180201T103636Z_C-20260522T153853Z.nc"
path2021 <- "/run/media/mmolinet/KER22/data_elise/raw/acooustic/LOCEAN_SOOP-BA_A_20210122T143044Z_MARIONDUFRESNE_FV02_EchointegrationAcoustic-18-38-70-120-200_END-20210307T041919Z_C-20260609T160140Z.nc"
path2023 <- "/run/media/mmolinet/KER22/data_elise/raw/acooustic/LOCEAN_SOOP-BA_A_20230123T103153Z_MARIONDUFRESNE_FV02_EchointegrationAcoustic-18-38-70-120-200_END-20230227T021804Z_C-20260728T105027Z.nc"

# Hyperparameters
year <- 2023
freq <- 200
thr_depth <- 165 #m

# Part I - filtering
ds <- nc_open(path2023)
Sv <- ncvar_get(ds, "Sv") # (depths, time, frequency)
print(dim(Sv))
time <- ncvar_get(ds, "time")
depth <- ncvar_get(ds, "depth")

# filter by channel
idx_freq <- which(ncvar_get(ds, "instrument_frequency")==freq)
Sv <- Sv[,,idx_freq] # (depths, time,)

# filter by latitude and longitude
lat <- ncvar_get(ds, "latitude") # (time,)
lon <- ncvar_get(ds, "longitude") # (time,)

idx <- which(lon > 40 & lat < -30)

n_total <- length(lon)
n_gardees <- length(idx)
n_supprimees <- n_total - n_gardees

cat("Latitudes totales :", n_total, "\n")
cat("Latitudes gardées :", n_gardees, "\n")
cat("Latitudes supprimées :", n_supprimees, "\n")
cat("Proportion gardée :", n_gardees / n_total, "\n")
cat("Proportion supprimée :", n_supprimees / n_total, "\n")


lat <- lat[idx]
lon <- lon[idx]
Sv <- Sv[,idx]
print(dim(Sv))
time <- time[idx]
# filter depth 
depth_idx <- which.min(abs(depth - thr_depth))
depth <- depth[1:depth_idx]
Sv <- Sv[1:depth_idx,]

# filter by diurnal period (day/night)
diurnal_period <- ncvar_get(ds, "day")
diurnal_period <- diurnal_period[idx] # filtrage lat/lon

idx_day <- which(diurnal_period == 3)
idx_night <- which(diurnal_period == 1)

print(dim(Sv))
print(idx_day)
Sv_day <- Sv[, idx_day]
Sv_night <- Sv[, idx_night]
print(dim(Sv_day))
print(dim(Sv_night))

time_day <- time[idx_day]
lat_day <- lat[idx_day]
lon_day <- lon[idx_day]

time_night <- time[idx_night]
lat_night <- lat[idx_night]
lon_night <- lon[idx_night]

nc_close(ds)

df_day <- list(
  profiles = t(Sv_day),
  lat = lat_day,
  lon = lon_day,
  depth = depth,
  time = time_day
)

df_night <- list(
  profiles = t(Sv_night),
  lat = lat_night,
  lon = lon_night,
  depth = depth,
  time = time_night
)

rm(Sv, Sv_day, Sv_night, ds)

# Save intermediary rds files
saveRDS(df_day, paste0("/run/media/mmolinet/KER22/MIO_internship_III/data_preprocessed/NASC/transect_2018_2022_2023/Sv_day_", freq, "kHz_", year, ".rds"))
saveRDS(df_night, paste0("/run/media/mmolinet/KER22/MIO_internship_III/data_preprocessed/NASC/transect_2018_2022_2023/Sv_night_", freq, "kHz_", year, ".rds"))

