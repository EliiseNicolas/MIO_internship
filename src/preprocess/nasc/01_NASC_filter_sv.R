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

# Paths
path2018 <- ""

# Hyperparameters
year <- 2018
path <- ""
freq <- 200
thr_depth <- 165 #m

# Part I - filtering
ds <- nc_open(path)
Sv <- ncvar_get(ds, "Sv") # (time, depths, frequency)
time <- ncvar_get(ds, "time")
depth <- ncvar_get(ds, "depth")

# filter by channel
idx_freq <- which(ncvar_get(ds, "instrument_frequency")==freq)
Sv <- Sv[,,idx_freq] # (time, depths,)

# filter by latitude and longitude
lat <- ncvar_get(ds, "latitude") # (time,)
lon <- ncvar_get(ds, "longitude") # (time,)

idx <- which(lat < -30 & lon > 40)

lat <- lat[idx]
lon <- lon[idx]

Sv <- Sv[idx,]

# filter depth 
depth_idx <- which.min(abs(depth - thr_depth))
depth <- depth[1:depth_idx]
Sv <- Sv[,1:depth_idx]

# filter by diurnal period (day/night)
diurnal_period <- ncvar_get(ds, "day")
diurnal_period <- diurnal_period[idx] # filtrage lat/lon

idx_day <- which(diurnal_period == 3)
idx_night <- which(diurnal_period == 1)

Sv_day <- Sv[idx_day, ]
Sv_night <- Sv[idx_night, ]

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

# Save intermediary rds files
saveRDS(df_day, paste0("/Sv_day_", year, ".rds"))
saveRDS(df_night, paste0("/Sv_night_", year, ".rds"))

