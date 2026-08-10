library(ncdf4)

path_2023_station <- "G:/données elise/raw/acoustic/Acoustique_obsaus_2023_-100_-50dB_50ESU_10m.nc"
path_2023_transect <- "G:/données elise/raw/acoustic/LOCEAN_SOOP-BA_A_20230123T103153Z_MARIONDUFRESNE_FV02_EchointegrationAcoustic-18-38-70-120-200_END-20230227T021804Z_C-20260728T105027Z.nc"

ds_st <- nc_open(path_2023_station)
lat <- ncvar_get(ds_st, "latitude")
lon <- ncvar_get(ds_st, "longitude")
time <- ncvar_get(ds_st, "time")
ds_st$dim$time$units
dates <- as.POSIXct(
  "1950-01-01 00:00:00",
  tz = "UTC"
) + time * 86400
print(c(min(lat), max(lat), min(lon), max(lon)))
print(range(dates))
nc_close(ds_st)

ds_tr <- nc_open(path_2023_transect)
lat <- ncvar_get(ds_tr, "latitude")
lon <- ncvar_get(ds_tr, "longitude")
time <- ncvar_get(ds_tr, "time")
ds_tr$dim$time$units
dates <- as.POSIXct(
  "1950-01-01 00:00:00",
  tz = "UTC"
) + time * 86400
print(c(min(lat), max(lat), min(lon), max(lon)))
print(range(dates))
nc_close(ds_tr)