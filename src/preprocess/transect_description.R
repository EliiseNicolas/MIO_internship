# Description 

# Compute number of days we have in a transect

# Libraries
library(ncdf4)

# Paths
path2018tr <- "/run/media/mmolinet/KER22/données elise/raw/acoustic/LOCEAN_SOOP-BA_A_20180105T121559Z_MARIONDUFRESNE_FV02_EchointegrationAcoustic-18-38-70-120-200_END-20180201T103636Z_C-20260522T153853Z.nc"
path2021tr <- "/run/media/mmolinet/KER22/données elise/raw/acoustic/LOCEAN_SOOP-BA_A_20210122T143044Z_MARIONDUFRESNE_FV02_EchointegrationAcoustic-18-38-70-120-200_END-20210307T041919Z_C-20260609T160140Z.nc"
path2023tr <- "/run/media/mmolinet/KER22/données elise/raw/acoustic/LOCEAN_SOOP-BA_A_20230123T103153Z_MARIONDUFRESNE_FV02_EchointegrationAcoustic-18-38-70-120-200_END-20230227T021804Z_C-20260728T105027Z.nc"
path2021st <- "/run/media/mmolinet/KER22/données elise/raw/acoustic/Acoustique_obsaus_2021_-100_-50dB_50ESU_10m.nc"
path2022st <- "/run/media/mmolinet/KER22/données elise/raw/acoustic/Acoustique_obsaus_2022_-100_-50dB_50ESU_10m.nc"
path2023st <- "/run/media/mmolinet/KER22/données elise/raw/acoustic/Acoustique_obsaus_2023_-100_-50dB_50ESU_10m.nc"


nb_days <- function(path, year){
  
  nc <- nc_open(path)
  
  # Extraire les variables
  lat <- ncvar_get(nc, "latitude")
  lon <- ncvar_get(nc, "longitude")
  time <- as.POSIXct("1950-01-01 00:00:00", tz="UTC") +
    ncvar_get(nc, "time") * 86400
  day <- ncvar_get(nc, "day")
  nc_close(nc)
  
  # créer un dataframe des pings
  df <- data.frame(
    time = time,
    latitude = lat,
    longitude = lon, 
    day = day
  )
  
  # filtrage spatial
  print(length(unique(as.Date(df$time))))
  df <- df[df$latitude < -30 & df$longitude > 40, ]
  print(length(unique(as.Date(df$time))))
  df <- df[df$day == 3,]
  print(length(unique(as.Date(df$time))))
}

# nb_days(path2018st, 2018)
nb_days(path2021st, 2021)
nb_days(path2022st, 2022)
nb_days(path2023st, 2023)


compare_days <- function(path_st, path_tr){
  ds_st <- nc_open(path_st)
  lat_st <- ncvar_get(ds_st, "latitude")
  lon_st <- ncvar_get(ds_st, "longitude")
  time_st <- as.POSIXct("1950-01-01 00:00:00", tz="UTC") +
    ncvar_get(ds_st, "time") * 86400
  day_st <- ncvar_get(ds_st, "day")
  nc_close(ds_st)
  
  ds_tr <- nc_open(path_tr)
  lat_tr <- ncvar_get(ds_tr, "latitude")
  lon_tr <- ncvar_get(ds_tr, "longitude")
  time_tr <- as.POSIXct("1950-01-01 00:00:00", tz="UTC") +
    ncvar_get(ds_tr, "time") * 86400
  # print(as.Date(time_tr))
  day_tr <- ncvar_get(ds_tr, "day")
  nc_close(ds_tr)
  
  # filtrage
  keep_st <- which(lat_st < -30 & lon_st > 40)
  dates_st <- unique(as.Date(time_st[keep_st]))
  
  keep_tr <- which(lat_tr < -30 & lon_tr > 40)
  dates_tr <- unique(as.Date(time_tr[keep_tr]))
  
  # display results
  common_days <- intersect(dates_st, dates_tr)
  
  cat("Nombre de jours uniques ST :", length(dates_st), "\n")
  cat("Nombre de jours uniques TR :", length(dates_tr), "\n")
  print(paste("Jours uniques ST :", dates_st, "\n"))
  print(paste("Jours uniques TR :", dates_tr, "\n"))
  cat("Nombre de jours communs :", length(common_days), "\n")
  
}
compare_days(path2021st, path2021tr)
compare_days(path2023st, path2023tr)
