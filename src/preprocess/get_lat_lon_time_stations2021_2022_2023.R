library(ncdf4)

path2021 <- "G:/données elise/raw/acoustic/Acoustique_obsaus_2021_-100_-50dB_50ESU_10m.nc"
path2022 <- "G:/données elise/raw/acoustic/Acoustique_obsaus_2022_-100_-50dB_50ESU_10m.nc"
path2023 <- "G:/données elise/raw/acoustic/Acoustique_obsaus_2023_-100_-50dB_50ESU_10m.nc"

l <- list(path2021, path2022, path2023)
table_coord <- data.frame()

for (p in l) {
  ds <- nc_open(p)

  lat <- ncvar_get(ds, "latitude")
  lon <- ncvar_get(ds, "longitude")
  time <- ncvar_get(ds, "time")
  
  nc_close(ds)
  
  # création du tableau
  tmp <- data.frame(
    lat = lat,
    lon = lon,
    time = time
  )

  print(length(lat))
  print(length(lon))
  print(length(time))
  
  # ajout aux données précédentes
  table_coord <- rbind(table_coord, tmp)
}
nrow(table_coord)
# creer netcdf
n <- nrow(table_coord)
dim_point <- ncdim_def(
  name = "point",
  units = "count",
  vals = 1:n
)

var_time <- ncvar_def(
  name = "time",
  units = "days since 1950-01-01 00:00:00 UTC",
  dim = list(dim_point),
  missval = -9999
)

var_lat <- ncvar_def(
  name = "latitude",
  units = "degrees_north",
  dim = list(dim_point),
  missval = -9999
)

var_lon <- ncvar_def(
  name = "longitude",
  units = "degrees_east",
  dim = list(dim_point),
  missval = -9999
)

# Créer le fichier
nc_out <- nc_create(
  "C:/Users/mmolinet/elisou_ta_stagiaire_pref/prepross/list_lat_lon_time_2021_2022_2023_stations.nc",
  list(var_lat, var_lon, var_time)
)
# Écrire les données
ncvar_put(nc_out, var_lat, table_coord$lat)
ncvar_put(nc_out, var_lon, table_coord$lon)
ncvar_put(nc_out, var_time, table_coord$time)

# Fermer
nc_close(nc_out)

# check 


# Ouvrir le nouveau fichier NetCDF
nc <- nc_open("C:/Users/mmolinet/elisou_ta_stagiaire_pref/prepross/list_lat_lon_time_2021_2022_2023_stations.nc")

# Lire les variables
lat <- ncvar_get(nc, "latitude")
lon <- ncvar_get(nc, "longitude")
time <- ncvar_get(nc, "time")
print(nc$var[["time"]])
nc_close(nc)

# Reconvertir le temps en date
date <- as.POSIXct(
  time * 86400,
  origin = "1950-01-01 00:00:00",
  tz = "UTC"
)

# Extraire l'année
year <- format(date, "%Y")

# Créer un dataframe
coord <- data.frame(
  lon = lon,
  lat = lat,
  year = year
)

head(coord)

all.equal(
  table_coord$lat,
  coord$lat
)

all.equal(
  table_coord$lon,
  coord$lon
)
all.equal(
  table_coord$time,
  as.vector(time)
)

summary(coord$lon)
summary(coord$lat)
table(coord$year)