# Description 

# On veut créer un dataset de ftle qui contient toutes les données ftle des années 2018 2021 2022 2023 aux dates des transects / fod_clusters
# A partir de la liste de dates et des coordonnées lon/lat des FOD on veut cropper les données de ftle et les rassembler dans un meme RDS 

rm(list=ls())
# Libraries 
library(ncdf4)

# Global variables
lat_fod <- "F:/data_elise/fod_elise_2018_2021_2022_2023/mclust_results_2018_2021_2022_2023/lat.rds"
lon_fod <- "F:/data_elise/fod_elise_2018_2021_2022_2023/mclust_results_2018_2021_2022_2023/lon.rds"
time_fod <- "F:/data_elise/fod_elise_2018_2021_2022_2023/mclust_results_2018_2021_2022_2023/time.rds"
folder_path <- "G:/data_elise/raw/ftle"

# Ouvrir lat lon time FOD
lat <- readRDS(lat_fod); lon <- readRDS(lon_fod); time <- readRDS(time_fod)
lat_min <- min(lat); lat_max <- max(lat)
lon_min <- min(lon); lon_max <- max(lon)

# Convertir les dates en str
dates <- unique(format(time, "%Y-%m-%d"))
print(dates)

# lire les fichiers ftle correspondants aux dates de fod
files <- list.files(
  path = folder_path,
  pattern = "\\.nc$",
  full.names = TRUE
)
print(files)

# On garde uniquement les fichiers correspondant aux dates
files <- files[
  grepl(
    paste(dates, collapse = "|"),
    basename(files)
  )
]

file_dates <- unique(sub(".*?(\\d{4}-\\d{2}-\\d{2}).*", "\\1", basename(files)))
dates_sans_fichier <- setdiff(dates, file_dates)
cat("Dates FOD sans fichier FTLE correspondant :", length(dates_sans_fichier), "/", length(dates), "\n")
print(dates_sans_fichier)
# "2018-02-02" "2018-02-03" "2018-02-04" "2018-02-05" "2018-02-06" "2018-02-07" "2018-02-08" "2018-02-09" "2018-02-10" "2018-02-11"
# [11] "2018-02-12" "2018-02-13" "2018-02-14" "2018-02-15" "2018-02-16" "2018-02-17" "2018-02-18" "2018-02-19" "2018-02-20" "2018-02-21"
# [21] "2018-02-22" "2018-02-23" "2018-02-24" "2018-02-25" "2018-02-26" "2018-02-27" "2018-02-28" "2018-03-01" "2018-03-02" "2018-03-03"
# [31] "2021-01-09" "2021-01-10" "2021-01-11" "2021-01-12" "2022-01-09" "2022-01-10" "2022-01-11" "2022-01-12" "2022-01-13" "2022-01-14"
# [41] "2022-01-15" "2022-01-16" "2022-01-17" "2022-01-18" "2022-01-19" "2022-01-20" "2022-01-21" "2022-01-22" "2022-01-23" "2022-01-24"
# [51] "2022-01-25" "2022-01-26" "2022-01-27" "2022-01-28" "2022-01-29" "2022-01-30" "2022-01-31" "2023-01-09" "2023-01-10" "2023-01-11"
# [61] "2023-01-12" "2023-01-13" "2023-01-14" "2023-01-15" "2023-01-16" "2023-01-17" "2023-01-18" "2023-01-19" "2023-01-20" "2023-01-21"
# [71] "2023-01-22" "2023-01-23" "2023-01-24" "2023-01-25" "2023-02-10" "2023-02-22" "2023-02-24" "2023-02-26" "2023-02-27" "2023-02-28"
# [81] "2023-03-01" "2023-03-02" "2023-03-03" pas focément un pb pcq pas forcément de nasc a cette periode

# Ouvrir chaque fichier, cropper aux lon/lat min/max 
ftle <- ""

# Ouvrir chaque fichier, cropper aux lon/lat min/max
# ============================================================
# Determination de la grille de sortie a partir du premier fichier
# (necessaire pour dimensionner les arrays AVANT la boucle - dans
# la version pigments, lon_count/lat_count etaient utilises avant
# d'etre definis, ce qui plantait).
# ============================================================

f_ref <- nc_open(files[1])

lon_ftle <- ncvar_get(f_ref, "lon")
lat_ftle <- ncvar_get(f_ref, "lat")

lon_idx_ref <- which(lon_ftle >= lon_min & lon_ftle <= lon_max)
lat_idx_ref <- which(lat_ftle >= lat_min & lat_ftle <= lat_max)

range(lon_ftle[lon_idx_ref])
range(lat_ftle[lat_idx_ref])

nc_close(f_ref)

n_lon <- length(lon_idx_ref)
n_lat <- length(lat_idx_ref)
n_dates <- length(files)

ftle <- list()
ftle$lon <- NULL
ftle$lat <- NULL
ftle$date <- as.Date(character(n_dates))
ftle$ftle <- array(NA_real_, dim = c(n_dates, n_lon, n_lat))

# ============================================================
# Boucle sur tous les fichiers
# ============================================================

for (i in seq_along(files)) {
  
  f <- files[i]
  ds <- nc_open(f)
  
  # ---------------------------------------------------
  # DATE (format YYYY-MM-DD dans le nom de fichier, pas YYYYMMDD)
  # ---------------------------------------------------
  
  file_date <- sub(".*?(\\d{4}-\\d{2}-\\d{2}).*", "\\1", basename(f))
  date <- as.Date(file_date, format = "%Y-%m-%d")
  ftle$date[i] <- date
  
  # ---------------------------------------------------
  # COORDONNEES
  # ---------------------------------------------------
  
  nc_lon <- ncvar_get(ds, "lon")
  nc_lat <- ncvar_get(ds, "lat")
  
  lon_idx <- which(nc_lon >= lon_min & nc_lon <= lon_max)
  lat_idx <- which(nc_lat >= lat_min & nc_lat <= lat_max)
  
  lon_start <- min(lon_idx); lon_count <- length(lon_idx)
  lat_start <- min(lat_idx); lat_count <- length(lat_idx)
  
  lon_crop <- nc_lon[lon_idx]
  lat_crop <- nc_lat[lat_idx]
  
  if (lon_count != n_lon || lat_count != n_lat) {
    cat("  ATTENTION dimensions differentes pour", file_date, ": ", lon_count, "x", lat_count,
        "(attendu", n_lon, "x", n_lat, ") -> fichier ignore\n")
    nc_close(ds)
    next
  }
  
  # Stocker les coordonnees une seule fois
  if (is.null(ftle$lon)) {
    ftle$lon <- lon_crop
    ftle$lat <- lat_crop
  }
  
  # ---------------------------------------------------
  # FTLE
  # ---------------------------------------------------
  # dims du fichier : (lon, lat, time), time toujours de taille 1
  # -> on lit un seul pas de temps et on droppe cette dimension.
  
  ftle_slab <- ncvar_get(
    ds, "FTLE",
    start = c(lon_start, lat_start, 1),
    count = c(lon_count, lat_count, 1)
  )
  # ftle_slab est deja 2D (lon x lat) car ncvar_get droppe les dims
  # de taille 1 par defaut -> pas besoin de drop manuel, mais on
  # verifie quand meme au cas ou (collapse = FALSE ferait un array 3D).
  if (length(dim(ftle_slab)) == 3) ftle_slab <- ftle_slab[, , 1]
  
  ftle$ftle[i, , ] <- ftle_slab
  
  nc_close(ds)
  
  cat(i, "/", length(files), "-", format(date, "%Y-%m-%d"), "\n")
}

str(ftle)

saveRDS(ftle, "G:/data_elise/processed/ftle_colocated_transect/ftle_2018_2021_2022_2023_cropped.rds")

n_dates_ok <- sum(!is.na(ftle$ftle[, 1, 1]))
cat(sprintf("Dates FTLE effectivement chargees : %d / %d\n", n_dates_ok, n_dates))