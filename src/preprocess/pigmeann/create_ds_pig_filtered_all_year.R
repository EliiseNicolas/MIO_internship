# Description 

# On veut créer un dataset de pigments qui contient toutes les données pigmentaires des années 2018 2021 2022 2023 aux dates des transects / fod_clusters
# A partir de la liste de dates et des coordonnées lon/lat des FOD on veut cropper les données de pigments et les rassembler dans un meme RDS 

# Libraries 
library(ncdf4)

# Global variables
lat_fod <- "F:/data_elise/fod_elise_2018_2021_2022_2023/mclust_results_2018_2021_2022_2023/lat.rds"
lon_fod <- "F:/data_elise/fod_elise_2018_2021_2022_2023/mclust_results_2018_2021_2022_2023/lon.rds"
time_fod <- "F:/data_elise/fod_elise_2018_2021_2022_2023/mclust_results_2018_2021_2022_2023/time.rds"
folder_path <- "G:/data_elise/raw/PIGMeANN/daily"

# Ouvrir lat lon time FOD
lat <- readRDS(lat_fod); lon <- readRDS(lon_fod); time <- readRDS(time_fod)
lat_min <- min(lat); lat_max <- max(lat)
lon_min <- min(lon); lon_max <- max(lon)

# Convertir les dates en str
dates <- unique(format(time, "%Y%m%d"))
print(dates)

# lire les fichiers pigments correspondants aux dates de fod
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

file_dates <- unique(
  sub(".*?(\\d{8}).*", "\\1", basename(files))
)
dates_sans_fichier <- setdiff(dates, file_dates)
print(dates_sans_fichier)
print(files)
# dates sans fichier : [1] "20230109" "20230110" "20230111" "20230112" "20230113" "20230114" "20230115" "20230116" "20230117" "20230118" "20230119" "20230120" "20230121" 
# "20230122" "20230123" "20230124" "20230125" "20230301" "20230302" "20230303"

# Ouvrir chaque fichier, cropper aux lon/lat min/max 
pigments <- ""

f <- nc_open(files[1])
# get lat/lon
lon_pig <- ncvar_get(f, "lon")
lat_pig <- ncvar_get(f, "lat")

lon_idx <- which(lon_pig >= lon_min & lon_pig <= lon_max)
lat_idx <- which(lat_pig >= lat_min & lat_pig <= lat_max)

range(lon_pig[lon_idx])
range(lat_pig[lat_idx])

# extraire chaque variable 
var_names <- names(f$var)

var_names <- var_names[
  grepl("^c_cond_", var_names) |
    grepl("^use_", var_names) |
    var_names == "in_domain"
]

print(var_names)

pig_names <- c(
  "Chla", "Per", "But", "Fuco", "Hex",
  "Allo", "Zea", "Chlb", "DvChla"
)

n_dates <- length(files)

# Dimensions de la grille
n_lon <- lon_count
n_lat <- lat_count

pigments <- list()

pigments$lon <- NULL
pigments$lat <- NULL
pigments$date <- as.Date(character(n_dates))

# Créer les tableaux 3D
for (pig in pig_names) {
  pigments[[paste0("c_cond_", pig)]] <- array(
    NA_real_,
    dim = c(n_dates, n_lon, n_lat)
  )
}


for (i in seq_along(files)) {
  
  f <- files[i]
  
  ds <- nc_open(f)
  
  # ---------------------------------------------------
  # DATE
  # ---------------------------------------------------
  
  file_date <- sub(
    ".*?(\\d{8}).*",
    "\\1",
    basename(f)
  )
  
  date <- as.Date(
    file_date,
    format = "%Y%m%d"
  )
  
  pigments$date[i] <- date
  
  
  # ---------------------------------------------------
  # COORDONNEES
  # ---------------------------------------------------
  
  nc_lon <- ncvar_get(ds, "lon")
  nc_lat <- ncvar_get(ds, "lat")
  
  lon_idx <- which(
    nc_lon >= lon_min &
      nc_lon <= lon_max
  )
  
  lat_idx <- which(
    nc_lat >= lat_min &
      nc_lat <= lat_max
  )
  
  lon_start <- min(lon_idx)
  lon_count <- length(lon_idx)
  
  lat_start <- min(lat_idx)
  lat_count <- length(lat_idx)
  
  lon_crop <- nc_lon[lon_idx]
  lat_crop <- nc_lat[lat_idx]
  
  
  # Stocker les coordonnées une seule fois
  if (is.null(pigments$lon)) {
    pigments$lon <- lon_crop
    pigments$lat <- lat_crop
  }
  
  
  # ---------------------------------------------------
  # IN_DOMAIN
  # ---------------------------------------------------
  
  in_domain <- ncvar_get(
    ds,
    "in_domain",
    start = c(
      lon_start,
      lat_start,
      1
    ),
    count = c(
      lon_count,
      lat_count,
      1
    )
  )
  
  
  # ---------------------------------------------------
  # PIGMENTS
  # ---------------------------------------------------
  
  for (pig in pig_names) {
    
    c_cond <- ncvar_get(
      ds,
      paste0("c_cond_", pig),
      start = c(
        lon_start,
        lat_start,
        1
      ),
      count = c(
        lon_count,
        lat_count,
        1
      )
    )
    
    use <- ncvar_get(
      ds,
      paste0("use_", pig),
      start = c(
        lon_start,
        lat_start,
        1
      ),
      count = c(
        lon_count,
        lat_count,
        1
      )
    )
    
    
    # -------------------------------------------------
    # FILTRAGE 1%-99%
    # -------------------------------------------------
    
    q <- quantile(
      c_cond,
      c(0.01, 0.99),
      na.rm = TRUE
    )
    
    c_cond[
      c_cond < q[1] |
        c_cond > q[2]
    ] <- NA
    
    
    # -------------------------------------------------
    # FILTRAGE USE
    # -------------------------------------------------
    
    c_cond[use == 0] <- 0
    
    
    # -------------------------------------------------
    # FILTRAGE IN_DOMAIN
    # -------------------------------------------------
    
    c_cond[in_domain == 0] <- NA
    
    
    # -------------------------------------------------
    # STOCKAGE
    # -------------------------------------------------
    
    pigments[[paste0("c_cond_", pig)]][
      i, , 
    ] <- c_cond
  }
  
  nc_close(ds)
  
  cat(
    i, "/", length(files),
    "-", format(date, "%Y-%m-%d"),
    "\n"
  )
}

str(pigments)

saveRDS(pigments, "F:/data_elise/pigmeann/pigments_2018_2021_2022_2023_crop.rds")
