library(ncdf4)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
rm(list = ls())

folder_path <- "/mnt/KER22/data_elise/raw/PIGMeANN/daily"
year <- 2018

# ---------------------------------------------------------
# PARAMETRES
# ---------------------------------------------------------

# Taille du voisinage autour de la station
# p = 1 : pixel le plus proche
# p = 3 : fenêtre 3 x 3
# p = 5 : fenêtre 5 x 5
# p = 10 : fenêtre 10 x 10
p <- 3

list_pigs <- c("Chla", "Per", "But", "Fuco", "Hex","Allo", "Zea", "Chlb", "DvChla")


# ---------------------------------------------------------
# NASC
# ---------------------------------------------------------
nasc_path <- "~/Documents/stage_MIO/pt_III/data_preprocessed/NASC/transect_2018_2022_2023/NASC_mean_pig_grid_by_year_2018_2021_2023_day_200kHz.rds"
nasc_ds <- readRDS(nasc_path)
print(head(nasc_ds$time))
dates_all_str <- format(as.Date(nasc_ds$time), "%Y%m%d")
dates_unique_str <- (unique(dates_all_str))

print(length(dates_unique_str)) # 75

# ---------------------------------------------------------
# FICHIERS PIGMeANN
# ---------------------------------------------------------

files <- list.files(
  path = folder_path,
  pattern = "\\.nc$",
  full.names = TRUE
)
print(files)

# On garde uniquement les fichiers correspondant aux dates
files <- files[
  grepl(
    paste(dates_unique_str, collapse = "|"),
    basename(files)
  )
]

print(length(files)) # 75 aussi donc OK
print(basename(files))


# ---------------------------------------------------------
# FONCTION : INDICES D'UNE FENETRE
# ---------------------------------------------------------

get_window <- function(index, n, p) {
  
  r <- floor(p / 2)
  
  i_min <- max(1, index - r)
  i_max <- min(n, index + r)
  
  return(i_min:i_max)
}

# ---------------------------------------------------------
# EXTRACTION PIGMeANN
# ---------------------------------------------------------

# Initialisation dataframe
list_pigs = c(
  "Chla", "Per", "But", "Fuco", "Hex",
  "Allo", "Zea", "Chlb", "DvChla"
)

pig_cond <- data.frame(
  time = nasc_ds$time,
  lat_sv = nasc_ds$lat,
  lon_sv = nasc_ds$lon,
  lat_pig = NA_real_,
  lon_pig = NA_real_
)

# Ajouter les colonnes pigments
pig_cond[list_pigs] <- NA_real_


for(date_i in dates_unique_str){
  
  # indice ESU à cette date
  ind <- which(dates_all_str == date_i)
  
  if(length(ind) == 0)
    next
  
  # ouvrir le fichier du jour i
  filename <- list.files(
    folder_path,
    pattern = paste0(".*", date_i, "\\.nc$"),
    full.names = TRUE
  )
  
  if(length(filename) == 0)
    next
  
  ds <- nc_open(filename[1])
  
  # récupérer les coordonnées PIGMeANN du jour
  lat_pig <- ds$dim$lat$vals
  lon_pig <- ds$dim$lon$vals
  
  nlat <- length(lat_pig)
  nlon <- length(lon_pig)
  
  # trouver les pixels PIGMeANN les plus proches des stations
  idx_lon <- sapply(
    nasc_ds$lon[ind],
    function(x) which.min(abs(lon_pig - x))
  )
  
  idx_lat <- sapply(
    nasc_ds$lat[ind],
    function(x) which.min(abs(lat_pig - x))
  )
  
  # enregistrer les coordonnées PIGMeANN associées
  pig_cond$lat_pig[ind] <- lat_pig[idx_lat]
  pig_cond$lon_pig[ind] <- lon_pig[idx_lon]
  
  
  # =====================================================
  # FENETRE LAT/LON A LIRE
  # =====================================================
  
  r <- (p - 1) / 2
  
  lon_start <- max(1, min(idx_lon) - r)
  lon_end   <- min(nlon, max(idx_lon) + r)
  
  lat_start <- max(1, min(idx_lat) - r)
  lat_end   <- min(nlat, max(idx_lat) + r)
  
  lon_count <- lon_end - lon_start + 1
  lat_count <- lat_end - lat_start + 1
  
  
  for (pig in list_pigs){
    
    # ===================================================
    # LIRE UNIQUEMENT LA FENETRE LAT/LON
    # ===================================================
    
    c_cond <- ncvar_get(
      ds,
      paste0("c_cond_", pig),
      start = c(lon_start, lat_start, 1),
      count = c(lon_count, lat_count, 1)
    )
    
    use <- ncvar_get(
      ds,
      paste0("use_", pig),
      start = c(lon_start, lat_start, 1),
      count = c(lon_count, lat_count, 1)
    )
    
    in_domain <- ncvar_get(
      ds,
      "in_domain",
      start = c(lon_start, lat_start, 1),
      count = c(lon_count, lat_count, 1)
    )
    
    
    # ===================================================
    # FILTRAGE
    # ===================================================
    
    q <- quantile(
      c_cond,
      c(0.01,0.99),
      na.rm = TRUE
    )
    
    c_cond[
      c_cond < q[1] |
        c_cond > q[2]
    ] <- NA
    
    c_cond[use == 0] <- 0
    c_cond[in_domain == 0] <- NA
    
    
    # ===================================================
    # EXTRACTION AU VOISINAGE DU POINT
    # ===================================================
    
    for(j in seq_along(ind)){
      
      i <- ind[j]
      
      # coordonnées de la station dans la fenêtre lue
      ilon <- idx_lon[j] - lon_start + 1
      ilat <- idx_lat[j] - lat_start + 1
      
      r <- (p - 1) / 2
      
      lon_win <- max(
        1,
        ilon - r
      ):min(
        nrow(c_cond),
        ilon + r
      )
      
      lat_win <- max(
        1,
        ilat - r
      ):min(
        ncol(c_cond),
        ilat + r
      )
      
      pig_cond[i, pig] <- mean(
        c_cond[lon_win, lat_win],
        na.rm = TRUE
      )
    }
    
    # libérer les petites matrices avant le pigment suivant
    rm(c_cond, use, in_domain)
    gc()
  }
  
  nc_close(ds)
  rm(ds)
  gc()
}

saveRDS(
  pig_cond,
  file = paste0(
    "pig_mask",
    p*p,
    "_1d_",
    year,
    ".rds"
  )
)

str(pig_cond)

pig_vars <- c(
  "Chla",
  "Fuco",
  "But",
  "Per",
  "Hex",
  "Allo",
  "Chlb",
  "Zea", 
  "DvChla"
)

p1 <- pig_cond %>%
  select(any_of(pig_vars)) %>%
  pivot_longer(
    everything(),
    names_to = "pigment",
    values_to = "value"
  ) %>%
  filter(is.finite(value)) %>%
  ggplot(aes(x = value)) +
  geom_histogram(bins = 50) +
  facet_wrap(~pigment, scales = "free") +
  theme_bw() +
  scale_x_continuous(n.breaks = 3)


# ---------------------------------------------------------
# Valeurs manquantes
# ---------------------------------------------------------

nb_missing <- sum(!is.finite(pig_cond[[pig_vars[1]]]))

pct_missing <- nb_missing / nrow(pig_cond) * 100

p2 <- ggplot() +
  annotate(
    "text",
    x = 0,
    y = 0,
    label = paste0(
      "Missing values :\n",
      nb_missing,
      " (",
      round(pct_missing, 2),
      "%)",
      "\n\nConsidered values :\n",
      nrow(pig_cond) - nb_missing
    ),
    size = 4
  ) +
  theme_void()


# ---------------------------------------------------------
# Affichage côte à côte
# ---------------------------------------------------------

p1 + p2 +
  plot_layout(widths = c(3, 1)) +
  plot_annotation(
    title = paste0(
      "Distribution of pigment variables on 2018, 2021, 2023 transect day data"
    ),
    subtitle = "Missing values are excluded from histograms"
  )
