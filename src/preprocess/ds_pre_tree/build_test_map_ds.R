# ============================================================
# Extraction FTLE, pigments et FOD pour une date précise
# Grille de référence = grille la plus fine (pigments, 1080x720)
# ============================================================

rm(list = ls())

# ------------------------------------------------------------
# Date cible
# ------------------------------------------------------------

date_cible <- as.Date("2023-01-26")

# ------------------------------------------------------------
# Paths -- ATTENTION : chemins NATIFS (haute résolution),
# pas les fichiers "grids_custom" agrégés à 200km
# ------------------------------------------------------------

path_ftle <- "F:/data_elise/ftle/ftle_2018_2021_2022_2023_cropped.rds"   # 901 x 600
path_pig  <- "F:/data_elise/pigmeann/pigments_2018_2021_2022_2023_crop.rds"  # 1080 x 720

path_fod_clusters <- "F:/data_elise/fod_elise_2018_2021_2022_2023/mclust_results_2018_2021_2022_2023/transitions_upgraded/cluster_transition_map_renamed.rds"
path_fod_lon  <- "F:/data_elise/fod_elise_2018_2021_2022_2023/mclust_results_2018_2021_2022_2023/lon.rds"
path_fod_lat  <- "F:/data_elise/fod_elise_2018_2021_2022_2023/mclust_results_2018_2021_2022_2023/lat.rds"
path_fod_time <- "F:/data_elise/fod_elise_2018_2021_2022_2023/mclust_results_2018_2021_2022_2023/time.rds"

# ------------------------------------------------------------
# Chargement
# ------------------------------------------------------------

ftle_grid <- readRDS(path_ftle)
pigs_grid <- readRDS(path_pig)

fod_clusters <- readRDS(path_fod_clusters)
fod_lon  <- readRDS(path_fod_lon)
fod_lat  <- readRDS(path_fod_lat)
fod_time <- readRDS(path_fod_time)
fod_date <- as.Date(fod_time)

list_pigs <- c("Chla", "Per", "But", "Fuco", "Hex", "Allo", "Zea", "Chlb", "DvChla")

# ------------------------------------------------------------
# Grille de référence = grille PIGMENTS (la plus fine : 1080x720)
# ------------------------------------------------------------

grid_lon <- pigs_grid$lon   # 1080
grid_lat <- pigs_grid$lat   # 720

# ------------------------------------------------------------
# Indices de reprojection (nearest neighbor -- duplique les valeurs
# des grilles plus grossières sur la grille de référence)
# ------------------------------------------------------------

idx_lon_ftle <- sapply(grid_lon, function(x) which.min(abs(ftle_grid$lon - x)))
idx_lat_ftle <- sapply(grid_lat, function(x) which.min(abs(ftle_grid$lat - x)))

idx_lon_fod <- sapply(grid_lon, function(x) which.min(abs(fod_lon - x)))
idx_lat_fod <- sapply(grid_lat, function(x) which.min(abs(fod_lat - x)))

# ------------------------------------------------------------
# FTLE : extraire le jour, puis reprojeter sur la grille pigments
# ------------------------------------------------------------

idx_date_ftle <- which(ftle_grid$date == date_cible)
if (length(idx_date_ftle) == 0) stop("Date absente du dataset FTLE")

ftle_day_native <- ftle_grid$ftle[idx_date_ftle, , ]        # matrix [901, 600]
ftle_day <- ftle_day_native[idx_lon_ftle, idx_lat_ftle]      # matrix [1080, 720] (dupliqué)

# ------------------------------------------------------------
# PIGMENTS : déjà sur la grille de référence, pas de reprojection
# ------------------------------------------------------------

idx_date_pig <- which(pigs_grid$date == date_cible)
if (length(idx_date_pig) == 0) stop("Date absente du dataset pigments")

pig_day <- list()
for (p in list_pigs) {
  pig_day[[p]] <- pigs_grid[[paste0("c_cond_", p)]][idx_date_pig, , ]  # matrix [1080, 720]
}

# ------------------------------------------------------------
# FOD : extraire le jour, puis reprojeter sur la grille pigments
# ------------------------------------------------------------

idx_date_fod <- which(fod_date == date_cible)
if (length(idx_date_fod) == 0) stop("Date absente du dataset FOD")

fod_day_native <- fod_clusters[, , idx_date_fod]     # matrix [541, 361]
fod_day <- fod_day_native[idx_lon_fod, idx_lat_fod]   # matrix [1080, 720] (dupliqué)

# ------------------------------------------------------------
# Assemblage final -- tout est déjà sur la grille de référence
# ------------------------------------------------------------

day_ds <- list(
  date = date_cible,
  lon  = grid_lon,
  lat  = grid_lat,
  ftle = ftle_day,
  pig  = pig_day,
  fod  = fod_day
)

str(day_ds)

# ------------------------------------------------------------
# Sauvegarde
# ------------------------------------------------------------

saveRDS(
  day_ds,
  paste0("F:/data_elise/ds_day_ftle_pig_fod/ds_ftle_pig_fod_", format(date_cible, "%Y%m%d"), ".rds")
)

# ============================================================
# Plot d'une variable pigment sur carte, pour une date donnée
# ============================================================


# ------------------------------------------------------------
# Chargement
# ------------------------------------------------------------

pigs_grid <- readRDS(path_pig)

idx_date <- which(pigs_grid$date == date_cible)
pigment    <- "Chla"
if (length(idx_date) == 0) stop("Date absente du dataset pigments")

# ------------------------------------------------------------
# Extraction du pigment pour cette date
# ------------------------------------------------------------

pig_day <- pigs_grid[[paste0("c_cond_", pigment)]][idx_date, , ]  # matrix [lon, lat]

# ------------------------------------------------------------
# Data.frame pour ggplot
# ------------------------------------------------------------

df_plot <- expand.grid(lon = pigs_grid$lon, lat = pigs_grid$lat)
df_plot$value <- as.vector(pig_day)

pct_valid <- 100 * sum(!is.na(df_plot$value)) / nrow(df_plot)

# ------------------------------------------------------------
# Carte
# ------------------------------------------------------------

ggplot(df_plot, aes(x = lon, y = lat, fill = value)) +
  geom_raster() +
  scale_fill_viridis_c(na.value = "grey90") +
  coord_quickmap() +
  theme_bw() +
  labs(
    title = paste(pigment, "-", format(date_cible, "%Y-%m-%d")),
    subtitle = paste0(round(pct_valid, 1), "% de pixels valides"),
    x = "Longitude", y = "Latitude",
    fill = pigment
  )