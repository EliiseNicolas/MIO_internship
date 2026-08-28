# Description
# 
# A partir d'une résolution en km, on veut regridder la grille de pigment (en °) selon cette résolution
# On veut ensuite exporter en RDS les coordonnées lat/lon correspondant

path_pig <- pigments <- c("Chla", "Per", "But", "Fuco", "Hex", "Allo", "Zea", "Chlb", "DvChla")
path_pig <- "F:/data_elise/pigmeann/pigments_2018_2021_2022_2023_crop.rds"
pigs <- readRDS(path_pig)
str(pigs)
lat <- pigs$lat
lon <- pigs$lon
grids <- list(
  c(lon_km = 1500, lat_km = 1000),
  c(lon_km = 1000, lat_km = 700),
  c(lon_km = 500,  lat_km = 500),
  c(lon_km = 200,  lat_km = 200)
)

build_grid <- function(lat, lon, lon_res_km, lat_res_km, lat0 = NULL) {
  
  if (is.null(lat0)) lat0 <- mean(range(lat))
  
  km_per_deg_lon <- 111.32 * cos(lat0 * pi / 180)
  km_per_deg_lat <- 111.32
  
  lon_res_deg <- lon_res_km / km_per_deg_lon
  lat_res_deg <- lat_res_km / km_per_deg_lat
  
  lon_breaks <- seq(floor(min(lon) / lon_res_deg) * lon_res_deg, ceiling(max(lon) / lon_res_deg) * lon_res_deg, by = lon_res_deg)
  lat_breaks <- seq(floor(min(lat) / lat_res_deg) * lat_res_deg, ceiling(max(lat) / lat_res_deg) * lat_res_deg, by = lat_res_deg)
  
  list(
    lon_breaks = lon_breaks, lat_breaks = lat_breaks,
    lon_centers = head(lon_breaks, -1) + lon_res_deg / 2,
    lat_centers = head(lat_breaks, -1) + lat_res_deg / 2,
    lon_res_deg = lon_res_deg, lat_res_deg = lat_res_deg,
    lon_res_km = lon_res_km, lat_res_km = lat_res_km, lat0 = lat0
  )
}

# ============================================================
# Construction des 4 grilles + apercu rapide du nombre de cellules
# ============================================================

grids_list <- lapply(grids, function(res) build_grid(lat = lat, lon = lon, lon_res_km = res["lon_km"], lat_res_km = res["lat_km"]))
names(grids_list) <- sapply(grids, function(res) sprintf("lon%d_lat%d", res["lon_km"], res["lat_km"]))
print(grids_list)

grid_summary <- data.frame(
  grille = names(grids_list),
  lon_res_km = sapply(grids, function(r) r["lon_km"]),
  lat_res_km = sapply(grids, function(r) r["lat_km"]),
  n_cells_lon = sapply(grids_list, function(g) length(g$lon_centers)),
  n_cells_lat = sapply(grids_list, function(g) length(g$lat_centers))
)
grid_summary$n_cells_total <- grid_summary$n_cells_lon * grid_summary$n_cells_lat
print(grid_summary)

# ============================================================
# Sauvegarde des grilles personnalisees au format attendu par le
# pipeline NASC existant (liste $lon/$lat = vecteurs de coordonnees
# des centres de cellule, meme structure que le pig_grid original).
# ============================================================

grid_output_dir <- "F:/data_elise/sv_cropped/grids_custom"
dir.create(grid_output_dir, recursive = TRUE, showWarnings = FALSE)

for (nm in names(grids_list)) {
  g <- grids_list[[nm]]
  pig_grid_custom <- list(lon = g$lon_centers, lat = g$lat_centers)
  
  out_path <- file.path(grid_output_dir, paste0("pigmeann_grid_", nm, ".rds"))
  saveRDS(pig_grid_custom, out_path)
  cat("Sauvegarde :", out_path, "(", length(g$lon_centers), "x", length(g$lat_centers), "cellules )\n")
}