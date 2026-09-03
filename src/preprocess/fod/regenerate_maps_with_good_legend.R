save_path <- "F:/data_elise/fod_elise_2018_2021_2022_2023/mclust_results_2018_2021_2022_2023/transitions_upgraded"

# ------------------------------------------------------------
# Chargement des variables nécessaires
# ------------------------------------------------------------

lon      <- readRDS(file.path(save_path, "lon.rds"))
lat      <- readRDS(file.path(save_path, "lat.rds"))
time_bis <- readRDS(file.path(save_path, "time.rds"))

cluster_transition_map_new <- readRDS(file.path(save_path, "cluster_transition_map_renamed.rds"))

stopifnot(dim(cluster_transition_map_new)[1] == length(lon))
stopifnot(dim(cluster_transition_map_new)[2] == length(lat))
stopifnot(dim(cluster_transition_map_new)[3] == length(time_bis))

out_dir <- file.path(save_path, "maps_regenerated_1")
dir.create(out_dir, showWarnings = FALSE)

# ============================================================
# COULEURS DES CLUSTERS ET DES TRANSITIONS
# ============================================================

cluster_cols <- c(
  "1" = "#2166AC", "2" = "#67A9CF", "3" = "#1A9850",
  "4" = "#A6D96A", "5" = "#FDAE61", "6" = "#D73027"
)

transition_cols <- c(
  "7"  = "#3B73B9", "8"  = "#3FA7B5", "9"  = "#45B97C",
  "10" = "#C46A00", "11" = "#E85D04", "12" = "#C9184A",
  "13" = "#8F1D3F"
)

fod_cols <- c(cluster_cols, transition_cols)

legend_codes  <- c(1, 7, 2, 8, 9, 3, 10, 4, 11, 12, 5, 13, 6)
legend_labels <- c(
  "C1", "T1-2", "C2", "T1-3", "T2-3", "C3", "T3-4",
  "C4", "T4-6", "T4-5", "C5", "T5-6", "C6"
)

# Correspondance code numerique -> libelle (manquait dans la version precedente)
fod_label_map  <- setNames(legend_labels, as.character(legend_codes))
fod_cols_named <- setNames(fod_cols[as.character(legend_codes)], legend_labels)

# ------------------------------------------------------------
# Fonctions de formatage des axes en degres
# ------------------------------------------------------------

label_lon <- function(x) paste0(abs(x), ifelse(x >= 0, "\u00b0E", "\u00b0W"))
label_lat <- function(x) paste0(abs(x), ifelse(x >= 0, "\u00b0N", "\u00b0S"))

# ------------------------------------------------------------
# Boucle sur les dates : une carte ggplot carree par date
# ------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

world <- ne_countries(scale = "medium", returnclass = "sf")

for (i in seq_along(time_bis)) {
  
  map_slice <- cluster_transition_map_new[, , i]   # [lon, lat]
  
  map_df <- expand.grid(lon = lon, lat = lat) %>%
    mutate(
      code    = as.vector(map_slice),
      fod_lbl = factor(fod_label_map[as.character(code)], levels = legend_labels)
    )
  
  p <- ggplot() +
    geom_raster(data = map_df, aes(x = lon, y = lat, fill = fod_lbl)) +
    geom_sf(data = world, fill = "grey85", color = "grey50", linewidth = 0.2) +
    coord_sf(xlim = range(lon), ylim = range(lat)) +
    scale_fill_manual(
      values = fod_cols_named, name = "FOD",
      na.value = "white", drop = FALSE,
      guide = guide_legend(reverse = TRUE)
    ) +
    labs(
      x = "Longitude", y = "Latitude",
      title = paste("Functional Oceanographic Domains (clusters) -",
                    format(time_bis[i], "%Y-%m-%d"))
    ) +
    theme_minimal()
  
  print(p)
  
  ggsave(
    filename = sprintf(file.path(out_dir, "cluster_transition_map_renamed_%s.png"),
                       format(time_bis[i], "%Y%m%d")),
    plot = p, width = 8, height = 8, units = "in", dpi = 300
  )
}

cat("Cartes régénérées dans :", out_dir, "\n")

# ============================================================
# TABLEAU DE VERIFICATION (compte les 13 classes + NA)
# ============================================================

tab <- table(cluster_transition_map_new, useNA = "ifany")
print(tab)
cat("Nombre de classes non-NA :", length(tab) - any(is.na(names(tab))), "\n")