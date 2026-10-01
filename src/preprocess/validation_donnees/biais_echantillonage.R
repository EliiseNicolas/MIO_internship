library(sf)
library(dplyr)
library(ggplot2)
library(rnaturalearth)   # pour le trait de côte
library(patchwork)       # pour juxtaposer les cartes

options(scipen = 999) 

# Data validation : biais echantillonage 
path <- "F:/data_elise/ds_NASC_pig_ftle_fod/ds_NASC_per_esu_all/new_ratio_ds_NASC_per_esu_pig_ftle_fod_2018_2021_2022_2023_transect_38kHz_mask9.rds"

ds <- readRDS(path)
str(ds)
ds <- ds[ds$day == 3, ]
# Carte des points de NASC avec grilles de résolution (20x20km ou 1500 x 1500 km)


# 1. Points NASC valides -> sf
pts <- ds %>%
  filter(!is.na(nasc), !is.na(lat_nasc), !is.na(lon_nasc)) %>%
  st_as_sf(coords = c("lon_nasc", "lat_nasc"), crs = 4326)

# 2. Projection Lambert azimutale équivalente centrée sur la zone d'étude
lon0 <- mean(st_coordinates(pts)[, 1])
lat0 <- mean(st_coordinates(pts)[, 2])
crs_laea <- sprintf("+proj=laea +lat_0=%.3f +lon_0=%.3f +datum=WGS84 +units=m",
                    lat0, lon0)

xy <- st_coordinates(st_transform(pts, crs_laea))

# 3. Comptage par cellule pour une résolution donnée (en km)
grille_comptage <- function(xy, res_km) {
  res <- res_km * 1000
  data.frame(ix = floor(xy[, 1] / res),
             iy = floor(xy[, 2] / res)) %>%
    count(ix, iy, name = "n") %>%
    mutate(x = (ix + 0.5) * res,   # centre de la cellule
           y = (iy + 0.5) * res,
           res_km = res_km)
}

sf_use_s2(FALSE)

# Terre découpée sur la zone d'étude puis projetée (évite les géométries
# invalides quand on projette le monde entier en LAEA)
b <- st_bbox(pts)
land_zone <- land %>%
  st_crop(c(xmin = b[["xmin"]] - 20, ymin = b[["ymin"]] - 20,
            xmax = b[["xmax"]] + 20, ymax = b[["ymax"]] + 20)) %>%
  st_union() %>%
  st_transform(crs_laea) %>%
  st_make_valid()

# Grille complète (cellules vides incluses), cellules terrestres exclues
grille_complete <- function(xy, res_km, land_p = land_zone) {
  res <- res_km * 1000
  occ <- grille_comptage(xy, res_km)
  
  full <- expand.grid(ix = seq(min(occ$ix), max(occ$ix)),
                      iy = seq(min(occ$iy), max(occ$iy))) %>%
    left_join(select(occ, ix, iy, n), by = c("ix", "iy")) %>%
    mutate(n = coalesce(n, 0L),
           x = (ix + 0.5) * res,
           y = (iy + 0.5) * res,
           res_km = res_km)
  
  # une cellule est "à terre" si son centre est sur un continent/île
  centres <- st_as_sf(full, coords = c("x", "y"), crs = crs_laea)
  full$terre <- lengths(st_intersects(centres, land_p)) > 0
  filter(full, !terre)
}

carte_densite <- function(res_km, marge_km = 0) {
  g <- grille_complete(xy, res_km)
  pct_vide <- 100 * mean(g$n == 0)
  res <- res_km * 1000
  m <- marge_km * 1000
  
  # limites = bord extérieur des cellules de la grille (+ marge éventuelle)
  xlim <- range(g$x) + c(-1, 1) * (res / 2 + m)
  ylim <- range(g$y) + c(-1, 1) * (res / 2 + m)
  
  ggplot() +
    geom_tile(data = g, aes(x = x, y = y, fill = if_else(n == 0, NA, n)),
              width = res, height = res,
              colour = if (res_km >= 100) "grey70" else NA, linewidth = 0.1) +
    geom_sf(data = land_zone, fill = "grey80", colour = "grey50", linewidth = 0.2) +
    scale_fill_viridis_c(trans = "log10", na.value = "white",
                         name = "Nb points\nNASC") +
    coord_sf(crs = crs_laea, xlim = xlim, ylim = ylim, expand = FALSE) +
    labs(title = sprintf("Densité d'échantillonnage NASC - grille %g x %g km",
                         res_km, res_km),
         subtitle = sprintf("%.1f %% de cellules vides",
                            pct_vide),
         x = NULL, y = NULL) +
    theme_bw()
}

p <- carte_densite(20)

ggsave("C:/Users/mmolinet/elisou_ta_stagiaire_pref/MIO_internship_III/src/preprocess/validation_donnees/densite_nasc_20km.png", p, width = 20, height = 18, units = "cm", dpi = 300)

# Histogramme par campagne
library(dplyr)
library(ggplot2)

d <- ds %>%
  filter(!is.na(nasc), !is.na(time_nasc)) %>%
  mutate(annee = factor(format(time_nasc, "%Y")),
         mois  = factor(format(time_nasc, "%m"), levels = sprintf("%02d", 1:12),
                        labels = month.abb))

# 1. Nombre de données par année
p <- ggplot(d, aes(x = annee)) +
  geom_bar(fill = "steelblue") +
  scale_y_continuous(labels = scales::label_number(big.mark = " ")) +
  labs(x = "Année", y = "Nombre de points NASC", title = "Nombre de données de NASC par Année") +
  theme_bw() 
ggsave("C:/Users/mmolinet/elisou_ta_stagiaire_pref/MIO_internship_III/src/preprocess/validation_donnees/densite_nasc_annee.png", p, width = 20, height = 18, units = "cm", dpi = 300)

# 2. Nombre de données par mois, une facette par année
p <- ggplot(d, aes(x = mois)) +
  geom_bar(fill = "steelblue") +
  facet_wrap(~ annee) +
  scale_x_discrete(drop = TRUE) +
  labs(x = "Mois", y = "Nombre de points NASC", title = "Nombre de données de NASC par mois et Année") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("C:/Users/mmolinet/elisou_ta_stagiaire_pref/MIO_internship_III/src/preprocess/validation_donnees/densite_nasc_mois_et_annee.png", p, width = 20, height = 18, units = "cm", dpi = 300)
