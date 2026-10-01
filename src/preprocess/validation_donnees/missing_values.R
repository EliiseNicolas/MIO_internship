library(sf)
library(dplyr)
library(ggplot2)
library(rnaturalearth)   # pour le trait de côte
library(patchwork)       # pour juxtaposer les cartes

options(scipen = 999) 

# Data validation : données manquantes
path <- "F:/data_elise/ds_NASC_pig_ftle_fod/ds_NASC_per_esu_all/new_ratio_ds_NASC_per_esu_pig_ftle_fod_2018_2021_2022_2023_transect_38kHz_mask9.rds"

ds <- readRDS(path)
str(ds)
ds <- ds[ds$day == 3, ]

# Cartes de valeurs manquantes par région 
library(tidyr)
library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)
library(rnaturalearth)

sf_use_s2(FALSE)

# Données avec coordonnées NASC valides
d <- ds %>%
  filter(!is.na(lat_nasc), !is.na(lon_nasc)) %>%
  mutate(fod = na_if(fod, "NA"))

pts <- st_as_sf(d, coords = c("lon_nasc", "lat_nasc"), crs = 4326)

# Projection équivalente centrée sur la zone
lon0 <- mean(d$lon_nasc)
lat0 <- mean(d$lat_nasc)
crs_laea <- sprintf("+proj=laea +lat_0=%.3f +lon_0=%.3f +datum=WGS84 +units=m",
                    lat0, lon0)

xy_all <- st_coordinates(st_transform(pts, crs_laea))

# Trait de côte découpé sur la zone puis projeté
land <- ne_countries(scale = "medium", returnclass = "sf")
b <- st_bbox(pts)
land_zone <- land %>%
  st_crop(c(xmin = b[["xmin"]] - 20, ymin = b[["ymin"]] - 20,
            xmax = b[["xmax"]] + 20, ymax = b[["ymax"]] + 20)) %>%
  st_union() %>%
  st_transform(crs_laea) %>%
  st_make_valid()


vars <- c("nasc", "fod", "Chla", "ftle")

d <- ds %>%
  filter(!is.na(lat_nasc), !is.na(lon_nasc)) %>%
  mutate(fod = na_if(fod, "NA"))

xy_all <- d %>%
  st_as_sf(coords = c("lon_nasc", "lat_nasc"), crs = 4326) %>%
  st_transform(crs_laea) %>%
  st_coordinates()

# Comptage des manquants par cellule et par variable
grille_na <- function(res_km) {
  res <- res_km * 1000
  d %>%
    transmute(ix = floor(xy_all[, 1] / res),
              iy = floor(xy_all[, 2] / res),
              across(all_of(vars), is.na)) %>%          # is.na() détecte aussi NaN
    pivot_longer(all_of(vars), names_to = "variable", values_to = "manquant") %>%
    group_by(variable, ix, iy) %>%
    summarise(n_tot = n(), n_na = sum(manquant), .groups = "drop") %>%
    mutate(pct_na = 100 * n_na / n_tot,
           x = (ix + 0.5) * res,
           y = (iy + 0.5) * res)
}

carte_na <- function(res_km, marge_km = 0) {
  g <- grille_na(res_km)
  res <- res_km * 1000
  m <- marge_km * 1000
  
  # % global de manquants dans le titre de chaque facette
  etiq <- g %>%
    group_by(variable) %>%
    summarise(p = 100 * sum(n_na) / sum(n_tot)) %>%
    mutate(lab = sprintf("%s (%.1f %% manquants)", variable, p))
  g$variable <- factor(g$variable, levels = vars,
                       labels = etiq$lab[match(vars, etiq$variable)])
  
  ggplot() +
    geom_tile(data = g, aes(x = x, y = y, fill = if_else(n_na == 0, NA, n_na)),
              width = res, height = res) +
    geom_sf(data = land_zone, fill = "grey80", colour = "grey50", linewidth = 0.2) +
    scale_fill_viridis_c(option = "magma", trans = "log10", na.value = "grey90",
                         name = "Nb valeurs\nmanquantes",
                         labels = scales::label_number(big.mark = " ")) +
    coord_sf(crs = crs_laea,
             xlim = range(g$x) + c(-1, 1) * (res / 2 + m),
             ylim = range(g$y) + c(-1, 1) * (res / 2 + m),
             expand = FALSE) +
    facet_wrap(~ variable) +
    labs(title = sprintf("Valeurs manquantes par cellule - grille %g x %g km",
                         res_km, res_km),
         x = NULL, y = NULL) +
    theme_bw()
}
cols <- scales::hue_pal()(4)   # couleurs par défaut, dans l'ordre de vars
couleurs <- c(nasc = cols[1], fod = cols[3], Chla = cols[2], ftle = cols[4])
p_na <- carte_na(20)

ggsave("C:/Users/mmolinet/elisou_ta_stagiaire_pref/MIO_internship_III/src/preprocess/validation_donnees/na_20km.png", p_na, width = 25, height = 22, units = "cm", dpi = 300)

# Histogramme de valeurs manquantes par campagne/mois
# 1. Identification des campagnes
na_long <- d %>%
  filter(!is.na(time_nasc)) %>%
  transmute(annee   = format(time_nasc, "%Y"),
            mois_num = as.integer(format(time_nasc, "%m")),
            across(all_of(vars), is.na)) %>%
  pivot_longer(all_of(vars), names_to = "variable", values_to = "manquant") %>%
  mutate(variable = factor(variable, levels = vars))

# Mois présents dans l'ensemble du jeu de données (les 3 mois d'échantillonnage)
mois_presents <- sort(unique(na_long$mois_num))
na_long <- na_long %>%
  mutate(mois = factor(month.abb[mois_num], levels = month.abb[mois_presents]))

# 1. Nombre de valeurs manquantes par année (tous mois confondus)
na_an <- na_long %>%
  group_by(annee, variable) %>%
  summarise(n_na = sum(manquant), .groups = "drop")

p_an <- ggplot(na_an, aes(x = annee, y = n_na, fill = variable)) +
  geom_col(position = "dodge") +
  scale_y_continuous(labels = scales::label_number(big.mark = " ")) +
  labs(x = "Année", y = "Nombre de valeurs manquantes", fill = "Variable") +
  theme_bw() + scale_fill_manual(values = couleurs)
p_an
ggsave("C:/Users/mmolinet/elisou_ta_stagiaire_pref/MIO_internship_III/src/preprocess/validation_donnees/hist_miss_vals_annee.png", p_an, width = 25, height = 22, units = "cm", dpi = 300)

# 2. Par année et par mois : les 3 mois toujours sur l'axe x
na_an_mois <- na_long %>%
  group_by(annee, mois, variable) %>%
  summarise(n_na = sum(manquant), .groups = "drop")

p_an_mois <- ggplot(na_an_mois, aes(x = mois, y = n_na, fill = variable)) +
  geom_col(position = position_dodge(preserve = "single")) +
  facet_wrap(~ annee) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(labels = scales::label_number(big.mark = " ")) +
  labs(x = "Mois", y = "Nombre de valeurs manquantes", fill = "Variable") +
  theme_bw() + scale_fill_manual(values = couleurs)
p_an_mois
ggsave("C:/Users/mmolinet/elisou_ta_stagiaire_pref/MIO_internship_III/src/preprocess/validation_donnees/hist_miss_vals_mois_annee.png", p_an_mois, width = 25, height = 22, units = "cm", dpi = 300)


# variables manquantes conjointement 
install.packages(c("naniar", "UpSetR"))
library(naniar)

ds_na <- ds %>% mutate(fod = na_if(fod, "NA"))

gg_miss_upset(ds_na,
              nsets = n_var_miss(ds_na),   # toutes les variables ayant au moins un NA
              nintersects = 40,            # nb max de combinaisons affichées
              order.by = "freq")
png("C:/Users/mmolinet/elisou_ta_stagiaire_pref/MIO_internship_III/src/preprocess/validation_donnees/upset_na.png", width = 30, height = 18, units = "cm", res = 300)
gg_miss_upset(ds_na, nsets = n_var_miss(ds_na), nintersects = 40)
dev.off()