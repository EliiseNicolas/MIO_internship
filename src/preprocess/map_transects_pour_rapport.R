datas <- readRDS("F:/data_elise/ds_NASC_pig_ftle_fod/ds_NASC_per_esu_all/new_ratio_ds_NASC_per_esu_pig_ftle_fod_2018_2021_2022_2023_transect_38kHz_mask9.rds")
str(datas)
## ---------------------------------------------------------------
## Plot des transects colorés par année (2018, 2021, 2022, 2023)
## avec fond de carte (côtes)
## ---------------------------------------------------------------

## Packages nécessaires
# install.packages(c("ggplot2", "sf", "rnaturalearth", "rnaturalearthdata", "dplyr"))

library(ggplot2)
library(sf)
library(rnaturalearth)
library(dplyr)
## -----------------------------------------------------------------
## 1. Préparation des données
## -----------------------------------------------------------------
## Adaptez "df" au nom réel de votre data.frame si besoin
# df <- read.csv("votre_fichier.csv")

df <- datas %>%
  mutate(
    annee = format(as.Date(time_nasc), "%Y"),
    annee = factor(annee, levels = c("2018", "2021", "2022", "2023"))
  ) %>%
  filter(!is.na(annee), !is.na(lon_nasc), !is.na(lat_nasc))


## -----------------------------------------------------------------
## 2. Fond de carte (côtes) - téléchargé via rnaturalearth
## -----------------------------------------------------------------
world <- ne_countries(scale = "medium", returnclass = "sf")

## Emprise spatiale des données (+ marge de 1°)
marge <- 1
xlim <- range(df$lon_nasc, na.rm = TRUE) + c(-marge, marge)
ylim <- range(df$lat_nasc, na.rm = TRUE) + c(-marge, marge)

## -----------------------------------------------------------------
## 2bis. Fond de bathymétrie - téléchargé via marmap (NOAA ETOPO1)
## -----------------------------------------------------------------
# install.packages("marmap")
library(marmap)

bathy <- getNOAA.bathy(
  lon1 = xlim[1], lon2 = xlim[2],
  lat1 = ylim[1], lat2 = ylim[2],
  resolution = 4,      # en minutes ; diminuer (ex: 1) pour plus de détail, plus lent
  keep = TRUE          # met en cache le fichier localement (dossier de travail)
)

## Conversion en data.frame pour ggplot (uniquement les profondeurs < 0 = mer)
bathy_df <- as.xyz(bathy)
names(bathy_df) <- c("lon", "lat", "depth")
bathy_df <- bathy_df %>% filter(depth <= 0)

## Palette de bleus pour la bathymétrie (plus profond = plus foncé)
pal_bathy <- colorRampPalette(c("#08306b", "#2171b5", "#6baed6", "#c6dbef", "#f7fbff"))

## -----------------------------------------------------------------
## 3. Carte des transects colorés par année
## -----------------------------------------------------------------
p <- ggplot() +
  # fond de bathymétrie
  geom_raster(data = bathy_df, aes(x = lon, y = lat, fill = depth)) +
  scale_fill_gradientn(
    colours = pal_bathy(100),
    name = "Profondeur (m)"
  ) +
  # isobathe côtière (0 m) pour bien marquer la côte
  geom_contour(
    data = bathy_df, aes(x = lon, y = lat, z = depth),
    breaks = 0, color = "grey30", linewidth = 0.3
  ) +
  # fond de carte (terres émergées)
  geom_sf(data = world, fill = "grey85", color = "grey40", linewidth = 0.2, inherit.aes = FALSE) +
  # transects (points reliés par le déplacement du bateau)
  geom_point(
    data = df,
    aes(x = lon_nasc, y = lat_nasc, color = annee),
    size = 0.6, alpha = 0.8
  ) +
  scale_color_manual(
    values = c(
      "2018" = "#ffdd00",
      "2021" = "#d95f02",
      "2022" = "#e6194B",
      "2023" = "#000000"
    ),
    name = "Année"
  ) +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  labs(
    title = "Transects du Marion Dufresnes II par année de campagne océanographique",
    x = "Longitude", y = "Latitude"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_line(color = "grey90", linewidth = 0.2),
    legend.position = "right"
  ) +
  # deux légendes de remplissage (profondeur) et couleur (année)
  guides(
    fill = guide_colorbar(order = 1),
    color = guide_legend(order = 2, override.aes = list(size = 2))
  )

print(p)
## Pour enregistrer la figure :
ggsave("C:/Users/mmolinet/elisou_ta_stagiaire_pref/MIO_internship_III/figures/figures_add_report/transects_par_annee_2.png", p, width = 8, height = 7, dpi = 300)

# ## -----------------------------------------------------------------
# ## 4. Variante : un panel par année (facet)
# ## -----------------------------------------------------------------
# p_facet <- ggplot() +
#   geom_sf(data = world, fill = "grey85", color = "grey40", linewidth = 0.2) +
#   geom_point(
#     data = df,
#     aes(x = lon_nasc, y = lat_nasc, color = annee),
#     size = 0.6, alpha = 0.8, show.legend = FALSE
#   ) +
#   scale_color_manual(
#     values = c(
#       "2018" = "#1b9e77",
#       "2021" = "#d95f02",
#       "2022" = "#7570b3",
#       "2023" = "#e7298a"
#     )
#   ) +
#   coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
#   facet_wrap(~ annee) +
#   labs(title = "Transects par campagne (facettes)", x = "Longitude", y = "Latitude") +
#   theme_minimal(base_size = 11) +
#   theme(panel.background = element_rect(fill = "aliceblue"))
# 
# print(p_facet)