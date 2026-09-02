# ============================================================
# Communautes phytoplanctoniques (PCA + k-means) a partir des
# ratios pigment/Chla du RDS unique deja aligne (all_ds), au lieu
# des donnees ESU/NASC ponctuelles du script d'origine.
#
# 1 - normalisation des ratios pigmentaires (proportions)
# 2 - PCA sur ces proportions
# 3 - k-means sur les scores PCA -> communautes phytoplanctoniques
# 4 - distributions des pigments par communaute (+ coloration PFT)
# 5 - distribution spatiale des communautes (carte raster complete)
#
# ADAPTATIONS PAR RAPPORT AU SCRIPT D'ORIGINE :
#  - Source : grille complete (pixel x date, ~777 600 pixels x
#    N dates) au lieu de quelques centaines de points ESU le long
#    de transects acoustiques -> N potentiellement ENORME. On
#    calcule donc les positions valides directement sur les arrays
#    (masque booleen + indexation matricielle), SANS jamais
#    materialiser un data.frame de la grille complete avant filtrage
#    (un expand.grid() + filter() classique serait bien trop gros).
#  - Noms des variables : ceux produits par le pipeline d'extraction
#    (suffixe "_chla" en minuscule + "chla_total"), pas les noms
#    "Xxx_Chla" du fichier NASC d'origine.
#  - PCA/k-means ajustes sur un SOUS-ECHANTILLON (sample_n_fit) pour
#    rester tractable, puis projection de TOUTE la grille valide sur
#    ce resultat (memes formules center/scale/rotation, puis
#    affectation au centroide k-means le plus proche).
#  - Pas de colonne NASC dans all_ds -> la partie "on observe le NASC
#    le long de ces distributions" mentionnee dans la description
#    d'origine n'est pas reproduite ici.
#  - Carte finale : raster complet (geom_raster) au lieu de points
#    ESU (geom_point) ; deux modes possibles (communaute dominante
#    par pixel sur toute la periode, ou carte d'une date precise).
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(factoextra)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
# ------------------------------------------------------------
# Options
# ------------------------------------------------------------

path_all_ds <- "F:/data_elise/prediction_ds/ds_ftle_pig_fod_ALL_DATES.rds"
n_cluster   <- 5
out_dir <- "C:/Users/mmolinet/elisou_ta_stagiaire_pref/MIO_internship_III/figures/commu_phyto_cluster_pigs"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
# La grille complete (pixel x date) peut representer des dizaines de
# millions de lignes valides -- bien plus que les 403 lignes du
# script d'origine. On ajuste la PCA/k-means sur un sous-echantillon
# pour rester tractable (prcomp scale bien, mais kmeans(nstart=100)
# sur des millions de lignes serait tres lent), puis on projette
# TOUTES les lignes valides dans cet espace pour la carte finale.
# Mets sample_n_fit <- Inf pour utiliser toutes les donnees (a ne
# faire que si ta RAM et ton temps de calcul le permettent).
sample_n_fit <- 200000
set.seed(123)

pig_vars <- c(
  "chla_total", "per_chla", "but_chla", "fuco_chla", "hex_chla",
  "allo_chla", "zea_chla", "chlb_chla", "dvchla_chla"
)

# ============================================================
# 0. Chargement + positions valides (tous ratios pigmentaires finis)
# ============================================================

all_ds <- readRDS(path_all_ds)

dates <- all_ds$date
lons  <- all_ds$lon
lats  <- all_ds$lat

# Masque booleen [n_date, n_lon, n_lat] : TRUE la ou TOUS les ratios
# pigmentaires sont renseignes. Calcule directement sur les arrays
# (economie de memoire majeure vs. construire d'abord la grille
# complete en data.frame puis filtrer).
valid_mask <- Reduce(`&`, lapply(pig_vars, function(v) is.finite(all_ds$pig[[v]])))
valid_idx  <- which(valid_mask, arr.ind = TRUE)   # matrice [n_valid, 3] : date_idx, lon_idx, lat_idx

data_pca <- data.frame(
  date = dates[valid_idx[, 1]],
  lon  = lons[valid_idx[, 2]],
  lat  = lats[valid_idx[, 3]]
)

for (v in pig_vars) {
  data_pca[[v]] <- all_ds$pig[[v]][valid_idx]
}

nrow(data_pca)

# ============================================================
# 1 - normalisation des ratios pigmentaires (proportions)
# ============================================================

data_pca <- data_pca %>%
  mutate(
    pig_sum = rowSums(across(all_of(pig_vars)))
  ) %>%
  filter(pig_sum > 0) %>%
  mutate(
    across(all_of(pig_vars), ~ .x / pig_sum)
  )

nrow(data_pca)

# ------------------------------------------------------------
# Sous-echantillon pour ajuster la PCA/k-means
# ------------------------------------------------------------

if (is.finite(sample_n_fit) && nrow(data_pca) > sample_n_fit) {
  fit_idx <- sample(nrow(data_pca), sample_n_fit)
} else {
  fit_idx <- seq_len(nrow(data_pca))
}
data_fit <- data_pca[fit_idx, ]

# ------------------------------------------------------------ 2 - PCA
pca <- prcomp(
  data_fit %>% select(all_of(pig_vars)),
  center = TRUE,
  scale. = TRUE
)

summary(pca)

fviz_eig(pca, addlabels = TRUE)

fviz_pca_var(pca, col.var = "contrib", repel = TRUE)

pc_scores_fit <- as.data.frame(pca$x)
data_fit <- bind_cols(data_fit, pc_scores_fit)

pca_scores <- data_fit %>%
  select(PC1, PC2, PC3, PC4, PC5, PC6, PC7)

# ------------------------------------------------------------
# selection du nombre de cluster pour le kmeans (elbow method)
# ------------------------------------------------------------

set.seed(123)

wss <- sapply(1:10, function(k) {
  kmeans(pca_scores, centers = k, nstart = 100)$tot.withinss
})

wss_df <- data.frame(k = 1:10, WSS = wss)

ggplot(wss_df, aes(x = k, y = WSS)) +
  geom_line() +
  geom_point(size = 3) +
  labs(
    x = "Nombre de clusters",
    y = "Within-cluster sum of squares",
    title = "Méthode du coude"
  ) +
  theme_minimal()

# cluster (n_cluster)

km <- kmeans(pca_scores, centers = n_cluster, nstart = 100)

data_fit$community <- factor(km$cluster)

# ============================================================
# 3 - Projection de TOUTE la grille valide dans l'espace PCA, puis
#     affectation au centroide k-means le plus proche (pas seulement
#     l'echantillon utilise pour ajuster la PCA/k-means)
# ============================================================

pig_mat_all   <- as.matrix(data_pca[, pig_vars])
pc_scores_all <- scale(pig_mat_all, center = pca$center, scale = pca$scale) %*% pca$rotation

assign_cluster <- function(pc_mat, centers) {
  # distance euclidienne a chaque centroide k-means, sur les memes
  # composantes que celles utilisees par kmeans() (ici PC1-PC7)
  pc_mat_sub <- pc_mat[, colnames(centers), drop = FALSE]
  d <- sapply(seq_len(nrow(centers)), function(k) {
    rowSums((pc_mat_sub - matrix(centers[k, ], nrow(pc_mat_sub), ncol(centers), byrow = TRUE))^2)
  })
  max.col(-d)   # indice de la distance minimale par ligne
}

data_pca$community <- factor(assign_cluster(pc_scores_all, km$centers))

# ============================================================
# 4 - Distributions des pigments par communaute (sur l'echantillon
#     d'ajustement, comme dans le script d'origine)
# ============================================================

data_long <- data_fit %>%
  select(community, all_of(pig_vars)) %>%
  pivot_longer(
    cols = all_of(pig_vars),
    names_to = "pigment",
    values_to = "proportion"
  )

# Palette commune aux deux graphiques
community_colors <- c(
  "1" = "#1b9e77",
  "2" = "#d95f02",
  "3" = "#7570b3",
  "4" = "#e7298a",
  "5" = "#66a61e"
)

p <- ggplot(data_long, aes(x = community, y = proportion, color = community)) +
  stat_summary(fun = mean, geom = "point", size = 2.5) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
               geom = "errorbar", width = 0.3, linewidth = 0.6) +
  facet_wrap(~ pigment, scales = "free_y") +
  scale_color_manual(values = community_colors, name = "Communauté") +
  labs(
    x = "Communauté",
    y = "Proportion du pigment (moyenne ± 1 écart-type)",
    title = "Distribution des pigments au sein des communautés"
  ) +
  theme_minimal()

ggsave(
  file.path(out_dir, "distrib_pigments_in_communautes_kmeans_k5_PFT.png"),
  plot = p,
  width = 14,
  height = 8,
  units = "in",
  dpi = 300
)

# ---- coloration des distrib par Phytoplankton Functional Type (PFT) ----

data_long <- data_fit %>%
  select(community, all_of(pig_vars)) %>%
  pivot_longer(
    cols = all_of(pig_vars),
    names_to = "pigment",
    values_to = "proportion"
  ) %>%
  mutate(
    PFT = case_when(
      pigment %in% c("dvchla_chla", "zea_chla")          ~ "Picocyanobacteria",
      pigment %in% c("allo_chla", "hex_chla", "but_chla") ~ "Flagellates",
      pigment == "fuco_chla"                              ~ "Diatoms",
      TRUE                                                 ~ "Other"
    )
  )

p_pigments <- ggplot(
  data_long,
  aes(x = pigment, y = proportion, color = PFT)
) +
  stat_summary(fun = mean, geom = "point", size = 2.5,
               position = position_dodge(width = 0.3)) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
               geom = "errorbar", width = 0.3, linewidth = 0.6,
               position = position_dodge(width = 0.3)) +
  facet_wrap(~ community, ncol = 1, axes = "all_x") +
  scale_color_manual(
    values = c(
      "Picocyanobacteria" = "#E69F00",
      "Flagellates"       = "#56B4E9",
      "Diatoms"           = "#009E73",
      "Other"             = "grey40"
    )
  ) +
  labs(
    x = "Pigment",
    y = "Proportion (moyenne ± 1 écart-type)",
    color = "Functional group",
    title = "Pigment distributions within phytoplankton communities"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_pigments

ggsave(
  file.path(out_dir, "distribution_pigments_communautes_kmeans_k5_PFT.png"),
  plot = p_pigments,
  width = 8,
  height = 12,
  units = "in",
  dpi = 300
)

# ============================================================
# 5 - Distribution spatiale des communautes (grille complete)
# ============================================================
# Deux modes possibles :
#  - "dominant"    : communaute la plus frequente par pixel sur
#                    toute la periode (vision "province" stable)
#  - "single_date" : carte d'une seule date precise
# La grille etant desormais complete (et non des points ESU
# eparses), on utilise geom_raster plutot que geom_point.

map_mode <- "single_date"
# map_date <- dates[1]   # utilise seulement si map_mode == "single_date"
for (map_date in dates) {
  map_date <- as.Date(map_date)
  
  if (map_mode == "dominant") {
    map_data <- data_pca %>%
      dplyr::count(lon, lat, community) %>%
      dplyr::group_by(lon, lat) %>%
      dplyr::slice_max(n, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup()
    map_title <- "Communauté phytoplanctonique dominante par pixel (toutes dates confondues)"
    date_str <- "all_dates"
  } else {
    map_data <- data_pca %>% dplyr::filter(.data$date == map_date)
    map_title <- paste("Communauté phytoplanctonique -", format(map_date, "%Y-%m-%d"))
    date_str <- format(map_date, "%Y-%m-%d")
  }
  
  world <- ne_countries(scale = "medium", returnclass = "sf")
  
  p_map <- ggplot() +
    geom_raster(data = map_data, aes(x = lon, y = lat, fill = community)) +
    geom_sf(data = world, fill = "grey85", color = "grey50", linewidth = 0.2) +
    coord_sf(xlim = range(lons), ylim = range(lats)) +
    scale_fill_manual(values = community_colors, name = "Communauté") +
    labs(x = "Longitude", y = "Latitude", title = map_title) +
    theme_minimal()
  
  print(p_map)
  
  ggsave(
    file.path(out_dir, paste0(date_str, "_map_pigments_communautes_kmeans_k5_PFT.png")),
    plot = p_map, width = 8, height = 8, units = "in", dpi = 300
  )
}



# p_map_purity <- ggplot() +
#   geom_raster(data = map_data, aes(x = lon, y = lat, fill = purity)) +
#   geom_sf(data = world, fill = "grey85", color = "grey50", linewidth = 0.2) +
#   coord_sf(xlim = range(lons), ylim = range(lats)) +
#   scale_fill_viridis_c(limits = c(1 / n_cluster, 1), name = "Pureté") +
#   labs(
#     x = "Longitude",
#     y = "Latitude",
#     title = "Pureté de la communauté dominante (proportion de jours ou elle est observee)"
#   ) +
#   theme_minimal()
# 
# p_map_purity