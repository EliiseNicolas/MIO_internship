# # Description 
# 
# from pigment concentrations along 2018-2021-2023 obsaustral oceanographic campain (day only)
# 1 - we normalize igment concentration
# 2 - we do a PCA on pigment concentrations
# 3 - we clusterize PCA scores to get groups of phytoplanctonic communities
# 4 - we observe NASC along those distribs

# libraries
library(dplyr)
library(tidyr)
library(ggplot2)

# Global variables
path <- "~/Documents/stage_MIO/pt_III/data_preprocessed/pigmeann/transect/pigs_transect_day_2018_2021_2023/pig_mask9_1d_2018_2021_2023.rds"
pig_vars <- c(
  "Chla", "Per", "But", "Fuco", "Hex",
  "Allo", "Zea", "Chlb", "DvChla"
)

pig <- readRDS(path)
str(pig)

# 1 - pigment concentration normalization
# garder seulement les données ou il y a des pigments partout
data_pca <- pig %>%
  filter(
    if_all(
      all_of(pig_vars),
      is.finite
    )
  )
nrow(data_pca) # 403

# normalisation
data_pca <- data_pca %>%
  mutate(
    pig_sum = rowSums(
      across(all_of(pig_vars))
    )
  ) %>%
  filter(pig_sum > 0) %>%
  mutate(
    across(
      all_of(pig_vars),
      ~ .x / pig_sum
    )
  )

# 2 - PCA
pca <- prcomp(
  data_pca %>%
    select(all_of(pig_vars)),
  center = TRUE,
  scale. = TRUE
)

summary(pca)

library(factoextra)

fviz_eig(
  pca,
  addlabels = TRUE
)

fviz_pca_var(
  pca,
  col.var = "contrib",
  repel = TRUE
)

data_pca <- data_pca %>%
  mutate(
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    PC3 = pca$x[, 3]
  )

set.seed(123)

pca_scores <- data_pca %>%
  select(PC1, PC2, PC3)

# selection du nombre de cluster 
set.seed(123)

wss <- sapply(1:10, function(k) {
  kmeans(
    pca_scores,
    centers = k,
    nstart = 100
  )$tot.withinss
})

wss_df <- data.frame(
  k = 1:10,
  WSS = wss
)

ggplot(wss_df, aes(x = k, y = WSS)) +
  geom_line() +
  geom_point(size = 3) +
  labs(
    x = "Nombre de clusters",
    y = "Within-cluster sum of squares",
    title = "Méthode du coude"
  ) +
  theme_minimal()

# cluster 5

km <- kmeans(
  pca_scores,
  centers = 5,
  nstart = 100
)

data_pca$community <- factor(km$cluster)

data_long <- data_pca %>%
  select(community, all_of(pig_vars)) %>%
  pivot_longer(
    cols = all_of(pig_vars),
    names_to = "pigment",
    values_to = "proportion"
  )

ggplot(data_long, aes(x = community, y = proportion)) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.7
  ) +
  geom_jitter(
    width = 0.15,
    alpha = 0.5,
    size = 1.5
  ) +
  facet_wrap(
    ~ pigment,
    scales = "free_y"
  ) +
  labs(
    x = "Communauté",
    y = "Proportion du pigment",
    title = "Distribution des pigments au sein des communautés"
  ) +
  theme_minimal()

ggplot(data_long, aes(x = pigment, y = proportion)) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.7
  ) +
  geom_jitter(
    width = 0.15,
    alpha = 0.5,
    size = 1.5
  ) +
  facet_wrap(~ community, ncol = 1) +
  labs(
    x = "Pigment",
    y = "Proportion du pigment",
    title = "Distribution des pigments par communauté"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# coloration des distrib par PFT
data_long <- data_pca %>%
  select(community, all_of(pig_vars)) %>%
  pivot_longer(
    cols = all_of(pig_vars),
    names_to = "pigment",
    values_to = "proportion"
  ) %>%
  mutate(
    PFT = case_when(
      pigment %in% c("DvChla", "Zea") ~ "Picocyanobacteria",
      pigment %in% c("Allo", "Hex", "But") ~ "Flagellates",
      pigment == "Fuco" ~ "Diatoms",
      TRUE ~ "Other"
    )
  )
p_pigments <- ggplot(
  data_long,
  aes(
    x = pigment,
    y = proportion,
    fill = PFT
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.7
  ) +
  geom_jitter(
    aes(color = PFT),
    width = 0.15,
    alpha = 0.5,
    size = 1.5
  ) +
  facet_wrap(
    ~ community,
    ncol = 1,
    axes = "all_x"
  ) +
  scale_fill_manual(
    values = c(
      "Picocyanobacteria" = "#E69F00",
      "Flagellates" = "#56B4E9",
      "Diatoms" = "#009E73"
    )
  ) +
  scale_color_manual(
    values = c(
      "Picocyanobacteria" = "#E69F00",
      "Flagellates" = "#56B4E9",
      "Diatoms" = "#009E73"
    )
  ) +
  labs(
    x = "Pigment",
    y = "Proportion",
    fill = "Functional group",
    color = "Functional group",
    title = "Pigment distributions within phytoplankton communities"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

p_pigments

# p_pigments <- ggplot(
#   data_long,
#   aes(x = pigment, y = proportion)
# ) +
#   geom_boxplot(
#     outlier.shape = NA,
#     width = 0.7
#   ) +
#   geom_jitter(
#     width = 0.15,
#     alpha = 0.5,
#     size = 1.5
#   ) +
#   facet_wrap(
#     ~ community,
#     ncol = 1,
#     axes = "all_x"
#   ) +
#   labs(
#     x = "Pigment",
#     y = "Proportion",
#     title = "Distribution des pigments par communauté"
#   ) +
#   theme_minimal() +
#   theme(
#     axis.text.x = element_text(
#       angle = 45,
#       hjust = 1
#     )
#   )
# 
# p_pigments

ggsave(
  "~/Documents/stage_MIO/distribution_pigments_communautes_kmeans_k5_PFT.png",
  plot = p_pigments,
  width = 8,
  height = 12,
  units = "in",
  dpi = 300
)


# distribution spatiale des clusters 
install.packages("rnaturalearth")
library(sf)
library(rnaturalearth)
library(ggplot2)

# Fond de carte
world <- ne_countries(
  scale = "medium",
  returnclass = "sf"
)

# Carte des communautés
p_map <- ggplot() +
  
  # continents
  geom_sf(
    data = world,
    fill = "grey85",
    color = "grey50",
    linewidth = 0.2
  ) +
  
  # observations
  geom_point(
    data = data_pca,
    aes(
      x = lon_sv,
      y = lat_sv,
      color = community
    ),
    size = 3,
    alpha = 0.8
  ) +
  
  labs(
    x = "Longitude",
    y = "Latitude",
    color = "Communauté",
    title = "Distribution spatiale des communautés phytoplanctoniques"
  ) +
  
  coord_sf() +
  
  theme_minimal()

p_map

# regarder nasc dans les clusters
nasc_ds <- readRDS("~/Documents/stage_MIO/pt_III/data_preprocessed/NASC/transect_2018_2022_2023/NASC_mean_pig_grid_by_year_2018_2021_2023_day_200kHz.rds")
str(nasc_ds)

# joindre le NASC aux communautés
data_cluster <- data_pca %>%
  select(
    time,
    lat_sv,
    lon_sv,
    community
  ) %>%
  left_join(
    nasc_ds %>%
      select(
        time,
        lat,
        lon,
        NASC
      ),
    by = c(
      "time" = "time",
      "lat_sv" = "lat",
      "lon_sv" = "lon"
    )
  )

str(data_cluster)

# visualiser ditrib de NASC au sein des clusters
prop.table(table(data_cluster$community, data_cluster$NASC), margin = 1)

ggplot(data_cluster, aes(x = log(NASC))) +
  geom_histogram(bins = 30) +
  facet_wrap(~ community, ncol = 1) +
  labs(
    x = "log(NASC)",
    y = "Effectif",
    title="log(NASC) distribution per cluster of pigments concentration PCA", 
    subtitle = "2018, 2021, 2023 transect day dataset"
  ) +
  theme_minimal()
