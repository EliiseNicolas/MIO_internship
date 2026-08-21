# Description 

# On veut étudier les tendances des pigments sur les années 2018 2021 2023 au sein des FOD
# On veut utiliser toutes les données pigmentaires au sein d'un cluster, pas seulement les points co-localisés avec le nasc
# 
# On veut avoir les plots suivants : 
# - Proportion des pigments pour chaque cluster pour chaque année 
# - évolution des moyenne de concentration des pigments pour chaque année et chaque fod
# - évolution de la proportion des pigments a la chlorophylle totale pour chaque année et chaque cluster de fod

# libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(emmeans)
library(multcomp)
library(multcompView)
library(purrr)

# Global variables
path_pig <- "F:/data_elise/pigmeann/pigments_2018_2021_2022_2023_crop.rds"
path_fod <- "F:/data_elise/fod_elise_2018_2021_2022_2023/mclust_results_2018_2021_2022_2023/cluster_map.rds"
path_time <- "F:/data_elise/fod_elise_2018_2021_2022_2023/mclust_results_2018_2021_2022_2023/time.rds"
path_lat <- "F:/data_elise/fod_elise_2018_2021_2022_2023/mclust_results_2018_2021_2022_2023/lat.rds"
path_lon <- "F:/data_elise/fod_elise_2018_2021_2022_2023/mclust_results_2018_2021_2022_2023/lon.rds"

pigments <- c(
  "Fuco",
  "Perid",
  "19but",
  "19hex",
  "Allo",
  "Zeax",
  "Lut",
  "Chlb",
  "Chla"
)

# ouvrir les ds 
pig <- readRDS(path_pig)
fod <- readRDS(path_fod)
time <- readRDS(path_time); lat <- readRDS(path_lat); lon <- readRDS(path_lon)

# Statistiques de pigments dans chaque fod

# garder uniquement les dates communes 
fod_dates <- as.Date(time)
pig_dates <- as.Date(pig$date)

dates_communes <- intersect(fod_dates, pig_dates)

length(fod_dates)
length(pig_dates)
length(dates_communes) # OK 196 dates

idx_fod <- which(fod_dates %in% dates_communes)
idx_pig <- which(pig_dates %in% dates_communes)

library(dplyr)
library(tidyr)
library(ggplot2)

# ============================================================
# 1. Correspondance spatiale : grille pigment -> grille FOD
# ============================================================

# Pour chaque longitude pigmentaire,
# trouver la longitude FOD la plus proche
idx_lon <- sapply(pig$lon, function(x) {
  which.min(abs(lon - x))
})

# Pour chaque latitude pigmentaire,
# trouver la latitude FOD la plus proche
idx_lat <- sapply(pig$lat, function(x) {
  which.min(abs(lat - x))
})


# ============================================================
# 2. Créer le tableau pigment + FOD
# ============================================================

pig_fod_list <- vector("list", length(dates_communes))

for (i in seq_along(dates_communes)) {
  
  # indices correspondants
  i_pig <- date_match$i_pig[i]
  i_fod <- date_match$i_fod[i]
  
  # ----------------------------------------------------------
  # Carte FOD de la date
  # ----------------------------------------------------------
  
  fod_i <- fod[, , i_fod]
  
  # Adapter la carte FOD à la grille pigmentaire
  # nearest neighbour
  fod_pig <- fod_i[idx_lon, idx_lat]
  
  
  # ----------------------------------------------------------
  # Créer le dataframe
  # ----------------------------------------------------------
  
  pig_fod_list[[i]] <- data.frame(
    
    date = pig$date[i_pig],
    
    lon = rep(pig$lon, times = length(pig$lat)),
    lat = rep(pig$lat, each = length(pig$lon)),
    
    Chla = as.vector(pig$c_cond_Chla[i_pig, , ]),
    Perid = as.vector(pig$c_cond_Per[i_pig, , ]),
    But = as.vector(pig$c_cond_But[i_pig, , ]),
    Fuco = as.vector(pig$c_cond_Fuco[i_pig, , ]),
    Hex = as.vector(pig$c_cond_Hex[i_pig, , ]),
    Allo = as.vector(pig$c_cond_Allo[i_pig, , ]),
    Zeax = as.vector(pig$c_cond_Zea[i_pig, , ]),
    Chlb = as.vector(pig$c_cond_Chlb[i_pig, , ]),
    DvChla = as.vector(pig$c_cond_DvChla[i_pig, , ]),
    
    # FOD associé
    fod = as.vector(fod_pig)
  )
}

# Combiner toutes les dates
pig_fod <- bind_rows(pig_fod_list)

str(pig_fod)
pig_fod <- pig_fod %>%
  filter(
    !is.na(Chla),
    !is.na(fod)
  ) %>%
  mutate(
    year = format(date, "%Y"),
    fod = factor(fod)
  )

table(pig_fod$year, pig_fod$fod)

pig_long <- pig_fod %>%
  pivot_longer(
    cols = c(
      Chla,
      Perid,
      But,
      Fuco,
      Hex,
      Allo,
      Zeax,
      Chlb,
      DvChla
    ),
    names_to = "pigment",
    values_to = "concentration"
  )

ggplot(
  pig_long,
  aes(x = fod, y = concentration)
) +
  geom_boxplot(
    na.rm = TRUE
  ) +
  facet_grid(
    pigment ~ year,
    scales = "free_y"
  ) +
  labs(
    x = "FOD cluster",
    y = "Pigment concentration"
  ) +
  theme_classic()
