# =====================================================================
# 00_config.R -- Configuration centrale du pipeline NASC
# =====================================================================
# Tous les autres scripts commencent par : source("R/00_config.R")
# Modifier ICI les chemins, fréquences, variables et grilles à tester,
# rien d'autre ne devrait avoir besoin d'être changé ailleurs.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(FNN)
  library(ranger)
  library(xgboost)
  library(rpart)
  library(rpart.plot)
  library(patchwork)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- Fréquences et période diurne ------------------------------------------
FREQS          <- c(38, 120)   # étendre à c(38, 70, 120, 200) si besoin
DIURNAL_PERIOD <- 3            # 3 = jour, 1 = nuit
DP_LABEL       <- "day"
N_CV           <- 10           # nb de folds pour l'entraînement final (naive ET blocked)

# ---- Chemins ----------------------------------------------------------------
# A adapter à l'environnement d'exécution (ici recopié depuis les scripts
# d'origine -- chemins Windows locaux).
PATH_TEMPLATE <- function(freq) {
  paste0(
    "F:/data_elise/ds_NASC_pig_ftle_fod/ds_NASC_per_esu_all/",
    "ds_NASC_per_esu_pig_ftle_fod_2018_2021_2022_2023_transect_",
    freq, "kHz_mask9.rds"
  )
}
PATH_GRID_DAY <- "F:/data_elise/ds_day_ftle_pig_fod/ds_ftle_pig_fod_20230126.rds"

OUTPUT_ROOT <- "outputs_pipeline"
dir.create(OUTPUT_ROOT, showWarnings = FALSE, recursive = TRUE)

path_out <- function(...) file.path(OUTPUT_ROOT, ...)

# ---- Variables ---------------------------------------------------------------
VARS_NUM <- c(
  "NASC", "per_ratio_chla", "but_ratio_chla", "fuco_ratio_chla",
  "hex_ratio_chla", "allo_ratio_chla", "zea_ratio_chla",
  "chlb_ratio_chla", "total_chla", "ftle"
)

COVARIATES_NUM <- setdiff(VARS_NUM, "NASC")
COVARIATES_ALL <- c(COVARIATES_NUM, "fod")
RESPONSE_VAR   <- "NASC"

# ---- Schémas de blocage spatio-temporel à tester ------------------------------
# Un buffer proportionnel à la taille de cellule est utilisé par défaut :
# ADAPTER ces valeurs si le variogramme des résidus indique une portée
# d'autocorrélation différente (cf. section "métadonnées / diagnostics").
SPATIAL_RESOLUTIONS <- list(
  list(label = "1500x1000km", lon_km = 1500, lat_km = 1000, buffer_km = 100),
  list(label = "200x200km",   lon_km = 200,  lat_km = 200,  buffer_km = 20),
  list(label = "20x20km",     lon_km = 20,   lat_km = 20,   buffer_km = 5)
)
TEMPORAL_RESOLUTIONS <- list(
  list(label = "1j", block_days = 1, buffer_days = 1)
)

# Nombre de folds pour le blocage : ADAPTATIF plutôt que fixe. Une
# résolution fine (20x20km) peut avoir des centaines de blocs -- n'en
# tirer que N_CV=10 jetterait trop d'information et donnerait une
# variance inter-fold peu fiable. Une résolution grossière (1500x1000km)
# peut n'avoir que quelques blocs -- viser un nombre fixe de folds
# identique à toutes les résolutions n'aurait pas de sens. On prend donc
# une fraction des blocs disponibles, bornée par un min et un max.
BLOCK_MIN_BLOCK_N        <- 50    # nb minimal d'observations dans un bloc pour l'utiliser
BLOCK_MAX_FOLDS_FRACTION <- 0.3   # fraction des blocs disponibles utilisée comme folds
BLOCK_MIN_FOLDS          <- 5     # plancher (même si peu de blocs dispo)
BLOCK_MAX_FOLDS_ABS      <- 30    # plafond (même si beaucoup de blocs dispo)

# Nombre de folds pour le naive RS 80/20 : reste fixe (ce n'est pas un
# nombre de "blocs disponibles", juste un nombre de répétitions Monte-Carlo).
NAIVE_N_FOLDS <- N_CV

set.seed(42)
