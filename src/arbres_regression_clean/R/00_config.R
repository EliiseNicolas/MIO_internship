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
  library(randomForestSRC)   # nécessaire pour 14_run_rfsrc_reconstruction.R
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
    "new_ratio_ds_NASC_per_esu_pig_ftle_fod_2018_2021_2022_2023_transect_",
    freq, "kHz_mask9.rds"
  )
}

# ---- Grille de prédiction spatiale ------------------------------------------
# UN SEUL fichier désormais (après correction du script de génération,
# cf. data_generation/generate_ds_ftle_pig_fod_all_dates.R) : toutes les
# dates -- y compris 2023-01-26 -- avec pigments bruts ET versions
# normalisées déjà calculées (chla_total, chla_totpig, ..., en minuscules),
# et `ftle`/`pig`/`fod` tous alignés sur le même ordre de dimensions
# [date, lon, lat]. ADAPTER le chemin ci-dessous à ton fichier réel.
PATH_GRID_ALL_DATES <- "F:/data_elise/prediction_ds/ds_ftle_pig_fod_ALL_DATES.rds"

# Date utilisée pour les cartes "mono-date" (12_run_grid_prediction.R) --
# extraite directement du fichier ci-dessus, plus besoin d'un fichier séparé.
TARGET_DATE_SINGLE <- as.Date("2023-01-26")

# Table de correspondance entre les noms (minuscules) des champs déjà
# normalisés dans le fichier et les noms utilisés dans COVARIATES_NUM
# (voir plus bas) -- PAS de formule à deviner, ces champs sont déjà
# calculés dans le fichier source.
MULTIDATE_PIG_NAME_MAP <- c(
  chla_total    = "Chla_total",
  chla_totpig   = "Chla_totpig",
  per_totpig    = "Per_totpig",
  but_totpig    = "But_totpig",
  fuco_totpig   = "Fuco_totpig",
  hex_totpig    = "Hex_totpig",
  allo_totpig   = "Allo_totpig",
  zea_totpig    = "Zea_totpig",
  chlb_totpig   = "Chlb_totpig",
  dvchla_totpig = "DvChla_totpig"
)

OUTPUT_ROOT <- "outputs_pipeline"
dir.create(OUTPUT_ROOT, showWarnings = FALSE, recursive = TRUE)

path_out <- function(...) file.path(OUTPUT_ROOT, ...)

# ---- Variables ---------------------------------------------------------------
# Structure de données mise à jour : les prédicteurs sont désormais
# fournis directement (déjà normalisés par rapport au pigment total),
# plus besoin de calculer des ratios à la main.
#   - ftle              : directement dans la table
#   - fod               : facteur (chaîne "NA" à convertir en NA)
#   - Chla_total        : Chla normalisé (variable distincte de Chla_totpig)
#   - <pigment>_totpig  : chaque pigment normalisé par le pigment total
#     (Chla_totpig, Per_totpig, But_totpig, Fuco_totpig, Hex_totpig,
#      Allo_totpig, Zea_totpig, Chlb_totpig, DvChla_totpig)
RESPONSE_VAR   <- "NASC"   # construit à partir de la colonne source `nasc` (log10)

COVARIATES_NUM <- c(
  "ftle", "Chla_total",
  "Chla_totpig", "Per_totpig", "But_totpig", "Fuco_totpig",
  "Hex_totpig", "Allo_totpig", "Zea_totpig", "Chlb_totpig", "DvChla_totpig"
)
COVARIATES_ALL <- c(COVARIATES_NUM, "fod")

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
# nombre de "blocs disponibles", juste un nombre de répétitions/folds).
NAIVE_N_FOLDS <- N_CV

# Méthode de CV pour le naive RS 80/20 :
#   "kfold_repeated" (recommandé) : répétitions de K-fold CV classique
#     (partition disjointe, K dérivé de frac_train). Voir la justification
#     complète dans 02_folds.R.
# Méthode de CV pour le naive RS 80/20 :
#   "kfold_repeated" (recommandé pour un nouveau tuning) : répétitions de
#     K-fold CV classique (partition disjointe, K dérivé de frac_train).
#     Voir la justification complète dans 02_folds.R.
#   "monte_carlo" : tirages aléatoires indépendants (chevauchement des
#     tests possible) -- ACTIF ici car c'est la méthode déjà utilisée
#     pour le tuning existant (10_run_tuning.R) : il faut la même méthode
#     pour le training/les scripts d'analyse, sinon les folds ne
#     correspondent plus aux hyperparamètres déjà tunés.
NAIVE_CV_METHOD <- "monte_carlo"

# ---- Échelles partagées entre plots (axes Y / colorbars identiques) --------
# Pour permettre une comparaison visuelle DIRECTE entre modèles, schémas
# et fréquences (même axe = même valeur d'un plot à l'autre).
#   "global"   : UNE seule échelle pour tout, toutes fréquences confondues.
#     ATTENTION : 38 kHz et 120 kHz sont deux fréquences acoustiques
#     distinctes, dont le NASC peut avoir une magnitude intrinsèquement
#     différente (pas juste un effet du modèle) -- forcer la même échelle
#     conflate performance du modèle et différence physique entre
#     fréquences. Défendable seulement si tu veux explicitement comparer
#     les deux fréquences sur un pied d'égalité visuel.
#   "per_freq" : une échelle par fréquence (modèles/schémas partagés à
#     l'intérieur d'une fréquence, mais 38 kHz et 120 kHz peuvent différer).
#     Plus défendable scientifiquement si les deux fréquences ont des
#     gammes de NASC différentes -- à envisager si les plots "global"
#     donnent l'impression qu'une fréquence "écrase" l'autre visuellement.
SHARED_SCALE_SCOPE <- "global"

# ---- Prédiction sur la grille multi-date (133 jours) ------------------------
# 133 dates x plusieurs schémas x 2 modèles peut vite représenter des
# milliers de cartes. Par défaut on ne produit les cartes datées que pour
# le schéma naive (le plus rapide/simple) -- étends cette liste si tu
# veux aussi les 133 cartes pour un ou plusieurs schémas bloqués.
MULTIDATE_PREDICTION_SCHEMES <- c("naive_RS_80_20")

set.seed(42)
