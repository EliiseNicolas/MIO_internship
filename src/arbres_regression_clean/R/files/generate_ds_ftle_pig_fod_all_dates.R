# ============================================================
# Extraction FTLE, pigments et FOD pour TOUTES les dates
# -> un seul objet / un seul fichier RDS en sortie
# Version entièrement vectorisée : aucune boucle sur les dates.
#
# CORRECTIF (par rapport à la version précédente) :
#   `fod_all` avait un ordre de dimensions [lon, lat, date],
#   incohérent avec `ftle_all` et `pig_all` qui ont [date, lon, lat].
#   Cette incohérence obligeait tout code consommateur à traiter `fod`
#   comme un cas particulier (risque d'erreur d'indexation). On
#   applique maintenant aperm() pour ramener `fod_all` à
#   [date, lon, lat], comme les deux autres variables.
#
# DEUX POINTS NON MODIFIÉS ICI (choix scientifiques à confirmer,
# pas des bugs de structure) -- voir les commentaires marqués
# "A CONFIRMER" plus bas :
#   - chla_total est une copie exacte de Chla (DvChla exclue) ;
#   - le dénominateur des ratios *_totpig exclut Chla de la somme,
#      y compris pour calculer chla_totpig lui-même.
# ============================================================

rm(list = ls())

# ------------------------------------------------------------
# Options
# ------------------------------------------------------------

out_file <- "F:/data_elise/prediction_ds/ds_ftle_pig_fod_ALL_DATES.rds"

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

path_ftle <- "F:/data_elise/ftle/ftle_2018_2021_2022_2023_cropped.rds"
path_pig  <- "F:/data_elise/pigmeann/pigments_2018_2021_2022_2023_crop.rds"

path_fod_clusters <- "F:/data_elise/fod_elise_2018_2021_2022_2023/mclust_results_2018_2021_2022_2023/transitions_upgraded/cluster_transition_map_renamed.rds"
path_fod_lon  <- "F:/data_elise/fod_elise_2018_2021_2022_2023/mclust_results_2018_2021_2022_2023/lon.rds"
path_fod_lat  <- "F:/data_elise/fod_elise_2018_2021_2022_2023/mclust_results_2018_2021_2022_2023/lat.rds"
path_fod_time <- "F:/data_elise/fod_elise_2018_2021_2022_2023/mclust_results_2018_2021_2022_2023/time.rds"

# ------------------------------------------------------------
# Chargement
# ------------------------------------------------------------

ftle_grid <- readRDS(path_ftle)
pigs_grid <- readRDS(path_pig)

fod_clusters <- readRDS(path_fod_clusters)
fod_lon  <- readRDS(path_fod_lon)
fod_lat  <- readRDS(path_fod_lat)
fod_time <- readRDS(path_fod_time)
fod_date <- as.Date(fod_time)

list_pigs <- c("Chla", "Per", "But", "Fuco", "Hex", "Allo", "Zea", "Chlb", "DvChla")

# ------------------------------------------------------------
# Grille de référence = pigments (1080 x 720)
# ------------------------------------------------------------

grid_lon <- pigs_grid$lon
grid_lat <- pigs_grid$lat

# ------------------------------------------------------------
# Indices de reprojection spatiale (nearest neighbor), 1 seule fois
# ------------------------------------------------------------

idx_lon_ftle <- vapply(grid_lon, function(x) which.min(abs(ftle_grid$lon - x)), integer(1))
idx_lat_ftle <- vapply(grid_lat, function(x) which.min(abs(ftle_grid$lat - x)), integer(1))

idx_lon_fod <- vapply(grid_lon, function(x) which.min(abs(fod_lon - x)), integer(1))
idx_lat_fod <- vapply(grid_lat, function(x) which.min(abs(fod_lat - x)), integer(1))

# ------------------------------------------------------------
# Dates communes aux trois jeux de données
# ------------------------------------------------------------

dates_communes <- sort(as.Date(Reduce(intersect, list(ftle_grid$date, pigs_grid$date, fod_date))))
n_date <- length(dates_communes)
cat(n_date, "dates communes.\n")

idx_date_ftle <- match(dates_communes, ftle_grid$date)
idx_date_pig  <- match(dates_communes, pigs_grid$date)
idx_date_fod  <- match(dates_communes, fod_date)

# ------------------------------------------------------------
# Extraction + reprojection en UNE SEULE opération d'indexation
# par variable (produit cartésien des indices sur les 3 dimensions)
# ------------------------------------------------------------

# --- FTLE : natif [ndate_tot, 901, 600] -> [n_date, 1080, 720] ---
ftle_all <- ftle_grid$ftle[idx_date_ftle, idx_lon_ftle, idx_lat_ftle]

# --- FOD : natif [541, 361, ndate_tot] -> extraction [1080, 720, n_date] ---
# PUIS reordonnancement en [n_date, 1080, 720] via aperm(), pour que
# fod_all ait EXACTEMENT le même ordre de dimensions que ftle_all et
# les pigments ci-dessous (date en premier). C'est le correctif
# principal de cette version du script.
fod_all_raw <- fod_clusters[idx_lon_fod, idx_lat_fod, idx_date_fod]   # [1080, 720, n_date]
fod_all     <- aperm(fod_all_raw, c(3, 1, 2))                        # [n_date, 1080, 720]
rm(fod_all_raw)

# --- Pigments : déjà sur la grille de référence, juste sélection des dates ---
pig_all <- lapply(list_pigs, function(p) {
  pigs_grid[[paste0("c_cond_", p)]][idx_date_pig, , ]   # [n_date, 1080, 720]
})
names(pig_all) <- list_pigs

# ------------------------------------------------------------
# chla_total = Chla seule (DvChla volontairement exclue)
# A CONFIRMER : ceci duplique exactement pig_all$Chla (même valeurs).
# Si c'est bien voulu (garder la concentration brute de Chla à côté
# des ratios, sous un nom dédié), rien à changer. Si ce n'était pas
# l'intention, remplace la ligne ci-dessous par la formule voulue
# (ex. Chla / (Chla + DvChla) pour une "chlorophylle a totale" au
# sens phytoplanctonique classique).
# ------------------------------------------------------------

pig_all[["chla_total"]] <- pig_all[["Chla"]]

# ------------------------------------------------------------
# Somme des concentrations de tous les pigments SAUF Chla
# (calculee une seule fois, reutilisee comme denominateur commun)
#
# A CONFIRMER : ce dénominateur exclut Chla de la somme pour TOUS les
# ratios *_totpig, y compris chla_totpig (= Chla / somme des 8 AUTRES
# pigments, pas Chla / somme des 9 pigments). Si l'intention est une
# vraie fraction du pigment total (incluant Chla dans la somme pour
# tout le monde), remplace `pigs_sans_chla` par `list_pigs` complet
# ci-dessous (ligne commentée en alternative).
# ------------------------------------------------------------

pigs_sans_chla <- setdiff(list_pigs, "Chla")

pig_sum_others <- Reduce(`+`, pig_all[pigs_sans_chla])   # somme elementwise, meme dims [n_date, 1080, 720]

# Alternative (total incluant Chla, decommenter si c'est l'intention) :
# pig_sum_others <- Reduce(`+`, pig_all[list_pigs])

# ------------------------------------------------------------
# Ratio de chaque pigment (les 9 pigments bruts, Chla incluse) sur
# la somme des pigments hors Chla, ajoutes dans pig_all avec un
# suffixe "_totpig" (ex : chla_totpig, but_totpig, fuco_totpig, ...).
# Vectorise via lapply : une division element par element par
# pigment, aucune boucle sur les dates/pixels.
# ------------------------------------------------------------

ratio_list <- lapply(list_pigs, function(p) pig_all[[p]] / pig_sum_others)
names(ratio_list) <- paste0(tolower(list_pigs), "_totpig")

pig_all <- c(pig_all, ratio_list)

# ------------------------------------------------------------
# Assemblage final : un seul objet pour toutes les dates
#
# `ftle`, `pig` ET `fod` ont maintenant TOUS le même ordre de
# dimensions [n_date, 1080, 720] -- plus de cas particulier pour fod
# côté consommateur du fichier.
# ------------------------------------------------------------

all_ds <- list(
  date = dates_communes,   # vecteur de n_date dates
  lon  = grid_lon,
  lat  = grid_lat,
  ftle = ftle_all,         # array [n_date, 1080, 720]
  pig  = pig_all,          # liste : 9 pigments bruts + chla_total + ratios *_totpig
  fod  = fod_all           # array [n_date, 1080, 720]  (corrige, etait [1080,720,n_date])
)

str(all_ds, max.level = 2)
cat("Taille en mémoire :", format(object.size(all_ds), units = "GB"), "\n")

# ------------------------------------------------------------
# Sauvegarde (un seul fichier)
# ------------------------------------------------------------

saveRDS(all_ds, out_file)
cat("Sauvegardé dans :", out_file, "\n")
