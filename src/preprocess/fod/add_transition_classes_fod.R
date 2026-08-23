# Description
# 
# A partir des données de FOD sauvegardées, On veut restructurer la classe de transition 0
# en plusieurs classes de transition 

# Libraries

# Global variables
save_path <-"F:/data_elise/fod_elise_2018_2021_2022_2023/mclust_results_2018_2021_2022_2023"

# Chargement des résultats
cluster_soft <- readRDS(file.path(save_path, "cluster_soft.rds"))
cluster_prob <- readRDS(file.path(save_path, "cluster_probabilities.rds"))
cluster_map <- readRDS(file.path(save_path, "cluster_map.rds"))
mask_common <- readRDS(file.path(save_path, "mask_common.rds"))
lon <- readRDS(file.path(save_path, "lon.rds"))
lat <- readRDS(file.path(save_path, "lat.rds"))
time_bis <- readRDS(file.path(save_path, "time.rds"))
nclust <- readRDS(file.path(save_path, "nclust.rds"))
seuil <- readRDS(file.path(save_path, "softclass_threshold.rds"))

nlon <- length(lon)
nlat <- length(lat)
ntime <- length(time_bis)

# Identification des profils du cluster 0
idx_transition <- which(cluster_soft == 0)
prob_transition <- cluster_prob[idx_transition, , drop = FALSE]

# Deux clusters ayant les probabilités les plus élevées pour chaque profil
top2 <- t(apply(prob_transition, 1, function(x) order(x, decreasing = TRUE)[1:2]))
colnames(top2) <- c("cluster1", "cluster2")

# ----------------------------------------- Définir les transitions
# Codes utilisés pour les nouvelles classes
# 0 = transition non identifiée
# 1-6 = clusters Mclust
# 7 = transition 1-2
# 8 = transition 1-4
# 9 = transition 2-4
# 10 = transition 3-4
# 11 = transition 3-5
# 12 = transition 3-6
# 13 = transition 5-6

transition_codes <- c("1-2" = 7, "1-4" = 8, "2-4" = 9, "3-4" = 10, "3-5" = 11, "3-6" = 12, "5-6" = 13)
transition_labels <- c("0" = "Transition non identifiée", "1" = "Cluster 1", "2" = "Cluster 2", "3" = "Cluster 3", "4" = "Cluster 4", "5" = "Cluster 5", "6" = "Cluster 6", "7" = "Transition 1-2", "8" = "Transition 1-4", "9" = "Transition 2-4", "10" = "Transition 3-4", "11" = "Transition 3-5", "12" = "Transition 3-6", "13" = "Transition 5-6")

# ----------------------------------------- Créer les nouvelles classes de transition
cluster_transition <- cluster_soft

for (i in seq_len(nrow(top2))) {
  cl1 <- top2[i, 1]
  cl2 <- top2[i, 2]
  pair <- paste(sort(c(cl1, cl2)), collapse = "-")
  if (pair %in% names(transition_codes)) cluster_transition[idx_transition[i]] <- transition_codes[pair]
}

# ---------------------------------------- Infos sur les transitions
transition_info <- data.frame(profile = idx_transition, cluster_1 = top2[, 1], cluster_2 = top2[, 2], prob_1 = prob_transition[cbind(seq_len(nrow(prob_transition)), top2[, 1])], prob_2 = prob_transition[cbind(seq_len(nrow(prob_transition)), top2[, 2])])
transition_info$pair <- apply(transition_info[, c("cluster_1", "cluster_2")], 1, function(x) paste(sort(x), collapse = "-"))
transition_info$transition_cluster <- cluster_transition[transition_info$profile]

# ---------------------------------------- Remettre sur grille spatiale et temporelle
cluster_transition_flat <- rep(NA_integer_, length(mask_common))
cluster_transition_flat[mask_common] <- cluster_transition
cluster_transition_map <- array(cluster_transition_flat, dim = c(nlon, nlat, ntime))

# ============================================================
# PROBABILITES SUR LA GRILLE
# ============================================================

prob_maps <- vector("list", nclust)

for (cl in seq_len(nclust)) {
  prob_flat <- rep(NA_real_, length(mask_common))
  prob_flat[mask_common] <- cluster_prob[, cl]
  prob_maps[[cl]] <- array(prob_flat, dim = c(nlon, nlat, ntime))
}

names(prob_maps) <- paste0("cluster_", seq_len(nclust))

# ============================================================
# VERIFICATIONS
# ============================================================

print(table(cluster_soft))
print(table(cluster_transition))
print(data.frame(cluster = names(table(cluster_transition)), label = transition_labels[names(table(cluster_transition))], n = as.vector(table(cluster_transition))))

# ============================================================
# SAUVEGARDE
# ============================================================

saveRDS(cluster_transition, file = file.path(save_path, "cluster_transition.rds"))
saveRDS(cluster_transition_flat, file = file.path(save_path, "cluster_transition_flat.rds"))
saveRDS(cluster_transition_map, file = file.path(save_path, "cluster_transition_map.rds"))
saveRDS(transition_codes, file = file.path(save_path, "transition_codes.rds"))
saveRDS(transition_info, file = file.path(save_path, "transition_information.rds"))
saveRDS(prob_maps, file = file.path(save_path, "cluster_probability_maps.rds"))
# Correspondance anciens clusters -> nouveaux clusters
cluster_rename <- c("1" = 1, "2" = 2, "3" = 4, "4" = 3, "5" = 6, "6" = 5)

# Correspondance anciens codes de transition -> nouveaux codes
# Anciennes transitions :
# 7  = 1-2  -> nouveau 1-2
# 8  = 1-4  -> nouveau 1-3
# 9  = 2-4  -> nouveau 2-3
# 10 = 3-4  -> nouveau 3-4
# 11 = 3-5  -> nouveau 4-6
# 12 = 3-6  -> nouveau 4-5
# 13 = 5-6  -> nouveau 5-6

transition_rename <- c("7" = 7, "8" = 8, "9" = 9, "10" = 10, "11" = 11, "12" = 12, "13" = 13)

# Nouveau vecteur de classification
cluster_transition_new <- cluster_transition

# Renommer les clusters 1 à 6
for (old in 1:6) {
  cluster_transition_new[cluster_transition == old] <- cluster_rename[as.character(old)]
}

# Les transitions gardent leurs codes numériques
for (old in 7:13) {
  cluster_transition_new[cluster_transition == old] <- transition_rename[as.character(old)]
}

# Les 0 restent 0
cluster_transition_new[cluster_transition == 0] <- 0

# ============================================================
# REMISE SUR LA GRILLE
# ============================================================

cluster_transition_flat_new <- rep(NA_integer_, length(mask_common))
cluster_transition_flat_new[mask_common] <- cluster_transition_new

cluster_transition_map_new <- array(cluster_transition_flat_new, dim = c(length(lon), length(lat), length(time_bis)))

# ============================================================
# NOUVELLE NOMENCLATURE DES TRANSITIONS
# ============================================================

transition_labels_new <- c("1-2", "1-3", "2-3", "3-4", "4-6", "4-5", "5-6")

transition_summary_new <- data.frame(code = 7:13, transition = transition_labels_new, n = as.integer(table(factor(cluster_transition_new, levels = 7:13))))

print(transition_summary_new)

# ============================================================
# SAUVEGARDE
# ============================================================

saveRDS(cluster_transition_new, file = file.path(save_path, "cluster_transition_renamed.rds"))
saveRDS(cluster_transition_flat_new, file = file.path(save_path, "cluster_transition_flat_renamed.rds"))
saveRDS(cluster_transition_map_new, file = file.path(save_path, "cluster_transition_map_renamed.rds"))
saveRDS(cluster_rename, file = file.path(save_path, "cluster_rename.rds"))
saveRDS(transition_summary_new, file = file.path(save_path, "transition_summary_renamed.rds"))

# ============================================================
# CARTES AVEC DEGRADE LATITUDINAL DES CLUSTERS
# Ordre Sud -> Nord : 1 -> 2 -> 4 -> 3 -> 6 -> 5
# ============================================================

library(fields)

# ============================================================
# CARTES AVEC NOUVELLE NOMENCLATURE
# Ordre Sud -> Nord :
# C1 -> T1-2 -> C2 -> T2-3 -> C3 -> T3-4 ->
# C4 -> T4-6 -> T4-5 -> C5 -> T5-6 -> C6
# ============================================================

library(fields)

# ============================================================
# COULEURS DES CLUSTERS ET DES TRANSITIONS
# ============================================================

cluster_cols <- c(
  "1" = "#2166AC",
  "2" = "#67A9CF",
  "3" = "#1A9850",
  "4" = "#A6D96A",
  "5" = "#FDAE61",
  "6" = "#D73027"
)

transition_cols <- c(
  "7"  = "#3B73B9",  # 1-2
  "8"  = "#3FA7B5",  # 1-3
  "9"  = "#45B97C",  # 2-3
  "10" = "#C46A00",  # 3-4
  "11" = "#E85D04",  # 4-6
  "12" = "#C9184A",  # 4-5
  "13" = "#8F1D3F"   # 5-6
)

cols <- c(
  "grey80",
  cluster_cols[c("1", "2", "3", "4", "5", "6")],
  transition_cols
)

breaks <- seq(-0.5, 13.5, 1)

# ============================================================
# ORDRE GEOGRAPHIQUE SUD -> NORD
# ============================================================

legend_codes <- c(
  1, 7, 2, 9, 3, 10, 4, 11, 12, 5, 13, 6
)

legend_labels <- c(
  "C1", "T1-2",
  "C2", "T2-3",
  "C3", "T3-4",
  "C4", "T4-6",
  "T4-5", "C5",
  "T5-6", "C6"
)

legend_colors <- cols[legend_codes + 1]

# ============================================================
# CARTES
# ============================================================

for (i in seq_along(time_bis)) {
  
  png(
    sprintf(
      file.path(
        save_path,
        "cluster_transition_map_renamed_%s.png"
      ),
      format(time_bis[i], "%Y%m%d")
    ),
    width = 2400,
    height = 1800,
    res = 300
  )
  
  # Marge basse augmentée pour descendre la légende
  par(
    mar = c(12, 4, 4, 2),
    xpd = NA
  )
  
  # ==========================================================
  # CARTE
  # ==========================================================
  
  image(
    x = lon,
    y = lat,
    z = cluster_transition_map_new[, , i],
    col = cols,
    breaks = breaks,
    xlab = "Longitude",
    ylab = "Latitude"
  )
  
  title(
    main = paste(
      "Functional Oceanographic Domains (clusters)",
      format(time_bis[i], "%Y-%m-%d")
    ),
    cex.main = 1.4
  )
  
  # ==========================================================
  # LEGENDE HORIZONTALE
  # ==========================================================
  
  legend_positions <- seq(
    min(lon) + 3,
    max(lon) - 3,
    length.out = length(legend_codes)
  )
  
  # Position des carrés
  legend_y <- min(lat) - 12
  
  # Position des labels
  legend_label_y <- min(lat) - 17
  
  # Carrés de couleur
  points(
    x = legend_positions,
    y = rep(legend_y, length(legend_positions)),
    pch = 15,
    cex = 2.4,
    col = legend_colors
  )
  
  # Labels
  text(
    x = legend_positions,
    y = rep(legend_label_y, length(legend_positions)),
    labels = legend_labels,
    cex = 1.05,
    font = 2
  )
  
  # Indication du gradient
  text(
    x = mean(range(lon)),
    y = min(lat) - 23,
    labels = "South → North",
    cex = 1.1,
    font = 2
  )
  
  # Fermeture du fichier PNG
  dev.off()
}

# ============================================================
# TABLEAU FINAL DES TRANSITIONS
# ============================================================

transition_summary_new <- data.frame(
  code = 7:13,
  transition = c(
    "1-2",
    "1-3",
    "2-3",
    "3-4",
    "4-6",
    "4-5",
    "5-6"
  ),
  n = as.integer(
    table(
      factor(
        cluster_transition_new,
        levels = 7:13
      )
    )
  )
)

print(transition_summary_new)

