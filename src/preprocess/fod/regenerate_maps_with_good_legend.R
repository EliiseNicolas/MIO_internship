
save_path <- "F:/data_elise/fod_elise_2018_2021_2022_2023/mclust_results_2018_2021_2022_2023/transitions_upgraded"

# ------------------------------------------------------------
# Chargement des variables nécessaires
# ------------------------------------------------------------

lon      <- readRDS(file.path(save_path, "lon.rds"))
lat      <- readRDS(file.path(save_path, "lat.rds"))
time_bis <- readRDS(file.path(save_path, "time.rds"))

cluster_transition_map_new <- readRDS(file.path(save_path, "cluster_transition_map_renamed.rds"))
# Si vous préférez repartir du vecteur + du mask plutôt que de la map déjà reconstruite :
# cluster_transition_new <- readRDS(file.path(save_path, "cluster_transition_renamed.rds"))
# mask_common            <- readRDS(file.path(save_path, "mask_common.rds"))
# cluster_transition_flat_new <- rep(NA_integer_, length(mask_common))
# cluster_transition_flat_new[mask_common] <- cluster_transition_new
# cluster_transition_map_new <- array(cluster_transition_flat_new,
#                                      dim = c(length(lon), length(lat), length(time_bis)))

# Vérification rapide des dimensions
stopifnot(dim(cluster_transition_map_new)[1] == length(lon))
stopifnot(dim(cluster_transition_map_new)[2] == length(lat))
stopifnot(dim(cluster_transition_map_new)[3] == length(time_bis))

# ------------------------------------------------------------
# Dossier de sortie pour les nouvelles cartes
# ------------------------------------------------------------

out_dir <- file.path(save_path, "maps_regenerated")
dir.create(out_dir, showWarnings = FALSE)

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
# ORDRE GEOGRAPHIQUE SUD -> NORD (LEGENDE CORRIGEE)
# Ordre : C1 -> T1-2 -> C2 -> T1-3 -> T2-3 -> C3 -> T3-4 ->
#         C4 -> T4-6 -> T4-5 -> C5 -> T5-6 -> C6
# (le code 8 = T1-3 avait été omis dans la version précédente)
# ============================================================

legend_codes <- c(
  1, 7, 2, 8, 9, 3, 10, 4, 11, 12, 5, 13, 6
)

legend_labels <- c(
  "C1", "T1-2",
  "C2", "T1-3", "T2-3",
  "C3", "T3-4",
  "C4", "T4-6",
  "T4-5", "C5",
  "T5-6", "C6"
)

# Vérification : tous les codes de transition_cols/cluster_cols sont bien présents
stopifnot(all(1:6 %in% legend_codes))
stopifnot(all(7:13 %in% legend_codes))

legend_colors <- cols[legend_codes + 1]

# ============================================================
# GENERATION DES CARTES
# ============================================================

for (i in seq_along(time_bis)) {
  
  png(
    sprintf(
      file.path(out_dir, "cluster_transition_map_renamed_%s.png"),
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
  
  # ----------------------------------------------------------
  # CARTE
  # ----------------------------------------------------------
  
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
  
  # ----------------------------------------------------------
  # LEGENDE HORIZONTALE
  # ----------------------------------------------------------
  
  legend_positions <- seq(
    min(lon) + 3,
    max(lon) - 3,
    length.out = length(legend_codes)
  )
  
  legend_y <- min(lat) - 12
  legend_label_y <- min(lat) - 17
  
  points(
    x = legend_positions,
    y = rep(legend_y, length(legend_positions)),
    pch = 15,
    cex = 2.4,
    col = legend_colors
  )
  
  text(
    x = legend_positions,
    y = rep(legend_label_y, length(legend_positions)),
    labels = legend_labels,
    cex = 1.05,
    font = 2
  )
  
  text(
    x = mean(range(lon)),
    y = min(lat) - 23,
    labels = "South → North",
    cex = 1.1,
    font = 2
  )
  
  dev.off()
}

cat("Cartes régénérées dans :", out_dir, "\n")

# ============================================================
# TABLEAU DE VERIFICATION (compte les 13 classes + NA)
# ============================================================

tab <- table(cluster_transition_map_new, useNA = "ifany")
print(tab)
cat("Nombre de classes non-NA :", length(tab) - any(is.na(names(tab))), "\n")