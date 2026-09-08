# ============================================================
# Profils moyens de température / salinité par cluster
# (clusters 1-6 + classes de transition 7-13)
# A partir des objets déjà sauvegardés par FOD_all.Rmd et
# add_transition_classes_fod.R
# ============================================================

# ---- Chemins ----
save_path <- "F:/data_elise/fod_elise_2018_2021_2022_2023/mclust_results_2018_2021_2022_2023"

# ---- Chargement des données ----
depth        <- readRDS(file.path(save_path, "depth.rds"))
temp_rec     <- readRDS(file.path(save_path, "temperature_profiles.rds"))  # (n_profils x n_depth)
sal_rec      <- readRDS(file.path(save_path, "salinity_profiles.rds"))     # (n_profils x n_depth)

# classification avec transitions renommées, alignée sur les profils
# (même longueur que nrow(temp_rec)/nrow(sal_rec))
cluster_transition_new <- readRDS(file.path(save_path, "transitions_upgraded/cluster_transition_renamed.rds"))

stopifnot(length(cluster_transition_new) == nrow(temp_rec))

# ============================================================
# TEMPLATE COULEURS / LEGENDE
# ============================================================

cluster_cols <- c(
  "1" = "#2166AC", "2" = "#67A9CF", "3" = "#1A9850",
  "4" = "#A6D96A", "5" = "#FDAE61", "6" = "#D73027"
)
transition_cols <- c(
  "7"  = "#3B73B9", "8"  = "#3FA7B5", "9"  = "#45B97C",
  "10" = "#C46A00", "11" = "#E85D04", "12" = "#C9184A",
  "13" = "#8F1D3F"
)
fod_cols <- c(cluster_cols, transition_cols)

legend_codes  <- c(1, 7, 2, 8, 9, 3, 10, 4, 11, 12, 5, 13, 6)
legend_labels <- c(
  "C1", "T1-2", "C2", "T1-3", "T2-3", "C3", "T3-4",
  "C4", "T4-6", "T4-5", "C5", "T5-6", "C6"
)

# ============================================================
# CALCUL DES STATISTIQUES PAR CLASSE (moyenne + IQ)
# ============================================================

mean_temp_cl <- list()
q1_temp_cl   <- list()
q3_temp_cl   <- list()

mean_sal_cl <- list()
q1_sal_cl   <- list()
q3_sal_cl   <- list()

for (code in legend_codes) {
  
  key <- as.character(code)
  ind <- cluster_transition_new == code
  
  if (sum(ind) == 0) next  # classe absente des données -> on saute
  
  mean_temp_cl[[key]] <- colMeans(temp_rec[ind, , drop = FALSE], na.rm = TRUE)
  q1_temp_cl[[key]]   <- apply(temp_rec[ind, , drop = FALSE], 2, quantile, probs = 0.25, na.rm = TRUE)
  q3_temp_cl[[key]]   <- apply(temp_rec[ind, , drop = FALSE], 2, quantile, probs = 0.75, na.rm = TRUE)
  
  mean_sal_cl[[key]] <- colMeans(sal_rec[ind, , drop = FALSE], na.rm = TRUE)
  q1_sal_cl[[key]]   <- apply(sal_rec[ind, , drop = FALSE], 2, quantile, probs = 0.25, na.rm = TRUE)
  q3_sal_cl[[key]]   <- apply(sal_rec[ind, , drop = FALSE], 2, quantile, probs = 0.75, na.rm = TRUE)
}

# codes réellement présents, dans l'ordre du template
present_codes  <- legend_codes[as.character(legend_codes) %in% names(mean_temp_cl)]
present_labels <- legend_labels[legend_codes %in% present_codes]

# ============================================================
# FONCTION DE PLOT GENERIQUE
# ============================================================

plot_profiles_cluster <- function(mean_cl, q1_cl, q3_cl, depth,
                                  codes, labels, cols,
                                  xlab, out_file,
                                  width = 2000, height = 2800, res = 300) {
  
  png(out_file, width = width, height = height, res = res)
  
  xr <- range(unlist(q1_cl[as.character(codes)]),
              unlist(q3_cl[as.character(codes)]))
  
  plot(
    NULL,
    xlim = xr,
    ylim = rev(range(depth)),
    xlab = xlab,
    ylab = "Depth (m)"
  )
  
  for (code in codes) {
    
    key <- as.character(code)
    col <- cols[key]
    col_fill <- adjustcolor(col, alpha.f = 0.25)
    
    polygon(
      x = c(q1_cl[[key]], rev(q3_cl[[key]])),
      y = c(depth, rev(depth)),
      col = col_fill,
      border = NA
    )
    
    lines(
      mean_cl[[key]],
      depth,
      col = col,
      lwd = 3
    )
  }
  
  legend(
    "bottomright",
    legend = labels,
    col = cols[as.character(codes)],
    lwd = 3,
    bty = "n",
    cex = 0.85,
    ncol = 2
  )
  
  dev.off()
}

# ============================================================
# TRACE TEMPERATURE
# ============================================================

plot_profiles_cluster(
  mean_cl  = mean_temp_cl,
  q1_cl    = q1_temp_cl,
  q3_cl    = q3_temp_cl,
  depth    = depth,
  codes    = present_codes,
  labels   = present_labels,
  cols     = fod_cols,
  xlab     = "Temperature (°C)",
  out_file = file.path(save_path, "temperature_clusters_transitions.png")
)

# ============================================================
# TRACE SALINITE
# ============================================================

plot_profiles_cluster(
  mean_cl  = mean_sal_cl,
  q1_cl    = q1_sal_cl,
  q3_cl    = q3_sal_cl,
  depth    = depth,
  codes    = present_codes,
  labels   = present_labels,
  cols     = fod_cols,
  xlab     = "Salinity (PSU)",
  out_file = file.path(save_path, "salinity_clusters_transitions.png")
)

# ============================================================
# SAUVEGARDE DES STATISTIQUES (optionnel, pour réutilisation)
# ============================================================

saveRDS(mean_temp_cl, file = file.path(save_path, "mean_temperature_clusters_transitions.rds"))
saveRDS(q1_temp_cl,   file = file.path(save_path, "q1_temperature_clusters_transitions.rds"))
saveRDS(q3_temp_cl,   file = file.path(save_path, "q3_temperature_clusters_transitions.rds"))

saveRDS(mean_sal_cl, file = file.path(save_path, "mean_salinity_clusters_transitions.rds"))
saveRDS(q1_sal_cl,   file = file.path(save_path, "q1_salinity_clusters_transitions.rds"))
saveRDS(q3_sal_cl,   file = file.path(save_path, "q3_salinity_clusters_transitions.rds"))