# ============================================================
# Regrillage des pigments sur les grilles personnalisees
# ------------------------------------------------------------
# Reutilise les memes fichiers de grille que pour l'acoustique
# (pigmeann_grid_*.rds), pour garantir que pigments et acoustique
# tombent exactement sur les memes cellules spatiales.
# ============================================================

library(dplyr)
library(ggplot2)

# ============================================================
# Helper : orientation de l'array (time, lon, lat)
# ------------------------------------------------------------
# pigs$c_cond_* peut etre stocke en (time, lon, lat) ou (time, lat, lon)
# selon le pipeline d'origine. On detecte et corrige automatiquement
# pour eviter un desalignement silencieux avec pigs$lon / pigs$lat.
# ============================================================

fix_lonlat_order <- function(x, pigs) {
  d <- dim(x)
  n_lon <- length(pigs$lon)
  n_lat <- length(pigs$lat)
  if (d[2] == n_lon && d[3] == n_lat) return(x)
  if (d[2] == n_lat && d[3] == n_lon) return(aperm(x, c(1, 3, 2)))
  stop("Dimensions de l'array incompatibles avec pigs$lon / pigs$lat.")
}

# --- Attribution au point de grille le plus proche (identique au pipeline NASC) ---
assign_cell <- function(coord, grid_vals) {
  grid_sorted <- sort(unique(grid_vals))
  idx <- findInterval(coord, grid_sorted, all.inside = TRUE)
  d_before <- abs(coord - grid_sorted[idx])
  d_after  <- abs(coord - grid_sorted[pmin(idx + 1, length(grid_sorted))])
  idx_final <- ifelse(d_after < d_before, idx + 1, idx)
  grid_sorted[idx_final]
}

# --- Moyenne (arithmetique, pas de log) d'une matrice [lon x lat] regroupee par cellule ---
regrid_matrix_mean <- function(mat2d, lon_group, lat_group) {
  valid <- !is.na(mat2d)
  mat0 <- mat2d
  mat0[!valid] <- 0
  
  sum_lon <- rowsum(mat0, group = lon_group)
  count_lon <- rowsum(valid * 1, group = lon_group)
  
  sum_final <- t(rowsum(t(sum_lon), group = lat_group))
  count_final <- t(rowsum(t(count_lon), group = lat_group))
  
  mean_final <- sum_final / count_final
  mean_final[count_final == 0] <- NA
  list(mean = mean_final, n_valid = count_final)
}

# --- Regrillage d'une variable pigment, toutes dates ---
regrid_pigment_variable <- function(mat, lon_group, lat_group, n_new_lon, n_new_lat) {
  n_time <- dim(mat)[1]
  mean_arr <- array(NA_real_, dim = c(n_time, n_new_lon, n_new_lat))
  n_valid_arr <- array(0L, dim = c(n_time, n_new_lon, n_new_lat))
  
  for (d in seq_len(n_time)) {
    res <- regrid_matrix_mean(mat[d, , ], lon_group, lat_group)
    mean_arr[d, , ] <- res$mean
    n_valid_arr[d, , ] <- res$n_valid
  }
  list(mean = mean_arr, n_valid = n_valid_arr)
}

# --- Regrillage complet de tous les pigments sur une grille + stats ---
regrid_pigments_to_grid <- function(pigs, grid, pigments, grid_name) {
  
  lon_levels <- as.character(sort(unique(grid$lon)))
  lat_levels <- as.character(sort(unique(grid$lat)))
  lon_group <- factor(as.character(assign_cell(pigs$lon, grid$lon)), levels = lon_levels)
  lat_group <- factor(as.character(assign_cell(pigs$lat, grid$lat)), levels = lat_levels)
  
  n_new_lon <- nlevels(lon_group)
  n_new_lat <- nlevels(lat_group)
  n_time <- length(pigs$date)
  
  pigs_new <- list(date = pigs$date, lon = as.numeric(lon_levels), lat = as.numeric(lat_levels))
  coverage_rows <- list()
  daily_rows <- list()
  
  for (pig in pigments) {
    cat("  Regrid pigment :", pig, "\n")
    mat <- fix_lonlat_order(pigs[[paste0("c_cond_", pig)]], pigs)
    
    res <- regrid_pigment_variable(mat, lon_group, lat_group, n_new_lon, n_new_lat)
    pigs_new[[paste0("c_cond_", pig)]] <- res$mean
    pigs_new[[paste0("n_valid_", pig)]] <- res$n_valid
    
    # (a) Total par cellule, toutes dates confondues : diagnostic spatial pur.
    n_valid_total_by_cell <- as.vector(apply(res$n_valid, c(2, 3), sum))
    
    # (b) Par cellule ET par jour, seulement la ou une moyenne a ete calculee.
    n_valid_flat <- as.vector(res$n_valid)
    n_valid_used <- n_valid_flat[n_valid_flat > 0]
    
    coverage_rows[[pig]] <- data.frame(
      grid = grid_name, pigment = pig,
      n_cells = length(n_valid_total_by_cell),
      n_cells_jamais_touchees = sum(n_valid_total_by_cell == 0),
      mean_n_valid_total_cell = mean(n_valid_total_by_cell),
      sd_n_valid_total_cell = sd(n_valid_total_by_cell),
      cv_n_valid_total_cell = sd(n_valid_total_by_cell) / mean(n_valid_total_by_cell),
      n_cell_jours_valides = length(n_valid_used),
      mean_n_valid_par_jour = mean(n_valid_used),
      sd_n_valid_par_jour = sd(n_valid_used),
      cv_n_valid_par_jour = sd(n_valid_used) / mean(n_valid_used),
      min_n_valid_par_jour = min(n_valid_used), max_n_valid_par_jour = max(n_valid_used)
    )
    
    # (c) Nouveau : moyenne du nombre de pixels utilises par cellule, PAR DATE
    #     (uniquement sur les cellules touchees ce jour-la), pour visualiser
    #     l'evolution temporelle de la couverture nuageuse post-regrillage.
    mean_n_valid_by_date <- numeric(n_time)
    n_cells_touched_by_date <- integer(n_time)
    for (d in seq_len(n_time)) {
      touched <- res$n_valid[d, , ][res$n_valid[d, , ] > 0]
      mean_n_valid_by_date[d] <- if (length(touched) > 0) mean(touched) else NA
      n_cells_touched_by_date[d] <- length(touched)
    }
    
    daily_rows[[pig]] <- data.frame(
      grid = grid_name, pigment = pig, date = pigs$date,
      mean_n_valid_cell = mean_n_valid_by_date, n_cells_touched = n_cells_touched_by_date
    )
  }
  
  list(pigs_new = pigs_new, coverage = bind_rows(coverage_rows), daily = bind_rows(daily_rows))
}

# ============================================================
# BOUCLE SUR TOUTES LES GRILLES
# ============================================================

grid_output_dir <- "F:/data_elise/sv_cropped/grids_custom"
pig_output_dir  <- "F:/data_elise/pigmeann/grids_custom"
dir.create(pig_output_dir, recursive = TRUE, showWarnings = FALSE)
path_out_figures <- "C:/Users/mmolinet/elisou_ta_stagiaire_pref/MIO_internship_III/figures/diag_regridding/regridding_pig"
grid_paths <- list.files(grid_output_dir, pattern = "^pigmeann_grid_.*\\.rds$", full.names = TRUE)

path_pig <- "F:/data_elise/pigmeann/pigments_2018_2021_2022_2023_crop.rds"
pigs <- readRDS(path_pig)
pigments <- c("Chla", "Per", "But", "Fuco", "Hex", "Allo", "Zea", "Chlb", "DvChla")

coverage_summary_all <- list()
daily_summary_all <- list()

for (path_pig_grid in grid_paths) {
  grid_name <- tools::file_path_sans_ext(basename(path_pig_grid))
  cat("Grille :", grid_name, "\n")
  grid <- readRDS(path_pig_grid)
  
  out <- regrid_pigments_to_grid(pigs, grid, pigments, grid_name)
  
  saveRDS(out$pigs_new, file.path(pig_output_dir, paste0("pigments_", grid_name, ".rds")))
  coverage_summary_all[[grid_name]] <- out$coverage
  daily_summary_all[[grid_name]] <- out$daily
}

coverage_summary_df <- bind_rows(coverage_summary_all)
daily_summary_df <- bind_rows(daily_summary_all)

print(coverage_summary_df)
write.csv(coverage_summary_df, file.path(pig_output_dir, "coverage_summary_pigments_grids.csv"), row.names = FALSE)
write.csv(daily_summary_df, file.path(pig_output_dir, "daily_summary_pigments_grids.csv"), row.names = FALSE)


# ============================================================
# HISTOGRAMME : nombre moyen de pixels utilises par cellule, par date
# ------------------------------------------------------------
# Un plot par grille (facette par pigment), meme style que les
# histogrammes "ESU par jour" du pipeline NASC.
# ============================================================

for (gname in unique(daily_summary_df$grid)) {
  
  df_grid <- daily_summary_df %>% filter(grid == gname)
  mean_global <- mean(df_grid$mean_n_valid_cell, na.rm = TRUE)
  
  p_daily <- ggplot(df_grid, aes(x = date, y = mean_n_valid_cell)) +
    geom_col(fill = "steelblue") +
    geom_hline(yintercept = mean_global, linetype = "dashed", color = "red", linewidth = 0.5) +
    facet_wrap(~pigment, ncol = 3) +
    labs(
      title = paste("Nombre moyen de pixels utilisés par cellule (post-regrillage) -", gname),
      subtitle = paste("Moyenne globale, toutes dates/pigments confondus =", round(mean_global, 1)),
      x = "Date", y = "Nb moyen de pixels valides par cellule"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6))
  
  print(p_daily)
  
  ggsave(
    file.path(path_out_figures, paste0("diag_mean_n_valid_par_date_", gname, ".png")),
    p_daily, width = 14, height = 10, dpi = 300
  )
}

# ============================================================
# PLOT : CV du nombre de points utilises par jour, par grille et pigment
# ============================================================

p_cv <- ggplot(coverage_summary_df, aes(x = grid, y = cv_n_valid_par_jour, fill = pigment)) +
  geom_col(position = "dodge") +
  labs(
    title = "Variabilité du nombre de pixels utilisés par moyenne journalière (pigments)",
    subtitle = "CV = écart-type / moyenne, calculé uniquement sur les cellule-jours où une moyenne a été produite",
    x = "Grille", y = "Coefficient de variation"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_cv)

ggsave(file.path(path_out_figures, "diag_cv_n_valid_pigments.png"), p_cv, width = 10, height = 6, dpi = 300)