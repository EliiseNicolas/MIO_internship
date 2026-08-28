# ============================================================
# Regrillage du FTLE sur les grilles personnalisees
# ------------------------------------------------------------
# Meme logique que le regrillage pigments, mais une seule variable
# (pas de boucle "pigments"). La resolution native du FTLE est
# differente de celle des pigments/grilles cibles, mais ca n'a pas
# d'impact : assign_cell rattache chaque pixel FTLE au point le
# plus proche de la grille cible, quelle que soit sa resolution
# d'origine.
# ============================================================

library(dplyr)
library(ggplot2)

# --- Attribution au point de grille le plus proche (identique au pipeline NASC/pigments) ---
assign_cell <- function(coord, grid_vals) {
  grid_sorted <- sort(unique(grid_vals))
  idx <- findInterval(coord, grid_sorted, all.inside = TRUE)
  d_before <- abs(coord - grid_sorted[idx])
  d_after  <- abs(coord - grid_sorted[pmin(idx + 1, length(grid_sorted))])
  idx_final <- ifelse(d_after < d_before, idx + 1, idx)
  grid_sorted[idx_final]
}

# --- Moyenne (arithmetique) d'une matrice [lon x lat] regroupee par cellule ---
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

# --- Regrillage complet du FTLE sur une grille + stats ---
regrid_ftle_to_grid <- function(ftle, grid, grid_name) {
  
  lon_levels <- as.character(sort(unique(grid$lon)))
  lat_levels <- as.character(sort(unique(grid$lat)))
  lon_group <- factor(as.character(assign_cell(ftle$lon, grid$lon)), levels = lon_levels)
  lat_group <- factor(as.character(assign_cell(ftle$lat, grid$lat)), levels = lat_levels)
  
  n_new_lon <- nlevels(lon_group)
  n_new_lat <- nlevels(lat_group)
  n_time <- length(ftle$date)
  
  mean_arr <- array(NA_real_, dim = c(n_time, n_new_lon, n_new_lat))
  n_valid_arr <- array(0L, dim = c(n_time, n_new_lon, n_new_lat))
  
  for (d in seq_len(n_time)) {
    res <- regrid_matrix_mean(ftle$ftle[d, , ], lon_group, lat_group)
    mean_arr[d, , ] <- res$mean
    n_valid_arr[d, , ] <- res$n_valid
  }
  
  ftle_new <- list(
    date = ftle$date, lon = as.numeric(lon_levels), lat = as.numeric(lat_levels),
    ftle = mean_arr, n_valid = n_valid_arr
  )
  
  # (a) Total par cellule, toutes dates confondues : diagnostic spatial pur.
  n_valid_total_by_cell <- as.vector(apply(n_valid_arr, c(2, 3), sum))
  
  # (b) Par cellule ET par jour, seulement la ou une moyenne a ete calculee.
  n_valid_flat <- as.vector(n_valid_arr)
  n_valid_used <- n_valid_flat[n_valid_flat > 0]
  
  coverage <- data.frame(
    grid = grid_name, variable = "ftle",
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
  
  # (c) Nombre moyen de pixels utilises par cellule, PAR DATE (cellules touchees seulement).
  mean_n_valid_by_date <- numeric(n_time)
  n_cells_touched_by_date <- integer(n_time)
  for (d in seq_len(n_time)) {
    touched <- n_valid_arr[d, , ][n_valid_arr[d, , ] > 0]
    mean_n_valid_by_date[d] <- if (length(touched) > 0) mean(touched) else NA
    n_cells_touched_by_date[d] <- length(touched)
  }
  
  daily <- data.frame(
    grid = grid_name, variable = "ftle", date = ftle$date,
    mean_n_valid_cell = mean_n_valid_by_date, n_cells_touched = n_cells_touched_by_date
  )
  
  list(ftle_new = ftle_new, coverage = coverage, daily = daily)
}

# ============================================================
# BOUCLE SUR TOUTES LES GRILLES
# ============================================================

grid_output_dir <- "F:/data_elise/sv_cropped/grids_custom"
ftle_output_dir <- "F:/data_elise/ftle/grids_custom"
dir.create(ftle_output_dir, recursive = TRUE, showWarnings = FALSE)
path_out_figures <- "C:/Users/mmolinet/elisou_ta_stagiaire_pref/MIO_internship_III/figures/diag_regridding/regridding_ftle"
dir.create(path_out_figures, recursive = TRUE, showWarnings = FALSE)
grid_paths <- list.files(grid_output_dir, pattern = "^pigmeann_grid_.*\\.rds$", full.names = TRUE)

path_ftle <- "F:/data_elise/ftle/ftle_2018_2021_2022_2023_cropped.rds"
ftle <- readRDS(path_ftle)

coverage_summary_all <- list()
daily_summary_all <- list()

for (path_grid in grid_paths) {
  grid_name <- tools::file_path_sans_ext(basename(path_grid))
  cat("Grille :", grid_name, "\n")
  grid <- readRDS(path_grid)
  
  out <- regrid_ftle_to_grid(ftle, grid, grid_name)
  
  saveRDS(out$ftle_new, file.path(ftle_output_dir, paste0("ftle_", grid_name, ".rds")))
  coverage_summary_all[[grid_name]] <- out$coverage
  daily_summary_all[[grid_name]] <- out$daily
}

coverage_summary_df <- bind_rows(coverage_summary_all)
daily_summary_df <- bind_rows(daily_summary_all)

print(coverage_summary_df)
write.csv(coverage_summary_df, file.path(ftle_output_dir, "coverage_summary_ftle_grids.csv"), row.names = FALSE)
write.csv(daily_summary_df, file.path(ftle_output_dir, "daily_summary_ftle_grids.csv"), row.names = FALSE)

# ============================================================
# HISTOGRAMME : nombre moyen de pixels utilises par cellule, par date
# ============================================================

for (gname in unique(daily_summary_df$grid)) {
  
  df_grid <- daily_summary_df %>% filter(grid == gname)
  mean_global <- mean(df_grid$mean_n_valid_cell, na.rm = TRUE)
  
  p_daily <- ggplot(df_grid, aes(x = date, y = mean_n_valid_cell)) +
    geom_col(fill = "darkorange") +
    geom_hline(yintercept = mean_global, linetype = "dashed", color = "red", linewidth = 0.5) +
    labs(
      title = paste("Nombre moyen de pixels FTLE utilisés par cellule (post-regrillage) -", gname),
      subtitle = paste("Moyenne globale, toutes dates confondues =", round(mean_global, 1)),
      x = "Date", y = "Nb moyen de pixels valides par cellule"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6))
  
  print(p_daily)
  
  ggsave(
    file.path(path_out_figures, paste0("diag_mean_n_valid_par_date_", gname, ".png")),
    p_daily, width = 12, height = 6, dpi = 300
  )
}

# ============================================================
# PLOT : CV du nombre de points utilises par jour, par grille
# ============================================================

p_cv <- ggplot(coverage_summary_df, aes(x = grid, y = cv_n_valid_par_jour)) +
  geom_col(fill = "darkorange") +
  labs(
    title = "Variabilité du nombre de pixels utilisés par moyenne journalière (FTLE)",
    subtitle = "CV = écart-type / moyenne, calculé uniquement sur les cellule-jours où une moyenne a été produite",
    x = "Grille", y = "Coefficient de variation"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_cv)

ggsave(file.path(path_out_figures, "diag_cv_n_valid_ftle.png"), p_cv, width = 8, height = 6, dpi = 300)