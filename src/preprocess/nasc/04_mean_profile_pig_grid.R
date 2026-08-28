# Description

# 01_NASC_filter_sv.R
# From obsaustral ncdf transect data 2018, 2021, 2023 : 
#   1) filter by channel, keep only 200kHz data
#   2) cut data on depth dimension to remove most NA
#   3) filter lon >40°E and lat <-30°N
#   4) filter day/night 
# => save intermediary rds file : 2018_day.rds, 2018_night.rds, 2021_day.rds,...

# 02_NASC_concat_sv_diurnal_period.R
#   5) Concatenate each year of data
# => save intermediary rds file : sv_2018_2021_2022_2023.rds

# 03_NASC_check_pigmeann_grid.R
#   6) get pigmeann grid (check if pigmeann grid consistent on every year of data)

# 04_NASC_mean_profile_pig_grid.R
#   7) mean profile by grid, date and diurnal period
# => save intermediary rds file : mean_pig_grid_2018_2021_2023.rds, ...

# 05_NASC_mean_pig_grid
#   8) compute NASC for each mean profile 
# => save intermediary rds file : nasc_mean_pig_grid_2018_2021_2023_day.rds,..



# packages
library(ggplot2)
library(patchwork)

# Global Variables
rm(list=ls())
freqs <- c(18 , 38, 70, 120, 200) #
grid_output_dir <- "F:/data_elise/sv_cropped/grids_custom"
grid_paths <- list.files(grid_output_dir, pattern = "^pigmeann_grid_.*\\.rds$", full.names = TRUE)

fig_base_dir <- "C:/Users/mmolinet/elisou_ta_stagiaire_pref/MIO_internship_III/figures/Sv_profiles"

# ============================================================
# Helper : dossier + libelle de periode selon idx (1=night, 2=day_and_night, 3=day)
# ============================================================
period_label <- function(idx) c("night", "day_and_night", "day")[idx]

for (freq in freqs) {
  path_sv <- paste0("F:/data_elise/sv_cropped/sv_cropped_all_years/Sv_2018_2021_2022_2023_", freq, "kHz.rds")
  sv <- readRDS(path_sv)
  
  for (path_pig_grid in grid_paths) {
    
    grid_name <- tools::file_path_sans_ext(basename(path_pig_grid))  # ex: "pigmeann_grid_lon1500_lat1000"
    pig_grid <- readRDS(path_pig_grid)
  
  
  
  
  ####################################################################### profil moyen par grille
  # --- Résolution de la grille ---
  lon_res <- diff(pig_grid$lon)  # pas en longitude (°)
  lat_res <- diff(pig_grid$lat)  # pas en latitude (°)
  
  # vérifier si la grille est régulière
  summary(lon_res)
  summary(lat_res)
  
  # pas moyen (au cas où il y a de petites variations d'arrondi)
  dlon <- mean(lon_res)
  dlat <- mean(lat_res)
  
  cat("Résolution grille : ", round(dlon, 4), "° lon x ", round(dlat, 4), "° lat\n")
  
  
  # -------------- profil moyen par grid pigmeann et par jour
  sv_date <- as.Date(sv$time, origin = "1950-01-01")  # à adapter si origine différente
  table(sv_date)
  sv_date <- format(sv_date, "%Y-%m-%d")
  sv_day <- sv$day
  
  # Hidtogramme de nombre de profils par jour
  for (idx in c(1, 2, 3)){
    if (idx != 2) {
      daily_counts <- as.data.frame(table(sv_date[sv_day==idx]))
    }
    else {
      daily_counts <- as.data.frame(table(sv_date))
    }
    
    colnames(daily_counts) <- c("date", "n_profils")
    days <- c("night", "day and night", "day")
    
    mean_profils <- mean(daily_counts$n_profils, na.rm = TRUE)
    
    p <- ggplot(daily_counts, aes(x = date, y = n_profils)) +
      geom_col() +
      geom_hline(yintercept = mean_profils, linetype = "dashed", color = "red", linewidth = 0.5) +
      labs(
        x = "Date", y = "ESU number", title = "ESU number per day",
        subtitle = paste(
          "Transect datas from 2018, 2021, 2022, 2023 OC at", freq, "kHz,", days[idx],
          "- grid", grid_name,
          "\nMean ESU per day =", round(mean_profils, 1)
        )
      ) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
    
    print(p)
    
    period_dir <- file.path(fig_base_dir, "nb_profiles_per_day", period_label(idx))
    dir.create(period_dir, recursive = TRUE, showWarnings = FALSE)
    ggsave(
      filename = file.path(period_dir, paste0("nb_profiles_per_day_", grid_name, "_", freq, "kHz.png")),
      plot = p, width = 10, height = 6, dpi = 300, units = "in"
    )
  }
  
  
  # --- 1. Attribution à la cellule de grille (comme avant) ---
  assign_cell <- function(coord, grid_vals) {
    grid_sorted <- sort(unique(grid_vals))
    idx <- findInterval(coord, grid_sorted, all.inside = TRUE)
    d_before <- abs(coord - grid_sorted[idx])
    d_after  <- abs(coord - grid_sorted[pmin(idx + 1, length(grid_sorted))])
    idx_final <- ifelse(d_after < d_before, idx + 1, idx)
    grid_sorted[idx_final]
  }
  
  lon_cell <- assign_cell(sv$lon, pig_grid$lon)
  lat_cell <- assign_cell(sv$lat, pig_grid$lat)
  
  lon_cell <- assign_cell(sv$lon, pig_grid$lon)
  lat_cell <- assign_cell(sv$lat, pig_grid$lat)
  
  # ============================================================
  # DIAGNOSTIC 1 : transect vs grille (localisation spatiale)
  # ------------------------------------------------------------
  # pig_grid ne contient que les centres de cellule -> on reconstruit
  # les bordures (a mi-chemin entre 2 centres consecutifs) pour tracer
  # le quadrillage.
  # ============================================================
  
  lon_centers_g <- sort(unique(pig_grid$lon))
  lat_centers_g <- sort(unique(pig_grid$lat))
  
  lon_edges <- c(
    lon_centers_g[1] - diff(lon_centers_g)[1] / 2,
    head(lon_centers_g, -1) + diff(lon_centers_g) / 2,
    tail(lon_centers_g, 1) + tail(diff(lon_centers_g), 1) / 2
  )
  lat_edges <- c(
    lat_centers_g[1] - diff(lat_centers_g)[1] / 2,
    head(lat_centers_g, -1) + diff(lat_centers_g) / 2,
    tail(lat_centers_g, 1) + tail(diff(lat_centers_g), 1) / 2
  )
  
  p_transect_grid <- ggplot() +
    geom_point(aes(x = lon, y = lat), data = data.frame(lon = sv$lon, lat = sv$lat), size = 0.2, alpha = 0.15, color = "grey40") +
    geom_vline(xintercept = lon_edges, color = "steelblue", linewidth = 0.3) +
    geom_hline(yintercept = lat_edges, color = "steelblue", linewidth = 0.3) +
    geom_point(aes(x = lon, y = lat), data = expand.grid(lon = lon_centers_g, lat = lat_centers_g), color = "firebrick", size = 1) +
    labs(title = paste("Transect vs grille -", grid_name, "-", freq, "kHz"), x = "Longitude", y = "Latitude") +
    coord_fixed() +
    theme_minimal()
  
  print(p_transect_grid)
  
  ggsave(
    filename = paste0("F:/data_elise/sv_cropped/grids_custom/diag_transect_grid_", grid_name, "_", freq, "kHz.png"),
    plot = p_transect_grid, width = 10, height = 8, dpi = 300, units = "in"
  )
  
  # ============================================================
  # DIAGNOSTIC 2 : repartition du nombre de points PAR CELLULE
  # ------------------------------------------------------------
  # Agrege sur toutes les dates/jours confondus (contrairement a
  # n_profils qui est par cellule x date x jour) -> montre si les
  # points bruts se concentrent sur une seule cellule de grille ou
  # sont bien repartis spatialement.
  # ============================================================
  
  cell_id_all <- paste(lon_cell, lat_cell, sep = "|")
  counts_per_cell <- as.data.frame(table(cell_id_all))
  colnames(counts_per_cell) <- c("cell_id", "n_points")
  
  n_cells_total <- length(pig_grid$lon) * length(pig_grid$lat)
  n_cells_occupied <- nrow(counts_per_cell)
  cv_points <- sd(counts_per_cell$n_points) / mean(counts_per_cell$n_points)
  
  cat(sprintf(
    "Grille %s (%d kHz) : %d/%d cellules occupees (%.1f%%)\n",
    grid_name, freq, n_cells_occupied, n_cells_total, 100 * n_cells_occupied / n_cells_total
  ))
  cat(sprintf(
    "  Points/cellule : moyenne=%.0f, ecart-type=%.0f, variance=%.0f, CV=%.2f, min=%d, max=%d\n",
    mean(counts_per_cell$n_points), sd(counts_per_cell$n_points), var(counts_per_cell$n_points),
    cv_points, min(counts_per_cell$n_points), max(counts_per_cell$n_points)
  ))
  
  p_hist_counts <- ggplot(counts_per_cell, aes(x = n_points)) +
    geom_histogram(bins = 30, fill = "steelblue", color = "white") +
    labs(
      title = paste("Points par cellule -", grid_name, "-", freq, "kHz"),
      subtitle = sprintf("%d/%d cellules occupees | moyenne=%.0f | CV=%.2f (0=uniforme, >1=tres concentre)",
                         n_cells_occupied, n_cells_total, mean(counts_per_cell$n_points), cv_points),
      x = "Nombre de pings par cellule", y = "Nombre de cellules"
    ) +
    theme_minimal()
  
  print(p_hist_counts)
  
  ggsave(
    filename = paste0("F:/data_elise/sv_cropped/grids_custom/diag_hist_points_par_cellule_", grid_name, "_", freq, "kHz.png"),
    plot = p_hist_counts, width = 8, height = 6, dpi = 300, units = "in"
  )
  
  
  # clé de groupement incluant maintenant l'année
  cell_id_day <- paste(
    lon_cell,
    lat_cell,
    sv_date,
    sv$day,
    sep = "|"
  )
  
  # --- 2. Conversion dB -> linéaire ---
  sv_linear <- 10^(sv$profiles / 10)
  
  # --- 3. Moyenne par cellule ET par année, robuste aux NA ---
  sum_by_cell   <- rowsum(sv_linear, group = cell_id_day, na.rm = TRUE)
  count_by_cell <- rowsum(
    matrix(as.numeric(!is.na(sv_linear)), 
           nrow = nrow(sv_linear), 
           ncol = ncol(sv_linear)),
    group = cell_id_day
  )
  
  # Moyenne linéaire
  mean_profiles_linear <- sum_by_cell / count_by_cell
  mean_profiles_linear[count_by_cell == 0] <- NA
  
  # --- 4. Reconversion en dB ---
  mean_profiles_db <- 10 * log10(mean_profiles_linear)
  
  # --- 5. récupérer pour chaque profil moyen : lon,lat et time moyen + nb de profils qui ont servi à calculer la moyenne.
  group_levels <- unique(cell_id_day)
  
  group_info <- do.call(rbind, strsplit(group_levels, "\\|"))
  colnames(group_info) <- c(
    "lon",
    "lat",
    "date",
    "day"
  )
  
  head(group_info)
  
  # Temps moyen
  mean_time_by_cell <- tapply(
    as.numeric(sv$time),
    cell_id_day,
    mean,
    na.rm = TRUE
  )
  
  # Reconversion en POSIXct
  mean_time_by_cell <- as.POSIXct(
    mean_time_by_cell,
    origin = "1970-01-01",
    tz = "UTC"
  )
  
  print(mean_time_by_cell)
  
  # Objet final
  mean_pig_grid_by_day <- list(
    profiles  = mean_profiles_db,
    lon       = as.numeric(group_info[, 1]),
    lat       = as.numeric(group_info[, 2]),
    date      = as.Date(group_info[, 3]),
    day = as.numeric(group_info[, "day"]),
    time      =  mean_time_by_cell,
    depth     = sv$depth,
    n_profils = as.vector(table(cell_id_day))
  )
  
  str(mean_pig_grid_by_day, max.level = 1) # verif du dataset final
  
  # --- Sauvegarde ---
  
  # fichier global avec toutes les années (time moyen conservé)
  saveRDS(
    mean_pig_grid_by_day,
    paste0("F:/data_elise/sv_cropped/sv_cropped_all_years/mean_Sv_2018_2021_2022_2023_", grid_name, "_", freq, "kHz.rds")
  )
  
  # --- Histogramme du nombre d'ESU post moyennage
  mean_pig_grid_df <- data.frame(
    lon       = mean_pig_grid_by_day$lon,
    lat       = mean_pig_grid_by_day$lat,
    date      = mean_pig_grid_by_day$date,
    day       = mean_pig_grid_by_day$day,
    time      = mean_pig_grid_by_day$time,
    n_profils = mean_pig_grid_by_day$n_profils
  )
  
  # Ajouter les profondeurs comme colonnes
  profiles_df <- as.data.frame(mean_pig_grid_by_day$profiles)
  
  colnames(profiles_df) <- paste0(
    "depth_",
    mean_pig_grid_by_day$depth
  )
  
  mean_pig_grid_df <- cbind(
    mean_pig_grid_df,
    profiles_df
  )
  
  str(mean_pig_grid_df)
  
  for(idx in c(1, 2, 3)){
    if (idx != 2){
      df_day <- mean_pig_grid_df[mean_pig_grid_df$day == idx,]
    }
    else{
      df_day <- mean_pig_grid_df
    }
    
    depth_cols <- grep("^depth_", names(df_day), value = TRUE)
    
    df_day$valid_profile <- apply(df_day[, depth_cols], 1, function(x) any(!is.na(x)))
    
    daily_counts <- as.data.frame(table(df_day$date[df_day$valid_profile]))
    colnames(daily_counts) <- c("date", "n_ESU")
    
    mean_ESU <- mean(daily_counts$n_ESU, na.rm = TRUE)
    
    p <- ggplot(daily_counts, aes(x = date, y = n_ESU)) +
      geom_col() +
      geom_hline(yintercept = mean_ESU, linetype = "dashed", color = "red", linewidth = 0.5) +
      labs(
        x = "Date", y = "ESU number",
        title = paste("ESU number per day after computation of mean (by pigment grid", grid_name , ", date and diurnal period)"),
        subtitle = paste(
          "Transect datas from 2018, 2021, 2022, 2023 OC at", freq, "kHz,", days[idx],
          "\nMean ESU per day =", round(mean_ESU, 1)
        )
      ) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
    
    print(p)
    
    period_dir <- file.path(fig_base_dir, "nb_profiles_per_day", period_label(idx))
    dir.create(period_dir, recursive = TRUE, showWarnings = FALSE)
    ggsave(
      filename = file.path(period_dir, paste0("post_mean_ESU_per_day_", grid_name, "_", freq, "kHz.png")),
      plot = p, width = 10, height = 6, dpi = 300, units = "in"
    )
  }
  }
}
 