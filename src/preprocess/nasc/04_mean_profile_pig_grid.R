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
freqs <- c(18, 38, 70, 120, 200)
for (freq in freqs){
  path_sv <-paste0("F:/data_elise/sv_cropped/sv_cropped_all_years/Sv_2018_2021_2023_", freq, "kHz.rds")
  path_pig_grid <- "F:/data_elise/sv_cropped/pigmeann_grid.rds"
  
  pig_grid <- readRDS(path_pig_grid)
  sv <- readRDS(path_sv)
  
  
  
  
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
  
  
  # -------------- profil moyen par grid pigmeann et par année
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
    
    # Renommer les colonnes
    colnames(daily_counts) <- c("date", "n_profils")
    
    days <- c("night", "day and night", "day")
    
    # Histogramme
    mean_profils <- mean(
      daily_counts$n_profils,
      na.rm = TRUE
    )
    
    # Histogramme
    p <- ggplot(
      daily_counts,
      aes(x = date, y = n_profils)
    ) +
      geom_col() +
      geom_hline(
        yintercept = mean_profils,
        linetype = "dashed",
        color = "red",
        linewidth = 0.5
      ) +
      labs(
        x = "Date",
        y = "ESU number",
        title = "ESU number per day",
        subtitle = paste(
          "Transect datas from 2018, 2021, 2023 OC at",
          freq,
          "kHz,",
          days[idx],
          "\nMean ESU per day =",
          round(mean_profils, 1)
        )
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(
          angle = 90,
          vjust = 0.5
        )
      )
    
    print(p)
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
    paste0("F:/data_elise/sv_cropped/sv_cropped_all_years/mean_Sv_pig_grid_by_date_2018_2021_2023_", freq, "kHz.rds")
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
    
    depth_cols <- grep(
      "^depth_",
      names(df_day),
      value = TRUE
    )
    
    df_day$valid_profile <- apply(
      df_day[, depth_cols],
      1,
      function(x) any(!is.na(x))
    )
    
    daily_counts <- as.data.frame(
      table(
        df_day$date[df_day$valid_profile]
      )
    )
    
    colnames(daily_counts) <- c(
      "date",
      "n_ESU"
    )
    
    mean_ESU <- mean(daily_counts$n_ESU, na.rm = TRUE)
    
    p <- ggplot(daily_counts, aes(x = date, y = n_ESU)) +
      geom_col() +
      geom_hline(
        yintercept = mean_ESU,
        linetype = "dashed",
        color = "red",
        linewidth = 0.5
      ) +
      labs(
        x = "Date",
        y = "ESU number",
        title = "ESU number per day after computation of mean (by pigment grid, date and diurnal period)",
        subtitle = paste(
          "Transect datas from 2018, 2021, 2023 OC at",
          freq,
          "kHz,",
          days[idx],
          "\nMean ESU per day =",
          round(mean_ESU, 1)
        )
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(
          angle = 90,
          vjust = 0.5
        )
      )
    
    print(p)
  }
}








# # ----------- Plot des effets du moyennage
# 
# # ------------------------------------------------------------
# # 1. Variance par cellule × jour et par profondeur
# #    sur les valeurs linéaires
# # ------------------------------------------------------------
# 
# sum_sq_by_cell <- rowsum(
#   sv_linear^2,
#   group = cell_id_day,
#   na.rm = TRUE
# )
# 
# var_profiles_linear <- (
#   sum_sq_by_cell -
#     count_by_cell * mean_profiles_linear^2
# ) / (count_by_cell - 1)
# 
# # Pas de variance si moins de 2 profils
# var_profiles_linear[count_by_cell <= 1] <- NA
# 
# 
# # ------------------------------------------------------------
# # 2. Data frame récapitulatif par cellule × jour
# # ------------------------------------------------------------
# 
# n_profils_df <- data.frame(
#   cell_id   = rownames(mean_profiles_db),
#   n_profils = mean_pig_grid_by_day$n_profils,
#   date      = mean_pig_grid_by_day$date
# )
# 
# # Vérification
# head(n_profils_df)
# 
# 
# # ------------------------------------------------------------
# # 3. Statistiques sur le nombre de profils
# # ------------------------------------------------------------
# 
# mean_n <- mean(
#   n_profils_df$n_profils,
#   na.rm = TRUE
# )
# 
# sd_n <- sd(
#   n_profils_df$n_profils,
#   na.rm = TRUE
# )
# 
# 
# # ------------------------------------------------------------
# # 4. Statistiques sur la variance
# # ------------------------------------------------------------
# 
# mean_variance <- mean(
#   var_profiles_linear,
#   na.rm = TRUE
# )
# 
# sd_variance <- sd(
#   var_profiles_linear,
#   na.rm = TRUE
# )
# 
# # PLOTS
# 
# p1 <- ggplot(
#   n_profils_df,
#   aes(x = n_profils)
# ) +
#   geom_histogram(
#     bins = 30,
#     fill = "steelblue",
#     color = "white"
#   ) +
#   geom_vline(
#     xintercept = mean_n,
#     color = "red",
#     linetype = "dashed",
#     linewidth = 0.6
#   ) +
#   labs(
#     title = "Nombre de profils par cellule et par jour",
#     subtitle = paste0(
#       "Moyenne = ", round(mean_n, 1),
#       " ± ", round(sd_n, 1)
#     ),
#     x = "Nombre de profils",
#     y = "Nombre de cellules × jours"
#   ) +
#   theme_minimal()
# 
# # ------------------------------------------------------------
# # Data frame pour le plot de variance
# # ------------------------------------------------------------
# 
# # ------------------------------------------------------------
# # Data frame pour le plot de variance
# # ------------------------------------------------------------
# 
# variance_values <- as.vector(var_profiles_linear)
# 
# var_df <- data.frame(
#   variance = variance_values
# )
# 
# # Vérification
# class(var_df)
# str(var_df)
# dim(var_df)
# 
# # Retirer les NA et les valeurs <= 0
# var_df <- var_df[
#   is.finite(var_df$variance) &
#     var_df$variance > 0,
#   ,
#   drop = FALSE
# ]
# 
# p2 <- ggplot(
#   var_df,
#   aes(x = variance)
# ) +
#   geom_histogram(
#     bins = 40,
#     fill = "darkorange",
#     color = "white"
#   ) +
#   geom_vline(
#     xintercept = mean_variance,
#     color = "red",
#     linetype = "dashed",
#     linewidth = 0.6
#   ) +
#   scale_x_log10() +
#   labs(
#     title = "Variance du moyennage",
#     subtitle = paste0(
#       "Moyenne = ", signif(mean_variance, 3)
#     ),
#     x = "Variance (échelle linéaire, log10)",
#     y = "Fréquence"
#   ) +
#   theme_minimal()
# 
# p1 / p2


# ################################################################### plots des transects sur la grille pigmeann
# # construire un data.frame propre à partir des coordonnées réelles
# df_sv <- data.frame(
#   lon = sv$lon,
#   lat = sv$lat, 
#   time = sv$time
# )
# 
# 
# # retirer les NA éventuels
# df_sv <- df_sv[!is.na(df_sv$lon) & !is.na(df_sv$lat), ]
# df_sv$date <- as.Date(df_sv$time, origin = "1950-01-01")
# df_sv$year <- factor(format(df_sv$date, "%Y"))
# 
# 
# # plots 
# ggplot() +
#   
#   # grille PIGMeANN
#   geom_vline(
#     xintercept = pig_grid$lon,
#     color = "grey80",
#     linewidth = 0.15
#   ) +
#   geom_hline(
#     yintercept = pig_grid$lat,
#     color = "grey80",
#     linewidth = 0.15
#   ) +
#   
#   # points Sv colorés par année
#   geom_point(
#     data = df_sv,
#     aes(x = lon, y = lat, color = year),
#     alpha = 0.3,
#     size = 0.4
#   ) +
#   
#   scale_color_manual(
#     values = c("2018" = "#1b9e77", "2021" = "#d95f02", "2023" = "#7570b3")
#   ) +
#   
#   guides(color = guide_legend(override.aes = list(alpha = 1, size = 2))) +
#   
#   coord_fixed() +
#   
#   labs(
#     title = paste0("Points Sv sur la grille PIGMeANN - ", freq, "kHz"),
#     x = "Longitude (°E)",
#     y = "Latitude (°)",
#     color = "Année"
#   ) +
#   
#   theme_minimal()
 