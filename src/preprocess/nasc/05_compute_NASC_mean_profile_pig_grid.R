# Description

# from mean profile per pigemann grid point, compute NASC
# Save NASC in rds file containing 4 columns : time, lat lon NASC

# rm(list=ls())

# library

# global variables
freqs <- c(18, 38, 70, 120, 200)

######################################################################################## Code pour les Sv moyennés
# Code pour le NASC sur chaque profil moyen par grille pigment/date/période diurne
for (freq in freqs){
  
  path_mean_profile <- paste0("F:/data_elise/sv_cropped/sv_cropped_all_years/mean_Sv_pig_grid_by_date_2018_2021_2023_", freq, "kHz.rds")
  
  # Ouvrir le fichier
  mean_profiles <- readRDS(path_mean_profile)
  Sv <- mean_profiles$profiles
  str(mean_profiles)
  
  # ------------------------------------------------------------
  # Diagnostic des NA
  # ------------------------------------------------------------
  
  # Profils entièrement NA
  pct_all_na <- 100 * sum(rowSums(is.na(Sv)) == ncol(Sv)) / nrow(Sv)
  cat("Profils entièrement NA :", round(pct_all_na, 2), "%\n")
  
  # Profils avec >20% de NA
  pct_more_20 <- 100 * sum(rowSums(is.na(Sv)) > 0.2 * ncol(Sv)) / nrow(Sv)
  cat("Profils avec >20% de NA :", round(pct_more_20, 2), "%\n")
  
  # ------------------------------------------------------------
  # Retirer les profils qui ont >50% de NA
  # ------------------------------------------------------------
  
  mask_keep <- rowSums(is.na(Sv)) <= 0.5 * ncol(Sv)
  
  cat("Profils supprimés :", sum(!mask_keep), "sur", nrow(Sv), "\n")
  
  mean_profiles$profiles <- mean_profiles$profiles[mask_keep, , drop = FALSE]
  mean_profiles$lon <- mean_profiles$lon[mask_keep]
  mean_profiles$lat <- mean_profiles$lat[mask_keep]
  mean_profiles$date <- mean_profiles$date[mask_keep]
  mean_profiles$day <- mean_profiles$day[mask_keep]
  mean_profiles$time <- mean_profiles$time[mask_keep]
  mean_profiles$n_profils <- mean_profiles$n_profils[mask_keep]
  
  Sv <- mean_profiles$profiles
  
  str(mean_profiles)
  
  # ------------------------------------------------------------
  # Où sont les NA restants ? 
  # ------------------------------------------------------------
  
  years <- c("2018", "2021", "2023")
  
  mean_profiles$year <- format(as.Date(mean_profiles$date), "%Y")
  
  na_by_year_depth <- list()
  
  for (year in years){
    
    idx <- mean_profiles$year == year
    profiles_year <- mean_profiles$profiles[idx, , drop = FALSE]
    
    # Nombre de NA à chaque profondeur
    n_na <- colSums(is.na(profiles_year))
    
    # Pourcentage de profils avec NA
    pct_na <- 100 * n_na / nrow(profiles_year)
    
    na_by_year_depth[[year]] <- data.frame(depth = mean_profiles$depth, n_NA = n_na, pct_NA = pct_na, n_profiles = nrow(profiles_year), year = year)
  }
  
  na_df <- do.call(rbind, na_by_year_depth)
  rownames(na_df) <- NULL
  
  # Plot des NA par profondeur
  p_na <- ggplot(na_df, aes(x = depth, y = pct_NA, color = year)) +
    geom_line() +
    geom_point(size = 0.5) +
    labs(x = "Depth (m)", y = "Profiles with NA (%)", color = "Year", title = "Remaining NA by depth", subtitle = paste("After removing profiles with >50% NA -", freq, "kHz")) +
    theme_minimal()

  print(p_na)
  
  # ------------------------------------------------------------
  # Interpoler les NA restants + extrapoler début/fin
  # ------------------------------------------------------------
  
  cat("NA avant interpolation :", sum(is.na(mean_profiles$profiles)), "\n")
  
  depth <- mean_profiles$depth
  
  for (i in seq_len(nrow(mean_profiles$profiles))){
    
    y <- mean_profiles$profiles[i, ]
    ok <- !is.na(y)
    
    if (sum(ok) >= 2){
      mean_profiles$profiles[i, !ok] <- approx(x = depth[ok], y = y[ok], xout = depth[!ok], method = "linear", rule = 2)$y
    }
  }
  
  Sv <- mean_profiles$profiles
  
  cat("NA après interpolation :", sum(is.na(Sv)), "\n")
  
  # ------------------------------------------------------------
  # Calcul du NASC
  # ------------------------------------------------------------
  
  sv <- 10^(Sv / 10)
  print(dim(sv))
  
  int <- rowSums(sv, na.rm = TRUE)
  print(length(int))
  
  depth_step <- mean(diff(mean_profiles$depth))
  print(depth_step)
  
  sa <- int * depth_step
  
  NASC <- 4 * pi * 1852^2 * sa
  
  print(summary(NASC))
  
  # ------------------------------------------------------------
  # Data frame final
  # ------------------------------------------------------------
  
  nasc_df <- data.frame(time = mean_profiles$time, lat = mean_profiles$lat, lon = mean_profiles$lon, day = mean_profiles$day, NASC = NASC, n_profils = mean_profiles$n_profils)
  
  print(length(nasc_df$NASC))
  str(nasc_df)
  
  # ------------------------------------------------------------
  # Sauvegarde
  # ------------------------------------------------------------
  
  saveRDS(nasc_df, paste0("F:/data_elise/NASC/NASC_pig_mean/mean_Sv_pig_grid_by_date_2018_2021_2023_", freq, "kHz.rds"))
  
  
  # ------------------------------------------------------------
  # Plot du NASC au cours du temps par année
  # ------------------------------------------------------------
  
  years <- c("2018", "2021", "2023")
  nasc_df$year <- format(as.Date(nasc_df$time), "%Y")
  
  for (year in years){
    
    df_year <- nasc_df[nasc_df$year == year, ]
    
    p <- ggplot(df_year, aes(x = time, y = log(NASC))) +
      geom_point(size = 0.2) +
      labs(x = "Time", y = "log(NASC)", title = paste("NASC over time in", year), subtitle = paste("Transect Datas,", freq, "kHz"), caption = "Profiles with >50% NAs removed, linear interpolation + extrapolation") +
      theme_minimal() +
      theme(plot.caption = element_text(hjust = 0))
    
    print(p)
  }
}





################################################################################## Code pour le NASC sur chaque ESU
for (freq in c(18, 38, 70, 120, 200)){
  path_esu <- paste0("F:/data_elise/sv_cropped/sv_cropped_all_years/Sv_2018_2021_2023_", freq, "kHz.rds")
  
  # open file
  esu <- readRDS(path_esu)
  str(esu)
  
  # diagnostic de NASC
  print(sum(rowSums(is.na(esu$profiles)) == ncol(esu$profiles))/nrow(esu$profiles)*100)
  # 200kHz : 1.22% des données entièrement NA
  print(sum(rowSums(is.na(esu$profiles)) > 0.2*ncol(esu$profiles))/nrow(esu$profiles)*100)
  # 200kHz : 25% ont + de 20% de NA
  
  # retirer les profiles qui ont >50% NA
 
  mask_keep <- rowSums(is.na(esu$profiles)) <= 0.5*ncol(esu$profiles)
  cat(
    "Profils supprimés :",
    sum(!mask_keep),
    "sur",
    nrow(esu$profiles),
    "\n"
  )
  
  esu$profiles <- esu$profiles[mask_keep, , drop = FALSE]
  esu$lat <- esu$lat[mask_keep, drop = FALSE]
  esu$lon <- esu$lon[mask_keep, drop = FALSE]
  esu$time <- esu$time[mask_keep, drop = FALSE]
  esu$day <- esu$day[mask_keep, drop = FALSE]
  
  str(esu)
  # ----------------------------------------------- Ou sont les NAs restants ? (quels profondeurs)
  years <- c("2018", "2021", "2023")
  
  # année de chaque profil
  esu$year <- format(as.Date(esu$time), "%Y")
  
  na_by_year_depth <- list()
  
  for (year in years) {
    
    # Profils de l'année
    idx <- esu$year == year
    
    profiles_year <- esu$profiles[idx, , drop = FALSE]
    
    # Nombre de NA à chaque profondeur
    n_na <- colSums(is.na(profiles_year))
    
    # Pourcentage de profils avec NA
    pct_na <- 100 * n_na / nrow(profiles_year)
    
    na_by_year_depth[[year]] <- data.frame(
      depth = esu$depth,
      n_NA = n_na,
      pct_NA = pct_na,
      n_profiles = nrow(profiles_year)
    )
  }
  
  na_df <- do.call(
    rbind,
    lapply(names(na_by_year_depth), function(year) {
      
      df <- na_by_year_depth[[year]]
      df$year <- year
      
      df
    })
  )
  
  rownames(na_df) <- NULL
  
  ggplot(
    na_df,
    aes(x = depth, y = pct_NA, color = year)
  ) +
    geom_line() +
    geom_point(size = 0.5) +
    labs(
      x = "Depth (m)",
      y = "Profiles with NA (%)",
      color = "Year",
      title = "Remaining NA by depth",
      subtitle = paste(
        "After removing profiles with >50% NA -",
        freq, "kHz"
      )
    ) +
    theme_minimal()
  
  #------------------------------------------------  Interpoler les NA restant et extrapoler sur les profondeur début/fin
  cat("NA avant interpolation :", sum(is.na(esu$profiles)), "\n")
  # 200kHz : 22740 NA
  
  depth <- esu$depth
  
  for (i in 1:nrow(esu$profiles)) {
    
    y <- esu$profiles[i, ]
    
    # Indices des valeurs valides
    ok <- !is.na(y)
    
    # Interpolation uniquement s'il y a au moins 2 valeurs valides
    if (sum(ok) >= 2) {
      
      esu$profiles[i, is.na(y)] <- approx(
        x = depth[ok],
        y = y[ok],
        xout = depth[is.na(y)],
        method = "linear",
        rule = 2
      )$y
    }
  }
  
  # Nombre de NA après interpolation
  cat("NA après interpolation :", sum(is.na(esu$profiles)), "\n")
  # 200kHz : 4652 NA
  
  # compute nasc 
  sv <- 10^(esu$profiles/10) # linear sv
  print(dim(sv)) # (time, depth)
  int <- rowSums(sv, na.rm = TRUE)
  print(length(int)) # (time,)
  depth_step <- mean(diff(esu$depth))
  print(depth_step) # résolution de 2m
  sa <- int * depth_step # (time,)
  NASC <- 4 * pi * 1852**2 * sa # (time,) 
  print(NASC)
  str(esu)
  nasc_df <- data.frame(
    time = esu$time,
    lat = esu$lat,
    lon = esu$lon,
    day = esu$day,
    NASC = NASC
  )
  print(length(nasc_df$NASC))
  str(nasc_df)
  saveRDS(nasc_df, paste0("F:/data_elise/NASC/NASC_all_ESU/NASC_per_ESU_2018_2021_2023_", freq, "kHz.rds"))
  
  
  # plot de NASC over time
  years <- c("2018", "2021", "2023")
  nasc_df$year <- format(as.Date(nasc_df$time), "%Y")
  sum(nasc_df$NASC == 0, na.rm = TRUE)
  
  sum(nasc_df$NASC == 0, na.rm = TRUE)
  
  sum(
    rowSums(is.na(esu$profiles)) == ncol(esu$profiles)
  )
  for (year in years) {
    
    df_year <- nasc_df[nasc_df$year == year, ]
    
    p <- ggplot(
      df_year,
      aes(x = time, y = log(NASC))
    ) +
      geom_point(size = 0.2) +
      labs(
        x = "Time",
        y = "log(NASC)",
        title = paste("NASC over time in", year),
        subtitle = paste(
          "Transect Datas,", freq, "kHz"
        ),
        caption = "Profiles with >50% NAs removed, linear interpolation + extrapolation"
      ) +
      theme_minimal() +
      theme(
        plot.caption = element_text(hjust = 0)
      )
    
    print(p)
  }
}




# Note sur les NASC par ESU et moyennés

# Les NASCs calculés sur la moyenne ont l'avantage de présenter moins de valeurs manquantes.
# Les NASCs calculés sur chaque ESU independante peuvent être intégrés sur des profondeurs différentes selon les valeurs manquantes.
# Si on retire les profils ayant >20% de NA, alors on supprime presque toutes les données 2021...
# J'ai donc supprimé tous les profiles ayant >50% de NA et ai fait une interpolation linéaire. Les NAs en fin de profil ne peuvent pas être complétés.
# 
# Une solution serait de prendre la dernière valeur et de l'appliquer sur le reste des profondeurs manquantes. c'est mieux que de laisser vide.
