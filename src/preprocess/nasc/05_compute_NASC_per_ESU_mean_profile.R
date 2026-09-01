# global variables
freqs <- c(18, 38, 70, 120, 200)

fig_output_dir <- "C:/Users/mmolinet/elisou_ta_stagiaire_pref/MIO_internship_III/figures/diag_NASC"
dir.create(fig_output_dir, recursive = TRUE, showWarnings = FALSE)

library(ggplot2)

for (freq in freqs){
  
  path_mean_profile <- paste0("F:/data_elise/sv_cropped/sv_cropped_all_years/Sv_2018_2021_2022_2023_", freq, "kHz.rds") # par ESU
  
  # Ouvrir le fichier
  mean_profiles <- readRDS(path_mean_profile)
  Sv <- mean_profiles$profiles
  str(mean_profiles)
  print(unique(format(as.Date(mean_profiles$time), "%Y")))
  
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
  mean_profiles$lon      <- mean_profiles$lon[mask_keep]
  mean_profiles$lat      <- mean_profiles$lat[mask_keep]
  mean_profiles$time     <- mean_profiles$time[mask_keep]
  mean_profiles$day      <- mean_profiles$day[mask_keep]
  
  Sv <- mean_profiles$profiles
  
  str(mean_profiles)
  
  # ------------------------------------------------------------
  # Où sont les NA restants ?
  # ------------------------------------------------------------
  
  years <- c("2018", "2021", "2022", "2023")
  
  # date dérivée de time (le champ 'date' n'existe plus dans la nouvelle structure)
  mean_profiles$year <- format(as.Date(mean_profiles$time), "%Y")
  
  na_by_year_depth <- list()
  
  for (year in years){
    
    idx <- mean_profiles$year == year
    profiles_year <- mean_profiles$profiles[idx, , drop = FALSE]
    
    if (nrow(profiles_year) == 0) next
    
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
    labs(
      x = "Depth (m)", y = "Profiles with NA (%)", color = "Year",
      title = "Remaining NA by depth",
      subtitle = paste("After removing profiles with >50% NA -", freq, "kHz ")
    ) +
    theme_minimal()
  
  print(p_na)
  
  ggsave(
    filename = file.path(fig_output_dir, paste0("diag_NA_by_depth_per_esu_", freq, "kHz.png")),
    plot = p_na, width = 10, height = 6, dpi = 300, units = "in"
  )
  
  # ------------------------------------------------------------
  # Interpoler les NA restants + extrapoler début/fin
  # ------------------------------------------------------------
  
  cat("NA avant interpolation :", sum(is.na(mean_profiles$profiles)), "\n")
  
  depth <- mean_profiles$depth
  
  for (r in seq_len(nrow(mean_profiles$profiles))){
    
    y <- mean_profiles$profiles[r, ]
    ok <- !is.na(y)
    
    if (sum(ok) >= 2){
      mean_profiles$profiles[r, !ok] <- approx(x = depth[ok], y = y[ok], xout = depth[!ok], method = "linear", rule = 2)$y
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
  # n_profils retiré : n'existe plus dans la nouvelle structure
  # (chaque ligne est un ESU individuel, plus un profil moyen agrégé)
  
  nasc_df <- data.frame(
    time = mean_profiles$time,
    lat  = mean_profiles$lat,
    lon  = mean_profiles$lon,
    day  = mean_profiles$day,
    NASC = NASC
  )
  
  print(length(nasc_df$NASC))
  str(nasc_df)
  
  # ------------------------------------------------------------
  # Sauvegarde
  # ------------------------------------------------------------
  
  saveRDS(nasc_df, paste0("F:/data_elise/NASC/NASC_all_esu/NASC_per_esu_2018_2022_2021_2023_", freq, "kHz.rds"))
  
  # ------------------------------------------------------------
  # Plot du NASC au cours du temps par année
  # ------------------------------------------------------------
  
  nasc_df$year <- format(as.Date(nasc_df$time), "%Y")
  
  for (year in years){
    
    df_year <- nasc_df[nasc_df$year == year, ]
    
    if (nrow(df_year) == 0) next
    
    p <- ggplot(df_year, aes(x = time, y = log(NASC))) +
      geom_point(size = 1) +
      labs(
        x = "Time", y = "log(NASC)",
        title = paste("NASC over time in", year),
        subtitle = paste("Transect Datas,", freq, "kHz "),
        caption = "Profiles with >50% NAs removed, linear interpolation + extrapolation"
      ) +
      theme_minimal() +
      theme(plot.caption = element_text(hjust = 0))
    
    print(p)
    
    ggsave(
      filename = file.path(fig_output_dir, paste0("diag_NASC_time_per_esu_", freq, "kHz_", year, ".png")),
      plot = p, width = 10, height = 6, dpi = 300, units = "in"
    )
  }
}