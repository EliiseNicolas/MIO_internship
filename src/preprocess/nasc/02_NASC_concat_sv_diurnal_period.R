# Description
# 01_NASC_filter_sv
        # From obsaustral ncdf transect data 2018, 2021, 2023 : 
        #   1) filter by channel, keep only 200kHz data
        #   2) cut data on depth dimension to remove most NA
        #   3) filter lon >40°E and lat <-30°N
        #   4) filter day/night 
        # => save intermediary rds file : 2018_day.rds, 2018_night.rds, 2021_day.rds,...

# 02_NASC_concat_sv_diurnal_period
  #   5) Concatenate each year of same period (i.e. 2018_day, 2021_day, 2023_day)
  # => save intermediary rds file : 2018_2021_2022_2023_day.rds, 2018_2021_2022_2023_night.rds

#   6) get pigmeann grid (check if pigmeann grid consistent on every year of data)
#   7) mean profile by grid (for day and night data)
# => save intermediary rds file : mean_pig_grid_2018_2021_2023_day.rds, ...

#   8) compute NASC for each mean profile 
# => save intermediary rds file : nasc_mean_pig_grid_2018_2021_2023_day.rds,..

# Packages
library(dplyr)

# Global Variables
rm(list=ls())

freqs <- c(38, 70, 120, 200)
for (freq in freqs){
  path2018 <- paste0("F:/data_elise/sv_cropped/sv_cropped_per_year/", freq, "kHz/Sv_2018_", freq, "kHz.rds")
  path2021 <- paste0("F:/data_elise/sv_cropped/sv_cropped_per_year/", freq, "kHz/Sv_2021_", freq, "kHz.rds")
  path2023 <- paste0("F:/data_elise/sv_cropped/sv_cropped_per_year/", freq, "kHz/Sv_2023_", freq, "kHz.rds")
  
  # Read data
  df2018 <- readRDS(path2018)
  df2021 <- readRDS(path2021)
  df2023 <- readRDS(path2023)
  
  # Check that depth grid is identical
  stopifnot(
    identical(df2018$depth, df2021$depth),
    identical(df2018$depth, df2023$depth)
  )
  
  # Concatenate Sv profiles
  Sv_all <- rbind(
    df2018$profiles,
    df2021$profiles,
    df2023$profiles
  ) #(n_profiles, depth)
  
  # Concatenate metadata
  metadata_all <- bind_rows(
    data.frame(
      lat = df2018$lat,
      lon = df2018$lon,
      time = df2018$time,
      day = df2018$day
    ),
    data.frame(
      lat = df2021$lat,
      lon = df2021$lon,
      time = df2021$time,
      day = df2021$day
    ),
    data.frame(
      lat = df2023$lat,
      lon = df2023$lon,
      time = df2023$time,
      day = df2023$day
    )
  )
  
  # Check dimensions
  stopifnot(
    nrow(Sv_all) == nrow(metadata_all)
  )
  
  # Depth vector
  depth <- df2018$depth
  
  Sv_concat <- list(
    profiles = Sv_all,
    lat = metadata_all$lat,
    lon = metadata_all$lon,
    time = metadata_all$time,
    day = metadata_all$day,
    depth = depth
  )
  
  str(Sv_concat)
  
  saveRDS(
    Sv_concat,
    file = paste0(
      "F:/data_elise/sv_cropped/sv_cropped_all_years/Sv_2018_2021_2023_", freq, "kHz.rds"
    )
  )
}

