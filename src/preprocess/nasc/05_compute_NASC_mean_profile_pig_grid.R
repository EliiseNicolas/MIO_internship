# Description

# from mean profile per pigemann grid point, compute NASC
# Save NASC in rds file containing 4 columns : time, lat lon NASC

# rm(list=ls())

# library

# global variables
path_mean_profile <- "~/Documents/stage_MIO/pt_III/data_preprocessed/NASC/transect_2018_2022_2023/mean_pig_grid_by_year_2018_2021_2023_day_200kHz.rds"

# open file
mean_profiles <- readRDS(path_mean_profile)
Sv <- mean_profiles$profiles

# compute nasc 
sv <- 10^(Sv/10) # linear sv
print(dim(sv))
int <- rowSums(sv, na.rm = TRUE)
print(length(int))
sa <- int *1 # car on a un step de profondeur de 1m 
NASC <- 4 * pi * 1852**2 * sa # integration sur la profondeur

time_sv <- as.POSIXct(
  mean_profiles$time * 86400,
  origin = "1950-01-01",
  tz = "UTC"
)
print(head(time_sv))

nasc_df <- data.frame(
  time = time_sv,
  lat = mean_profiles$lat,
  lon = mean_profiles$lon,
  NASC = NASC,
  n_profils <- mean_profiles$n_profils
)
print(length(mean_profiles$NASC))
saveRDS(nasc_df, "~/Documents/stage_MIO/pt_III/data_preprocessed/NASC/transect_2018_2022_2023/NASC_mean_pig_grid_by_year_2018_2021_2023_day_200kHz.rds")
