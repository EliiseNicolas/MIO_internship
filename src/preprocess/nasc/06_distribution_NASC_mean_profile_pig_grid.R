# Description 

# Compute NASC distribution for each year of data

# Global variable
path_nasc <- "~/Documents/stage_MIO/pt_III/data_preprocessed/NASC/transect_2018_2022_2023/NASC_mean_pig_grid_by_year_2018_2021_2023_day_200kHz.rds"

# Open nasc ds
nasc_ds <- readRDS(path_nasc)

# Global distribution (2018, 2021, 2023)
hist(
  log(nasc_ds$NASC),
  breaks = 200,
  main = "Distribution du log(NASC) données transect 200kHz 2018-2021-2023",
  xlab = "log(NASC)",
  ylab = "Fréquence"
)

# distribution par année
year <- as.numeric(format(nasc_ds$time, "%Y"))
years <- c(2018, 2021, 2023)

for (y in years){
  mask <- year == y
  hist(
    log(nasc_ds$NASC[mask]),
    breaks = 200,
    main = paste("Distribution du log(NASC) données transect 200kHz", y),
    xlab = "log(NASC)",
    ylab = "Fréquence"
  )
}
