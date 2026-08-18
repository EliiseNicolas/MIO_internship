# Description 

# Compute NASC distribution for each year of data

# Global variable
path_nasc <- "F:/data_elise/NASC/120kHz/NASC_mean_pig_grid_by_year_2018_2021_2023_day_120kHz.rds"

# Open nasc ds
nasc_ds <- readRDS(path_nasc)
x <- log(nasc_ds$NASC[nasc_ds$NASC > 0])
moyenne <- mean(x, na.rm = TRUE)
print(moyenne)
# Global distribution (2018, 2021, 2023)
hist(
  log(nasc_ds$NASC),
  breaks = 200,
  main = "Distribution du log(NASC) données transect 120kHz day 2018-2021-2023",
  xlab = "log(NASC)",
  ylab = "Fréquence"
)
abline(v = moyenne, col = "red", lwd = 2)

# distribution par année
year <- as.numeric(format(nasc_ds$time, "%Y"))
years <- c(2018, 2021, 2023)

for (y in years){
  mask <- year == y & nasc_ds$NASC > 0
  x <- log(nasc_ds$NASC[mask])
  moyenne <- mean(x, na.rm = TRUE)
  
  hist(
    x,
    breaks = 200,
    main = paste("Distribution du log(NASC) données transect 120kHz day", y),
    xlab = "log(NASC)",
    ylab = "Fréquence"
  )
  
  # Trait vertical rouge = moyenne
  abline(v = moyenne, col = "red", lwd = 2)
}
