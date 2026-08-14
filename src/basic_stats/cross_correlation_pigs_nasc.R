# Description 

# cross correlation of NASC and pigs values on 2018, 2021, 2023 day transect data

# librairies
library(dplyr)
library(tidyr)
library(ggplot2)

# global variables
path_nasc <- "~/Documents/stage_MIO/pt_III/data_preprocessed/NASC/transect_2018_2022_2023/NASC_mean_pig_grid_by_year_2018_2021_2023_day_200kHz.rds"
path_pigs <- "~/Documents/stage_MIO/pt_III/data_preprocessed/pigmeann/transect/pigs_transect_day_2018_2021_2023/pig_mask9_1d_2018_2021_2023.rds"

pig_vars <- c(
  "Chla", "Per", "But", "Fuco", "Hex",
  "Allo", "Zea", "Chlb", "DvChla"
)

# open ds
nasc_ds <- readRDS(path_nasc)
pig_ds <- readRDS(path_pigs)

# cross correlation
str(nasc_ds)
str(pig_ds)

data_corr <- nasc_ds %>%
  select(time, lat, lon, NASC) %>%
  left_join(
    pig_ds %>%
      select(time, lat_sv, lon_sv, all_of(pig_vars)),
    by = c(
      "time" = "time",
      "lat"  = "lat_sv",
      "lon"  = "lon_sv"
    )
  )
sum(is.finite(pig_ds$Chla))
sum(is.finite(nasc_ds$NASC))
sum(is.finite(data_corr$Chla))
sum(is.finite(data_corr$NASC))

cor_results <- sapply(
  pig_vars,
  function(pig) {
    
    data_tmp <- data_corr %>%
      select(NASC, all_of(pig)) %>%
      filter(
        is.finite(NASC),
        is.finite(.data[[pig]])
      )
    
    cor(
      data_tmp$NASC,
      data_tmp[[pig]],
      method = "pearson"
    )
  }
)

cor_results


cor_df <- data.frame(
  pigment = names(cor_results),
  correlation = as.numeric(cor_results)
)
ggplot(cor_df, aes(x = pigment, y = "NASC", fill = correlation)) +
  geom_tile() +
  geom_text(
    aes(label = round(correlation, 2)),
    size = 5
  ) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-1, 1)
  ) +
  theme_bw() +
  labs(
    title = "Correlation between NASC and pigments (2018, 2021, 2023 transect day data)",
    x = "Pigment",
    y = NULL,
    fill = "Pearson r"
  )

