# ============================================================
# DESCRIPTION DU DATASET
# ============================================================

# We want to understand the distribution of each entry variable
# used in the regression tree / random forest,
# and quantify missing values for each variable.

# ------------------------------------------------------------
# Libraries
# ------------------------------------------------------------

rm(list = ls())

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)


# ------------------------------------------------------------
# Path and dataset
# ------------------------------------------------------------

freq <- 120
path <- paste0("F:/data_elise/ds_NASC_pig_ftle_fod/ds_NASC_mean_pig_grid_all/ds_NASC_per_esu_pig_ftle_fod_2018_2021_2022_2023_transect_120kHz_lon1500_lat1000.rds")

datas <- readRDS(path)

print(nrow(datas))
str(datas)

# filter datas day/night
diurnal_period <- 3 # 3 : day, 1: night
dp <- "day"
ds <- datas[datas$day == diurnal_period,]

# ============================================================
# MISSING VALUES
# ============================================================
# Variables utilisées pour le modèle
print(names(datas))
model_vars <- setdiff(
  names(datas),
  c("time_nasc", "lat_nasc", "lon_nasc", "day" ,"lat_fod", "lon_fod", "lat_pig", "lon_pig", "Chla", "Per", "But", "Fuco","Hex","Allo", "Zea","Chlb","DvChla", "lat_ftle", "lon_ftle")
)
print(vars)

# Nombre de valeurs manquantes
nb_missing <- sapply(
  ds[model_vars],
  function(x) sum(!is.finite(x))
)

# Pourcentage de valeurs manquantes
pct_missing <- nb_missing / nrow(ds) * 100

missing_df <- data.frame(
  variable = model_vars,
  missing = nb_missing,
  percentage = pct_missing,
  available = nrow(ds) - nb_missing
)

print(missing_df)


# ============================================================
# FOD
# ============================================================

# FOD est maintenant directement une variable numérique
# et non plus des variables one-hot fod0 ... fod6.

# On transforme FOD en facteur pour obtenir une distribution
# par cluster.

ds_fod <- ds %>%
  mutate(
    fod_cluster = as.character(fod)
  )

# Remplacer les valeurs manquantes par "NA"
ds_fod$fod_cluster[
  is.na(ds_fod$fod_cluster)
] <- "NA"

# Distribution
fod_distribution <- ds_fod %>%
  count(fod_cluster) %>%
  mutate(
    percentage = n / sum(n) * 100
  )

print(fod_distribution)


# ------------------------------------------------------------
# Plot FOD
# ------------------------------------------------------------

ggplot(
  fod_distribution,
  aes(
    x = fod_cluster,
    y = percentage
  )
) +
  geom_col() +
  labs(
    x = "FOD cluster",
    y = "Percentage (%)",
    title = paste("FOD distribution - 2018, 2021, 2023 - Transect -", dp, "-", freq, "kHz")
  ) +
  theme_bw()


# ============================================================
# FTLE
# ============================================================

ggplot(
  ds,
  aes(
    x = ftle,
    y = after_stat(count / sum(count) * 100)
  )
) +
  geom_histogram(
    bins = 200
  ) +
  theme_bw() +
  labs(
    x = "FTLE",
    y = "Percentage (%)",
    title = paste("FTLE distribution - 2018, 2021, 2023 - Transect -", dp, "-", freq, "kHz")
  )


# ============================================================
# PIGMENTS
# ============================================================

pig_vars <- c(
  "Chla",
  "Per",
  "But",
  "Fuco",
  "Hex",
  "Allo",
  "Zea",
  "Chlb",
  "DvChla"
)

pig_vars <- c(
  "Chla_Chla",
  "Per_Chla",
  "But_Chla",
  "Fuco_Chla",
  "Hex_Chla",
  "Allo_Chla",
  "Zea_Chla",
  "Chlb_Chla",
  "DvChla_Chla"
)

pig_vars <- c(
  "Chla_total",
  "Per_total",
  "But_total",
  "Fuco_total",
  "Hex_total",
  "Allo_total",
  "Zea_total",
  "Chlb_total",
  "DvChla_total"
)
# ------------------------------------------------------------
# Long format
# ------------------------------------------------------------

pig_long <- ds %>%
  select(all_of(pig_vars)) %>%
  pivot_longer(
    cols = everything(),
    names_to = "pigment",
    values_to = "value"
  ) %>%
  filter(
    is.finite(value)
  )


# ------------------------------------------------------------
# Histograms
# ------------------------------------------------------------

p1 <- ggplot(
  pig_long,
  aes(x = value)
) +
  geom_histogram(
    bins = 50
  ) +
  facet_wrap(
    ~ pigment,
    scales = "free"
  ) +
  theme_bw() +
  scale_x_continuous(
    n.breaks = 3
  ) +
  labs(
    x = "Concentration",
    y = "Count"
  )


# ------------------------------------------------------------
# Missing values pigments
# ------------------------------------------------------------

# Nombre total de lignes
n_total <- nrow(ds)

# Nombre de lignes avec les 9 pigments disponibles
n_complete <- sum(
  complete.cases(ds[pig_vars])
)

# Nombre de lignes avec au moins une valeur manquante
n_missing <- n_total - n_complete

# Pourcentage de lignes avec au moins une valeur manquante
pct_missing <- n_missing / n_total * 100

# Pourcentage de lignes complètes
pct_complete <- n_complete / n_total * 100


# ------------------------------------------------------------
# Text panel
# ------------------------------------------------------------

p2 <- ggplot() +
  annotate(
    "text",
    x = 0,
    y = 0,
    label = paste0(
      "Total rows :\n",
      n_total,
      
      "\n\nRows with missing values :\n",
      n_missing,
      " (",
      round(pct_missing, 2),
      "%)",
      
      "\n\nRows with no missing values :\n",
      n_complete,
      " (",
      round(pct_complete, 2),
      "%)"
    ),
    size = 4
  ) +
  theme_void()


# ------------------------------------------------------------
# Combined pigment plot
# ------------------------------------------------------------

p1 + p2 +
  plot_layout(
    widths = c(3, 1)
  ) +
  plot_annotation(
    title = paste("Distribution of pigment variables - 2018, 2021, 2023 - Transect -", dp, "-", freq, "kHz"),
    subtitle = "Missing values are excluded from histograms"
  )


# ============================================================
# DAYS WITH PIGMENT DATA
# ============================================================

idx <- is.finite(ds$Chla)

dates_chla <- unique(
  as.Date(ds$time_nasc[idx])
)

print(dates_chla)

cat(
  "Number of days with Chla data :",
  length(dates_chla),
  "\n"
)


# ============================================================
# NASC
# ============================================================
print(range(ds$nasc))
q <- quantile(ds$nasc, probs = c(0.05, 0.95), na.rm = TRUE)
ds_filtered <- ds |>
  dplyr::filter(nasc >= q[1], nasc <= q[2])

# ds_filtered$nasc <- log10(ds_filtered$nasc)
print(range(ds_filtered$nasc))
ggplot(
  ds_filtered,
  aes(
    x = nasc,
    y = after_stat(count / sum(count) * 100)
  )
) +
  geom_histogram(
    bins = 200
  ) +
  # scale_x_log10(
  #   labels = scales::label_number()
  # ) +
  theme_bw() +
  labs(
    x = "NASC",
    y = "Percentage (%)",
    title = paste("NASC distribution - 2018, 2021, 2023 - Transect -", dp, "-", freq, "kHz")
  )
