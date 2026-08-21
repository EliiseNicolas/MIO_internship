# Description

# 01_NASC_filter_sv.R
# From obsaustral ncdf transect data 2018, 2021, 2023 : 
#   1) filter by channel, keep only 200kHz data
#   2) cut data on depth dimension to remove most NA
#   3) filter lon >40°E and lat <-30°N
#   4) filter day/night 
# => save intermediary rds file : 2018_day.rds, 2018_night.rds, 2021_day.rds,...

# 02_NASC_concat_sv_diurnal_period.R
#   5) Concatenate each year of same period (i.e. 2018_day, 2021_day, 2023_day)
# => save intermediary rds file : 2018_2021_2022_2023_day.rds, 2018_2021_2022_2023_night.rds

# 03_NASC_check_pigmeann_grid.R
#   6) get pigmeann grid (check if pigmeann grid consistent on every year of data)

# 04_NASC_mean_profile_pig_grid.R
#   7) mean profile by grid (for day and night data)
# => save intermediary rds file : mean_pig_grid_2018_2021_2023_day.rds, ...

# 05_NASC_mean_pig_grid
#   8) compute NASC for each mean profile 
# => save intermediary rds file : nasc_mean_pig_grid_2018_2021_2023_day.rds,..

# rm(list=ls())

# packages
library(ggplot2)
library(patchwork)

# Paths
path_pig_grid <- "F:/data_elise/NASC/pigmeann_grid.rds"
diurnal_period <- "night"
path_sv <-paste0("F:/data_elise/sv_cropped/120kHz/Sv_2018_2021_2023_", diurnal_period, "_120kHz.rds")




pig_grid <- readRDS(path_pig_grid)
sv <- readRDS(path_sv)

################################################################### plots des transects sur la grille pigmeann
# construire un data.frame propre à partir des coordonnées réelles
df_sv <- data.frame(
  lon = sv$lon,
  lat = sv$lat, 
  time = sv$time
)


# retirer les NA éventuels
df_sv <- df_sv[!is.na(df_sv$lon) & !is.na(df_sv$lat), ]
df_sv$date <- as.Date(df_sv$time, origin = "1950-01-01")
df_sv$year <- factor(format(df_sv$date, "%Y"))


# plots 
ggplot() +
  
  # grille PIGMeANN
  geom_vline(
    xintercept = pig_grid$lon,
    color = "grey80",
    linewidth = 0.15
  ) +
  geom_hline(
    yintercept = pig_grid$lat,
    color = "grey80",
    linewidth = 0.15
  ) +
  
  # points Sv colorés par année
  geom_point(
    data = df_sv,
    aes(x = lon, y = lat, color = year),
    alpha = 0.3,
    size = 0.4
  ) +
  
  scale_color_manual(
    values = c("2018" = "#1b9e77", "2021" = "#d95f02", "2023" = "#7570b3")
  ) +
  
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2))) +
  
  coord_fixed() +
  
  labs(
    title = paste0("Points Sv sur la grille PIGMeANN - ", diurnal_period, " 120kHz"),
    x = "Longitude (°E)",
    y = "Latitude (°)",
    color = "Année"
  ) +
  
  theme_minimal()


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
sv_day <- format(sv_date, "%Y-%m-%d")
# Hidtogramme de nombre de profils par jour 
daily_counts <- as.data.frame(table(sv_date))

# Renommer les colonnes
colnames(daily_counts) <- c("date", "n_profils")

# Histogramme
ggplot(daily_counts, aes(x = date, y = n_profils)) +
  geom_col() +
  labs(
    x = "Date",
    y = "Nombre de profils",
    title = "Nombre de profils par jour",
    subtitle = "Données transect 2018, 2021, 2023, 120kHz, day"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5)
  )

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
  sv_day,
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

head(group_info)

# Temps moyen
mean_time_by_cell <- tapply(
  sv$time,
  cell_id_day,
  mean,
  na.rm = TRUE
)


# Objet final
mean_pig_grid_by_day <- list(
  profiles  = mean_profiles_db,
  lon       = as.numeric(group_info[, 1]),
  lat       = as.numeric(group_info[, 2]),
  date      = as.Date(group_info[, 3]),
  time      =  as.POSIXct(
    mean_time_by_cell[group_levels] * 86400,
    origin = "1950-01-01",
    tz = "UTC"
  ),
  depth     = sv$depth,
  n_profils = as.vector(table(cell_id_day))
)

str(mean_pig_grid_by_day, max.level = 1) # verif du dataset final

# --- Sauvegarde ---

# fichier global avec toutes les années (time moyen conservé)
saveRDS(
  mean_pig_grid_by_day,
  paste0("F:/data_elise/NASC/120kHz/mean_sv_profile_pig_grid_by_year_2018_2021_2023_", diurnal_period, "_120kHz_v2.rds", ".rds")
)


# ----------- Plot des effets du moyennage

# ------------------------------------------------------------
# 1. Variance par cellule × jour et par profondeur
#    sur les valeurs linéaires
# ------------------------------------------------------------

sum_sq_by_cell <- rowsum(
  sv_linear^2,
  group = cell_id_day,
  na.rm = TRUE
)

var_profiles_linear <- (
  sum_sq_by_cell -
    count_by_cell * mean_profiles_linear^2
) / (count_by_cell - 1)

# Pas de variance si moins de 2 profils
var_profiles_linear[count_by_cell <= 1] <- NA


# ------------------------------------------------------------
# 2. Data frame récapitulatif par cellule × jour
# ------------------------------------------------------------

n_profils_df <- data.frame(
  cell_id   = rownames(mean_profiles_db),
  n_profils = mean_pig_grid_by_day$n_profils,
  date      = mean_pig_grid_by_day$date
)

# Vérification
head(n_profils_df)


# ------------------------------------------------------------
# 3. Statistiques sur le nombre de profils
# ------------------------------------------------------------

mean_n <- mean(
  n_profils_df$n_profils,
  na.rm = TRUE
)

sd_n <- sd(
  n_profils_df$n_profils,
  na.rm = TRUE
)


# ------------------------------------------------------------
# 4. Statistiques sur la variance
# ------------------------------------------------------------

mean_variance <- mean(
  var_profiles_linear,
  na.rm = TRUE
)

sd_variance <- sd(
  var_profiles_linear,
  na.rm = TRUE
)

# PLOTS

p1 <- ggplot(
  n_profils_df,
  aes(x = n_profils)
) +
  geom_histogram(
    bins = 30,
    fill = "steelblue",
    color = "white"
  ) +
  geom_vline(
    xintercept = mean_n,
    color = "red",
    linetype = "dashed",
    linewidth = 0.6
  ) +
  labs(
    title = "Nombre de profils par cellule et par jour",
    subtitle = paste0(
      "Moyenne = ", round(mean_n, 1),
      " ± ", round(sd_n, 1)
    ),
    x = "Nombre de profils",
    y = "Nombre de cellules × jours"
  ) +
  theme_minimal()

# ------------------------------------------------------------
# Data frame pour le plot de variance
# ------------------------------------------------------------

# ------------------------------------------------------------
# Data frame pour le plot de variance
# ------------------------------------------------------------

variance_values <- as.vector(var_profiles_linear)

var_df <- data.frame(
  variance = variance_values
)

# Vérification
class(var_df)
str(var_df)
dim(var_df)

# Retirer les NA et les valeurs <= 0
var_df <- var_df[
  is.finite(var_df$variance) &
    var_df$variance > 0,
  ,
  drop = FALSE
]

p2 <- ggplot(
  var_df,
  aes(x = variance)
) +
  geom_histogram(
    bins = 40,
    fill = "darkorange",
    color = "white"
  ) +
  geom_vline(
    xintercept = mean_variance,
    color = "red",
    linetype = "dashed",
    linewidth = 0.6
  ) +
  scale_x_log10() +
  labs(
    title = "Variance du moyennage",
    subtitle = paste0(
      "Moyenne = ", signif(mean_variance, 3)
    ),
    x = "Variance (échelle linéaire, log10)",
    y = "Fréquence"
  ) +
  theme_minimal()

p1 / p2


# # --- 1. Variance par cellule (× année) et par profondeur, sur les valeurs linéaires ---
# sum_sq_by_cell <- rowsum(sv_linear^2, group = cell_id_day, na.rm = TRUE)
# 
# var_profiles_linear <- (sum_sq_by_cell - count_by_cell * mean_profiles_linear^2) /
#   (count_by_cell - 1)
# var_profiles_linear[count_by_cell <= 1] <- NA
# 
# # --- 2. Data frame récapitulatif par cellule × année ---
# n_profils_df <- data.frame(
#   cell_id   = rownames(mean_profiles_db),
#   n_profils = mean_pig_grid_by_year$n_profils,
#   year      = mean_years   # dérivé du time moyen, calculé précédemment
# )
# 
# mean_n <- mean(n_profils_df$n_profils, na.rm = TRUE)
# sd_n   <- sd(n_profils_df$n_profils, na.rm = TRUE)
# 
# mean_variance <- mean(var_profiles_linear, na.rm = TRUE)
# sd_variance   <- sd(var_profiles_linear, na.rm = TRUE)
# 
# # --- 3. Plots ---
# library(ggplot2)
# library(patchwork)
# 
# # Plot A : distribution du nombre de profils par cellule, facetté par année
# p1 <- ggplot(n_profils_df, aes(x = n_profils)) +
#   geom_histogram(bins = 30, fill = "steelblue", color = "white") +
#   geom_vline(xintercept = mean_n, color = "red", linetype = "dashed", linewidth = 0.6) +
#   annotate(
#     "text",
#     x = mean_n, y = Inf,
#     label = paste0("moyenne = ", round(mean_n, 1), " ± ", round(sd_n, 1)),
#     vjust = 2, hjust = -0.05, color = "red", size = 3.5
#   ) +
#   facet_wrap(~ year) +
#   labs(
#     title = "Nombre de profils par point de grille (par année)",
#     x = "Nombre de profils",
#     y = "Nombre de cellules"
#   ) +
#   theme_minimal()
# 
# # Plot B : distribution de la variance par cellule × année × profondeur
# var_df <- data.frame(
#   variance = as.vector(var_profiles_linear),
#   year     = rep(mean_years, times = ncol(var_profiles_linear))
# )
# var_df <- var_df[!is.na(var_df$variance), ]
# 
# p2 <- ggplot(var_df, aes(x = variance)) +
#   geom_histogram(bins = 40, fill = "darkorange", color = "white") +
#   geom_vline(xintercept = mean_variance, color = "red", linetype = "dashed", linewidth = 0.6) +
#   annotate(
#     "text",
#     x = mean_variance, y = Inf,
#     label = paste0("moyenne = ", signif(mean_variance, 3)),
#     vjust = 2, hjust = -0.05, color = "red", size = 3.5
#   ) +
#   scale_x_log10() +
#   facet_wrap(~ year) +
#   labs(
#     title = "Variance du moyennage par année (échelle linéaire, cellule × profondeur)",
#     x = "Variance (échelle log10)",
#     y = "Fréquence"
#   ) +
#   theme_minimal()
# 
# p1 / p2   # empilés verticalement (plus lisible avec le facettage)

 