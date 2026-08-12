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

# packages
library(ggplot2)

# Paths
path_pig_grid <- "~/Documents/stage_MIO/pt_III/data_preprocessed/NASC/transect_2018_2022_2023/pigmeann_grid.rds"
diurnal_period <- "day"
path_sv <- paste0("/mnt/KER22/MIO_internship_III/data_preprocessed/NASC/transect_2018_2022_2023/Sv_2018_2021_2023_", diurnal_period,"_200kHz.rds")




pig_grid <- readRDS(path_pig_grid)
sv <- readRDS(path_sv)

# j'avais oublié le filtrage des profondeurs de surface pour sv sur mes scipts précédents
idx_depth <- which(sv$depth > 15)
sv$profiles <- sv$profiles[, idx_depth]
sv$depth <- sv$depth[idx_depth]

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

# plot 
dev.off()      # ferme le device courant
# si erreur "cannot shut down device 1", répéter jusqu'à ce que ça échoue proprement
graphics.off() # ferme tous les devices
# --- Plot coloré par année ---
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
    title = paste0("Points Sv sur la grille PIGMeANN - ", diurnal_period),
    x = "Longitude (°E)",
    y = "Latitude (°)",
    color = "Année"
  ) +
  
  theme_minimal()

# surface d'une grille
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
sv_year <- format(sv_date, "%Y")

table(sv_year)

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
cell_id_year <- paste(lon_cell, lat_cell, sv_year, sep = "_")


# --- 2. Conversion dB -> linéaire ---
sv_linear <- 10^(sv$profiles / 10)

# --- 3. Moyenne par cellule ET par année, robuste aux NA ---
sum_by_cell   <- rowsum(sv_linear, group = cell_id_year, na.rm = TRUE)
not_na_mat    <- !is.na(sv_linear)
count_by_cell <- rowsum(not_na_mat * 1, group = cell_id_year)

mean_profiles_linear <- sum_by_cell / count_by_cell
mean_profiles_linear[count_by_cell == 0] <- NA

# --- 4. Reconversion en dB ---
mean_profiles_db <- 10 * log10(mean_profiles_linear)

# --- 5. Reconstruire l'objet avec lon / lat / time moyen ---
cell_names   <- rownames(mean_profiles_db)
cell_coords  <- do.call(rbind, strsplit(cell_names, "_"))

# temps moyen par cellule (même groupement que pour les profils)
mean_time_by_cell <- tapply(sv$time, cell_id_year, mean, na.rm = TRUE)

mean_pig_grid_by_year <- list(
  profiles  = mean_profiles_db,               # matrice n_cellules(x année) x n_depth, en dB
  lon       = as.numeric(cell_coords[, 1]),
  lat       = as.numeric(cell_coords[, 2]),
  time      = as.numeric(mean_time_by_cell[cell_names]),  # temps moyen (même unité que sv$time)
  depth     = sv$depth,
  n_profils = as.vector(table(cell_id_year)[cell_names])
)

str(mean_pig_grid_by_year, max.level = 1)

# vérification : reconvertir en date pour contrôle visuel
mean_dates <- as.Date(mean_pig_grid_by_year$time, origin = "1950-01-01")  # adapte l'origine si besoin
table(format(mean_dates, "%Y"))  # vérifier la répartition par année déduite du temps moyen

# --- Sauvegarde ---

# fichier global avec toutes les années (time moyen conservé)
saveRDS(
  mean_pig_grid_by_year,
  paste0("/home/elise/Documents/stage_MIO/pt_III/MIO_internship_III/data_preprocessed/NASC/transect_2018_2022_2023/mean_pig_grid_by_year_2018_2021_2023_", diurnal_period, ".rds")
)

mean_dates <- as.Date(mean_pig_grid_by_year$time, origin = "1950-01-01")  # adapte l'origine si besoin
mean_years <- format(mean_dates, "%Y")

# ----------- Plot des effets du moyennage
# --- 1. Variance par cellule (× année) et par profondeur, sur les valeurs linéaires ---
sum_sq_by_cell <- rowsum(sv_linear^2, group = cell_id_year, na.rm = TRUE)

var_profiles_linear <- (sum_sq_by_cell - count_by_cell * mean_profiles_linear^2) /
  (count_by_cell - 1)
var_profiles_linear[count_by_cell <= 1] <- NA

# --- 2. Data frame récapitulatif par cellule × année ---
n_profils_df <- data.frame(
  cell_id   = rownames(mean_profiles_db),
  n_profils = mean_pig_grid_by_year$n_profils,
  year      = mean_years   # dérivé du time moyen, calculé précédemment
)

mean_n <- mean(n_profils_df$n_profils, na.rm = TRUE)
sd_n   <- sd(n_profils_df$n_profils, na.rm = TRUE)

mean_variance <- mean(var_profiles_linear, na.rm = TRUE)
sd_variance   <- sd(var_profiles_linear, na.rm = TRUE)

# --- 3. Plots ---
library(ggplot2)
library(patchwork)

# Plot A : distribution du nombre de profils par cellule, facetté par année
p1 <- ggplot(n_profils_df, aes(x = n_profils)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white") +
  geom_vline(xintercept = mean_n, color = "red", linetype = "dashed", linewidth = 0.6) +
  annotate(
    "text",
    x = mean_n, y = Inf,
    label = paste0("moyenne = ", round(mean_n, 1), " ± ", round(sd_n, 1)),
    vjust = 2, hjust = -0.05, color = "red", size = 3.5
  ) +
  facet_wrap(~ year) +
  labs(
    title = "Nombre de profils par point de grille (par année)",
    x = "Nombre de profils",
    y = "Nombre de cellules"
  ) +
  theme_minimal()

# Plot B : distribution de la variance par cellule × année × profondeur
var_df <- data.frame(
  variance = as.vector(var_profiles_linear),
  year     = rep(mean_years, times = ncol(var_profiles_linear))
)
var_df <- var_df[!is.na(var_df$variance), ]

p2 <- ggplot(var_df, aes(x = variance)) +
  geom_histogram(bins = 40, fill = "darkorange", color = "white") +
  geom_vline(xintercept = mean_variance, color = "red", linetype = "dashed", linewidth = 0.6) +
  annotate(
    "text",
    x = mean_variance, y = Inf,
    label = paste0("moyenne = ", signif(mean_variance, 3)),
    vjust = 2, hjust = -0.05, color = "red", size = 3.5
  ) +
  scale_x_log10() +
  facet_wrap(~ year) +
  labs(
    title = "Variance du moyennage par année (échelle linéaire, cellule × profondeur)",
    x = "Variance (échelle log10)",
    y = "Fréquence"
  ) +
  theme_minimal()

p1 / p2   # empilés verticalement (plus lisible avec le facettage)
# 
# assign_cell <- function(coord, grid_vals) {
#   grid_sorted <- sort(unique(grid_vals))  # gère l'ordre décroissant automatiquement
#   
#   idx <- findInterval(coord, grid_sorted, all.inside = TRUE)
#   d_before <- abs(coord - grid_sorted[idx])
#   d_after  <- abs(coord - grid_sorted[pmin(idx + 1, length(grid_sorted))])
#   idx_final <- ifelse(d_after < d_before, idx + 1, idx)
#   
#   grid_sorted[idx_final]
# }
# 
# lon_cell <- assign_cell(sv$lon, pig_grid$lon)
# lat_cell <- assign_cell(sv$lat, pig_grid$lat)
# 
# cell_id <- paste(lon_cell, lat_cell, sep = "_")
# 
# # --- 2. Conversion dB -> linéaire avant moyennage ---
# sv_linear <- 10^(sv$profiles / 10)   # matrice n_profils x n_depth
# 
# # --- 3. Moyenne par cellule de grille, pour chaque profondeur ---
# library(dplyr)
# 
# # calcul de la moyenne linéaire par cellule (colonne = profondeur)
# # somme des valeurs non-NA par cellule et par profondeur
# sum_by_cell <- rowsum(sv_linear, group = cell_id, na.rm = TRUE)
# 
# # comptage des valeurs non-NA par cellule et par profondeur
# not_na_mat <- !is.na(sv_linear)
# count_by_cell <- rowsum(not_na_mat * 1, group = cell_id)
# 
# # moyenne = somme / comptage (élément par élément)
# mean_profiles_linear <- sum_by_cell / count_by_cell
# mean_profiles_linear[count_by_cell == 0] <- NA
# 
# # --- 4. Reconversion en dB ---
# mean_profiles_db <- 10 * log10(mean_profiles_linear)
# 
# # --- 5. Reconstruire un objet propre avec coordonnées de grille ---
# 
# # extraire lon/lat depuis les noms de ligne (rownames = cell_id = "lon_lat")
# cell_names <- rownames(mean_profiles_db)
# cell_coords <- do.call(rbind, strsplit(cell_names, "_"))
# 
# mean_pig_grid <- list(
#   profiles  = mean_profiles_db,                    # matrice n_cellules x n_depth (en dB)
#   lon       = as.numeric(cell_coords[, 1]),
#   lat       = as.numeric(cell_coords[, 2]),
#   depth     = sv$depth,
#   n_profils = as.vector(table(cell_id)[cell_names])
# )
# 
# str(mean_pig_grid, max.level = 1)
# 
# graphics.off()
# 
# ggplot() +
#   
#   # grille PIGMeANN
#   geom_vline(
#     xintercept = pig_grid$lon,
#     color = "grey85",
#     linewidth = 0.15
#   ) +
#   geom_hline(
#     yintercept = pig_grid$lat,
#     color = "grey85",
#     linewidth = 0.15
#   ) +
#   
#   # points correspondant aux profils moyens par cellule
#   geom_point(
#     data = data.frame(
#       lon = mean_pig_grid$lon,
#       lat = mean_pig_grid$lat,
#       n_profils = mean_pig_grid$n_profils
#     ),
#     aes(x = lon, y = lat, color = n_profils),
#     size = 1
#   ) +
#   
#   scale_color_viridis_c(name = "N profils") +
#   
#   coord_fixed() +
#   
#   labs(
#     title = paste0("Profils moyens par cellule PIGMeANN - ", diurnal_period),
#     x = "Longitude (°E)",
#     y = "Latitude (°)"
#   ) +
#   
#   theme_minimal()
# 
# # diagnostic effet moyennage profil seon grid pigmeann
# 
# # --- 1. Variance par cellule et par profondeur (sur les valeurs linéaires) ---
# 
# # somme des carrés par cellule (nécessaire pour calculer la variance)
# sum_sq_by_cell <- rowsum(sv_linear^2, group = cell_id, na.rm = TRUE)
# 
# # variance = (somme des carrés - n * moyenne^2) / (n - 1)
# # on réutilise sum_by_cell, count_by_cell, mean_profiles_linear déjà calculés
# var_profiles_linear <- (sum_sq_by_cell - count_by_cell * mean_profiles_linear^2) /
#   (count_by_cell - 1)
# 
# # les cellules avec 1 seul profil ont une variance non définie -> NA
# var_profiles_linear[count_by_cell <= 1] <- NA
# 
# # --- 2. Résumé du nombre de profils par cellule ---
# n_profils_df <- data.frame(
#   cell_id   = names(mean_pig_grid$n_profils) %||% rownames(mean_profiles_db),
#   n_profils = mean_pig_grid$n_profils
# )
# 
# mean_n <- mean(n_profils_df$n_profils, na.rm = TRUE)
# sd_n   <- sd(n_profils_df$n_profils, na.rm = TRUE)
# 
# cat("Nombre moyen de profils par cellule :", round(mean_n, 1),
#     "± (écart-type)", round(sd_n, 1), "\n")
# 
# # --- 3. Variance moyenne du moyennage (moyenne sur toutes profondeurs & cellules) ---
# mean_variance <- mean(var_profiles_linear, na.rm = TRUE)
# sd_variance   <- sd(var_profiles_linear, na.rm = TRUE)
# 
# cat("Variance moyenne (échelle linéaire) du moyennage :",
#     signif(mean_variance, 3), "±", signif(sd_variance, 3), "\n")
# 
# # --- 4. Plots ---
# library(ggplot2)
# library(patchwork)  # pour assembler 2 plots côte à côte (install.packages("patchwork") si besoin)
# 
# # Plot A : distribution du nombre de profils par cellule
# p1 <- ggplot(n_profils_df, aes(x = n_profils)) +
#   geom_histogram(bins = 30, fill = "steelblue", color = "white") +
#   geom_vline(xintercept = mean_n, color = "red", linetype = "dashed", linewidth = 0.6) +
#   annotate(
#     "text",
#     x = mean_n, y = Inf,
#     label = paste0("moyenne = ", round(mean_n, 1), " ± ", round(sd_n, 1)),
#     vjust = 2, hjust = -0.05, color = "red", size = 3.5
#   ) +
#   labs(
#     title = "Nombre de profils par point de grille",
#     x = "Nombre de profils",
#     y = "Nombre de cellules"
#   ) +
#   theme_minimal()
# 
# # Plot B : distribution de la variance par cellule/profondeur (échelle linéaire)
# var_vec <- as.vector(var_profiles_linear)
# var_vec <- var_vec[!is.na(var_vec)]
# 
# p2 <- ggplot(data.frame(variance = var_vec), aes(x = variance)) +
#   geom_histogram(bins = 40, fill = "darkorange", color = "white") +
#   geom_vline(xintercept = mean_variance, color = "red", linetype = "dashed", linewidth = 0.6) +
#   annotate(
#     "text",
#     x = mean_variance, y = Inf,
#     label = paste0("moyenne = ", signif(mean_variance, 3)),
#     vjust = 2, hjust = -0.05, color = "red", size = 3.5
#   ) +
#   scale_x_log10() +   # la variance en linéaire a souvent une distribution très asymétrique
#   labs(
#     title = "Variance du moyennage (échelle linéaire, par cellule × profondeur)",
#     x = "Variance (échelle log10)",
#     y = "Fréquence"
#   ) +
#   theme_minimal()
# 
# p1 + p2
 