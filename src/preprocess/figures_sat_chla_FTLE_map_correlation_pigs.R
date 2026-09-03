library(dplyr)
library(tidyr)
library(ggplot2)
library(factoextra)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(lubridate)
library(patchwork)

# ------------------------------------------------------------
# Options
# ------------------------------------------------------------
path_all_ds <- "F:/data_elise/prediction_ds/ds_ftle_pig_fod_ALL_DATES.rds"
out_dir     <- "C:/Users/mmolinet/elisou_ta_stagiaire_pref/MIO_internship_III/figures/figures_add_report"

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

all_ds <- readRDS(path_all_ds)
str(all_ds)

date <- "2023-01-26"

# Stride used to subsample the (lon, lat) grid when building the pigment
# correlation matrices (all 133 dates x 1080 x 720 pixels is ~100M rows per
# variable -> far too much for cor() in memory). Increase to go faster /
# use less RAM, decrease (min 1) to use more of the data.
subsample_stride <- 1

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
world <- ne_countries(scale = "medium", returnclass = "sf")

lon <- all_ds$lon
lat <- all_ds$lat

# Build a long data.frame (lon, lat, value) from a [lon, lat] matrix,
# matching the column-major flattening order (lon varies fastest).
mat_to_df <- function(mat2d, lon, lat) {
  data.frame(
    lon   = rep(lon, times = length(lat)),
    lat   = rep(lat, each  = length(lon)),
    value = as.vector(mat2d)
  )
}

plot_map <- function(df, title, fill_label, palette = "viridis") {
  ggplot() +
    geom_raster(data = df, aes(x = lon, y = lat, fill = value)) +
    geom_sf(data = world, fill = "grey85", color = "grey40", linewidth = 0.2) +
    coord_sf(xlim = range(lon, na.rm = TRUE),
             ylim = range(lat, na.rm = TRUE),
             expand = FALSE) +
    scale_fill_viridis_c(option = palette, na.value = NA, name = fill_label) +
    labs(title = title, x = "Longitude", y = "Latitude") +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank())
}

# ------------------------------------------------------------
# 1) Carte de Chla a la date
# ------------------------------------------------------------
date_idx <- which(all_ds$date == as.Date(date))
if (length(date_idx) == 0) stop("Date not found in all_ds$date: ", date)

chla_mat <- all_ds$pig$Chla[date_idx, , ]           # [lon, lat]
df_chla  <- mat_to_df(chla_mat, lon, lat)

p_chla <- plot_map(
  df_chla,
  title      = paste("Carte de concentration de Chlorophylle A -", date),
  fill_label = expression(Chl*italic(a)~(mg~m^-3)),
  palette    = "viridis"
)
print(p_chla)
ggsave(file.path(out_dir, paste0("map_chla_", date, ".png")),
       p_chla, width = 8, height = 6, dpi = 300)

# ------------------------------------------------------------
# 2) Carte de FTLE a la date
# ------------------------------------------------------------
ftle_mat <- all_ds$ftle[date_idx, , ]               # [lon, lat]
df_ftle  <- mat_to_df(ftle_mat, lon, lat)

p_ftle <- plot_map(
  df_ftle,
  title      = paste("Carte de FTLE -", date),
  fill_label = expression(FTLE~(d^-1)),
  palette    = "magma"
)
print(p_ftle)
ggsave(file.path(out_dir, paste0("map_ftle_", date, ".png")),
       p_ftle, width = 8, height = 6, dpi = 300)

# ------------------------------------------------------------
# Flatten pigments (subsampled grid) for correlation matrices
# ------------------------------------------------------------
pig_names <- names(all_ds$pig)
# garder seulement les concentrations brutes
pig_names <- pig_names[!grepl("_chla$", pig_names)]
pig_names <- pig_names[!grepl("_total$", pig_names)]
print((pig_names))

lon_idx <- seq(1, length(lon), by = subsample_stride)
lat_idx <- seq(1, length(lat), by = subsample_stride)

# Matrix with one column per pigment variable, flattened in the same
# (date, lon, lat) order for every variable.
pig_mat <- sapply(pig_names, function(nm) {
  as.vector(all_ds$pig[[nm]][, lon_idx, lat_idx])
})
colnames(pig_mat) <- pig_names

n_date  <- length(all_ds$date)
n_lon_s <- length(lon_idx)
n_lat_s <- length(lat_idx)

year_vec <- rep(year(all_ds$date), times = n_lon_s * n_lat_s)

# ------------------------------------------------------------
# 3) Matrice de correlation des pigments, toutes dates confondues
# ------------------------------------------------------------
cor_all <- cor(pig_mat, use = "pairwise.complete.obs")

df_cor_all <- as.data.frame(as.table(cor_all))
names(df_cor_all) <- c("var1", "var2", "corr")

p_cor_all <- ggplot(df_cor_all, aes(var1, var2, fill = corr)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", corr)), size = 2.5) +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                       midpoint = 0, limits = c(-1, 1), name = "r") +
  labs(title = "Correlation des pigments - toutes dates confondues",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_cor_all)
ggsave(file.path(out_dir, "corr_pigments_all_dates.png"),
       p_cor_all, width = 9, height = 8, dpi = 300)

# ------------------------------------------------------------
# 4) Matrice de correlation des pigments, par annee
# ------------------------------------------------------------
years <- sort(unique(year_vec))

# ------------------------------------------------------------
# Repartir de zéro pour éviter les objets pollués des essais précédents
# ------------------------------------------------------------
rm(list = c("df_cor_year", "plots_list", "p_cor_year")[
  c("df_cor_year", "plots_list", "p_cor_year") %in% ls()
])

# ------------------------------------------------------------
# 1) Matrices de corrélation par année
# ------------------------------------------------------------
df_cor_year <- do.call(rbind, lapply(years, function(y) {
  cor_y <- cor(pig_mat[year_vec == y, , drop = FALSE],
               use = "pairwise.complete.obs")
  d <- as.data.frame(as.table(cor_y))
  names(d) <- c("var1", "var2", "corr")
  d$year <- as.character(y)
  d
}))

# ------------------------------------------------------------
# 2) Matrice toutes années confondues -> label combiné
# ------------------------------------------------------------
all_label <- paste(years, collapse = "-")   # "2018-2021-2022-2023"

cor_all_y <- cor(pig_mat, use = "pairwise.complete.obs")
df_cor_all_y <- as.data.frame(as.table(cor_all_y))
names(df_cor_all_y) <- c("var1", "var2", "corr")
df_cor_all_y$year <- all_label

df_cor_year <- rbind(df_cor_year, df_cor_all_y)

df_cor_year$year <- factor(df_cor_year$year,
                           levels = c(all_label, as.character(years)))

# ------------------------------------------------------------
# 3) Fonction pour construire un seul panel
# ------------------------------------------------------------
library(patchwork)

make_corr_plot <- function(df, title) {
  ggplot(df, aes(var1, var2, fill = corr)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", corr)), size = 2) +
    scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                         midpoint = 0, limits = c(-1, 1), name = "r") +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = 9) +
    theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 6),
          axis.text.y  = element_text(size = 6),
          plot.title   = element_text(size = 10, hjust = 0.5))
}

# ------------------------------------------------------------
# 4) Un ggplot par panel (reconstruit maintenant que les labels sont bons)
# ------------------------------------------------------------
panel_levels <- levels(df_cor_year$year)
years_chr    <- setdiff(panel_levels, all_label)

plots_list <- lapply(panel_levels, function(lv) {
  make_corr_plot(df_cor_year[df_cor_year$year == lv, ], title = lv)
})
names(plots_list) <- panel_levels

stopifnot(length(years_chr) == 4)

# ------------------------------------------------------------
# 5) Assemblage patchwork
# ------------------------------------------------------------
design <- "
#AA#
BBCC
DDEE
"

p_cor_year <- wrap_plots(
  plots_list[[all_label]],
  plots_list[[years_chr[1]]],
  plots_list[[years_chr[2]]],
  plots_list[[years_chr[3]]],
  plots_list[[years_chr[4]]],
  design = design,
  guides = "collect"
) +
  plot_annotation(title = "Correlation des concentrations pigments")

print(p_cor_year)

ggsave(file.path(out_dir, "corr_pigments_by_year.png"),
       p_cor_year, width = 12, height = 10, dpi = 300)


message("4 figures saved to: ", out_dir)


############## correlation nasc pigments

freq <- 120
nasc_ds <- readRDS(paste0(
  "F:/data_elise/ds_NASC_pig_ftle_fod/ds_NASC_per_esu_all/ds_NASC_per_esu_pig_ftle_fod_2018_2021_2022_2023_transect_",
  freq, "kHz_mask9.rds"
))
str(nasc_ds)

# ------------------------------------------------------------
# Sélection des variables se terminant par "_Chla"
# ------------------------------------------------------------
pig_chla_vars <- grep("_Chla$", names(nasc_ds), value = TRUE)
print(pig_chla_vars)

# "Chla_Chla" est un ratio de Chla sur lui-même (toujours = 1) : à exclure,
# sinon la facette correspondante n'aura aucun intérêt
pig_chla_vars <- setdiff(pig_chla_vars, "Chla_Chla")

# ------------------------------------------------------------
# Passage en format long
# ------------------------------------------------------------
df_long <- nasc_ds %>%
  select(nasc, all_of(pig_chla_vars)) %>%
  pivot_longer(cols = all_of(pig_chla_vars),
               names_to = "pigment", values_to = "value") %>%
  filter(is.finite(nasc), is.finite(value), nasc > 0)   # nasc > 0 pour le log

# ------------------------------------------------------------
# Coefficient de corrélation (Pearson) par pigment
# ------------------------------------------------------------
cor_df <- df_long %>%
  group_by(pigment) %>%
  summarise(r = cor(nasc, value, use = "pairwise.complete.obs"),
            n = n(),
            .groups = "drop") %>%
  mutate(label = sprintf("r = %.2f (n = %s)", r, format(n, big.mark = " ")))

# positions des labels : coin haut-gauche de chaque facette
label_pos <- df_long %>%
  group_by(pigment) %>%
  summarise(x = min(value, na.rm = TRUE),
            y = max(nasc, na.rm = TRUE),
            .groups = "drop") %>%
  left_join(cor_df, by = "pigment")

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------
p_nasc_pig <- ggplot(df_long, aes(x = value, y = nasc)) +
  geom_point(alpha = 0.1, size = 0.5, color = "steelblue") +
  geom_smooth(method = "lm", color = "firebrick", se = TRUE, linewidth = 0.7) +
  geom_text(data = label_pos, aes(x = x, y = y, label = label),
            hjust = 0, vjust = 1, size = 3, inherit.aes = FALSE) +
  facet_wrap(~ pigment, scales = "free_x", ncol = 3) +
  scale_y_log10() +
  labs(title = "Covariation du NASC avec les ratios pigmentaires (2018, 2021, 2022, 2023)",
       x = "Ratio pigmentaire", y = "NASC (echelle log10)") +
  theme_minimal(base_size = 10) +
  theme(strip.text = element_text(face = "bold"))

print(p_nasc_pig)

ggsave(file.path(out_dir, "nasc_vs_pigments_chla.png"),
       p_nasc_pig, width = 12, height = 10, dpi = 300)

# ==============================================================
# 1) NASC vs FTLE (variable continue)
# ==============================================================
df_ftle <- nasc_ds %>%
  select(nasc, ftle) %>%
  filter(is.finite(nasc), is.finite(ftle), nasc > 0)

r_ftle <- cor(df_ftle$nasc, df_ftle$ftle, use = "pairwise.complete.obs")
label_ftle <- sprintf("r = %.2f (n = %s)", r_ftle,
                      format(nrow(df_ftle), big.mark = " "))

p_nasc_ftle <- ggplot(df_ftle, aes(x = ftle, y = nasc)) +
  geom_hex(bins = 40) +
  scale_fill_viridis_c(trans = "log10", name = "count") +
  geom_smooth(method = "lm", color = "firebrick", se = FALSE, linewidth = 0.7) +
  annotate("text", x = min(df_ftle$ftle, na.rm = TRUE),
           y = max(df_ftle$nasc, na.rm = TRUE),
           label = label_ftle, hjust = 0, vjust = 1, size = 3.5, color = "black") +
  scale_y_log10() +
  labs(title = "Covariation du NASC avec le FTLE(2018, 2021, 2022, 2023)",
       x = "FTLE", y = "NASC (echelle log10)") +
  theme_minimal(base_size = 11)

print(p_nasc_ftle)

ggsave(file.path(out_dir, "nasc_vs_ftle.png"),
       p_nasc_ftle, width = 8, height = 6, dpi = 300)

# ==============================================================
# 2) NASC vs FOD (variable categorielle)
# ==============================================================
# ==============================================================
# NASC vs FOD, avec couleurs/labels alignes sur les clusters
# ==============================================================
cluster_cols <- c(
  "1" = "#2166AC",
  "2" = "#67A9CF",
  "3" = "#1A9850",
  "4" = "#A6D96A",
  "5" = "#FDAE61",
  "6" = "#D73027"
)

transition_cols <- c(
  "7"  = "#3B73B9",  # 1-2
  "8"  = "#3FA7B5",  # 1-3
  "9"  = "#45B97C",  # 2-3
  "10" = "#C46A00",  # 3-4
  "11" = "#E85D04",  # 4-6
  "12" = "#C9184A",  # 4-5
  "13" = "#8F1D3F"   # 5-6
)
fod_cols <- c(cluster_cols, transition_cols)

# ------------------------------------------------------------
# Renommage des codes FOD en libelles lisibles (clusters = C.,
# transitions = T.-.), regroupes par cluster de depart (toutes les
# transitions issues de C1, puis C2, etc.)
# ------------------------------------------------------------

legend_codes <- c(
  1, 7, 8, 2, 9, 3, 10, 4, 12, 11, 5, 13, 6
)

legend_labels <- c(
  "C1", "T1-2", "T1-3",
  "C2", "T2-3",
  "C3", "T3-4",
  "C4", "T4-5", "T4-6",
  "C5", "T5-6", "C6"
)

fod_label_map <- setNames(legend_labels, as.character(legend_codes))
fod_cols_named <- setNames(fod_cols[as.character(legend_codes)], legend_labels)
relabel_fod <- function(fod_vec) {
  fod_chr <- trimws(as.character(fod_vec))
  
  # NA réel ou chaîne "NA" -> niveau explicite "NA"
  fod_chr[is.na(fod_chr) | fod_chr == "NA"] <- "NA"
  
  lbl <- fod_label_map[fod_chr]
  
  # Si un code n'est pas dans le mapping, on conserve le code
  unknown <- is.na(lbl)
  lbl[unknown] <- fod_chr[unknown]
  
  factor(
    lbl,
    levels = c(
      legend_labels,
      sort(setdiff(unique(lbl), legend_labels))
    )
  )
}
df_fod <- nasc_ds %>%
  select(nasc, fod) %>%
  filter(is.finite(nasc), nasc > 0) %>%
  mutate(fod_label = relabel_fod(fod)) %>%
  filter(fod_label != "NA")   # on retire la categorie NA du plot

# effectifs par categorie, pour les afficher sous les boites
n_by_fod <- df_fod %>%
  count(fod_label) %>%
  mutate(label = paste0("n=", format(n, big.mark = " ")))

p_nasc_fod <- ggplot(df_fod, aes(x = fod_label, y = nasc, fill = fod_label)) +
  geom_boxplot(outlier.alpha = 0.1, outlier.size = 0.5) +
  geom_text(data = n_by_fod,
            aes(x = fod_label, y = min(df_fod$nasc, na.rm = TRUE), label = label),
            inherit.aes = FALSE, vjust = 1.5, size = 3) +
  scale_y_log10() +
  scale_fill_manual(values = fod_cols_named, guide = "none",
                    na.value = "grey70") +   # au cas ou un code inconnu resterait
  labs(title = "Covariation du NASC avec le FOD",
       x = "FOD", y = "log10(NASC)") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_nasc_fod)

ggsave(file.path(out_dir, "nasc_vs_fod.png"),
       p_nasc_fod, width = 9, height = 6, dpi = 300)

