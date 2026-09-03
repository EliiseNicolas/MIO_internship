# ============================================================
# Communautes phytoplanctoniques (PCA + k-means) a partir des
# ratios pigment/(somme des pigments hors Chla) du RDS unique
# deja aligne (all_ds).
#
# ADAPTATION : les ratios pigmentaires ne sont plus calcules sur
# Chla (suffixe "_chla") mais sur la somme de tous les pigments
# SAUF Chla (suffixe "_totpig"). "chla_total" (= Chla seule) est
# conserve comme variable a part entiere. "chla_totpig" (ratio de
# Chla elle-meme sur la somme des autres) existe dans all_ds mais
# n'est pas utilise ici : il ferait doublon avec "chla_total".
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(factoextra)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

# ------------------------------------------------------------
# Options
# ------------------------------------------------------------

path_all_ds <- "F:/data_elise/prediction_ds/ds_ftle_pig_fod_ALL_DATES.rds"
n_cluster   <- 6
out_dir <- "C:/Users/mmolinet/elisou_ta_stagiaire_pref/MIO_internship_III/figures/commu_phyto_cluster_pigs"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

sample_n_fit <- 200000
set.seed(123)

# ------------------------------------------------------------
# Variables utilisees pour la PCA : Chla en concentration totale +
# les 8 ratios pigment/(somme des pigments hors Chla)
# ------------------------------------------------------------
pig_vars <- c(
  "chla_total", "per_totpig", "but_totpig", "fuco_totpig", "hex_totpig",
  "allo_totpig", "zea_totpig", "chlb_totpig", "dvchla_totpig"
)

# ============================================================
# 0. Chargement + positions valides (tous ratios pigmentaires finis)
# ============================================================

all_ds <- readRDS(path_all_ds)
str(all_ds)
dates <- all_ds$date
lons  <- all_ds$lon
lats  <- all_ds$lat

valid_mask <- Reduce(`&`, lapply(pig_vars, function(v) is.finite(all_ds$pig[[v]])))
valid_idx  <- which(valid_mask, arr.ind = TRUE)   # [n_valid, 3] : date_idx, lon_idx, lat_idx

data_pca <- data.frame(
  date = dates[valid_idx[, 1]],
  lon  = lons[valid_idx[, 2]],
  lat  = lats[valid_idx[, 3]]
)

for (v in pig_vars) {
  data_pca[[v]] <- all_ds$pig[[v]][valid_idx]
}

nrow(data_pca)

# ============================================================
# 1 - normalisation des ratios pigmentaires (proportions)
# ============================================================

data_pca <- data_pca %>%
  mutate(
    pig_sum = rowSums(across(all_of(pig_vars)))
  ) %>%
  filter(pig_sum > 0) %>%
  mutate(
    across(all_of(pig_vars), ~ .x / pig_sum)
  )

nrow(data_pca)

# ------------------------------------------------------------
# Sous-echantillon pour ajuster la PCA/k-means
# ------------------------------------------------------------

if (is.finite(sample_n_fit) && nrow(data_pca) > sample_n_fit) {
  fit_idx <- sample(nrow(data_pca), sample_n_fit)
} else {
  fit_idx <- seq_len(nrow(data_pca))
}
data_fit <- data_pca[fit_idx, ]

# ------------------------------------------------------------ 2 - PCA
pca <- prcomp(
  data_fit %>% select(all_of(pig_vars)),
  center = TRUE,
  scale. = TRUE
)

summary(pca)
fviz_eig(pca, addlabels = TRUE)
fviz_pca_var(pca, col.var = "contrib", repel = TRUE)

pc_scores_fit <- as.data.frame(pca$x)
data_fit <- bind_cols(data_fit, pc_scores_fit)

pca_scores <- data_fit %>%
  select(PC1, PC2, PC3, PC4, PC5, PC6, PC7)

# ------------------------------------------------------------
# selection du nombre de cluster pour le kmeans (elbow method)
# ------------------------------------------------------------

set.seed(123)

wss <- sapply(1:10, function(k) {
  kmeans(pca_scores, centers = k, nstart = 100)$tot.withinss
})

wss_df <- data.frame(k = 1:10, WSS = wss)

p <- ggplot(wss_df, aes(x = k, y = WSS)) +
  geom_line() +
  geom_point(size = 3) +
  labs(x = "Nombre de clusters", y = "Within-cluster sum of squares",
       title = "Methode du coude") +
  theme_minimal()
print(p)
ggsave(file.path(out_dir, "elbow_method_kmeans.png"), plot = p,
       width = 8, height = 8, units = "in", dpi = 300)

km <- kmeans(pca_scores, centers = n_cluster, nstart = 100)
data_fit$community <- factor(km$cluster)

# ============================================================
# 3 - Projection de TOUTE la grille valide dans l'espace PCA
# ============================================================

pig_mat_all   <- as.matrix(data_pca[, pig_vars])
pc_scores_all <- scale(pig_mat_all, center = pca$center, scale = pca$scale) %*% pca$rotation

assign_cluster <- function(pc_mat, centers) {
  pc_mat_sub <- pc_mat[, colnames(centers), drop = FALSE]
  d <- sapply(seq_len(nrow(centers)), function(k) {
    rowSums((pc_mat_sub - matrix(centers[k, ], nrow(pc_mat_sub), ncol(centers), byrow = TRUE))^2)
  })
  max.col(-d)
}

data_pca$community <- factor(assign_cluster(pc_scores_all, km$centers))

# ============================================================
# 4 - Distributions des pigments par communaute
# ============================================================

data_long <- data_fit %>%
  select(community, all_of(pig_vars)) %>%
  pivot_longer(cols = all_of(pig_vars), names_to = "pigment", values_to = "proportion")

community_colors <- c(
  "1" = "#1b9e77", "2" = "#d95f02", "3" = "#7570b3",
  "4" = "#e7298a", "5" = "#66a61e", "6" = "#e6ab02"
)

p <- ggplot(data_long, aes(x = community, y = proportion, color = community)) +
  stat_summary(fun = mean, geom = "point", size = 2.5) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
               geom = "errorbar", width = 0.3, linewidth = 0.6) +
  facet_wrap(~ pigment, scales = "free_y") +
  scale_color_manual(values = community_colors, name = "Communaute") +
  labs(x = "Communaute", y = "Proportion du pigment (moyenne +/- 1 ecart-type)",
       title = "Distribution des pigments au sein des communautes") +
  theme_minimal()

ggsave(file.path(out_dir, "distrib_pigments_in_communautes_kmeans_k6_PFT.png"),
       plot = p, width = 14, height = 8, units = "in", dpi = 300)

library(rstatix)
library(multcompView)

# ============================================================
# 7 - Tests statistiques : differences de pigments entre communautes
# ============================================================

kw_results <- data_long %>%
  group_by(pigment) %>%
  rstatix::kruskal_test(proportion ~ community) %>%
  rstatix::adjust_pvalue(method = "BH") %>%
  rstatix::add_significance("p.adj")
print(kw_results, n = Inf)

kw_effsize <- data_long %>%
  group_by(pigment) %>%
  rstatix::kruskal_effsize(proportion ~ community)
print(kw_effsize, n = Inf)

dunn_results <- data_long %>%
  group_by(pigment) %>%
  rstatix::dunn_test(proportion ~ community, p.adjust.method = "BH")
print(dunn_results, n = Inf)

write.csv(kw_results,   file.path(out_dir, "stats_kruskal_wallis_pigments.csv"), row.names = FALSE)
write.csv(dunn_results, file.path(out_dir, "stats_dunn_posthoc_pigments.csv"),   row.names = FALSE)

get_cld <- function(dunn_df) {
  pvals <- dunn_df$p.adj
  names(pvals) <- paste(dunn_df$group1, dunn_df$group2, sep = "-")
  cld <- multcompView::multcompLetters(pvals)$Letters
  data.frame(community = names(cld), cld = cld, row.names = NULL)
}

cld_by_pigment <- dunn_results %>%
  group_by(pigment) %>%
  group_modify(~ get_cld(.x)) %>%
  ungroup() %>%
  mutate(community = factor(community, levels = levels(data_long$community)))

label_pos <- data_long %>%
  group_by(pigment, community) %>%
  summarise(y_pos = mean(proportion) + sd(proportion), .groups = "drop")

cld_plot_data <- cld_by_pigment %>%
  left_join(label_pos, by = c("pigment", "community"))

p <- ggplot(data_long, aes(x = community, y = proportion, color = community)) +
  stat_summary(fun = mean, geom = "point", size = 2.5) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
               geom = "errorbar", width = 0.3, linewidth = 0.6) +
  geom_text(data = cld_plot_data, aes(x = community, y = y_pos, label = cld),
            color = "black", vjust = -0.6, size = 3.5, fontface = "bold", inherit.aes = FALSE) +
  facet_wrap(~ pigment, scales = "free_y") +
  scale_color_manual(values = community_colors, name = "Communaute") +
  labs(x = "Communaute", y = "Proportion du pigment (moyenne +/- 1 ecart-type)",
       title = "Distribution des pigments au sein des communautes",
       subtitle = "Lettres differentes = difference significative (test de Dunn, p < 0.05, ajuste BH)") +
  theme_minimal()

print(p)

ggsave(file.path(out_dir, "distrib_pigments_in_communautes_kmeans_k6_PFT_stats.png"),
       plot = p, width = 14, height = 8, units = "in", dpi = 300)

# ------------------------------------------------------------
# Coloration des distributions par Phytoplankton Functional Type
# (PFT) -- noms de pigments mis a jour avec le suffixe "_totpig"
# ------------------------------------------------------------

data_long <- data_fit %>%
  select(community, all_of(pig_vars)) %>%
  pivot_longer(cols = all_of(pig_vars), names_to = "pigment", values_to = "proportion") %>%
  mutate(
    PFT = case_when(
      pigment %in% c("dvchla_totpig", "zea_totpig")            ~ "Picocyanobacteria",
      pigment %in% c("allo_totpig", "hex_totpig", "but_totpig") ~ "Flagellates",
      pigment == "fuco_totpig"                                  ~ "Diatoms",
      TRUE                                                      ~ "Other"
    )
  )

p_pigments <- ggplot(data_long, aes(x = pigment, y = proportion, color = PFT)) +
  stat_summary(fun = mean, geom = "point", size = 2.5, position = position_dodge(width = 0.3)) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
               geom = "errorbar", width = 0.3, linewidth = 0.6, position = position_dodge(width = 0.3)) +
  facet_wrap(~ community, ncol = 1, axes = "all_x") +
  scale_color_manual(values = c(
    "Picocyanobacteria" = "#E69F00", "Flagellates" = "#56B4E9",
    "Diatoms" = "#009E73", "Other" = "grey40"
  )) +
  labs(x = "Pigment", y = "Proportion (moyenne +/- 1 ecart-type)", color = "Functional group",
       title = "Pigment distributions within phytoplankton communities") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_pigments

ggsave(file.path(out_dir, "distribution_pigments_communautes_kmeans_k6_PFT.png"),
       plot = p_pigments, width = 8, height = 12, units = "in", dpi = 300)

# ============================================================
# 8 - Distribution des pigments par communaute + tests entre
#     pigments au sein de chaque communaute
# ============================================================

data_long <- data_fit %>%
  select(community, all_of(pig_vars)) %>%
  pivot_longer(cols = all_of(pig_vars), names_to = "pigment", values_to = "proportion") %>%
  mutate(
    PFT = case_when(
      pigment %in% c("dvchla_totpig", "zea_totpig")            ~ "Picocyanobacteria",
      pigment %in% c("allo_totpig", "hex_totpig", "but_totpig") ~ "Flagellates",
      pigment == "fuco_totpig"                                  ~ "Diatoms",
      TRUE                                                      ~ "Other"
    )
  )

kw_pigment_results <- data_long %>%
  group_by(community) %>%
  rstatix::kruskal_test(proportion ~ pigment) %>%
  rstatix::adjust_pvalue(method = "BH") %>%
  rstatix::add_significance("p.adj")
print(kw_pigment_results, n = Inf)

kw_pigment_effsize <- data_long %>%
  group_by(community) %>%
  rstatix::kruskal_effsize(proportion ~ pigment)
print(kw_pigment_effsize, n = Inf)

dunn_pigment_results <- data_long %>%
  group_by(community) %>%
  rstatix::dunn_test(proportion ~ pigment, p.adjust.method = "BH")
print(dunn_pigment_results, n = Inf)

write.csv(kw_pigment_results,
          file.path(out_dir, "stats_kruskal_wallis_pigments_within_community.csv"), row.names = FALSE)
write.csv(dunn_pigment_results,
          file.path(out_dir, "stats_dunn_pigments_within_community.csv"), row.names = FALSE)

get_cld_pigments <- function(dunn_df) {
  pvals <- dunn_df$p.adj
  names(pvals) <- paste(dunn_df$group1, dunn_df$group2, sep = "-")
  cld <- multcompView::multcompLetters(pvals)$Letters
  data.frame(pigment = names(cld), cld = cld, row.names = NULL)
}

cld_pigment_data <- dunn_pigment_results %>%
  group_by(community) %>%
  group_modify(~ get_cld_pigments(.x)) %>%
  ungroup()

label_pos_pigment <- data_long %>%
  group_by(community, pigment) %>%
  summarise(mean_prop = mean(proportion, na.rm = TRUE),
            sd_prop = sd(proportion, na.rm = TRUE),
            y_pos = mean_prop + sd_prop, .groups = "drop")

cld_pigment_plot_data <- cld_pigment_data %>%
  left_join(label_pos_pigment, by = c("community", "pigment"))

p_pigments <- ggplot(data_long, aes(x = pigment, y = proportion, color = PFT)) +
  stat_summary(fun = mean, geom = "point", size = 2.5, position = position_dodge(width = 0.3)) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
               geom = "errorbar", width = 0.3, linewidth = 0.6, position = position_dodge(width = 0.3)) +
  geom_text(data = cld_pigment_plot_data, aes(x = pigment, y = y_pos, label = cld),
            color = "black", vjust = -0.6, size = 3.5, fontface = "bold", inherit.aes = FALSE) +
  facet_wrap(~ community, ncol = 1, scales = "free_y") +
  scale_color_manual(values = c(
    "Picocyanobacteria" = "#E69F00", "Flagellates" = "#56B4E9",
    "Diatoms" = "#009E73", "Other" = "grey40"
  )) +
  labs(x = "Pigment", y = "Proportion (moyenne +/- 1 ecart-type)", color = "Functional group",
       title = "Pigment distributions within phytoplankton communities",
       subtitle = "Lettres differentes = difference significative (Kruskal-Wallis + Dunn, p < 0.05, correction BH)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_pigments)

ggsave(file.path(out_dir, "distribution_pigments_communautes_kmeans_k6_PFT_stats.png"),
       plot = p_pigments, width = 8, height = 12, units = "in", dpi = 300)

# ============================================================
# 5 - Distribution spatiale des communautes (grille complete)
# ============================================================

map_mode <- "dominant" #"single_date"

for (map_date in dates[1]) {
  map_date <- as.Date(map_date)
  
  if (map_mode == "dominant") {
    map_data <- data_pca %>%
      dplyr::count(lon, lat, community) %>%
      dplyr::group_by(lon, lat) %>%
      dplyr::slice_max(n, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup()
    map_title <- "Communaute phytoplanctonique dominante par pixel (toutes dates confondues)"
    date_str <- "all_dates"
  } else {
    map_data <- data_pca %>% dplyr::filter(.data$date == map_date)
    map_title <- paste("Communaute phytoplanctonique -", format(map_date, "%Y-%m-%d"))
    date_str <- format(map_date, "%Y-%m-%d")
  }
  
  world <- ne_countries(scale = "medium", returnclass = "sf")
  
  p_map <- ggplot() +
    geom_raster(data = map_data, aes(x = lon, y = lat, fill = community)) +
    geom_sf(data = world, fill = "grey85", color = "grey50", linewidth = 0.2) +
    coord_sf(xlim = range(lons), ylim = range(lats)) +
    scale_fill_manual(values = community_colors, name = "Communaute") +
    labs(x = "Longitude", y = "Latitude", title = map_title) +
    theme_minimal()
  
  print(p_map)
  
  ggsave(file.path(out_dir, paste0(date_str, "_map_pigments_communautes_kmeans_k6_PFT.png")),
         plot = p_map, width = 8, height = 8, units = "in", dpi = 300)
}

# ============================================================
# 6 - Carte de purete de la communaute dominante (toutes dates)
# ============================================================

map_data_purity <- data_pca %>%
  dplyr::count(lon, lat, community, name = "n") %>%
  dplyr::group_by(lon, lat) %>%
  dplyr::mutate(total = sum(n)) %>%
  dplyr::slice_max(n, n = 1, with_ties = FALSE) %>%
  dplyr::mutate(purity = n / total) %>%
  dplyr::ungroup()

world <- ne_countries(scale = "medium", returnclass = "sf")

p_map_purity <- ggplot() +
  geom_raster(data = map_data_purity, aes(x = lon, y = lat, fill = community, alpha = purity)) +
  geom_sf(data = world, fill = "grey85", color = "grey50", linewidth = 0.2) +
  coord_sf(xlim = range(lons), ylim = range(lats)) +
  scale_fill_manual(values = community_colors, name = "Communaute\ndominante") +
  scale_alpha_continuous(range = c(0.25, 1), limits = c(1 / n_cluster, 1), name = "Purete") +
  labs(x = "Longitude", y = "Latitude",
       title = "Purete de la communaute dominante (proportion de dates ou elle est observee)") +
  theme_minimal()

print(p_map_purity)

ggsave(file.path(out_dir, "map_purity_communautes_kmeans_k6_PFT.png"),
       plot = p_map_purity, width = 8, height = 8, units = "in", dpi = 300)

# ============================================================
# Sauvegarde des clusters sous forme de matrice [date, lon, lat]
# ============================================================

cluster_array <- array(
  NA_integer_,
  dim = dim(valid_mask),
  dimnames = list(
    date = as.character(dates),
    lon  = as.character(lons),
    lat  = as.character(lats)
  )
)

cluster_array[valid_idx] <- as.integer(data_pca$community)

dim(cluster_array)
table(cluster_array, useNA = "ifany")

out_dir_clusters <- "F:/data_elise/commu_phyto"
if (!dir.exists(out_dir_clusters)) dir.create(out_dir_clusters, recursive = TRUE)

saveRDS(cluster_array, file.path(out_dir_clusters, "PCA_kmeans_clusters_array.rds"))
saveRDS(lons,  file.path(out_dir_clusters, "lons_PCA_kmeans_clusters_array.rds"))
saveRDS(lats,  file.path(out_dir_clusters, "lats_PCA_kmeans_clusters_array.rds"))
saveRDS(dates, file.path(out_dir_clusters, "dates_PCA_kmeans_clusters_array.rds"))

# ============================================================
# Distrib NASC par cluster de pigment (communaute phytoplanctonique)
# ============================================================

compute_interannual_stats_generic <- function(dat, value_col, entity_col,
                                              years_keep = c("2018","2021","2022","2023"),
                                              log_transform = TRUE, y_limits = NULL) {
  dat_sub <- dat %>%
    dplyr::filter(year %in% years_keep) %>%
    dplyr::mutate(val = if (log_transform) log10(.data[[value_col]]) else .data[[value_col]])
  
  daily_means <- dat_sub %>%
    dplyr::group_by(.data[[entity_col]], date, year) %>%
    dplyr::summarise(val_day = mean(val, na.rm = TRUE), .groups = "drop")
  
  summary_stats <- daily_means %>%
    dplyr::group_by(.data[[entity_col]], year) %>%
    dplyr::summarise(n_days = dplyr::n(), m = mean(val_day, na.rm = TRUE),
                     var_v = var(val_day, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(
      var_v  = tidyr::replace_na(var_v, 0), sd_v = sqrt(var_v),
      mean_c = if (log_transform) 10^m else m,
      ymin   = if (log_transform) 10^(m - sd_v) else m - sd_v,
      ymax   = if (log_transform) 10^(m + sd_v) else m + sd_v
    )
  
  letters_by_entity <- daily_means %>%
    dplyr::group_by(.data[[entity_col]]) %>%
    dplyr::group_split() %>%
    purrr::map_dfr(function(dat_e) {
      this_e <- as.character(unique(dat_e[[entity_col]]))
      if (dplyr::n_distinct(dat_e$year) < 2) {
        out_e <- data.frame(year = as.character(unique(dat_e$year)), letter = "a")
      } else {
        model   <- lm(val_day ~ year, data = dat_e)
        emm     <- emmeans::emmeans(model, ~ year)
        cld_out <- multcomp::cld(emm, Letters = letters, adjust = "sidak", sort = FALSE)
        out_e <- data.frame(year = as.character(cld_out$year), letter = trimws(cld_out$.group))
      }
      out_e[[entity_col]] <- this_e
      out_e
    })
  
  out <- summary_stats %>%
    dplyr::mutate(!!entity_col := as.character(.data[[entity_col]]), year = as.character(year)) %>%
    dplyr::left_join(letters_by_entity, by = c("year", entity_col))
  
  if (log_transform) {
    out <- out %>% dplyr::mutate(y = 10^(log10(ymax) + 1))
  } else {
    top    <- if (!is.null(y_limits)) y_limits[2] else max(dat_sub[[value_col]], na.rm = TRUE)
    bottom <- if (!is.null(y_limits)) y_limits[1] else min(dat_sub[[value_col]], na.rm = TRUE)
    out <- out %>% dplyr::mutate(y = bottom + 0.94 * (top - bottom))
  }
  out
}

plot_interannual_generic <- function(dat, value_col, entity_col, label,
                                     save_dir = out_dir, save = TRUE,
                                     years_keep = c("2018","2021","2022","2023"),
                                     log_transform = TRUE, y_limits = NULL) {
  stats_df <- compute_interannual_stats_generic(dat, value_col, entity_col, years_keep, log_transform, y_limits)
  
  dat_all <- dat %>%
    dplyr::filter(year %in% years_keep) %>%
    dplyr::mutate(!!entity_col := as.character(.data[[entity_col]]), year = as.character(year)) %>%
    dplyr::left_join(stats_df %>% dplyr::select(year, dplyr::all_of(entity_col), mean_c), by = c("year", entity_col))
  
  dat_all  <- dat_all  %>% dplyr::mutate(year = factor(year, levels = years_keep))
  stats_df <- stats_df %>% dplyr::mutate(year = factor(year, levels = years_keep))
  
  p <- ggplot(dat_all, aes(x = year, y = .data[[value_col]])) +
    geom_violin(trim = TRUE, scale = "width", alpha = 0.8, linewidth = 0.3, fill = "grey80") +
    geom_errorbar(data = stats_df, aes(x = year, y = mean_c, ymin = ymin, ymax = ymax),
                  inherit.aes = FALSE, width = 0.12, linewidth = 0.5, color = "black") +
    geom_crossbar(data = stats_df, aes(x = year, y = mean_c, ymin = mean_c, ymax = mean_c),
                  inherit.aes = FALSE, width = 0.35, color = "black", linewidth = 0.4, fatten = 1) +
    geom_text(data = stats_df, aes(x = year, y = y, label = letter),
              inherit.aes = FALSE, size = 3.5, fontface = "bold") +
    facet_wrap(as.formula(paste("~", entity_col))) +
    labs(x = "Year", y = label, title = paste(label, "distribution across years"),
         caption = if (log_transform) {
           "Tiret = moyenne geometrique ; barre = +/- ecart-type (log10) ; lettres = groupes Sidak (p < 0.05)"
         } else {
           "Tiret = moyenne ; barre = +/- ecart-type ; lettres = groupes Sidak (p < 0.05)"
         }) +
    theme_classic() +
    theme(axis.title = element_text(size = 12), axis.text = element_text(size = 11),
          axis.text.x = element_text(angle = 45, hjust = 1), strip.text = element_text(size = 12, face = "bold"),
          plot.title = element_text(size = 14, face = "bold"), panel.spacing = unit(1.2, "lines"))
  
  if (log_transform) {
    p <- p + scale_y_log10(expand = expansion(mult = c(0.05, 0.45)))
  } else if (!is.null(y_limits)) {
    p <- p + coord_cartesian(ylim = y_limits, clip = "off")
  }
  
  if (save) {
    dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
    fname <- gsub("[^A-Za-z0-9_]", "_", label)
    ggsave(file.path(save_dir, paste0("violin_interannual_", fname, ".png")), plot = p, width = 12, height = 8, dpi = 300)
    write.csv(stats_df, file.path(save_dir, paste0("stats_interannual_", fname, ".csv")), row.names = FALSE)
  }
  p
}

# Necessite pca, km, pig_vars, community_colors (script PCA/k-means).
years_keep <- c("2018", "2021", "2022", "2023")
freq <- 38
nasc_pig <- readRDS(paste0(
  "F:/data_elise/ds_NASC_pig_ftle_fod/ds_NASC_per_esu_all/new_ratio_ds_NASC_per_esu_pig_ftle_fod_2018_2021_2022_2023_transect_",
  freq, "kHz_mask9.rds"
))
str(nasc_pig)
# ------------------------------------------------------------
# Correspondance entre les noms all_ds (pig_vars, cles) et les
# colonnes du dataset NASC (valeurs) : ce dernier a ete regenere
# avec la meme nomenclature "_totpig" (cf. adaptation du dataset
# NASC : Chla_total + ratio pigment/somme(hors Chla)).
# ------------------------------------------------------------
pig_map <- c(
  chla_total    = "Chla_total",
  per_totpig    = "Per_totpig",
  but_totpig    = "But_totpig",
  fuco_totpig   = "Fuco_totpig",
  hex_totpig    = "Hex_totpig",
  allo_totpig   = "Allo_totpig",
  zea_totpig    = "Zea_totpig",
  chlb_totpig   = "Chlb_totpig",
  dvchla_totpig = "DvChla_totpig"
)

stopifnot(all(pig_map %in% names(nasc_pig)))
stopifnot(identical(names(pig_map), pig_vars))  # meme ordre que la PCA

nasc_for_pca <- nasc_pig[, unname(pig_map)]
colnames(nasc_for_pca) <- names(pig_map)

valid_rows <- rowSums(is.finite(as.matrix(nasc_for_pca))) == length(pig_vars)
pig_sum    <- rep(NA_real_, nrow(nasc_for_pca))
pig_sum[valid_rows] <- rowSums(nasc_for_pca[valid_rows, ])

keep <- valid_rows & !is.na(pig_sum) & (pig_sum > 0)

nasc_prop <- as.matrix(nasc_for_pca[keep, ]) / pig_sum[keep]

pc_scores_nasc <- scale(nasc_prop, center = pca$center, scale = pca$scale) %*% pca$rotation

nasc_pig$community <- NA_character_
nasc_pig$community[keep] <- as.character(assign_cluster(pc_scores_nasc, km$centers))
nasc_pig$community <- factor(nasc_pig$community, levels = names(community_colors))

cat("NASC avec communaute assignee :", sum(!is.na(nasc_pig$community)), "/", nrow(nasc_pig), "\n")

nasc_long_comm <- nasc_pig %>%
  dplyr::mutate(date = as.Date(time_nasc), year = format(time_nasc, "%Y"), variable = "NASC") %>%
  dplyr::filter(is.finite(nasc), nasc > 0, year %in% years_keep, !is.na(community))

compute_stats_by_group <- function(dat, value_col, entity_col, entity_val, group_col,
                                   years_keep = c("2018","2021","2022","2023"),
                                   log_transform = TRUE, y_limits = NULL) {
  dat_sub <- dat %>%
    dplyr::filter(.data[[entity_col]] == entity_val, year %in% years_keep) %>%
    dplyr::mutate(val = if (log_transform) log10(.data[[value_col]]) else .data[[value_col]])
  
  daily_means <- dat_sub %>%
    dplyr::group_by(date, year, .data[[group_col]]) %>%
    dplyr::summarise(val_day = mean(val, na.rm = TRUE), .groups = "drop")
  
  summary_stats <- daily_means %>%
    dplyr::group_by(year, .data[[group_col]]) %>%
    dplyr::summarise(n_days = dplyr::n(), m = mean(val_day, na.rm = TRUE),
                     var_v = var(val_day, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(
      var_v  = tidyr::replace_na(var_v, 0), sd_v = sqrt(var_v),
      mean_c = if (log_transform) 10^m else m,
      ymin   = if (log_transform) 10^(m - sd_v) else m - sd_v,
      ymax   = if (log_transform) 10^(m + sd_v) else m + sd_v
    )
  
  letters_by_group <- daily_means %>%
    dplyr::group_by(.data[[group_col]]) %>%
    dplyr::group_split() %>%
    purrr::map_dfr(function(dat_g) {
      this_g <- as.character(unique(dat_g[[group_col]]))
      if (dplyr::n_distinct(dat_g$year) < 2) {
        res <- data.frame(year = as.character(unique(dat_g$year)), letter = "a")
      } else {
        model   <- lm(val_day ~ year, data = dat_g)
        emm     <- emmeans::emmeans(model, ~ year)
        cld_out <- multcomp::cld(emm, Letters = letters, adjust = "sidak", sort = FALSE)
        res <- data.frame(year = as.character(cld_out$year), letter = trimws(cld_out$.group))
      }
      res[[group_col]] <- this_g
      res
    })
  
  out <- summary_stats %>%
    dplyr::mutate(!!group_col := as.character(.data[[group_col]]), year = as.character(year)) %>%
    dplyr::left_join(letters_by_group, by = c("year", group_col))
  
  if (log_transform) {
    out <- out %>% dplyr::mutate(y = 10^(log10(ymax) + 1.5))
  } else {
    range_by_group <- dat_sub %>%
      dplyr::group_by(.data[[group_col]]) %>%
      dplyr::summarise(val_min = min(.data[[value_col]], na.rm = TRUE),
                       val_max = max(.data[[value_col]], na.rm = TRUE), .groups = "drop") %>%
      dplyr::mutate(!!group_col := as.character(.data[[group_col]]))
    out <- out %>%
      dplyr::left_join(range_by_group, by = group_col) %>%
      dplyr::mutate(y = if (!is.null(y_limits)) y_limits[1] + 0.94 * (y_limits[2] - y_limits[1])
                    else val_max + 0.03 * (val_max - val_min))
  }
  out
}

plot_by_group <- function(dat, value_col, entity_col, entity_val, label, group_col, palette,
                          save_dir = out_dir, save = TRUE,
                          years_keep = c("2018","2021","2022","2023"),
                          log_transform = TRUE, y_limits = NULL) {
  stats_df <- compute_stats_by_group(dat, value_col, entity_col, entity_val, group_col,
                                     years_keep, log_transform, y_limits)
  
  dat_sub <- dat %>%
    dplyr::filter(.data[[entity_col]] == entity_val, year %in% years_keep) %>%
    dplyr::mutate(!!group_col := as.character(.data[[group_col]]), year = as.character(year)) %>%
    dplyr::left_join(stats_df %>% dplyr::select(year, dplyr::all_of(group_col), mean_c),
                     by = c("year", group_col))
  
  dat_sub  <- dat_sub  %>% dplyr::mutate(year = factor(year, levels = years_keep),
                                         !!group_col := droplevels(factor(.data[[group_col]])))
  stats_df <- stats_df %>% dplyr::mutate(year = factor(year, levels = years_keep),
                                         !!group_col := factor(.data[[group_col]], levels = levels(dat_sub[[group_col]])))
  
  strip_colors <- palette[levels(dat_sub[[group_col]])]
  
  p <- ggplot(dat_sub, aes(x = year, y = .data[[value_col]])) +
    geom_violin(trim = TRUE, scale = "width", alpha = 0.8, linewidth = 0.3, fill = "grey75") +
    geom_errorbar(data = stats_df, aes(x = year, y = mean_c, ymin = ymin, ymax = ymax),
                  inherit.aes = FALSE, width = 0.12, linewidth = 0.5, color = "black") +
    geom_crossbar(data = stats_df, aes(x = year, y = mean_c, ymin = mean_c, ymax = mean_c),
                  inherit.aes = FALSE, width = 0.35, color = "black", linewidth = 0.2) +
    geom_text(data = stats_df, aes(x = year, y = y, label = letter),
              inherit.aes = FALSE, size = 3.5, fontface = "bold") +
    ggh4x::facet_wrap2(
      vars(!!rlang::sym(group_col)), scales = "free_y",
      strip = ggh4x::strip_themed(
        background_x = ggh4x::elem_list_rect(fill = strip_colors),
        text_x       = ggh4x::elem_list_text(colour = "white", face = "bold")
      )
    ) +
    labs(x = "Year", y = label,
         title = paste(label, "distribution within", group_col, "clusters across years"),
         caption = "Tiret = moyenne ; barre = +/- ecart-type ; lettres = groupes Sidak (p < 0.05)") +
    theme_classic() +
    theme(axis.title = element_text(size = 12), axis.text = element_text(size = 11),
          axis.text.x = element_text(angle = 45, hjust = 1), strip.text = element_text(size = 12, face = "bold"),
          plot.title = element_text(size = 14, face = "bold"), panel.spacing = unit(1.2, "lines"))
  
  if (log_transform) {
    p <- p + scale_y_log10(expand = expansion(mult = c(0.05, 0.45)))
  } else if (!is.null(y_limits)) {
    p <- p + coord_cartesian(ylim = y_limits, clip = "off")
  }
  
  if (save) {
    dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
    ggsave(file.path(save_dir, paste0("violin_", label, ".png")), plot = p, width = 10, height = 6.5, dpi = 300)
    write.csv(stats_df %>% dplyr::mutate(variable = label),
              file.path(save_dir, paste0("stats_", label, ".csv")), row.names = FALSE)
  }
  p
}

plot_by_group(nasc_long_comm, value_col = "nasc", entity_col = "variable",
              entity_val = "NASC", label = paste0("NASC (", freq, " kHz) par communaute phytoplanctonique"),
              group_col = "community", palette = community_colors,
              years_keep = years_keep, log_transform = TRUE)

# ============================================================
# NASC par communaute, par annee (comparaison entre communautes
# au sein de chaque annee)
# ============================================================

compute_stats_by_year <- function(dat, value_col, entity_col, entity_val, group_col,
                                  years_keep = c("2018","2021","2022","2023"),
                                  log_transform = TRUE, y_limits = NULL) {
  dat_sub <- dat %>%
    dplyr::filter(.data[[entity_col]] == entity_val, year %in% years_keep) %>%
    dplyr::mutate(val = if (log_transform) log10(.data[[value_col]]) else .data[[value_col]])
  
  daily_means <- dat_sub %>%
    dplyr::group_by(date, year, .data[[group_col]]) %>%
    dplyr::summarise(val_day = mean(val, na.rm = TRUE), .groups = "drop")
  
  summary_stats <- daily_means %>%
    dplyr::group_by(year, .data[[group_col]]) %>%
    dplyr::summarise(n_days = dplyr::n(), m = mean(val_day, na.rm = TRUE),
                     var_v = var(val_day, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(
      var_v  = tidyr::replace_na(var_v, 0), sd_v = sqrt(var_v),
      mean_c = if (log_transform) 10^m else m,
      ymin   = if (log_transform) 10^(m - sd_v) else m - sd_v,
      ymax   = if (log_transform) 10^(m + sd_v) else m + sd_v
    )
  
  letters_by_year <- daily_means %>%
    dplyr::group_by(year) %>%
    dplyr::group_split() %>%
    purrr::map_dfr(function(dat_y) {
      this_y <- as.character(unique(dat_y$year))
      if (dplyr::n_distinct(dat_y[[group_col]]) < 2) {
        res <- data.frame(group_val = as.character(unique(dat_y[[group_col]])), letter = "a")
        names(res)[1] <- group_col
      } else {
        form_lm <- stats::as.formula(paste("val_day ~", group_col))
        model   <- lm(form_lm, data = dat_y)
        emm     <- emmeans::emmeans(model, as.formula(paste("~", group_col)))
        cld_out <- multcomp::cld(emm, Letters = letters, adjust = "sidak", sort = FALSE)
        res <- data.frame(group_val = as.character(cld_out[[group_col]]), letter = trimws(cld_out$.group))
        names(res)[1] <- group_col
      }
      res$year <- this_y
      res
    })
  
  out <- summary_stats %>%
    dplyr::mutate(!!group_col := as.character(.data[[group_col]]), year = as.character(year)) %>%
    dplyr::left_join(letters_by_year, by = c("year", group_col))
  
  if (log_transform) {
    out <- out %>% dplyr::mutate(y = 10^(log10(ymax) + 1.5))
  } else {
    range_by_year <- dat_sub %>%
      dplyr::group_by(year) %>%
      dplyr::summarise(val_min = min(.data[[value_col]], na.rm = TRUE),
                       val_max = max(.data[[value_col]], na.rm = TRUE), .groups = "drop") %>%
      dplyr::mutate(year = as.character(year))
    out <- out %>%
      dplyr::left_join(range_by_year, by = "year") %>%
      dplyr::mutate(y = if (!is.null(y_limits)) y_limits[1] + 0.94 * (y_limits[2] - y_limits[1])
                    else val_max + 0.03 * (val_max - val_min))
  }
  out
}

plot_by_year <- function(dat, value_col, entity_col, entity_val, label, group_col, palette,
                         save_dir = out_dir, save = TRUE,
                         years_keep = c("2018","2021","2022","2023"),
                         log_transform = TRUE, y_limits = NULL) {
  stats_df <- compute_stats_by_year(dat, value_col, entity_col, entity_val, group_col,
                                    years_keep, log_transform, y_limits)
  
  dat_sub <- dat %>%
    dplyr::filter(.data[[entity_col]] == entity_val, year %in% years_keep) %>%
    dplyr::mutate(!!group_col := as.character(.data[[group_col]]), year = as.character(year)) %>%
    dplyr::left_join(stats_df %>% dplyr::select(year, dplyr::all_of(group_col), mean_c),
                     by = c("year", group_col))
  
  group_levels <- names(palette)
  dat_sub  <- dat_sub  %>% dplyr::mutate(!!group_col := factor(.data[[group_col]], levels = group_levels))
  stats_df <- stats_df %>% dplyr::mutate(!!group_col := factor(.data[[group_col]], levels = group_levels))
  
  dat_sub  <- dat_sub  %>% dplyr::mutate(year = factor(year, levels = years_keep))
  stats_df <- stats_df %>% dplyr::mutate(year = factor(year, levels = years_keep))
  
  p <- ggplot(dat_sub, aes(x = .data[[group_col]], y = .data[[value_col]], fill = .data[[group_col]])) +
    geom_violin(trim = TRUE, scale = "width", alpha = 0.8, linewidth = 0.3) +
    geom_errorbar(data = stats_df, aes(x = .data[[group_col]], y = mean_c, ymin = ymin, ymax = ymax),
                  inherit.aes = FALSE, width = 0.12, linewidth = 0.5, color = "black") +
    geom_crossbar(data = stats_df, aes(x = .data[[group_col]], y = mean_c, ymin = mean_c, ymax = mean_c),
                  inherit.aes = FALSE, width = 0.35, color = "black", linewidth = 0.2) +
    geom_text(data = stats_df, aes(x = .data[[group_col]], y = y, label = letter),
              inherit.aes = FALSE, size = 3.5, fontface = "bold") +
    scale_fill_manual(values = palette, name = "Communaute") +
    facet_wrap(~ year) +
    labs(x = "Communaute", y = label,
         title = paste(label, "par communaute phytoplanctonique, par annee"),
         caption = "Tiret = moyenne ; barre = +/- ecart-type ; lettres = groupes Sidak (p < 0.05), comparaison entre communautes au sein de chaque annee") +
    theme_classic() +
    theme(axis.title = element_text(size = 12), axis.text = element_text(size = 11),
          strip.text = element_text(size = 12, face = "bold"),
          plot.title = element_text(size = 14, face = "bold"), panel.spacing = unit(1.2, "lines"))
  
  if (log_transform) {
    p <- p + scale_y_log10(expand = expansion(mult = c(0.05, 0.45)))
  } else if (!is.null(y_limits)) {
    p <- p + coord_cartesian(ylim = y_limits, clip = "off")
  }
  
  if (save) {
    dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
    ggsave(file.path(save_dir, paste0("violin_", label, ".png")), plot = p, width = 11, height = 7, dpi = 300)
    write.csv(stats_df %>% dplyr::mutate(variable = label),
              file.path(save_dir, paste0("stats_", label, ".csv")), row.names = FALSE)
  }
  p
}

plot_by_year(nasc_long_comm, value_col = "nasc", entity_col = "variable",
             entity_val = "NASC", label = paste0("NASC (", freq, " kHz) par communaute phytoplanctonique par annee"),
             group_col = "community", palette = community_colors,
             years_keep = years_keep, log_transform = TRUE)

# ============================================================
# Effectifs NASC par communaute phytoplanctonique, par annee
# ============================================================

compute_n_nasc_by_community <- function(dat, value_col, entity_col, entity_val, group_col,
                                        years_keep = c("2018","2021","2022","2023")) {
  dat %>%
    dplyr::filter(.data[[entity_col]] == entity_val, year %in% years_keep,
                  is.finite(.data[[value_col]])) %>%
    dplyr::count(year, .data[[group_col]], name = "n_obs")
}

plot_n_nasc_by_community <- function(dat, value_col, entity_col, entity_val, label, group_col, palette,
                                     save_dir = out_dir, save = TRUE,
                                     years_keep = c("2018","2021","2022","2023")) {
  
  n_df <- compute_n_nasc_by_community(dat, value_col, entity_col, entity_val, group_col, years_keep) %>%
    dplyr::mutate(
      year = factor(year, levels = years_keep),
      !!group_col := factor(.data[[group_col]], levels = names(palette))
    )
  
  p <- ggplot(n_df, aes(x = year, y = n_obs, fill = .data[[group_col]])) +
    geom_col(position = position_dodge2(preserve = "single"), color = "black", linewidth = 0.2) +
    geom_text(aes(label = n_obs), position = position_dodge2(width = 0.9, preserve = "single"),
              vjust = -0.4, size = 3) +
    scale_fill_manual(values = palette, name = "Communaute") +
    labs(
      x = "Year", y = "Nombre de mesures NASC",
      title = paste(label, "- effectifs NASC par communaute phytoplanctonique, par annee")
    ) +
    theme_classic() +
    theme(
      axis.title  = element_text(size = 12),
      axis.text   = element_text(size = 11),
      plot.title  = element_text(size = 13, face = "bold")
    )
  
  
  dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(save_dir, paste0("n_nasc_", label, ".png")), plot = p, width = 9, height = 5.5, dpi = 300)
  write.csv(n_df, file.path(save_dir, paste0("n_nasc_", label, ".csv")), row.names = FALSE)

  
  p
}

plot_n_nasc_by_community(nasc_long_comm, value_col = "nasc", entity_col = "variable",
                         entity_val = "NASC",
                         label = paste0("NASC (", freq, " kHz) par communaute phytoplanctonique"),
                         group_col = "community", palette = community_colors,
                         years_keep = years_keep)
