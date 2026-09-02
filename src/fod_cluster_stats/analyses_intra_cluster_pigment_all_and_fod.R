# ============================================================
# Etude des tendances des pigments (2018, 2021, 2023) au sein
# des clusters FOD - VERSION CORRIGEE
#
# On utilise TOUTES les donnees pigmentaires d'un cluster,
# pas seulement les points co-localises avec le nasc.
#
# 3 graphiques demandes :
#  1. Proportion des pigments pour chaque cluster, chaque annee
#  2. Evolution des moyennes de concentration des pigments,
#     par annee et par cluster FOD
#  3. Evolution de la proportion des pigments par rapport a la
#     chlorophylle totale, par annee et par cluster FOD
# ============================================================

# libraries -----------------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(emmeans)
library(multcomp)
library(multcompView)
library(purrr)
library(tidytext) 
# Global variables ------------------------------------------------------------
path_pig  <- "F:/data_elise/pigmeann/pigments_2018_2021_2022_2023_crop.rds"
path_fod  <- "F:/data_elise/fod_elise_2018_2021_2022_2023/mclust_results_2018_2021_2022_2023/transitions_upgraded/cluster_transition_map_renamed.rds"
path_time <- "F:/data_elise/fod_elise_2018_2021_2022_2023/mclust_results_2018_2021_2022_2023/time.rds"
path_lat  <- "F:/data_elise/fod_elise_2018_2021_2022_2023/mclust_results_2018_2021_2022_2023/lat.rds"
path_lon  <- "F:/data_elise/fod_elise_2018_2021_2022_2023/mclust_results_2018_2021_2022_2023/lon.rds"
out_dir <- "C:/Users/mmolinet/elisou_ta_stagiaire_pref/MIO_internship_III/figures/distribution_per_fod_cluster/violin_plots_colors"
# ouvrir les donnees ----------------------------------------------------------
pig  <- readRDS(path_pig)
str(pig)
fod  <- readRDS(path_fod)
time <- readRDS(path_time)
lat  <- readRDS(path_lat)
lon  <- readRDS(path_lon)

# Palette de couleurs pour les clusters et transitions FOD -------------------
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
# Palette combinee, utilisable pour tous les codes FOD (clusters + transitions)
fod_cols <- c(cluster_cols, transition_cols)

# # ============================================================
# # 1. Correspondance temporelle : dates FOD <-> dates pigments
# # ============================================================
# # BUG CORRIGE : intersect() renvoie un simple vecteur de dates,
# # pas un objet avec des colonnes i_pig / i_fod. On construit ici
# # une vraie table de correspondance avec les indices des deux
# # tableaux via une jointure.
# 
# fod_dates <- as.Date(time)
# pig_dates <- as.Date(pig$date)
# 
# date_match <- inner_join(
#   data.frame(i_fod = seq_along(fod_dates), date = fod_dates),
#   data.frame(i_pig = seq_along(pig_dates), date = pig_dates),
#   by = "date"
# )
# 
# cat("Dates FOD     :", length(fod_dates), "\n")
# cat("Dates pigment :", length(pig_dates), "\n")
# cat("Dates communes:", nrow(date_match), "\n")   # remplace length(dates_communes)
# 
# # ============================================================
# # 2. Correspondance spatiale : grille pigment -> grille FOD
# # ============================================================
# 
# # Pour chaque longitude pigmentaire, trouver la longitude FOD la plus proche
# idx_lon <- sapply(pig$lon, function(x) which.min(abs(lon - x)))
# 
# # Pour chaque latitude pigmentaire, trouver la latitude FOD la plus proche
# idx_lat <- sapply(pig$lat, function(x) which.min(abs(lat - x)))
# 
# # ============================================================
# # 3. Construire le tableau pigment + FOD pour toutes les dates
# #    communes -- VERSION VECTORISEE (sans boucle for)
# # ============================================================
# 
# n_lon  <- length(pig$lon)
# n_lat  <- length(pig$lat)
# n_date <- nrow(date_match)
# 
# # ------------------------------------------------------------
# # 3.1. Regriddage de FOD sur la grille pigments, pour TOUTES les
# #      dates communes d'un coup (indexation par array 3D)
# # ------------------------------------------------------------
# # fod a pour dimensions [lon_fod, lat_fod, date_fod] ; on sélectionne
# # les indices lon/lat regrillés (idx_lon, idx_lat) et les dates
# # communes (date_match$i_fod), dans cet ordre -> résultat de
# # dimensions [n_lon, n_lat, n_date], déjà aligné sur la grille pigments.
# 
# fod_regridded <- fod[idx_lon, idx_lat, date_match$i_fod]
# 
# # ------------------------------------------------------------
# # 3.2. Extraction des pigments pour les dates communes, en gardant
# #      la dimension [date, lon, lat] d'origine de pig$c_cond_XXX
# # ------------------------------------------------------------
# 
# extract_pigment <- function(pig_array, i_pig_vec) {
#   # pig_array : [date, lon, lat] -> on garde seulement les dates communes
#   pig_array[i_pig_vec, , , drop = FALSE]
# }
# 
# pig_Chla   <- extract_pigment(pig$c_cond_Chla,   date_match$i_pig)
# pig_Per    <- extract_pigment(pig$c_cond_Per,    date_match$i_pig)
# pig_But    <- extract_pigment(pig$c_cond_But,    date_match$i_pig)
# pig_Fuco   <- extract_pigment(pig$c_cond_Fuco,   date_match$i_pig)
# pig_Hex    <- extract_pigment(pig$c_cond_Hex,    date_match$i_pig)
# pig_Allo   <- extract_pigment(pig$c_cond_Allo,   date_match$i_pig)
# pig_Zea    <- extract_pigment(pig$c_cond_Zea,    date_match$i_pig)
# pig_Chlb   <- extract_pigment(pig$c_cond_Chlb,   date_match$i_pig)
# pig_DvChla <- extract_pigment(pig$c_cond_DvChla, date_match$i_pig)
# 
# # ------------------------------------------------------------
# # 3.3. Réorganiser tout en [lon, lat, date] pour matcher fod_regridded,
# #      puis aplatir en un seul data.frame long
# # ------------------------------------------------------------
# # aperm() transpose [date, lon, lat] -> [lon, lat, date], cohérent
# # avec fod_regridded ; as.vector() aplatit ensuite en respectant
# # l'ordre lon (le plus vite), puis lat, puis date -- exactement
# # l'ordre produit par expand.grid ci-dessous.
# 
# to_long_vector <- function(arr) as.vector(aperm(arr, c(2, 3, 1)))
# 
# pig_fod <- expand.grid(
#   lon_idx  = seq_len(n_lon),
#   lat_idx  = seq_len(n_lat),
#   date_idx = seq_len(n_date)
# ) %>%
#   mutate(
#     lon    = pig$lon[lon_idx],
#     lat    = pig$lat[lat_idx],
#     date   = pig$date[date_match$i_pig[date_idx]],
#     # BUG CORRIGE : "year = format(date, '%Y')" etait ecrit comme une
#     # instruction isolee AVANT la definition de plot_pigment(), hors de
#     # tout mutate(). Ca ne creait pas de colonne dans pig_fod / pig_long :
#     # ca creait un objet global "year" (longueur 1, voire un conflit avec
#     # la fonction base::date()), sans rapport avec les lignes du tableau.
#     # ggplot cherchait alors "year" dans l'environnement global plutot que
#     # dans les donnees -> mismatch de longueur avec les ~12.4M de lignes
#     # ("Aesthetics must be either length 1 or the same as the data").
#     # Fix : creer year ICI, dans le mutate(), comme vraie colonne alignee
#     # ligne par ligne avec le reste du tableau.
#     year   = format(date, "%Y"),
#     
#     Chla   = to_long_vector(pig_Chla),
#     Perid  = to_long_vector(pig_Per),
#     But    = to_long_vector(pig_But),
#     Fuco   = to_long_vector(pig_Fuco),
#     Hex    = to_long_vector(pig_Hex),
#     Allo   = to_long_vector(pig_Allo),
#     Zeax   = to_long_vector(pig_Zea),
#     Chlb   = to_long_vector(pig_Chlb),
#     DvChla = to_long_vector(pig_DvChla),
#     
#     fod    = factor(as.vector(fod_regridded))
#   ) %>%
#   dplyr::select(-lon_idx, -lat_idx, -date_idx)
# 
# str(pig_fod)
# saveRDS(pig_fod, "C:/Users/mmolinet/elisou_ta_stagiaire_pref/MIO_internship_III/src/analyses_intra_cluster/pig_fod_associated")

pig_fod <- readRDS("F:/data_elise/distrib_intra_cluster_fod/pig_fod_associated")
# ============================================================
# 4. Format long : une ligne par pixel x pigment
# ============================================================

pigments_plot <- c("Chla", "Perid", "But", "Fuco", "Hex", "Allo", "Zeax", "Chlb", "DvChla")

pig_long <- pig_fod %>%
  pivot_longer(
    cols      = all_of(pigments_plot),
    names_to  = "pigment",
    values_to = "concentration"
  ) %>%
  filter(concentration > 0)   # log10 impossible sur des valeurs <= 0

# ============================================================
# 5. Distribution des concentrations par annee x FOD (boxplots)
# ============================================================
# statistiques 


compute_pigment_stats <- function(pigment_name, dat = pig_long,
                                  years_keep = c("2018", "2021", "2022", "2023")) {
  
  dat_pig <- dat %>%
    dplyr::filter(pigment == pigment_name, year %in% years_keep)
  
  # ---- Agrégation par DATE x FOD d'abord (une ligne = une date) ----
  # réduit la pseudo-réplication liée à l'autocorrélation spatiale des pixels
  daily_means <- dat_pig %>%
    dplyr::mutate(log_conc = log10(concentration)) %>%
    dplyr::group_by(date, year, fod) %>%
    dplyr::summarise(log_conc_day = mean(log_conc, na.rm = TRUE), .groups = "drop")
  
  # ---- Résumé par année (basé sur les moyennes journalières, pas les pixels) ----
  summary_stats <- daily_means %>%
    dplyr::group_by(year, fod) %>%
    dplyr::summarise(
      n_days   = n(),
      log_mean = mean(log_conc_day, na.rm = TRUE),
      var_log  = var(log_conc_day,  na.rm = TRUE),
      .groups  = "drop"
    ) %>%
    dplyr::mutate(
      var_log = tidyr::replace_na(var_log, 0),  # BUGFIX: var() vaut NA si n_days == 1
      # CORRECTIF (2) : la légende annonce "+/- écart-type", donc on utilise
      # sqrt(var_log) et non var_log brute (qui n'a pas la même échelle et
      # produisait des barres d'erreur incohérentes avec le texte affiché).
      sd_log = sqrt(var_log),
      mean_c = 10^log_mean,
      ymin   = 10^(log_mean - sd_log),
      ymax   = 10^(log_mean + sd_log)
    )
  
  # ---- Test sur les moyennes journalières (1 valeur par date, pas par pixel) ----
  letters_by_fod <- daily_means %>%
    dplyr::group_by(fod) %>%
    dplyr::group_split() %>%
    purrr::map_dfr(function(dat_fod) {
      this_fod <- as.character(unique(dat_fod$fod))
      
      if (dplyr::n_distinct(dat_fod$year) < 2) {
        return(data.frame(fod = this_fod, year = as.character(unique(dat_fod$year)), letter = "a"))
      }
      
      model   <- lm(log_conc_day ~ year, data = dat_fod)
      emm     <- emmeans::emmeans(model, ~ year)
      cld_out <- multcomp::cld(emm, Letters = letters, adjust = "sidak", sort = FALSE)
      
      data.frame(
        fod    = this_fod,
        year   = as.character(cld_out$year),
        letter = trimws(cld_out$.group)
      )
    })
  
  summary_stats %>%
    dplyr::mutate(fod = as.character(fod), year = as.character(year)) %>%  # BUGFIX: cohérence des types avant le join
    dplyr::left_join(letters_by_fod, by = c("year", "fod")) %>%
    # Décalage FIXE en espace log10 (+1.5 dex) au-dessus de la barre
    # d'erreur : donne un écart visuel constant quelle que soit l'échelle
    # y du panel. Comme le facet_wrap est en échelle fixe (pas free_y),
    # un décalage global est valable pour toutes les facettes.
    dplyr::mutate(y = 10^(log10(ymax) + 1.5))
}


plot_pigment <- function(pigment_name, dat = pig_long, save_dir = out_dir,
                         save = TRUE, years_keep = c("2018", "2021", "2022", "2023")) {
  
  stats_df <- compute_pigment_stats(pigment_name, dat, years_keep)
  
  dat_pig <- dat %>%
    dplyr::filter(pigment == pigment_name, year %in% years_keep) %>%
    dplyr::mutate(fod = as.character(fod), year = as.character(year)) %>%  # BUGFIX: cohérence des types avant le join
    dplyr::left_join(stats_df %>% dplyr::select(year, fod, mean_c), by = c("year", "fod"))
  
  # ---- Ordre chronologique sur l'axe x (plus de tri par concentration) ----
  dat_pig  <- dat_pig  %>% dplyr::mutate(year = factor(year, levels = years_keep))
  stats_df <- stats_df %>% dplyr::mutate(year = factor(year, levels = years_keep))
  
  # ---- Facettes FOD triées par ordre numérique croissant (1, 2, 3...) ----
  # BUGFIX: `fod` est une chaîne de caractères, donc facet_wrap() les
  # triait par ordre alphabétique ("1","10","11","12","13","2","3"...).
  # On force un facteur avec les niveaux triés numériquement ; les
  # valeurs non numériques (ex. "NA") sont placées en dernier.
  fod_levels_sorted <- unique(c(dat_pig$fod, stats_df$fod))
  fod_levels_sorted <- fod_levels_sorted[order(suppressWarnings(as.numeric(fod_levels_sorted)), na.last = TRUE)]
  dat_pig  <- dat_pig  %>% dplyr::mutate(fod = factor(fod, levels = fod_levels_sorted))
  stats_df <- stats_df %>% dplyr::mutate(fod = factor(fod, levels = fod_levels_sorted))
  
  p <- ggplot(dat_pig, aes(x = year, y = concentration, fill = fod)) +
    geom_violin(trim = TRUE, scale = "width", alpha = 0.8, linewidth = 0.3) +
    scale_fill_manual(values = fod_cols, name = "FOD") +
    geom_errorbar(
      data = stats_df,
      aes(x = year, y = mean_c, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE, width = 0.12, linewidth = 0.5, color = "black"
    ) +
    geom_crossbar(
      data = stats_df,
      aes(x = year, y = mean_c, ymin = mean_c, ymax = mean_c),
      inherit.aes = FALSE, width = 0.35, color = "black", linewidth = 0.2
    ) +
    geom_text(
      data = stats_df,
      aes(x = year, y = y, label = letter),
      inherit.aes = FALSE, size = 3.5, fontface = "bold"
    ) +
    # BUGFIX: axe y partagé (fixe) entre toutes les facettes FOD --
    # plus de scales = "free_y" -- pour pouvoir comparer directement
    # les niveaux de concentration entre clusters.
    facet_wrap(~ fod) +
    scale_y_log10(expand = expansion(mult = c(0.05, 0.45))) +  # marge du haut augmentée (était 0.3)
    labs(
      x = "Year", y = pigment_name,
      title = paste(pigment_name, "distribution within FOD clusters across years"),
      caption = paste0(
        "Tiret = moyenne géométrique ; barre = ± écart-type (log10) ; ",
        "lettres = groupes Šidák (p < 0.05) ; années en ordre chronologique ; "
      )
    ) +
    theme_classic() +
    theme(
      axis.title  = element_text(size = 12),
      axis.text   = element_text(size = 11),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text  = element_text(size = 12, face = "bold"),
      plot.title  = element_text(size = 14, face = "bold"),
      panel.spacing = unit(1.2, "lines")
    )
  
  if (save) {
    dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
    ggsave(
      filename = file.path(save_dir, paste0("violin_", pigment_name, ".png")),
      plot = p, width = 10, height = 6.5, dpi = 300
    )
    write.csv(
      stats_df %>% dplyr::mutate(pigment = pigment_name),
      file.path(save_dir, paste0("stats_", pigment_name, ".csv")),
      row.names = FALSE
    )
  }
  
  p
}


for (p_name in unique(pig_long$pigment)) {
  plot_pigment(p_name)
}

########################################################## plot interanual

n_date_all <- length(pig$date)
n_lon_all  <- length(pig$lon)
n_lat_all  <- length(pig$lat)

# Les arrays c_cond_XXX ont pour dimensions [date, lon, lat] : le premier
# indice (date) varie le plus vite dans l'ordre de stockage R (column-major).
# expand.grid() varie aussi son premier argument le plus vite -> on met
# date_idx en premier pour que l'ordre corresponde exactement à
# as.vector(array).
pig_long_all <- expand.grid(
  date_idx = seq_len(n_date_all),
  lon_idx  = seq_len(n_lon_all),
  lat_idx  = seq_len(n_lat_all)
) %>%
  dplyr::mutate(
    date   = pig$date[date_idx],
    year   = format(date, "%Y"),
    lon    = pig$lon[lon_idx],
    lat    = pig$lat[lat_idx],
    Chla   = as.vector(pig$c_cond_Chla),
    Perid  = as.vector(pig$c_cond_Per),
    But    = as.vector(pig$c_cond_But),
    Fuco   = as.vector(pig$c_cond_Fuco),
    Hex    = as.vector(pig$c_cond_Hex),
    Allo   = as.vector(pig$c_cond_Allo),
    Zeax   = as.vector(pig$c_cond_Zea),
    Chlb   = as.vector(pig$c_cond_Chlb),
    DvChla = as.vector(pig$c_cond_DvChla)
  ) %>%
  dplyr::select(-date_idx, -lon_idx, -lat_idx) %>%
  tidyr::pivot_longer(
    cols      = all_of(pigments_plot),
    names_to  = "pigment",
    values_to = "concentration"
  ) %>%
  dplyr::filter(concentration > 0)   # log10 impossible sur des valeurs <= 0

# ============================================================
# 7. Distribution interannuelle par pigment, TOUS FOD confondus
#    (1 seul plot, 1 facette = 1 pigment)
# ============================================================

compute_interannual_stats <- function(dat = pig_long_all,
                                      years_keep = c("2018", "2021", "2022", "2023")) {
  
  dat_filt <- dat %>%
    dplyr::filter(year %in% years_keep)
  
  # ---- Agrégation par DATE x PIGMENT d'abord (tous FOD confondus) ----
  daily_means <- dat_filt %>%
    dplyr::mutate(log_conc = log10(concentration)) %>%
    dplyr::group_by(pigment, date, year) %>%
    dplyr::summarise(log_conc_day = mean(log_conc, na.rm = TRUE), .groups = "drop")
  
  # ---- Résumé par année x pigment (basé sur les moyennes journalières) ----
  summary_stats <- daily_means %>%
    dplyr::group_by(pigment, year) %>%
    dplyr::summarise(
      n_days   = n(),
      log_mean = mean(log_conc_day, na.rm = TRUE),
      var_log  = var(log_conc_day,  na.rm = TRUE),
      .groups  = "drop"
    ) %>%
    dplyr::mutate(
      var_log = tidyr::replace_na(var_log, 0),
      sd_log  = sqrt(var_log),
      mean_c  = 10^log_mean,
      ymin    = 10^(log_mean - sd_log),
      ymax    = 10^(log_mean + sd_log)
    )
  
  # ---- Test sur les moyennes journalières, par pigment ----
  letters_by_pigment <- daily_means %>%
    dplyr::group_by(pigment) %>%
    dplyr::group_split() %>%
    purrr::map_dfr(function(dat_p) {
      this_pig <- as.character(unique(dat_p$pigment))
      
      if (dplyr::n_distinct(dat_p$year) < 2) {
        return(data.frame(pigment = this_pig, year = as.character(unique(dat_p$year)), letter = "a"))
      }
      
      model   <- lm(log_conc_day ~ year, data = dat_p)
      emm     <- emmeans::emmeans(model, ~ year)
      cld_out <- multcomp::cld(emm, Letters = letters, adjust = "sidak", sort = FALSE)
      
      data.frame(
        pigment = this_pig,
        year    = as.character(cld_out$year),
        letter  = trimws(cld_out$.group)
      )
    })
  
  summary_stats %>%
    dplyr::mutate(pigment = as.character(pigment), year = as.character(year)) %>%
    dplyr::left_join(letters_by_pigment, by = c("year", "pigment")) %>%
    dplyr::mutate(y = 10^(log10(ymax) + 1.5))  # même logique de décalage log10 que compute_pigment_stats
}


plot_interannual <- function(dat = pig_long_all, save_dir = out_dir, save = TRUE,
                             years_keep = c("2018", "2021", "2022", "2023")) {
  
  stats_df <- compute_interannual_stats(dat, years_keep)
  
  dat_all <- dat %>%
    dplyr::filter(year %in% years_keep) %>%
    dplyr::mutate(pigment = as.character(pigment), year = as.character(year)) %>%
    dplyr::left_join(stats_df %>% dplyr::select(year, pigment, mean_c), by = c("year", "pigment"))
  
  # ---- Ordre chronologique sur l'axe x ----
  dat_all  <- dat_all  %>% dplyr::mutate(year = factor(year, levels = years_keep))
  stats_df <- stats_df %>% dplyr::mutate(year = factor(year, levels = years_keep))
  
  p <- ggplot(dat_all, aes(x = year, y = concentration)) +
    geom_violin(trim = TRUE, scale = "width", alpha = 0.8, linewidth = 0.3, fill = "grey80") +
    geom_errorbar(
      data = stats_df,
      aes(x = year, y = mean_c, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE, width = 0.12, linewidth = 0.5, color = "black"
    ) +
    geom_crossbar(
      data = stats_df,
      aes(x = year, y = mean_c, ymin = mean_c, ymax = mean_c),
      inherit.aes = FALSE, width = 0.35, color = "black", linewidth = 0.4, fatten = 1
    ) +
    geom_text(
      data = stats_df,
      aes(x = year, y = y, label = letter),
      inherit.aes = FALSE, size = 3.5, fontface = "bold"
    ) +
    # BUGFIX: axe y partagé (fixe) entre toutes les facettes -- plus de
    # scales = "free_y" -- pour pouvoir comparer directement les niveaux
    # de concentration entre pigments.
    facet_wrap(~ pigment) +
    scale_y_log10(expand = expansion(mult = c(0.05, 0.45))) +
    labs(
      x = "Year", y = "Concentration",
      title = "Pigment distribution across years",
      caption = paste0(
        "Tiret = moyenne géométrique ; barre = ± écart-type (log10) ; ",
        "lettres = groupes Šidák (p < 0.05) ; années en ordre chronologique ; "
      )
    ) +
    theme_classic() +
    theme(
      axis.title  = element_text(size = 12),
      axis.text   = element_text(size = 11),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text  = element_text(size = 12, face = "bold"),
      plot.title  = element_text(size = 14, face = "bold"),
      panel.spacing = unit(1.2, "lines")
    )
  
  if (save) {
    dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
    ggsave(
      filename = file.path(save_dir, "violin_interannual_all_pigments.png"),
      plot = p, width = 12, height = 8, dpi = 300
    )
    write.csv(
      stats_df,
      file.path(save_dir, "stats_interannual_all_pigments.csv"),
      row.names = FALSE
    )
  }
  
  p
}

plot_interannual()