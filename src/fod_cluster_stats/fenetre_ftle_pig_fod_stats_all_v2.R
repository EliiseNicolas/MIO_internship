# ============================================================
# Distributions des RATIOS pigment/(somme des pigments hors Chla),
# du FTLE, de la Chla absolue, de la distribution des clusters FOD
# eux-memes, et du NASC -- au sein / a travers les clusters FOD,
# a partir du RDS unique (all_ds) deja aligne.
#
# NOUVEAUTES vs version precedente :
#   1) Distribution interannuelle des clusters FOD (couverture par
#      cluster, par annee)
#   2) Distribution de la Chla absolue (chla_total), interannuelle
#      ET intra-cluster FOD
#   3) Pour le NASC : version "facette = annee, x = FOD" pour
#      comparer les clusters entre eux au sein d'une meme annee
#   4) Pour les distributions INTERANNUELLES (FTLE, FOD, NASC, Chla,
#      ratios) : les lettres de significativite sont maintenant
#      positionnees juste sous le titre de chaque sous-graphe
#      (position fixe y = Inf, independante de l'echelle des
#      donnees), comme c'est deja le cas pour les figures par
#      cluster FOD (pigments).
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(emmeans)
library(multcomp)
library(multcompView)
library(purrr)
library(ggh4x)

# ------------------------------------------------------------
# Chemins
# ------------------------------------------------------------

path_all_ds <- "F:/data_elise/prediction_ds/ds_ftle_pig_fod_ALL_DATES.rds"
out_dir <- "C:/Users/mmolinet/elisou_ta_stagiaire_pref/MIO_internship_III/figures/distribution_per_fod_cluster/violin_plots_ratio_ftle"

years_keep <- c("2018", "2021", "2022", "2023")

# ------------------------------------------------------------
# Palette de couleurs FOD (identique au script pigments)
# ------------------------------------------------------------

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
  fod_chr[is.na(fod_chr) | fod_chr == "NA"] <- "NA"
  lbl <- fod_label_map[fod_chr]
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

# ============================================================
# 1. Chargement du dataset unique deja aligne
# ============================================================

all_ds <- readRDS(path_all_ds)
str(all_ds)
dates <- all_ds$date
lons  <- all_ds$lon
lats  <- all_ds$lat

n_date <- length(dates)
n_lon  <- length(lons)
n_lat  <- length(lats)

# fod est stocke en [lon, lat, date] -> on le remet en [date, lon, lat]
fod_arr <- aperm(all_ds$fod, c(3, 1, 2))

# ------------------------------------------------------------
# Variables ratio pigment/(somme des pigments hors Chla)
# (suffixe "_totpig") a etudier
# ------------------------------------------------------------
ratio_vars <- grep("_totpig$", names(all_ds$pig), value = TRUE)

# Variables etudiees, conformement a la structure du dataset :
# - les 8 ratios pigment/(somme des pigments hors Chla), suffixe
#   "_totpig" : per_totpig, but_totpig, fuco_totpig, hex_totpig,
#   allo_totpig, zea_totpig, chlb_totpig, dvchla_totpig
# - la Chla absolue : chla_total
# Les concentrations brutes des autres pigments (Chla, Per, But,
# Fuco, Hex, Allo, Zea, Chlb, DvChla) ne sont PAS etudiees ici.

cat("Ratios trouves :", paste(ratio_vars, collapse = ", "), "\n")

# ------------------------------------------------------------
# NOUVEAU (point 2) : variable Chla absolue
# ATTENTION : nom suppose "chla_total" d'apres les commentaires du
# script d'origine -- verifie ce nom dans names(all_ds$pig) et
# adapte si besoin.
# ------------------------------------------------------------
chla_var <- "chla_total"
if (!chla_var %in% names(all_ds$pig)) {
  warning(paste0("'", chla_var, "' introuvable dans all_ds$pig -- verifie le nom exact ",
                 "du champ de Chla absolue (names(all_ds$pig)) et corrige la variable chla_var."))
}

# ============================================================
# 2. Table maitre au format long (1 ligne = 1 pixel x 1 date)
# ============================================================

master_long <- expand.grid(
  date_idx = seq_len(n_date),
  lon_idx  = seq_len(n_lon),
  lat_idx  = seq_len(n_lat)
) %>%
  dplyr::mutate(
    date = dates[date_idx],
    year = format(date, "%Y"),
    lon  = lons[lon_idx],
    lat  = lats[lat_idx],
    fod  = factor(as.vector(fod_arr)),
    ftle = as.vector(all_ds$ftle)
  )

# extraction des ratios ET de la Chla absolue
vars_to_extract <- unique(c(ratio_vars, chla_var))
for (v in vars_to_extract) {
  if (v %in% names(all_ds$pig)) {
    master_long[[v]] <- as.vector(all_ds$pig[[v]])
  }
}

master_long <- master_long %>% dplyr::select(-date_idx, -lon_idx, -lat_idx)

# ------------------------------------------------------------
# Format long des ratios
# ------------------------------------------------------------

ratio_long <- master_long %>%
  dplyr::select(date, year, lon, lat, fod, dplyr::all_of(ratio_vars)) %>%
  tidyr::pivot_longer(
    cols      = dplyr::all_of(ratio_vars),
    names_to  = "ratio_name",
    values_to = "ratio"
  ) %>%
  dplyr::filter(is.finite(ratio), ratio >= 0)

# ------------------------------------------------------------
# Table FTLE seule
# ------------------------------------------------------------

ftle_long <- master_long %>%
  dplyr::select(date, year, lon, lat, fod, ftle) %>%
  dplyr::filter(is.finite(ftle)) %>%
  dplyr::mutate(variable = "FTLE")

# ------------------------------------------------------------
# NOUVEAU (point 2) : table Chla absolue seule
# ------------------------------------------------------------

chla_long <- master_long %>%
  dplyr::select(date, year, lon, lat, fod, dplyr::all_of(chla_var)) %>%
  dplyr::rename(chla = dplyr::all_of(chla_var)) %>%
  dplyr::filter(is.finite(chla), chla > 0) %>%
  dplyr::mutate(variable = "Chla")

# ------------------------------------------------------------
# NOUVEAU (point 1) : table brute pour la distribution des FOD
# eux-memes (independante de toute variable continue)
# ------------------------------------------------------------

fod_distrib_base <- master_long %>%
  dplyr::filter(year %in% years_keep, !is.na(fod), as.character(fod) != "NA") %>%
  dplyr::mutate(fod = as.character(fod))

# ============================================================
# 3. Fonctions generiques -- par cluster FOD x annee
#    (facette = FOD, x = annee)  -- INCHANGE
# ============================================================

compute_stats_by_fod <- function(dat, value_col, entity_col, entity_val,
                                 years_keep = c("2018","2021","2022","2023"),
                                 log_transform = TRUE,
                                 y_limits = NULL) {
  
  dat_sub <- dat %>%
    dplyr::filter(.data[[entity_col]] == entity_val, year %in% years_keep) %>%
    dplyr::mutate(val = if (log_transform) log10(.data[[value_col]]) else .data[[value_col]])
  
  daily_means <- dat_sub %>%
    dplyr::group_by(date, year, fod) %>%
    dplyr::summarise(val_day = mean(val, na.rm = TRUE), .groups = "drop")
  
  summary_stats <- daily_means %>%
    dplyr::group_by(year, fod) %>%
    dplyr::summarise(
      n_days = dplyr::n(),
      m      = mean(val_day, na.rm = TRUE),
      var_v  = var(val_day,  na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      var_v  = tidyr::replace_na(var_v, 0),
      sd_v   = sqrt(var_v),
      mean_c = if (log_transform) 10^m else m,
      ymin   = if (log_transform) 10^(m - sd_v) else m - sd_v,
      ymax   = if (log_transform) 10^(m + sd_v) else m + sd_v
    )
  
  letters_by_fod <- daily_means %>%
    dplyr::group_by(fod) %>%
    dplyr::group_split() %>%
    purrr::map_dfr(function(dat_fod) {
      this_fod <- as.character(unique(dat_fod$fod))
      if (dplyr::n_distinct(dat_fod$year) < 2) {
        return(data.frame(fod = this_fod, year = as.character(unique(dat_fod$year)), letter = "a"))
      }
      model   <- lm(val_day ~ year, data = dat_fod)
      emm     <- emmeans::emmeans(model, ~ year)
      cld_out <- multcomp::cld(emm, Letters = letters, adjust = "sidak", sort = FALSE)
      data.frame(fod = this_fod, year = as.character(cld_out$year), letter = trimws(cld_out$.group))
    })
  
  out <- summary_stats %>%
    dplyr::mutate(fod = as.character(fod), year = as.character(year)) %>%
    dplyr::left_join(letters_by_fod, by = c("year", "fod"))
  
  out
}

plot_by_fod <- function(dat, value_col, entity_col, entity_val, label,
                        save_dir = out_dir, save = TRUE,
                        years_keep = c("2018","2021","2022","2023"),
                        log_transform = TRUE,
                        y_limits = NULL) {
  
  stats_df <- compute_stats_by_fod(dat, value_col, entity_col, entity_val, years_keep, log_transform, y_limits)
  
  dat_sub <- dat %>%
    dplyr::filter(.data[[entity_col]] == entity_val, year %in% years_keep) %>%
    dplyr::mutate(fod = as.character(fod), year = as.character(year)) %>%
    dplyr::left_join(stats_df %>% dplyr::select(year, fod, mean_c), by = c("year", "fod"))
  
  dat_sub  <- dat_sub  %>% dplyr::mutate(year = factor(year, levels = years_keep))
  stats_df <- stats_df %>% dplyr::mutate(year = factor(year, levels = years_keep))
  
  dat_sub  <- dat_sub  %>% dplyr::mutate(fod = relabel_fod(fod))
  stats_df <- stats_df %>% dplyr::mutate(fod = relabel_fod(fod))
  
  dat_sub  <- dat_sub  %>% dplyr::mutate(fod = droplevels(fod))
  stats_df <- stats_df %>% dplyr::mutate(fod = droplevels(fod))
  
  strip_colors <- fod_cols_named[levels(dat_sub$fod)]
  strip_colors[is.na(strip_colors) | levels(dat_sub$fod) == "NA"] <- "grey50"
  
  p <- ggplot(dat_sub, aes(x = year, y = .data[[value_col]])) +
    geom_violin(trim = TRUE, scale = "width", alpha = 0.8,
                linewidth = 0.3, fill = "grey75") +
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
      aes(x = year, y = Inf, label = letter),
      inherit.aes = FALSE, vjust = 1.4, size = 3.5, fontface = "bold"
    ) +
    ggh4x::facet_wrap2(
      ~ fod, scales = "free_y",
      strip = ggh4x::strip_themed(
        background_x = ggh4x::elem_list_rect(fill = strip_colors),
        text_x       = ggh4x::elem_list_text(colour = "white", face = "bold")
      )
    ) +
    labs(
      x = "Year", y = label,
      title = paste(label, "distribution within FOD clusters across years"),
      caption = "Tiret = moyenne ; barre = +/- ecart-type ; lettres = groupes Sidak (p < 0.05)"
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

# ============================================================
# 4. Fonctions generiques -- distribution interannuelle, tous FOD
#    confondus (facette = variable/entite, x = annee)
#
#    MODIF (point 4) : les lettres de significativite sont
#    positionnees a y = Inf (juste sous le titre du sous-graphe),
#    au lieu d'une position calculee a partir de ymax -- ce qui
#    les rend visibles et bien placees quelle que soit l'echelle
#    des donnees, comme pour les figures par cluster FOD.
# ============================================================

compute_interannual_stats_generic <- function(dat, value_col, entity_col,
                                              years_keep = c("2018","2021","2022","2023"),
                                              log_transform = TRUE,
                                              y_limits = NULL) {
  
  dat_sub <- dat %>%
    dplyr::filter(year %in% years_keep) %>%
    dplyr::mutate(val = if (log_transform) log10(.data[[value_col]]) else .data[[value_col]])
  
  daily_means <- dat_sub %>%
    dplyr::group_by(.data[[entity_col]], date, year) %>%
    dplyr::summarise(val_day = mean(val, na.rm = TRUE), .groups = "drop")
  
  summary_stats <- daily_means %>%
    dplyr::group_by(.data[[entity_col]], year) %>%
    dplyr::summarise(
      n_days = dplyr::n(),
      m      = mean(val_day, na.rm = TRUE),
      var_v  = var(val_day,  na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      var_v  = tidyr::replace_na(var_v, 0),
      sd_v   = sqrt(var_v),
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
  
  summary_stats %>%
    dplyr::mutate(!!entity_col := as.character(.data[[entity_col]]), year = as.character(year)) %>%
    dplyr::left_join(letters_by_entity, by = c("year", entity_col))
}

plot_interannual_generic <- function(dat, value_col, entity_col, label,
                                     save_dir = out_dir, save = TRUE,
                                     years_keep = c("2018","2021","2022","2023"),
                                     log_transform = TRUE,
                                     y_limits = NULL) {
  
  stats_df <- compute_interannual_stats_generic(dat, value_col, entity_col, years_keep, log_transform, y_limits)
  
  dat_all <- dat %>%
    dplyr::filter(year %in% years_keep) %>%
    dplyr::mutate(!!entity_col := as.character(.data[[entity_col]]), year = as.character(year)) %>%
    dplyr::left_join(stats_df %>% dplyr::select(year, dplyr::all_of(entity_col), mean_c), by = c("year", entity_col))
  
  dat_all  <- dat_all  %>% dplyr::mutate(year = factor(year, levels = years_keep))
  stats_df <- stats_df %>% dplyr::mutate(year = factor(year, levels = years_keep))
  
  p <- ggplot(dat_all, aes(x = year, y = .data[[value_col]])) +
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
      aes(x = year, y = Inf, label = letter),
      inherit.aes = FALSE, vjust = 1.4, size = 3.5, fontface = "bold"
    ) +
    facet_wrap(as.formula(paste("~", entity_col))) +
    labs(
      x = "Year", y = label,
      title = paste(label, "distribution across years"),
      caption = if (log_transform) {
        "Tiret = moyenne geometrique ; barre = +/- ecart-type (log10) ; lettres = groupes Sidak (p < 0.05)"
      } else {
        "Tiret = moyenne ; barre = +/- ecart-type ; lettres = groupes Sidak (p < 0.05)"
      }
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
  
  if (log_transform) {
    p <- p + scale_y_log10(expand = expansion(mult = c(0.05, 0.45)))
  } else if (!is.null(y_limits)) {
    p <- p + coord_cartesian(ylim = y_limits, clip = "off")
  } else {
    p <- p + scale_y_continuous(expand = expansion(mult = c(0.05, 0.45)))
  }
  
  if (save) {
    dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
    fname <- gsub("[^A-Za-z0-9_]", "_", label)
    ggsave(file.path(save_dir, paste0("violin_interannual_", fname, ".png")), plot = p, width = 12, height = 8, dpi = 300)
    write.csv(stats_df, file.path(save_dir, paste0("stats_interannual_", fname, ".csv")), row.names = FALSE)
  }
  
  p
}

# ============================================================
# 4bis. NOUVEAU (point 3) -- comparaison des clusters FOD AU SEIN
#       de chaque annee : facette = annee, x = FOD.
#       Utilise pour le NASC (voir section 8), reutilisable pour
#       d'autres variables (ftle, chla, ratios...).
# ============================================================

compute_stats_by_year <- function(dat, value_col, entity_col, entity_val,
                                  years_keep = c("2018","2021","2022","2023"),
                                  log_transform = TRUE) {
  
  dat_sub <- dat %>%
    dplyr::filter(.data[[entity_col]] == entity_val, year %in% years_keep,
                  !is.na(fod), as.character(fod) != "NA") %>%
    dplyr::mutate(fod = as.character(fod),
                  val = if (log_transform) log10(.data[[value_col]]) else .data[[value_col]])
  
  daily_means <- dat_sub %>%
    dplyr::group_by(date, year, fod) %>%
    dplyr::summarise(val_day = mean(val, na.rm = TRUE), .groups = "drop")
  
  summary_stats <- daily_means %>%
    dplyr::group_by(year, fod) %>%
    dplyr::summarise(
      n_days = dplyr::n(), m = mean(val_day, na.rm = TRUE), var_v = var(val_day, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      var_v  = tidyr::replace_na(var_v, 0),
      sd_v   = sqrt(var_v),
      mean_c = if (log_transform) 10^m else m,
      ymin   = if (log_transform) 10^(m - sd_v) else m - sd_v,
      ymax   = if (log_transform) 10^(m + sd_v) else m + sd_v
    )
  
  # lettres = comparaison des clusters FOD ENTRE EUX, separement pour
  # chaque annee (donc l'interpretation est "au sein de cette annee")
  letters_by_year <- daily_means %>%
    dplyr::group_by(year) %>%
    dplyr::group_split() %>%
    purrr::map_dfr(function(dat_y) {
      this_year <- unique(dat_y$year)
      if (dplyr::n_distinct(dat_y$fod) < 2) {
        return(data.frame(year = this_year, fod = unique(dat_y$fod), letter = "a"))
      }
      model   <- lm(val_day ~ fod, data = dat_y)
      emm     <- emmeans::emmeans(model, ~ fod)
      cld_out <- multcomp::cld(emm, Letters = letters, adjust = "sidak", sort = FALSE)
      data.frame(year = this_year, fod = as.character(cld_out$fod), letter = trimws(cld_out$.group))
    })
  
  summary_stats %>%
    dplyr::mutate(fod = as.character(fod), year = as.character(year)) %>%
    dplyr::left_join(letters_by_year, by = c("year", "fod"))
}

plot_by_year <- function(dat, value_col, entity_col, entity_val, label,
                         save_dir = out_dir, save = TRUE,
                         years_keep = c("2018","2021","2022","2023"),
                         log_transform = TRUE) {
  
  stats_df <- compute_stats_by_year(dat, value_col, entity_col, entity_val, years_keep, log_transform)
  
  dat_sub <- dat %>%
    dplyr::filter(.data[[entity_col]] == entity_val, year %in% years_keep,
                  !is.na(fod), as.character(fod) != "NA") %>%
    dplyr::mutate(fod = as.character(fod), year = as.character(year)) %>%
    dplyr::left_join(stats_df %>% dplyr::select(year, fod, mean_c), by = c("year", "fod"))
  
  dat_sub  <- dat_sub  %>% dplyr::mutate(fod = relabel_fod(fod)) %>% dplyr::mutate(fod = droplevels(fod))
  stats_df <- stats_df %>% dplyr::mutate(fod = relabel_fod(fod)) %>% dplyr::mutate(fod = droplevels(fod))
  
  dat_sub  <- dat_sub  %>% dplyr::mutate(year = factor(year, levels = years_keep))
  stats_df <- stats_df %>% dplyr::mutate(year = factor(year, levels = years_keep))
  
  p <- ggplot(dat_sub, aes(x = fod, y = .data[[value_col]], fill = fod)) +
    geom_violin(trim = TRUE, scale = "width", alpha = 0.8, linewidth = 0.3) +
    scale_fill_manual(values = fod_cols_named, guide = "none") +
    geom_errorbar(
      data = stats_df, aes(x = fod, y = mean_c, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE, width = 0.12, linewidth = 0.5, color = "black"
    ) +
    geom_crossbar(
      data = stats_df, aes(x = fod, y = mean_c, ymin = mean_c, ymax = mean_c),
      inherit.aes = FALSE, width = 0.35, color = "black", linewidth = 0.2
    ) +
    geom_text(
      data = stats_df, aes(x = fod, y = Inf, label = letter),
      inherit.aes = FALSE, vjust = 1.4, size = 3.5, fontface = "bold"
    ) +
    facet_wrap(~ year, scales = "free_y") +
    labs(
      x = "FOD", y = label,
      title = paste(label, "- comparaison des clusters FOD, par annee"),
      caption = "Tiret = moyenne ; barre = +/- ecart-type ; lettres = groupes Sidak (p < 0.05), comparaison inter-cluster au sein de chaque annee"
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
  
  if (log_transform) {
    p <- p + scale_y_log10(expand = expansion(mult = c(0.05, 0.45)))
  } else {
    p <- p + scale_y_continuous(expand = expansion(mult = c(0.05, 0.45)))
  }
  
  if (save) {
    dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
    fname <- gsub("[^A-Za-z0-9_]", "_", label)
    ggsave(file.path(save_dir, paste0("violin_by_year_", fname, ".png")), plot = p, width = 11, height = 7, dpi = 300)
    write.csv(stats_df, file.path(save_dir, paste0("stats_by_year_", fname, ".csv")), row.names = FALSE)
  }
  
  p
}

# ============================================================
# 5. Diagrammes du nombre de donnees effectives (jours)
#    -- INCHANGE
# ============================================================

compute_n_days_by_fod <- function(dat, value_col, entity_col, entity_val,
                                  years_keep = c("2018","2021","2022","2023")) {
  dat %>%
    dplyr::filter(.data[[entity_col]] == entity_val, year %in% years_keep,
                  is.finite(.data[[value_col]])) %>%
    dplyr::distinct(date, year, fod) %>%
    dplyr::group_by(year, fod) %>%
    dplyr::summarise(n_days = dplyr::n(), .groups = "drop")
}

plot_n_days_by_fod <- function(dat, value_col, entity_col, entity_val, label,
                               save_dir = out_dir, save = TRUE,
                               years_keep = c("2018","2021","2022","2023")) {
  
  n_days_df <- compute_n_days_by_fod(dat, value_col, entity_col, entity_val, years_keep) %>%
    dplyr::mutate(year = factor(year, levels = years_keep), fod = relabel_fod(fod))
  
  p <- ggplot(n_days_df, aes(x = year, y = n_days, fill = fod)) +
    geom_col(position = position_dodge2(preserve = "single"), color = "black", linewidth = 0.2) +
    scale_fill_manual(values = fod_cols_named, name = "FOD") +
    labs(
      x = "Year", y = "Nombre de jours effectifs",
      title = paste(label, "- jours utilises pour la moyenne, par cluster FOD")
    ) +
    theme_classic() +
    theme(
      axis.title  = element_text(size = 12),
      axis.text   = element_text(size = 11),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title  = element_text(size = 13, face = "bold")
    )
  
  if (save) {
    dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
    ggsave(file.path(save_dir, paste0("ndays_", label, ".png")), plot = p, width = 9, height = 5.5, dpi = 300)
    write.csv(n_days_df, file.path(save_dir, paste0("ndays_", label, ".csv")), row.names = FALSE)
  }
  
  p
}

compute_n_days_interannual <- function(dat, value_col, entity_col,
                                       years_keep = c("2018","2021","2022","2023")) {
  dat %>%
    dplyr::filter(year %in% years_keep, is.finite(.data[[value_col]])) %>%
    dplyr::distinct(.data[[entity_col]], date, year) %>%
    dplyr::group_by(.data[[entity_col]], year) %>%
    dplyr::summarise(n_days = dplyr::n(), .groups = "drop")
}

plot_n_days_interannual <- function(dat, value_col, entity_col, label,
                                    save_dir = out_dir, save = TRUE,
                                    years_keep = c("2018","2021","2022","2023")) {
  
  n_days_df <- compute_n_days_interannual(dat, value_col, entity_col, years_keep) %>%
    dplyr::mutate(year = factor(year, levels = years_keep))
  
  p <- ggplot(n_days_df, aes(x = year, y = n_days)) +
    geom_col(fill = "grey60", color = "black", linewidth = 0.2) +
    facet_wrap(as.formula(paste("~", entity_col))) +
    labs(
      x = "Year", y = "Nombre de jours effectifs",
      title = paste(label, "- jours utilises pour la moyenne interannuelle (tous FOD confondus)")
    ) +
    theme_classic() +
    theme(
      axis.title  = element_text(size = 12),
      axis.text   = element_text(size = 11),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text  = element_text(size = 12, face = "bold"),
      plot.title  = element_text(size = 13, face = "bold"),
      panel.spacing = unit(1.2, "lines")
    )
  
  if (save) {
    dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
    fname <- gsub("[^A-Za-z0-9_]", "_", label)
    ggsave(file.path(save_dir, paste0("ndays_interannual_", fname, ".png")), plot = p, width = 11, height = 7, dpi = 300)
    write.csv(n_days_df, file.path(save_dir, paste0("ndays_interannual_", fname, ".csv")), row.names = FALSE)
  }
  
  p
}

# ============================================================
# 5bis. NOUVEAU (point 1) -- distribution interannuelle des
#       clusters FOD eux-memes (couverture en % de jours,
#       par cluster, comparee entre annees).
#       Lettres = test de proportions par paires (FDR), comme
#       pour un test de Sidak mais adapte a des donnees de
#       comptage/proportion plutot qu'a une moyenne continue.
# ============================================================

compute_fod_distrib_stats <- function(dat, years_keep = c("2018","2021","2022","2023")) {
  
  counts <- dat %>%
    dplyr::group_by(fod, year) %>%
    dplyr::summarise(n_days = dplyr::n_distinct(date), .groups = "drop")
  
  totals_per_year <- counts %>%
    dplyr::group_by(year) %>%
    dplyr::summarise(tot = sum(n_days), .groups = "drop")
  
  counts <- counts %>%
    dplyr::left_join(totals_per_year, by = "year") %>%
    dplyr::mutate(pct = 100 * n_days / tot)
  
  letters_by_fod <- counts %>%
    dplyr::group_by(fod) %>%
    dplyr::group_split() %>%
    purrr::map_dfr(function(d) {
      this_fod <- unique(d$fod)
      if (nrow(d) < 2) {
        return(data.frame(fod = this_fod, year = d$year, letter = "a"))
      }
      pw   <- pairwise.prop.test(d$n_days, d$tot, p.adjust.method = "fdr")
      pmat <- pw$p.value
      pv <- c(); nm <- c()
      for (i in seq_len(nrow(pmat))) for (j in seq_len(ncol(pmat))) {
        if (!is.na(pmat[i, j])) {
          nm <- c(nm, paste(rownames(pmat)[i], colnames(pmat)[j], sep = "-"))
          pv <- c(pv, pmat[i, j])
        }
      }
      names(pv) <- nm
      lt <- multcompView::multcompLetters(pv)$Letters
      data.frame(fod = this_fod, year = names(lt), letter = as.character(lt))
    })
  
  counts %>% dplyr::left_join(letters_by_fod, by = c("fod", "year"))
}

plot_fod_distrib <- function(dat, years_keep = c("2018","2021","2022","2023"),
                             save_dir = out_dir, save = TRUE) {
  
  stats_df <- compute_fod_distrib_stats(dat, years_keep) %>%
    dplyr::mutate(year = factor(year, levels = years_keep), fod = relabel_fod(fod)) %>%
    dplyr::mutate(fod = droplevels(fod))
  
  p <- ggplot(stats_df, aes(x = year, y = pct)) +
    geom_col(fill = "grey65", color = "black", linewidth = 0.2, width = 0.7) +
    geom_text(aes(y = Inf, label = letter), vjust = 1.4, size = 3.5, fontface = "bold") +
    facet_wrap(~ fod, scales = "free_y") +
    labs(
      x = "Year", y = "% des jours (par cluster FOD)",
      title = "Distribution interannuelle des clusters FOD",
      caption = "Lettres = groupes de proportions significativement differentes (test par paires, FDR, p < 0.05)"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.45))) +
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
    ggsave(file.path(save_dir, "distribution_FOD_interannuelle.png"), plot = p, width = 11, height = 7, dpi = 300)
    write.csv(stats_df, file.path(save_dir, "stats_distribution_FOD_interannuelle.csv"), row.names = FALSE)
  }
  
  p
}

# ============================================================
# 6. Execution -- ratios pigment/(somme des pigments hors Chla)
#    Echelle LINEAIRE, ylim adapte a chaque ratio
# ============================================================

# for (rv in ratio_vars) {
#   rv_vals  <- ratio_long$ratio[ratio_long$ratio_name == rv & ratio_long$year %in% years_keep]
#   stats_rv <- compute_stats_by_fod(ratio_long, value_col = "ratio", entity_col = "ratio_name",
#                                    entity_val = rv, years_keep = years_keep, log_transform = FALSE)
#   
#   rv_upper <- max(
#     stats::quantile(rv_vals, probs = 0.999, na.rm = TRUE),
#     max(stats_rv$ymax, na.rm = TRUE)
#   )
#   y_lim_rv <- c(0, rv_upper * 1.25)
#   
#   plot_by_fod(ratio_long, value_col = "ratio", entity_col = "ratio_name",
#               entity_val = rv, label = rv, years_keep = years_keep,
#               log_transform = FALSE, y_limits = y_lim_rv)
# }

# 6.2 Distribution interannuelle, tous FOD confondus, un plot facette par ratio
plot_interannual_generic(ratio_long, value_col = "ratio", entity_col = "ratio_name",
                         label = "ratios_pigment_totpig", years_keep = years_keep,
                         log_transform = FALSE)
plot_n_days_interannual(ratio_long, value_col = "ratio", entity_col = "ratio_name",
                        label = "ratios_pigment_totpig", years_keep = years_keep)

# ============================================================
# 7. Execution -- FTLE (echelle lineaire, axe libre par facette)
# ============================================================

plot_by_fod(ftle_long, value_col = "ftle", entity_col = "variable",
            entity_val = "FTLE", label = "FTLE", years_keep = years_keep, log_transform = FALSE)
plot_n_days_by_fod(ftle_long, value_col = "ftle", entity_col = "variable",
                   entity_val = "FTLE", label = "FTLE", years_keep = years_keep)

plot_interannual_generic(ftle_long, value_col = "ftle", entity_col = "variable",
                         label = "FTLE", years_keep = years_keep, log_transform = FALSE)
plot_n_days_interannual(ftle_long, value_col = "ftle", entity_col = "variable",
                        label = "FTLE", years_keep = years_keep)

# ============================================================
# 7bis. NOUVEAU (point 2) -- Execution Chla absolue (chla_total)
#       Echelle log (usuelle pour de la concentration en pigment)
# ============================================================

plot_by_fod(chla_long, value_col = "chla", entity_col = "variable",
            entity_val = "Chla", label = "Chla_a", years_keep = years_keep,
            log_transform = TRUE)
plot_n_days_by_fod(chla_long, value_col = "chla", entity_col = "variable",
                   entity_val = "Chla", label = "Chla_a", years_keep = years_keep)

plot_interannual_generic(chla_long, value_col = "chla", entity_col = "variable",
                         label = "Chla_a", years_keep = years_keep, log_transform = TRUE)
plot_n_days_interannual(chla_long, value_col = "chla", entity_col = "variable",
                        label = "Chla_a", years_keep = years_keep)

# ============================================================
# 7ter. NOUVEAU (point 1) -- Execution distribution interannuelle
#       des clusters FOD
# ============================================================

plot_fod_distrib(fod_distrib_base, years_keep = years_keep)

# ============================================================
# 8. Execution -- NASC
# ============================================================

freq <- 120
nasc_ds <- readRDS(paste0(
  "F:/data_elise/ds_NASC_pig_ftle_fod/ds_NASC_per_esu_all/new_ratio_ds_NASC_per_esu_pig_ftle_fod_2018_2021_2022_2023_transect_",
  freq, "kHz_mask9.rds"
))

nasc_long <- nasc_ds %>%
  dplyr::mutate(
    date = as.Date(time_nasc),
    year = format(time_nasc, "%Y"),
    fod = trimws(as.character(fod)),
    fod = ifelse(is.na(fod) | fod == "NA", NA_character_, fod),
    variable = "NASC"
  ) %>%
  dplyr::filter(
    is.finite(nasc),
    nasc > 0,
    year %in% years_keep
  )

cat("NASC valides :", nrow(nasc_long), "/", nrow(nasc_ds), "\n")
print(unique(nasc_long$fod))

plot_by_fod(nasc_long, value_col = "nasc", entity_col = "variable",
            entity_val = "NASC", label = paste0("NASC (", freq, " kHz)"), years_keep = years_keep,
            log_transform = TRUE)
plot_n_days_by_fod(nasc_long, value_col = "nasc", entity_col = "variable",
                   entity_val = "NASC", label = paste0("NASC (", freq, " kHz)"), years_keep = years_keep)

plot_interannual_generic(nasc_long, value_col = "nasc", entity_col = "variable",
                         label = paste0("NASC (", freq, " kHz)"), years_keep = years_keep, log_transform = TRUE)
plot_n_days_interannual(nasc_long, value_col = "nasc", entity_col = "variable",
                        label = paste0("NASC (", freq, " kHz)"), years_keep = years_keep)

# NOUVEAU (point 3) : comparaison des clusters FOD entre eux, au
# sein de chaque annee -- pour voir si certains clusters ont plus
# de NASC que d'autres
plot_by_year(nasc_long, value_col = "nasc", entity_col = "variable",
             entity_val = "NASC", label = paste0("NASC (", freq, " kHz)"),
             years_keep = years_keep, log_transform = TRUE)

sort(unique(nasc_ds$fod))