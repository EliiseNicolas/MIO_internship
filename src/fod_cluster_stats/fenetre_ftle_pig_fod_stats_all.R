# ============================================================
# Distributions des RATIOS pigment/Chla et du FTLE au sein des
# clusters FOD, a partir du RDS unique (all_ds) deja aligne
# (memes dates, meme grille pour ftle / pig / fod).
#
# Reprend exactement la meme logique statistique et graphique que
# le script pigments d'origine (violon + moyenne +/- variance +
# lettres Sidak), generalisee pour s'appliquer :
#  - aux ratios pigment/Chla (colonnes *_chla dans all_ds$pig),
#    par annee x cluster FOD, puis en distribution interannuelle
#    (echelle LINEAIRE, ylim adapte par ratio -- cf section 6)
#  - au FTLE, par annee x cluster FOD, puis en distribution
#    interannuelle (echelle lineaire, axe libre par facette)
#
# Etiquettes de facette (strip) colorees par FOD plutot que le
# violon lui-meme ; NA en gris ; lettres des tests Sidak calees
# pres du haut visible du panneau (94% de la hauteur), sans jamais
# depasser.
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

# ============================================================
# 1. Chargement du dataset unique deja aligne
# ============================================================

all_ds <- readRDS(path_all_ds)

dates <- all_ds$date
lons  <- all_ds$lon
lats  <- all_ds$lat

n_date <- length(dates)
n_lon  <- length(lons)
n_lat  <- length(lats)

# fod est stocke en [lon, lat, date] -> on le remet en [date, lon, lat]
# pour etre coherent avec ftle et les arrays de pig (deja [date, lon, lat])
fod_arr <- aperm(all_ds$fod, c(3, 1, 2))

# Variables ratio pigment/Chla (suffixe _chla) a etudier
ratio_vars <- grep("_chla$", names(all_ds$pig), value = TRUE)

# chla_chla vaut toujours 1 (Chla / Chla) : sans interet, exclu
ratio_vars <- setdiff(ratio_vars, "chla_chla")

cat("Ratios trouves :", paste(ratio_vars, collapse = ", "), "\n")

# ============================================================
# 2. Table maitre au format long (1 ligne = 1 pixel x 1 date)
# ============================================================
# date_idx en premier : c'est la 1ere dimension des arrays (varie le
# plus vite en column-major R), pour matcher exactement as.vector().

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

for (v in ratio_vars) {
  master_long[[v]] <- as.vector(all_ds$pig[[v]])
}

master_long <- master_long %>% dplyr::select(-date_idx, -lon_idx, -lat_idx)

# ------------------------------------------------------------
# Format long des ratios (equivalent de pig_long)
# ------------------------------------------------------------
# Plus de contrainte "> 0" : on ne passe plus par log10, donc 0 est
# une valeur valide pour un ratio. On garde uniquement les valeurs
# finies (exclut Inf/NaN issus d'une division par 0/NA).

ratio_long <- master_long %>%
  dplyr::select(date, year, lon, lat, fod, dplyr::all_of(ratio_vars)) %>%
  tidyr::pivot_longer(
    cols      = dplyr::all_of(ratio_vars),
    names_to  = "ratio_name",
    values_to = "ratio"
  ) %>%
  dplyr::filter(is.finite(ratio), ratio >= 0)

# ------------------------------------------------------------
# Table FTLE seule (une seule variable, pas de pivot necessaire)
# ------------------------------------------------------------

ftle_long <- master_long %>%
  dplyr::select(date, year, lon, lat, fod, ftle) %>%
  dplyr::filter(is.finite(ftle)) %>%
  dplyr::mutate(variable = "FTLE")   # colonne factice pour reutiliser les memes fonctions generiques

# ============================================================
# 3. Fonctions generiques -- par cluster FOD x annee
# ============================================================

compute_stats_by_fod <- function(dat, value_col, entity_col, entity_val,
                                 years_keep = c("2018","2021","2022","2023"),
                                 log_transform = TRUE,
                                 y_limits = NULL) {
  
  dat_sub <- dat %>%
    dplyr::filter(.data[[entity_col]] == entity_val, year %in% years_keep) %>%
    dplyr::mutate(val = if (log_transform) log10(.data[[value_col]]) else .data[[value_col]])
  
  # ---- Moyenne par DATE x FOD d'abord (reduit la pseudo-replication) ----
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
      var_v  = tidyr::replace_na(var_v, 0),   # var() vaut NA si n_days == 1
      sd_v   = sqrt(var_v),                   # ecart-type, coherent avec la legende ("+/- ecart-type")
      mean_c = if (log_transform) 10^m else m,
      ymin   = if (log_transform) 10^(m - sd_v) else m - sd_v,
      ymax   = if (log_transform) 10^(m + sd_v) else m + sd_v
    )
  
  # ---- Test sur les moyennes journalieres ----
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
  
  # Position des lettres : le plus haut possible sans jamais depasser.
  #  - log_transform : +1.5 dex au-dessus du haut de la barre d'erreur
  #    (multiplicatif, scale-invariant).
  #  - lineaire : calee a 94% de la hauteur visible du panneau.
  #    * y_limits fixe (coord_cartesian, ex. ratios) -> le haut visible
  #      vaut exactement y_limits[2].
  #    * echelle libre (ex. FTLE, scales="free_y") -> le haut visible
  #      depend du max reel des donnees de CETTE facette.
  if (log_transform) {
    out <- out %>% dplyr::mutate(y = 10^(log10(ymax) + 1.5))
  } else {
    range_by_fod <- dat_sub %>%
      dplyr::group_by(fod) %>%
      dplyr::summarise(
        val_min = min(.data[[value_col]], na.rm = TRUE),
        val_max = max(.data[[value_col]], na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::mutate(fod = as.character(fod))
    
    out <- out %>%
      dplyr::left_join(range_by_fod, by = "fod") %>%
      dplyr::mutate(
        y = if (!is.null(y_limits)) {
          y_limits[1] + 0.94 * (y_limits[2] - y_limits[1])
        } else {
          val_max + 0.03 * (val_max - val_min)
        }
      )
  }
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
      aes(x = year, y = y, label = letter),
      inherit.aes = FALSE, size = 3.5, fontface = "bold"
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
    # coord_cartesian recadre l'affichage sans retirer de donnees
    # (contrairement a scale_y_continuous(limits=...)), donc les
    # violons/stats hors bornes sont juste coupes a l'affichage,
    # pas supprimes du calcul.
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
#    confondus (equivalent de compute_interannual_stats / plot_interannual)
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
      aes(x = year, y = y, label = letter),
      inherit.aes = FALSE, size = 3.5, fontface = "bold"
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
# 5. Diagrammes du nombre de donnees effectives (jours) utilisees
#    pour calculer chaque moyenne -- PAS de statistiques (pas de
#    barre d'erreur, pas de lettres). Meme granularite de comptage
#    que daily_means dans compute_stats_by_fod / compute_interannual_
#    stats_generic (1 jour = 1 date ou au moins 1 pixel valide),
#    pour visualiser directement l'effectif derriere chaque violon.
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
# 6. Execution -- ratios pigment/Chla
#    Echelle LINEAIRE, ylim adapte a chaque ratio (99.9e percentile
#    et sommet des barres d'erreur, avec marge de 25%), pour rendre
#    lisibles les ratios a faible concentration sans jamais couper
#    une barre d'erreur.
# ============================================================

for (rv in ratio_vars) {
  rv_vals  <- ratio_long$ratio[ratio_long$ratio_name == rv & ratio_long$year %in% years_keep]
  stats_rv <- compute_stats_by_fod(ratio_long, value_col = "ratio", entity_col = "ratio_name",
                                   entity_val = rv, years_keep = years_keep, log_transform = FALSE)
  
  rv_upper <- max(
    stats::quantile(rv_vals, probs = 0.999, na.rm = TRUE),
    max(stats_rv$ymax, na.rm = TRUE)
  )
  y_lim_rv <- c(0, rv_upper * 1.25)   # marge genereuse pour laisser respirer violon + lettres
  
  plot_by_fod(ratio_long, value_col = "ratio", entity_col = "ratio_name",
              entity_val = rv, label = rv, years_keep = years_keep,
              log_transform = FALSE, y_limits = y_lim_rv)
  # plot_n_days_by_fod(ratio_long, value_col = "ratio", entity_col = "ratio_name",
  #                    entity_val = rv, label = rv, years_keep = years_keep)
}

# 6.2 Distribution interannuelle, tous FOD confondus, un plot facette par ratio
plot_interannual_generic(ratio_long, value_col = "ratio", entity_col = "ratio_name",
                         label = "ratios_pigment_chla", years_keep = years_keep,
                         log_transform = FALSE)
plot_n_days_interannual(ratio_long, value_col = "ratio", entity_col = "ratio_name",
                        label = "ratios_pigment_chla", years_keep = years_keep)

# ============================================================
# 7. Execution -- FTLE (echelle lineaire, axe libre par facette)
# ============================================================

# 7.1 Par cluster FOD x annee (echelle lineaire)
plot_by_fod(ftle_long, value_col = "ftle", entity_col = "variable",
            entity_val = "FTLE", label = "FTLE", years_keep = years_keep, log_transform = FALSE)
plot_n_days_by_fod(ftle_long, value_col = "ftle", entity_col = "variable",
                   entity_val = "FTLE", label = "FTLE", years_keep = years_keep)

# 7.2 Distribution interannuelle, tous FOD confondus (echelle lineaire)
plot_interannual_generic(ftle_long, value_col = "ftle", entity_col = "variable",
                         label = "FTLE", years_keep = years_keep, log_transform = FALSE)
plot_n_days_interannual(ftle_long, value_col = "ftle", entity_col = "variable",
                        label = "FTLE", years_keep = years_keep)

# ============================================================
# 8. Execution -- NASC
# ============================================================

freq <- 38
nasc_ds <- readRDS(paste0(
  "F:/data_elise/ds_NASC_pig_ftle_fod/ds_NASC_per_esu_all/ds_NASC_per_esu_pig_ftle_fod_2018_2021_2022_2023_transect_",
  freq, "kHz_mask9.rds"
))

# ------------------------------------------------------------
# Mise en forme -- meme structure que ratio_long / ftle_long
# (colonnes date, year, fod, valeur + colonne "variable" factice
# pour reutiliser plot_by_fod / plot_interannual_generic tel quel)
# ------------------------------------------------------------

nasc_long <- nasc_ds %>%
  dplyr::mutate(
    date = as.Date(time_nasc),
    year = format(time_nasc, "%Y"),
    
    # Suppression des espaces avant/après les codes FOD
    fod = trimws(as.character(fod)),
    
    # "NA" et NA réel traités comme NA
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
# 8.1 Par cluster FOD x annee (echelle log, comme les pigments d'origine)
plot_by_fod(nasc_long, value_col = "nasc", entity_col = "variable",
            entity_val = "NASC", label = "NASC", years_keep = years_keep,
            log_transform = TRUE)
plot_n_days_by_fod(nasc_long, value_col = "nasc", entity_col = "variable",
                   entity_val = "NASC", label = "NASC", years_keep = years_keep)

# 8.2 Distribution interannuelle, tous FOD confondus (echelle log)
plot_interannual_generic(nasc_long, value_col = "nasc", entity_col = "variable",
                         label = "NASC", years_keep = years_keep, log_transform = TRUE)
plot_n_days_interannual(nasc_long, value_col = "nasc", entity_col = "variable",
                        label = "NASC", years_keep = years_keep)

sort(unique(nasc_ds$fod))
