library(dplyr)
library(tidyr)
library(ggplot2)
library(ranger)
options(scipen = 999) 

# Data validation : données manquantes
path <- "F:/data_elise/ds_NASC_pig_ftle_fod/ds_NASC_per_esu_all/new_ratio_ds_NASC_per_esu_pig_ftle_fod_2018_2021_2022_2023_transect_38kHz_mask9.rds"

ds <- readRDS(path)
str(ds)
ds <- ds[ds$day == 3, ]

# Cartes de valeurs manquantes par région 
library(tidyr)
library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)
library(rnaturalearth)
set.seed(42)

output_dir <- "C:/Users/mmolinet/elisou_ta_stagiaire_pref/MIO_internship_III/src/preprocess/validation_donnees"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Sauvegarde d'un graphique (affiché aussi à l'écran)
sauver_plot <- function(p, nom, width = 25, height = 18) {
  print(p)
  ggsave(file.path(output_dir, paste0(nom, ".png")), p,
         width = width, height = height, units = "cm", dpi = 300)
}

# Sauvegarde d'un tableau (format CSV français : ";" et virgule décimale)
sauver_tab <- function(tab, nom) {
  write.csv2(tab, file.path(output_dir, paste0(nom, ".csv")), row.names = FALSE)
}
# =====================================================================
# 1. Préparation
# =====================================================================

pig_ratio   <- c("Per_totpig", "But_totpig", "Fuco_totpig", "Hex_totpig",
                 "Allo_totpig", "Zea_totpig", "Chlb_totpig", "DvChla_totpig")
covars      <- c("ftle", "Chla", "total_pig", pig_ratio)   # fod exclue : c'est la partition
covars_phys <- c("ftle")                                   # jeu sans pigments (moins de NA)

dat <- ds %>%
  filter(!is.na(nasc), !is.na(lat_nasc), !is.na(lon_nasc), !is.na(time_nasc)) %>%
  mutate(zone     = factor(na_if(fod, "NA")),
         annee    = factor(format(time_nasc, "%Y")),
         bloc     = as.character(as.Date(time_nasc)),   # blocs journaliers pour la CV
         log_nasc = log10(nasc + 1))

# Contrôles : effectifs par zone, et croisement années x zones
table(dat$zone, useNA = "ifany")
table(dat$annee, dat$zone, useNA = "ifany")

# =====================================================================
# 2. Décalage univarié (covariables + NASC)
# =====================================================================

long <- dat %>%
  select(annee, zone, all_of(covars), log_nasc) %>%
  pivot_longer(c(all_of(covars), log_nasc),
               names_to = "variable", values_to = "valeur") %>%
  filter(is.finite(valeur)) %>%
  mutate(variable = factor(variable, levels = c("log_nasc", covars)))

# Densités superposées
densites <- function(part) {
  ggplot(filter(long, !is.na(.data[[part]])),
         aes(x = valeur, colour = .data[[part]])) +
    geom_density(linewidth = 0.5) +
    facet_wrap(~ variable, scales = "free") +
    labs(x = NULL, y = "Densité", colour = part) +
    theme_bw()
}

# Coefficient de recouvrement de deux densités (1 = identiques, 0 = disjointes)
ovl <- function(a, b, n = 512) {
  rng <- range(c(a, b))
  if (diff(rng) == 0) return(1)
  da <- density(a, from = rng[1], to = rng[2], n = n)
  db <- density(b, from = rng[1], to = rng[2], n = n)
  sum(pmin(da$y, db$y)) * diff(da$x[1:2])
}

# Tailles d'effet par variable et par modalité (modalité vs reste)
decalage_univ <- function(part) {
  long %>%
    group_by(variable) %>%
    group_modify(function(df, key) {
      df <- df[!is.na(df[[part]]), ]
      g  <- droplevels(df[[part]])
      bind_rows(lapply(levels(g), function(m) {
        a <- df$valeur[g == m]
        b <- df$valeur[g != m]
        if (length(a) < 20 || length(b) < 20) return(NULL)
        tibble(modalite       = m,
               n              = length(a),
               D_ks           = unname(suppressWarnings(ks.test(a, b)$statistic)),
               recouvrement   = ovl(a, b),
               pct_hors_plage = 100 * mean(a < min(b) | a > max(b)))
      }))
    }) %>%
    ungroup()
}

# Carte de chaleur d'une métrique
heat <- function(tab, metrique, titre) {
  ggplot(tab, aes(x = modalite, y = variable, fill = .data[[metrique]])) +
    geom_tile(colour = "white") +
    geom_text(aes(label = round(.data[[metrique]], 2)), size = 3) +
    scale_fill_viridis_c(name = metrique) +
    labs(title = titre, x = NULL, y = NULL) +
    theme_minimal()
}

densites("annee")
densites("zone")

u_an   <- decalage_univ("annee")
u_zone <- decalage_univ("zone")

heat(u_an,   "D_ks",           "D de Kolmogorov-Smirnov (année vs autres années)")
heat(u_an,   "recouvrement",   "Recouvrement des densités (année vs autres années)")
heat(u_an,   "pct_hors_plage", "% d'observations hors de la plage des autres années")
heat(u_zone, "D_ks",           "D de Kolmogorov-Smirnov (zone vs autres zones)")
heat(u_zone, "recouvrement",   "Recouvrement des densités (zone vs autres zones)")
heat(u_zone, "pct_hors_plage", "% d'observations hors de la plage des autres zones")

# =====================================================================
# 3. Validation adversariale (décalage multivarié)
# =====================================================================

# AUC par la formule des rangs (Mann-Whitney)
auc <- function(y, p) {
  r  <- rank(p)
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

validation_adv <- function(part, vars, n_max = 20000, k = 5, num.trees = 300) {
  d0 <- dat %>%
    select(all_of(c(part, "bloc", vars))) %>%
    filter(!is.na(.data[[part]]), if_all(all_of(vars), is.finite))
  mods <- levels(droplevels(d0[[part]]))
  
  res <- lapply(mods, function(m) {
    # classe 1 = la modalité, classe 0 = le reste ; sous-échantillonnage équilibré
    d <- d0 %>%
      mutate(y = factor(as.integer(.data[[part]] == m))) %>%
      group_by(y) %>%
      slice_sample(n = n_max) %>%
      ungroup()
    if (n_distinct(d$y) < 2) return(NULL)
    
    # CV par blocs journaliers
    blocs  <- unique(d$bloc)
    d$fold <- sample(rep_len(1:k, length(blocs)))[match(d$bloc, blocs)]
    
    p <- numeric(nrow(d)); imp <- 0
    for (f in 1:k) {
      tr <- d$fold != f
      rf <- ranger(y ~ ., data = d[tr, c("y", vars)], probability = TRUE,
                   num.trees = num.trees, importance = "permutation")
      p[!tr] <- predict(rf, d[!tr, vars])$predictions[, "1"]
      imp <- imp + rf$variable.importance / k
    }
    
    list(auc = tibble(modalite = m, n = nrow(d),
                      auc = auc(as.integer(as.character(d$y)), p)),
         imp = tibble(modalite = m, variable = names(imp), importance = imp))
  })
  
  res <- Filter(Negate(is.null), res)
  list(auc = bind_rows(lapply(res, `[[`, "auc")),
       imp = bind_rows(lapply(res, `[[`, "imp")))
}

adv <- list(
  annee_phys = validation_adv("annee", covars_phys),
  annee_tout = validation_adv("annee", covars),
  zone_phys  = validation_adv("zone",  covars_phys),
  zone_tout  = validation_adv("zone",  covars)
)

# AUC par modalité
auc_tab <- bind_rows(lapply(adv, `[[`, "auc"), .id = "analyse")
auc_tab

ggplot(auc_tab, aes(x = modalite, y = auc)) +
  geom_col(fill = "steelblue") +
  geom_hline(yintercept = 0.5, linetype = 2) +
  facet_wrap(~ analyse, scales = "free_x") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = NULL, y = "AUC adversariale (modalité vs reste)") +
  theme_bw()

# Variables portant le décalage
heat(adv$annee_tout$imp, "importance", "Importance (permutation) - années")
heat(adv$zone_tout$imp,  "importance", "Importance (permutation) - zones")

# ---- Section 2 : décalage univarié ----

sauver_plot(densites("annee"), "densites_annee", width = 30, height = 22)
sauver_plot(densites("zone"),  "densites_zone",  width = 30, height = 22)

u_an   <- decalage_univ("annee")
u_zone <- decalage_univ("zone")
sauver_tab(u_an,   "univarie_annee")
sauver_tab(u_zone, "univarie_zone")

sauver_plot(heat(u_an, "D_ks", "D de Kolmogorov-Smirnov (année vs autres années)"),
            "ks_annee")
sauver_plot(heat(u_an, "recouvrement", "Recouvrement des densités (année vs autres années)"),
            "recouvrement_annee")
sauver_plot(heat(u_an, "pct_hors_plage", "% d'observations hors de la plage des autres années"),
            "hors_plage_annee")
sauver_plot(heat(u_zone, "D_ks", "D de Kolmogorov-Smirnov (zone vs autres zones)"),
            "ks_zone")
sauver_plot(heat(u_zone, "recouvrement", "Recouvrement des densités (zone vs autres zones)"),
            "recouvrement_zone")
sauver_plot(heat(u_zone, "pct_hors_plage", "% d'observations hors de la plage des autres zones"),
            "hors_plage_zone")

# ---- Section 3 : validation adversariale ----

auc_tab <- bind_rows(lapply(adv, `[[`, "auc"), .id = "analyse")
imp_tab <- bind_rows(lapply(adv, `[[`, "imp"), .id = "analyse")
sauver_tab(auc_tab, "adversarial_auc")
sauver_tab(imp_tab, "adversarial_importance")

p_auc <- ggplot(auc_tab, aes(x = modalite, y = auc)) +
  geom_col(fill = "steelblue") +
  geom_hline(yintercept = 0.5, linetype = 2) +
  facet_wrap(~ analyse, scales = "free_x") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = NULL, y = "AUC adversariale (modalité vs reste)") +
  theme_bw()
sauver_plot(p_auc, "adversarial_auc")

sauver_plot(heat(adv$annee_tout$imp, "importance", "Importance (permutation) - années"),
            "adversarial_importance_annee")
sauver_plot(heat(adv$zone_tout$imp, "importance", "Importance (permutation) - zones"),
            "adversarial_importance_zone")

# Objets complets, pour recharger sans tout recalculer (readRDS)
saveRDS(list(u_an = u_an, u_zone = u_zone, adv = adv),
        file.path(output_dir, "resultats_decalage.rds"))