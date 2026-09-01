# =====================================================================
# Variogrammes empiriques des RÉSIDUS d'un random forest (NASC ~ fod +
# ftle + pigments) -- version révisée
# =====================================================================
# Corrections apportées suite au diagnostic précédent :
#  - seed FIXE dans ranger() -- sans ça, deux exécutions donnent des
#    résidus OOF différents (randomisation interne de la forêt), donc des
#    variogrammes différents, sans que ce soit un vrai changement de signal.
#  - RÉPÉTITION (n_repeats) du calcul de résidus OOF avec différents seeds,
#    pour distinguer un vrai signal d'un artefact d'une seule exécution.
#  - RÉSOLUTION TEMPORELLE COMPLÈTE (pas de round à la journée) -- tes ESU
#    sont espacées de 3-4 secondes, un arrondi à la journée (as.Date)
#    écrasait complètement cette échelle et aplatissait le variogramme.
#  - BOUNDARIES MULTI-ÉCHELLES (fines près de zéro, grossières au loin) --
#    un seul `width` fixe ne peut pas couvrir à la fois la seconde et les
#    ~60 jours d'une campagne.
#  - SOUS-ÉCHANTILLONNAGE pour le calcul du variogramme (pas pour l'entraî-
#    nement du modèle) -- des boundaries fines sur ~70k points donnent des
#    milliards de paires proches ; un sous-échantillon de quelques milliers
#    de points donne déjà largement assez de paires par classe tout en
#    restant rapide à calculer, y compris répété plusieurs fois.
#
# Résidus utilisés : OUT-OF-FOLD d'une CV ALÉATOIRE classique (PAS la CV
# bloquée spatiale/temporelle) -- toujours pour la même raison qu'avant :
# ne pas mesurer la dépendance avec un outil qui l'a déjà supprimée.
#
# install.packages(c("gstat", "sp"))

library(dplyr)
library(ranger)
library(gstat)
library(sp)
library(ggplot2)
library(purrr)
freq <- 18
diurnal_period <- 3 # 3 : day, 1: night
dp <- "day"
path_ds <- paste0("F:/data_elise/ds_NASC_pig_ftle_fod/ds_NASC_per_esu_all/ds_NASC_per_esu_pig_ftle_fod_2018_2021_2022_2023_transect_", freq, "kHz_mask9.rds")
datas <- readRDS(path_ds)

# datas <- datas[datas$day == diurnal_period, ]

q <- quantile(datas$nasc, probs = c(0.05, 0.95), na.rm = TRUE)
datas <- datas |> dplyr::filter(nasc >= q[1], nasc <= q[2])
datas$nasc <- log10(datas$nasc)

datas$fod <- as.factor(datas$fod)
datas$fod[datas$fod == "NA"] <- NA
datas$fod <- droplevels(datas$fod)

df <- data.frame(
  NASC            = datas$nasc,
  year            = format(datas$time_nasc, "%Y"),
  time            = datas$time_nasc,
  day             = datas$day,
  lat             = datas$lat_nasc,
  lon             = datas$lon_nasc,
  fod             = datas$fod,
  ftle            = datas$ftle,
  total_chla      = datas$Chla,
  per_ratio_chla  = datas$Per_Chla,
  but_ratio_chla  = datas$But_Chla,
  fuco_ratio_chla = datas$Fuco_Chla,
  hex_ratio_chla  = datas$Hex_Chla,
  allo_ratio_chla = datas$Allo_Chla,
  zea_ratio_chla  = datas$Zea_Chla,
  chlb_ratio_chla = datas$Chlb_Chla
)

vars_num <- c("NASC", "per_ratio_chla", "but_ratio_chla", "fuco_ratio_chla",
              "hex_ratio_chla", "allo_ratio_chla", "zea_ratio_chla",
              "chlb_ratio_chla", "total_chla", "ftle")

df <- df |>
  dplyr::filter(
    if_all(all_of(setdiff(vars_num, "NASC")), ~ !is.na(.)),
    !is.na(fod), fod != "NA"
  )

covariates_num <- setdiff(vars_num, "NASC")
covariates_all <- c(covariates_num, "fod")
response_var   <- "NASC"

cat("Nombre d'observations après filtrage :", nrow(df), "\n")
lon0 <- mean(df$lon, na.rm = TRUE)
lat0 <- mean(df$lat, na.rm = TRUE)
km_per_deg_lat <- 110.574
km_per_deg_lon <- 111.320 * cos(lat0 * pi / 180)
df$x_km <- (df$lon - lon0) * km_per_deg_lon
df$y_km <- (df$lat - lat0) * km_per_deg_lat
# ---------------------------------------------------------------------
# 0. Suppose déjà en mémoire : df, covariates_all, response_var
# ---------------------------------------------------------------------
# Ce script est volontairement AUTONOME : il ne dépend pas de
# results_spatial (qui suppose l'autre script déjà exécuté en entier avec
# son tuning). C'est logique ici -- le variogramme des résidus sert
# normalement à INFORMER le choix des blocs/buffers en amont, pas l'inverse.
# Des hyperparamètres RF par défaut raisonnables suffisent pour ce
# diagnostic (le but est de voir la structure spatiale/temporelle du
# résidu, pas d'avoir le modèle le plus optimisé possible).
#
# Si tu as déjà lancé l'autre script et que `results_spatial` existe en
# mémoire, décommente la ligne suivante pour réutiliser ses hyperparamètres :
# params <- results_spatial$best_params

params <- list(mtry = 3, min.node.size = 5, num.trees = 300)

n_repeats       <- 3      # nombre de répétitions pour juger de la stabilité
variogram_n_sub <- 8000   # taille du sous-échantillon utilisé pour les variogrammes

# ---------------------------------------------------------------------
# 1. RÉSIDUS OUT-OF-FOLD -- CV ALÉATOIRE, RÉPÉTÉE n_repeats FOIS
# ---------------------------------------------------------------------
compute_oof_residuals <- function(seed, data, response, covs, params, k = 10) {
  set.seed(seed)
  n <- nrow(data)
  fold_id_random <- sample(rep(1:k, length.out = n))
  oof_pred <- rep(NA_real_, n)
  
  for (kk in seq_len(k)) {
    train_idx <- which(fold_id_random != kk)
    test_idx  <- which(fold_id_random == kk)
    
    train_df <- data[train_idx, c(response, covs)]
    train_df <- train_df[stats::complete.cases(train_df), ]
    
    model <- ranger(
      stats::as.formula(paste(response, "~ .")), data = train_df,
      mtry = params$mtry, min.node.size = params$min.node.size,
      num.trees = params$num.trees,
      num.threads = max(1, parallel::detectCores() - 1),
      seed = seed   # <-- reproductibilité : sans ça, la randomisation interne
      #     de ranger() change les résultats à chaque exécution
    )
    
    oof_pred[test_idx] <- predict(model, data[test_idx, covs])$predictions
  }
  
  data[[response]] - oof_pred
}

residual_matrix <- map(seq_len(n_repeats), function(r) {
  cat("Répétition", r, "/", n_repeats, "...\n")
  compute_oof_residuals(seed = 100 + r, data = df, response = response_var,
                        covs = covariates_all, params = params)
})

for (r in seq_len(n_repeats)) {
  cat("Répétition", r, "-- résidu : moyenne =",
      round(mean(residual_matrix[[r]], na.rm = TRUE), 4),
      " | sd =", round(sd(residual_matrix[[r]], na.rm = TRUE), 4), "\n")
}

# ---------------------------------------------------------------------
# 1bis. CARTE DES RÉSIDUS COLORÉE PAR VALEUR -- vérifier l'hypothèse de
# sous-zones distinctes (suggérée par la forme en bosse/creux du
# variogramme spatial, incompatible avec un continuum spatial homogène)
# ---------------------------------------------------------------------
# On utilise la répétition 1 (les 3 étant quasi identiques comme vu sur les
# variogrammes) et TOUTES les données (pas le sous-échantillon utilisé pour
# les variogrammes) -- la carte bénéficie de la résolution complète.
df_map <- df %>%
  mutate(residual = residual_matrix[[1]]) %>%
  filter(!is.na(residual), !is.na(lon), !is.na(lat))

ggplot(df_map, aes(x = lon, y = lat, color = residual)) +
  geom_point(size = 0.6, alpha = 0.7) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                        name = "Résidu\n(log10 NASC)") +
  coord_quickmap() +
  labs(title = "Carte des résidus du random forest",
       subtitle = "Bleu = sous-estimation, rouge = sur-estimation -- cherche des zones cohérentes",
       x = "Longitude", y = "Latitude") +
  theme_minimal()

# Complément utile : moyenne du résidu par bloc spatial grossier (ex. grille
# de 50x50 km) -- si des blocs entiers ressortent nettement positifs ou
# négatifs, ça confirme l'hypothèse de sous-zones plutôt qu'un bruit spatial
# homogène.
df_map_grid <- df_map %>%
  mutate(block_x = floor(x_km / 50), block_y = floor(y_km / 50)) %>%
  group_by(block_x, block_y) %>%
  summarise(mean_residual = mean(residual), n = n(),
            lon = mean(lon), lat = mean(lat), .groups = "drop") %>%
  filter(n >= 20)   # écarte les blocs avec trop peu de points pour être fiables

ggplot(df_map_grid, aes(x = lon, y = lat, color = mean_residual, size = n)) +
  geom_point(alpha = 0.85) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                        name = "Résidu moyen\npar bloc (50km)") +
  coord_quickmap() +
  labs(title = "Résidu moyen par bloc spatial (50x50 km)",
       subtitle = "Des blocs entiers nettement colorés = signe de sous-zones distinctes",
       x = "Longitude", y = "Latitude", size = "n points") +
  theme_minimal()

# ---------------------------------------------------------------------
# 2. FONCTION GÉNÉRIQUE : variogramme empirique avec boundaries multi-échelles
# ---------------------------------------------------------------------
compute_variogram <- function(coords_df, residual, boundaries, coord_names) {
  d <- coords_df
  d$residual <- residual
  d <- d[stats::complete.cases(d[, c("residual", coord_names)]), ]
  coordinates(d) <- stats::as.formula(paste("~", paste(coord_names, collapse = "+")))
  variogram(residual ~ 1, data = d, boundaries = boundaries)
}

# ---------------------------------------------------------------------
# 3. VARIOGRAMME SPATIAL DES RÉSIDUS -- boundaries fines près de zéro
# ---------------------------------------------------------------------
# À 3-4 secondes d'intervalle et une vitesse de navigation typique, deux
# ESU consécutives sont espacées de quelques dizaines de mètres seulement
# -- des boundaries fines (50-100 m) près de zéro sont nécessaires pour
# voir une éventuelle décorrélation à cette échelle, avant d'élargir vers
# le km puis la dizaine de km.
spatial_boundaries <- c(
  seq(0, 2, by = 0.05),      # 0 à 2 km, pas de 50 m
  seq(2.5, 20, by = 0.5),    # 2 à 20 km, pas de 500 m
  seq(22, 100, by = 2),      # 20 à 100 km, pas de 2 km
  seq(105, 500, by = 10)     # 100 à 500 km, pas de 10 km
)

set.seed(42)
spatial_variograms <- map(seq_len(n_repeats), function(r) {
  sub_idx <- sample(seq_len(nrow(df)), min(variogram_n_sub, nrow(df)))
  vg <- compute_variogram(
    coords_df  = df[sub_idx, c("x_km", "y_km")],
    residual   = residual_matrix[[r]][sub_idx],
    boundaries = spatial_boundaries,
    coord_names = c("x_km", "y_km")
  )
  vg$repeat_id <- paste0("rep_", r)
  vg
})

spatial_variograms_df <- bind_rows(map(spatial_variograms, as.data.frame))

cat("\nVariogramme empirique spatial (résidus), toutes répétitions :\n")
print(head(spatial_variograms_df, 20))

ggplot(spatial_variograms_df, aes(x = dist, y = gamma, color = repeat_id)) +
  geom_point(alpha = 0.6) +
  geom_line(alpha = 0.4) +
  labs(title = "Variogramme empirique des résidus -- spatial",
       subtitle = paste0(n_repeats, " répétitions (sous-échantillon de ", variogram_n_sub, " points)"),
       x = "Distance (km)", y = expression(gamma(h)), color = "Répétition") +
  theme_minimal()

# Ajustement d'un modèle théorique sur CHAQUE répétition, pour comparer
# les portées obtenues -- si elles varient beaucoup, la portée n'est pas
# fiable et il vaut mieux se fier à la lecture visuelle du graphique global.
spatial_fits <- map_dfr(seq_along(spatial_variograms), function(r) {
  vg <- spatial_variograms[[r]]
  vgm_init <- vgm(psill = var(vg$gamma, na.rm = TRUE) * 2,
                  model = "Exp", range = 50, nugget = min(vg$gamma, na.rm = TRUE))
  fit <- tryCatch(fit.variogram(vg, model = vgm_init), error = function(e) NULL)
  if (is.null(fit)) return(tibble(repeat_id = paste0("rep_", r), range_km = NA_real_))
  tibble(repeat_id = paste0("rep_", r),
         nugget = fit$psill[1], sill = fit$psill[2], range_km = fit$range[2])
})

cat("\nPortées spatiales ajustées, par répétition :\n")
print(spatial_fits)
spatial_fits <- spatial_fits %>% mutate(nugget_ratio = nugget / (nugget + sill))
cat("\nRatio nugget / (nugget + palier) -- proportion de variance NON structurée :\n")
print(spatial_fits %>% select(repeat_id, nugget_ratio))
if (mean(spatial_fits$nugget_ratio, na.rm = TRUE) > 0.9) {
  cat("\n[!] Plus de 90% de la variance résiduelle spatiale est du bruit pur (nugget).\n",
      "    La portée ci-dessus décrit une composante quasi négligeable -- ne pas la\n",
      "    présenter comme un chiffre précis. La dépendance spatiale résiduelle est\n",
      "    globalement négligeable une fois les covariables prises en compte.\n")
}
cat("-> Si ces valeurs varient fortement d'une répétition à l'autre,\n",
    "   préfère une lecture visuelle du graphique poolé plutôt qu'un chiffre unique.\n")

# ---------------------------------------------------------------------
# 3bis. MODÈLE À EFFET DE TROU ("Hol") -- alternative au modèle exponentiel
# ---------------------------------------------------------------------
# Le variogramme empirique observé (bosse vers 100-250 km, creux vers
# 300-320 km, remontée après) n'est PAS une forme monotone croissante --
# aucun modèle classique (Exp/Sph/Gau) ne peut s'y ajuster correctement,
# d'où le "No convergence" systématique plus haut. Le modèle "Hol" (hole
# effect, basé sur une fonction sinus amortie) est fait pour ce genre de
# motif oscillant, souvent le signe de sous-zones spatiales distinctes
# plutôt que d'un continuum homogène.
spatial_fits_hole <- map_dfr(seq_along(spatial_variograms), function(r) {
  vg <- spatial_variograms[[r]]
  vgm_init_hole <- vgm(psill = var(vg$gamma, na.rm = TRUE),
                       model = "Hol", range = 300,
                       nugget = min(vg$gamma, na.rm = TRUE))
  fit <- tryCatch(fit.variogram(vg, model = vgm_init_hole), error = function(e) NULL)
  if (is.null(fit)) return(tibble(repeat_id = paste0("rep_", r), range_km = NA_real_))
  tibble(repeat_id = paste0("rep_", r),
         nugget = fit$psill[1], sill = fit$psill[2], range_km = fit$range[2])
})

cat("\nModèle 'Hol' (effet de trou) -- portées ajustées par répétition :\n")
print(spatial_fits_hole)
cat("-> `range_km` ici correspond à la pseudo-période de l'oscillation,\n",
    "   pas à une portée d'indépendance classique -- interprète-la comme\n",
    "   la distance caractéristique séparant les zones de résidus similaires,\n",
    "   pas comme un seuil au-delà duquel l'indépendance est acquise.\n")

# Superpose le modèle Hole ajusté (répétition 1) sur le variogramme poolé,
# pour vérifier visuellement la qualité de l'ajustement.
vgm_hole_r1 <- vgm(psill = spatial_fits_hole$sill[1], model = "Hol",
                   range = spatial_fits_hole$range_km[1], nugget = spatial_fits_hole$nugget[1])
plot(spatial_variograms[[1]], vgm_hole_r1,
     main = "Variogramme spatial -- ajustement modèle 'Hol' (répétition 1)")

# ---------------------------------------------------------------------
# 3ter. VARIOGRAMME SPATIAL RESTREINT AUX PAIRES INTRA-TRANSECT
# ---------------------------------------------------------------------
# La carte des résidus révèle que les données ne forment PAS un champ
# spatial continu, mais une poignée de transects distincts et séparés
# (lignes fines, grands vides entre elles). Le variogramme global (3) mêle
# donc des paires INTRA-transect (vraie dépendance spatiale continue) et
# des paires INTER-transects (quelques distances arbitraires entre zones
# d'étude différentes, sans rapport avec une dépendance spatiale continue)
# -- c'est probablement la cause du motif bosse/creux, pas un vrai effet
# de trou. On restreint donc le calcul aux paires DANS UN MÊME transect
# (approximé ici par le jour de campagne).
#
# Sous-échantillonnage par groupe (max_points_per_group) pour rester
# calculable : le nombre de paires croît en n² par transect.
max_points_per_group <- 800

compute_within_group_variogram <- function(data, group_var, coord_names,
                                           residual_col, boundaries,
                                           max_points_per_group = 800, seed = 1) {
  set.seed(seed)
  bin_mids <- (boundaries[-1] + boundaries[-length(boundaries)]) / 2
  n_bins   <- length(bin_mids)
  sum_sq   <- numeric(n_bins)
  sum_np   <- numeric(n_bins)
  
  groups <- unique(data[[group_var]])
  for (g in groups) {
    sub <- data[data[[group_var]] == g, ]
    if (nrow(sub) < 2) next
    if (nrow(sub) > max_points_per_group) {
      sub <- sub[sample(seq_len(nrow(sub)), max_points_per_group), ]
    }
    
    coords <- as.matrix(sub[, coord_names])
    resid  <- sub[[residual_col]]
    d      <- as.matrix(dist(coords))
    idx    <- which(upper.tri(d), arr.ind = TRUE)
    dists  <- d[idx]
    sqdiff <- (resid[idx[, 1]] - resid[idx[, 2]])^2
    
    bin_id <- findInterval(dists, boundaries)
    valid  <- bin_id >= 1 & bin_id <= n_bins & dists <= max(boundaries)
    
    agg <- tapply(sqdiff[valid], bin_id[valid], sum)
    cnt <- tapply(sqdiff[valid], bin_id[valid], length)
    b_idx <- as.integer(names(agg))
    sum_sq[b_idx] <- sum_sq[b_idx] + as.numeric(agg)
    sum_np[b_idx] <- sum_np[b_idx] + as.numeric(cnt)
  }
  
  tibble(dist = bin_mids, gamma = sum_sq / (2 * sum_np), np = sum_np) %>%
    filter(np > 0)
}

df_transect <- df %>%
  mutate(residual = residual_matrix[[1]], date = as.Date(time)) %>%
  filter(!is.na(residual), !is.na(x_km), !is.na(y_km))

vgm_within_transect <- compute_within_group_variogram(
  data = df_transect, group_var = "date", coord_names = c("x_km", "y_km"),
  residual_col = "residual", boundaries = spatial_boundaries,
  max_points_per_group = max_points_per_group
)

cat("\nVariogramme spatial INTRA-transect (paires du même transect uniquement) :\n")
print(vgm_within_transect)

ggplot(vgm_within_transect, aes(x = dist, y = gamma)) +
  geom_point(color = "darkred") +
  geom_line(color = "darkred", alpha = 0.5) +
  labs(title = "Variogramme spatial -- paires INTRA-transect uniquement",
       subtitle = "Exclut les comparaisons entre transects séparés (pas de sens en continu)",
       x = "Distance (km)", y = expression(gamma(h))) +
  theme_minimal()

# Ajustement d'un modèle exponentiel classique sur ce variogramme
# "propre" -- devrait converger correctement, sans le motif bosse/creux
# artificiel du variogramme global.
vgm_object_transect <- vgm_within_transect %>%
  mutate(dir.hor = 0, dir.ver = 0, id = "var1") %>%
  select(np, dist, gamma, dir.hor, dir.ver, id)
class(vgm_object_transect) <- c("gstatVariogram", "data.frame")

vgm_init_transect <- vgm(psill = var(vgm_within_transect$gamma, na.rm = TRUE),
                         model = "Exp", range = 30,
                         nugget = min(vgm_within_transect$gamma, na.rm = TRUE))
vgm_fit_transect <- fit.variogram(vgm_object_transect, model = vgm_init_transect)
cat("\nModèle ajusté -- variogramme intra-transect :\n")
print(vgm_fit_transect)
cat("-> Portée spatiale intra-transect estimée :",
    round(vgm_fit_transect$range[2], 1), "km\n",
    "   (c'est la portée la plus fiable pour dimensionner ton buffer spatial,\n",
    "    car elle exclut les comparaisons entre transects sans rapport spatial continu)\n")

# ---------------------------------------------------------------------
# 4. VARIOGRAMME TEMPOREL DES RÉSIDUS -- résolution complète, échelle HEURES
# ---------------------------------------------------------------------
# time_num en HEURES FRACTIONNAIRES (pas de round à la journée) -- garde la
# résolution native de tes données (3-4 secondes) tout en restant lisible
# sur un axe log (0,05h à 1440h, soit ~3 min à 60 jours).
# Boundaries : très fines les premières minutes, puis élargies
# jusqu'à ~1440 heures (60 jours, durée d'une campagne).
temporal_boundaries <- c(
  seq(0, 1, length.out = 20),      # 0 à 1h, résolution ~3 min
  seq(1.2, 24, by = 1.2),          # 1h à 1 jour, pas de ~1h12
  seq(26, 240, by = 5),            # 1 à 10 jours, pas de 5h
  seq(250, 1440, by = 24)          # 10 à 60 jours, pas de 1 jour (24h)
)

set.seed(43)
temporal_variograms <- map(seq_len(n_repeats), function(r) {
  df_t <- df %>%
    mutate(time_num = as.numeric(difftime(time, min(time, na.rm = TRUE), units = "hours")),
           dummy_y  = 0)
  sub_idx <- sample(seq_len(nrow(df_t)), min(variogram_n_sub, nrow(df_t)))
  vg <- compute_variogram(
    coords_df  = df_t[sub_idx, c("time_num", "dummy_y")],
    residual   = residual_matrix[[r]][sub_idx],
    boundaries = temporal_boundaries,
    coord_names = c("time_num", "dummy_y")
  )
  vg$repeat_id <- paste0("rep_", r)
  vg
})

temporal_variograms_df <- bind_rows(map(temporal_variograms, as.data.frame))

cat("\nVariogramme empirique temporel (résidus), toutes répétitions :\n")
print(head(temporal_variograms_df, 20))

ggplot(temporal_variograms_df, aes(x = dist, y = gamma, color = repeat_id)) +
  geom_point(alpha = 0.6) +
  geom_line(alpha = 0.4) +
  scale_x_log10() +
  labs(title = "Variogramme empirique des résidus -- temporel",
       subtitle = paste0(n_repeats, " répétitions -- axe des x en échelle log (heures)"),
       x = "Écart temporel (heures, échelle log)", y = expression(gamma(h)), color = "Répétition") +
  theme_minimal()

temporal_fits <- map_dfr(seq_along(temporal_variograms), function(r) {
  vg <- temporal_variograms[[r]]
  vg <- vg[vg$np >= 30, ]   # écarte les classes peu fiables avant ajustement
  vgm_init <- vgm(psill = var(vg$gamma, na.rm = TRUE) * 2,
                  model = "Exp", range = 120, nugget = min(vg$gamma, na.rm = TRUE))
  fit <- tryCatch(fit.variogram(vg, model = vgm_init), error = function(e) NULL)
  if (is.null(fit)) return(tibble(repeat_id = paste0("rep_", r), range_hours = NA_real_))
  tibble(repeat_id = paste0("rep_", r),
         nugget = fit$psill[1], sill = fit$psill[2], range_hours = fit$range[2])
})

cat("\nPortées temporelles ajustées, par répétition :\n")
print(temporal_fits)
temporal_fits <- temporal_fits %>% mutate(nugget_ratio = nugget / (nugget + sill))
cat("\nRatio nugget / (nugget + palier) -- proportion de variance NON structurée :\n")
print(temporal_fits %>% select(repeat_id, nugget_ratio))
if (mean(temporal_fits$nugget_ratio, na.rm = TRUE) > 0.9) {
  cat("\n[!] Plus de 90% de la variance résiduelle temporelle est du bruit pur (nugget).\n",
      "    La portée ci-dessus décrit une composante quasi négligeable -- ne pas la\n",
      "    présenter comme un chiffre précis. La dépendance temporelle résiduelle est\n",
      "    globalement négligeable une fois les covariables prises en compte.\n")
}
cat("(Pour convertir en jours : diviser range_hours par 24.)\n")

# ---------------------------------------------------------------------
# 4bis. VARIOGRAMME TEMPOREL RESTREINT AUX PAIRES INTRA-TRANSECT
# ---------------------------------------------------------------------
# Même correction que pour le spatial (section 3ter) : si deux transects
# différents ont eu lieu à quelques heures/jours d'écart pendant une même
# campagne, le variogramme temporel global (section 4) les compare comme
# un seul continuum alors que ce sont deux endroits différents -- d'où
# probablement les pics erratiques au-delà de 5-10h. On restreint donc aux
# paires DANS UN MÊME transect (même fonction que pour le spatial, juste
# avec le temps comme coordonnée).
df_transect_t <- df %>%
  mutate(residual = residual_matrix[[1]], date = as.Date(time),
         time_num = as.numeric(difftime(time, min(time, na.rm = TRUE), units = "hours")),
         dummy_y  = 0) %>%
  filter(!is.na(residual), !is.na(time_num))

vgm_within_transect_t <- compute_within_group_variogram(
  data = df_transect_t, group_var = "date", coord_names = c("time_num", "dummy_y"),
  residual_col = "residual", boundaries = temporal_boundaries,
  max_points_per_group = max_points_per_group
)

cat("\nVariogramme temporel INTRA-transect (paires du même transect uniquement) :\n")
print(vgm_within_transect_t)

ggplot(vgm_within_transect_t, aes(x = dist, y = gamma)) +
  geom_point(color = "darkorange") +
  geom_line(color = "darkorange", alpha = 0.5) +
  scale_x_log10() +
  labs(title = "Variogramme temporel -- paires INTRA-transect uniquement",
       subtitle = "Exclut les comparaisons entre transects séparés dans le temps",
       x = "Écart temporel (heures, échelle log)", y = expression(gamma(h))) +
  theme_minimal()

# Ajustement d'un modèle exponentiel classique sur ce variogramme "propre"
vgm_object_transect_t <- vgm_within_transect_t %>%
  filter(np >= 30) %>%   # écarte les classes peu fiables avant ajustement
  mutate(dir.hor = 0, dir.ver = 0, id = "var1") %>%
  select(np, dist, gamma, dir.hor, dir.ver, id)
class(vgm_object_transect_t) <- c("gstatVariogram", "data.frame")

vgm_init_transect_t <- vgm(psill = var(vgm_within_transect_t$gamma, na.rm = TRUE),
                           model = "Exp", range = 10,
                           nugget = min(vgm_within_transect_t$gamma, na.rm = TRUE))
vgm_fit_transect_t <- fit.variogram(vgm_object_transect_t, model = vgm_init_transect_t)
cat("\nModèle ajusté -- variogramme temporel intra-transect :\n")
print(vgm_fit_transect_t)
cat("-> Portée temporelle intra-transect estimée :",
    round(vgm_fit_transect_t$range[2], 1), "heures (soit",
    round(vgm_fit_transect_t$range[2] / 24, 2), "jours)\n",
    "   (c'est la portée la plus fiable pour dimensionner ton buffer temporel,\n",
    "    car elle exclut les comparaisons entre transects distincts dans le temps)\n")

print(vgm_within_transect$np)
print(vgm_within_transect_t[vgm_within_transect_t$dist > 15, ])
# ---------------------------------------------------------------------
# 5. LECTURE DES RÉSULTATS
# ---------------------------------------------------------------------
# - Regarde D'ABORD les graphiques poolés (toutes répétitions superposées)
#   avant de te fier à un chiffre de portée unique. Une forme cohérente
#   entre répétitions = signal réel ; une forme qui change du tout au tout
#   = pas assez de signal pour conclure à une portée précise.
# - L'axe log sur le graphique temporel permet de voir à la fois l'échelle
#   fine (minutes/heures, near-field) et l'échelle large (jours) sur le
#   même graphique -- regarde si gamma monte tôt (dépendance à courte
#   échelle uniquement) ou continue de monter sur toute la plage (structure
#   à plus grande échelle encore présente).
# - N'oublie pas la mise en garde précédente : le variogramme temporel
#   intra-campagne mélange partiellement temps et espace, puisque le
#   bateau se déplace continûment pendant l'acquisition.