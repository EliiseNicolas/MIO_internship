# =====================================================================
# 02_folds.R -- construction des schémas de validation croisée
# =====================================================================
# Toutes les fonctions prennent le data.frame `df` en paramètre explicite
# (pas de variable globale) afin d'être réutilisables telles quelles pour
# le tuning ET pour l'entraînement final -- avec les MÊMES folds.

rmse_fn <- function(obs, pred) sqrt(mean((obs - pred)^2, na.rm = TRUE))

# ---------------------------------------------------------------------
# A. NAIVE RANDOM SPLIT (RS) 80/20 -- deux methodes possibles :
#
#   "kfold_repeated" (DEFAUT recommande) : repetitions de K-fold CV
#     classique (partition disjointe), K derive de frac_train (K=5 pour
#     frac_train=0.8, chaque fold test = 1/K = 20%). Chaque repetition
#     est une VRAIE partition -- chaque observation testee exactement
#     une fois par repetition, sans chevauchement entre folds d'une
#     meme repetition. Avec n_folds=10 et K=5, ca donne 2 repetitions
#     de 5-fold. Preferable au Monte-Carlo : couverture complete et
#     disjointe, pas de sous-estimation de la variance inter-fold liee
#     au chevauchement des tests, et coherent avec la logique de
#     partition des schemas bloques (comparaison naive vs bloque plus
#     symetrique).
#
#   "monte_carlo" (ancien defaut, gardee pour comparaison/compatibilite) :
#     n_folds tirages aleatoires INDEPENDANTS de 80% du jeu de donnees
#     pour le train a chaque fois. Les jeux de test se CHEVAUCHENT d'un
#     tirage a l'autre (pas une partition) -- a eviter par defaut, mais
#     utile pour retrouver exactement le comportement des runs precedents.
# ---------------------------------------------------------------------
build_naive_folds <- function(df, n_folds = NAIVE_N_FOLDS, frac_train = 0.8,
                               method = NAIVE_CV_METHOD, seed = 123) {
  switch(method,
    "kfold_repeated" = build_naive_folds_kfold_repeated(df, n_folds, frac_train, seed),
    "monte_carlo"    = build_naive_folds_monte_carlo(df, n_folds, frac_train, seed),
    stop("method inconnue : '", method, "' (attendu : 'kfold_repeated' ou 'monte_carlo')")
  )
}

build_naive_folds_kfold_repeated <- function(df, n_folds, frac_train, seed) {
  n <- nrow(df)
  k_per_repeat <- round(1 / (1 - frac_train))   # ex: frac_train=0.8 -> k=5
  n_repeats    <- ceiling(n_folds / k_per_repeat)

  set.seed(seed)
  folds <- list()
  fold_counter <- 0

  for (r in seq_len(n_repeats)) {
    perm  <- sample(seq_len(n))                 # nouvelle permutation a chaque repetition
    block <- cut(seq_len(n), breaks = k_per_repeat, labels = FALSE)  # k blocs contigus dans l'ordre permute

    for (k in seq_len(k_per_repeat)) {
      fold_counter <- fold_counter + 1
      if (fold_counter > n_folds) break

      test_idx  <- perm[block == k]
      train_idx <- setdiff(seq_len(n), test_idx)
      folds[[fold_counter]] <- list(train = train_idx, test = test_idx)
    }
  }

  names(folds) <- paste0("rs_", seq_along(folds))
  list(data = df, folds = folds, scheme = "naive_RS_80_20",
       cv_method = "kfold_repeated", k_per_repeat = k_per_repeat, n_repeats = n_repeats)
}

build_naive_folds_monte_carlo <- function(df, n_folds, frac_train, seed) {
  set.seed(seed)
  n <- nrow(df)
  folds <- map(seq_len(n_folds), function(k) {
    train_idx <- sample(seq_len(n), size = round(frac_train * n))
    test_idx  <- setdiff(seq_len(n), train_idx)
    list(train = train_idx, test = test_idx)
  })
  names(folds) <- paste0("rs_", seq_len(n_folds))
  list(data = df, folds = folds, scheme = "naive_RS_80_20", cv_method = "monte_carlo")
}

# ---------------------------------------------------------------------
# B. BLOCAGE SPATIAL (grille en km + buffer autour du bloc test)
# ---------------------------------------------------------------------
assign_spatial_block <- function(data, cellsize_x, cellsize_y) {
  data %>% mutate(
    block_x       = floor(x_km / cellsize_x),
    block_y       = floor(y_km / cellsize_y),
    spatial_block = paste(block_x, block_y, sep = "_")
  )
}

make_spatial_fold_buffered <- function(data, block_id, buffer_km) {
  test_idx            <- which(data$spatial_block == block_id)
  candidate_train_idx <- which(data$spatial_block != block_id)
  test_mat      <- as.matrix(data[test_idx, c("x_km", "y_km")])
  candidate_mat <- as.matrix(data[candidate_train_idx, c("x_km", "y_km")])
  dist_to_test  <- FNN::knnx.dist(data = test_mat, query = candidate_mat, k = 1)[, 1]
  train_idx     <- candidate_train_idx[dist_to_test > buffer_km]
  list(train = train_idx, test = test_idx)
}

# ---------------------------------------------------------------------
# Selectionne n_target blocs parmi blocks_used, en garantissant une
# couverture du bas, du milieu ET du haut de la plage de LATITUDE des
# blocs eligibles -- au lieu d'un tirage uniforme qui peut, par hasard,
# sur-representer une bande de latitude et sous-representer les autres.
# Important pour la comparaison inter-folds : deux blocs a la meme
# distance (km) du train peuvent avoir une difficulte tres differente
# selon qu'ils sont au centre ou au bord de la plage de latitude
# (extrapolation vers un regime moyen jamais vu, cf. discussion sur la
# stationnarite de 1er ordre / le gradient latitudinal).
# ---------------------------------------------------------------------
stratify_blocks_by_latitude <- function(data_blocked, blocks_used, n_target, seed) {
  if (length(blocks_used) <= n_target) return(blocks_used)

  block_lat <- data_blocked %>%
    filter(spatial_block %in% blocks_used) %>%
    group_by(spatial_block) %>%
    summarise(mean_lat = mean(lat, na.rm = TRUE), .groups = "drop")

  breaks <- stats::quantile(block_lat$mean_lat, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE, type = 1)
  breaks[1] <- breaks[1] - 1e-6   # garantit l'inclusion du minimum dans le 1er intervalle
  block_lat$lat_stratum <- cut(block_lat$mean_lat, breaks = breaks,
                                labels = c("bas", "milieu", "haut"), include.lowest = TRUE)

  set.seed(seed)
  strata   <- c("bas", "milieu", "haut")
  quota    <- rep(n_target %/% 3, 3)
  reste    <- n_target %% 3
  if (reste > 0) quota[seq_len(reste)] <- quota[seq_len(reste)] + 1  # distribue le reste

  selected <- character(0)
  for (i in seq_along(strata)) {
    pool <- block_lat$spatial_block[block_lat$lat_stratum == strata[i]]
    take <- min(quota[i], length(pool))
    if (take > 0) selected <- c(selected, sample(pool, take))
  }

  # Si un stratum manque de blocs pour atteindre son quota, complete
  # depuis le reste du pool disponible (garantit qu'on atteint n_target
  # meme si la repartition bas/milieu/haut est tres inegale).
  if (length(selected) < n_target) {
    reste_pool <- setdiff(blocks_used, selected)
    manquant   <- n_target - length(selected)
    if (length(reste_pool) > 0) {
      selected <- c(selected, sample(reste_pool, min(manquant, length(reste_pool))))
    }
  }

  selected
}

build_spatial_folds <- function(df, cellsize_x, cellsize_y, buffer_km,
                                 min_block_n         = BLOCK_MIN_BLOCK_N,
                                 max_folds_fraction  = BLOCK_MAX_FOLDS_FRACTION,
                                 min_folds           = BLOCK_MIN_FOLDS,
                                 max_folds_abs       = BLOCK_MAX_FOLDS_ABS,
                                 seed = 42) {
  data_blocked <- assign_spatial_block(df, cellsize_x, cellsize_y)
  blocks_count <- data_blocked %>% count(spatial_block)
  blocks_used  <- blocks_count %>% filter(n >= min_block_n) %>% pull(spatial_block)

  # Nombre de folds ADAPTATIF : une fraction des blocs disponibles,
  # borné par [min_folds, max_folds_abs, nb de blocs disponibles].
  n_target <- round(length(blocks_used) * max_folds_fraction)
  n_target <- max(min_folds, min(n_target, max_folds_abs, length(blocks_used)))

  set.seed(seed)
  if (length(blocks_used) > n_target) {
    blocks_used <- stratify_blocks_by_latitude(data_blocked, blocks_used, n_target, seed)
  }
  cat(sprintf(
    "  %sx%skm : %d blocs disponibles (>= %d obs) -> %d folds retenus, stratifies par latitude (cible adaptative : %d)\n",
    cellsize_x, cellsize_y, length(blocks_count %>% filter(n >= min_block_n) %>% pull(spatial_block)),
    min_block_n, length(blocks_used), n_target
  ))

  folds <- map(blocks_used, ~ make_spatial_fold_buffered(data_blocked, .x, buffer_km))
  names(folds) <- paste0("s_", blocks_used)

  valid <- map_lgl(folds, ~ length(.x$train) > 0 && length(.x$test) > 0)
  if (any(!valid)) cat("  [!] ", sum(!valid), "bloc(s) spatial(aux) ignore(s) : train ou test vide apres buffer\n")
  folds <- folds[valid]

  list(
    data = data_blocked, folds = folds, scheme = "spatial_block",
    cellsize_x = cellsize_x, cellsize_y = cellsize_y, buffer_km = buffer_km,
    n_blocks_total = nrow(blocks_count), n_blocks_used = length(folds),
    n_folds_target = n_target
  )
}

# ---------------------------------------------------------------------
# C. BLOCAGE TEMPOREL (blocs de N jours + buffer)
# ---------------------------------------------------------------------
fast_min_gap_days <- function(train_dates, test_dates) {
  train_num   <- as.numeric(as.Date(train_dates))
  test_sorted <- sort(unique(as.numeric(as.Date(test_dates))))
  idx       <- findInterval(train_num, test_sorted)
  left_val  <- test_sorted[pmax(idx, 1)]
  right_val <- test_sorted[pmin(idx + 1, length(test_sorted))]
  pmin(abs(train_num - left_val), abs(train_num - right_val))
}

assign_temporal_block <- function(data, block_days) {
  t0 <- min(as.numeric(as.Date(data$time)))
  data %>% mutate(temporal_block = floor((as.numeric(as.Date(time)) - t0) / block_days))
}

make_temporal_fold <- function(data, block_id, buffer_days) {
  test_idx            <- which(data$temporal_block == block_id)
  candidate_train_idx <- which(data$temporal_block != block_id)
  gaps      <- fast_min_gap_days(data$time[candidate_train_idx], data$time[test_idx])
  train_idx <- candidate_train_idx[gaps > buffer_days]
  list(train = train_idx, test = test_idx)
}

build_temporal_folds <- function(df, block_days, buffer_days,
                                  min_block_n         = BLOCK_MIN_BLOCK_N,
                                  max_folds_fraction  = BLOCK_MAX_FOLDS_FRACTION,
                                  min_folds           = BLOCK_MIN_FOLDS,
                                  max_folds_abs       = BLOCK_MAX_FOLDS_ABS,
                                  seed = 42) {
  data_blocked <- assign_temporal_block(df, block_days)
  blocks_count <- data_blocked %>% count(temporal_block)
  blocks_used  <- blocks_count %>% filter(n >= min_block_n) %>% pull(temporal_block)

  # Nombre de folds ADAPTATIF, comme pour le blocage spatial (voir plus haut)
  n_target <- round(length(blocks_used) * max_folds_fraction)
  n_target <- max(min_folds, min(n_target, max_folds_abs, length(blocks_used)))

  set.seed(seed)
  if (length(blocks_used) > n_target) {
    blocks_used <- sample(blocks_used, n_target)
  }
  cat(sprintf(
    "  blocs de %d j : %d blocs disponibles (>= %d obs) -> %d folds retenus (cible adaptative : %d)\n",
    block_days, length(blocks_count %>% filter(n >= min_block_n) %>% pull(temporal_block)),
    min_block_n, length(blocks_used), n_target
  ))

  folds <- map(blocks_used, ~ make_temporal_fold(data_blocked, .x, buffer_days))
  names(folds) <- paste0("t_", blocks_used)

  valid <- map_lgl(folds, ~ length(.x$train) > 0 && length(.x$test) > 0)
  if (any(!valid)) cat("  [!] ", sum(!valid), "bloc(s) temporel(s) ignore(s) : train ou test vide apres buffer\n")
  folds <- folds[valid]

  list(
    data = data_blocked, folds = folds, scheme = "temporal_block",
    block_days = block_days, buffer_days = buffer_days,
    n_blocks_total = nrow(blocks_count), n_blocks_used = length(folds),
    n_folds_target = n_target
  )
}

# ---------------------------------------------------------------------
# D. Fonction "meta" : construit TOUS les schémas demandés pour un `df`.
#    -> Appelée UNE SEULE FOIS par fréquence, puis les mêmes objets
#    `scheme` (mêmes folds) sont réutilisés pour le tuning ET
#    l'entraînement final : c'est ce qui garantit la cohérence entre
#    les deux étapes.
# ---------------------------------------------------------------------
build_all_schemes <- function(df) {
  schemes <- list()
  schemes[["naive_RS_80_20"]] <- build_naive_folds(df)

  for (r in SPATIAL_RESOLUTIONS) {
    schemes[[paste0("blocked_spatial_", r$label)]] <-
      build_spatial_folds(df, r$lon_km, r$lat_km, r$buffer_km)
  }
  for (r in TEMPORAL_RESOLUTIONS) {
    schemes[[paste0("blocked_temporal_", r$label)]] <-
      build_temporal_folds(df, r$block_days, r$buffer_days)
  }
  schemes
}
