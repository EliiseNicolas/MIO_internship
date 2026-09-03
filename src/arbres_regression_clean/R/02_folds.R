# =====================================================================
# 02_folds.R -- construction des schémas de validation croisée
# =====================================================================
# Toutes les fonctions prennent le data.frame `df` en paramètre explicite
# (pas de variable globale) afin d'être réutilisables telles quelles pour
# le tuning ET pour l'entraînement final -- avec les MÊMES folds.

rmse_fn <- function(obs, pred) sqrt(mean((obs - pred)^2, na.rm = TRUE))

# ---------------------------------------------------------------------
# A. NAIVE RANDOM SPLIT (RS) 80/20, répété N_CV fois (Monte-Carlo CV)
#    -- c'est la définition standard d'un "RS 80/20" répété : à chaque
#    itération on retire un nouveau tirage aléatoire de 80% du jeu de
#    données pour l'entraînement. Les folds ne sont PAS une partition
#    (contrairement à une k-fold CV classique).
# ---------------------------------------------------------------------
build_naive_folds <- function(df, n_folds = NAIVE_N_FOLDS, frac_train = 0.8, seed = 123) {
  set.seed(seed)
  n <- nrow(df)
  folds <- map(seq_len(n_folds), function(k) {
    train_idx <- sample(seq_len(n), size = round(frac_train * n))
    test_idx  <- setdiff(seq_len(n), train_idx)
    list(train = train_idx, test = test_idx)
  })
  names(folds) <- paste0("rs_", seq_len(n_folds))
  list(data = df, folds = folds, scheme = "naive_RS_80_20")
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
    blocks_used <- sample(blocks_used, n_target)
  }
  cat(sprintf(
    "  %sx%skm : %d blocs disponibles (>= %d obs) -> %d folds retenus (cible adaptative : %d)\n",
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
