# =====================================================================
# 17_run_map_comparison.R -- SCRIPT 7 : comparaison des cartes produites
# =====================================================================
# Compare quantitativement toutes les cartes de prediction produites
# pour TARGET_DATE_SINGLE par 12_run_grid_prediction.R (RF/XGB, naive +
# schemas bloques, avec/sans NA) et 14_run_rfsrc_reconstruction.R
# (randomForestSRC, imputation complete) -- pas juste visuellement, mais
# via : correlation entre cartes, RMSE entre paires de cartes (sur les
# pixels communs), et cartes de difference pour TOUTES les paires
# possibles (N*(N-1)/2 avec N couches -- peut representer plusieurs
# dizaines/centaines de fichiers, ranges dans un sous-dossier all_pairs/).
#
# PREREQUIS : avoir lance 12_run_grid_prediction.R ET
# 14_run_rfsrc_reconstruction.R (ce script lit leurs sorties .rds, ne
# refait AUCUNE prediction lui-meme).
#
# Sorties, sous outputs_pipeline/map_comparison/<freq>kHz/ :
#   - correlation_matrix.csv, correlation_heatmap.png
#   - pairwise_rmse.csv (RMSE + n pixels communs, TOUTES les paires)
#   - all_pairs/diff_<layer1>_MOINS_<layer2>.png (TOUTES les combinaisons)

source("R/00_config.R")
source("R/01_data_prep.R")
source("R/05_plots.R")

prediction_dir <- path_out("predictions")
rfsrc_dir      <- path_out("rfsrc_reconstruction")
out_root       <- path_out("map_comparison")
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------
# Chargement de toutes les cartes deja produites (aucune reprediction)
# ---------------------------------------------------------------------
predictions_path <- file.path(prediction_dir, "predictions_all_combined.rds")
if (!file.exists(predictions_path)) {
  stop("Fichier introuvable : ", predictions_path,
       " -- lance d'abord 12_run_grid_prediction.R.")
}
predictions_all <- readRDS(predictions_path)

rfsrc_layers <- purrr::map_dfr(FREQS, function(f) {
  p <- file.path(rfsrc_dir, paste0(f, "kHz"), "reconstruction_grid.rds")
  if (!file.exists(p)) {
    cat("  [!] reconstruction rfsrc introuvable pour", f, "kHz -- ignoree (lance 14_run_rfsrc_reconstruction.R).\n")
    return(NULL)
  }
  readRDS(p)
})

all_layers <- bind_rows(predictions_all, rfsrc_layers)
cat(sprintf("Cartes chargees : %d lignes, %d couches (layer_id) distinctes\n",
            nrow(all_layers), length(unique(all_layers$layer_id))))

# ---------------------------------------------------------------------
# Fonctions de comparaison
# ---------------------------------------------------------------------
pairwise_rmse_table <- function(wide_df, layer_cols) {
  pairs <- combn(layer_cols, 2, simplify = FALSE)
  purrr::map_dfr(pairs, function(p) {
    sub <- wide_df[, p]
    complete <- stats::complete.cases(sub)
    n_common <- sum(complete)
    if (n_common == 0) {
      return(tibble(layer_1 = p[1], layer_2 = p[2], n_common_pixels = 0,
                     rmse = NA_real_, correlation = NA_real_))
    }
    diff <- sub[complete, 1] - sub[complete, 2]
    tibble(
      layer_1 = p[1], layer_2 = p[2],
      n_common_pixels = n_common,
      rmse = sqrt(mean(diff^2)),
      correlation = cor(sub[complete, 1], sub[complete, 2])
    )
  })
}

for (freq in FREQS) {

  cat("\n============================================================\n")
  cat("MAP COMPARISON -- FREQUENCE :", freq, "kHz\n")
  cat("============================================================\n")

  out_dir <- file.path(out_root, paste0(freq, "kHz"))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  layers_freq <- all_layers %>% filter(freq == !!freq)
  if (nrow(layers_freq) == 0) {
    cat("  [!] aucune carte pour cette frequence -- saute.\n")
    next
  }

  layer_ids <- unique(layers_freq$layer_id)
  cat(sprintf("  %d couches : %s\n", length(layer_ids), paste(layer_ids, collapse = ", ")))

  # ---- Passage en format large : une colonne par couche, jointe sur lon/lat ----
  wide_df <- layers_freq %>%
    select(lon, lat, layer_id, NASC_pred) %>%
    distinct(lon, lat, layer_id, .keep_all = TRUE) %>%
    pivot_wider(names_from = layer_id, values_from = NASC_pred)

  layer_cols <- setdiff(names(wide_df), c("lon", "lat"))
  if (length(layer_cols) < 2) {
    cat("  [!] moins de 2 couches disponibles -- pas de comparaison possible.\n")
    next
  }

  # ---- Correlation entre TOUTES les paires de cartes ----
  cor_matrix <- cor(wide_df[, layer_cols], use = "pairwise.complete.obs")
  write.csv(cor_matrix, file.path(out_dir, "correlation_matrix.csv"))

  p_cor <- plot_correlation_heatmap(cor_matrix, subtitle = sprintf("%d kHz - %s", freq, format(TARGET_DATE_SINGLE)))
  ggsave(file.path(out_dir, "correlation_heatmap.png"), p_cor,
         width = max(6, length(layer_cols) * 0.9), height = max(5, length(layer_cols) * 0.8), dpi = 150)

  # ---- RMSE entre TOUTES les paires (sur pixels communs) ----
  rmse_table <- pairwise_rmse_table(as.data.frame(wide_df[, layer_cols]), layer_cols)
  write.csv(rmse_table, file.path(out_dir, "pairwise_rmse.csv"), row.names = FALSE)
  cat("  -> correlation_matrix.csv, correlation_heatmap.png, pairwise_rmse.csv\n")

  # ---- Cartes de difference pour TOUTES les combinaisons possibles ----
  # Echelle symetrique PARTAGEE entre TOUTES les cartes de difference
  # (calculee a partir de toutes les paires en meme temps, avant de
  # tracer) -- pour qu'une meme intensite de rouge/bleu represente le
  # meme ecart quelle que soit la paire affichee. ATTENTION AU VOLUME :
  # avec N couches, ca fait N*(N-1)/2 cartes (ex. 16 couches -> 120
  # cartes) -- rangees dans un sous-dossier all_pairs/ pour ne pas
  # polluer le dossier principal.
  diff_dir <- file.path(out_dir, "all_pairs")
  dir.create(diff_dir, showWarnings = FALSE, recursive = TRUE)

  all_pairs <- combn(layer_cols, 2, simplify = FALSE)
  cat(sprintf("  %d paires possibles (%d couches) -> cartes de difference dans all_pairs/\n",
              length(all_pairs), length(layer_cols)))

  # -- Phase A : calcule tous les diffs d'abord (pour l'echelle commune) --
  diffs_computed <- purrr::map(all_pairs, function(p) {
    list(
      a = p[1], b = p[2],
      diff_df = wide_df %>% transmute(lon, lat, diff = .data[[p[1]]] - .data[[p[2]]])
    )
  })

  all_diff_vals <- unlist(lapply(diffs_computed, function(x) x$diff_df$diff))
  shared_diff_limit <- max(abs(all_diff_vals), na.rm = TRUE)
  shared_diff_limits <- c(-shared_diff_limit, shared_diff_limit)

  # -- Phase B : trace avec l'echelle commune --
  for (d in diffs_computed) {
    pair_label <- paste0(d$a, "_MOINS_", d$b)
    p_diff <- plot_map_difference(
      d$diff_df, diff_col = "diff",
      title = paste0("Difference : ", d$a, " moins ", d$b),
      subtitle = sprintf("%d kHz - %s", freq, format(TARGET_DATE_SINGLE)),
      limits = shared_diff_limits
    )
    ggsave(file.path(diff_dir, paste0("diff_", pair_label, ".png")), p_diff, width = 8, height = 6, dpi = 150)
  }
  cat(sprintf("  -> %d cartes de difference enregistrees dans %s\n", length(diffs_computed), diff_dir))
}

cat("\nComparaison de cartes terminee. Resultats dans :", normalizePath(out_root), "\n")
