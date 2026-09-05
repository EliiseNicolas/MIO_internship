# =====================================================================
# 18_run_noise_robustness_test.R -- SCRIPT 8 : robustesse au bruit
# =====================================================================
# Repond a l'axe "robustesse a la perturbation" (equivalent tabulaire de
# la robustesse a la perturbation de texture en vision, cf. discussion
# ViT/CNN) : pour chaque freq x modele (CART/RF/XGB) x schema, on
# entraine UNE FOIS par fold (hyperparametres deja tunes par
# 10_run_tuning.R) puis on reevalue sur des copies de plus en plus
# bruitees du meme jeu de test (bruit gaussien sur les covariables
# numeriques, amplitude exprimee en fraction de l'ecart-type de chaque
# covariable -- voir 08_robustness.R pour le detail).
#
# Question posee : quel modele se degrade le moins vite quand les
# covariables sont perturbees ? Et est-ce coherent entre naive et
# blocage (un modele robuste au bruit l'est-il aussi au changement de
# domaine spatial/temporel) ?
#
# ECHELLE Y PARTAGEE (cf. SHARED_SCALE_SCOPE, 00_config.R) : le RMSE est
# la MEME unite pour CART/RF/XGB (contrairement a l'importance des
# variables) -- partager l'axe Y entre modeles est donc parfaitement
# legitime ici. Structure en deux phases : Phase 1 calcule tout (par
# freq) et sauve les CSV ; Phase 2 calcule les echelles partagees a
# partir de TOUTES les frequences/modeles/schemas puis trace.
#
# Sorties, sous outputs_pipeline/noise_robustness/<freq>kHz/ :
#   - noise_robustness_raw.csv (toutes les combinaisons fold x niveau x modele x schema)
#   - noise_robustness_<model>.png (RMSE vs bruit, facette par schema)
#   - noise_robustness_relative_<model>.png (degradation % vs bruit)
#   - noise_robustness_all_models.png (comparaison directe CART/RF/XGB, schema naive)

source("R/00_config.R")
source("R/01_data_prep.R")
source("R/02_folds.R")
source("R/03_models.R")
source("R/04_diagnostics.R")
source("R/05_plots.R")
source("R/08_robustness.R")

tuning_dir <- path_out("tuning")
out_root   <- path_out("noise_robustness")
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

NOISE_LEVELS <- c(0, 0.1, 0.25, 0.5, 1, 2)
MODELS <- c("cart", "rf", "xgb")

ALL_SCHEMES <- c(
  "naive_RS_80_20",
  paste0("blocked_spatial_", map_chr(SPATIAL_RESOLUTIONS, "label")),
  paste0("blocked_temporal_", map_chr(TEMPORAL_RESOLUTIONS, "label"))
)

# ---------------------------------------------------------------------
# PHASE 1 : calcule tout, par frequence, sauve les CSV
# ---------------------------------------------------------------------
noise_by_freq <- list()

for (freq in FREQS) {

  cat("\n============================================================\n")
  cat("NOISE ROBUSTNESS TEST (PHASE 1) -- FREQUENCE :", freq, "kHz\n")
  cat("============================================================\n")

  out_dir <- file.path(out_root, paste0(freq, "kHz"))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  prep_complete <- load_and_clean(freq, drop_na_numeric = TRUE)   # pour CART/RF
  prep_xgb      <- load_and_clean(freq, drop_na_numeric = FALSE)  # pour XGB
  fod_levels    <- prep_xgb$fod_levels

  schemes_rf  <- build_all_schemes(prep_complete$df)
  schemes_xgb <- build_all_schemes(prep_xgb$df)

  results <- list()

  for (model in MODELS) {
    backend <- if (model == "xgb") make_backend("xgb", fod_levels = fod_levels) else make_backend(model)
    schemes <- if (model == "xgb") schemes_xgb else schemes_rf

    for (scheme_name in ALL_SCHEMES) {

      tuning_path <- file.path(tuning_dir, sprintf("%s_%dkHz_%s.rds", model, freq, scheme_name))
      if (!file.exists(tuning_path)) {
        cat("  [!] tuning introuvable pour", model, scheme_name, "-- saute.\n")
        next
      }
      params <- readRDS(tuning_path)$best_params

      cat(" ", model, "-", scheme_name, "\n")
      res <- run_noise_robustness_scheme(
        schemes[[scheme_name]], params, backend,
        noise_levels = NOISE_LEVELS, label = paste(model, scheme_name)
      )
      results[[length(results) + 1]] <- res %>% mutate(model = model, scheme = scheme_name)
    }
  }

  noise_all <- bind_rows(results)
  write.csv(noise_all, file.path(out_dir, "noise_robustness_raw.csv"), row.names = FALSE)
  noise_by_freq[[as.character(freq)]] <- list(out_dir = out_dir, noise_all = noise_all)
}

cat("\nPhase 1 terminee.\n")

# ---------------------------------------------------------------------
# PHASE 2 : echelles partagees + generation des plots
# ---------------------------------------------------------------------
cat("\n============================================================\n")
cat("GENERATION DES PLOTS (PHASE 2) -- echelles partagees\n")
cat("============================================================\n")

scope_key <- function(f) if (SHARED_SCALE_SCOPE == "per_freq") as.character(f) else "ALL"
scope_keys <- unique(vapply(FREQS, scope_key, character(1)))

for (sk in scope_keys) {

  freqs_in_scope <- FREQS[vapply(FREQS, function(f) scope_key(f) == sk, logical(1))]
  noise_scope <- bind_rows(lapply(noise_by_freq[as.character(freqs_in_scope)], `[[`, "noise_all"))

  range_rmse_abs <- range(noise_scope$rmse_test, na.rm = TRUE)
  range_rmse_abs[1] <- max(0, range_rmse_abs[1])

  degrad_scope <- noise_scope %>%
    group_by(model, scheme, noise_level) %>%
    summarise(mean_rmse = mean(rmse_test), .groups = "drop") %>%
    group_by(model, scheme) %>%
    mutate(rmse_baseline = mean_rmse[noise_level == min(noise_level)],
           degradation_pct = 100 * (mean_rmse - rmse_baseline) / rmse_baseline) %>%
    ungroup()
  range_degrad <- range(degrad_scope$degradation_pct, na.rm = TRUE)
  range_degrad[1] <- min(0, range_degrad[1])  # inclut toujours 0 (baseline)

  cat(sprintf("Echelles partagees (%s) : RMSE=[%.3f,%.3f] degradation%%=[%.1f,%.1f]\n",
              sk, range_rmse_abs[1], range_rmse_abs[2], range_degrad[1], range_degrad[2]))

  for (freq in freqs_in_scope) {
    entry <- noise_by_freq[[as.character(freq)]]
    out_dir <- entry$out_dir
    noise_all <- entry$noise_all

    for (m in unique(noise_all$model)) {
      sub <- noise_all %>% filter(model == m)

      p_abs <- plot_noise_robustness(sub, group_col = "scheme",
                                      subtitle = sprintf("%s - %d kHz", toupper(m), freq), ylim = range_rmse_abs)
      ggsave(file.path(out_dir, sprintf("noise_robustness_%s.png", m)), p_abs, width = 9, height = 6, dpi = 150)

      p_rel <- plot_noise_robustness_relative(sub, group_col = "scheme",
                                               subtitle = sprintf("%s - %d kHz", toupper(m), freq), ylim = range_degrad)
      ggsave(file.path(out_dir, sprintf("noise_robustness_relative_%s.png", m)), p_rel, width = 9, height = 6, dpi = 150)
    }

    naive_only <- noise_all %>% filter(scheme == "naive_RS_80_20")
    if (nrow(naive_only) > 0) {
      p_compare <- plot_noise_robustness(naive_only, group_col = "model",
                                          subtitle = sprintf("Schema naive - %d kHz", freq), ylim = range_rmse_abs)
      ggsave(file.path(out_dir, "noise_robustness_all_models.png"), p_compare, width = 8, height = 6, dpi = 150)
    }

    cat("  -> resultats enregistres dans", out_dir, "\n")
  }
}

cat("\nTest de robustesse au bruit termine. Resultats dans :", normalizePath(out_root), "\n")
