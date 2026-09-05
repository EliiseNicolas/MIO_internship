# =====================================================================
# 15_run_transfer_tuning_test.R -- SCRIPT 5 : cout de ne PAS re-tuner
# =====================================================================
# Repond a l'axe "quelle perte de performance si on choisit une strategie
# moins couteuse" (equivalent du sondage lineaire vs reglage fin complet
# en deep learning) : applique-t-on les hyperparametres tunes en NAIVE
# TELS QUELS sur les schemas bloques (rapide, un seul tuning), ou
# re-tune-t-on separement pour chaque schema (plus couteux, cf.
# 10_run_tuning.R et sa justification) ?
#
# Pour chaque freq x modele (CART/RF/XGB) x schema BLOQUE :
#   - "sans_retuning" : hyperparametres tunes sur NAIVE, appliques direct
#     sur les folds du schema bloque
#   - "avec_retuning" : hyperparametres tunes SPECIFIQUEMENT pour ce
#     schema bloque (deja calcules par 10_run_tuning.R)
# Compare le RMSE (moyenne inter-fold) des deux approches, sur les MEMES
# folds -- seul le choix des hyperparametres differe.
#
# Sorties, sous outputs_pipeline/transfer_tuning_test/ :
#   - transfer_tuning_comparison.csv
#   - transfer_tuning_comparison.png

source("R/00_config.R")
source("R/01_data_prep.R")
source("R/02_folds.R")
source("R/03_models.R")
source("R/04_diagnostics.R")
source("R/05_plots.R")

tuning_dir <- path_out("tuning")
out_dir    <- path_out("transfer_tuning_test")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

BLOCKED_SCHEMES <- c(
  paste0("blocked_spatial_", map_chr(SPATIAL_RESOLUTIONS, "label")),
  paste0("blocked_temporal_", map_chr(TEMPORAL_RESOLUTIONS, "label"))
)
MODELS <- c("cart", "rf", "xgb")

results <- list()

for (freq in FREQS) {

  cat("\n============================================================\n")
  cat("TRANSFER TUNING TEST -- FREQUENCE :", freq, "kHz\n")
  cat("============================================================\n")

  prep_complete <- load_and_clean(freq, drop_na_numeric = TRUE)   # pour CART/RF
  prep_xgb      <- load_and_clean(freq, drop_na_numeric = FALSE)  # pour XGB
  fod_levels    <- prep_xgb$fod_levels

  schemes_rf  <- build_all_schemes(prep_complete$df)
  schemes_xgb <- build_all_schemes(prep_xgb$df)

  for (model in MODELS) {
    backend <- if (model == "xgb") make_backend("xgb", fod_levels = fod_levels) else make_backend(model)
    schemes <- if (model == "xgb") schemes_xgb else schemes_rf

    naive_tuning_path <- file.path(tuning_dir, sprintf("%s_%dkHz_naive_RS_80_20.rds", model, freq))
    if (!file.exists(naive_tuning_path)) {
      cat("  [!] tuning naive introuvable pour", model, "--", naive_tuning_path, "-- saute.\n")
      next
    }
    naive_params <- readRDS(naive_tuning_path)$best_params

    for (scheme_name in BLOCKED_SCHEMES) {

      own_tuning_path <- file.path(tuning_dir, sprintf("%s_%dkHz_%s.rds", model, freq, scheme_name))
      if (!file.exists(own_tuning_path)) {
        cat("  [!] tuning introuvable pour", model, scheme_name, "-- saute.\n")
        next
      }
      own_params <- readRDS(own_tuning_path)$best_params

      cat(" ", model, "-", scheme_name, "\n")
      scheme_blocked <- schemes[[scheme_name]]

      cv_naive_params <- run_cv_scheme(scheme_blocked, naive_params, backend, label = "sans_retuning")
      cv_own_params   <- run_cv_scheme(scheme_blocked, own_params,   backend, label = "avec_retuning")

      rmse_naive <- mean(cv_naive_params$metrics$rmse_test)
      rmse_own   <- mean(cv_own_params$metrics$rmse_test)

      results[[length(results) + 1]] <- tibble(
        freq = freq, model = model, scheme = scheme_name,
        rmse_sans_retuning = rmse_naive,
        rmse_avec_retuning = rmse_own,
        perte_relative_pct = 100 * (rmse_naive - rmse_own) / rmse_own
      )
    }
  }
}

comparison <- bind_rows(results)
write.csv(comparison, file.path(out_dir, "transfer_tuning_comparison.csv"), row.names = FALSE)

comparison_long <- comparison %>%
  pivot_longer(cols = c(rmse_sans_retuning, rmse_avec_retuning),
               names_to = "strategie", values_to = "rmse")

p <- ggplot(comparison_long, aes(x = scheme, y = rmse, fill = strategie)) +
  geom_col(position = "dodge") +
  facet_grid(model ~ freq, scales = "free_y", labeller = label_both) +
  labs(title = "Cout de ne pas re-tuner : hyperparametres naive vs re-tunes, appliques en blocage",
       subtitle = "Memes folds dans les deux cas -- seul le choix des hyperparametres differe",
       x = "Schema bloque", y = "RMSE test (moyenne inter-fold)", fill = NULL) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
ggsave(file.path(out_dir, "transfer_tuning_comparison.png"), p, width = 12, height = 8, dpi = 150)

cat("\n=== Resume (perte relative en %, positif = moins bon sans retuning) ===\n")
print(comparison %>% select(freq, model, scheme, perte_relative_pct))

cat("\nTermine. Resultats dans :", normalizePath(out_dir), "\n")
