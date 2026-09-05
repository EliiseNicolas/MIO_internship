# =====================================================================
# 16_run_cross_scheme_analysis.R -- SCRIPT 6 : diagnostics inter-schemas
# =====================================================================
# Regroupe 4 axes de comparaison au-dela du RMSE global (cf. discussion
# ViT/CNN), a partir des sorties DEJA enregistrees par 11_run_training.R
# (pas de reentrainement, juste de l'agregation + du tracer) :
#
# 1) IMPORTANCE NAIVE vs BLOQUE -- le modele regarde-t-il la bonne chose,
#    ou exploite-t-il un raccourci d'autocorrelation ? Une variable
#    dont l'importance s'effondre entre naive et blocage suggere qu'elle
#    servait surtout de proxy de proximite spatio-temporelle, pas d'un
#    vrai signal ecologique.
#
# 2) CALIBRATION -- le NASC predit suit-il la bonne moyenne par tranche
#    de valeur, pas juste correle, par modele x schema.
#
# 3) RESIDU vs LATITUDE -- diagnostic direct d'un gradient latitudinal
#    non capture par les covariables. Une courbe lissee plate -> le
#    gradient est bien absorbe. Une pente ou une forme systematique ->
#    biais residuel lie a la latitude, potentiellement pire en blocage
#    qu'en naive (la CV naive le masque partiellement).
#
# 4) RMSE vs DISTANCE test->train -- fuite spatio-temporelle : si le
#    RMSE augmente nettement avec la distance, une partie du score
#    naive est un artefact de proximite plutot qu'un vrai pouvoir
#    predictif.
#
# Sorties, sous outputs_pipeline/cross_scheme_analysis/<freq>kHz/ :
#   - importance_comparison.png
#   - calibration_<model>.png (une par modele, facette par schema)
#   - residual_vs_latitude.png
#   - rmse_vs_geo_distance.png
#   - rmse_vs_covariate_distance.png

source("R/00_config.R")
source("R/01_data_prep.R")
source("R/05_plots.R")
source("R/07_cross_scheme_utils.R")

training_dir <- path_out("training")
out_root     <- path_out("cross_scheme_analysis")
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

configs <- list_training_configs(training_dir)
cat(sprintf("Configurations trouvees : %d (freq x modele x schema)\n", nrow(configs)))

for (freq in FREQS) {

  cat("\n============================================================\n")
  cat("CROSS-SCHEME ANALYSIS -- FREQUENCE :", freq, "kHz\n")
  cat("============================================================\n")

  out_dir <- file.path(out_root, paste0(freq, "kHz"))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  configs_freq <- configs %>% filter(freq == !!freq)
  if (nrow(configs_freq) == 0) {
    cat("  [!] aucune configuration pour cette frequence -- saute.\n")
    next
  }

  # ---- 1) Importance naive vs bloque ----
  importance_all <- load_tagged_csv(configs_freq, "importance_all.csv")
  if (nrow(importance_all) > 0) {
    p_imp <- plot_importance_comparison(importance_all, subtitle = sprintf("%d kHz", freq))
    ggsave(file.path(out_dir, "importance_comparison.png"), p_imp, width = 11, height = 7, dpi = 150)
    cat("  -> importance_comparison.png\n")
  }

  # ---- 2) Calibration, une figure par modele (facette par schema) ----
  obs_pred_all <- load_tagged_csv(configs_freq, "obs_pred_all.csv")
  if (nrow(obs_pred_all) > 0) {
    for (m in unique(obs_pred_all$model)) {
      sub <- obs_pred_all %>% filter(model == m)
      calib_by_scheme <- sub %>%
        group_by(scheme) %>%
        group_modify(~ {
          binned <- .x %>% mutate(bin = dplyr::ntile(pred, 10)) %>%
            group_by(bin) %>%
            summarise(mean_obs = mean(obs), mean_pred = mean(pred), n = dplyr::n(), .groups = "drop")
          binned
        }) %>% ungroup()

      p_cal <- ggplot(calib_by_scheme, aes(x = mean_pred, y = mean_obs)) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
        geom_line(color = "steelblue", alpha = 0.6) +
        geom_point(aes(size = n), color = "steelblue") +
        facet_wrap(~scheme) +
        labs(title = sprintf("Calibration - %s - %d kHz", toupper(m), freq),
             x = "NASC predit moyen (log10), par bin", y = "NASC observe moyen (log10), par bin", size = "n obs") +
        theme_bw()
      ggsave(file.path(out_dir, sprintf("calibration_%s.png", m)), p_cal, width = 11, height = 7, dpi = 150)
      cat("  -> calibration_", m, ".png\n", sep = "")
    }

    # ---- 3) Residu vs latitude, par schema (facette par modele) ----
    # Diagnostic direct du gradient latitudinal : une courbe lissee plate
    # -> le gradient est bien absorbe par les covariables. Une pente ou
    # une forme systematique -> biais residuel lie a la latitude, non
    # capture par le modele (potentiellement pire en blocage qu'en naive).
    p_lat <- plot_residual_vs_latitude(obs_pred_all, subtitle = sprintf("%d kHz", freq))
    ggsave(file.path(out_dir, "residual_vs_latitude.png"), p_lat, width = 11, height = 7, dpi = 150)
    cat("  -> residual_vs_latitude.png\n")
  }

  # ---- 4) RMSE vs distance test->train (fuite spatio-temporelle) ----
  metrics_all <- load_tagged_csv(configs_freq, "metrics_par_fold.csv")
  if (nrow(metrics_all) > 0) {
    p_geo <- plot_rmse_vs_distance(metrics_all, "mean_geo_dist_km", subtitle = sprintf("%d kHz", freq))
    ggsave(file.path(out_dir, "rmse_vs_geo_distance.png"), p_geo, width = 9, height = 6, dpi = 150)

    p_cov <- plot_rmse_vs_distance(metrics_all, "mean_covariate_dist", subtitle = sprintf("%d kHz", freq))
    ggsave(file.path(out_dir, "rmse_vs_covariate_distance.png"), p_cov, width = 9, height = 6, dpi = 150)
    cat("  -> rmse_vs_geo_distance.png, rmse_vs_covariate_distance.png\n")
  }
}

cat("\nAnalyse inter-schemas terminee. Resultats dans :", normalizePath(out_root), "\n")
