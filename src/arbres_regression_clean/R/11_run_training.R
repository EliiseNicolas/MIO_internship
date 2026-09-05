# =====================================================================
# 11_run_training.R -- SCRIPT 2 : entraînement CART / RF / XGB
# =====================================================================
# Pour chaque fréquence (38, 120 kHz), pour CART, RF et XGB, pour :
#   - naive RS 80/20 (NAIVE_N_FOLDS tirages, fixé à N_CV = 10)
#   - blocage spatial 1500x1000km, 200x200km, 20x20km (nb de folds
#     ADAPTATIF, cf. 00_config.R / 02_folds.R -- pas forcément 10)
#   - blocage temporel 1j (nb de folds adaptatif)
# on charge les meilleurs hyperparamètres trouvés en 10_run_tuning.R et
# on entraîne un modèle par fold, en sauvegardant TOUS les diagnostics
# et plots demandés. CART et RF partagent le même schéma/folds
# (`schemes_rf`, données complètes) ; XGB utilise `schemes_xgb`
# (données avec NA préservés).
#
# STRUCTURE EN DEUX PHASES (pour des axes Y / colorbars identiques
# entre modèles/schémas/fréquences, cf. SHARED_SCALE_SCOPE dans
# 00_config.R) :
#   Phase 1 : entraîne chaque combinaison (freq x modele x schema), sauve
#     immédiatement les modèles/CSV sur disque, et ne garde en mémoire
#     que les tables LEGERES (metrics, obs_pred, importance, learning
#     curve) -- pas les modèles eux-mêmes, ni les data.frame complets des
#     schémas -- pour ne pas exploser la mémoire.
#   Phase 2 : calcule les échelles partagées (RMSE, R², variance,
#     obs/pred, résidus -- PAS l'importance, volontairement : Gain XGB et
#     importance par permutation RF ne sont pas la même unité, forcer le
#     même axe serait trompeur, pas juste esthétique -- voir README) à
#     partir de TOUTES les combinaisons collectées, puis reproduit tous
#     les plots avec ces échelles communes.
#
# Sorties, sous outputs_pipeline/training/<freq>kHz/<model>/<schema>/ :
#   - models.rds, metrics_par_fold.csv, obs_pred_all.csv, importance_all.csv
#   - learning_curve_summary.csv
#   - tous les plots (voir liste dans generate_scale_dependent_plots ci-dessous)

source("R/00_config.R")
source("R/01_data_prep.R")
source("R/02_folds.R")
source("R/03_models.R")
source("R/04_diagnostics.R")
source("R/05_plots.R")

training_dir <- path_out("training")
tuning_dir   <- path_out("tuning")
dir.create(training_dir, showWarnings = FALSE, recursive = TRUE)

MODEL_SCHEMES <- c("naive_RS_80_20",
                   paste0("blocked_spatial_", map_chr(SPATIAL_RESOLUTIONS, "label")),
                   paste0("blocked_temporal_", map_chr(TEMPORAL_RESOLUTIONS, "label")))
MODELS <- c("cart", "rf", "xgb")

# Si TRUE (defaut) : une combinaison (freq/modele/schema) deja entrainee
# avec succes lors d'un run precedent (tous ses fichiers de sortie
# presents sur disque) est RECHARGEE au lieu d'etre reentrainee -- utile
# pour reprendre apres une interruption/erreur sans perdre le travail
# deja fait. Mettre a FALSE pour forcer un reentrainement complet (ex.
# apres avoir change les hyperparametres tunes).
SKIP_EXISTING_TRAINING <- TRUE

REQUIRED_OUTPUT_FILES <- c(
  "models.rds", "metrics_par_fold.csv", "obs_pred_all.csv",
  "importance_all.csv", "learning_curve_summary.csv"
)

# ---------------------------------------------------------------------
# PHASE 1 : entraînement + sauvegarde immédiate des objets lourds +
# collecte des tables légères pour la phase 2
# ---------------------------------------------------------------------
all_runs <- list()

for (freq in FREQS) {

  cat("\n============================================================\n")
  cat("ENTRAINEMENT (PHASE 1) -- FREQUENCE :", freq, "kHz\n")
  cat("============================================================\n")

  prep_complete <- load_and_clean(freq, drop_na_numeric = TRUE)   # pour CART/RF
  prep_xgb      <- load_and_clean(freq, drop_na_numeric = FALSE)  # pour XGB
  fod_levels    <- prep_xgb$fod_levels

  schemes_rf  <- build_all_schemes(prep_complete$df)
  schemes_xgb <- build_all_schemes(prep_xgb$df)

  for (scheme_name in MODEL_SCHEMES) {
    for (model in MODELS) {

      out_dir <- file.path(training_dir, paste0(freq, "kHz"), model, scheme_name)

      # ---- Deja fait ? Recharge depuis le disque, ne reentraine pas ----
      already_done <- SKIP_EXISTING_TRAINING &&
        all(file.exists(file.path(out_dir, REQUIRED_OUTPUT_FILES)))

      if (already_done) {
        cat("  [deja fait]", model, "-", scheme_name, "-- rechargement depuis disque\n")
        metrics_df    <- read.csv(file.path(out_dir, "metrics_par_fold.csv"))
        obs_pred_df   <- read.csv(file.path(out_dir, "obs_pred_all.csv"))
        importance_df <- read.csv(file.path(out_dir, "importance_all.csv"))
        lc_summary_df <- read.csv(file.path(out_dir, "learning_curve_summary.csv"))

        all_runs[[length(all_runs) + 1]] <- list(
          freq = freq, model = model, scheme = scheme_name, out_dir = out_dir,
          n_folds = nrow(metrics_df),
          metrics = metrics_df, obs_pred = obs_pred_df,
          importance = importance_df, lc_summary = lc_summary_df
        )
        next
      }

      cat(" ", model, "-", scheme_name, "\n")

      tuning_path <- file.path(tuning_dir, sprintf("%s_%dkHz_%s.rds", model, freq, scheme_name))
      if (!file.exists(tuning_path)) {
        cat("  [!] tuning introuvable --", tuning_path, "-- saute.\n")
        next
      }
      params <- readRDS(tuning_path)$best_params

      if (model == "xgb") {
        backend <- make_backend("xgb", fod_levels = fod_levels)
        scheme  <- schemes_xgb[[scheme_name]]
      } else {
        # CART et RF partagent le meme schema/folds (donnees completes)
        backend <- make_backend(model)
        scheme  <- schemes_rf[[scheme_name]]
      }

      cv_res <- run_cv_scheme(scheme, params, backend, label = scheme_name)
      lc     <- compute_learning_curve(scheme, params, backend)

      dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

      # -- sauvegarde immediate des objets lourds (modeles) --
      saveRDS(cv_res$models, file.path(out_dir, "models.rds"))
      write.csv(cv_res$metrics,     file.path(out_dir, "metrics_par_fold.csv"), row.names = FALSE)
      write.csv(cv_res$obs_pred,    file.path(out_dir, "obs_pred_all.csv"), row.names = FALSE)
      write.csv(cv_res$importance,  file.path(out_dir, "importance_all.csv"), row.names = FALSE)
      write.csv(lc$summary,         file.path(out_dir, "learning_curve_summary.csv"), row.names = FALSE)

      # -- plots qui n'ont PAS besoin d'echelle partagee : generes tout de
      #    suite (cartes de reference spatiale, categorique) --
      p_points <- plot_points_colored_by_fold(scheme, subtitle = sprintf("%s - %d kHz - %s", toupper(model), freq, scheme_name))
      ggsave(file.path(out_dir, "16_carte_points_par_fold.png"), p_points, width = 8, height = 6, dpi = 150)
      if (scheme$scheme != "naive_RS_80_20") {
        p_buffer <- plot_spatial_train_test(scheme, subtitle = sprintf("%s - %d kHz - %s", toupper(model), freq, scheme_name))
        ggsave(file.path(out_dir, "17_carte_train_test_buffer.png"), p_buffer, width = 10, height = 8, dpi = 150)
      }
      # -- plots de metadonnees (distances/distributions) : pas prioritaires
      #    pour la comparabilite inter-config (chaque config est inspectee
      #    individuellement), generes tout de suite aussi --
      p_dist   <- plot_distances_per_fold(cv_res$metrics, subtitle = scheme_name)
      ggsave(file.path(out_dir, "11_distances_test_train.png"), p_dist, width = 9, height = 5, dpi = 150)
      p_cov    <- plot_covariate_stats(cv_res$covariate_stats, subtitle = scheme_name)
      ggsave(file.path(out_dir, "12_covariable_stats.png"), p_cov, width = 10, height = 8, dpi = 150)
      p_covvar <- plot_covariate_variance(cv_res$covariate_stats, subtitle = scheme_name)
      ggsave(file.path(out_dir, "13_covariable_variance.png"), p_covvar, width = 10, height = 8, dpi = 150)
      p_fod    <- plot_fod_distribution(cv_res$fod_dist, subtitle = scheme_name)
      ggsave(file.path(out_dir, "14_fod_distribution.png"), p_fod, width = 9, height = 7, dpi = 150)
      p_distr  <- plot_distributions_per_fold(cv_res$numeric_dist, subtitle = scheme_name)
      ggsave(file.path(out_dir, "15_distributions_par_fold.png"), p_distr,
             width = 12, height = max(6, 1.2 * length(scheme$folds)), dpi = 150)

      # -- collecte LEGERE pour la phase 2 (echelles partagees) --
      all_runs[[length(all_runs) + 1]] <- list(
        freq = freq, model = model, scheme = scheme_name, out_dir = out_dir,
        n_folds = length(scheme$folds),
        metrics    = cv_res$metrics,
        obs_pred   = cv_res$obs_pred,
        importance = cv_res$importance,
        lc_summary = lc$summary
      )
    }
  }
}

cat(sprintf("\nPhase 1 terminee : %d configurations entrainees.\n", length(all_runs)))

# ---------------------------------------------------------------------
# PHASE 2 : calcul des echelles partagees, puis generation de tous les
# plots "sensibles a l'echelle" avec ces bornes communes
# ---------------------------------------------------------------------
cat("\n============================================================\n")
cat("GENERATION DES PLOTS (PHASE 2) -- echelles partagees\n")
cat("============================================================\n")

scope_key <- function(run) if (SHARED_SCALE_SCOPE == "per_freq") as.character(run$freq) else "ALL"

# -- regroupe les runs par cle de portee (soit tout ensemble, soit par freq) --
scope_keys <- unique(vapply(all_runs, scope_key, character(1)))

for (sk in scope_keys) {

  runs_scope <- Filter(function(r) scope_key(r) == sk, all_runs)
  scope_label <- if (SHARED_SCALE_SCOPE == "per_freq") sprintf(" (%s kHz)", sk) else " (toutes frequences)"

  metrics_all_scope    <- bind_rows(lapply(runs_scope, `[[`, "metrics"))
  obs_pred_all_scope   <- bind_rows(lapply(runs_scope, `[[`, "obs_pred"))
  lc_all_scope         <- bind_rows(lapply(runs_scope, `[[`, "lc_summary"))

  # -- echelles partagees (memes unites -> comparaison directe legitime) --
  range_rmse <- range(c(metrics_all_scope$rmse_test, metrics_all_scope$rmse_train,
                         lc_all_scope$mean_rmse_train + lc_all_scope$sd_rmse_train,
                         lc_all_scope$mean_rmse_test  + lc_all_scope$sd_rmse_test,
                         lc_all_scope$mean_rmse_train - lc_all_scope$sd_rmse_train,
                         lc_all_scope$mean_rmse_test  - lc_all_scope$sd_rmse_test),
                       na.rm = TRUE)
  range_rmse[1] <- max(0, range_rmse[1])  # RMSE >= 0

  range_r2   <- range(metrics_all_scope$r2_test, na.rm = TRUE)
  range_var  <- range(metrics_all_scope$var_intra_fold_test, na.rm = TRUE)

  range_obs_pred <- range(c(obs_pred_all_scope$obs, obs_pred_all_scope$pred), na.rm = TRUE)
  range_residual <- range(obs_pred_all_scope$residual, na.rm = TRUE)
  range_residual <- c(-max(abs(range_residual)), max(abs(range_residual)))  # symetrique autour de 0
  range_abs_residual <- c(0, max(abs(obs_pred_all_scope$residual), na.rm = TRUE))

  cat(sprintf("Echelles partagees%s : RMSE=[%.3f,%.3f] R2=[%.3f,%.3f] obs/pred=[%.3f,%.3f] residu=[%.3f,%.3f]\n",
              scope_label, range_rmse[1], range_rmse[2], range_r2[1], range_r2[2],
              range_obs_pred[1], range_obs_pred[2], range_residual[1], range_residual[2]))

  # -- Importance : PARTAGEE PAR MODELE (memes unites au sein d'un
  #    modele donne, entre schemas/frequences), PAS ENTRE MODELES
  #    differents (Gain XGB != importance par permutation RF/CART) --
  range_importance_by_model <- runs_scope %>%
    purrr::map_dfr(~ mutate(.x$importance, model = .x$model)) %>%
    group_by(model) %>%
    summarise(min_imp = min(importance, na.rm = TRUE), max_imp = max(importance, na.rm = TRUE), .groups = "drop")

  for (run in runs_scope) {

    label <- sprintf("%s - %d kHz - %s", toupper(run$model), run$freq, run$scheme)
    save_plot <- function(p, name, w = 8, h = 6) {
      ggsave(file.path(run$out_dir, name), p, width = w, height = h, dpi = 150)
    }

    save_plot(plot_learning_curve(run$lc_summary, subtitle = label, ylim = range_rmse), "01_learning_curve.png")
    save_plot(plot_obs_vs_pred_scatter(run$obs_pred, subtitle = label, axis_limits = range_obs_pred),
              "02_obs_vs_pred_scatter.png")
    save_plot(plot_obs_vs_pred_hist(run$obs_pred, subtitle = label), "03_obs_vs_pred_hist.png")
    save_plot(plot_obs_vs_pred_by_fold(run$obs_pred, subtitle = label, axis_limits = range_obs_pred),
              "04_obs_vs_pred_by_fold.png", 10, 8)
    save_plot(plot_metric_bar(run$metrics, "rmse_test", subtitle = label, ylim = range_rmse), "05_rmse_par_fold.png")
    save_plot(plot_metric_bar(run$metrics, "r2_test", fill = "darkgreen", hline0 = TRUE, subtitle = label, ylim = range_r2),
              "06_r2_par_fold.png")

    imp_range_model <- range_importance_by_model %>% filter(model == run$model)
    ylim_imp <- if (nrow(imp_range_model) > 0) c(imp_range_model$min_imp, imp_range_model$max_imp) else NULL
    save_plot(plot_importance_mean_sd(run$importance, subtitle = label, ylim = ylim_imp), "07_importance_variables.png")

    save_plot(plot_residual_map(run$obs_pred, subtitle = label, color_limits = range_residual), "08_carte_residus.png")
    save_plot(plot_abs_residual_map(run$obs_pred, subtitle = label, color_limits = range_abs_residual),
              "09_carte_residus_abs.png")
    save_plot(plot_metric_bar(run$metrics, "var_intra_fold_test", fill = "orange", subtitle = label, ylim = range_var),
              "10_variance_intra_fold.png")

    global_rmse <- rmse_fn(run$obs_pred$obs, run$obs_pred$pred)
    global_r2   <- 1 - sum((run$obs_pred$obs - run$obs_pred$pred)^2) /
      sum((run$obs_pred$obs - mean(run$obs_pred$obs))^2)
    writeLines(
      sprintf("RMSE globale (out-of-fold) = %.4f\nR2 global (out-of-fold) = %.4f\nNombre de folds = %d",
              global_rmse, global_r2, run$n_folds),
      file.path(run$out_dir, "resume_global.txt")
    )
    cat(sprintf("  -> %s : RMSE=%.4f R2=%.4f (%d folds)\n", label, global_rmse, global_r2, run$n_folds))
  }
}

cat("\nEntrainement termine. Resultats dans :", normalizePath(training_dir), "\n")
