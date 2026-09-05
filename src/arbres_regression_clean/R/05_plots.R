# =====================================================================
# 05_plots.R -- toutes les fonctions de tracé, indépendantes du modèle
# =====================================================================
# Convention : chaque fonction retourne un objet ggplot (ou patchwork),
# à sauvegarder ensuite avec ggsave(). Aucune fonction n'appelle print()
# ni ggsave() elle-même -- ça reste au script appelant.

theme_pipeline <- theme_bw()

# ---------------------------------------------------------------------
# 1. Courbe d'apprentissage
# ---------------------------------------------------------------------
plot_learning_curve <- function(lc_summary, subtitle = "", ylim = NULL) {
  p <- ggplot(lc_summary, aes(x = fraction)) +
    geom_ribbon(aes(ymin = mean_rmse_train - sd_rmse_train, ymax = mean_rmse_train + sd_rmse_train),
                fill = "steelblue", alpha = 0.15) +
    geom_ribbon(aes(ymin = mean_rmse_test - sd_rmse_test, ymax = mean_rmse_test + sd_rmse_test),
                fill = "firebrick", alpha = 0.15) +
    geom_line(aes(y = mean_rmse_train, color = "Train")) +
    geom_point(aes(y = mean_rmse_train, color = "Train")) +
    geom_line(aes(y = mean_rmse_test, color = "Test")) +
    geom_point(aes(y = mean_rmse_test, color = "Test")) +
    scale_color_manual(values = c(Train = "steelblue", Test = "firebrick")) +
    labs(title = "Courbe d'apprentissage", subtitle = subtitle,
         x = "Fraction du train utilisée", y = "RMSE", color = NULL) +
    theme_pipeline
  if (!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim)
  p
}

# ---------------------------------------------------------------------
# 2. Courbe de validation (RMSE en fonction d'un hyperparamètre) --
#    utile pendant le tuning pour visualiser le compromis biais/variance
# ---------------------------------------------------------------------
plot_validation_curve <- function(tuning_summary, x_var, facet_var = NULL, subtitle = "") {
  p <- ggplot(tuning_summary, aes(x = .data[[x_var]], y = mean_rmse_test)) +
    geom_line(aes(group = 1), color = "steelblue") +
    geom_point() +
    geom_errorbar(aes(ymin = mean_rmse_test - sd_rmse_test, ymax = mean_rmse_test + sd_rmse_test), width = 0.1) +
    labs(title = paste0("Tuning : RMSE test vs ", x_var), subtitle = subtitle,
         x = x_var, y = "RMSE test (moyenne inter-fold)") +
    theme_pipeline
  if (!is.null(facet_var)) p <- p + facet_wrap(vars(.data[[facet_var]]), labeller = label_both)
  p
}

# ---------------------------------------------------------------------
# 3. Observé vs prédit : nuage de points (coloré par fold) + histogramme
#    des résidus, et version facettée par fold
# ---------------------------------------------------------------------
plot_obs_vs_pred_scatter <- function(obs_pred_all, subtitle = "", axis_limits = NULL) {
  rmse_g <- rmse_fn(obs_pred_all$obs, obs_pred_all$pred)
  r2_g   <- 1 - sum((obs_pred_all$obs - obs_pred_all$pred)^2) /
    sum((obs_pred_all$obs - mean(obs_pred_all$obs))^2)

  if (is.null(axis_limits)) axis_limits <- range(c(obs_pred_all$obs, obs_pred_all$pred), na.rm = TRUE)

  ggplot(obs_pred_all, aes(x = obs, y = pred, color = factor(fold_id))) +
    geom_point(alpha = 0.4, size = 1.2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
    coord_cartesian(xlim = axis_limits, ylim = axis_limits) +
    annotate("label", x = axis_limits[1], y = axis_limits[2],
             hjust = 0, vjust = 1, label = sprintf("RMSE = %.3f\nR\u00b2 = %.3f", rmse_g, r2_g),
             fill = "white", alpha = 0.85, size = 3.5) +
    labs(title = "NASC observe vs predit (out-of-fold, tous folds poolés)", subtitle = subtitle,
         x = "NASC observe (log10)", y = "NASC predit (log10)", color = "Fold") +
    theme_pipeline
}

plot_obs_vs_pred_hist <- function(obs_pred_all, subtitle = "") {
  long <- bind_rows(
    tibble(value = obs_pred_all$obs,  type = "Observe"),
    tibble(value = obs_pred_all$pred, type = "Predit")
  )
  ggplot(long, aes(x = value, fill = type)) +
    geom_histogram(alpha = 0.5, position = "identity", bins = 40) +
    labs(title = "Distribution NASC observe vs predit", subtitle = subtitle,
         x = "NASC (log10)", y = "Effectif", fill = NULL) +
    theme_pipeline
}

plot_obs_vs_pred_by_fold <- function(obs_pred_all, subtitle = "", axis_limits = NULL) {
  p <- ggplot(obs_pred_all, aes(x = obs, y = pred)) +
    geom_point(alpha = 0.4, size = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    facet_wrap(~fold_id) +
    labs(title = "NASC observe vs predit, par fold", subtitle = subtitle,
         x = "NASC observe (log10)", y = "NASC predit (log10)") +
    theme_pipeline
  if (!is.null(axis_limits)) p <- p + coord_cartesian(xlim = axis_limits, ylim = axis_limits)
  p
}

# ---------------------------------------------------------------------
# 4. Importance des variables : moyenne +/- écart-type inter-fold
# ---------------------------------------------------------------------
plot_importance_mean_sd <- function(importance_all, subtitle = "", ylim = NULL) {
  summary_imp <- importance_all %>%
    group_by(variable) %>%
    summarise(mean_imp = mean(importance, na.rm = TRUE),
              sd_imp   = sd(importance, na.rm = TRUE), .groups = "drop")

  p <- ggplot(summary_imp, aes(x = reorder(variable, mean_imp), y = mean_imp)) +
    geom_col(fill = "darkorange") +
    geom_errorbar(aes(ymin = pmax(mean_imp - sd_imp, 0), ymax = mean_imp + sd_imp), width = 0.2) +
    labs(title = "Importance moyenne des variables (+/- ecart-type inter-fold)", subtitle = subtitle,
         x = NULL, y = "Importance") +
    theme_pipeline
  if (!is.null(ylim)) p <- p + coord_flip(ylim = ylim) else p <- p + coord_flip()
  p
}

# ---------------------------------------------------------------------
# 5. Carte des résidus (out-of-fold)
# ---------------------------------------------------------------------
plot_residual_map <- function(obs_pred_all, subtitle = "", color_limits = NULL) {
  ggplot(obs_pred_all, aes(x = lon, y = lat, color = residual)) +
    geom_point(size = 1.2, alpha = 0.7) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = color_limits) +
    coord_quickmap() +
    labs(title = "Residus spatialises (observe - predit), out-of-fold", subtitle = subtitle,
         x = "Longitude", y = "Latitude", color = "Residu\n(log10 NASC)") +
    theme_pipeline
}

plot_abs_residual_map <- function(obs_pred_all, subtitle = "", color_limits = NULL) {
  ggplot(obs_pred_all, aes(x = lon, y = lat, color = abs(residual))) +
    geom_point(size = 1.2, alpha = 0.7) +
    scale_color_viridis_c(option = "magma", direction = -1, limits = color_limits) +
    coord_quickmap() +
    labs(title = "Amplitude des erreurs de prediction, spatialisee", subtitle = subtitle,
         x = "Longitude", y = "Latitude", color = "|Residu|\n(log10 NASC)") +
    theme_pipeline
}

# ---------------------------------------------------------------------
# 6. Carte de prédiction sur grille (avec ou sans NA)
# ---------------------------------------------------------------------
plot_prediction_map <- function(grid_df, title, subtitle = "", limits = NULL) {
  ggplot(grid_df, aes(x = lon, y = lat, fill = NASC_pred)) +
    geom_raster() +
    scale_fill_viridis_c(limits = limits) +
    coord_quickmap() +
    theme_pipeline +
    labs(title = title, subtitle = subtitle, x = "Longitude", y = "Latitude", fill = "log10(NASC)")
}

# ---------------------------------------------------------------------
# 7. Métriques par fold : RMSE, R², variance intra-fold
# ---------------------------------------------------------------------
plot_metric_bar <- function(metrics_all, metric, fill = "steelblue", hline0 = FALSE, subtitle = "", ylim = NULL) {
  p <- ggplot(metrics_all, aes(x = fold_id, y = .data[[metric]])) +
    geom_col(fill = fill) +
    labs(title = paste0(metric, " par fold"), subtitle = subtitle, x = "Fold", y = metric) +
    theme_pipeline + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  if (hline0) p <- p + geom_hline(yintercept = 0, linetype = "dashed", color = "red")
  if (!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim)
  p
}

# ---------------------------------------------------------------------
# 8. Distance géographique / covariables, par fold
# ---------------------------------------------------------------------
plot_distances_per_fold <- function(metrics_all, subtitle = "") {
  dist_long <- metrics_all %>%
    select(fold_id, mean_geo_dist_km, sd_geo_dist_km,
           mean_covariate_dist, sd_covariate_dist) %>%
    pivot_longer(
      cols = c(mean_geo_dist_km, mean_covariate_dist),
      names_to = "metric", values_to = "mean_value"
    ) %>%
    mutate(sd_value = if_else(metric == "mean_geo_dist_km", sd_geo_dist_km, sd_covariate_dist)) %>%
    select(fold_id, metric, mean_value, sd_value)

  ggplot(dist_long, aes(x = fold_id, y = mean_value, fill = metric)) +
    geom_col(show.legend = FALSE) +
    geom_errorbar(aes(ymin = pmax(mean_value - sd_value, 0), ymax = mean_value + sd_value), width = 0.2) +
    facet_wrap(~metric, scales = "free_y", labeller = as_labeller(c(
      mean_geo_dist_km    = "Distance geographique (km), test -> train le plus proche",
      mean_covariate_dist = "Distance covariables (standardisee), test -> train le plus proche"
    ))) +
    labs(title = "Distance test -> train, par fold (moyenne +/- ecart-type intra-fold)",
         subtitle = subtitle, x = "Fold", y = NULL) +
    theme_pipeline + theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# ---------------------------------------------------------------------
# 9. Variance/moyenne des covariables numériques, train vs test, par fold
# ---------------------------------------------------------------------
plot_covariate_stats <- function(covariate_stats_all, subtitle = "") {
  ggplot(covariate_stats_all, aes(x = fold_id, y = mean, color = set)) +
    geom_point(position = position_dodge(width = 0.4)) +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2, position = position_dodge(width = 0.4)) +
    facet_wrap(~variable, scales = "free_y") +
    labs(title = "Moyenne +/- ecart-type de chaque covariable, par fold (train vs test)", subtitle = subtitle,
         x = "Fold", y = NULL, color = NULL) +
    theme_pipeline + theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

plot_covariate_variance <- function(covariate_stats_all, subtitle = "") {
  ggplot(covariate_stats_all, aes(x = fold_id, y = var, fill = set)) +
    geom_col(position = "dodge") +
    facet_wrap(~variable, scales = "free_y") +
    labs(title = "Variance de chaque covariable, par fold (train vs test)", subtitle = subtitle,
         x = "Fold", y = "Variance", fill = NULL) +
    theme_pipeline + theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# ---------------------------------------------------------------------
# 10. Distribution complète de chaque covariable, train vs test, par fold
# ---------------------------------------------------------------------
plot_distributions_per_fold <- function(numeric_dist_all, subtitle = "") {
  ggplot(numeric_dist_all, aes(x = value, fill = set, color = set)) +
    geom_density(alpha = 0.3) +
    facet_grid(fold_id ~ variable, scales = "free") +
    labs(title = "Distribution des covariables numeriques, train vs test, par fold", subtitle = subtitle,
         x = NULL, y = "Densite", fill = NULL, color = NULL) +
    theme_pipeline +
    theme(strip.text.y = element_text(angle = 0), axis.text = element_blank())
}

plot_fod_distribution <- function(fod_dist_all, subtitle = "") {
  ggplot(fod_dist_all, aes(x = fod, y = prop, fill = set)) +
    geom_col(position = "dodge") +
    facet_wrap(~fold_id) +
    labs(title = "Distribution des clusters FOD, train vs test, par fold", subtitle = subtitle,
         x = "Cluster FOD", y = "Proportion", fill = NULL) +
    theme_pipeline + theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# ---------------------------------------------------------------------
# 11. Localisation spatiale train/test (utile surtout pour le blocage)
# ---------------------------------------------------------------------
plot_spatial_train_test <- function(scheme, subtitle = "") {
  data <- scheme$data
  data$row_id <- seq_len(nrow(data))

  status_by_fold <- imap_dfr(scheme$folds, function(fold, fid) {
    df_status <- data %>% select(lon, lat)
    df_status$status  <- "Exclu par le buffer"
    df_status$status[fold$train] <- "Train conserve"
    df_status$status[fold$test]  <- "Test (bloc isole)"
    df_status$fold_id <- fid
    df_status
  })

  ggplot(status_by_fold, aes(x = lon, y = lat, color = status)) +
    geom_point(size = 0.5, alpha = 0.6) +
    scale_color_manual(values = c(
      "Test (bloc isole)"   = "red",
      "Train conserve"      = "steelblue",
      "Exclu par le buffer" = "grey80"
    )) +
    coord_quickmap() +
    facet_wrap(~fold_id) +
    labs(title = "Localisation spatiale train / test / exclu-par-buffer, par fold", subtitle = subtitle,
         x = "Longitude", y = "Latitude", color = NULL) +
    theme_pipeline + theme(legend.position = "bottom")
}

plot_points_colored_by_fold <- function(scheme, subtitle = "") {
  data <- scheme$data
  data$row_id <- seq_len(nrow(data))
  fold_lookup <- imap_dfr(scheme$folds, ~ tibble(row_id = .x$test, fold_id = .y))
  data <- data %>% left_join(fold_lookup, by = "row_id") %>%
    mutate(fold_id = ifelse(is.na(fold_id), "Non utilise", fold_id))

  ggplot(data, aes(x = lon, y = lat, color = fold_id)) +
    geom_point(size = 0.8, alpha = 0.7) +
    coord_quickmap() +
    labs(title = "Points colores par fold (bloc test)", subtitle = subtitle,
         x = "Longitude", y = "Latitude", color = "Fold") +
    theme_pipeline
}

# =====================================================================
# COMPARAISONS INTER-SCHEMAS / INTER-MODELES (axes de comparaison type
# "au-dela de l'accuracy globale", cf. discussion ViT/CNN)
# =====================================================================

# ---------------------------------------------------------------------
# 12. Importance des variables : naive vs blocage (le classement change-t-il ?)
# Detecte si une variable domine seulement grace a l'autocorrelation
# spatio-temporelle exploitee par la CV naive (raccourci, "Clever Hans"),
# et disparait une fois cette fuite supprimee par le blocage.
# ---------------------------------------------------------------------
plot_importance_comparison <- function(importance_all, subtitle = "") {
  summary_imp <- importance_all %>%
    group_by(model, scheme, variable) %>%
    summarise(mean_imp = mean(importance, na.rm = TRUE), .groups = "drop")

  ggplot(summary_imp, aes(x = reorder(variable, mean_imp), y = mean_imp, fill = scheme)) +
    geom_col(position = "dodge") +
    coord_flip() +
    facet_wrap(~model, scales = "free_x") +
    labs(title = "Importance des variables : naive vs schemas bloques", subtitle = subtitle,
         x = NULL, y = "Importance (moyenne inter-fold)", fill = "Schema") +
    theme_pipeline + theme(legend.position = "bottom")
}

# ---------------------------------------------------------------------
# 13. Courbe de calibration : le NASC predit suit-il la bonne moyenne
# par tranche de valeur, pas seulement correle (utile pour juger si le
# modele "regarde la bonne chose" de facon quantitative, sans biais
# systematique par gamme de valeurs).
# ---------------------------------------------------------------------
plot_calibration_curve <- function(obs_pred, n_bins = 10, subtitle = "") {
  binned <- obs_pred %>%
    mutate(bin = dplyr::ntile(pred, n_bins)) %>%
    group_by(bin) %>%
    summarise(mean_obs = mean(obs), mean_pred = mean(pred), n = dplyr::n(), .groups = "drop")

  ggplot(binned, aes(x = mean_pred, y = mean_obs)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    geom_line(color = "steelblue", alpha = 0.6) +
    geom_point(aes(size = n), color = "steelblue") +
    labs(title = sprintf("Courbe de calibration (%d bins par prediction)", n_bins), subtitle = subtitle,
         x = "NASC predit moyen (log10), par bin", y = "NASC observe moyen (log10), par bin", size = "n obs") +
    theme_pipeline
}

# ---------------------------------------------------------------------
# 14. RMSE vs distance test->train, par fold (diagnostic de fuite) :
# si le RMSE augmente nettement avec la distance, une partie du score
# "naive" est un artefact de proximite spatio-temporelle plutot qu'un
# vrai pouvoir predictif -- pas juste une observation qualitative, un
# vrai test quantitatif de la fuite.
# ---------------------------------------------------------------------
plot_rmse_vs_distance <- function(metrics_all, distance_col = "mean_geo_dist_km", subtitle = "") {
  ggplot(metrics_all, aes(x = .data[[distance_col]], y = rmse_test, color = scheme, shape = model)) +
    geom_point(size = 2, alpha = 0.8) +
    geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "black", linetype = "dashed", linewidth = 0.5) +
    labs(title = paste("RMSE test vs", distance_col, ", par fold"), subtitle = subtitle,
         x = distance_col, y = "RMSE test", color = "Schema", shape = "Modele") +
    theme_pipeline
}

# ---------------------------------------------------------------------
# 15. Carte de difference entre deux sorties de prediction (modeles,
# schemas, ou strategies de gestion du NA differentes)
# ---------------------------------------------------------------------
plot_map_difference <- function(grid_df, diff_col = "diff", title = "", subtitle = "", limits = NULL) {
  ggplot(grid_df, aes(x = lon, y = lat, fill = .data[[diff_col]])) +
    geom_raster() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = limits) +
    coord_quickmap() +
    theme_pipeline +
    labs(title = title, subtitle = subtitle, x = "Longitude", y = "Latitude", fill = "Difference\n(log10 NASC)")
}

# ---------------------------------------------------------------------
# 16. Heatmap de correlation entre toutes les cartes de sortie produites
# (une carte = une combinaison freq/schema/modele/gestion-du-NA)
# ---------------------------------------------------------------------
plot_correlation_heatmap <- function(cor_matrix, subtitle = "") {
  cor_long <- as.data.frame(cor_matrix) %>%
    tibble::rownames_to_column("layer_1") %>%
    pivot_longer(-layer_1, names_to = "layer_2", values_to = "correlation")

  ggplot(cor_long, aes(x = layer_1, y = layer_2, fill = correlation)) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.2f", correlation)), size = 2.5) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-1, 1)) +
    labs(title = "Correlation entre les cartes de prediction", subtitle = subtitle, x = NULL, y = NULL) +
    theme_pipeline + theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# ---------------------------------------------------------------------
# 17. Robustesse au bruit gaussien : RMSE (et sa degradation relative)
# en fonction du niveau de bruit ajoute aux covariables du test.
# Equivalent tabulaire d'une courbe de robustesse a la perturbation de
# texture en vision.
# ---------------------------------------------------------------------
plot_noise_robustness <- function(noise_df, group_col = "scheme", subtitle = "", ylim = NULL) {
  summary_df <- noise_df %>%
    group_by(.data[[group_col]], noise_level) %>%
    summarise(mean_rmse = mean(rmse_test), sd_rmse = sd(rmse_test), .groups = "drop")

  p <- ggplot(summary_df, aes(x = noise_level, y = mean_rmse, color = .data[[group_col]])) +
    geom_ribbon(aes(ymin = mean_rmse - sd_rmse, ymax = mean_rmse + sd_rmse, fill = .data[[group_col]]),
                alpha = 0.1, color = NA) +
    geom_line() +
    geom_point() +
    labs(title = "Robustesse au bruit gaussien (covariables du test)", subtitle = subtitle,
         x = "Niveau de bruit (fraction de l'ecart-type de la covariable)",
         y = "RMSE test (moyenne inter-fold)", color = NULL, fill = NULL) +
    theme_pipeline
  if (!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim)
  p
}

plot_noise_robustness_relative <- function(noise_df, group_col = "scheme", subtitle = "", ylim = NULL) {
  summary_df <- noise_df %>%
    group_by(.data[[group_col]], noise_level) %>%
    summarise(mean_rmse = mean(rmse_test), .groups = "drop") %>%
    group_by(.data[[group_col]]) %>%
    mutate(rmse_baseline = mean_rmse[noise_level == min(noise_level)],
           degradation_pct = 100 * (mean_rmse - rmse_baseline) / rmse_baseline) %>%
    ungroup()

  p <- ggplot(summary_df, aes(x = noise_level, y = degradation_pct, color = .data[[group_col]])) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_line() +
    geom_point() +
    labs(title = "Degradation relative du RMSE sous bruit croissant", subtitle = subtitle,
         x = "Niveau de bruit (fraction de l'ecart-type de la covariable)",
         y = "Degradation du RMSE (%, vs niveau 0)", color = NULL) +
    theme_pipeline
  if (!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim)
  p
}

# ---------------------------------------------------------------------
# 18. Residu vs latitude, par schema/modele -- diagnostic direct d'un
# gradient latitudinal non capture par les covariables : si la courbe
# lissee (LOESS) reste plate autour de 0, le gradient est bien absorbe ;
# si elle a une pente ou une forme systematique, il reste un biais lie a
# la latitude non modelise -- potentiellement pire en blocage qu'en
# naive (la CV naive le masque partiellement via la fuite spatiale).
# ---------------------------------------------------------------------
plot_residual_vs_latitude <- function(obs_pred_all, subtitle = "") {
  ggplot(obs_pred_all, aes(x = lat, y = residual, color = scheme)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    geom_point(alpha = 0.12, size = 0.5) +
    geom_smooth(method = "loess", se = TRUE, linewidth = 0.8, formula = y ~ x) +
    facet_wrap(~model) +
    labs(title = "Residu (observe - predit) vs latitude, par schema", subtitle = subtitle,
         x = "Latitude", y = "Residu (log10 NASC)", color = "Schema") +
    theme_pipeline
}

# ---------------------------------------------------------------------
# 19. Variogramme empirique (semi-variance vs distance) -- pour juger si
# un buffer choisi "a dire d'expert" est trop court/trop long par
# rapport a la portee reelle d'autocorrelation spatiale/temporelle des
# residus DETENDANCES (cf. 09_variogram.R pour le detrending).
# ---------------------------------------------------------------------
plot_variogram <- function(variogram_df, current_buffer = NULL, subtitle = "", x_label = "Distance (km)") {
  p <- ggplot(variogram_df, aes(x = lag_mid, y = semivariance)) +
    geom_point(aes(size = n_pairs), color = "steelblue", alpha = 0.8) +
    geom_line(color = "steelblue", alpha = 0.5) +
    labs(title = "Variogramme empirique des residus detendances", subtitle = subtitle,
         x = x_label, y = "Semi-variance", size = "n paires") +
    theme_pipeline
  if (!is.null(current_buffer)) {
    p <- p + geom_vline(xintercept = current_buffer, linetype = "dashed", color = "firebrick") +
      annotate("text", x = current_buffer, y = max(variogram_df$semivariance, na.rm = TRUE),
               label = paste0("buffer actuel = ", current_buffer), color = "firebrick",
               hjust = -0.05, vjust = 1, size = 3.2)
  }
  p
}

