# =====================================================================
# generate_covariate_distance_significance_table.R -- test de
# correlation RMSE test <-> decalage train/test, COVARIABLE PAR
# COVARIABLE (pas juste la distance euclidienne agregee), par fold.
# =====================================================================
# Reutilise covariate_stats.csv et metrics_par_fold.csv deja sur disque
# (aucun reentrainement, aucune reprediction).
#
# Pour chaque covariable, on calcule un decalage standardise entre
# train et test, par fold : (moyenne_test - moyenne_train) / ecart_type_train
# (mesure type "d de Cohen" -- decalage exprime en nombre d'ecarts-types
# du train). On teste ensuite la correlation de Pearson entre la valeur
# ABSOLUE de ce decalage et le RMSE de test, a travers les folds d'une
# meme combinaison freq x modele x schema.
#
# Objectif : identifier si UNE covariable en particulier porte le
# decalage train/test (et donc potentiellement la difficulte
# d'extrapolation), plutot que de se limiter a une distance euclidienne
# agregee sur toutes les covariables (qui peut masquer l'effet d'une
# seule variable derriere le bruit des autres).
#
# A lancer depuis la racine de nasc_pipeline/ :
#   setwd("chemin/vers/arbres_regression_clean")
#   source("generate_covariate_distance_significance_table.R")
#
# Sortie : covariate_distance_significance_table.tex (longtable, a
# inclure via \input{...} -- necessite \usepackage{longtable, booktabs}).

source("R/00_config.R")
source("R/07_cross_scheme_utils.R")

SCHEME_LABELS <- c(
  "naive_RS_80_20" = "Naive RS 80/20",
  "blocked_spatial_1500x1000km" = "Bloc 1500$\\times$1000\\,km",
  "blocked_spatial_200x200km"   = "Bloc 200$\\times$200\\,km",
  "blocked_spatial_20x20km"     = "Bloc 20$\\times$20\\,km",
  "blocked_temporal_1j"         = "Bloc temporel 1j"
)

MIN_FOLDS_FOR_TEST <- 4

# Restriction demandee : 38 kHz uniquement, XGB, un tableau SEPARE par schema
TARGET_FREQ     <- 38
MODELS_TO_TEST  <- c("rf")
SCHEMES_TO_TEST <- c("blocked_spatial_1500x1000km", "blocked_spatial_200x200km")

configs <- list_training_configs() %>%
  dplyr::filter(freq == TARGET_FREQ, model %in% MODELS_TO_TEST, scheme %in% SCHEMES_TO_TEST)

metrics_all  <- load_tagged_csv(configs, "metrics_par_fold.csv")
covstats_all <- load_tagged_csv(configs, "covariate_stats.csv")

cat(sprintf("Configurations chargees : %d combinaisons freq x modele x schema\n", nrow(configs)))

if (!all(c("n_train", "n_test") %in% names(metrics_all))) {
  stop(
    "Colonnes 'n_train'/'n_test' introuvables dans metrics_par_fold.csv (colonnes presentes : ",
    paste(names(metrics_all), collapse = ", "), "). ",
    "Le test de Welch en a besoin -- verifie le nom exact de la colonne d'effectif test ",
    "dans 04_diagnostics.R::compute_fold_diagnostics() et adapte n_by_fold ci-dessous."
  )
}

# ---------------------------------------------------------------------
# Decalage standardise train/test, par covariable et par fold
# ---------------------------------------------------------------------
cov_wide <- covstats_all %>%
  dplyr::select(freq, model, scheme, fold_id, variable, set, mean, sd) %>%
  tidyr::pivot_wider(names_from = set, values_from = c(mean, sd)) %>%
  dplyr::mutate(
    mean_shift     = (mean_test - mean_train) / sd_train,
    abs_mean_shift = abs(mean_shift)
  ) %>%
  dplyr::filter(is.finite(abs_mean_shift))   # ecarte les cas sd_train = 0 (rares)

joined <- cov_wide %>%
  dplyr::left_join(
    metrics_all %>% dplyr::select(freq, model, scheme, fold_id, rmse_test),
    by = c("freq", "model", "scheme", "fold_id")
  )

# ---------------------------------------------------------------------
# Test de correlation, pour chaque combinaison freq x modele x schema x
# covariable
# ---------------------------------------------------------------------
run_cor_test <- function(df) {
  n_folds <- nrow(df)
  if (n_folds < MIN_FOLDS_FOR_TEST || length(unique(df$abs_mean_shift)) < 2) {
    return(tibble(n_folds = n_folds, r = NA_real_, p_value = NA_real_))
  }
  test <- suppressWarnings(stats::cor.test(df$rmse_test, df$abs_mean_shift, method = "pearson"))
  tibble(n_folds = n_folds, r = unname(test$estimate), p_value = test$p.value)
}

results <- joined %>%
  dplyr::group_by(freq, model, scheme, variable) %>%
  dplyr::group_modify(~ run_cor_test(.x)) %>%
  dplyr::ungroup()

sig_stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else if (p < 0.10) "." else ""
}
results$sig <- sapply(results$p_value, sig_stars)

write.csv(results, "covariate_distance_significance_results.csv", row.names = FALSE)

# ---------------------------------------------------------------------
# Generation d'UN TABLEAU LATEX SEPARE PAR SCHEMA (fichiers .tex distincts)
# ---------------------------------------------------------------------
write_scheme_table <- function(scheme, out_file) {
  scheme_label <- SCHEME_LABELS[[scheme]]

  lines <- c(
    "\\begin{table}[htbp]",
    "\\centering",
    sprintf(
      "\\caption{Test de corr\\'elation de Pearson entre RMSE de test et d\\'ecalage standardis\\'e train/test (type $d$ de Cohen), covariable par covariable -- %s, %d kHz, %s}",
      scheme_label, TARGET_FREQ, toupper(paste(MODELS_TO_TEST, collapse = "/"))
    ),
    sprintf("\\label{tab:covariate_distance_significance_%s}",
            gsub("[^a-zA-Z0-9]", "_", scheme)),
    "\\small",
    "\\begin{tabular}{lrrl}",
    "\\toprule",
    "Covariable & $n$ folds & $r$ (Pearson) & $p$-valeur \\\\",
    "\\midrule"
  )

  sub <- results %>%
    dplyr::filter(freq == TARGET_FREQ, scheme == !!scheme, model %in% MODELS_TO_TEST) %>%
    dplyr::arrange(dplyr::desc(abs(r)))

  for (i in seq_len(nrow(sub))) {
    row <- sub[i, ]
    r_cell <- if (is.na(row$r)) "n/a" else sprintf("%.3f%s", row$r, row$sig)
    p_cell <- if (is.na(row$p_value)) "n/a" else sprintf("%.3f", row$p_value)
    lines <- c(lines, sprintf(
      "\\texttt{%s} & %d & %s & %s \\\\",
      row$variable, row$n_folds, r_cell, p_cell
    ))
  }

  lines <- c(lines,
    "\\bottomrule",
    "\\end{tabular}",
    "\\vspace{0.3em}",
    "{\\footnotesize Significativit\\'e : *** $p<0{,}001$ ; ** $p<0{,}01$ ; * $p<0{,}05$ ; . $p<0{,}10$. Tri\\'e par $|r|$ d\\'ecroissant.}",
    "\\end{table}"
  )

  writeLines(lines, out_file)
  cat("Ecrit :", out_file, "\n")
}

write_scheme_table("blocked_spatial_1500x1000km", "covariate_distance_significance_1500x1000km.tex")
write_scheme_table("blocked_spatial_200x200km",   "covariate_distance_significance_200x200km.tex")

cat("\nInclusion : \\input{covariate_distance_significance_1500x1000km.tex} (necessite \\usepackage{booktabs})\n")
cat("            \\input{covariate_distance_significance_200x200km.tex}\n")

# ---------------------------------------------------------------------
# Tableau DETAILLE PAR FOLD : une ligne par fold, une colonne par
# covariable -- chaque cellule combine (1) le decalage standardise
# train/test (type d de Cohen) ET (2) sa PROPRE significativite,
# obtenue par un test t de Welch comparant DIRECTEMENT la distribution
# de cette covariable entre train et test, DANS CE FOLD PRECIS (pas une
# correlation a travers les folds -- un test independant par cellule).
#
# Le test de Welch pour deux echantillons independants se calcule
# directement a partir des moyennes/ecarts-types/effectifs deja
# disponibles (covariate_stats.csv + n_train/n_test de
# metrics_par_fold.csv), sans avoir besoin des observations individuelles :
#   t  = (mean_test - mean_train) / sqrt(sd_test^2/n_test + sd_train^2/n_train)
#   df = Welch-Satterthwaite (approx.)
# ---------------------------------------------------------------------
welch_t_test_from_summary <- function(mean1, sd1, n1, mean2, sd2, n2) {
  se2 <- sd1^2 / n1 + sd2^2 / n2
  if (!is.finite(se2) || se2 <= 0) return(list(t = NA_real_, p = NA_real_))
  t_stat <- (mean1 - mean2) / sqrt(se2)
  df <- se2^2 / ( (sd1^2/n1)^2/(n1-1) + (sd2^2/n2)^2/(n2-1) )
  if (!is.finite(df) || df <= 0) return(list(t = t_stat, p = NA_real_))
  p_val <- 2 * stats::pt(-abs(t_stat), df)
  list(t = t_stat, p = p_val)
}

write_fold_detail_table <- function(scheme, out_file) {
  scheme_label <- SCHEME_LABELS[[scheme]]

  # n_train/n_test par fold (constants pour toutes les covariables d'un
  # meme fold) -- necessaires pour le test de Welch
  n_by_fold <- metrics_all %>%
    dplyr::filter(freq == TARGET_FREQ, scheme == !!scheme, model %in% MODELS_TO_TEST) %>%
    dplyr::select(fold_id, n_train, n_test)

  cov_detail <- cov_wide %>%
    dplyr::filter(freq == TARGET_FREQ, scheme == !!scheme, model %in% MODELS_TO_TEST) %>%
    dplyr::left_join(n_by_fold, by = "fold_id") %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      welch = list(welch_t_test_from_summary(mean_test, sd_test, n_test, mean_train, sd_train, n_train)),
      t_stat = welch$t, p_value = welch$p
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(sig = sapply(p_value, sig_stars))

  rmse_by_fold <- metrics_all %>%
    dplyr::filter(freq == TARGET_FREQ, scheme == !!scheme, model %in% MODELS_TO_TEST) %>%
    dplyr::select(fold_id, rmse_test)

  # Cellule = "decalage(etoiles)" -- ex. "1.23**"
  cell_wide <- cov_detail %>%
    dplyr::mutate(cell = ifelse(is.na(mean_shift), "n/a", sprintf("%.2f%s", mean_shift, sig))) %>%
    dplyr::select(fold_id, variable, cell) %>%
    tidyr::pivot_wider(names_from = variable, values_from = cell) %>%
    dplyr::left_join(rmse_by_fold, by = "fold_id") %>%
    dplyr::arrange(fold_id)

  cov_cols <- setdiff(names(cell_wide), c("fold_id", "rmse_test"))
  n_cov <- length(cov_cols)

  lines <- c(
    "\\begin{table}[htbp]",
    "\\centering",
    sprintf(
      "\\caption{D\\'ecalage standardis\\'e train/test (type $d$ de Cohen) par covariable, par fold -- test t de Welch (train vs test, dans chaque fold) -- %s, %d kHz, %s}",
      scheme_label, TARGET_FREQ, toupper(paste(MODELS_TO_TEST, collapse = "/"))
    ),
    sprintf("\\label{tab:covariate_fold_detail_%s}", gsub("[^a-zA-Z0-9]", "_", scheme)),
    "\\resizebox{\\textwidth}{!}{%",
    "\\tiny",
    paste0("\\begin{tabular}{l", paste(rep("r", n_cov), collapse = ""), "r}"),
    "\\toprule",
    paste0("Fold & ", paste0("\\texttt{", cov_cols, "}", collapse = " & "), " & RMSE \\\\"),
    "\\midrule"
  )

  for (i in seq_len(nrow(cell_wide))) {
    row <- cell_wide[i, ]
    lines <- c(lines, sprintf(
      "%s & %s & %.3f \\\\",
      row$fold_id, paste(unlist(row[cov_cols]), collapse = " & "), row$rmse_test
    ))
  }

  lines <- c(lines,
    "\\bottomrule",
    "\\end{tabular}",
    "}",
    "\\vspace{0.3em}",
    "{\\footnotesize D\\'ecalage = (moyenne\\_test - moyenne\\_train) / \\'ecart-type\\_train. Chaque valeur est test\\'ee INDIVIDUELLEMENT (test t de Welch, train vs test, DANS CE FOLD) : *** $p<0{,}001$ ; ** $p<0{,}01$ ; * $p<0{,}05$ ; . $p<0{,}10$ ; sans \\'etoile : $p \\geq 0{,}10$.}",
    "\\end{table}"
  )

  writeLines(lines, out_file)
  cat("Ecrit :", out_file, "\n")
}

write_fold_detail_table("blocked_spatial_1500x1000km", "covariate_fold_detail_1500x1000km.tex")
write_fold_detail_table("blocked_spatial_200x200km",   "covariate_fold_detail_200x200km.tex")

cat("\nInclusion : \\input{covariate_fold_detail_1500x1000km.tex} (necessite \\usepackage{graphicx} pour \\resizebox)\n")
cat("            \\input{covariate_fold_detail_200x200km.tex}\n")

cat("\n=== Resume (covariables significatives a p<0.05, toutes combinaisons confondues) ===\n")
print(results %>% dplyr::filter(!is.na(p_value), p_value < 0.05) %>%
        dplyr::arrange(dplyr::desc(abs(r))) %>%
        dplyr::select(freq, model, scheme, variable, n_folds, r, p_value))
