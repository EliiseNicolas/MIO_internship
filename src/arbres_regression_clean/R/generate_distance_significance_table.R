# =====================================================================
# generate_distance_significance_table.R -- test de correlation
# RMSE test <-> distance geographique moyenne au train, par fold,
# pour chaque combinaison freq x modele x schema deja entrainee.
# =====================================================================
# Reutilise uniquement les metrics_par_fold.csv deja sur disque (aucun
# reentrainement, aucune reprediction) -- test de Pearson entre
# rmse_test et mean_geo_dist_km, a travers les folds d'une meme
# combinaison. Un schema avec trop peu de folds (< 4) pour qu'un test
# de correlation ait un sens est marque "n/a" plutot que teste.
#
# A lancer depuis la racine de nasc_pipeline/ :
#   setwd("chemin/vers/arbres_regression_clean")
#   source("generate_distance_significance_table.R")
#
# Sortie : distance_significance_table.tex (a inclure via \input{...},
# necessite \usepackage{booktabs}).

source("R/00_config.R")
source("R/07_cross_scheme_utils.R")

SCHEME_LABELS <- c(
  "naive_RS_80_20" = "Naive RS 80/20",
  "blocked_spatial_1500x1000km" = "Bloc 1500$\\times$1000\\,km",
  "blocked_spatial_200x200km"   = "Bloc 200$\\times$200\\,km",
  "blocked_spatial_20x20km"     = "Bloc 20$\\times$20\\,km",
  "blocked_temporal_1j"         = "Bloc temporel 1j"
)

MIN_FOLDS_FOR_TEST <- 4   # en dessous, un test de correlation n'a pas grand sens

# Restriction demandee : XGB uniquement, schemas 1500x1000km et 200x200km
MODELS_TO_TEST  <- c("xgb")
SCHEMES_TO_TEST <- c("blocked_spatial_1500x1000km", "blocked_spatial_200x200km")

configs <- list_training_configs() %>%
  dplyr::filter(model %in% MODELS_TO_TEST, scheme %in% SCHEMES_TO_TEST)
metrics_all <- load_tagged_csv(configs, "metrics_par_fold.csv")

cat(sprintf("Configurations chargees : %d lignes (freq x modele x schema x fold)\n", nrow(metrics_all)))

# ---------------------------------------------------------------------
# Test de correlation, pour chaque combinaison freq x modele x schema
# ---------------------------------------------------------------------
results <- metrics_all %>%
  group_by(freq, model, scheme) %>%
  group_modify(~ {
    n_folds <- nrow(.x)
    if (n_folds < MIN_FOLDS_FOR_TEST || length(unique(.x$mean_geo_dist_km)) < 2) {
      return(tibble(n_folds = n_folds, r = NA_real_, p_value = NA_real_))
    }
    test <- suppressWarnings(stats::cor.test(.x$rmse_test, .x$mean_geo_dist_km, method = "pearson"))
    tibble(n_folds = n_folds, r = unname(test$estimate), p_value = test$p.value)
  }) %>%
  ungroup()

sig_stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else if (p < 0.10) "." else ""
}
results$sig <- sapply(results$p_value, sig_stars)

write.csv(results, "distance_significance_results.csv", row.names = FALSE)

# ---------------------------------------------------------------------
# Generation du tableau LaTeX
# ---------------------------------------------------------------------
lines <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Test de corr\\'elation de Pearson entre RMSE de test et distance g\\'eographique moyenne au train, par fold, pour chaque combinaison fr\\'equence $\\times$ mod\\`ele $\\times$ sch\\'ema}",
  "\\label{tab:distance_significance}",
  "\\small",
  "\\begin{tabular}{llrrrl}",
  "\\toprule",
  "Fr\\'equence & Sch\\'ema & Mod\\`ele & $n$ folds & $r$ (Pearson) & $p$-valeur \\\\",
  "\\midrule"
)

for (freq in FREQS) {
  sub_freq <- results %>% dplyr::filter(freq == !!freq) %>% arrange(scheme, model)
  for (scheme in names(SCHEME_LABELS)) {
    sub_scheme <- sub_freq %>% dplyr::filter(scheme == !!scheme)
    if (nrow(sub_scheme) == 0) next
    for (i in seq_len(nrow(sub_scheme))) {
      row <- sub_scheme[i, ]
      freq_cell   <- if (i == 1 && scheme == names(SCHEME_LABELS)[1]) paste0(freq, " kHz") else ""
      scheme_cell <- if (i == 1) SCHEME_LABELS[[scheme]] else ""
      r_cell   <- if (is.na(row$r)) "n/a" else sprintf("%.3f%s", row$r, row$sig)
      p_cell   <- if (is.na(row$p_value)) "n/a" else sprintf("%.3f", row$p_value)
      lines <- c(lines, sprintf(
        "%s & %s & %s & %d & %s & %s \\\\",
        freq_cell, scheme_cell, toupper(row$model), row$n_folds, r_cell, p_cell
      ))
    }
  }
  lines <- c(lines, "\\midrule")
}
lines[length(lines)] <- "\\bottomrule"

lines <- c(lines,
  "\\end{tabular}",
  "\\vspace{0.3em}",
  "{\\footnotesize Significativit\\'e : *** $p<0{,}001$ ; ** $p<0{,}01$ ; * $p<0{,}05$ ; . $p<0{,}10$. ``n/a'' : moins de 4 folds disponibles, test non r\\'ealis\\'e.}",
  "\\end{table}"
)

writeLines(lines, "distance_significance_table.tex")
cat("\nEcrit : distance_significance_table.tex\n")
cat("Inclusion : \\input{distance_significance_table.tex} (necessite \\usepackage{booktabs})\n")

cat("\n=== Resume (combinaisons significatives a p<0.05) ===\n")
print(results %>% dplyr::filter(!is.na(p_value), p_value < 0.05) %>%
        select(freq, model, scheme, n_folds, r, p_value))
