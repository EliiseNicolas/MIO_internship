# =====================================================================
# generate_tuning_summary_table.R -- tableau LaTeX des meilleurs
# hyperparametres, a partir des .rds deja produits par 10_run_tuning.R
# =====================================================================
# Ne relit AUCUN plot, AUCUN gros fichier -- juste les .rds de tuning
# (30 fichiers, quelques Ko chacun). Genere un fichier .tex par modele,
# avec un tableau (une ligne par freq x schema) : hyperparametres
# retenus + RMSE test minimal obtenu pendant le tuning.
#
# A lancer depuis la racine de nasc_pipeline/ :
#   setwd("chemin/vers/arbres_regression_clean")
#   source("generate_tuning_summary_table.R")
#
# Sorties : tuning_summary_cart.tex, tuning_summary_rf.tex,
# tuning_summary_xgb.tex -- a inclure directement dans le rapport via
# \input{tuning_summary_cart.tex} (necessite \usepackage{booktabs}).

source("R/00_config.R")

tuning_dir <- path_out("tuning")

SCHEMES <- c(
  "naive_RS_80_20",
  paste0("blocked_spatial_", map_chr(SPATIAL_RESOLUTIONS, "label")),
  paste0("blocked_temporal_", map_chr(TEMPORAL_RESOLUTIONS, "label"))
)
SCHEME_LABELS <- c(
  "naive_RS_80_20" = "Naive RS 80/20",
  "blocked_spatial_1500x1000km" = "Bloc 1500$\\times$1000\\,km",
  "blocked_spatial_200x200km"   = "Bloc 200$\\times$200\\,km",
  "blocked_spatial_20x20km"     = "Bloc 20$\\times$20\\,km",
  "blocked_temporal_1j"         = "Bloc temporel 1j"
)

# ---------------------------------------------------------------------
# Charge tous les .rds d'un modele donne, extrait best_params + RMSE min
# ---------------------------------------------------------------------
build_summary_table <- function(model) {
  rows <- list()
  for (freq in FREQS) {
    for (scheme in SCHEMES) {
      path <- file.path(tuning_dir, sprintf("%s_%dkHz_%s.rds", model, freq, scheme))
      if (!file.exists(path)) {
        cat("[!] manquant :", path, "\n")
        next
      }
      tuning <- readRDS(path)
      best_rmse <- min(tuning$tuning_results$mean_rmse_test, na.rm = TRUE)
      rows[[length(rows) + 1]] <- c(
        list(freq = freq, scheme = scheme, rmse = best_rmse),
        tuning$best_params
      )
    }
  }
  bind_rows(lapply(rows, as_tibble))
}

# ---------------------------------------------------------------------
# Genere le code LaTeX (tabularx) pour un modele donne
# ---------------------------------------------------------------------
write_latex_table <- function(df, model, param_cols, param_labels, out_file) {
  lines <- character(0)
  lines <- c(lines,
    "\\begin{table}[htbp]",
    "\\centering",
    sprintf("\\caption{Meilleurs hyperparam\\`etres retenus lors du tuning -- %s}", toupper(model)),
    sprintf("\\label{tab:tuning_summary_%s}", model),
    "\\small",
    paste0("\\begin{tabularx}{\\textwidth}{l l ", paste(rep("r", length(param_cols)), collapse = " "), " r}"),
    "\\toprule",
    paste0("Fr\\'equence & Sch\\'ema & ", paste(param_labels, collapse = " & "), " & RMSE test \\\\"),
    "\\midrule"
  )

  for (freq in FREQS) {
    sub <- df %>% filter(freq == !!freq)
    for (i in seq_len(nrow(sub))) {
      row <- sub[i, ]
      scheme_label <- SCHEME_LABELS[[row$scheme]]
      param_vals <- sapply(param_cols, function(p) {
        v <- row[[p]]
        if (is.numeric(v)) formatC(v, digits = 4, format = "g") else as.character(v)
      })
      freq_cell <- if (i == 1) paste0(freq, " kHz") else ""
      lines <- c(lines, paste0(
        freq_cell, " & ", scheme_label, " & ",
        paste(param_vals, collapse = " & "), " & ",
        formatC(row$rmse, digits = 4, format = "f"), " \\\\"
      ))
    }
    lines <- c(lines, "\\midrule")
  }
  lines[length(lines)] <- "\\bottomrule"  # remplace le dernier \midrule par \bottomrule

  lines <- c(lines, "\\end{tabularx}", "\\end{table}")

  writeLines(lines, out_file)
  cat("Ecrit :", out_file, "\n")
}

# ---------------------------------------------------------------------
# Un tableau par modele (colonnes d'hyperparametres differentes)
# ---------------------------------------------------------------------
summary_cart <- build_summary_table("cart")
write_latex_table(summary_cart, "cart",
                   param_cols = c("cp", "minsplit", "maxdepth"),
                   param_labels = c("\\texttt{cp}", "\\texttt{minsplit}", "\\texttt{maxdepth}"),
                   out_file = "tuning_summary_cart.tex")

summary_rf <- build_summary_table("rf")
write_latex_table(summary_rf, "rf",
                   param_cols = c("mtry", "min.node.size", "num.trees"),
                   param_labels = c("\\texttt{mtry}", "\\texttt{min.node.size}", "\\texttt{num.trees}"),
                   out_file = "tuning_summary_rf.tex")

summary_xgb <- build_summary_table("xgb")
write_latex_table(summary_xgb, "xgb",
                   param_cols = c("max_depth", "eta", "min_child_weight", "nrounds"),
                   param_labels = c("\\texttt{max\\_depth}", "\\texttt{eta}", "\\texttt{min\\_child\\_weight}", "\\texttt{nrounds}"),
                   out_file = "tuning_summary_xgb.tex")

cat("\nTermine. 3 fichiers .tex generes dans :", getwd(), "\n")
cat("Inclusion dans le rapport : \\input{tuning_summary_cart.tex} (necessite \\usepackage{booktabs, tabularx})\n")
