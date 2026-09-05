# =====================================================================
# 07_cross_scheme_utils.R -- comparer schémas/modèles entre eux
# =====================================================================
# Utilisé par 15_run_transfer_tuning_test.R, 16_run_cross_scheme_analysis.R
# et 17_run_map_comparison.R : lit les CSV déjà enregistrés par
# 11_run_training.R (metrics_par_fold.csv, obs_pred_all.csv,
# importance_all.csv) sous outputs_pipeline/training/<freq>kHz/<model>/<schema>/,
# et les combine en ajoutant des colonnes freq/model/scheme (absentes des
# CSV individuels, qui ne connaissent que leur propre configuration).

# ---------------------------------------------------------------------
# Enumere toutes les configurations (freq, modele, schema) pour
# lesquelles 11_run_training.R a produit une sortie.
# ---------------------------------------------------------------------
list_training_configs <- function(training_dir = path_out("training"),
                                   models = c("cart", "rf", "xgb")) {
  if (!dir.exists(training_dir)) {
    stop("Dossier introuvable : ", training_dir, " -- as-tu bien lance 11_run_training.R ?")
  }
  freq_dirs <- list.dirs(training_dir, recursive = FALSE)

  configs <- purrr::map_dfr(freq_dirs, function(fd) {
    freq <- as.integer(gsub("kHz$", "", basename(fd)))
    purrr::map_dfr(models, function(m) {
      model_dir <- file.path(fd, m)
      if (!dir.exists(model_dir)) return(NULL)
      scheme_dirs <- list.dirs(model_dir, recursive = FALSE)
      tibble(freq = freq, model = m, scheme = basename(scheme_dirs), dir = scheme_dirs)
    })
  })

  if (nrow(configs) == 0) {
    stop("Aucune configuration trouvee sous ", training_dir, " -- verifie que 11_run_training.R s'est bien execute.")
  }
  configs
}

# ---------------------------------------------------------------------
# Charge un fichier CSV donne (ex. "metrics_par_fold.csv") pour TOUTES
# les configurations enumerees par list_training_configs(), en ajoutant
# les colonnes freq/model/scheme. Ignore silencieusement les
# configurations ou le fichier n'existe pas (avec un avis).
# ---------------------------------------------------------------------
load_tagged_csv <- function(configs, filename) {
  purrr::pmap_dfr(configs, function(freq, model, scheme, dir) {
    path <- file.path(dir, filename)
    if (!file.exists(path)) {
      warning("Fichier manquant, ignore : ", path)
      return(NULL)
    }
    read.csv(path) %>% mutate(freq = freq, model = model, scheme = scheme)
  })
}
