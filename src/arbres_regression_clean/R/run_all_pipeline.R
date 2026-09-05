# =====================================================================
# run_all_pipeline.R -- lance TOUTE la pipeline, dans l'ordre
# =====================================================================
# A placer à la racine de nasc_pipeline/ (à côté du dossier R/).
# A lancer via RStudio : Tools -> Background Jobs -> Start Background Job
# (recommandé pour une exécution longue/de nuit -- libère la console et
# garantit l'ordre, cf. discussion précédente) -- ou directement avec
# source("run_all_pipeline.R") dans la console si tu es disponible pour
# suivre l'exécution.
#
# IMPORTANT : le répertoire de travail (getwd()) doit être le dossier
# nasc_pipeline (celui qui contient R/) AVANT de lancer ce script.
#
# Chaque étape peut être désactivée individuellement ci-dessous (utile
# pour relancer seulement la suite après une interruption, sans refaire
# ce qui a déjà tourné).

# ---------------------------------------------------------------------
# ETAPES A EXECUTER (mettre à FALSE pour sauter une étape déjà faite)
# ---------------------------------------------------------------------
RUN_TUNING               <- FALSE  # 10_run_tuning.R -- deja fait, on saute
RUN_TRAINING             <- TRUE   # 11_run_training.R
RUN_PREDICTION_SINGLE    <- TRUE   # 12_run_grid_prediction.R      (1 date)
RUN_PREDICTION_MULTIDATE <- TRUE   # 13_run_grid_prediction_multidate.R (133 dates -- LONG)
RUN_RFSRC_RECONSTRUCTION <- TRUE   # 14_run_rfsrc_reconstruction.R
RUN_TRANSFER_TUNING_TEST <- TRUE   # 15_run_transfer_tuning_test.R (cout de ne pas re-tuner)
RUN_CROSS_SCHEME_ANALYSIS <- TRUE  # 16_run_cross_scheme_analysis.R (importance/calibration/fuite)
RUN_MAP_COMPARISON       <- TRUE   # 17_run_map_comparison.R (necessite 12_ ET 14_ deja lances)
RUN_NOISE_ROBUSTNESS_TEST <- TRUE  # 18_run_noise_robustness_test.R (robustesse au bruit gaussien)
RUN_VARIOGRAM_ANALYSIS   <- TRUE   # 19_run_variogram_analysis.R (variogrammes, aide au choix des buffers)

# ---------------------------------------------------------------------
# Utilitaires de log
# ---------------------------------------------------------------------
t_start_global <- Sys.time()

log_step <- function(msg) {
  cat("\n", strrep("=", 70), "\n", sep = "")
  cat(msg, "--", format(Sys.time()), "\n")
  cat(strrep("=", 70), "\n\n", sep = "")
}

run_step <- function(condition, script_path, step_name) {
  if (!isTRUE(condition)) {
    log_step(paste0("[SAUTE] ", step_name, " (desactive en haut du script)"))
    return(invisible(NULL))
  }
  log_step(paste0("DEBUT : ", step_name))
  t0 <- Sys.time()
  source(script_path)
  dt <- difftime(Sys.time(), t0, units = "mins")
  log_step(sprintf("FIN : %s (duree : %.1f min)", step_name, as.numeric(dt)))
}

# ---------------------------------------------------------------------
# Verification du repertoire de travail AVANT de commencer
# ---------------------------------------------------------------------
if (!dir.exists("R") || !file.exists("R/00_config.R")) {
  stop(
    "Le dossier 'R/' (ou 'R/00_config.R') est introuvable depuis le ",
    "repertoire de travail actuel (", getwd(), "). ",
    "Fais setwd('chemin/vers/nasc_pipeline') avant de lancer ce script."
  )
}

# ---------------------------------------------------------------------
# Execution, dans l'ordre
# ---------------------------------------------------------------------
run_step(RUN_TUNING,               "R/10_run_tuning.R",                  "SCRIPT 1 - Tuning CART/RF/XGB")

# Verification avant training : on a besoin des fichiers de tuning
n_tuning_expected <- 2 * 5 * 3   # 2 freqs x 5 schemas x 3 modeles (CART/RF/XGB)
n_tuning_found <- length(list.files("outputs_pipeline/tuning", pattern = "\\.rds$"))
if (RUN_TRAINING && n_tuning_found < n_tuning_expected) {
  stop(sprintf(
    paste0(
      "Seulement %d/%d fichiers de tuning trouves dans outputs_pipeline/tuning/. ",
      "Le training a besoin de TOUS les fichiers de tuning -- verifie que ",
      "10_run_tuning.R s'est bien termine sans erreur avant de continuer ",
      "(ou mets RUN_TUNING <- TRUE pour le relancer)."
    ),
    n_tuning_found, n_tuning_expected
  ))
}

run_step(RUN_TRAINING,             "R/11_run_training.R",                "SCRIPT 2 - Entrainement CART/RF/XGB")

# Verification avant prediction : on a besoin des modeles entraines
if ((RUN_PREDICTION_SINGLE || RUN_PREDICTION_MULTIDATE || RUN_RFSRC_RECONSTRUCTION) &&
    !dir.exists("outputs_pipeline/training")) {
  stop(
    "Aucun modele entraine trouve dans outputs_pipeline/training/. ",
    "Lance 11_run_training.R (RUN_TRAINING <- TRUE) avant les etapes de prediction."
  )
}

run_step(RUN_PREDICTION_SINGLE,    "R/12_run_grid_prediction.R",         "SCRIPT 3 - Prediction grille (1 date)")
run_step(RUN_PREDICTION_MULTIDATE, "R/13_run_grid_prediction_multidate.R", "SCRIPT 3bis - Prediction grille (133 dates)")
run_step(RUN_RFSRC_RECONSTRUCTION, "R/14_run_rfsrc_reconstruction.R",    "SCRIPT 4 - Reconstruction randomForestSRC")

# Verification avant analyse inter-schemas / comparaison de tuning : on
# a besoin des CSV de diagnostics produits par 11_run_training.R
if ((RUN_TRANSFER_TUNING_TEST || RUN_CROSS_SCHEME_ANALYSIS) && !dir.exists("outputs_pipeline/training")) {
  stop(
    "Aucune sortie de training trouvee dans outputs_pipeline/training/. ",
    "Lance 11_run_training.R avant les etapes d'analyse inter-schemas."
  )
}

run_step(RUN_TRANSFER_TUNING_TEST,  "R/15_run_transfer_tuning_test.R",   "SCRIPT 5 - Cout de ne pas re-tuner")
run_step(RUN_CROSS_SCHEME_ANALYSIS, "R/16_run_cross_scheme_analysis.R",  "SCRIPT 6 - Diagnostics inter-schemas")

# Verification avant comparaison de cartes : on a besoin des grilles
# numeriques produites par 12_ ET 14_ (pas juste les PNG)
if (RUN_MAP_COMPARISON &&
    (!file.exists("outputs_pipeline/predictions/predictions_all_combined.rds"))) {
  stop(
    "Fichier outputs_pipeline/predictions/predictions_all_combined.rds introuvable. ",
    "Lance 12_run_grid_prediction.R (RUN_PREDICTION_SINGLE <- TRUE) avant la comparaison de cartes."
  )
}

run_step(RUN_MAP_COMPARISON, "R/17_run_map_comparison.R", "SCRIPT 7 - Comparaison des cartes de sortie")

# RUN_NOISE_ROBUSTNESS_TEST a besoin des memes fichiers de tuning que
# RUN_TRANSFER_TUNING_TEST -- meme verification.
if (RUN_NOISE_ROBUSTNESS_TEST && n_tuning_found < n_tuning_expected) {
  stop(
    "Fichiers de tuning incomplets -- lance 10_run_tuning.R avant ",
    "18_run_noise_robustness_test.R."
  )
}
run_step(RUN_NOISE_ROBUSTNESS_TEST, "R/18_run_noise_robustness_test.R", "SCRIPT 8 - Robustesse au bruit gaussien")

# RUN_VARIOGRAM_ANALYSIS n'a besoin d'aucun modele ni tuning -- calcul
# purement descriptif sur les donnees d'entrainement, independant.
run_step(RUN_VARIOGRAM_ANALYSIS, "R/19_run_variogram_analysis.R", "SCRIPT 9 - Variogrammes empiriques")

# ---------------------------------------------------------------------
# Bilan final
# ---------------------------------------------------------------------
dt_total <- difftime(Sys.time(), t_start_global, units = "hours")
log_step(sprintf("PIPELINE COMPLETE TERMINEE -- duree totale : %.2f heures", as.numeric(dt_total)))
cat("Resultats dans :", normalizePath("outputs_pipeline"), "\n")
