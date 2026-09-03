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
RUN_TUNING               <- TRUE   # 10_run_tuning.R
RUN_TRAINING             <- TRUE   # 11_run_training.R
RUN_PREDICTION_SINGLE    <- TRUE   # 12_run_grid_prediction.R      (1 date)
RUN_PREDICTION_MULTIDATE <- TRUE   # 13_run_grid_prediction_multidate.R (133 dates -- LONG)
RUN_RFSRC_RECONSTRUCTION <- TRUE   # 14_run_rfsrc_reconstruction.R

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

# ---------------------------------------------------------------------
# Bilan final
# ---------------------------------------------------------------------
dt_total <- difftime(Sys.time(), t_start_global, units = "hours")
log_step(sprintf("PIPELINE COMPLETE TERMINEE -- duree totale : %.2f heures", as.numeric(dt_total)))
cat("Resultats dans :", normalizePath("outputs_pipeline"), "\n")
