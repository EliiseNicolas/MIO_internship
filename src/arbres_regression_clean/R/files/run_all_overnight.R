# =====================================================================
# run_all_overnight.R -- enchaîne tuning PUIS training, dans l'ordre
# =====================================================================
# A lancer via RStudio : Tools -> Background Jobs -> Start Background Job
# (choisir ce fichier). Ca tourne dans un process R séparé, ta console
# reste libre, et 11_run_training.R ne démarre qu'une fois 10_run_tuning.R
# TERMINE (même process = pas de risque de lire des fichiers de tuning
# pas encore écrits).
#
# IMPORTANT : le "Répertoire de travail" du Background Job doit être le
# dossier nasc_pipeline (celui qui contient R/) -- RStudio propose un
# champ "Working Directory" dans la fenêtre de lancement du job, choisis
# la même valeur que ton setwd() habituel.

cat("=== DEBUT TUNING ===", format(Sys.time()), "\n")
source("R/10_run_tuning.R")
cat("=== FIN TUNING ===", format(Sys.time()), "\n")

cat("=== DEBUT TRAINING ===", format(Sys.time()), "\n")
source("R/11_run_training.R")
cat("=== FIN TRAINING ===", format(Sys.time()), "\n")

cat("=== TOUT EST TERMINE ===", format(Sys.time()), "\n")
