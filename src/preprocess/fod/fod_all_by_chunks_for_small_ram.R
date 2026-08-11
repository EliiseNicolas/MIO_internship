# ============================================================
# FOD TEMPERATURE + SALINITE
# Version adaptée à ~20 Go de RAM
#
# Données :
#   661 longitude
#   241 latitude
#   19 profondeurs
#   216 dates
#
# Nombre de profils :
#   661 * 241 * 216 = 34 354 776
#
# Principe :
#   - jamais thetao_flat complet en RAM
#   - jamais so_flat complet en RAM
#   - calcul par blocs temporels
#   - PCA exacte sur tous les profils valides
#   - Mclust ajusté sur un échantillon aléatoire
#   - prédiction Mclust sur tous les profils par blocs
#   - probabilités sauvegardées
#   - statistiques des profils sauvegardées
# ============================================================


# ============================================================
# 1. LIBRAIRIES
# ============================================================

library(splines)
library(mclust)
library(fields)

# Les librairies suivantes ne sont pas nécessaires ici :
# library(ncdf4)
# library(rpart)
# library(fda)
# library(zoo)
# library(pbapply)
# library(glue)


# ============================================================
# 2. PARAMETRES
# ============================================================

rm(list = ls())
gc()

# ------------------------------------------------------------
# Fichier RDS contenant thetao, so, time, lon, lat, depth
# ------------------------------------------------------------

path_data <- "/run/media/mmolinet/KER22/MIO_internship_III/data_preprocessed/concat_temp_sal/thetao_so_crop_2018_2023.rds"

# ------------------------------------------------------------
# Dossier de sortie
# ------------------------------------------------------------

output_dir <- "~/Elisou/MIO_internship_III/data_preprocessed/fod/FOD_2018_2021_2022_2023"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


# ------------------------------------------------------------
# Paramètres B-spline
# ------------------------------------------------------------

K <- 25
r <- 4
lambda_spline <- 0.25


# ------------------------------------------------------------
# Nombre d'harmoniques FPCA
# ------------------------------------------------------------

nharm <- 6


# ------------------------------------------------------------
# Paramètres Mclust
# ------------------------------------------------------------

G_range <- 4:6

# IMPORTANT :
# Mclust est ajusté sur un échantillon et non sur
# les ~34 millions de profils.
#
# 100 000 est raisonnable avec 20 Go de RAM.
#
# Tu peux essayer 200000 si le calcul passe bien.
gmm_sample_n <- 100000

# seuil pour créer la classe transition
seuil <- 0.75

set.seed(1234)


# ------------------------------------------------------------
# Taille des blocs temporels
# ------------------------------------------------------------

# 10 jours = environ 0.5 Go de données brutes
# pour thetao + so
#
# Si tu as un problème de RAM :
#   block_time <- 5
#
# Si tout fonctionne bien :
#   block_time <- 20
#
block_time <- 10


# ============================================================
# 3. CHARGEMENT DES DONNEES
# ============================================================

cat("\nChargement des données RDS...\n")

data_final <- readRDS(path_data)

cat("Données chargées.\n")

thetao <- data_final$thetao
so     <- data_final$so
time   <- data_final$time
lon    <- data_final$lon
lat    <- data_final$lat
depth  <- data_final$depth


# ------------------------------------------------------------
# Dimensions
# ------------------------------------------------------------

dim_ds <- dim(thetao)

nlon   <- dim_ds[1]
nlat   <- dim_ds[2]
ndepth <- dim_ds[3]
ntime  <- dim_ds[4]

n_spatial <- nlon * nlat
n_total   <- n_spatial * ntime


cat("\nDimensions :\n")
cat("nlon   =", nlon, "\n")
cat("nlat   =", nlat, "\n")
cat("ndepth =", ndepth, "\n")
cat("ntime  =", ntime, "\n")
cat("Profils totaux =", format(n_total, big.mark = " "), "\n")


# ============================================================
# 4. TEMPS
# ============================================================

time_bis <- as.POSIXct(
  time * 3600,
  origin = "1950-01-01",
  tz = "UTC"
)


# ============================================================
# 5. FONCTION POUR EXTRAIRE UN BLOC
# ============================================================

# Cette fonction transforme :
#
#   lon x lat x depth x time
#
# en :
#
#   (lon x lat x time) x depth
#
# mais uniquement pour le bloc demandé.

get_flat_block <- function(x, t_start, t_end) {
  
  block <- x[
    , ,
    ,
    t_start:t_end,
    drop = FALSE
  ]
  
  flat <- matrix(
    aperm(
      block,
      c(1, 2, 4, 3)
    ),
    ncol = ndepth
  )
  
  rm(block)
  
  return(flat)
}


# ============================================================
# 6. BASE B-SPLINE
# ============================================================

cat("\nConstruction de la base B-spline...\n")

phi <- bs(
  depth,
  df = K,
  degree = r - 1,
  intercept = TRUE
)

phi <- as.matrix(phi)

cat("Dimensions phi :", dim(phi), "\n")


# ------------------------------------------------------------
# Matrice de pénalisation
# ------------------------------------------------------------

D2 <- diff(
  diag(K),
  differences = 2
)

R <- t(D2) %*% D2


# ------------------------------------------------------------
# Résolution
# ------------------------------------------------------------

A <- crossprod(phi) +
  lambda_spline * R

cat("Rang de A :", qr(A)$rank, "\n")

B <- solve(
  A,
  t(phi)
)

# B = K x depth


# ============================================================
# 7. PASSAGE 1 :
#    CREATION DU MASQUE + MOYENNE DES COEFFICIENTS
# ============================================================

cat("\n============================================\n")
cat("PASSAGE 1 : masque + moyenne coefficients\n")
cat("============================================\n")


# ------------------------------------------------------------
# Masque global
#
# dimensions :
# nlon x nlat x ntime
# ------------------------------------------------------------

mask_common <- array(
  FALSE,
  dim = c(nlon, nlat, ntime)
)


# Somme des coefficients
# 2K = température + salinité

sum_coef <- numeric(2 * K)

n_valid <- 0


# ------------------------------------------------------------
# Boucle temporelle
# ------------------------------------------------------------

for (t_start in seq(1, ntime, by = block_time)) {
  
  t_end <- min(
    t_start + block_time - 1,
    ntime
  )
  
  cat(
    "Bloc",
    t_start,
    "-",
    t_end,
    "/",
    ntime,
    "\n"
  )
  
  
  # ----------------------------------------------------------
  # Lecture bloc température
  # ----------------------------------------------------------
  
  theta_flat <- get_flat_block(
    thetao,
    t_start,
    t_end
  )
  
  
  # ----------------------------------------------------------
  # Lecture bloc salinité
  # ----------------------------------------------------------
  
  so_flat <- get_flat_block(
    so,
    t_start,
    t_end
  )
  
  
  # ----------------------------------------------------------
  # Profils sans NA/NaN/Inf
  # ----------------------------------------------------------
  
  valid <- apply(
    is.finite(theta_flat),
    1,
    all
  ) &
    apply(
      is.finite(so_flat),
      1,
      all
    )
  
  
  # ----------------------------------------------------------
  # Sauvegarde du masque
  # ----------------------------------------------------------
  
  ntime_block <- t_end - t_start + 1
  
  mask_common[
    , ,
    t_start:t_end
  ] <- array(
    valid,
    dim = c(
      nlon,
      nlat,
      ntime_block
    )
  )
  
  
  n_valid_block <- sum(valid)
  
  n_valid <- n_valid +
    n_valid_block
  
  
  if (n_valid_block > 0) {
    
    # --------------------------------------------------------
    # Coefficients température
    # --------------------------------------------------------
    
    coef_theta_block <-
      theta_flat[valid, , drop = FALSE] %*%
      t(B)
    
    
    # --------------------------------------------------------
    # Coefficients salinité
    # --------------------------------------------------------
    
    coef_so_block <-
      so_flat[valid, , drop = FALSE] %*%
      t(B)
    
    
    # --------------------------------------------------------
    # Coefficients combinés
    # --------------------------------------------------------
    
    coef_biv_block <- cbind(
      coef_theta_block,
      coef_so_block
    )
    
    
    # --------------------------------------------------------
    # Somme
    # --------------------------------------------------------
    
    sum_coef <- sum_coef +
      colSums(coef_biv_block)
    
    
    # --------------------------------------------------------
    # Libération mémoire
    # --------------------------------------------------------
    
    rm(
      coef_theta_block,
      coef_so_block,
      coef_biv_block
    )
  }
  
  
  rm(
    theta_flat,
    so_flat,
    valid
  )
  
  gc()
}


# ------------------------------------------------------------
# Moyenne globale des coefficients
# ------------------------------------------------------------

alpha_mean <- sum_coef / n_valid

rm(sum_coef)

gc()


cat("\nNombre de profils valides :\n")
cat(
  format(
    n_valid,
    big.mark = " "
  ),
  "\n"
)

cat(
  "Proportion valide :",
  n_valid / n_total,
  "\n"
)


# ------------------------------------------------------------
# Sauvegarde masque
# ------------------------------------------------------------

valid_idx <- which(mask_common)

saveRDS(
  list(
    mask = mask_common,
    valid_idx = valid_idx,
    n_valid = n_valid
  ),
  file.path(
    output_dir,
    "FOD_mask_2021.rds"
  )
)


# ============================================================
# 8. PASSAGE 2 :
#    COVARIANCE EXACTE DES COEFFICIENTS
# ============================================================

cat("\n============================================\n")
cat("PASSAGE 2 : covariance FPCA\n")
cat("============================================\n")


cov_sum <- matrix(
  0,
  nrow = 2 * K,
  ncol = 2 * K
)


for (t_start in seq(1, ntime, by = block_time)) {
  
  t_end <- min(
    t_start + block_time - 1,
    ntime
  )
  
  cat(
    "Bloc",
    t_start,
    "-",
    t_end,
    "/",
    ntime,
    "\n"
  )
  
  
  # ----------------------------------------------------------
  # Masque du bloc
  # ----------------------------------------------------------
  
  mask_block <- as.vector(
    mask_common[
      , ,
      t_start:t_end
    ]
  )
  
  
  if (sum(mask_block) == 0) {
    next
  }
  
  
  # ----------------------------------------------------------
  # Données
  # ----------------------------------------------------------
  
  theta_flat <- get_flat_block(
    thetao,
    t_start,
    t_end
  )
  
  so_flat <- get_flat_block(
    so,
    t_start,
    t_end
  )
  
  
  # ----------------------------------------------------------
  # Coefficients
  # ----------------------------------------------------------
  
  coef_theta_block <-
    theta_flat[mask_block, , drop = FALSE] %*%
    t(B)
  
  coef_so_block <-
    so_flat[mask_block, , drop = FALSE] %*%
    t(B)
  
  
  coef_biv_block <- cbind(
    coef_theta_block,
    coef_so_block
  )
  
  
  # ----------------------------------------------------------
  # Centrage
  # ----------------------------------------------------------
  
  coef_biv_block <- sweep(
    coef_biv_block,
    2,
    alpha_mean,
    "-"
  )
  
  
  # ----------------------------------------------------------
  # Contribution covariance
  # ----------------------------------------------------------
  
  cov_sum <- cov_sum +
    crossprod(coef_biv_block)
  
  
  rm(
    theta_flat,
    so_flat,
    coef_theta_block,
    coef_so_block,
    coef_biv_block,
    mask_block
  )
  
  gc()
}


# ------------------------------------------------------------
# Matrice covariance
# ------------------------------------------------------------

V <- cov_sum /
  (n_valid - 1)

rm(cov_sum)

gc()


# ============================================================
# 9. DECOMPOSITION SPECTRALE / FPCA
# ============================================================

cat("\n============================================\n")
cat("FPCA\n")
cat("============================================\n")


eig <- eigen(
  V,
  symmetric = TRUE
)


ord <- order(
  eig$values,
  decreasing = TRUE
)


lambda_all <- eig$values[ord]

U_all <- eig$vectors[
  ,
  ord
]


# ------------------------------------------------------------
# Variance expliquée
# ------------------------------------------------------------

prop_var <- lambda_all /
  sum(lambda_all)

cum_var <- cumsum(
  prop_var
)


cat("\nVariance expliquée :\n")

print(
  data.frame(
    PC = 1:length(prop_var),
    Eigenvalue = lambda_all,
    PropVar = prop_var,
    CumVar = cum_var
  )
)


# ------------------------------------------------------------
# Garder les 6 premières harmoniques
# ------------------------------------------------------------

U <- U_all[
  ,
  1:nharm,
  drop = FALSE
]

lambda <- lambda_all[
  1:nharm
]


cat(
  "\nNombre d'harmoniques :",
  nharm,
  "\n"
)


# ------------------------------------------------------------
# Sauvegarde modèle FPCA
# ------------------------------------------------------------

saveRDS(
  list(
    phi = phi,
    B = B,
    alpha_mean = alpha_mean,
    U = U,
    lambda = lambda,
    lambda_all = lambda_all,
    prop_var = prop_var,
    cum_var = cum_var,
    K = K,
    r = r,
    lambda_spline = lambda_spline,
    nharm = nharm,
    depth = depth
  ),
  file.path(
    output_dir,
    "FOD_FPCA_model_2021.rds"
  )
)


rm(
  V,
  eig,
  U_all
)

gc()


# ============================================================
# 10. ECHANTILLON POUR MCLUST
# ============================================================

cat("\n============================================\n")
cat("Préparation échantillon Mclust\n")
cat("============================================\n")


sample_n <- min(
  gmm_sample_n,
  n_valid
)


# ------------------------------------------------------------
# Tirage aléatoire parmi les profils valides
# ------------------------------------------------------------

sample_idx <- sort(
  sample(
    valid_idx,
    size = sample_n,
    replace = FALSE
  )
)


cat(
  "Nombre de profils utilisés pour Mclust :",
  format(
    sample_n,
    big.mark = " "
  ),
  "\n"
)


# ------------------------------------------------------------
# Matrice scores échantillon
# ------------------------------------------------------------

sample_scores <- matrix(
  NA_real_,
  nrow = sample_n,
  ncol = nharm
)


sample_counter <- 0


# ============================================================
# PASSAGE 3 :
# calcul des scores FPCA de l'échantillon
# ============================================================

for (t_start in seq(1, ntime, by = block_time)) {
  
  t_end <- min(
    t_start + block_time - 1,
    ntime
  )
  
  global_start <-
    (t_start - 1) * n_spatial + 1
  
  global_end <-
    t_end * n_spatial
  
  
  # indices échantillon dans ce bloc
  
  sel <- sample_idx[
    sample_idx >= global_start &
      sample_idx <= global_end
  ]
  
  
  if (length(sel) == 0) {
    next
  }
  
  
  # indices locaux
  
  local_idx <-
    sel - global_start + 1
  
  
  theta_flat <- get_flat_block(
    thetao,
    t_start,
    t_end
  )
  
  so_flat <- get_flat_block(
    so,
    t_start,
    t_end
  )
  
  
  # coefficients
  
  coef_theta_block <-
    theta_flat[
      local_idx,
      ,
      drop = FALSE
    ] %*%
    t(B)
  
  coef_so_block <-
    so_flat[
      local_idx,
      ,
      drop = FALSE
    ] %*%
    t(B)
  
  
  coef_biv_block <- cbind(
    coef_theta_block,
    coef_so_block
  )
  
  
  # centrage
  
  coef_biv_block <- sweep(
    coef_biv_block,
    2,
    alpha_mean,
    "-"
  )
  
  
  # scores
  
  scores_block <-
    coef_biv_block %*%
    U
  
  
  n_block <- nrow(
    scores_block
  )
  
  
  sample_scores[
    sample_counter + seq_len(n_block),
  ] <- scores_block
  
  
  sample_counter <-
    sample_counter + n_block
  
  
  rm(
    theta_flat,
    so_flat,
    coef_theta_block,
    coef_so_block,
    coef_biv_block,
    scores_block
  )
  
  gc()
}


# Vérification

stopifnot(
  sample_counter == sample_n
)


# ============================================================
# 11. MCLUST
# ============================================================

cat("\n============================================\n")
cat("MCLUST\n")
cat("============================================\n")


gmm <- Mclust(
  sample_scores,
  G = G_range,
  modelNames = "VVV"
)


# ------------------------------------------------------------
# Résumé
# ------------------------------------------------------------

print(
  summary(gmm)
)


nclust <- gmm$G

cat(
  "\nNombre de clusters sélectionné :",
  nclust,
  "\n"
)


# ------------------------------------------------------------
# Sauvegarde du modèle
# ------------------------------------------------------------

saveRDS(
  gmm,
  file.path(
    output_dir,
    "mclust_2021_results_nharm6.rds"
  )
)


rm(sample_scores)

gc()


# ============================================================
# 12. ALLOCATION DES RESULTATS
# ============================================================

# cluster pour tous les profils
cluster_flat <- rep(
  NA_integer_,
  n_total
)


# probabilités uniquement pour les profils valides
#
# n_valid x nclust
#
# Cela évite de créer une matrice de probabilités
# contenant des lignes inutiles pour les profils invalides.

probabilities <- matrix(
  NA_real_,
  nrow = n_valid,
  ncol = nclust
)


# cluster "dur"
cluster_valid <- integer(
  n_valid
)


# cluster avec transition
cluster_soft_valid <- integer(
  n_valid
)


# probabilité maximale
max_prob_valid <- numeric(
  n_valid
)


# ============================================================
# 13. PASSAGE 4 :
#     CLASSIFICATION DE TOUS LES PROFILS
# ============================================================

cat("\n============================================\n")
cat("CLASSIFICATION DE TOUS LES PROFILS\n")
cat("============================================\n")


valid_counter <- 0


# ------------------------------------------------------------
# Pour statistiques de profils
# ------------------------------------------------------------

sum_temp <- matrix(
  0,
  nrow = nclust + 1,
  ncol = ndepth
)

sum_sal <- matrix(
  0,
  nrow = nclust + 1,
  ncol = ndepth
)

count_cluster <- integer(
  nclust + 1
)


for (t_start in seq(1, ntime, by = block_time)) {
  
  t_end <- min(
    t_start + block_time - 1,
    ntime
  )
  
  cat(
    "Bloc",
    t_start,
    "-",
    t_end,
    "/",
    ntime,
    "\n"
  )
  
  
  # ----------------------------------------------------------
  # Masque
  # ----------------------------------------------------------
  
  mask_block <- as.vector(
    mask_common[
      , ,
      t_start:t_end
    ]
  )
  
  
  n_valid_block <- sum(
    mask_block
  )
  
  
  if (n_valid_block == 0) {
    next
  }
  
  
  # ----------------------------------------------------------
  # Données
  # ----------------------------------------------------------
  
  theta_flat <- get_flat_block(
    thetao,
    t_start,
    t_end
  )
  
  so_flat <- get_flat_block(
    so,
    t_start,
    t_end
  )
  
  
  # ----------------------------------------------------------
  # Coefficients
  # ----------------------------------------------------------
  
  coef_theta_block <-
    theta_flat[
      mask_block,
      ,
      drop = FALSE
    ] %*%
    t(B)
  
  coef_so_block <-
    so_flat[
      mask_block,
      ,
      drop = FALSE
    ] %*%
    t(B)
  
  
  # ----------------------------------------------------------
  # Coefficients centrés
  # ----------------------------------------------------------
  
  coef_biv_block <- cbind(
    coef_theta_block,
    coef_so_block
  )
  
  coef_biv_block <- sweep(
    coef_biv_block,
    2,
    alpha_mean,
    "-"
  )
  
  
  # ----------------------------------------------------------
  # Scores FPCA
  # ----------------------------------------------------------
  
  scores_block <-
    coef_biv_block %*%
    U
  
  
  # ----------------------------------------------------------
  # Prediction Mclust
  # ----------------------------------------------------------
  
  prediction <- predict(
    gmm,
    newdata = scores_block
  )
  
  
  z_block <- prediction$z
  
  cluster_block <-
    prediction$classification
  
  
  max_prob_block <-
    apply(
      z_block,
      1,
      max
    )
  
  
  # ----------------------------------------------------------
  # Soft classes
  # ----------------------------------------------------------
  
  cluster_soft_block <-
    cluster_block
  
  cluster_soft_block[
    max_prob_block < seuil
  ] <- 0
  
  
  # ----------------------------------------------------------
  # Indices globaux des profils valides
  # ----------------------------------------------------------
  
  global_idx_block <- which(
    mask_block
  ) +
    (t_start - 1) * n_spatial
  
  
  # ----------------------------------------------------------
  # Remplissage probabilités
  # ----------------------------------------------------------
  
  pos <- valid_counter +
    seq_len(n_valid_block)
  
  
  probabilities[
    pos,
  ] <- z_block
  
  
  cluster_valid[
    pos
  ] <- cluster_block
  
  
  cluster_soft_valid[
    pos
  ] <- cluster_soft_block
  
  
  max_prob_valid[
    pos
  ] <- max_prob_block
  
  
  # ----------------------------------------------------------
  # Carte spatiale/temporelle
  # ----------------------------------------------------------
  
  cluster_flat[
    global_idx_block
  ] <- cluster_soft_block
  
  
  valid_counter <-
    valid_counter + n_valid_block
  
  
  # ==========================================================
  # RECONSTRUCTION DES PROFILS
  #
  # uniquement pour le bloc courant
  # ==========================================================
  
  temp_rec_block <-
    coef_theta_block %*%
    t(phi)
  
  
  sal_rec_block <-
    coef_so_block %*%
    t(phi)
  
  
  # ----------------------------------------------------------
  # Moyennes par cluster
  #
  # rowsum évite de créer un gros sous-tableau
  # pour chaque cluster
  # ----------------------------------------------------------
  
  temp_sum_block <- rowsum(
    temp_rec_block,
    group = cluster_soft_block,
    reorder = FALSE
  )
  
  sal_sum_block <- rowsum(
    sal_rec_block,
    group = cluster_soft_block,
    reorder = FALSE
  )
  
  
  temp_groups <- as.integer(
    rownames(temp_sum_block)
  )
  
  sal_groups <- as.integer(
    rownames(sal_sum_block)
  )
  
  
  sum_temp[
    temp_groups + 1,
  ] <-
    sum_temp[
      temp_groups + 1,
    ] +
    temp_sum_block
  
  
  sum_sal[
    sal_groups + 1,
  ] <-
    sum_sal[
      sal_groups + 1,
    ] +
    sal_sum_block
  
  
  count_block <- tabulate(
    cluster_soft_block + 1,
    nbins = nclust + 1
  )
  
  
  count_cluster <-
    count_cluster +
    count_block
  
  
  # ----------------------------------------------------------
  # Libération mémoire
  # ----------------------------------------------------------
  
  rm(
    theta_flat,
    so_flat,
    coef_theta_block,
    coef_so_block,
    coef_biv_block,
    scores_block,
    prediction,
    z_block,
    cluster_block,
    max_prob_block,
    cluster_soft_block,
    global_idx_block,
    temp_rec_block,
    sal_rec_block,
    temp_sum_block,
    sal_sum_block,
    mask_block
  )
  
  gc()
}


# ============================================================
# 14. MOYENNES EXACTES DES PROFILS
# ============================================================

mean_temp_cl <- vector(
  "list",
  nclust + 1
)

mean_sal_cl <- vector(
  "list",
  nclust + 1
)


for (cl in 0:nclust) {
  
  if (count_cluster[cl + 1] > 0) {
    
    mean_temp_cl[[cl + 1]] <-
      sum_temp[cl + 1, ] /
      count_cluster[cl + 1]
    
    mean_sal_cl[[cl + 1]] <-
      sum_sal[cl + 1, ] /
      count_cluster[cl + 1]
    
  } else {
    
    mean_temp_cl[[cl + 1]] <-
      rep(NA_real_, ndepth)
    
    mean_sal_cl[[cl + 1]] <-
      rep(NA_real_, ndepth)
  }
}


rm(
  sum_temp,
  sum_sal
)

gc()


# ============================================================
# 15. QUARTILES EXACTS
#
# IMPORTANT :
#
# On ne peut pas conserver temp_rec et sal_rec pour
# les 34 millions de profils.
#
# On calcule donc les quartiles profondeur par profondeur.
#
# Pour chaque profondeur :
#   - allocation des valeurs par cluster
#   - passage sur tous les blocs
#   - quantile exact
#
# Cela nécessite 19 passages supplémentaires.
#
# C'est la partie la plus lente du script mais elle reste
# compatible avec 20 Go de RAM.
# ============================================================

cat("\n============================================\n")
cat("CALCUL DES QUARTILES EXACTS\n")
cat("============================================\n")


q1_temp_cl <- vector(
  "list",
  nclust + 1
)

q3_temp_cl <- vector(
  "list",
  nclust + 1
)

q1_sal_cl <- vector(
  "list",
  nclust + 1
)

q3_sal_cl <- vector(
  "list",
  nclust + 1
)


for (d in seq_len(ndepth)) {
  
  cat(
    "\nProfondeur",
    d,
    "/",
    ndepth,
    ":",
    depth[d],
    "m\n"
  )
  
  
  # ----------------------------------------------------------
  # Allocation des valeurs
  #
  # La somme des tailles = n_valid
  # ----------------------------------------------------------
  
  temp_values <- lapply(
    0:nclust,
    function(cl) {
      numeric(
        count_cluster[cl + 1]
      )
    }
  )
  
  sal_values <- lapply(
    0:nclust,
    function(cl) {
      numeric(
        count_cluster[cl + 1]
      )
    }
  )
  
  
  filled <- integer(
    nclust + 1
  )
  
  
  # ----------------------------------------------------------
  # Parcours de tous les blocs
  # ----------------------------------------------------------
  
  for (t_start in seq(1, ntime, by = block_time)) {
    
    t_end <- min(
      t_start + block_time - 1,
      ntime
    )
    
    
    mask_block <- as.vector(
      mask_common[
        , ,
        t_start:t_end
      ]
    )
    
    
    n_valid_block <- sum(
      mask_block
    )
    
    
    if (n_valid_block == 0) {
      next
    }
    
    
    # --------------------------------------------------------
    # Données
    # --------------------------------------------------------
    
    theta_flat <- get_flat_block(
      thetao,
      t_start,
      t_end
    )
    
    so_flat <- get_flat_block(
      so,
      t_start,
      t_end
    )
    
    
    # --------------------------------------------------------
    # Coefficients
    # --------------------------------------------------------
    
    coef_theta_block <-
      theta_flat[
        mask_block,
        ,
        drop = FALSE
      ] %*%
      t(B)
    
    coef_so_block <-
      so_flat[
        mask_block,
        ,
        drop = FALSE
      ] %*%
      t(B)
    
    
    # --------------------------------------------------------
    # Reconstruction UNIQUEMENT à la profondeur d
    #
    # Beaucoup moins de mémoire que temp_rec complet
    # --------------------------------------------------------
    
    temp_profile <- as.vector(
      coef_theta_block %*%
        phi[d, ]
    )
    
    sal_profile <- as.vector(
      coef_so_block %*%
        phi[d, ]
    )
    
    
    # --------------------------------------------------------
    # Classes du bloc
    # --------------------------------------------------------
    
    global_idx_block <- which(
      mask_block
    ) +
      (t_start - 1) * n_spatial
    
    
    classes_block <- cluster_flat[
      global_idx_block
    ]
    
    
    # --------------------------------------------------------
    # Stockage par cluster
    # --------------------------------------------------------
    
    for (cl in 0:nclust) {
      
      ind <- which(
        classes_block == cl
      )
      
      
      if (length(ind) > 0) {
        
        start_pos <-
          filled[cl + 1] + 1
        
        end_pos <-
          filled[cl + 1] + length(ind)
        
        
        temp_values[[cl + 1]][
          start_pos:end_pos
        ] <-
          temp_profile[ind]
        
        
        sal_values[[cl + 1]][
          start_pos:end_pos
        ] <-
          sal_profile[ind]
        
        
        filled[cl + 1] <-
          end_pos
      }
    }
    
    
    rm(
      theta_flat,
      so_flat,
      coef_theta_block,
      coef_so_block,
      temp_profile,
      sal_profile,
      global_idx_block,
      classes_block,
      mask_block
    )
    
    gc()
  }
  
  
  # ----------------------------------------------------------
  # Quantiles exacts
  # ----------------------------------------------------------
  
  for (cl in 0:nclust) {
    
    if (
      count_cluster[cl + 1] > 0
    ) {
      
      q_temp <- quantile(
        temp_values[[cl + 1]],
        probs = c(
          0.25,
          0.75
        ),
        names = FALSE
      )
      
      
      q_sal <- quantile(
        sal_values[[cl + 1]],
        probs = c(
          0.25,
          0.75
        ),
        names = FALSE
      )
      
      
      q1_temp_cl[[cl + 1]][d] <-
        q_temp[1]
      
      q3_temp_cl[[cl + 1]][d] <-
        q_temp[2]
      
      
      q1_sal_cl[[cl + 1]][d] <-
        q_sal[1]
      
      q3_sal_cl[[cl + 1]][d] <-
        q_sal[2]
      
    } else {
      
      q1_temp_cl[[cl + 1]][d] <- NA
      q3_temp_cl[[cl + 1]][d] <- NA
      
      q1_sal_cl[[cl + 1]][d] <- NA
      q3_sal_cl[[cl + 1]][d] <- NA
    }
  }
  
  
  rm(
    temp_values,
    sal_values,
    filled
  )
  
  gc()
}


# ============================================================
# 16. SAUVEGARDE DES PROBABILITES
# ============================================================

cat("\nSauvegarde des probabilités...\n")


saveRDS(
  list(
    probabilities = probabilities,
    valid_idx = valid_idx,
    nclust = nclust,
    threshold = seuil,
    cluster = cluster_valid,
    cluster_soft = cluster_soft_valid,
    max_prob = max_prob_valid
  ),
  file.path(
    output_dir,
    "FOD_cluster_probabilities_2021.rds"
  ),
  compress = FALSE
)


# On peut maintenant supprimer les probabilités
# de la RAM si on n'en a plus besoin.

rm(
  probabilities,
  cluster_valid,
  cluster_soft_valid,
  max_prob_valid
)

gc()


# ============================================================
# 17. CLUSTER MAP
# ============================================================

cluster_map <- array(
  cluster_flat,
  dim = c(
    nlon,
    nlat,
    ntime
  )
)


saveRDS(
  cluster_map,
  file.path(
    output_dir,
    "FOD_cluster_map_2021.rds"
  )
)


# ============================================================
# 18. SAUVEGARDE COORDONNEES
# ============================================================

saveRDS(
  list(
    lon = lon,
    lat = lat,
    depth = depth,
    time = time_bis
  ),
  file.path(
    output_dir,
    "FOD_coordinates_2021.rds"
  )
)


# ============================================================
# 19. SAUVEGARDE PROFILS DE CLUSTERS
# ============================================================

saveRDS(
  list(
    nclust = nclust,
    count_cluster = count_cluster,
    
    mean_temp_cl = mean_temp_cl,
    q1_temp_cl = q1_temp_cl,
    q3_temp_cl = q3_temp_cl,
    
    mean_sal_cl = mean_sal_cl,
    q1_sal_cl = q1_sal_cl,
    q3_sal_cl = q3_sal_cl,
    
    depth = depth
  ),
  file.path(
    output_dir,
    "FOD_cluster_profiles_2021.rds"
  )
)


# ============================================================
# 20. PLOTS TEMPERATURE
# ============================================================

cols <- 1:(nclust + 1)

cols_fill <- adjustcolor(
  cols,
  alpha.f = 0.25
)


png(
  file.path(
    output_dir,
    "temperature_clusters.png"
  ),
  width = 1800,
  height = 2800,
  res = 300
)


plot(
  NULL,
  xlim = range(
    unlist(q1_temp_cl),
    unlist(q3_temp_cl),
    na.rm = TRUE
  ),
  ylim = rev(
    range(depth)
  ),
  xlab = "Temperature (°C)",
  ylab = "Depth (m)"
)


for (cl in 0:nclust) {
  
  if (
    all(
      is.na(
        q1_temp_cl[[cl + 1]]
      )
    )
  ) {
    next
  }
  
  
  polygon(
    x = c(
      q1_temp_cl[[cl + 1]],
      rev(q3_temp_cl[[cl + 1]])
    ),
    y = c(
      depth,
      rev(depth)
    ),
    col = cols_fill[cl + 1],
    border = NA
  )
  
  
  lines(
    mean_temp_cl[[cl + 1]],
    depth,
    col = cols[cl + 1],
    lwd = 3
  )
}


legend(
  "bottomright",
  legend = c(
    "Transition",
    paste(
      "Cluster",
      1:nclust
    )
  ),
  col = cols,
  lwd = 3,
  bty = "n"
)


dev.off()


# ============================================================
# 21. PLOTS SALINITE
# ============================================================

png(
  file.path(
    output_dir,
    "salinity_clusters.png"
  ),
  width = 1800,
  height = 2800,
  res = 300
)


plot(
  NULL,
  xlim = range(
    unlist(q1_sal_cl),
    unlist(q3_sal_cl),
    na.rm = TRUE
  ),
  ylim = rev(
    range(depth)
  ),
  xlab = "Salinity (PSU)",
  ylab = "Depth (m)"
)


for (cl in 0:nclust) {
  
  if (
    all(
      is.na(
        q1_sal_cl[[cl + 1]]
      )
    )
  ) {
    next
  }
  
  
  polygon(
    x = c(
      q1_sal_cl[[cl + 1]],
      rev(q3_sal_cl[[cl + 1]])
    ),
    y = c(
      depth,
      rev(depth)
    ),
    col = cols_fill[cl + 1],
    border = NA
  )
  
  
  lines(
    mean_sal_cl[[cl + 1]],
    depth,
    col = cols[cl + 1],
    lwd = 3
  )
}


legend(
  "bottomright",
  legend = c(
    "Transition",
    paste(
      "Cluster",
      1:nclust
    )
  ),
  col = cols,
  lwd = 3,
  bty = "n"
)


dev.off()


# ============================================================
# 22. CARTES JOUR PAR JOUR
# ============================================================

cat("\nCréation des cartes...\n")


for (i in seq_along(time_bis)) {
  
  cat(
    "Carte",
    i,
    "/",
    length(time_bis),
    "\n"
  )
  
  
  png(
    file.path(
      output_dir,
      sprintf(
        "cluster_map_%s.png",
        format(
          time_bis[i],
          "%Y%m%d"
        )
      )
    ),
    width = 2000,
    height = 1500,
    res = 300
  )
  
  
  image.plot(
    x = lon,
    y = lat,
    z = cluster_map[, , i],
    
    col = cols,
    
    breaks = seq(
      -0.5,
      nclust + 0.5,
      1
    ),
    
    legend.only = FALSE,
    
    axis.args = list(
      at = 0:nclust,
      labels = c(
        "Transition",
        1:nclust
      )
    )
  )
  
  
  title(
    main = paste(
      "FOD map -",
      format(
        time_bis[i],
        "%Y-%m-%d"
      )
    )
  )
  
  
  dev.off()
}


# ============================================================
# 23. SAUVEGARDE PARAMETRES GENERAUX
# ============================================================

saveRDS(
  list(
    nlon = nlon,
    nlat = nlat,
    ndepth = ndepth,
    ntime = ntime,
    n_total = n_total,
    n_valid = n_valid,
    
    K = K,
    r = r,
    lambda_spline = lambda_spline,
    
    nharm = nharm,
    
    G_range = G_range,
    nclust = nclust,
    
    gmm_sample_n = gmm_sample_n,
    
    seuil = seuil,
    
    block_time = block_time
  ),
  file.path(
    output_dir,
    "FOD_parameters_2021.rds"
  )
)


# ============================================================
# 24. NETTOYAGE
# ============================================================

rm(
  thetao,
  so,
  data_final,
  mask_common,
  valid_idx,
  cluster_flat,
  cluster_map
)

gc()


cat("\n")
cat("============================================\n")
cat("      TRAITEMENT FOD TERMINE\n")
cat("============================================\n")
cat("\n")
cat(
  "Résultats sauvegardés dans :",
  output_dir,
  "\n"
)