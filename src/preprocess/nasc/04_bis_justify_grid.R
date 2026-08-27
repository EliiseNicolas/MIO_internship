# ============================================================
# AUTOCORRELATION SPATIALE DES PIGMENTS
#
# Objectif :
#   Tester l'effet de la résolution spatiale sur :
#     1. la variance locale moyenne
#     2. le Moran's I global
#
# Pour chaque facteur :
#
#   grille originale
#        ↓
#   agrégation facteur x facteur
#        ↓
#   nouvelle grille
#        ↓
#   variance locale 3x3
#        ↓
#   Moran global
#
# Les calculs sont faits date par date et pigment par pigment
# afin de limiter la mémoire utilisée.
# ============================================================


# ============================================================
# 0 - Paramètres
# ============================================================

freq <- 200

path_pig <- paste0(
  "F:/data_elise/pigmeann/",
  "pigments_2018_2021_2022_2023_crop.rds"
)

pigs <- readRDS(path_pig)

# Facteurs à tester
facteurs <- c(1, 2, 4, 8)

# Exemple :
# facteur = 1 -> grille originale
# facteur = 2 -> moyenne 2x2
# facteur = 4 -> moyenne 4x4
# facteur = 8 -> moyenne 8x8


# ============================================================
# 1 - Identifier les pigments
# ============================================================

pig_names <- grep("^c_cond_", names(pigs), value = TRUE)

cat("Pigments détectés :", paste(pig_names, collapse = ", "), "\n")

npig <- length(pig_names)


# ============================================================
# 2 - Supprimer les dates complètement NA
# ============================================================

keep <- sapply(seq_along(pigs$date), function(i) {
  any(sapply(pig_names, function(p) {
    any(!is.na(pigs[[p]][i, , ]))
  }))
})

cat("Dates initiales :", length(pigs$date), "\n")
cat("Dates supprimées :", sum(!keep), "\n")
cat("Dates conservées :", sum(keep), "\n")

pigs_filtered <- pigs
pigs_filtered$date <- pigs$date[keep]

for (p in pig_names) {
  pigs_filtered[[p]] <- pigs[[p]][keep, , , drop = FALSE]
}


# ============================================================
# 3 - Dimensions
# ============================================================

date <- pigs_filtered$date
lon <- pigs_filtered$lon
lat <- pigs_filtered$lat

ndate <- length(date)
nx <- length(lon)
ny <- length(lat)

lon_res <- mean(abs(diff(lon)))
lat_res <- mean(abs(diff(lat)))

cat("Grille originale :", nx, "x", ny, "\n")
cat("Résolution originale :", round(lon_res, 4), "° x", round(lat_res, 4), "°\n")


# ============================================================
# 4 - Fonction d'agrégation spatiale
# ============================================================
#
# Exemple facteur = 4 :
#
# 4 x 4 pixels originaux
#          ↓
#       1 pixel
#
# La moyenne est calculée en ignorant les NA.
#
# Les pixels qui ne rentrent pas dans un bloc complet sont
# retirés sur le bord.
# ============================================================

aggregate_grid <- function(x, facteur) {
  nx <- nrow(x)
  ny <- ncol(x)
  
  # dimensions utilisables
  nx_new <- floor(nx / facteur)
  ny_new <- floor(ny / facteur)
  
  nx_use <- nx_new * facteur
  ny_use <- ny_new * facteur
  
  # retirer les pixels de bord si nécessaire
  x <- x[1:nx_use, 1:ny_use, drop = FALSE]
  
  # matrice de sortie
  out <- matrix(NA_real_, nrow = nx_new, ncol = ny_new)
  
  # Agrégation
  #
  # On boucle sur la nouvelle grille et non sur la grille
  # originale.
  
  for (i in seq_len(nx_new)) {
    x1 <- (i - 1) * facteur + 1
    x2 <- i * facteur
    
    for (j in seq_len(ny_new)) {
      y1 <- (j - 1) * facteur + 1
      y2 <- j * facteur
      
      block <- x[x1:x2, y1:y2, drop = FALSE]
      
      out[i, j] <- if (all(is.na(block))) {
        NA_real_
      } else {
        mean(block, na.rm = TRUE)
      }
    }
  }
  
  out
}


# ============================================================
# 5 - Coordonnées de la nouvelle grille
# ============================================================

aggregate_coordinates <- function(coord, facteur) {
  n_new <- floor(length(coord) / facteur)
  
  coord_use <- coord[1:(n_new * facteur)]
  
  sapply(seq_len(n_new), function(i) {
    mean(coord_use[((i - 1) * facteur + 1):(i * facteur)])
  })
}


# ============================================================
# 6 - Fonction de décalage
# ============================================================

shift_matrix <- function(x, dx, dy) {
  nx <- nrow(x)
  ny <- ncol(x)
  
  out <- matrix(NA_real_, nrow = nx, ncol = ny)
  
  x_from <- max(1, 1 - dx):min(nx, nx - dx)
  x_to <- max(1, 1 + dx):min(nx, nx + dx)
  
  y_from <- max(1, 1 - dy):min(ny, ny - dy)
  y_to <- max(1, 1 + dy):min(ny, ny + dy)
  
  out[x_to, y_to] <- x[x_from, y_from]
  
  out
}


# ============================================================
# 7 - Variance locale
# ============================================================
#
# Pour chaque pixel :
#
#       voisin voisin voisin
#       voisin  PIXEL  voisin
#       voisin voisin voisin
#
# La variance est calculée sur les 9 valeurs disponibles.
# ============================================================

local_variance <- function(x) {
  neighbours <- list(
    shift_matrix(x, -1, -1),
    shift_matrix(x, -1,  0),
    shift_matrix(x, -1,  1),
    shift_matrix(x,  0, -1),
    shift_matrix(x,  0,  1),
    shift_matrix(x,  1, -1),
    shift_matrix(x,  1,  0),
    shift_matrix(x,  1,  1)
  )
  
  # Nombre de valeurs disponibles
  n_valid <- !is.na(x)
  
  # Somme
  sum_values <- ifelse(is.na(x), 0, x)
  
  # Somme des carrés
  sum_sq <- ifelse(is.na(x), 0, x^2)
  
  for (v in neighbours) {
    valid <- !is.na(v)
    
    n_valid <- n_valid + valid
    sum_values <- sum_values + ifelse(valid, v, 0)
    sum_sq <- sum_sq + ifelse(valid, v^2, 0)
  }
  
  # Moyenne locale
  mean_local <- sum_values / n_valid
  
  # Variance échantillonnale
  variance_local <- (
    sum_sq - n_valid * mean_local^2
  ) / (n_valid - 1)
  
  # Au moins 2 valeurs nécessaires
  variance_local[n_valid < 2] <- NA_real_
  
  variance_local
}


# ============================================================
# 8 - Moran global
# ============================================================
#
# Poids :
#
#     1 1 1
#     1 X 1
#     1 1 1
#
# Chaque voisin a le même poids.
#
# Les paires contenant un NA sont ignorées.
# ============================================================

global_moran <- function(x) {
  valid <- !is.na(x)
  
  n <- sum(valid)
  
  if (n < 2) {
    return(NA_real_)
  }
  
  # moyenne globale
  mean_x <- mean(x, na.rm = TRUE)
  
  # écarts à la moyenne
  z <- x - mean_x
  
  # dénominateur
  denominator <- sum(z[valid]^2)
  
  if (denominator == 0 || is.na(denominator)) {
    return(NA_real_)
  }
  
  # 8 voisins
  neighbours <- list(
    shift_matrix(z, -1, -1),
    shift_matrix(z, -1,  0),
    shift_matrix(z, -1,  1),
    shift_matrix(z,  0, -1),
    shift_matrix(z,  0,  1),
    shift_matrix(z,  1, -1),
    shift_matrix(z,  1,  0),
    shift_matrix(z,  1,  1)
  )
  
  numerator <- 0
  W <- 0
  
  for (v in neighbours) {
    ok <- !is.na(z) & !is.na(v)
    
    numerator <- numerator + sum(z[ok] * v[ok])
    W <- W + sum(ok)
  }
  
  if (W == 0) {
    return(NA_real_)
  }
  
  # Moran's I
  I <- (n / W) * numerator / denominator
  
  I
}


# ============================================================
# 9 - Fonction principale pour UN facteur
# ============================================================

calculate_factor <- function(facteur, pigs, pig_names, date, lon, lat) {
  cat("\n")
  cat("============================================\n")
  cat("FACTEUR :", facteur, "\n")
  cat("============================================\n")
  
  # ----------------------------------------------------------
  # Nouvelle grille
  # ----------------------------------------------------------
  
  nx_new <- floor(length(lon) / facteur)
  ny_new <- floor(length(lat) / facteur)
  
  lon_new <- aggregate_coordinates(lon, facteur)
  lat_new <- aggregate_coordinates(lat, facteur)
  
  cat("Nouvelle grille :", nx_new, "x", ny_new, "\n")
  
  cat(
    "Nouvelle résolution :",
    round(mean(abs(diff(lon_new))), 4),
    "° x",
    round(mean(abs(diff(lat_new))), 4),
    "°\n"
  )
  
  # ----------------------------------------------------------
  # Résultat variance locale
  # ----------------------------------------------------------
  
  mean_local_variance <- array(
    NA_real_,
    dim = c(length(date), nx_new, ny_new)
  )
  
  # ----------------------------------------------------------
  # Moran par pigment
  # ----------------------------------------------------------
  
  moran_global_pigments <- matrix(
    NA_real_,
    nrow = length(date),
    ncol = length(pig_names),
    dimnames = list(as.character(date), pig_names)
  )
  
  # ----------------------------------------------------------
  # Moran moyen
  # ----------------------------------------------------------
  
  mean_moran_global <- rep(NA_real_, length(date))
  names(mean_moran_global) <- as.character(date)
  
  # ----------------------------------------------------------
  # Boucle temporelle
  # ----------------------------------------------------------
  
  for (d in seq_along(date)) {
    cat("\rDate ", d, "/", length(date), " : ", as.character(date[d]), "          ")
    
    # accumulateur variance
    variance_sum <- matrix(0, nx_new, ny_new)
    variance_n <- matrix(0, nx_new, ny_new)
    
    # --------------------------------------------------------
    # Boucle pigments
    # --------------------------------------------------------
    
    for (p in seq_along(pig_names)) {
      # ------------------------------------------------------
      # Extraction grille originale
      # ------------------------------------------------------
      
      x <- pigs[[pig_names[p]]][d, , , drop = TRUE]
      
      # ------------------------------------------------------
      # REGRIDDING
      # ------------------------------------------------------
      
      x_agg <- aggregate_grid(x, facteur)
      
      # ------------------------------------------------------
      # VARIANCE LOCALE
      # ------------------------------------------------------
      
      v <- local_variance(x_agg)
      
      ok <- !is.na(v)
      
      variance_sum[ok] <- variance_sum[ok] + v[ok]
      variance_n[ok] <- variance_n[ok] + 1
      
      # ------------------------------------------------------
      # MORAN GLOBAL
      # ------------------------------------------------------
      
      moran_global_pigments[d, p] <- global_moran(x_agg)
    }
    
    # --------------------------------------------------------
    # Moyenne variance sur les pigments
    # --------------------------------------------------------
    
    variance_mean <- variance_sum / variance_n
    variance_mean[variance_n == 0] <- NA_real_
    
    mean_local_variance[d, , , drop = TRUE] <- variance_mean
    
    # --------------------------------------------------------
    # Moyenne Moran sur les pigments
    # --------------------------------------------------------
    
    mean_moran_global[d] <- mean(
      moran_global_pigments[d, ],
      na.rm = TRUE
    )
  }
  
  cat("\n")
  
  # ----------------------------------------------------------
  # Retour
  # ----------------------------------------------------------
  
  list(
    facteur = facteur,
    lon = lon_new,
    lat = lat_new,
    date = date,
    resolution_lon = mean(abs(diff(lon_new))),
    resolution_lat = mean(abs(diff(lat_new))),
    mean_local_variance = mean_local_variance,
    moran_global_pigments = moran_global_pigments,
    mean_moran_global = mean_moran_global
  )
}


# ============================================================
# 10 - Lancer tous les facteurs
# ============================================================

results <- vector("list", length(facteurs))

names(results) <- paste0("facteur_", facteurs)

for (i in seq_along(facteurs)) {
  results[[i]] <- calculate_factor(
    facteur = facteurs[i],
    pigs = pigs_filtered,
    pig_names = pig_names,
    date = date,
    lon = lon,
    lat = lat
  )
}

cat("\n")
cat("============================================\n")
cat("TOUS LES CALCULS SONT TERMINÉS\n")
cat("============================================\n")

#### Extraire les résultats
# Facteur 1
results$facteur_1

# Facteur 4
results$facteur_4

# Moran moyen pour facteur 4
results$facteur_4$mean_moran_global

# Moran de Chla pour facteur 4
results$facteur_4$moran_global_pigments[, "c_cond_Chla"]

# Variance locale moyenne pour facteur 4
results$facteur_4$mean_local_variance

#### Comparer
comparison <- data.frame(
  
  facteur = facteurs,
  
  resolution_lon = sapply(
    results,
    function(x) x$resolution_lon
  ),
  
  resolution_lat = sapply(
    results,
    function(x) x$resolution_lat
  ),
  
  mean_moran = sapply(
    results,
    function(x) mean(
      x$mean_moran_global,
      na.rm = TRUE
    )
  ),
  
  sd_moran = sapply(
    results,
    function(x) sd(
      x$mean_moran_global,
      na.rm = TRUE
    )
  )
  
)

comparison
