library(ncdf4)

path     <- "C:/Users/mmolinet/elisou_ta_stagiaire_pref/data_preprocessed/concat_temp_sal/thetao_so_crop_2018_2021_2022_2023.nc"
out_path <- "C:/Users/mmolinet/elisou_ta_stagiaire_pref/data_preprocessed/fod/FOD_2018_2021_2022_2023"

target_date <- as.Date("2023-01-26")

ds <- nc_open(path)

# ------------------------------------------------------------
# 1) Lire uniquement le vecteur time (leger, 1D) pour trouver l'indice
# ------------------------------------------------------------
time     <- ncvar_get(ds, "time")
time_bis <- as.POSIXct(time * 3600, origin = "1950-01-01", tz = "UTC")

time_idx <- which(as.Date(time_bis) == target_date)

if (length(time_idx) == 0) {
  nc_close(ds)
  stop("Date non trouvee dans le fichier: ", target_date)
}
if (length(time_idx) > 1) {
  warning("Plusieurs pas de temps trouves pour cette date, on prend le premier.")
  time_idx <- time_idx[1]
}

# ------------------------------------------------------------
# 2) Dimensions spatiales (lon, lat, depth) - lecture complete car 1D et legere
# ------------------------------------------------------------
lon   <- ncvar_get(ds, "longitude")
lat   <- ncvar_get(ds, "latitude")
depth <- ncvar_get(ds, "depth")

nlon   <- length(lon)
nlat   <- length(lat)
ndepth <- length(depth)

# ------------------------------------------------------------
# 3) Lire seulement la tranche temporelle voulue pour thetao et so
#    dims des variables : (lon, lat, depth, time)
#    start/count : on prend tout en lon/lat/depth, 1 seul pas de temps
# ------------------------------------------------------------
start_vec <- c(1, 1, 1, time_idx)
count_vec <- c(nlon, nlat, ndepth, 1)

thetao <- ncvar_get(ds, "thetao", start = start_vec, count = count_vec)
so     <- ncvar_get(ds, "so",     start = start_vec, count = count_vec)

nc_close(ds)

# thetao_date et so_date ont maintenant les dimensions (lon, lat, depth)
dim(thetao)
dim(so)

# ---------------------------------------------------- applatir les mats de temp et sal
dim_ds <- dim(thetao)

nlon   <- dim_ds[1]
nlat   <- dim_ds[2]
ndepth   <- dim_ds[3]

thetao_flat <- matrix(aperm(thetao, c(1,2,3)), ncol = ndepth)
so_flat <- matrix(aperm(so, c(1,2,3)), ncol = ndepth)
dim(thetao_flat)

# -----------------------------------------------------filtrage nans
sum(rowSums(is.na(thetao_flat)) > 0)/nrow(thetao_flat)
sum(rowSums(is.na(so_flat)) > 0)/nrow(so_flat) 

# mask
mask_common <- !apply(is.na(thetao_flat), 1, any) &
  !apply(is.na(so_flat), 1, any)
sum(!mask_common)/length(mask_common)

# apply mask
thetao_flat_clean <-thetao_flat[mask_common,]
so_flat_clean <- so_flat[mask_common,]
sum(is.na(so_flat_clean))

# rm(thetao_flat, so_flat)
gc()
rm(so, thetao)


### B-spline 

K <- 25 # nombre de fonction de base
r <- 4 # ordre de la fonction de base
lambda_spline <- 0.25 # poids de la pénalisation

length(depth)
# base de Bsplines
phi <- bs(
  depth,
  df = K,
  degree = r-1,
  intercept = TRUE
)

phi <- as.matrix(phi)

# penalisation
D2 <- diff(diag(K), differences = 2)
R <- t(D2) %*% D2

# résolution 
A <- crossprod(phi) + lambda_spline * R
qr(A)$rank
B <- solve(A, t(phi))

coef_thetao <- thetao_flat_clean %*% t(B)
coef_so <- so_flat_clean %*% t(B)

# plot B-spline 
library(splines)
library(ggplot2)
library(patchwork)

library(splines)
library(ggplot2)
library(patchwork)

# ------------------------------------------------------------
# Choix d'un profil a illustrer (un pixel/une colonne d'eau)
# ------------------------------------------------------------
library(splines)
library(ggplot2)
library(patchwork)

# ------------------------------------------------------------
# Coordonnees lon/lat associees a chaque ligne de thetao_flat / so_flat
# (meme ordre que l'aplatissement : lon varie le plus vite, puis lat)
# ------------------------------------------------------------
coords_all   <- expand.grid(lon = lon, lat = lat)   # nrow = nlon * nlat, lon plus rapide
coords_clean <- coords_all[mask_common, ]           # meme masque que thetao_flat_clean / so_flat_clean

# ------------------------------------------------------------
# Choix d'un profil a illustrer (un pixel/une colonne d'eau)
# ------------------------------------------------------------
idx <- 1   # <-- modifiez l'indice pour tester un autre profil

lon_idx <- coords_clean$lon[idx]
lat_idx <- coords_clean$lat[idx]

profile_thetao <- thetao_flat_clean[idx, ]
profile_so     <- so_flat_clean[idx, ]

recon_thetao <- as.vector(phi %*% coef_thetao[idx, ])
recon_so     <- as.vector(phi %*% coef_so[idx, ])

rmse_thetao <- sqrt(mean((profile_thetao - recon_thetao)^2))
rmse_so     <- sqrt(mean((profile_so     - recon_so)^2))

# ------------------------------------------------------------
# Mise en forme longue : profil brut vs reconstruction B-spline
# ------------------------------------------------------------
df_temp <- data.frame(
  depth = rep(depth, 2),
  value = c(profile_thetao, recon_thetao),
  type  = rep(c("Profil brut", "Reconstruction B-spline"), each = length(depth))
)

df_sal <- data.frame(
  depth = rep(depth, 2),
  value = c(profile_so, recon_so),
  type  = rep(c("Profil brut", "Reconstruction B-spline"), each = length(depth))
)

# ------------------------------------------------------------
# Fonction pour construire un panel avec les 2 courbes + RMSE
# ------------------------------------------------------------
make_profile_plot <- function(df, var_title, x_label, rmse_val) {
  ggplot(df, aes(x = value, y = depth, color = type)) +
    geom_path(linewidth = 0.8) +
    geom_point(size = 0.9) +
    scale_y_reverse() +
    scale_color_manual(values = c("Profil brut" = "black",
                                  "Reconstruction B-spline" = "firebrick"),
                       name = NULL) +
    annotate("label", x = min(df$value, na.rm = TRUE),
             y = min(df$depth, na.rm = TRUE),
             label = sprintf("RMSE = %.3f", rmse_val),
             hjust = 0, vjust = 0, size = 3.5,
             fill = "white", label.size = 0.2) +
    labs(title = var_title, x = x_label, y = "Profondeur (m)") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")
}

p_temp <- make_profile_plot(df_temp, "Temperature", "Temperature (°C)", rmse_thetao)
p_sal  <- make_profile_plot(df_sal,  "Salinité",    "Salinité (PSU)",      rmse_so)

# ------------------------------------------------------------
# Assemblage final : Temp | Sal, cote a cote
# ------------------------------------------------------------
p_final <- (p_temp | p_sal) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = sprintf("Profil brut vs reconstruction B-spline (K=%d, degre=%d, lambda=%.2f) - 2023-01-26 - lon=%.2f, lat=%.2f",
                    K, r - 1, lambda_spline, lon_idx, lat_idx)
  ) &
  theme(legend.position = "bottom")

print(p_final)

ggsave(file.path(out_dir, sprintf("profil_bspline_temp_sal_lon%.2f_lat%.2f.png", lon_idx, lat_idx)),
       p_final, width = 10, height = 14, dpi = 300)

# PCA
coef_biv <- cbind(
  coef_thetao,
  coef_so
)
dim(coef_biv)

# centrage
alpha_mean <- colMeans(coef_biv)

Xc <- scale(
  coef_biv,
  center = alpha_mean,
  scale = FALSE
)

# covariance des coefficients
V <- crossprod(Xc) / (nrow(Xc)-1)

# décomposition spectrale
eig <- eigen(V, symmetric = TRUE)

# tri décroissant
ord <- order(eig$values, decreasing = TRUE)
lambda <- eig$values[ord] # val propres
U <- eig$vectors[, ord] # vect propres

prop_var <- lambda / sum(lambda) # explained var
cum_var <- cumsum(prop_var) # cum explained var
# nharm <- which(cum_var >= 0.90)[1] # n PC explaining 90% var

nharm <- 2
U <- U[,1:nharm]
lambda <- lambda[1:nharm]

# scores FPCA
scores <- Xc %*% U # (n, n_harm)

# reconstruction des coefficients de B spline
alpha_rec <- sweep(scores %*% t(U), 2, alpha_mean, "+") # reconstruction des scores de B splines
coef_temp_rec <- alpha_rec[,1:K] # (n_data, K=25)
coef_sal_rec  <- alpha_rec[,(K+1):(2*K)] # (n_data, K=25)

# plot pca 
library(splines)
library(ggplot2)
library(patchwork)

# ------------------------------------------------------------
# Coordonnees lon/lat associees a chaque ligne (deja calculees plus haut)
# ------------------------------------------------------------
# coords_all / coords_clean deja definis dans le script precedent

# ------------------------------------------------------------
# Choix d'un profil a illustrer
# ------------------------------------------------------------
idx <- 1   # <-- modifiez l'indice pour tester un autre profil

lon_idx <- coords_clean$lon[idx]
lat_idx <- coords_clean$lat[idx]

# profils bruts
profile_thetao <- thetao_flat_clean[idx, ]
profile_so     <- so_flat_clean[idx, ]

# reconstruction B-spline seule (deja calculee avant la PCA)
recon_bspline_thetao <- as.vector(phi %*% coef_thetao[idx, ])
recon_bspline_so     <- as.vector(phi %*% coef_so[idx, ])

# reconstruction B-spline + PCA (coefficients reconstruits apres reduction PCA)
recon_pca_thetao <- as.vector(phi %*% coef_temp_rec[idx, ])
recon_pca_so     <- as.vector(phi %*% coef_sal_rec[idx, ])

# RMSE : brut vs reconstruction PCA (celle qui nous interesse ici)
rmse_thetao <- sqrt(mean((profile_thetao - recon_pca_thetao)^2))
rmse_so     <- sqrt(mean((profile_so     - recon_pca_so)^2))

# ------------------------------------------------------------
# Mise en forme longue : 3 courbes par variable
# ------------------------------------------------------------
df_temp <- data.frame(
  depth = rep(depth, 3),
  value = c(profile_thetao, recon_bspline_thetao, recon_pca_thetao),
  type  = rep(c("Profil brut", "Reconstruction B-spline", "Reconstruction PCA"),
              each = length(depth))
)

df_sal <- data.frame(
  depth = rep(depth, 3),
  value = c(profile_so, recon_bspline_so, recon_pca_so),
  type  = rep(c("Profil brut", "Reconstruction B-spline", "Reconstruction PCA"),
              each = length(depth))
)

# ------------------------------------------------------------
# Fonction pour construire un panel avec les 3 courbes + RMSE (PCA)
# ------------------------------------------------------------
make_profile_plot <- function(df, var_title, x_label, rmse_val) {
  ggplot(df, aes(x = value, y = depth, color = type)) +
    geom_path(linewidth = 0.8) +
    geom_point(size = 0.9) +
    scale_y_reverse() +
    scale_color_manual(values = c("Profil brut"              = "black",
                                  "Reconstruction B-spline"  = "grey50",
                                  "Reconstruction PCA"       = "firebrick"),
                       name = NULL) +
    annotate("label", x = min(df$value, na.rm = TRUE),
             y = min(df$depth, na.rm = TRUE),
             label = sprintf("RMSE (PCA) = %.3f", rmse_val),
             hjust = 0, vjust = 0, size = 3.5,
             fill = "white", label.size = 0.2) +
    labs(title = var_title, x = x_label, y = "Profondeur (m)") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")
}

p_temp <- make_profile_plot(df_temp, "Temperature", "Temperature (deg C)", rmse_thetao)
p_sal  <- make_profile_plot(df_sal,  "Salinite",    "Salinite (PSU)",      rmse_so)

# ------------------------------------------------------------
# Assemblage final : Temp | Sal, cote a cote
# ------------------------------------------------------------
p_final <- (p_temp | p_sal) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = sprintf("Profil brut vs reconstruction B-spline / PCA (nharm=%d) - lon=%.2f, lat=%.2f",
                    nharm, lon_idx, lat_idx)
  ) &
  theme(legend.position = "bottom")

print(p_final)

ggsave(file.path(out_dir, sprintf("profil_pca_temp_sal_lon%.2f_lat%.2f_nharm%.2f.png", lon_idx, lat_idx, nharm)),
       p_final, width = 6, height = 12, dpi = 300)
