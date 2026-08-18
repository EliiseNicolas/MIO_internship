# Description 

# from netcdf temperature-salinity CMEMS product :
# 01_concat_temp_sal.R
    #   1) Filter thetao and so by latitude (<-30°N), longitude(>40°E) and depth (0<depth<500m)
    #   2) Flatten data and save into rds files of shape (n_depths, n_profiles) with meta data associated (array of shape (n_profiles, lat, lon, time)

# 02_functional_estimation_Bspline.R
#   3) Bspline estimation for every profile

# 03_fPCA.R
#   4) Keep only 6PCs
#   5) bivariate fpca


# Packages
library(fda)

# Hyperparameters
year <- 2018
K <- 25
lambda <- 0.25
order <- 4
path <- paste0("/Elisou/MIO_internship_III/data_preprocessed/concat_temp_sal/profiles_temp_sal_", year, ".rds")
path <- "~/Elisou/MIO_internship_III/data_preprocessed/concat_temp_sal/profiles_temp_sal_2018.rds"


# Bspline estimation

bspline_estim <- function(path, K, lambda){
  # load data 
  df <- readRDS(path)
  n_profiles <- nrow(df$thetao)
  print(n_profiles)
  # create Bspline basis
  basis <- create.bspline.basis(rangeval = range(df$depth), nbasis = K, norder = 4)

  # Bspline estimation thetao
  yfd_thetao <- Data2fd(
    argvals = as.numeric(df$depth),
    y = t(as.matrix(df$thetao)),
    basisobj = basis,
    lambda = lambda
  )
  coefs_thetao <- yfd_thetao$coefs
  thetao_reconstructed <- eval.fd(as.numeric(df$depth), yfd_thetao)

  # Bspline estimation so
  yfd_thetao <- Data2fd(
    argvals = as.numeric(df$depth),
    y = t(as.matrix(df$so)),
    basisobj = basis,
    lambda = lambda
  )
  coefs_so <- yfd_so$coefs
  so_reconstructed <- eval.fd(as.numeric(df$depth), yfd_so)


  # -------------------------
  # Plot thetao
  # -------------------------
  profile_id <- 1

  plot(
    df$thetao[profile_id, ],
    df$depth,
    type = "l",
    col = "red",
    lwd = 2,
    xlab = "Temperature (°C)",
    ylab = "Depth (m)",
    main = paste("Profile", profile_id, "- thetao"),
    ylim = rev(range(df$depth))
  )
  
  lines(
    thetao_reconstructed[, profile_id],
    df$depth,
    type = "l",
    col = "black",
    lwd = 2
  )
  
  legend(
    "topright",
    legend = c("Original", "Reconstructed"),
    col = c("red", "black"),
    lwd = 2
  )

  # -------------------------
  # Plot so
  # -------------------------

  plot(
    df$so[profile_id, ],
    df$depth,
    type = "l",
    col = "red",
    lwd = 2,
    xlab = "Salinity (PSU)",
    ylab = "Depth (m)",
    main = paste("Profile", profile_id, "- so"),
    ylim = rev(range(df$depth))
  )
  
  lines(
    so_reconstructed[, profile_id],
    df$depth,
    type = "l",
    col = "black",
    lwd = 2
  )
  
  legend(
    "topright",
    legend = c("Original", "Reconstructed"),
    col = c("red", "black"),
    lwd = 2
  )

  # Return results
  return(list(
    thetao_fd = yfd_thetao,
    so_fd = yfd_so,
    coefs_thetao = coefs_thetao,
    coefs_so = coefs_so,
    basis = basis,
    K=K,
    lambda=lambda
  ))
}

results <- bspline_estim(path, K, lambda)

# load data 
df <- readRDS(path)
n_profiles <- nrow(df$thetao)
print(n_profiles)


# create Bspline basis
K = 25
lambda = 0.25
basis <- create.bspline.basis(rangeval = range(df$depth), nbasis = K, norder = 4)

# Bspline estimation thetao
yfd_thetao <- Data2fd(
  argvals = as.numeric(df$depth),
  y = t(as.matrix(df$thetao[1:1000,])),
  basisobj = basis,
  lambda = lambda
)
coefs_thetao <- yfd_thetao$coefs
thetao_reconstructed <- eval.fd(as.numeric(df$depth), yfd_thetao)


profile_id <- 1

plot(
  df$thetao[profile_id, ],
  df$depth,
  type = "l",
  col = "red",
  lwd = 2,
  xlab = "Temperature (°C)",
  ylab = "Depth (m)",
  main = paste("Profile", profile_id, "- thetao"),
  ylim = rev(range(df$depth))
)

lines(
  thetao_reconstructed[, profile_id],
  df$depth,
  type = "l",
  col = "black",
  lwd = 2
)

legend(
  "topright",
  legend = c("Original", "Reconstructed"),
  col = c("red", "black"),
  lwd = 2
)
