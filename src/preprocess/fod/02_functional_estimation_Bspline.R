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


# rm(list=ls())

# Packages
library(fda)

# Hyperparameters
year <- 2018
K <- 10
lambda <- 0.05
order <- 4
path <- paste0("/run/media/mmolinet/KER22/MIO_internship_III/data_preprocessed/concat_temp_sal/profiles_temp_sal_", year, ".rds")

# Bspline estimation
bspline_estim <- function(path, K, lambda){
  # load data 
  df <- readRDS(path)
  n_profiles <- ncol(df$thetao)
  
  # create Bspline basis
  basis <- create.bspline.basis(rangeval = range(df$depth), nbasis = K, norder = 4)
  
  # Bspline estimation thetao
  yfd <- Data2fd(
    argvals = depth,
    y = t(df$thetao),
    basisobj = basis,
    lambda = lambda,
    fdnames = c("depth", "profiles", "thetao")
  )
  coefs_thetao <- yfd$coefs
  
  # Bspline estimation so
  yfd <- Data2fd(
    argvals = depth,
    y = t(df$so),
    basisobj = basis,
    lambda = lambda,
    fdnames = c("depth", "profiles", "so")
  )
  coefs_so <- array(NA, dim = c(K, n_profiles))
}

bspline_estim(path, K, lambda)
