## BIVARIATE FPCA --------------------------------------------------------------
## 
## 
## 
## Script: Marius Molinet + Nadège Fonvieille
## Last update: 07/2026
## 


# ---- Clean working space
graphics.off()
rm(list = ls())
cat('\f')
gc()


# ---- Packages
library(fda)


# ---- Working directory
#path_project <- "C:/Users/mmolinet/THESE MIO/Projects/Regional_FOD_NF"
path_project <- "E:/THESE MIO/Projects/Regional_FOD_NF"
Zone <- "Ker-Arg"
Researcher <- "AGM"
filepath <- file.path(path_project,"data","Regional_FOD",paste0("FOD_",Zone,"_",Researcher))


# ---- Load data
depth <- as.numeric(readRDS(file.path(filepath, "depth.RDS")))
lat_lon_dates <- readRDS(file.path(filepath, "lat_lon_dates.RDS"))
profiles <- readRDS(file.path(filepath, "profiles.RDS"))

variable_name <- c("thetao", "so")
n_profiles <- dim(profiles)[1]
nvar <- 2
min_depth <- min(depth)
max_depth <- max(depth)
n_depth <- length(depth)


# ---- B-spline basis projection 
K <- 10 # Number of basis functions
lambda <- 0.05 # Penalization parameter

# Cubic B-spline basis (degree 3)
basis_d3 <- create.bspline.basis(rangeval = range(depth), 
                                 nbasis = K, norder = 4)

# Evaluate basis functions at regular_depth
basismat <- eval.basis(depth, basis_d3)

# Transform profiles
fdnames <- as.list(c("depth", "profiles", variable_name))
coefs <- array(NA, dim = c(K, n_profiles, nvar))
for (i in 1:nvar){
  yfd <- Data2fd(argvals = depth, 
                 y = t(profiles[,,i]), 
                 basisobj = basis_d3, lambda = lambda,
                 fdnames = fdnames)
  coefs[,,i] <- yfd$coefs
  rm(yfd)
}
yfd <- fd(coef = coefs, basisobj = basis_d3, fdnames = fdnames)

# Visualization
s <- sample(1:dim(profiles)[1], 5)
par(mfrow = c(1,nvar), mar = c(4,4,1,1))
idxdepth <- n_depth:1
idxcoefs <- K:1
for(i in 1:nvar){
  irange <- range(profiles[s,,i], na.rm = T)
  matplot(t(profiles[s,,i]), -depth, pch = 1, 
          xlab = variable_name[i], las = 1, ylab = "", xlim = irange, ylim = c(-max_depth, -min_depth))
  matlines(basismat[,] %*% coefs[,s,i], -depth, lty = 1)
}


# Check dimensions
assertthat::are_equal(dim(yfd$coefs)[2], dim(profiles)[1])
assertthat::are_equal(dim(yfd$coefs)[1], K)
assertthat::are_equal(dim(yfd$coefs)[3], nvar)
# If false => there is a problem

# Save functional elements
fd.elements <- list(yfd = yfd, basis = basis_d3, nbasis = K, lambda = lambda, 
                    depth = depth, variable_name = variable_name)
saveRDS(fd.elements, file.path(filepath, "fd.elements.RDS"))




# ---- mfPCA
# Matrix X: merged coefficients
X_full <- NULL
for(i in 1:nvar){
  X_full <- cbind(X_full, t(coefs[,,i]))
}

# Remove NA profiles
NA_profiles <- which(is.na(apply(X_full, 1, sum)))
X <- X_full[-NA_profiles,]

# Number of observations
N <- nrow(X)
gc()

# Mean coefficients
alpha_chap <- apply(X, 2, mean, na.rm = T) 

# Matrix C: centered coefficients
C <- sweep(X, 2, alpha_chap, "-")
gc()

# -- Metric W by steps
# W_i: inner products of basis functions
basis_obj <- fd.elements$basis
W_i <- eval.penalty(basis_obj)
# By bloc
W <- matrix(0, K*nvar, K*nvar)
for (i in 1:nvar) {
  i0 <- i * K - K + 1
  i1 <- i * K
  W[i0:i1,i0:i1] <- W_i
}
# To ensure symmetry
W <- (W + t(W))/2
# Cholesky decomposition: W^1/2
Wdem <- chol(W)
# W^-1/2
Wdeminv <- solve(Wdem)

# -- Metric M by steps
# Sigma2 for each variable
sigma2_i <- NULL
for(i in 1:nvar){
  i0 <- i * K - K + 1
  i1 <- i * K
  V_ii <- 1/N * t(C[,i0:i1]) %*% C[,i0:i1] %*% W_i
  sigma2_i <- c(sigma2_i, sum(diag(V_ii)))
}
# Diagonal
M <- NULL
for(i in 1:nvar) {
  M <- c(M, rep(1/sigma2_i[i], K))
}
# Final matrix
M <- diag(M)
# M^1/2
Mdem <- sqrt(M)
# M^-1/2
Mdeminv <- solve(Mdem)

# --- VWM and eigen decomposition (V = 1/N * t(C) %*% C)
VWM <- 1/N * Mdem %*% Wdem %*% crossprod(C) %*% t(Wdem) %*% Mdem
gc()

mfpca <- eigen(VWM)

# Deal with complex numbers if needed
if(!is.null(Im(mfpca$values))){
  mfpca$values <- Re(mfpca$values)
  mfpca$vectors <- Re(mfpca$vectors)
}
# Deal with negative eigen values
mfpca$values <- abs(mfpca$values)
# Eigen vectors b_k
mfpca$vectnotWM <- mfpca$vectors
# W-normalized eigen vectors B_k
mfpca$vectors <- Mdeminv %*% Wdeminv %*% mfpca$vectnotWM
# Deformation associated to eigen function 
mfpca$axe <- sweep(mfpca$vectors, 2, sqrt(mfpca$values), "*")
# Principal components
mfpca$pc <- C %*% W %*% M %*% mfpca$vectors
# Percentage of inertia per fpca axes
mfpca$pval <- round(mfpca$values/sum(mfpca$values)*100,3)

# Save other elements
mfpca$alpha_chap <- alpha_chap
mfpca$W <- W
mfpca$M <- M

s <- sample(1:nrow(C), 100000)
plot(mfpca$pc[s,1],mfpca$pc[s,2])

# Save list
saveRDS(mfpca, file.path(filepath, "mfpca.RDS"))



