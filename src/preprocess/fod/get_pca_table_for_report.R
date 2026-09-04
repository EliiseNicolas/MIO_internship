save_path <- "F:/data_elise/fod_elise_2018_2021_2022_2023/mclust_results_2018_2021_2022_2023"

coef_thetao <- readRDS(file.path(save_path, "coef_thetao.rds"))
coef_so     <- readRDS(file.path(save_path, "coef_so.rds"))

coef_biv <- cbind(coef_thetao, coef_so)

alpha_mean <- colMeans(coef_biv)
Xc <- scale(coef_biv, center = alpha_mean, scale = FALSE)

V   <- crossprod(Xc) / (nrow(Xc) - 1)
eig <- eigen(V, symmetric = TRUE)

ord    <- order(eig$values, decreasing = TRUE)
lambda <- eig$values[ord]

prop_var <- lambda / sum(lambda)
cum_var  <- cumsum(prop_var)

# Table de variance expliquee
variance_table <- data.frame(
  PC          = seq_along(lambda),
  eigenvalue  = lambda,
  prop_var    = prop_var,
  cum_var     = cum_var
)

head(variance_table, 10)