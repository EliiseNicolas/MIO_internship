## FOD MODEL AND IN SITU PROJECTION --------------------------------------------
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
library(mclust)
library(PCDimension)
library(dplyr)
library(ggplot2)
library(viridis)

# ---- Working directory
#path_project <- "C:/Users/mmolinet/THESE MIO/Projects/Regional_FOD_NF"
path_project <- "E:/THESE MIO/Projects/Regional_FOD_NF"
Zone <- "Ker-Arg"
Researcher <- "AGM"
filepath <- file.path(path_project,"data","Regional_FOD",paste0("FOD_",Zone,"_",Researcher))


# ---- Load data 
lat_lon_dates <- readRDS(file.path(filepath, "lat_lon_dates.RDS"))
fd.elements <- readRDS(file.path(filepath, "fd.elements.RDS"))
mfpca <- readRDS(file.path(filepath, "mfpca.RDS"))

# Functional elements
K <- fd.elements$nbasis
N <- nrow(mfpca$pc)

# ---- Step 1: Choose number of pc
# Eigenvalues
par(mfrow = c(1,1))
barplot(mfpca$pval)
cumsum(mfpca$pval)

# Here, I arbitrarily choose 3 PC as it reaches 98 % of total inertia
# 2 would be enough
npc <- 3

# You can also choose with objective methods (ex: broken stick model)
# https://cran.r-project.org/web/packages/PCDimension/vignettes/PCDimension.pdf
# bsDimension(mfpca$value)

# ----- Step 2: decide number of groups
# MBC is sensitive to initialization, so we execute 10 times the model 
# and choose G following ICL & BIC criterion
# If ICL and BIC are not conclusive, better to follow ICL criterion
# Takes a while to run

icl <- NULL
bic <- NULL
for(r in 1:5){
  # gc()
  print(r)
  
  set.seed(r)
  if(N > 10000){s <- sample(1:N, 10000)}else{s <- N}
  pc_to_cluster <- mfpca$pc[s,1:npc]
  
  # I choose VVV model as it gives more flexibility
  MBC.ICL <- mclustICL(pc_to_cluster, G = 1:10, modelNames = "VVV")
  MBC.BIC <- mclustBIC(pc_to_cluster, G = 1:10, modelNames = "VVV")
  icl <- rbind(icl, MBC.ICL[1:10,1])
  bic <- rbind(bic, MBC.BIC[1:10,1])
}

par(mfrow = c(1,1), mar = c(3,3,1,1))
matplot(t(icl), type = "b", lty = 1, pch = 1)
matplot(t(bic), type = "b", lty = 1, pch = 1)

saveRDS(icl, file.path(filepath, "FOD_ICL.RDS"))
saveRDS(bic, file.path(filepath, "FOD_BIC.RDS"))

# Here, I choose 5 groups 
G <- 6

# ---- Step 3: MBC
# Similarly, as MBC is sensitive to initialization, so we execute 10 times the model 
# and verify if the clustering is stable (even better to run 100 times)
# Notes: the color does not matter, we want to check stability in cluster's shape
# Takes a while to run

par(mfrow = c(3,2))
for(r in 1:10){
  # gc()
  print(r)
  
  set.seed(r)
  if(N > 10000){s <- sample(1:N, 10000)}else{s = N}
  pc_to_cluster <- mfpca$pc[s,1:npc]
  
  # I fixed VVV model as it gives more flexibility
  MBC <- Mclust(pc_to_cluster, G = G, modelNames = "VVV")
  plot(MBC$data, col = MBC$classification)
}

# Clustering with all data
MBC <- Mclust(mfpca$pc[,1:npc], G = G, modelNames = "VVV")
gc()
if(N > 20000){s <- sample(1:N, 20000)}else{s <- N}
plot(mfpca$pc[s,1:2], col = MBC$classification[s])

# Save results
saveRDS(MBC, file.path(filepath, paste0("MBC_", G, ".RDS")))


# ---- Check MBC attribution
classif <- MBC$classification

# Probability of attribution
proba <- apply(MBC$z, 1, max)

# Profiles with low probability of attribution
threshold <- 0.75
classif[proba < threshold] <- 0

# match with profiles idx
lat_lon_dates_FOD <- lat_lon_dates
lat_lon_dates_FOD$FOD <- NA
lat_lon_dates_FOD$FOD[-NA_profiles] <- classif


# Day zone map
days_vec <- sort(unique(lat_lon_dates_FOD$dates))
day <- days_vec[27]
zone <- "Ker"
col_mbc <- c("grey", viridis(n = MBC$G))

idx <- which(lat_lon_dates_FOD$dates == day & lat_lon_dates_FOD$zone == zone)
p <- ggplot(lat_lon_dates_FOD[idx, ], aes(x = lon, y = lat, fill = as.factor(FOD))) +
  geom_raster() +
  scale_fill_manual(breaks = as.factor(c(0:G)), values = col_mbc) +
  theme_minimal()

plotname <- paste0("FOD_map_",zone,"_", day, ".png")
ggsave(filename = plotname,
       plot = p,
       path = path_figure,
       width = 10,
       height = 7,
       dpi = 500)

# ONLY ONE TIME AFTER CHECKING !
# Dont forget to recode value to match FOD ordrer, do it after visualistion !!

# left = former ; right = new
classif <- recode(MBC$classification, `1` = 2, `2` = 3, `3` = 1, `4` = 5, `5` = 6, `6` = 4)
# then check before attributing
MBC$classification <- classif
# if you fail :
#classif <- apply(MBC$z, 1, function(x){which.max(x)})
# finally save
saveRDS(MBC, file.path(from, paste0("MBC_", G, ".RDS")))



lat_lon_dates_FOD$FOD <- NA
lat_lon_dates_FOD$FOD[-NA_profiles] <- classif
lat_lon_dates_FOD$proba <- NA
lat_lon_dates_FOD$proba[-NA_profiles] <- proba

# Saving lat_lon_dates with FOD
saveRDS(lat_lon_dates_FOD, file.path(from, paste0("lat_lon_dates_FOD_G", G, ".RDS")))



# ---- SES data selection
from_ses <- file.path(path_project,"data","ses",paste0("SES_",Zone,"_",Researcher))
SES_data <- readRDS(file.path(from_ses, "SES_data.RDS"))

# Number of individuals
nind <- length(SES_data)

ses_variable_name <- names(SES_data[[1]])[-1]
nvar <- length(ses_variable_name)

# Gather all ses profiles
ses_metadata <- NULL
ses_profiles <- replicate(nvar, NULL)
names(ses_profiles) <- ses_variable_name
for(n in 1:nind){
  unique_data <- unique(SES_data[[n]]$meta_data$date)
  ndays <- length(unique_data)
  ses_metadata <- rbind(ses_metadata, SES_data[[n]]$meta_data) # we remove dive above plateau after
  for(i in 1:nvar){
    ses_profiles[[i]] <- cbind(ses_profiles[[i]], SES_data[[n]][[ses_variable_name[i]]])
  }
}

# Check dimensions
assertthat::are_equal(nrow(ses_metadata), ncol(ses_profiles[[1]]))
assertthat::are_equal(nrow(ses_metadata), ncol(ses_profiles[[2]]))

# Select profiles reaching max_depth
# This option implies the exclusion of all profiles not reaching max_depth
# Check consistency between variables
depth_reached <- colSums(!is.na(ses_profiles$TEMPERATURE))
assertthat::are_equal(which(!is.na(ses_profiles[[1]][max_depth,])), which(!is.na(ses_profiles[[2]][max_depth,])))
dives_to_keep <- which(!is.na(ses_profiles[[1]][max_depth,] + ses_profiles[[2]][max_depth,]))

cat("SES Profiles considered: ",  length(dives_to_keep)/ncol(ses_profiles[[1]])*100, " %")

profiles <- replicate(nvar, NULL)
names(profiles) <- ses_variable_name
for(i in 1:nvar){
  profiles[[i]] <- ses_profiles[[i]][min_depth:max_depth, dives_to_keep]
}

# selecting profiles with no NA
metadata <- ses_metadata[dives_to_keep,]

# Depth of SES profiles
ses_depth <- min_depth:max_depth

# Visualize along track profiles considered
par(mfrow = c(1,1), mar = c(4,4,1,1))
plot(ses_metadata$lon, ses_metadata$lat, pch = 16, xlab = "longitude", ylab = "latitude")
points(metadata$lon, metadata$lat, col = "red2", pch = 16, cex = 0.4)
legend("topright", legend = "Dives considered", bty = "n", col = "red2", pch = 16)



# ---- SES profiles projection
profiles.fd <- replicate(nvar, NULL)
names(profiles.fd) <- ses_variable_name
for(i in 1:nvar){
  profiles.fd[[i]] <- Data2fd(ses_depth, profiles[[i]], basis_obj, lambda = lambda)
}

# Visualization
s <- sample(1:ncol(profiles[[i]]), 5)
par(mfrow = c(1,nvar), mar = c(4,4,1,1))
for(i in 1:nvar){
  matplot(profiles[[i]][,s], -ses_depth, pch = 1, col = "grey80",
          xlab = ses_variable_name[i], las = 1, ylab = "")
  matlines(basismat %*% profiles.fd[[i]]$coefs[,s], -depth, 
           lty = 1, col = 1)
}

# Matrix Delta of new coefficients
Delta <- NULL
for(i in 1:nvar){
  Delta <-  cbind(Delta, t(profiles.fd[[i]]$coefs))
}
# Centered matrix
Delta_c <-  sweep(Delta, 2, alpha_chap, "-")
# Matrix P of Principal Coordinates
P <- Delta_c %*% W %*% M %*% mfpca$vectors

# Number of PC considered
npc <- ncol(MBC$data)
prediction <- predict(MBC, P[,c(1:npc)])

# check with following plots before recoding
prediction$classification <- recode(prediction$classification, `1` = 2, `2` = 3, `3` = 1, `4` = 5, `5` = 6, `6` = 4)


# Save data 
results <- list(ses_metadata = ses_metadata, metadata = metadata, 
                matrix_P = P, FOD_prediction = prediction)
saveRDS(results, file.path(from_ses, paste0("G",G,"_SES_prediction.RDS")))



# Plotting
SES_FOD <- ses_metadata
SES_FOD$FOD <- NA
SES_FOD$FOD[dives_to_keep] <- prediction$classification
SES_FOD$lon <- as.numeric(SES_FOD$lon)
SES_FOD$lat <- as.numeric(SES_FOD$lat)
SES_FOD$zone <- "Ker"
SES_FOD$zone[SES_FOD$lon < 0] <- "Arg"

days_vec <- sort(unique(lat_lon_dates_FOD$dates))
day <- days_vec[27]
zone <- "Ker"
col_mbc <- c("grey", viridis(n = MBC$G))

idx <- which(lat_lon_dates_FOD$dates == day & lat_lon_dates_FOD$zone == zone)
p <- ggplot(lat_lon_dates_FOD[idx, ], aes(x = lon, y = lat, fill = as.factor(FOD))) +
  geom_raster() +
  geom_point(data = SES_FOD[SES_FOD$zone == zone, ], aes(x = lon, y = lat, col = as.factor(FOD)), inherit.aes = F) +
  scale_fill_manual(breaks = as.factor(c(0:G)), values = col_mbc) +
  scale_color_manual(breaks = as.factor(c(0:G)), values = col_mbc) +
  theme_minimal()

plotname <- paste0("SES_FOD_map_",zone,"_", day, ".png")
ggsave(filename = plotname,
       plot = p,
       path = path_figure,
       width = 10,
       height = 7,
       dpi = 500)



