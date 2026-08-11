## LOAD CMEMS FILES AND RESHAPE PROFILES INTO MATRIX ---------------------------
## 
## 
## 
## Script: Marius Molinet
## Last update: 07/2026
## 


# ---- Clean working space
graphics.off()
rm(list = ls())
cat('\f')
gc()


# ---- Packages
library(ncdf4)
library(stringr)
library(assertthat)

# ---- Working directory
#path_project <- "C:/Users/mmolinet/THESE MIO/Projects/Regional_FOD_NF"
path_project <- "E:/THESE MIO/Projects/Regional_FOD_NF"
Zone <- "Ker-Arg"
Researcher <- "AGM"
from <- file.path(path_project,"data","TS_profiles",paste0("FOD_",Zone,"_",Researcher))
to <- file.path(path_project,"data","Regional_FOD",paste0("FOD_",Zone,"_",Researcher))

# ---- Get files information 

# Files in folder
filelist <- list.files(path = from, pattern = "nc$", full.names = TRUE)

# Select depth
min_depth <- 20
max_depth <- 300


# ---- Load and Reshape data
so_profiles_list <- list()
thetao_profiles_list <- list()
lat_lon_dates_list <- list()
n_profiles <- 0
for (file in filelist){
  
  cat(file, '\n')
  
  ncid <- nc_open(file)
  lon <- ncvar_get(ncid, "longitude")
  lat <- ncvar_get(ncid, "latitude")
  dates <- as.POSIXct(ncvar_get(ncid, 'time')*3600, origin = "1950-01-01", tz = "UTC")
  depth <- ncvar_get(ncid, "depth")
  so <- ncvar_get(ncid, "so")
  thetao <- ncvar_get(ncid, "thetao")
  
  # selecting depth of interest
  idx_depth <- which(depth >= min_depth & depth <= max_depth)
  new_depth <- depth[idx_depth]
  n_new_depth <- length(new_depth)
  
  # Reshape array 
  thetao <- thetao[,,idx_depth,]
  thetao_mat <- matrix(aperm(thetao, c(1,2,4,3)), ncol = n_new_depth)
  so <- so[,,idx_depth,]
  so_mat <- matrix(aperm(so, c(1,2,4,3)), ncol = n_new_depth)
  lat_lon_dates <- expand.grid(lon = lon, lat = lat, dates = dates)
  n_profiles <- n_profiles + nrow(thetao_mat)
  
  # save in list
  thetao_profiles_list[[length(thetao_profiles_list) + 1]] <- thetao_mat
  so_profiles_list[[length(so_profiles_list) + 1]] <- so_mat
  lat_lon_dates_list[[length(lat_lon_dates_list) + 1]] <- lat_lon_dates
  
  rm(thetao, thetao_mat, so, so_mat, lat_lon_dates)
}

# Convert list to matrix
profiles <- array(NA, dim = c(n_profiles, n_new_depth, 2))
profiles[,,1] <- do.call(rbind, thetao_profiles_list)
profiles[,,2] <- do.call(rbind, so_profiles_list)

# Adding zone to lat_lon_dates
lat_lon_dates <- do.call(rbind, lat_lon_dates_list)
lat_lon_dates$zone <- NA
lat_lon_dates$zone[lat_lon_dates$lon < 0] <- "Arg"
lat_lon_dates$zone[lat_lon_dates$lon > 0] <- "Ker"

# clear memory
rm(lat_lon_dates_list, so_profiles_list, thetao_profiles_list)

# Save extracted profiles and associated variables
saveRDS(new_depth, file.path(to, "depth.RDS"))
saveRDS(lat_lon_dates, file.path(to, "lat_lon_dates.RDS"))
saveRDS(profiles, file.path(to, "profiles.RDS"))


# # check plot
# library(ggplot2)
# library(viridis)
# 
# sst <- profiles[,1,1]
# lat_lon_dates$sst <- sst
# 
# idx <- which(lat_lon_dates$dates == dates[19] & lat_lon_dates$zone == "Arg")
# 
# ggplot(lat_lon_dates[idx, ], aes(x = lon, y = lat, fill = sst)) +
#   geom_raster() +
#   scale_fill_viridis()




