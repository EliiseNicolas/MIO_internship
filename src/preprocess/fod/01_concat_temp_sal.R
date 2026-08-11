# Description 
# from netcdf temperature-salinity CMEMS product :
#   1) Filter thetao and so by latitude (<-30°N), longitude(>40°E) and depth (0<depth<500m)
#   2) Flatten data and save into rds files of shape (n_depths, n_profiles) with meta data associated (array of shape (n_profiles, lat, lon, time)
#   3) Bspline estimation for every profile
#   4) Keep only 6PCs
#   5) bivariate fpca

library(ncdf4)

session <- "mmolinet"

path2018 <- paste0(
  "/run/media/", session,
  "/KER22/data_elise/raw/temperature_salinite/",
  "cmems_mod_glo_phy_my_0.083deg_P1D-m_thetao-so_",
  "40.00E-95.00E_60.00S-20.00S_0.49-5727.92m_",
  "2018-01-09-2018-03-03.nc"
)

path2021 <- paste0(
  "/run/media/", session,
  "/KER22/data_elise/raw/temperature_salinite/",
  "cmems_mod_glo_phy_my_0.083deg_P1D-m_thetao-so_",
  "40.00E-95.00E_60.00S-20.00S_0.49-5727.92m_",
  "2021-01-09-2021-03-03.nc"
)

path2022 <- paste0(
  "/run/media/", session,
  "/KER22/data_elise/raw/temperature_salinite/",
  "cmems_mod_glo_phy_my_0.083deg_P1D-m_thetao-so_",
  "40.00E-95.00E_60.00S-20.00S_0.49-5727.92m_",
  "2022-01-09-2022-03-03.nc"
)

path2023 <- paste0(
  "/run/media/", session,
  "/KER22/data_elise/raw/temperature_salinite/",
  "cmems_mod_glo_phy_my_0.083deg_P1D-m_thetao-so_",
  "40.00E-95.00E_60.00S-20.00S_0.49-5727.92m_",
  "2023-01-09-2023-03-03.nc"
)

paths <- c(
  path2018,
  path2021,
  path2022,
  path2023
)

years <- c(2018, 2021, 2022, 2023)


# ---- Function ---------------------------------------------------------

open_ds <- function(path, year, min_depth = 1, max_depth = 500) {
  
  cat("\nProcessing:", year, "\n")
  
  ds <- nc_open(path)
  
  # Coordinates
  lat_all <- ds$dim$latitude$vals
  lon_all <- ds$dim$longitude$vals
  depth_all <- ds$dim$depth$vals
  time <- ds$dim$time$vals
  
  # Spatial selection
  idx_lat <- which(lat_all < -30)
  idx_lon <- which(lon_all > 40)
  
  lat <- lat_all[idx_lat]
  lon <- lon_all[idx_lon]
  
  # Depth selection
  idx_depth <- which(
    depth_all >= min_depth &
      depth_all <= max_depth
  )
  
  depth <- depth_all[idx_depth]
  
  # Dimensions
  n_lon <- length(idx_lon)
  n_lat <- length(idx_lat)
  n_depth <- length(idx_depth)
  n_time <- length(time)
  
  cat(
    "Grid:",
    n_lon, "lon ×",
    n_lat, "lat ×",
    n_depth, "depth ×",
    n_time, "time\n"
  )
  
  # ---- Read thetao --------------------------------------------------
  
  thetao <- ncvar_get(
    ds,
    "thetao",
    start = c(
      min(idx_lon),
      min(idx_lat),
      min(idx_depth),
      1
    ),
    count = c(
      n_lon,
      n_lat,
      n_depth,
      n_time
    )
  )
  
  # ---- Read so ------------------------------------------------------
  
  so <- ncvar_get(
    ds,
    "so",
    start = c(
      min(idx_lon),
      min(idx_lat),
      min(idx_depth),
      1
    ),
    count = c(
      n_lon,
      n_lat,
      n_depth,
      n_time
    )
  )
  
  nc_close(ds)
  
  # ---- Dates --------------------------------------------------------
  
  dates <- as.POSIXct(
    time * 3600,
    origin = "1950-01-01",
    tz = "UTC"
  )
  
  # ---- Flatten thetao ----------------------------------------------
  
  thetao_mat <- matrix(
    aperm(thetao, c(1, 2, 4, 3)),
    ncol = n_depth
  )
  
  rm(thetao)
  gc()
  
  # ---- Flatten so ---------------------------------------------------
  
  so_mat <- matrix(
    aperm(so, c(1, 2, 4, 3)),
    ncol = n_depth
  )
  
  rm(so)
  gc()
  
  # ---- Coordinates --------------------------------------------------
  
  lat_lon_dates <- expand.grid(
    lon = lon,
    lat = lat,
    dates = dates,
    KEEP.OUT.ATTRS = FALSE
  )
  
  # ---- Checks -------------------------------------------------------
  
  stopifnot(
    nrow(thetao_mat) == nrow(so_mat),
    nrow(thetao_mat) == nrow(lat_lon_dates)
  )
  
  # ---- Return -------------------------------------------------------
  
  list(
    year = year,
    thetao = thetao_mat,
    so = so_mat,
    lat_lon_dates = lat_lon_dates,
    depth = depth
  )
}


# ---- Output directory ------------------------------------------------

out_dir <- "/run/media/mmolinet/KER22/MIO_internship_III/data_preprocessed/concat_temp_sal"


# ---- Process year by year --------------------------------------------

for (i in seq_along(paths)) {
  
  year <- years[i]
  
  ds <- open_ds(
    paths[i],
    year,
    min_depth = 1,
    max_depth = 500
  )
  
  # Save directly
  saveRDS(
    ds,
    file.path(
      out_dir,
      paste0("profiles_temp_sal_", year, ".rds")
    )
  )
  
  # Free memory
  rm(ds)
  gc()
  
  cat("Finished:", year, "\n")
}





