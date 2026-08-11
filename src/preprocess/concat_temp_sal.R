library(ncdf4)

path2018 <- "/run/media/elise/KER22/data_elise/raw/temperature_salinite/cmems_mod_glo_phy_my_0.083deg_P1D-m_thetao-so_40.00E-95.00E_60.00S-20.00S_0.49-5727.92m_2018-01-09-2018-03-03.nc"
path2021 <- "/run/media/elise/KER22/data_elise/raw/temperature_salinite/cmems_mod_glo_phy_my_0.083deg_P1D-m_thetao-so_40.00E-95.00E_60.00S-20.00S_0.49-5727.92m_2021-01-09-2021-03-03.nc"
path2022 <- "/run/media/elise/KER22/data_elise/raw/temperature_salinite/cmems_mod_glo_phy_my_0.083deg_P1D-m_thetao-so_40.00E-95.00E_60.00S-20.00S_0.49-5727.92m_2022-01-09-2022-03-03.nc"
path2023 <- "/run/media/elise/KER22/data_elise/raw/temperature_salinite/cmems_mod_glo_phy_my_0.083deg_P1D-m_thetao-so_40.00E-95.00E_60.00S-20.00S_0.49-5727.92m_2023-01-09-2023-03-03.nc"

open_ds <- function(path, year){
  ds <- nc_open(path)
  thetao <- ncvar_get(ds, "thetao")
  so <- ncvar_get(ds, "so")
  lat <- ds$dim$latitude$vals
  lon <- ds$dim$longitude$vals
  time <- ds$dim$time$vals
  
  nc_close(ds)
  
  df <- list(
    year = year,
    thetao = thetao,
    so = so,
    lat = lat,
    lon = lon,
    time = time
  )
  
  return(df)
}
ds2018 <- open_ds(path2018, 2018)
str(ds2018)