# Description 

# from NASC dataset (computed with 05_compute_nasc_mean_profile_pig_grid.R),
# NASC dataset of shape (n_ESU, 4) with 4 columns time, lon, lat, NASC
# 
# So from this dataset, we want to associate FTLE data that matches time, lat, long
# and create a new ds : nasc_ftle_pig_grid_day_200kHz.rds containing (n_ESU, 5) 
# with 5 columns time, lon, lat, NASC, FTLE

# global variables
path_nasc <- "~/Documents/stage_MIO/pt_III/data_preprocessed/NASC/transect_2018_2022_2023/NASC_mean_pig_grid_by_year_2018_2021_2023_day_200kHz.rds"
path_ftle_folder <- 