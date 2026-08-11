# Description

# 01_NASC_filter_sv
      # From obsaustral ncdf transect data 2018, 2021, 2023 : 
      #   1) filter by channel, keep only 200kHz data
      #   2) cut data on depth dimension to remove most NA
      #   3) filter lon >40°E and lat <-30°N
      #   4) filter day/night 
      # => save intermediary rds file : 2018_day.rds, 2018_night.rds, 2021_day.rds,...

# 02_NASC_concat_sv_diurnal_period
      #   5) Concatenate each year of same period (i.e. 2018_day, 2021_day, 2023_day)
      # => save intermediary rds file : 2018_2021_2022_2023_day.rds, 2018_2021_2022_2023_night.rds

# 03_NASC_check_pigmeann_grid
      #   6) get pigmeann grid (check if pigmeann grid consistent on every year of data)
      #   7) mean profile by grid (for day and night data)
      # => save intermediary rds file : mean_pig_grid_2018_2021_2023_day.rds, ...

#   8) compute NASC for each mean profile 
# => save intermediary rds file : nasc_mean_pig_grid_2018_2021_2023_day.rds,..
