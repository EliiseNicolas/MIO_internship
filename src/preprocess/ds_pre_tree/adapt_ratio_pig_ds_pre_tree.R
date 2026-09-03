# Libraries
library(dplyr)

# global var
for (freq in c(38, 70, 120, 200)){
  ds <- readRDS(paste0("F:/data_elise/ds_NASC_pig_ftle_fod/ds_NASC_per_esu_all/ds_NASC_per_esu_pig_ftle_fod_2018_2021_2022_2023_transect_", freq, "kHz_mask9.rds"))
  
  str(ds)
  
  
  # ------------------------------------------------------------
  # Variables de pigments bruts (sans "_Chla" ni "_total")
  # ------------------------------------------------------------
  list_pigs <- c("Chla", "Per", "But", "Fuco", "Hex", "Allo", "Zea", "Chlb", "DvChla")
  
  # ------------------------------------------------------------
  # 1) Supprimer toutes les colonnes se terminant par "_Chla" ou "_total"
  #    (y compris total_pig si vous voulez repartir propre, sinon retirez-le
  #    du pattern ci-dessous)
  # ------------------------------------------------------------
  cols_to_drop <- grep("_Chla$|_total$", names(ds), value = TRUE)
  print(cols_to_drop)   # verification avant suppression
  
  ds <- ds %>% select(-all_of(cols_to_drop))
  
  # ------------------------------------------------------------
  # 2) chla_total = Chla seule
  # ------------------------------------------------------------
  ds$Chla_total <- ds$Chla
  
  # ------------------------------------------------------------
  # 3) Somme des pigments SAUF Chla (calculee une seule fois)
  # ------------------------------------------------------------
  pigs_sans_chla <- setdiff(list_pigs, "Chla")
  
  sum_others <- Reduce(`+`, ds[pigs_sans_chla])
  
  # ------------------------------------------------------------
  # 4) Ratio de chaque pigment (Chla incluse) sur la somme des
  #    pigments hors Chla, vectorise sans boucle
  # ------------------------------------------------------------
  ratio_df <- as.data.frame(
    lapply(list_pigs, function(p) ds[[p]] / sum_others)
  )
  names(ratio_df) <- paste0(list_pigs, "_totpig")
  
  ds <- cbind(ds, ratio_df)
  
  str(ds)
  
  saveRDS(ds, paste0("F:/data_elise/ds_NASC_pig_ftle_fod/ds_NASC_per_esu_all/new_ratio_ds_NASC_per_esu_pig_ftle_fod_2018_2021_2022_2023_transect_", freq, "kHz_mask9.rds"))
  
}
