setwd("C:/Users/mmolinet/elisou_ta_stagiaire_pref/MIO_internship_III/src/arbres_regression_clean")
source("R/00_config.R")
source("R/05_plots.R")

# Charge les 5 importance_all.csv (un par schema), pour RF/38kHz uniquement
schemes <- c("naive_RS_80_20",
             paste0("blocked_spatial_", c("1500x1000km", "200x200km", "20x20km")),
             "blocked_temporal_1j")

imp_rf_38 <- purrr::map_dfr(schemes, function(sc) {
  path <- file.path("outputs_pipeline/training/38kHz/rf", sc, "importance_all.csv")
  if (!file.exists(path)) { cat("[!] manquant :", path, "\n"); return(NULL) }
  read.csv(path) %>% dplyr::mutate(model = "rf", scheme = sc)
})

p <- plot_importance_comparison(imp_rf_38, subtitle = "RF - 38 kHz - tous schemas de CV")
print(p)

# Pour la sauvegarder :
ggsave("importance_rf_38kHz_tous_schemas.png", p, width = 9, height = 6, dpi = 150)