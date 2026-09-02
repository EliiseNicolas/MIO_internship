# ============================================================
# Recherche de journées avec faible couverture nuageuse
# (= forte disponibilité de données pigments)
# ============================================================

rm(list = ls())

# ------------------------------------------------------------
# Seuil de disponibilité recherché
# ------------------------------------------------------------

seuil_valid <- 0.70

# ------------------------------------------------------------
# Path -- fichier pigments NATIF (haute résolution, 1080x720)
# ------------------------------------------------------------

path_pig <- "F:/data_elise/pigmeann/pigments_2018_2021_2022_2023_crop.rds"

pigs_grid <- readRDS(path_pig)

list_pigs <- c("Chla", "Per", "But", "Fuco", "Hex", "Allo", "Zea", "Chlb", "DvChla")

n_dates <- length(pigs_grid$date)
n_pixels <- length(pigs_grid$lon) * length(pigs_grid$lat)

cat("Nombre de dates disponibles :", n_dates, "\n")
cat("Nombre de pixels par carte :", n_pixels, "\n")

# ------------------------------------------------------------
# Calcul du % de pixels valides pour Chla, pour chaque date
# ------------------------------------------------------------
# Chla est utilisé comme proxy principal (les ratios pigment/Chla
# nécessitent de toute façon Chla non-NA). Adapter si besoin.

pct_valid_chla <- sapply(seq_len(n_dates), function(i) {
  day_slice <- pigs_grid$c_cond_Chla[i, , ]
  sum(!is.na(day_slice)) / n_pixels
})

results <- data.frame(
  date = pigs_grid$date,
  pct_valid_chla = pct_valid_chla
)

results <- results[order(-results$pct_valid_chla), ]

cat("\nTop 10 des meilleures journées (Chla) :\n")
print(head(results, 10))

# ------------------------------------------------------------
# Version stricte : tous les pigments valides (pas juste Chla)
# ------------------------------------------------------------

pct_valid_all_pigs <- sapply(seq_len(n_dates), function(i) {
  valid <- rep(TRUE, n_pixels)
  for (p in list_pigs) {
    day_slice <- pigs_grid[[paste0("c_cond_", p)]][i, , ]
    valid <- valid & !is.na(as.vector(day_slice))
  }
  sum(valid) / n_pixels
})

results$pct_valid_all_pigs <- pct_valid_all_pigs
results <- results[order(-results$pct_valid_all_pigs), ]

cat("\nTop 10 des meilleures journées (tous pigments) :\n")
print(head(results, 10))

# ------------------------------------------------------------
# Première date qui dépasse le seuil (tous pigments)
# ------------------------------------------------------------

date_ok <- results$date[results$pct_valid_all_pigs >= seuil_valid]

if (length(date_ok) > 0) {
  cat("\n", length(date_ok), "date(s) dépassent", seuil_valid * 100, "% (tous pigments) :\n")
  print(sort(date_ok))
} else {
  cat("\nAucune date ne dépasse", seuil_valid * 100, "% avec tous les pigments valides.\n")
  cat("Meilleure date disponible :", format(results$date[1], "%Y-%m-%d"),
      "(", round(100 * results$pct_valid_all_pigs[1], 1), "%)\n")
}

# ------------------------------------------------------------
# Plot de la disponibilité au cours du temps
# ------------------------------------------------------------

library(ggplot2)

ggplot(results, aes(x = date, y = pct_valid_all_pigs)) +
  geom_line() +
  geom_point(size = 0.8) +
  geom_hline(yintercept = seuil_valid, linetype = "dashed", color = "red") +
  theme_bw() +
  labs(
    x = "Date",
    y = "% de pixels valides (tous pigments)",
    title = "Disponibilité des données pigments au cours du temps",
    subtitle = paste("Seuil recherché :", seuil_valid * 100, "%")
  )