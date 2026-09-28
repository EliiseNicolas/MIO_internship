library(ggplot2)
library(dplyr)
library(zoo)  # pour rollmean/rollmedian

ds <- readRDS("/run/media/mmolinet/KER22/data_elise/processed/ds_NASC_pig_ftle_fod/ds_NASC_per_esu_all/ds_NASC_per_esu_pig_ftle_fod_2018_2021_2023_transect_38kHz_mask9.rds")

ds_day <- ds %>%
  filter(time_nasc >= as.POSIXct("2023-01-28 18:00:00", tz = "UTC"),
         time_nasc <  as.POSIXct("2023-01-29 00:00:00", tz = "UTC")) %>%
  arrange(time_nasc)

# --- Option 1 : moyenne mobile ROBUSTE (médiane glissante, moins sensible aux extrêmes) ---
ds_day <- ds_day %>%
  mutate(nasc_roll = rollmedian(nasc, k = 21, fill = NA, align = "center"))
# k = taille de la fenêtre (nombre de points), à ajuster selon la densité de tes données
# rollmedian est plus robuste aux outliers que rollmean

p1 <- ggplot(ds_day, aes(x = time_nasc)) +
  geom_point(aes(y = nasc), size = 0.4, alpha = 0.2, color = "grey20") +
  geom_line(aes(y = nasc_roll), color = "steelblue", linewidth = 1) +
  labs(
    title = "NASC observé le 2023-01-28 (18h-00h)",
    x = "Temps", y = "NASC"
  ) +
  scale_x_datetime(date_labels = "%H:%M") +
  theme_minimal()

print(p1)

# --- Option 2 : lissage LOESS (courbe continue, plus "esthétique") ---
p2 <- ggplot(ds_day, aes(x = time_nasc, y = nasc)) +
  geom_point(size = 0.4, alpha = 0.2, color = "grey20") +
  geom_smooth(method = "loess", span = 0.1, se = FALSE, color = "darkred", linewidth = 1) +
  labs(
    title = "NASC observé le 2023-01-28 (18h-00h)",
    x = "Temps", y = "NASC"
  ) +
  scale_x_datetime(date_labels = "%H:%M") +
  theme_minimal()

print(p2)

ggsave("~/Elisou/mio_final_rep/MIO_internship/figures/figures_add_report/nasc_2023-01-28_smooth_bis.png", p1, width = 10, height = 5, dpi = 300)
