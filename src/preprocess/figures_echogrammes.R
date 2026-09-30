df <- readRDS("F:/data_elise/sv_cropped/sv_cropped_per_year/120kHz/Sv_2023_120kHz.rds")

str(df)


library(ggplot2)
library(reshape2)
library(dplyr)

## -----------------------------------------------------------------
## 1. Choix de la journée à afficher
## -----------------------------------------------------------------
## Option A : par date (recommandé, plus explicite)
print(unique(as.Date(df$time)))
date_choisie <- as.Date("2023-01-28")
idx_jour <- which(as.Date(df$time) == date_choisie)

## Option B : par numéro de jour (décommenter si préféré)
# jour_choisi <- 3
# idx <- which(df$day == jour_choisi)

if (length(idx_jour) == 0) stop("Aucune donnée pour cette journée : vérifiez date_choisie / jour_choisi.")

## -----------------------------------------------------------------
## 1bis. Restriction à une plage horaire (ici 18h -> 00h)
## -----------------------------------------------------------------
heure_debut <- 18   # inclus
heure_fin   <- 24   # exclu (24 = minuit)

heure_num <- as.numeric(format(df$time[idx_jour], "%H")) +
  as.numeric(format(df$time[idx_jour], "%M")) / 60

idx <- idx_jour[heure_num >= heure_debut & heure_num < heure_fin]

if (length(idx) == 0) stop("Aucune donnée dans la plage horaire choisie.")

## -----------------------------------------------------------------
## 2. Extraction de la sous-matrice profils pour cette journée
## -----------------------------------------------------------------
mat       <- df$profiles[idx, , drop = FALSE]   # lignes = temps, colonnes = profondeur
time_sel  <- df$time[idx]
depth_vec <- df$depth

## Mise en forme longue pour ggplot
df_echo <- melt(mat)
names(df_echo) <- c("time_idx", "depth_idx", "Sv")
df_echo$time  <- time_sel[df_echo$time_idx]
df_echo$depth <- depth_vec[df_echo$depth_idx]

## Remplacer les valeurs non plausibles (souvent -999, -Inf...) par NA
df_echo$Sv[!is.finite(df_echo$Sv)] <- NA

## -----------------------------------------------------------------
## 3. Largeur/hauteur explicite des cases (corrige les bandes blanches
##    dues aux écarts irréguliers entre pings, sans inventer de données)
## -----------------------------------------------------------------
pas_temps  <- as.numeric(median(diff(sort(unique(df_echo$time))), na.rm = TRUE))
pas_depth  <- median(diff(sort(unique(depth_vec))), na.rm = TRUE)

## -----------------------------------------------------------------
## 4. Échogramme
## -----------------------------------------------------------------
sv_min <- -100  # à ajuster selon vos données (même échelle pour tous les jours)
sv_max <- -60

p <- ggplot(df_echo, aes(x = time, y = depth, fill = Sv)) +
  geom_tile(width = pas_temps, height = pas_depth) +
  scale_y_reverse(name = "Profondeur (m)") +
  scale_fill_gradientn(
    colours = c("#2c2c8f", "#1f9e89", "#7ad151", "#fde725", "#e31a1c"),
    name = "Sv (dB)",
    na.value = "white",
    limits = c(sv_min, sv_max),
    oob = scales::squish   # écrase les valeurs hors plage au lieu de les blanchir
  ) +
  scale_x_datetime(
    date_labels = "%H:%M",
    limits = range(time_sel),   # bornes prises directement dans les données (même fuseau horaire)
    expand = c(0, 0)
  ) +
  labs(
    title = paste0("Échogramme - 120 kHz - ", format(date_choisie, "%d/%m/%Y"),
                   " (", heure_debut, "h-", ifelse(heure_fin == 24, "00", heure_fin), "h)"),
    x = "Heure"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank())

print(p)
## Pour enregistrer la figure :
ggsave(paste0("C:/Users/mmolinet/elisou_ta_stagiaire_pref/MIO_internship_III/figures/figures_add_report/echogramme_120kHz_", date_choisie, ".png"), p, width = 12, height = 5, dpi = 300)
