# =====================================================================
# generate_fod_distribution_relabeled.R -- redessine la distribution
# des clusters FOD par fold (38 kHz, RF, schemas 1500x1000km et
# 200x200km) avec les VRAIS labels/couleurs de cluster, a partir de
# fod_dist.csv deja sur disque (aucun recalcul).
# =====================================================================
# A lancer depuis la racine de nasc_pipeline/ :
#   setwd("chemin/vers/arbres_regression_clean")
#   source("generate_fod_distribution_relabeled.R")
#
# Sorties : fod_distribution_relabeled_1500x1000km.png,
#           fod_distribution_relabeled_200x200km.png

source("R/00_config.R")
source("R/05_plots.R")

TARGET_FREQ <- 38
MODEL       <- "rf"
SCHEMES     <- c("blocked_spatial_1500x1000km", "blocked_spatial_200x200km")

cluster_cols <- c(
  "1" = "#2166AC",
  "2" = "#67A9CF",
  "3" = "#1A9850",
  "4" = "#A6D96A",
  "5" = "#FDAE61",
  "6" = "#D73027"
)
transition_cols <- c(
  "7"  = "#3B73B9",  # 1-2
  "8"  = "#3FA7B5",  # 1-3
  "9"  = "#45B97C",  # 2-3
  "10" = "#C46A00",  # 3-4
  "11" = "#E85D04",  # 4-6
  "12" = "#C9184A",  # 4-5
  "13" = "#8F1D3F"   # 5-6
)
fod_cols <- c(cluster_cols, transition_cols)

legend_codes  <- c(1, 7, 8, 2, 9, 3, 10, 4, 12, 11, 5, 13, 6)
legend_labels <- c(
  "C1", "T1-2", "T1-3",
  "C2", "T2-3",
  "C3", "T3-4",
  "C4", "T4-5", "T4-6",
  "C5", "T5-6", "C6"
)
names(legend_labels) <- as.character(legend_codes)

training_dir <- path_out("training")

plot_fod_relabeled <- function(scheme) {
  csv_path <- file.path(training_dir, paste0(TARGET_FREQ, "kHz"), MODEL, scheme, "fod_dist.csv")
  if (!file.exists(csv_path)) {
    cat("[!] introuvable :", csv_path, "\n")
    return(NULL)
  }
  fod_dist_all <- read.csv(csv_path, colClasses = c(fod = "character"))

  # Seul changement par rapport a plot_fod_distribution() (05_plots.R) :
  # le code numerique brut de fod est remplace par son label lisible
  # (C1, T1-2, ...), dans l'ordre de legend_codes -- tout le reste
  # (fill=set, geom_col dodge, facet_wrap, labs, theme) est identique.
  fod_dist_all <- fod_dist_all %>%
    dplyr::filter(fod %in% names(legend_labels)) %>%
    dplyr::mutate(fod = factor(legend_labels[fod], levels = legend_labels[as.character(legend_codes)]))

  ggplot(fod_dist_all, aes(x = fod, y = prop, fill = set)) +
    geom_col(position = "dodge") +
    facet_wrap(~fold_id) +
    labs(title = "Distribution des clusters FOD, train vs test, par fold",
         subtitle = sprintf("RF - %d kHz - %s", TARGET_FREQ, scheme),
         x = "Cluster FOD", y = "Proportion", fill = NULL) +
    theme_pipeline + theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

for (scheme in SCHEMES) {
  p <- plot_fod_relabeled(scheme)
  if (is.null(p)) next
  out_file <- paste0("fod_distribution_relabeled_", gsub("blocked_spatial_", "", scheme), ".png")
  ggsave(out_file, p, width = 12, height = 5, dpi = 150)
  cat("Ecrit :", out_file, "\n")
}
