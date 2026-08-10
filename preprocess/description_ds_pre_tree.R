# -----------  Description
# We want to understand what is the distribution of each entry variable to our tree, 
# and we want to know how many missing values we have per variable 

# -----------  Libraries
rm(list = ls())
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# -----------  Paths and open dataset
path <- "/run/media/mmolinet/KER22/données elise/elisou_ta_stagiaire_pref/prepross/regression_ds_2021_2022_2023_m100.rds"
ds <- readRDS(path)
print(nrow(ds))
str(ds)

# -----------  Distributions
# FOD
fod_vars <- paste0("fod", 0:6)

# rebuild one vector containing all fod number intead of 6 onehot vectors
ds_fod <- ds %>%
  mutate(
    fod_cluster = apply(select(., all_of(fod_vars)), 1, function(x) {
      if(all(is.na(x)) || sum(x) == 0) {
        return(NA)
      } else {
        return(which(x == 1) - 1)
      }
    })
  )

# replace NA by "NA" category
ds_fod$fod_cluster <- as.character(ds_fod$fod_cluster)
ds_fod$fod_cluster[is.na(ds_fod$fod_cluster)] <- "NA"
fod_distribution <- ds_fod %>%
  count(fod_cluster) %>%
  mutate(
    percentage = n / sum(n) * 100
  )

# histogram of distribution
ggplot(fod_distribution, aes(x = fod_cluster, y = percentage)) +
  geom_col() +
  labs(
    x = "FOD cluster",
    y = "Percentage (%)",
    title = "FOD Distribution on 2021, 2022, 2023 station datas"
  ) +
  theme_bw()

# FTLE
ggplot(ds, aes(x = ftle, y = after_stat(count / sum(count) * 100))) +
  geom_histogram(bins = 200) +
  theme_bw() +
  labs(
    x = "FTLE",
    y = "Pourcentage (%)",
    title = "FTLE Distribution on 2021, 2022, 2023 station datas"
  )

# PIGS
pig_vars <- c(
  "Chla",
  "Fuco",
  "But",
  "Per",
  "Hex",
  "Allo",
  "Chlb",
  "Zea", 
  "DvChla"
)

p1 <- ds %>%
  select(any_of(pig_vars)) %>%
  pivot_longer(
    everything(),
    names_to = "pigment",
    values_to = "value"
  ) %>%
  filter(is.finite(value)) %>%
  ggplot(aes(x = value)) +
  geom_histogram(bins = 50) +
  facet_wrap(~pigment, scales="free") +
  theme_bw()+
  scale_x_continuous(n.breaks = 3)

nb_missing <- sum(!is.finite(ds[[pig_vars[1]]]))
pct_missing <- nb_missing / nrow(ds) * 100

p2 <- ggplot() +
  annotate(
    "text",
    x = 0,
    y = 0,
    label = paste0("Missing values :\n", nb_missing, " (", round(pct_missing, 2), "%)", "\n \n Considered values :\n", nrow(ds)-nb_missing),
    size = 4
  ) +
  theme_void()

# affichage côte à côte
p1 + p2 +
  plot_layout(widths = c(3, 1)) +
  plot_annotation(
    title = "Distribution of pigment variables on 2021, 2022, 2023 station datas",
    subtitle = "Missing values are excluded from histograms"
  )


# nombre de jours où on a des données
idx <- which(!is.na(ds$Chla))
print(unique(as.Date(ds$time)[idx]))
# ----------- NASC
# NASC
ggplot(ds, aes(x = nasc, y = after_stat(count / sum(count) * 100))) +
  geom_histogram(bins = 200) +
  scale_x_log10() +
  theme_bw() +
  labs(
    x = "NASC (log10 scale)",
    y = "Percentage (%)",
    title = "NASC Distribution (200kHz) on 2021, 2022, 2023 station data"
  )
