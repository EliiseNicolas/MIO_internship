# Description 

# Distribution of pigmeann data on the OBSAUSTRAL Transect
# Quel pigment varie ? lesquels ont une distribution nulle ?
# Quels pigments covarient très fortements ?

# librairies

# variable globales 
path_pig <- "~/Documents/stage_MIO/pt_III/data_preprocessed/pigmeann/transect/pigs_transect_day_2018_2021_2023/pig_mask9_1d_2018_2021_2023.rds"
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

# Open pig data
pigs <- readRDS(path_pig)

# plot distribution all year
p1 <- pig_cond %>%
  select(any_of(pig_vars)) %>%
  pivot_longer(
    everything(),
    names_to = "pigment",
    values_to = "value"
  ) %>%
  filter(is.finite(value)) %>%
  ggplot(aes(x = value)) +
  geom_histogram(bins = 50) +
  facet_wrap(~pigment, scales = "free") +
  theme_bw() +
  scale_x_continuous(n.breaks = 3)


nb_missing <- sum(!is.finite(pig_cond[[pig_vars[1]]]))

pct_missing <- nb_missing / nrow(pig_cond) * 100

p2 <- ggplot() +
  annotate(
    "text",
    x = 0,
    y = 0,
    label = paste0(
      "Missing values :\n",
      nb_missing,
      " (",
      round(pct_missing, 2),
      "%)",
      "\n\nConsidered values :\n",
      nrow(pig_cond) - nb_missing
    ),
    size = 4
  ) +
  theme_void()


p1 + p2 +
  plot_layout(widths = c(3, 1)) +
  plot_annotation(
    title = paste0(
      "Distribution of pigment variables on 2018, 2021, 2023 transect day data"
    ),
    subtitle = "Missing values are excluded from histograms"
  )
  
# plot distrib per year 
years <- format(as.Date(pigs$time), "%Y")
years <- format(as.Date(pigs$time), "%Y")

for (year in c("2018", "2021", "2023")) {
  
  year_mask <- which(years == year)
  
  pig_year <- pig_cond[year_mask, ]
  
  p1 <- pig_year %>%
    select(any_of(pig_vars)) %>%
    pivot_longer(
      everything(),
      names_to = "pigment",
      values_to = "value"
    ) %>%
    filter(is.finite(value)) %>%
    ggplot(aes(x = value)) +
    geom_histogram(bins = 50) +
    facet_wrap(~pigment, scales = "free") +
    theme_bw() +
    scale_x_continuous(n.breaks = 3)
  
  nb_missing <- sum(!is.finite(pig_year[[pig_vars[1]]]))
  
  pct_missing <- nb_missing / nrow(pig_year) * 100
  
  p2 <- ggplot() +
    annotate(
      "text",
      x = 0,
      y = 0,
      label = paste0(
        "Missing values :\n",
        nb_missing,
        " (",
        round(pct_missing, 2),
        "%)",
        "\n\nConsidered values :\n",
        nrow(pig_year) - nb_missing
      ),
      size = 4
    ) +
    theme_void()
  
  print(
    p1 + p2 +
      plot_layout(widths = c(3, 1)) +
      plot_annotation(
        title = paste0(
          "Distribution of pigment variables - ", year
        ),
        subtitle = "Missing values are excluded from histograms"
      )
  )
}

# correlation des pigments 2018 2021 2023
# Variables pigmentaires
pig_data <- pig_cond %>%
  select(any_of(pig_vars))

# Matrice de corrélation
cor_all <- cor(
  pig_data,
  use = "pairwise.complete.obs",
  method = "pearson"
)

cor_all

library(reshape2)

cor_all_long <- melt(cor_all)

ggplot(cor_all_long, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  geom_text(aes(label = round(value, 2)), size = 3) +
  scale_fill_gradient2(
    limits = c(-1, 1),
    midpoint = 0
  ) +
  theme_bw() +
  labs(
    title = "Correlation matrix of pigment variables (2018-2021-2023 transect day data)",
    x = NULL,
    y = NULL,
    fill = "Correlation"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


library(patchwork)

years_to_plot <- c("2018", "2021", "2023")

plots <- list()

# Toutes les années
cor_all <- cor(
  pig_cond %>% select(any_of(pig_vars)),
  use = "pairwise.complete.obs",
  method = "pearson"
)

plots[["All years"]] <- ggplot(
  melt(cor_all),
  aes(Var1, Var2, fill = value)
) +
  geom_tile() +
  geom_text(aes(label = round(value, 2)), size = 2.5) +
  scale_fill_gradient2(
    limits = c(-1, 1),
    midpoint = 0
  ) +
  theme_bw() +
  labs(
    title = "All years",
    x = NULL,
    y = NULL,
    fill = "r"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


# Une matrice par année
for (year in years_to_plot) {
  
  year_mask <- years == year
  
  cor_year <- cor(
    pig_cond[year_mask, ] %>% select(any_of(pig_vars)),
    use = "pairwise.complete.obs",
    method = "pearson"
  )
  
  plots[[year]] <- ggplot(
    melt(cor_year),
    aes(Var1, Var2, fill = value)
  ) +
    geom_tile() +
    geom_text(aes(label = round(value, 2)), size = 2.5) +
    scale_fill_gradient2(
      limits = c(-1, 1),
      midpoint = 0
    ) +
    theme_bw() +
    labs(
      title = year,
      x = NULL,
      y = NULL,
      fill = "r"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

# Affichage
wrap_plots(plots, ncol = 2)