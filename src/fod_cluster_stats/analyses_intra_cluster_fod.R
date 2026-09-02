# description 

# From the final dataset NASC - FTLE - FOD - PIGMENTS
# We want to analyse the distribution of NASC, PIGMENTS and FTLE inside a FOD

# Libraries 
library(dplyr)
library(tidyr)
library(ggplot2)
library(multcompView)

# Global variables 
path <- "F:/data_elise/ds_NASC_pig_ftle_fod/ds_NASC_pig_ftle_fod_2018_2021_2023_transect_120kHz_day_mask9.rds"
ds <- readRDS(path)
fod_clusters <- sort(unique(ds$fod))
str(ds)
pigment_vars <- c(
  "Chla", "Per", "But", "Fuco", "Hex",
  "Allo", "Zea", "Chlb", "DvChla"
)

#################################################### Distributions simples des FTLE, pigments et NASC
# Keep only variables of interest
ds_long <- ds %>%
  select(fod, nasc, ftle, all_of(pigment_vars)) %>%
  pivot_longer(
    cols = all_of(pigment_vars),
    names_to = "pigment",
    values_to = "value"
  )

# Dsitribution générale des pigments
pigment_summary <- ds_long %>%
  group_by(fod, pigment) %>%
  summarise(
    n = sum(!is.na(value)),
    mean = mean(value, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    min = min(value, na.rm = TRUE),
    max = max(value, na.rm = TRUE),
    .groups = "drop"
  )

pigment_summary

ggplot(ds_long, aes(x = factor(fod), y = value)) +
  geom_violin(
    aes(group = factor(fod)),
    na.rm = TRUE,
    scale = "width"
  ) +
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA,
    na.rm = TRUE
  ) +
  facet_wrap(~ pigment, scales = "free_y") +
  labs(
    x = "FOD cluster",
    y = "Pigment concentration"
  ) +
  theme_bw()

# Dsitribution du NASC
ggplot(ds, aes(x = factor(fod), y = log(nasc))) +
  geom_violin(
    aes(group = factor(fod)),
    scale = "width",
    na.rm = TRUE
  ) +
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA,
    na.rm = TRUE
  ) +
  labs(
    x = "FOD cluster",
    y = "log(NASC)"
  ) +
  theme_bw()

# Distribution du NASC
ggplot(ds, aes(x = factor(fod), y = ftle)) +
  geom_violin(
    aes(group = factor(fod)),
    scale = "width",
    na.rm = TRUE
  ) +
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA,
    na.rm = TRUE
  ) +
  labs(
    x = "FOD cluster",
    y = "FTLE"
  ) +
  theme_bw()


# Distribution individuelle des pigments par cluster
pigments_long <- ds %>%
  select(fod, all_of(pigment_vars)) %>%
  pivot_longer(
    cols = all_of(pigment_vars),
    names_to = "pigment",
    values_to = "concentration"
  ) %>%
  filter(!is.na(concentration))

for (cl in sort(unique(pigments_long$fod))) {
  
  dat_cl <- pigments_long %>%
    filter(fod == cl)
  
  p <- ggplot(dat_cl, aes(x = pigment, y = concentration)) +
    geom_boxplot(
      na.rm = TRUE,
      outlier.shape = NA
    ) +
    geom_jitter(
      width = 0.15,
      alpha = 0.3,
      size = 1,
      na.rm = TRUE
    ) +
    labs(
      title = paste("Pigment concentrations - FOD cluster", cl),
      x = "Pigment",
      y = "Concentration"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  print(p)
}



# ============================================================
# NASC selon les clusters FOD
# ANOVA + Tukey + groupes statistiques + boxplot
# ============================================================

# -----------------------------
# Packages
# -----------------------------

library(dplyr)
library(tidyr)
library(ggplot2)
library(multcompView)

# -----------------------------
# 1. Préparation des données
# -----------------------------

nasc_fod <- ds %>%
  filter(
    !is.na(fod),
    !is.na(nasc),
    nasc >= 0
  ) %>%
  mutate(
    fod = factor(fod),
    log_nasc = log10(nasc + 1)
  )

# Vérification
str(nasc_fod)

# Nombre d'observations par cluster
table(nasc_fod$fod)


# -----------------------------
# 2. ANOVA sur NASC transformé
# -----------------------------

anova_model <- aov(
  log_nasc ~ fod,
  data = nasc_fod
)

summary(anova_model)


# -----------------------------
# 3. Test post-hoc de Tukey
# -----------------------------

tukey <- TukeyHSD(anova_model)

# Résultats Tukey pour fod
tukey_df <- as.data.frame(tukey$fod)

# Conserver les noms des comparaisons
tukey_df$comparison <- rownames(tukey_df)

# Séparer les deux groupes
tukey_df <- tukey_df %>%
  tidyr::separate(
    comparison,
    into = c("group1", "group2"),
    sep = "-"
  ) %>%
  arrange(`p adj`)

# Afficher toutes les comparaisons
tukey_df


# -----------------------------
# 4. Groupes statistiques
# -----------------------------

letters <- multcompLetters4(
  anova_model,
  tukey
)

letters_df <- data.frame(
  fod = names(letters$fod$Letters),
  letter = letters$fod$Letters
)

letters_df


# -----------------------------
# 5. Ordre des clusters
#    selon les groupes statistiques
# -----------------------------

letters_df <- letters_df %>%
  mutate(
    letter_order = match(letter, sort(unique(letter)))
  ) %>%
  arrange(letter_order, fod)

order_fod <- letters_df$fod


# Appliquer l'ordre sur l'axe X
nasc_fod <- nasc_fod %>%
  mutate(
    fod = factor(fod, levels = order_fod)
  )

letters_df <- letters_df %>%
  mutate(
    fod = factor(fod, levels = order_fod)
  )


# -----------------------------
# 6. Position des lettres
# -----------------------------

y_pos <- nasc_fod %>%
  group_by(fod) %>%
  summarise(
    y = max(nasc, na.rm = TRUE) * 1.3,
    .groups = "drop"
  )

letters_df <- letters_df %>%
  left_join(
    y_pos,
    by = "fod"
  )


# -----------------------------
# 7. Boxplot NASC
# -----------------------------

p <- ggplot(
  nasc_fod,
  aes(
    x = fod,
    y = nasc
  )
) +
  
  # Boxplots
  geom_boxplot(
    outlier.shape = NA,
    width = 0.65
  ) +
  
  # Points individuels
  geom_jitter(
    width = 0.15,
    alpha = 0.12,
    size = 0.8
  ) +
  
  # Lettres statistiques
  geom_text(
    data = letters_df,
    aes(
      x = fod,
      y = y,
      label = letter
    ),
    inherit.aes = FALSE,
    size = 5,
    fontface = "bold"
  ) +
  
  # Échelle logarithmique
  scale_y_log10() +
  
  # Labels
  labs(
    x = "FOD cluster",
    y = "NASC",
    title = "NASC distribution across FOD clusters"
  ) +
  
  # Thème
  theme_classic() +
  
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11),
    plot.title = element_text(
      size = 14,
      face = "bold"
    )
  )

p


############################################################# FTLE

# ============================================================
# FTLE selon les clusters FOD
# ANOVA + Tukey + groupes statistiques + boxplot
# ============================================================

# -----------------------------
# Packages
# -----------------------------

library(dplyr)
library(tidyr)
library(ggplot2)
library(multcompView)


# -----------------------------
# 1. Préparation des données
# -----------------------------

ftle_fod <- ds %>%
  filter(
    !is.na(fod),
    !is.na(ftle),
    is.finite(ftle)
  ) %>%
  mutate(
    fod = factor(fod)
  )

# Vérification
str(ftle_fod)

# Nombre d'observations par cluster
table(ftle_fod$fod)


# -----------------------------
# 2. ANOVA sur FTLE
# -----------------------------

anova_model <- aov(
  ftle ~ fod,
  data = ftle_fod
)

summary(anova_model)


# -----------------------------
# 3. Test post-hoc de Tukey
# -----------------------------

tukey <- TukeyHSD(anova_model)

# Résultats Tukey pour fod
tukey_df <- as.data.frame(tukey$fod)

# Conserver les noms des comparaisons
tukey_df$comparison <- rownames(tukey_df)

# Séparer les deux clusters
tukey_df <- tukey_df %>%
  tidyr::separate(
    comparison,
    into = c("group1", "group2"),
    sep = "-"
  ) %>%
  arrange(`p adj`)

# Afficher toutes les comparaisons
tukey_df


# -----------------------------
# 4. Groupes statistiques
# -----------------------------

letters <- multcompLetters4(
  anova_model,
  tukey
)

letters_df <- data.frame(
  fod = names(letters$fod$Letters),
  letter = letters$fod$Letters
)

letters_df


# -----------------------------
# 5. Ordre des clusters
#    selon les groupes statistiques
# -----------------------------

letters_df <- letters_df %>%
  mutate(
    letter_order = match(
      letter,
      sort(unique(letter))
    )
  ) %>%
  arrange(
    letter_order,
    fod
  )

order_fod <- letters_df$fod


# Appliquer l'ordre sur l'axe X
ftle_fod <- ftle_fod %>%
  mutate(
    fod = factor(
      fod,
      levels = order_fod
    )
  )

letters_df <- letters_df %>%
  mutate(
    fod = factor(
      fod,
      levels = order_fod
    )
  )


# -----------------------------
# 6. Position des lettres
# -----------------------------

y_pos <- ftle_fod %>%
  group_by(fod) %>%
  summarise(
    y = max(ftle, na.rm = TRUE) * 1.10,
    .groups = "drop"
  )

letters_df <- letters_df %>%
  left_join(
    y_pos,
    by = "fod"
  )


# -----------------------------
# 7. Boxplot FTLE
# -----------------------------

p <- ggplot(
  ftle_fod,
  aes(
    x = fod,
    y = ftle
  )
) +
  
  # Boxplots
  geom_boxplot(
    outlier.shape = NA,
    width = 0.65
  ) +
  
  # Points individuels
  geom_jitter(
    width = 0.15,
    alpha = 0.12,
    size = 0.8
  ) +
  
  # Lettres statistiques
  geom_text(
    data = letters_df,
    aes(
      x = fod,
      y = y,
      label = letter
    ),
    inherit.aes = FALSE,
    size = 5,
    fontface = "bold"
  ) +
  
  # Labels
  labs(
    x = "FOD cluster",
    y = "FTLE",
    title = "FTLE distribution across FOD clusters"
  ) +
  
  # Thème
  theme_classic() +
  
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11),
    plot.title = element_text(
      size = 14,
      face = "bold"
    )
  )

p

###################################################################### PIGMENTS

# ============================================================
# PIGMENTS selon les clusters FOD
# ANOVA + Tukey + groupes statistiques + boxplots
# ============================================================

# -----------------------------
# Packages
# -----------------------------

library(dplyr)
library(tidyr)
library(ggplot2)
library(multcompView)


# ============================================================
# 1. Liste des pigments
# ============================================================

pigment_vars <- c(
  "Chla",
  "Per",
  "But",
  "Fuco",
  "Hex",
  "Allo",
  "Zea",
  "Chlb",
  "DvChla"
)


# ============================================================
# 2. Fonction d'analyse pour un pigment
# ============================================================

analyse_pigment <- function(pigment_name, data) {
  
  # -----------------------------
  # Préparation des données
  # -----------------------------
  
  dat <- data %>%
    select(
      fod,
      value = all_of(pigment_name)
    ) %>%
    filter(
      !is.na(fod),
      !is.na(value),
      is.finite(value)
    ) %>%
    mutate(
      fod = factor(fod)
    )
  
  
  # -----------------------------
  # ANOVA
  # -----------------------------
  
  model <- aov(
    value ~ fod,
    data = dat
  )
  
  
  # -----------------------------
  # Tukey
  # -----------------------------
  
  tukey <- TukeyHSD(model)
  
  
  # -----------------------------
  # Lettres statistiques
  # -----------------------------
  
  letters <- multcompLetters4(
    model,
    tukey
  )
  
  letters_df <- data.frame(
    fod = names(letters$fod$Letters),
    letter = letters$fod$Letters
  )
  
  
  # -----------------------------
  # Ordre selon les groupes
  # statistiques
  # -----------------------------
  
  letters_df <- letters_df %>%
    mutate(
      letter_order = match(
        letter,
        sort(unique(letter))
      )
    ) %>%
    arrange(
      letter_order,
      fod
    )
  
  order_fod <- letters_df$fod
  
  
  # -----------------------------
  # Appliquer l'ordre
  # -----------------------------
  
  dat <- dat %>%
    mutate(
      fod = factor(
        fod,
        levels = order_fod
      )
    )
  
  letters_df <- letters_df %>%
    mutate(
      fod = factor(
        fod,
        levels = order_fod
      )
    )
  
  
  # -----------------------------
  # Position des lettres
  # -----------------------------
  
  y_pos <- dat %>%
    group_by(fod) %>%
    summarise(
      y = max(value, na.rm = TRUE) * 1.10,
      .groups = "drop"
    )
  
  letters_df <- letters_df %>%
    left_join(
      y_pos,
      by = "fod"
    )
  
  
  # -----------------------------
  # Graphique
  # -----------------------------
  
  p <- ggplot(
    dat,
    aes(
      x = fod,
      y = value
    )
  ) +
    
    # Boxplots
    geom_boxplot(
      outlier.shape = NA,
      width = 0.65
    ) +
    
    # Points individuels
    geom_jitter(
      width = 0.15,
      alpha = 0.12,
      size = 0.8
    ) +
    
    # Lettres statistiques
    geom_text(
      data = letters_df,
      aes(
        x = fod,
        y = y,
        label = letter
      ),
      inherit.aes = FALSE,
      size = 5,
      fontface = "bold"
    ) +
    
    labs(
      x = "FOD cluster",
      y = pigment_name,
      title = paste(
        pigment_name,
        "distribution across FOD clusters"
      )
    ) +
    
    theme_classic() +
    
    theme(
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 11),
      plot.title = element_text(
        size = 14,
        face = "bold"
      )
    )
  
  
  # -----------------------------
  # Retourner tous les résultats
  # -----------------------------
  
  list(
    data = dat,
    model = model,
    anova = summary(model),
    tukey = tukey,
    letters = letters_df,
    plot = p
  )
}


# ============================================================
# 3. Lancer l'analyse pour tous les pigments
# ============================================================

pigment_results <- lapply(
  pigment_vars,
  analyse_pigment,
  data = ds
)

names(pigment_results) <- pigment_vars


# ============================================================
# 4. Afficher les graphiques
# ============================================================

pigment_results$Chla$plot
pigment_results$Per$plot
pigment_results$But$plot
pigment_results$Fuco$plot
pigment_results$Hex$plot
pigment_results$Allo$plot
pigment_results$Zea$plot
pigment_results$Chlb$plot
pigment_results$DvChla$plot


##################################################################### PROPORTION DES PIGMENTS
# ============================================================
# PROPORTION DES PIGMENTS DANS LE POOL PIGMENTAIRE TOTAL
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)

# -----------------------------
# 1. Liste des pigments
# -----------------------------

pigment_vars <- c(
  "Chla",
  "Per",
  "But",
  "Fuco",
  "Hex",
  "Allo",
  "Zea",
  "Chlb",
  "DvChla"
)


# -----------------------------
# 2. Sélection des données
# -----------------------------

pigment_prop <- ds %>%
  select(
    fod,
    all_of(pigment_vars)
  ) %>%
  filter(
    !is.na(fod)
  )


# -----------------------------
# 3. Garder uniquement les
#    observations complètes
# -----------------------------

pigment_prop <- pigment_prop %>%
  filter(
    if_all(
      all_of(pigment_vars),
      ~ !is.na(.x) & is.finite(.x)
    )
  )


# -----------------------------
# 4. Calcul du pigment total
# -----------------------------

pigment_prop <- pigment_prop %>%
  mutate(
    pigment_total = rowSums(
      across(all_of(pigment_vars))
    )
  )


# -----------------------------
# 5. Passer en format long
# -----------------------------

pigment_prop_long <- pigment_prop %>%
  pivot_longer(
    cols = all_of(pigment_vars),
    names_to = "pigment",
    values_to = "value"
  ) %>%
  mutate(
    proportion = 100 * value / pigment_total
  )


# Vérification :
# la somme des proportions doit être 100 %
pigment_prop_long %>%
  group_by(fod) %>%
  summarise(
    mean_total = mean(pigment_total, na.rm = TRUE)
  )

pigment_prop_long %>%
  group_by(fod) %>%
  summarise(
    mean_proportion = mean(proportion, na.rm = TRUE)
  )

pigment_prop_long %>%
  group_by(fod) %>%
  summarise(
    min_sum = min(
      tapply(
        proportion,
        rep(
          1:(n() / length(pigment_vars)),
          each = length(pigment_vars)
        ),
        sum
      )
    ),
    max_sum = max(
      tapply(
        proportion,
        rep(
          1:(n() / length(pigment_vars)),
          each = length(pigment_vars)
        ),
        sum
      )
    )
  )

ggplot(
  pigment_prop_long,
  aes(
    x = factor(fod),
    y = proportion,
    fill = pigment
  )
) +
  stat_summary(
    fun = mean,
    geom = "bar",
    position = "fill"
  ) +
  scale_y_continuous(
    labels = scales::percent_format()
  ) +
  labs(
    x = "FOD cluster",
    y = "Relative pigment composition",
    fill = "Pigment",
    title = "Relative pigment composition across FOD clusters"
  ) +
  theme_classic()

################################################################### ANALYSE par Année et par cluster

# Vérifier les années disponibles
year <- format(ds$time, "%Y")
sort(unique(year))
ds$year <- year
# Nombre d'observations par année et cluster
ds %>%
  filter(!is.na(year), !is.na(fod)) %>%
  count(year, fod)

# NASC par année et FOD
nasc_year_fod <- ds %>%
  filter(
    !is.na(year),
    !is.na(fod),
    !is.na(nasc),
    is.finite(nasc),
    nasc >= 0
  ) %>%
  mutate(
    year = factor(year),
    fod = factor(fod),
    log_nasc = log10(nasc + 1)
  )

ggplot(
  nasc_year_fod,
  aes(
    x = year,
    y = nasc
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.65
  ) +
  geom_jitter(
    width = 0.15,
    alpha = 0.12,
    size = 0.7
  ) +
  facet_wrap(
    ~ fod,
    scales = "free_y"
  ) +
  scale_y_log10() +
  labs(
    x = "Year",
    y = "NASC",
    title = "NASC distribution within FOD clusters across years"
  ) +
  theme_classic()

# graphique des moyennes seulement 

nasc_summary <- nasc_year_fod %>%
  group_by(fod, year) %>%
  summarise(
    n = n(),
    mean = mean(nasc, na.rm = TRUE),
    median = median(nasc, na.rm = TRUE),
    sd = sd(nasc, na.rm = TRUE),
    se = sd / sqrt(n),
    .groups = "drop"
  )

nasc_summary

ggplot(
  nasc_summary,
  aes(
    x = year,
    y = mean,
    group = 1
  )
) +
  geom_line() +
  geom_point(size = 3) +
  facet_wrap(
    ~ fod,
    scales = "free_y"
  ) +
  scale_y_log10() +
  labs(
    x = "Year",
    y = "Mean NASC",
    title = "Temporal evolution of mean NASC within FOD clusters"
  ) +
  theme_classic()

# médiane 
ggplot(
  nasc_summary,
  aes(
    x = year,
    y = median,
    group = 1
  )
) +
  geom_line() +
  geom_point(size = 3) +
  facet_wrap(
    ~ fod,
    scales = "free_y"
  ) +
  scale_y_log10() +
  labs(
    x = "Year",
    y = "Median NASC",
    title = "Temporal evolution of median NASC within FOD clusters"
  ) +
  theme_classic()

# statistiques
# ============================================================
# NASC selon ANNEE et FOD
# ANOVA + Tukey + groupes statistiques par année dans chaque FOD
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(emmeans)
library(multcompView)
library(multcomp)

ds$year <- format(ds$time, "%Y")

# ============================================================
# 1. Préparation des données
# ============================================================

nasc_year_fod <- ds %>%
  filter(
    !is.na(year),
    !is.na(fod),
    !is.na(nasc),
    is.finite(nasc),
    nasc >= 0
  ) %>%
  mutate(
    year = factor(year),
    fod = factor(fod),
    log_nasc = log10(nasc + 1)
  )


# Vérification
table(
  nasc_year_fod$year,
  nasc_year_fod$fod
)


# ============================================================
# 2. ANOVA avec interaction année × FOD
# ============================================================

nasc_model <- aov(
  log_nasc ~ year * fod,
  data = nasc_year_fod
)

summary(nasc_model)


# ============================================================
# 3. EMMEANS
#    Comparaison des années à l'intérieur de chaque FOD
# ============================================================

nasc_emm <- emmeans(
  nasc_model,
  ~ year | fod
)


# ============================================================
# 4. Comparaisons de Tukey
# ============================================================

nasc_pairs <- pairs(
  nasc_emm,
  adjust = "tukey"
)

nasc_pairs_df <- as.data.frame(nasc_pairs)

nasc_pairs_df

# ============================================================
# 5. Groupes statistiques par année dans chaque FOD
# ============================================================

letters_list <- cld(
  nasc_emm,
  adjust = "tukey",
  Letters = letters,
  sort = FALSE
)
letters_list
letters_df <- as.data.frame(letters_list) %>%
  dplyr::select(
    fod,
    year,
    .group
  ) %>%
  dplyr::mutate(
    letter = gsub(" ", "", .group)
  ) %>%
  dplyr::select(
    fod,
    year,
    letter
  )

letters_df

# ============================================================
# 6. Position des lettres
# ============================================================

y_pos <- nasc_year_fod %>%
  group_by(fod) %>%
  summarise(
    y = max(nasc, na.rm = TRUE) * 1.20,
    .groups = "drop"
  )

letters_df <- letters_df %>%
  left_join(
    y_pos,
    by = "fod"
  )

# ============================================================
# 7. Boxplots NASC année × FOD
# ============================================================

p <- ggplot(
  nasc_year_fod,
  aes(
    x = year,
    y = nasc
  )
) +
  
  # Boxplots
  geom_boxplot(
    outlier.shape = NA,
    width = 0.65
  ) +
  
  # Points individuels
  geom_jitter(
    width = 0.15,
    alpha = 0.12,
    size = 0.7
  ) +
  
  # Lettres statistiques
  geom_text(
    data = letters_df,
    aes(
      x = year,
      y = y,
      label = letter
    ),
    inherit.aes = FALSE,
    size = 3,
    fontface = "bold"
  ) +
  
  # Un panneau par FOD
  facet_wrap(
    ~ fod,
    scales = "free_y"
  ) +
  
  # Échelle logarithmique
  scale_y_log10() +
  
  # Labels
  labs(
    x = "Year",
    y = "NASC",
    title = "NASC distribution within FOD clusters across years"
  ) +
  
  theme_classic() +
  
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11),
    strip.text = element_text(
      size = 12,
      face = "bold"
    ),
    plot.title = element_text(
      size = 14,
      face = "bold"
    )
  )

p

########################################################## Pigment et FTLE
# ============================================================
# PACKAGES
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(emmeans)
library(multcompView)
library(multcomp)


# ============================================================
# ANNEE
# ============================================================

ds$year <- format(ds$time, "%Y")

# ============================================================
# FONCTION D'ANALYSE
# FTLE ou PIGMENT selon ANNEE et FOD
# ============================================================

analyse_year_fod <- function(
    variable,
    data,
    y_label = variable
) {
  
  # ----------------------------------------------------------
  # 1. Préparation des données
  # ----------------------------------------------------------
  
  dat <- data %>%
    dplyr::select(
      year,
      fod,
      value = all_of(variable)
    ) %>%
    dplyr::filter(
      !is.na(year),
      !is.na(fod),
      !is.na(value),
      is.finite(value)
    ) %>%
    dplyr::mutate(
      year = factor(year),
      fod = factor(fod)
    )
  
  
  # ----------------------------------------------------------
  # 2. Vérification des effectifs
  # ----------------------------------------------------------
  
  cat("\n\n============================================\n")
  cat("VARIABLE :", variable, "\n")
  cat("============================================\n")
  
  print(
    table(
      dat$year,
      dat$fod
    )
  )
  
  
  # ----------------------------------------------------------
  # 3. ANOVA avec interaction année × FOD
  # ----------------------------------------------------------
  
  model <- aov(
    value ~ year * fod,
    data = dat
  )
  
  cat("\n--- ANOVA ---\n")
  print(summary(model))
  
  
  # ----------------------------------------------------------
  # 4. EMMEANS
  # Comparaison des années à l'intérieur de chaque FOD
  # ----------------------------------------------------------
  
  emm <- emmeans(
    model,
    ~ year | fod
  )
  
  
  # ----------------------------------------------------------
  # 5. Comparaisons de Tukey
  # ----------------------------------------------------------
  
  pairs_df <- pairs(
    emm,
    adjust = "tukey"
  ) %>%
    as.data.frame()
  
  
  cat("\n--- Comparaisons entre années ---\n")
  print(pairs_df)
  
  
  # ----------------------------------------------------------
  # 6. Groupes statistiques
  # ----------------------------------------------------------
  
  letters_list <- cld(
    emm,
    adjust = "tukey",
    Letters = letters,
    sort = FALSE
  )
  
  
  cat("\n--- Groupes statistiques ---\n")
  print(letters_list)
  
  
  # ----------------------------------------------------------
  # 7. Préparation des lettres
  # ----------------------------------------------------------
  
  letters_df <- as.data.frame(letters_list) %>%
    dplyr::select(
      fod,
      year,
      .group
    ) %>%
    dplyr::mutate(
      letter = gsub(" ", "", .group)
    ) %>%
    dplyr::select(
      fod,
      year,
      letter
    ) %>%
    dplyr::filter(
      !is.na(letter)
    )
  
  
  # ----------------------------------------------------------
  # 8. Position des lettres
  # ----------------------------------------------------------
  
  y_pos <- dat %>%
    dplyr::group_by(fod) %>%
    dplyr::summarise(
      y = max(value, na.rm = TRUE) * 1.20,
      .groups = "drop"
    )
  
  
  letters_df <- letters_df %>%
    dplyr::left_join(
      y_pos,
      by = "fod"
    )
  
  
  # ----------------------------------------------------------
  # 9. Graphique
  # ----------------------------------------------------------
  
  p <- ggplot(
    dat,
    aes(
      x = year,
      y = value
    )
  ) +
    
    # Boxplots
    geom_boxplot(
      outlier.shape = NA,
      width = 0.65
    ) +
    
    # Points individuels
    geom_jitter(
      width = 0.15,
      alpha = 0.12,
      size = 0.7
    ) +
    
    # Lettres statistiques
    geom_text(
      data = letters_df,
      aes(
        x = year,
        y = y,
        label = letter
      ),
      inherit.aes = FALSE,
      size = 3,
      fontface = "bold"
    ) +
    
    # Un panneau par FOD
    facet_wrap(
      ~ fod,
      scales = "free_y"
    ) +
    
    # Labels
    labs(
      x = "Year",
      y = y_label,
      title = paste(
        y_label,
        "distribution within FOD clusters across years"
      )
    ) +
    
    theme_classic() +
    
    theme(
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 11),
      strip.text = element_text(
        size = 12,
        face = "bold"
      ),
      plot.title = element_text(
        size = 14,
        face = "bold"
      )
    )
  
  
  # ----------------------------------------------------------
  # 10. Retourner tous les résultats
  # ----------------------------------------------------------
  
  return(
    list(
      data = dat,
      model = model,
      anova = summary(model),
      emm = emm,
      pairs = pairs_df,
      letters = letters_df,
      plot = p
    )
  )
}

# ============================================================
# FTLE
# ============================================================

ftle_year_fod <- analyse_year_fod(
  variable = "ftle",
  data = ds,
  y_label = "FTLE"
)


# ============================================================
# ANALYSE DES PIGMENTS SELON ANNEE ET FOD
# ANOVA + EMMEANS + TUKEY + GROUPES STATISTIQUES + BOXPLOTS
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(emmeans)
library(multcomp)
library(multcompView)
library(purrr)

# ============================================================
# 1. Année
# ============================================================

ds$year <- format(ds$time, "%Y")


# ============================================================
# 2. Liste des pigments
# ============================================================

pigments <- c(
  "Chla",
  "Chlb",
  "DvChla",
  "Fuco",
  "But",
  "Hex",
  "Allo",
  "Zea"
)


# ============================================================
# 3. Fonction d'analyse
# ============================================================

analyse_pigment <- function(pigment, data = ds) {
  
  cat("\n")
  cat("============================================\n")
  cat("VARIABLE :", pigment, "\n")
  cat("============================================\n")
  
  
  # ----------------------------------------------------------
  # Préparation des données
  # ----------------------------------------------------------
  
  dat <- data %>%
    filter(
      !is.na(year),
      !is.na(fod),
      !is.na(.data[[pigment]]),
      is.finite(.data[[pigment]]),
      .data[[pigment]] >= 0
    ) %>%
    mutate(
      year = factor(year),
      fod = factor(fod),
      pigment_value = .data[[pigment]]
    )
  
  
  # ----------------------------------------------------------
  # Vérification des effectifs
  # ----------------------------------------------------------
  
  cat("\n--- Effectifs ---\n")
  
  print(
    table(
      dat$year,
      dat$fod
    )
  )
  
  
  # ----------------------------------------------------------
  # ANOVA
  # ----------------------------------------------------------
  
  model <- aov(
    pigment_value ~ year * fod,
    data = dat
  )
  
  cat("\n--- ANOVA ---\n")
  print(summary(model))
  
  
  # ----------------------------------------------------------
  # EMMEANS
  # Comparaison des années dans chaque FOD
  # ----------------------------------------------------------
  
  emm <- emmeans(
    model,
    ~ year | fod
  )
  
  
  # ----------------------------------------------------------
  # Comparaisons de Tukey
  # ----------------------------------------------------------
  
  pairs_res <- pairs(
    emm,
    adjust = "tukey"
  )
  
  pairs_df <- as.data.frame(pairs_res)
  
  cat("\n--- Comparaisons entre années ---\n")
  print(pairs_df)
  
  
  # ----------------------------------------------------------
  # Groupes statistiques
  # ----------------------------------------------------------
  
  letters_list <- cld(
    emm,
    adjust = "tukey",
    Letters = letters,
    sort = FALSE
  )
  
  letters_raw <- as.data.frame(letters_list)
  
  cat("\n--- Groupes statistiques ---\n")
  print(letters_raw)
  
  
  # ----------------------------------------------------------
  # Nettoyage des lettres
  # ----------------------------------------------------------
  
  letters_df <- letters_raw %>%
    dplyr::filter(
      !is.na(emmean),
      !is.na(.group)
    ) %>%
    dplyr::select(
      fod,
      year,
      emmean,
      .group
    ) %>%
    dplyr::mutate(
      letter = gsub(" ", "", .group)
    ) %>%
    dplyr::select(
      fod,
      year,
      emmean,
      letter
    )
  
  
  # ----------------------------------------------------------
  # Position des lettres
  # ----------------------------------------------------------
  
  y_pos <- dat %>%
    group_by(fod) %>%
    summarise(
      y = max(
        pigment_value,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    mutate(
      y = ifelse(
        y > 0,
        y * 1.20,
        0.001
      )
    )
  
  
  letters_df <- letters_df %>%
    left_join(
      y_pos,
      by = "fod"
    )
  
  
  # ----------------------------------------------------------
  # Graphique
  # ----------------------------------------------------------
  
  p <- ggplot(
    dat,
    aes(
      x = year,
      y = pigment_value
    )
  ) +
    
    geom_boxplot(
      outlier.shape = NA,
      width = 0.65
    ) +
    
    geom_jitter(
      width = 0.15,
      alpha = 0.12,
      size = 0.7
    ) +
    
    geom_text(
      data = letters_df,
      aes(
        x = year,
        y = y,
        label = letter
      ),
      inherit.aes = FALSE,
      size = 3,
      fontface = "bold"
    ) +
    
    facet_wrap(
      ~ fod,
      scales = "free_y"
    ) +
    
    scale_y_log10() +
    
    labs(
      x = "Year",
      y = pigment,
      title = paste(
        pigment,
        "distribution within FOD clusters across years"
      )
    ) +
    
    theme_classic() +
    
    theme(
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 11),
      strip.text = element_text(
        size = 12,
        face = "bold"
      ),
      plot.title = element_text(
        size = 14,
        face = "bold"
      )
    )
  
  
  # ----------------------------------------------------------
  # Retour de la fonction
  # ----------------------------------------------------------
  
  return(
    list(
      data = dat,
      model = model,
      emm = emm,
      pairs = pairs_df,
      letters = letters_df,
      plot = p
    )
  )
}


# ============================================================
# 4. Lancer l'analyse pour tous les pigments
# ============================================================

pigment_year_fod <- setNames(
  lapply(
    pigments,
    analyse_pigment
  ),
  pigments
)


# ============================================================
# 5. Afficher les graphiques
# ============================================================

pigment_year_fod$Chla$plot
pigment_year_fod$Chlb$plot
pigment_year_fod$Fuco$plot
pigment_year_fod$But$plot
pigment_year_fod$Hex$plot
pigment_year_fod$Allo$plot
pigment_year_fod$Zea$plot
pigment_year_fod$DvChla$plot
