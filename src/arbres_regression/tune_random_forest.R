# ============================================================
# RANDOM FOREST TUNING
# + SPATIO-TEMPORAL LEAKAGE REMOVAL
# + LOYO
# + RANDOM 80/20
# + FOD_0 TRUE/FALSE
# + LOG TRUE/FALSE
# ============================================================

library(randomForest)
library(dplyr)
library(ggplot2)

set.seed(123)


# ============================================================
# 01 - PARAMETRES
# ============================================================

freq <- 18
diurnal_period <- 3
dp <- "day"

pigment_type <- "chla_ratio"

n_repeats <- 10
train_fraction <- 0.8

# Leakage
time_thr <- 3600

distance_thresholds <- c(
  0.25,
  0.5,
  0.75,
  1
)

# FOD / log
FOD_values <- c(FALSE, TRUE)
log_values <- c(FALSE, TRUE)

# Random Forest
ntree_values <- c(
  100,
  300,
  500
)

nodesize_values <- c(
  1,
  3,
  5,
  10
)

# Proportion du nombre de prédicteurs
# permettant de construire mtry
mtry_fraction <- c(
  0.5,
  1,
  1.5
)


# ============================================================
# 02 - PREPARATION DES DONNEES
# ============================================================

prepare_data <- function(datas,
                         pigment_type = "chla_ratio",
                         FOD_0 = FALSE,
                         log = TRUE,
                         diurnal_period = 3) {
  
  dat <- datas
  
  
  # ----------------------------------------------------------
  # Day / night
  # ----------------------------------------------------------
  
  dat <- dat[
    dat$day == diurnal_period,
  ]
  
  
  # ----------------------------------------------------------
  # NASC extreme values
  # ----------------------------------------------------------
  
  q <- quantile(
    dat$nasc,
    probs = c(0.05, 0.95),
    na.rm = TRUE
  )
  
  dat <- dat |>
    filter(
      nasc >= q[1],
      nasc <= q[2]
    )
  
  
  # ----------------------------------------------------------
  # Log NASC
  # ----------------------------------------------------------
  
  if (log) {
    dat$nasc <- log10(dat$nasc)
  }
  
  
  # ----------------------------------------------------------
  # FOD_0
  # ----------------------------------------------------------
  
  if (FOD_0) {
    
    fod_0 <- as.character(dat$fod)
    
    fod_num <- suppressWarnings(
      as.numeric(fod_0)
    )
    
    fod_0[!is.na(fod_num) & fod_num > 6] <- "0"
    
    dat$fod <- fod_0
  }
  
  
  # ----------------------------------------------------------
  # Base dataframe
  # ----------------------------------------------------------
  
  df <- data.frame(
    NASC = dat$nasc,
    year = format(dat$time_nasc, "%Y"),
    time = dat$time_nasc,
    fod = dat$fod,
    ftle = dat$ftle
  )
  
  
  # ----------------------------------------------------------
  # Pigments
  # ----------------------------------------------------------
  
  if (pigment_type == "chla_ratio") {
    
    vars_num <- c(
      "NASC",
      "per_ratio_chla",
      "but_ratio_chla",
      "fuco_ratio_chla",
      "hex_ratio_chla",
      "allo_ratio_chla",
      "zea_ratio_chla",
      "chlb_ratio_chla",
      "total_chla",
      "ftle"
    )
    
    df$total_chla <- dat$Chla
    df$per_ratio_chla <- dat$Per_Chla
    df$but_ratio_chla <- dat$But_Chla
    df$fuco_ratio_chla <- dat$Fuco_Chla
    df$hex_ratio_chla <- dat$Hex_Chla
    df$allo_ratio_chla <- dat$Allo_Chla
    df$zea_ratio_chla <- dat$Zea_Chla
    df$chlb_ratio_chla <- dat$Chlb_Chla
  }
  
  
  if (pigment_type == "total_ratio") {
    
    vars_num <- c(
      "NASC",
      "chla_ratio_total",
      "per_ratio_total",
      "but_ratio_total",
      "fuco_ratio_total",
      "hex_ratio_total",
      "allo_ratio_total",
      "zea_ratio_total",
      "chlb_ratio_total",
      "total_pig",
      "ftle"
    )
    
    df$chla_ratio_total <- dat$Chla_total
    df$per_ratio_total <- dat$Per_total
    df$but_ratio_total <- dat$But_total
    df$fuco_ratio_total <- dat$Fuco_total
    df$hex_ratio_total <- dat$Hex_total
    df$allo_ratio_total <- dat$Allo_total
    df$zea_ratio_total <- dat$Zea_total
    df$chlb_ratio_total <- dat$Chlb_total
    df$total_pig <- dat$total_pig
  }
  
  
  if (pigment_type == "conc") {
    
    vars_num <- c(
      "NASC",
      "chla",
      "per",
      "but",
      "fuco",
      "hex",
      "allo",
      "zea",
      "chlb",
      "total_pig",
      "ftle"
    )
    
    df$chla <- dat$Chla
    df$per <- dat$Per
    df$but <- dat$But
    df$fuco <- dat$Fuco
    df$hex <- dat$Hex
    df$allo <- dat$Allo
    df$zea <- dat$Zea
    df$chlb <- dat$Chlb
    df$total_pig <- dat$total_pig
  }
  
  
  # ----------------------------------------------------------
  # FOD factor
  # ----------------------------------------------------------
  
  df$fod <- as.factor(df$fod)
  
  df$fod[df$fod == "NA"] <- NA
  
  df$fod <- droplevels(df$fod)
  
  
  # ----------------------------------------------------------
  # Missing values
  # ----------------------------------------------------------
  
  df <- df |>
    filter(
      if_all(
        all_of(setdiff(vars_num, "NASC")),
        ~ !is.na(.)
      ),
      !is.na(fod)
    )
  
  
  # ----------------------------------------------------------
  # Retour
  # ----------------------------------------------------------
  
  list(
    df = df,
    vars_num = vars_num
  )
}


# ============================================================
# 03 - SUPPRESSION DU DATA LEAKAGE
# ============================================================

remove_close_pairs <- function(
    df_input,
    predictor_vars,
    dist_thr = 0.5,
    time_thr = 3600
) {
  
  df_work <- df_input
  
  n_initial <- nrow(df_work)
  
  removed_rows <- integer(0)
  
  n_iterations <- 0
  
  
  repeat {
    
    if (nrow(df_work) < 2) {
      break
    }
    
    
    # --------------------------------------------------------
    # Distance calculée UNIQUEMENT sur les prédicteurs
    # --------------------------------------------------------
    
    df_scaled <- scale(
      df_work[, predictor_vars, drop = FALSE]
    )
    
    dist_mat <- as.matrix(
      dist(df_scaled)
    )
    
    
    # --------------------------------------------------------
    # Différence temporelle
    # --------------------------------------------------------
    
    time_diff <- abs(
      outer(
        as.numeric(df_work$time),
        as.numeric(df_work$time),
        FUN = "-"
      )
    )
    
    
    # --------------------------------------------------------
    # Paires proches
    # --------------------------------------------------------
    
    pairs <- which(
      time_diff <= time_thr &
        dist_mat < dist_thr,
      arr.ind = TRUE
    )
    
    
    pairs <- pairs[
      pairs[, 1] < pairs[, 2],
      ,
      drop = FALSE
    ]
    
    
    n_pairs <- nrow(pairs)
    
    
    if (n_pairs == 0) {
      break
    }
    
    
    # --------------------------------------------------------
    # Suppression du deuxième élément de la première paire
    # --------------------------------------------------------
    
    pair <- pairs[1, ]
    
    remove_row <- pair[2]
    
    removed_rows <- c(
      removed_rows,
      remove_row
    )
    
    n_iterations <- n_iterations + 1
    
    df_work <- df_work[
      -remove_row,
      ,
      drop = FALSE
    ]
  }
  
  
  list(
    df = df_work,
    n_initial = n_initial,
    n_remaining = nrow(df_work),
    n_removed = n_initial - nrow(df_work),
    percent_removed =
      100 * (n_initial - nrow(df_work)) / n_initial,
    n_iterations = n_iterations,
    removed_rows = removed_rows
  )
}


# ============================================================
# 04 - METRIQUES
# ============================================================

calculate_metrics <- function(obs, pred) {
  
  RMSE <- sqrt(
    mean(
      (pred - obs)^2
    )
  )
  
  R2 <- 1 -
    sum(
      (obs - pred)^2
    ) /
    sum(
      (obs - mean(obs))^2
    )
  
  MAE <- mean(
    abs(pred - obs)
  )
  
  data.frame(
    RMSE = RMSE,
    R2 = R2,
    MAE = MAE
  )
}


# ============================================================
# 05 - RANDOM FOREST
# ============================================================

fit_rf <- function(
    train,
    test,
    vars_num,
    mtry,
    nodesize,
    ntree
) {
  
  
  # ----------------------------------------------------------
  # Scaling
  # ----------------------------------------------------------
  
  train_scaled <- scale(
    train[, vars_num]
  )
  
  center <- attr(
    train_scaled,
    "scaled:center"
  )
  
  scale_values <- attr(
    train_scaled,
    "scaled:scale"
  )
  
  test_scaled <- scale(
    test[, vars_num],
    center = center,
    scale = scale_values
  )
  
  
  ds_train <- train
  ds_test <- test
  
  ds_train[, vars_num] <- train_scaled
  ds_test[, vars_num] <- test_scaled
  
  
  # ----------------------------------------------------------
  # RF
  # ----------------------------------------------------------
  
  model <- randomForest(
    NASC ~ .,
    data = ds_train,
    ntree = ntree,
    mtry = mtry,
    nodesize = nodesize,
    importance = FALSE
  )
  
  
  # ----------------------------------------------------------
  # Prediction
  # ----------------------------------------------------------
  
  prediction <- predict(
    model,
    newdata = ds_test
  )
  
  observed <- ds_test$NASC
  
  
  calculate_metrics(
    observed,
    prediction
  )
}


# ============================================================
# 06 - RESULTATS
# ============================================================

results <- data.frame()

combination_id <- 0


# ============================================================
# 07 - BOUCLE FOD / LOG / LEAKAGE
# ============================================================

for (FOD_0 in FOD_values) {
  
  for (log in log_values) {
    
    
    cat("\n\n")
    cat("####################################################\n")
    cat("FOD_0 =", FOD_0,
        " | LOG =", log, "\n")
    cat("####################################################\n")
    
    
    # --------------------------------------------------------
    # Préparation
    # --------------------------------------------------------
    
    prepared <- prepare_data(
      datas = datas,
      pigment_type = pigment_type,
      FOD_0 = FOD_0,
      log = log,
      diurnal_period = diurnal_period
    )
    
    df_raw <- prepared$df
    vars_num <- prepared$vars_num
    
    predictor_vars <- setdiff(
      vars_num,
      "NASC"
    )
    
    
    cat(
      "Observations avant leakage :",
      nrow(df_raw),
      "\n"
    )
    
    
    # ========================================================
    # DIFFERENTS SEUILS DE LEAKAGE
    # ========================================================
    
    for (dist_thr in distance_thresholds) {
      
      
      cat("\n")
      cat("-----------------------------------------------\n")
      cat("Distance threshold =", dist_thr, "\n")
      cat("-----------------------------------------------\n")
      
      
      # ------------------------------------------------------
      # Suppression leakage
      # ------------------------------------------------------
      
      leakage <- remove_close_pairs(
        df_input = df_raw,
        predictor_vars = predictor_vars,
        dist_thr = dist_thr,
        time_thr = time_thr
      )
      
      df <- leakage$df
      
      
      cat(
        "Observations supprimées :",
        leakage$n_removed,
        "\n"
      )
      
      cat(
        "Observations restantes :",
        leakage$n_remaining,
        "\n"
      )
      
      cat(
        "Pourcentage supprimé :",
        round(
          leakage$percent_removed,
          1
        ),
        "%\n"
      )
      
      
      # ------------------------------------------------------
      # Nombre de prédicteurs
      # ------------------------------------------------------
      
      p <- length(
        predictor_vars
      )
      
      
      # ------------------------------------------------------
      # mtry
      # ------------------------------------------------------
      
      mtry_values <- unique(
        pmin(
          p,
          pmax(
            1,
            round(
              p * mtry_fraction
            )
          )
        )
      )
      
      
      # ------------------------------------------------------
      # Grid
      # ------------------------------------------------------
      
      grid <- expand.grid(
        mtry = mtry_values,
        nodesize = nodesize_values,
        ntree = ntree_values
      )
      
      
      cat(
        "Nombre de combinaisons :",
        nrow(grid),
        "\n"
      )
      
      
      # ======================================================
      # LOYO
      # ======================================================
      
      for (i in seq_len(nrow(grid))) {
        
        hp <- grid[i, ]
        
        years <- sort(
          unique(df$year)
        )
        
        loyo_metrics <- data.frame()
        
        
        for (yr in years) {
          
          train <- df[
            df$year != yr,
          ]
          
          test <- df[
            df$year == yr,
          ]
          
          
          if (
            nrow(train) < 10 ||
            nrow(test) < 2
          ) {
            next
          }
          
          
          metrics <- fit_rf(
            train = train,
            test = test,
            vars_num = vars_num,
            mtry = hp$mtry,
            nodesize = hp$nodesize,
            ntree = hp$ntree
          )
          
          
          metrics$year <- yr
          
          loyo_metrics <- rbind(
            loyo_metrics,
            metrics
          )
        }
        
        
        if (nrow(loyo_metrics) > 0) {
          
          combination_id <- combination_id + 1
          
          results <- rbind(
            results,
            data.frame(
              combination = combination_id,
              validation = "LOYO",
              FOD_0 = FOD_0,
              log = log,
              dist_thr = dist_thr,
              n_initial = leakage$n_initial,
              n_removed = leakage$n_removed,
              percent_removed = leakage$percent_removed,
              n_remaining = leakage$n_remaining,
              mtry = hp$mtry,
              nodesize = hp$nodesize,
              ntree = hp$ntree,
              RMSE = mean(loyo_metrics$RMSE),
              RMSE_sd = sd(loyo_metrics$RMSE),
              R2 = mean(loyo_metrics$R2),
              R2_sd = sd(loyo_metrics$R2),
              MAE = mean(loyo_metrics$MAE),
              MAE_sd = sd(loyo_metrics$MAE)
            )
          )
        }
      }
      
      
      # ======================================================
      # RANDOM SPLIT 80/20
      # ======================================================
      
      for (i in seq_len(nrow(grid))) {
        
        hp <- grid[i, ]
        
        rs_metrics <- data.frame()
        
        
        for (rep in seq_len(n_repeats)) {
          
          set.seed(
            1000 + rep
          )
          
          
          n_train <- floor(
            train_fraction *
              nrow(df)
          )
          
          train_index <- sample(
            seq_len(nrow(df)),
            size = n_train
          )
          
          test_index <- setdiff(
            seq_len(nrow(df)),
            train_index
          )
          
          
          train <- df[
            train_index,
          ]
          
          test <- df[
            test_index,
          ]
          
          
          metrics <- fit_rf(
            train = train,
            test = test,
            vars_num = vars_num,
            mtry = hp$mtry,
            nodesize = hp$nodesize,
            ntree = hp$ntree
          )
          
          
          rs_metrics <- rbind(
            rs_metrics,
            metrics
          )
        }
        
        
        combination_id <- combination_id + 1
        
        results <- rbind(
          results,
          data.frame(
            combination = combination_id,
            validation = "RS",
            FOD_0 = FOD_0,
            log = log,
            dist_thr = dist_thr,
            n_initial = leakage$n_initial,
            n_removed = leakage$n_removed,
            percent_removed = leakage$percent_removed,
            n_remaining = leakage$n_remaining,
            mtry = hp$mtry,
            nodesize = hp$nodesize,
            ntree = hp$ntree,
            RMSE = mean(rs_metrics$RMSE),
            RMSE_sd = sd(rs_metrics$RMSE),
            R2 = mean(rs_metrics$R2),
            R2_sd = sd(rs_metrics$R2),
            MAE = mean(rs_metrics$MAE),
            MAE_sd = sd(rs_metrics$MAE)
          )
        )
      }
    }
  }
}


# ============================================================
# 08 - MEILLEUR MODELE LOYO
# ============================================================

best_LOYO <- results |>
  filter(
    validation == "LOYO"
  ) |>
  arrange(
    RMSE
  ) |>
  slice(1)

cat("\n\n")
cat("====================================================\n")
cat("MEILLEUR MODELE LOYO\n")
cat("====================================================\n")

print(best_LOYO)


# ============================================================
# 09 - MEILLEUR MODELE RANDOM SPLIT
# ============================================================

best_RS <- results |>
  filter(
    validation == "RS"
  ) |>
  arrange(
    RMSE
  ) |>
  slice(1)

cat("\n\n")
cat("====================================================\n")
cat("MEILLEUR MODELE RANDOM SPLIT\n")
cat("====================================================\n")

print(best_RS)


# ============================================================
# 10 - TOP 10 LOYO
# ============================================================

cat("\n\n")
cat("====================================================\n")
cat("TOP 10 LOYO\n")
cat("====================================================\n")

results |>
  filter(
    validation == "LOYO"
  ) |>
  arrange(
    RMSE
  ) |>
  select(
    validation,
    FOD_0,
    log,
    dist_thr,
    n_removed,
    n_remaining,
    mtry,
    nodesize,
    ntree,
    RMSE,
    RMSE_sd,
    R2,
    MAE
  ) |>
  slice_head(
    n = 10
  ) |>
  print()


# ============================================================
# 11 - TOP 10 RANDOM SPLIT
# ============================================================

cat("\n\n")
cat("====================================================\n")
cat("TOP 10 RANDOM SPLIT\n")
cat("====================================================\n")

results |>
  filter(
    validation == "RS"
  ) |>
  arrange(
    RMSE
  ) |>
  select(
    validation,
    FOD_0,
    log,
    dist_thr,
    n_removed,
    n_remaining,
    mtry,
    nodesize,
    ntree,
    RMSE,
    RMSE_sd,
    R2,
    MAE
  ) |>
  slice_head(
    n = 10
  ) |>
  print()


# ============================================================
# 12 - COMPARAISON LOG / FOD / LEAKAGE
# ============================================================

summary_conditions <- results |>
  group_by(
    validation,
    FOD_0,
    log,
    dist_thr
  ) |>
  summarise(
    RMSE_mean = mean(RMSE),
    RMSE_sd = mean(RMSE_sd),
    R2_mean = mean(R2),
    MAE_mean = mean(MAE),
    n_removed = first(n_removed),
    n_remaining = first(n_remaining),
    .groups = "drop"
  ) |>
  arrange(
    validation,
    RMSE_mean
  )

print(summary_conditions)


# ============================================================
# 13 - PLOT RMSE EN FONCTION DU LEAKAGE
# ============================================================

ggplot(
  summary_conditions,
  aes(
    x = dist_thr,
    y = RMSE_mean,
    color = validation
  )
) +
  geom_line() +
  geom_point(size = 3) +
  facet_grid(
    FOD_0 ~ log,
    labeller = labeller(
      FOD_0 = function(x)
        paste("FOD_0 =", x),
      log = function(x)
        paste("log =", x)
    )
  ) +
  theme_bw() +
  labs(
    x = "Distance threshold",
    y = "RMSE moyen",
    color = "Validation",
    title = "Sensibilité du tuning RF au data leakage"
  )


# ============================================================
# 14 - SAUVEGARDE
# ============================================================

write.csv(
  results,
  "RF_tuning_LOYO_RS_FOD_LOG_leakage.csv",
  row.names = FALSE
)

write.csv(
  summary_conditions,
  "RF_tuning_summary_conditions.csv",
  row.names = FALSE
)