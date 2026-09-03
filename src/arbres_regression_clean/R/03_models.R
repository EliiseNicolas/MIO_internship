# =====================================================================
# 03_models.R -- interface commune CART / RF / XGB
# =====================================================================
# Chaque modèle expose la même signature :
#   fit_<model>(train_df, params, ...)      -> objet modèle
#   predict_<model>(model, newdata, ...)    -> vecteur de prédictions
#   importance_<model>(model)               -> tibble(variable, importance)
#
# Pour CART et RF, `newdata` doit être complet (pas de NA) -- filtré en
# amont par `prep_fold_data()` (voir 04_diagnostics.R).
# Pour XGB, on travaille avec des matrices de design (les NA sont
# préservés et gérés nativement par XGBoost).

# ---- Design matrix pour XGBoost (dummies manuels pour fod, NA préservés) ----
build_design_matrix <- function(data, covs_num, fod_levels) {
  mat_num <- as.matrix(data[, covs_num])
  fod_dummies <- sapply(fod_levels, function(lvl) as.numeric(data$fod == lvl))
  colnames(fod_dummies) <- paste0("fod_", trimws(fod_levels))
  cbind(mat_num, fod_dummies)
}

# =====================================================================
# CART (rpart)
# =====================================================================
# Hyperparamètres tunés : cp (complexité), minsplit, maxdepth.
fit_cart <- function(train_df, params, response = RESPONSE_VAR, covs = COVARIATES_ALL) {
  form <- stats::as.formula(paste(response, "~", paste(covs, collapse = " + ")))
  rpart::rpart(
    form, data = train_df[, c(response, covs)],
    method = "anova",
    control = rpart::rpart.control(
      cp       = params$cp,
      minsplit = params$minsplit,
      maxdepth = params$maxdepth
    )
  )
}
predict_cart <- function(model, newdata) predict(model, newdata = newdata)
importance_cart <- function(model) {
  vi <- model$variable.importance
  if (is.null(vi)) return(tibble(variable = character(), importance = numeric()))
  tibble(variable = names(vi), importance = as.numeric(vi))
}

# =====================================================================
# Random Forest (ranger)
# =====================================================================
# Hyperparamètres tunés : mtry, min.node.size, num.trees.
fit_rf <- function(train_df, params, response = RESPONSE_VAR, covs = COVARIATES_ALL) {
  form <- stats::as.formula(paste(response, "~", paste(covs, collapse = " + ")))
  ranger::ranger(
    form, data = train_df[, c(response, covs)],
    mtry = params$mtry, min.node.size = params$min.node.size,
    num.trees = params$num.trees, importance = "permutation",
    num.threads = max(1, parallel::detectCores() - 1)
  )
}
predict_rf <- function(model, newdata) predict(model, data = newdata)$predictions
importance_rf <- function(model) {
  imp <- ranger::importance(model)
  tibble(variable = names(imp), importance = as.numeric(imp))
}

# =====================================================================
# XGBoost
# =====================================================================
# Hyperparamètres tunés : max_depth, eta, min_child_weight, nrounds
# (nrounds optimal trouvé par early stopping pendant le tuning).
fit_xgb <- function(train_df, params, fod_levels,
                     covs_num = COVARIATES_NUM, response = RESPONSE_VAR) {
  X <- build_design_matrix(train_df, covs_num, fod_levels)
  y <- train_df[[response]]
  dtrain <- xgboost::xgb.DMatrix(data = X, label = y, missing = NA)
  xgb_params <- list(
    max_depth        = params$max_depth,
    eta              = params$eta,
    min_child_weight = params$min_child_weight,
    subsample        = params$subsample %||% 0.8,
    colsample_bytree = params$colsample_bytree %||% 0.8,
    objective        = "reg:squarederror"
  )
  xgboost::xgb.train(
    params = xgb_params, data = dtrain,
    nrounds = params$nrounds, verbose = 0
  )
}
predict_xgb <- function(model, newdata, fod_levels, covs_num = COVARIATES_NUM) {
  X <- build_design_matrix(newdata, covs_num, fod_levels)
  d <- xgboost::xgb.DMatrix(data = X, missing = NA)
  predict(model, d)
}
importance_xgb <- function(model) {
  imp <- xgboost::xgb.importance(model = model)
  tibble(variable = imp$Feature, importance = imp$Gain)
}

# =====================================================================
# randomForestSRC (rfsrc) -- gère nativement le manquant via imputation
# (na.action = "na.impute"), à l'entraînement ET à la prédiction. Utilisé
# pour le test de reconstruction de carte complète
# (14_run_rfsrc_reconstruction.R), sur les covariables NORMALISEES
# (COVARIATES_ALL -- ftle, Chla_total, *_totpig, fod), les mêmes que le
# pipeline principal CART/RF/XGB.
# Hyperparamètres tunés : mtry, nodesize, ntree.
# =====================================================================
fit_rfsrc <- function(train_df, params, response = RESPONSE_VAR, covs = COVARIATES_ALL) {
  form <- stats::as.formula(paste(response, "~", paste(covs, collapse = " + ")))
  randomForestSRC::rfsrc(
    form, data = train_df[, c(response, covs)],
    mtry = params$mtry, nodesize = params$nodesize, ntree = params$ntree,
    na.action = "na.impute", importance = TRUE
  )
}
predict_rfsrc <- function(model, newdata) {
  # na.action = "na.impute" est aussi utilisé en prédiction : les valeurs
  # manquantes de `newdata` sont imputées par le modèle avant prédiction
  # -> permet de prédire sur TOUS les pixels de la grille, même incomplets.
  predict(model, newdata = newdata, na.action = "na.impute")$predicted
}
importance_rfsrc <- function(model) {
  imp <- model$importance
  tibble(variable = names(imp), importance = as.numeric(imp))
}

# =====================================================================
# Backend générique : encapsule fit/predict/importance + indique si le
# modèle a besoin de données complètes (sans NA). Utilisé partout
# ailleurs (tuning, diagnostics, entraînement) pour écrire un seul code
# qui fonctionne pour les 4 familles de modèles.
# =====================================================================
make_backend <- function(model_type, fod_levels = NULL, covs = NULL) {
  switch(model_type,
    "cart" = list(
      model_type = "cart",
      fit        = fit_cart,
      predict    = function(model, newdata) predict_cart(model, newdata),
      importance = importance_cart,
      needs_complete = TRUE
    ),
    "rf" = list(
      model_type = "rf",
      fit        = fit_rf,
      predict    = function(model, newdata) predict_rf(model, newdata),
      importance = importance_rf,
      needs_complete = TRUE
    ),
    "xgb" = list(
      model_type = "xgb",
      fit        = function(train_df, params) fit_xgb(train_df, params, fod_levels = fod_levels),
      predict    = function(model, newdata) predict_xgb(model, newdata, fod_levels = fod_levels),
      importance = importance_xgb,
      needs_complete = FALSE
    ),
    "rfsrc" = list(
      model_type = "rfsrc",
      fit        = function(train_df, params) fit_rfsrc(train_df, params, covs = covs %||% COVARIATES_ALL),
      predict    = function(model, newdata) predict_rfsrc(model, newdata),
      importance = importance_rfsrc,
      needs_complete = FALSE   # gère le NA nativement via imputation
    ),
    stop("model_type inconnu : ", model_type, " (attendu : cart, rf, xgb, rfsrc)")
  )
}
