# =====================================================================
# debug_xgbcv_structure.R -- affiche la VRAIE structure de xgb.cv()
# =====================================================================
# A lancer juste apres avoir cree `prep_xgb` et `scheme` (comme dans
# test_fix_xgb_tuning.R). Bypass complet de tune_xgb_early_stopping()
# pour voir exactement ce que renvoie xgb.cv() sur ta version du package,
# sans aucune couche d'abstraction.

X <- build_design_matrix(scheme$data, COVARIATES_NUM, prep_xgb$fod_levels)
y <- scheme$data[[RESPONSE_VAR]]
dall <- xgboost::xgb.DMatrix(data = X, label = y, missing = NA)
cv_folds <- map(scheme$folds, "test")

cat("Version du package xgboost installee :", as.character(packageVersion("xgboost")), "\n\n")

cv_result <- xgboost::xgb.cv(
  params = list(max_depth = 3, eta = 0.1, min_child_weight = 1,
                subsample = 0.8, colsample_bytree = 0.8,
                objective = "reg:squarederror", eval_metric = "rmse"),
  data = dall, nrounds = 50,
  folds = cv_folds, early_stopping_rounds = 20, verbose = 0
)

cat("=== class(cv_result) ===\n")
print(class(cv_result))

cat("\n=== names(cv_result) ===\n")
print(names(cv_result))

cat("\n=== cv_result$best_iteration ===\n")
print(cv_result$best_iteration)

cat("\n=== class(cv_result$evaluation_log) ===\n")
print(class(cv_result$evaluation_log))

cat("\n=== names(cv_result$evaluation_log) ===\n")
print(names(cv_result$evaluation_log))

cat("\n=== head(cv_result$evaluation_log) ===\n")
print(head(cv_result$evaluation_log))

cat("\n=== str(cv_result, max.level = 1) ===\n")
str(cv_result, max.level = 1)
