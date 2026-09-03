# Pipeline NASC -- CART / RF / XGB, naive vs blocage spatio-temporel

Refonte de tes 4 scripts en un pipeline modulaire unique, pour éviter la
duplication (le nettoyage des données, la construction des folds, les
diagnostics et les plots étaient copiés-collés dans chaque script).

## Architecture

```
R/
  00_config.R              chemins, fréquences, variables, grilles de blocage
  01_data_prep.R           chargement + nettoyage (identique pour tous les modèles)
  02_folds.R               construction des schémas : naive RS 80/20, blocage spatial, blocage temporel
  03_models.R              interface commune fit/predict/importance pour CART, RF, XGB
  04_diagnostics.R         métriques + distances + distributions par fold (générique)
  05_plots.R               toutes les fonctions de tracé (génériques)
  06_tuning.R              fonction de tuning générique + grilles par défaut
  10_run_tuning.R          SCRIPT 1 : tuning CART/RF/XGB, par fréquence x par schéma
  11_run_training.R        SCRIPT 2 : entraînement 10 folds RF/XGB avec params tunés + tous les diagnostics
  12_run_grid_prediction.R SCRIPT 3 : prédiction sur grille = moyenne des 10 modèles, cartes avec/sans NA
```

À exécuter dans l'ordre : `10_run_tuning.R` -> `11_run_training.R` -> `12_run_grid_prediction.R`.
Chaque script `10_/11_/12_` source les fichiers `00_` à `06_` dont il a besoin.

**Adapter avant de lancer** : les chemins dans `00_config.R`
(`PATH_TEMPLATE`, `PATH_GRID_DAY`), qui pointent vers `F:/data_elise/...`
(chemins Windows locaux, non accessibles depuis cet environnement -- je
n'ai donc pas pu exécuter/tester ces scripts avec tes données réelles).
Vérifie aussi les packages requis : `dplyr, tidyr, purrr, tibble,
ggplot2, FNN, ranger, xgboost, rpart, rpart.plot, patchwork`.

---

## 1. Faut-il un tuning différent selon le type de split ?

**Oui, clairement.** Les hyperparamètres doivent être choisis avec le
**même schéma de validation** que celui utilisé ensuite pour évaluer /
déployer le modèle :

- Un modèle tuné en **CV naive** (split aléatoire) "voit" des points
  d'entraînement très proches en espace/temps des points de test
  (fuite d'information locale). Le tuning va donc naturellement
  sélectionner des hyperparamètres **peu régularisés** (RF : `mtry`
  élevé, `min.node.size` faible, beaucoup d'arbres ; XGB : `max_depth`
  élevé, `eta` plus grand) puisque le modèle peut se permettre
  d'"apprendre par cœur" du signal local qui ne généralisera pas.
- En **blocage spatial/temporel**, cette fuite est supprimée (buffer),
  donc le même niveau de complexité va généralement **sur-apprendre**
  et donner un RMSE de test bien plus mauvais. Le tuning en blocage va
  au contraire pousser vers des hyperparamètres **plus conservateurs**
  (plus de régularisation), ce qui est justement ce qu'on veut pour une
  vraie tâche d'extrapolation spatiale/temporelle.
- Idéalement, on va même plus loin : **retuner séparément pour chaque
  résolution de blocage** (20x20km / 200x200km / 1500x1000km / 1j),
  car le volume de train disponible par fold et la structure de
  corrélation spatiale ne sont pas les mêmes à ces échelles. C'est ce
  que fait `10_run_tuning.R` : une tuning run par (fréquence x modèle x
  schéma), soit 2 fréquences x 3 modèles x 5 schémas = 30 tunings.

Sans ça, le risque classique est de publier un R² qui a l'air
excellent (obtenu avec des hyperparamètres tunés en naive) alors qu'il
s'effondre complètement dès qu'on regarde la performance en blocage
spatial réel -- ce qui fausserait la comparaison "naive vs blocage" que
tu veux faire, en mélangeant l'effet du schéma de CV et l'effet du
tuning.

---

## 2. Autres plots utiles pour l'entraînement (en plus de la learning curve)

En plus de la courbe d'apprentissage (train/test RMSE vs fraction de
train), le pipeline produit par défaut :

- **Courbe/heatmap de validation** (`plot_validation_curve`) : RMSE
  test en fonction de chaque hyperparamètre, pendant le tuning --
  permet de voir si on est sous-régularisé (RMSE train << RMSE test) ou
  sur-régularisé (RMSE train ~ RMSE test mais les deux élevés).
- **RMSE et R² par fold** -- repère les folds "difficiles" (souvent
  ceux avec le moins de train après buffer, ou une zone géographique
  atypique).
- **Variance intra-fold du NASC observé (test)** -- un R² faible sur un
  fold à faible variance n'a pas le même sens qu'un R² faible sur un
  fold à forte variance.
- **Indice d'extrapolation** (`mean_covariate_dist` train-test /
  `mean_covariate_dist` intra-train) -- déjà calculé dans les
  métriques, à tracer par fold : détecte les folds où le test est
  "hors distribution" par rapport au train.
- **Stabilité de l'importance des variables** entre folds (déjà en
  moyenne +/- écart-type) -- une variable importante seulement dans 1
  fold sur 10 est suspecte (surapprentissage sur une particularité
  locale).
- Pour CART spécifiquement : envisage un `rpart.plot()` de l'arbre
  final (interprétabilité) -- pas inclus par défaut ici mais trivial à
  ajouter (`rpart.plot::rpart.plot(model)`).
- Optionnel si tu veux aller plus loin : **courbes de calibration**
  (moyenne prédite vs moyenne observée par décile de prédiction) et
  **PDP/ICE** sur les 2-3 variables les plus importantes, utiles pour
  vérifier que le modèle apprend une relation physiquement plausible
  et pas un artefact.

Tous ces plots sont déjà générés automatiquement par
`11_run_training.R` (voir la liste numérotée dans le dossier de sortie
de chaque configuration).

---

## 3. Ce que produit chaque script

### `10_run_tuning.R`
Pour chaque fréquence x schéma (naive, 3 résolutions spatiales, 1
résolution temporelle) : tuning CART (`cp`, `minsplit`, `maxdepth`), RF
(`mtry`, `min.node.size`, `num.trees`), XGB (`max_depth`, `eta`,
`min_child_weight`, `nrounds` via early stopping). Sauvegarde
`outputs_pipeline/tuning/<modele>_<freq>kHz_<schema>.rds` (+ .png).

### `11_run_training.R`
Pour chaque fréquence x CART/RF/XGB x schéma : charge les
hyperparamètres tunés, entraîne un modèle par fold, et sauvegarde
sous `outputs_pipeline/training/<freq>kHz/<modele>/<schema>/` :

> CART et RF partagent exactement le même schéma/folds (`schemes_rf`,
> données complètes -- CART a lui aussi besoin de données sans NA dans
> ce pipeline). XGB utilise `schemes_xgb` (NA préservés). Le nombre de
> folds pour le naive est fixe (`NAIVE_N_FOLDS = 10`) ; pour le blocage
> spatial/temporel il est **adaptatif** (cf. section 5 ci-dessous), donc
> peut différer d'une résolution à l'autre -- tout le code (plots,
> agrégations) boucle sur `names(scheme$folds)`, pas sur un nombre fixe,
> donc ça ne casse rien.

- `models.rds` (les 10 modèles), `metrics_par_fold.csv`,
  `obs_pred_all.csv`, `importance_all.csv`
- `01_learning_curve.png`
- `02/03/04_obs_vs_pred_*.png` (nuage de points, histogramme, par fold)
- `05/06_rmse_r2_par_fold.png`
- `07_importance_variables.png` (moyenne +/- écart-type)
- `08/09_carte_residus*.png`
- `10_variance_intra_fold.png`
- `11_distances_test_train.png` (géographique + covariables)
- `12/13_covariable_stats/variance.png` (train vs test, par fold)
- `14_fod_distribution.png`
- `15_distributions_par_fold.png` (toutes les covariables numériques)
- `16_carte_points_par_fold.png`
- `17_carte_train_test_buffer.png` (uniquement pour les schémas bloqués
  -- localisation spatiale test / train conservé / exclu par le buffer)

### `12_run_grid_prediction.R`
Pour chaque fréquence x schéma : recharge les 10 modèles RF/XGB,
prédit sur la grille du jour choisi, moyenne les 10 prédictions :

- RF : une carte, pixels complets uniquement (pas de gestion native du
  NA), avec le nombre de pixels prédits / total en sous-titre.
- XGB : deux cartes côte à côte (avec NA gérés nativement / sans NA),
  même échelle de couleur pour comparaison directe.
- Une carte de comparaison globale (facettée modèle x schéma) par
  fréquence.

Les points 1 ("réel vs prédit") et 2 ("importance moyenne") de ta liste
sont déjà couverts par `11_run_training.R` à partir des prédictions
out-of-fold (la grille de prédiction n'a pas de NASC observé, donc pas
de "réel" à comparer à cet endroit).

---

## 4bis. Nombre de folds adaptatif pour le blocage (mise à jour)

Le nombre de folds pour les schémas bloqués (spatial et temporel) n'est
plus fixé à 10 pour toutes les résolutions : il est maintenant calculé
comme une **fraction des blocs disponibles** (30% par défaut), bornée
par un plancher et un plafond :

```r
BLOCK_MAX_FOLDS_FRACTION <- 0.3   # fraction des blocs disponibles utilisée
BLOCK_MIN_FOLDS          <- 5     # plancher
BLOCK_MAX_FOLDS_ABS      <- 30    # plafond
```

Pourquoi : à 20x20km il peut y avoir des centaines de blocs -- n'en
tirer que 10 jetterait beaucoup d'information et donnerait une
variance inter-fold peu fiable ; à 1500x1000km il peut n'y avoir que
quelques blocs -- viser 10 n'aurait pas de sens. Avec cette règle, une
résolution fine se retrouve avec jusqu'à 30 folds, une résolution
grossière avec 5 folds minimum (ou moins si vraiment trop peu de blocs
ont >= `BLOCK_MIN_BLOCK_N` observations, auquel cas un message
`[!]`-warning s'affiche). Le naive RS 80/20 reste fixé à
`NAIVE_N_FOLDS = 10`, car ce n'est pas un nombre de blocs mais un
nombre de répétitions Monte-Carlo choisi arbitrairement.

Conséquence pour l'interprétation : ne compare pas directement "RMSE
moyen sur 10 folds (naive)" à "RMSE moyen sur 27 folds (20x20km)"
comme si c'était la même unité statistique -- regarde plutôt la
distribution (boxplot / barres d'erreur) et le nombre de folds affiché
dans chaque `resume_global.txt`.

## 4. Autres éléments auxquels tu n'avais pas pensé (proposition)

- **Comparaison directe des schémas** : un plot RMSE/R² par schéma
  (naive vs chaque résolution de blocage), barres d'erreur =
  écart-type inter-fold -- pour visualiser en un coup d'œil l'écart de
  performance "optimiste" (naive) vs "réaliste" (blocage). Facile à
  reconstituer à partir de `metrics_par_fold.csv` de chaque config.
- **Nombre de folds réellement utilisés vs cible** : pour les
  résolutions de blocage fines (20x20km), il se peut qu'il n'y ait pas
  assez de blocs avec >= 50 observations pour atteindre 10 folds --
  le pipeline te prévient dans la console (`[!] seulement N blocs
  disponibles...`) mais pense à vérifier ce message avant d'interpréter
  les résultats de cette résolution.
- **Choix du buffer** : les buffers par défaut (100km / 20km / 5km
  pour le spatial, 1 jour pour le temporel) sont repris de tes scripts
  d'origine mais sont un peu arbitraires -- un variogramme empirique
  des résidus (portée de l'autocorrélation spatiale/temporelle)
  donnerait une justification plus solide que de choisir le buffer "au
  jugé".
- **Cohérence CART inclus dans le tuning mais pas dans l'entraînement
  10 folds** : tu n'as pas demandé de CART en entraînement/prédiction
  finale (seulement en tuning), donc `11_` et `12_` ne concernent que
  RF et XGB. Si tu veux aussi un CART "final" comme baseline
  interprétable, il suffit de dupliquer le bloc RF dans
  `11_run_training.R` en remplaçant `make_backend("rf")` par
  `make_backend("cart")`.
