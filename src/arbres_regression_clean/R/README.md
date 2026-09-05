# Pipeline NASC -- CART / RF / XGB, naive vs blocage spatio-temporel

Refonte de tes 4 scripts en un pipeline modulaire unique, pour éviter la
duplication (le nettoyage des données, la construction des folds, les
diagnostics et les plots étaient copiés-collés dans chaque script).

## Architecture

```
data_generation/
  generate_ds_ftle_pig_fod_all_dates.R   script (corrigé) qui produit le fichier grille unique
R/
  00_config.R                  chemins, fréquences, variables, grilles de blocage
  01_data_prep.R                chargement + nettoyage (identique pour tous les modèles)
  02_folds.R                    construction des schémas : naive RS 80/20, blocage spatial, blocage temporel
  03_models.R                   interface commune fit/predict/importance pour CART, RF, XGB, randomForestSRC
  04_diagnostics.R               métriques + distances + distributions par fold (générique)
  05_plots.R                     toutes les fonctions de tracé (génériques)
  06_tuning.R                    fonction de tuning générique + grilles par défaut
  07_cross_scheme_utils.R        utilitaires pour agréger les sorties de 11_ à travers freq/modèle/schéma
  08_robustness.R                bruit gaussien sur les covariables du test + boucle par fold
  09_variogram.R                 détendançage + variogramme empirique (aide au choix des buffers)
  10_run_tuning.R                SCRIPT 1 : tuning CART/RF/XGB, par fréquence x par schéma
  11_run_training.R              SCRIPT 2 : entraînement CART/RF/XGB avec params tunés + tous les diagnostics
  12_run_grid_prediction.R       SCRIPT 3 : prédiction sur grille (une date), cartes avec/sans NA
  13_run_grid_prediction_multidate.R  SCRIPT 3bis : prédiction sur grille, TOUTES les dates (133 cartes)
  14_run_rfsrc_reconstruction.R  SCRIPT 4 : test de reconstruction de carte complète (randomForestSRC)
  15_run_transfer_tuning_test.R  SCRIPT 5 : coût de ne pas re-tuner (naive params réutilisés vs re-tunés)
  16_run_cross_scheme_analysis.R SCRIPT 6 : importance naive vs bloqué, calibration, résidu vs latitude, fuite
  17_run_map_comparison.R        SCRIPT 7 : comparaison quantitative des cartes produites (diff, corrélation, RMSE)
  18_run_noise_robustness_test.R SCRIPT 8 : robustesse au bruit gaussien sur les covariables
  19_run_variogram_analysis.R    SCRIPT 9 : variogrammes empiriques spatial/temporel des résidus détendancés
```

À exécuter dans l'ordre : `10_run_tuning.R` -> `11_run_training.R` -> `12_/13_run_grid_prediction*.R`
-> `15_/16_/18_` (ces derniers analysent des sorties déjà produites, ne réentraînent rien).
`14_run_rfsrc_reconstruction.R` est indépendant (pas besoin d'avoir lancé 10/11 avant), mais
`17_run_map_comparison.R` en a besoin (il compare ses cartes à celles de `12_`).
Chaque script `10_` à `18_` source les fichiers `00_` à `08_` dont il a besoin.

### Tout lancer d'un coup : `run_all_pipeline.R`

À la racine de `nasc_pipeline/` (à côté de `R/`), `run_all_pipeline.R` enchaîne
les 9 scripts dans l'ordre, avec logs horodatés, et vérifie entre chaque
étape que les fichiers nécessaires à l'étape suivante existent bien
(ex. arrête tout avec un message clair si le tuning ne s'est pas
terminé correctement avant de lancer le training). Chaque étape peut
être désactivée individuellement en haut du fichier (`RUN_TUNING`,
`RUN_TRAINING`, etc.) -- utile pour relancer seulement la suite après
une interruption.

Pour une exécution longue (nuit) sans bloquer ta console RStudio :
**Tools -> Background Jobs -> Start Background Job**, sélectionne
`run_all_pipeline.R`, et vérifie que le "Working Directory" du job est
bien le dossier `nasc_pipeline`. Sinon, `source("run_all_pipeline.R")`
directement dans la console fonctionne aussi (mais bloque la console
jusqu'à la fin).

**Adapter avant de lancer** : les chemins dans `00_config.R`
(`PATH_TEMPLATE`, `PATH_GRID_ALL_DATES`), qui pointent vers `F:/data_elise/...`
(chemins Windows locaux, non accessibles depuis cet environnement -- je
n'ai donc pas pu exécuter/tester ces scripts avec tes données réelles).
Vérifie aussi les packages requis : `dplyr, tidyr, purrr, tibble,
ggplot2, FNN, ranger, xgboost, rpart, rpart.plot, patchwork,
randomForestSRC`.

### Grille de prédiction : fichier unique, structure corrigée

Après ton retour sur la structure de `all_years`, j'ai identifié un vrai
bug structurel dans le script de génération : `fod` avait un ordre de
dimensions `[lon, lat, date]`, incohérent avec `ftle`/`pig` qui ont
`[date, lon, lat]`. Corrigé dans
`data_generation/generate_ds_ftle_pig_fod_all_dates.R` via `aperm()`.

Conséquence pour le pipeline : **un seul fichier grille**
(`PATH_GRID_ALL_DATES`) sert maintenant à tout -- plus besoin d'un
fichier séparé pour 2023-01-26, ni de deviner une formule de
normalisation (l'ancienne grille mono-date ne contenait que les
pigments bruts). `extract_grid_for_date()` / `extract_grid_for_date_raw()`
dans `01_data_prep.R` en extraient une date donnée, sans plus aucun cas
particulier pour `fod`.

Deux points **non corrigés silencieusement** dans le script de
génération (choix scientifiques, pas des bugs de structure) -- à
confirmer toi-même avant de faire confiance aux résultats :
1. `chla_total` est une copie exacte de `Chla` brute (DvChla exclue
   volontairement) -- ça duplique une variable déjà présente, sans info
   nouvelle en tant que prédicteur séparé.
2. Le dénominateur des ratios `*_totpig` exclut Chla de la somme, y
   compris pour calculer `chla_totpig` lui-même (`= Chla / somme des 8
   AUTRES pigments`, pas `/ somme des 9 pigments`). Une alternative
   (somme incluant Chla) est laissée en commentaire dans le script si
   ce n'était pas l'intention.

### Script 4 (nouveau) : test de reconstruction avec randomForestSRC

`14_run_rfsrc_reconstruction.R` répond à la demande "tester la
reconstruction d'une carte entière avec RF" : contrairement à RF
(`ranger`), qui doit jeter tout pixel incomplet, `randomForestSRC` gère
le manquant nativement via imputation (`na.action = "na.impute"`), à
l'entraînement ET à la prédiction. Le script entraîne un modèle sur les
**covariables normalisées** (`COVARIATES_ALL` -- `ftle`, `Chla_total`,
`*_totpig`, `fod`), les mêmes que le pipeline principal CART/RF/XGB, et
produit une carte de 2023-01-26 SANS AUCUN trou, avec un comptage
explicite du nombre de pixels "sauvés" par rapport à l'approche RF
classique (complete-case uniquement). C'est un script exploratoire
(pas de tuning complet, pas de CV multi-fold) ; utilise `tune_model()`
avec `make_backend("rfsrc")` si tu veux un vrai tuning avant de
conclure. Le backend `rfsrc` (`03_models.R`) accepte aussi un jeu de
covariables personnalisé (`make_backend("rfsrc", covs = ...)`) si tu
veux un jour comparer avec les pigments bruts (`RAW_COVARIATES_ALL`).

ggplot2, FNN, ranger, xgboost, rpart, rpart.plot, patchwork`.

### Structure de données (mise à jour)

Le jeu d'entraînement (`ds`) a changé de format : les prédicteurs sont
désormais fournis **déjà normalisés**, plus besoin de calculer des
ratios `pigment / Chla` à la main comme dans les scripts d'origine.

- `NASC` (réponse) : dérivé de la colonne source `nasc` (log10, après
  filtre 5e/95e percentile) -- **inchangé**.
- `lat`, `lon`, `time` : dérivés de `lat_nasc`, `lon_nasc`, `time_nasc`
  -- **inchangé**. (Les colonnes `*_fod`, `*_pig`, `*_ftle` sont les
  positions des capteurs sources et ne sont PAS utilisées comme
  coordonnées du point NASC.)
- `fod` : facteur, `"NA"` (chaîne) convertie en vrai `NA` -- inchangé.
- `ftle` : inchangé.
- `Chla_total` : colonne fournie telle quelle.
- `Chla_totpig, Per_totpig, But_totpig, Fuco_totpig, Hex_totpig,
  Allo_totpig, Zea_totpig, Chlb_totpig, DvChla_totpig` : colonnes
  fournies telles quelles (chaque pigment normalisé par le pigment
  total).

`COVARIATES_NUM` dans `00_config.R` a été mis à jour en conséquence.

(La partie "grille de prédiction" de cette note a été remplacée par la
section "Grille de prédiction : fichier unique, structure corrigée"
plus haut, suite à la clarification de la structure de `all_years`.)


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

### Scripts 5, 6, 7 : au-delà du RMSE global

Transposition des axes de comparaison utilisés dans un projet DL
(ViT vs CNN) à CART/RF/XGB x naive/spatial/temporel :

| Axe DL (ViT vs CNN) | Équivalent ici | Script |
|---|---|---|
| Perf en régime de données limitées | Courbe d'apprentissage par schéma (dispo. par fold très réduite en blocage fin 20x20km) | déjà dans `11_` (`01_learning_curve.png`) |
| Robustesse à l'occlusion / perturbation | Bruit gaussien croissant sur les covariables du test (fraction de l'écart-type) | `18_run_noise_robustness_test.R` |
| Coût d'un transfert moins cher (linear probing vs fine-tuning) | Hyperparamètres tunés en naive réutilisés tels quels sur les schémas bloqués, vs re-tunés | `15_run_transfer_tuning_test.R` |
| Le modèle regarde-t-il la bonne région pour la bonne raison ? | Importance des variables : naive vs bloqué (une variable qui s'effondre = raccourci d'autocorrélation) | `16_run_cross_scheme_analysis.R` |
| Biais architecturaux utiles après changement de domaine | CART/RF/XGB comparés sur le changement naive -> bloqué (le "domain shift" ici) | `16_` (RMSE vs distance) + `17_` (cartes) |
| Fuite d'autocorrélation (spécifique à ce projet) | RMSE vs distance test->train, par fold | `16_run_cross_scheme_analysis.R` |
| Calibration (spécifique à ce projet) | Moyenne prédite vs observée par bin de prédiction | `16_run_cross_scheme_analysis.R` |
| Couverture de reconstruction (spécifique à ce projet) | % de pixels prédits selon la stratégie de gestion du NA | `14_run_rfsrc_reconstruction.R` + `17_` |

**`15_run_transfer_tuning_test.R`** : pour chaque modèle x schéma bloqué,
compare le RMSE avec les hyperparamètres tunés sur naive réutilisés
tels quels, contre le RMSE avec les hyperparamètres re-tunés
spécifiquement pour ce schéma (déjà calculés par `10_`). Mêmes folds
dans les deux cas -- seul le choix des hyperparamètres change. Répond
à la question "combien perd-on à ne PAS re-tuner par schéma", posée
plus tôt dans la conversation.

**`16_run_cross_scheme_analysis.R`** : n'entraîne rien, agrège juste les
CSV déjà produits par `11_` (`metrics_par_fold.csv`, `obs_pred_all.csv`,
`importance_all.csv`) à travers tous les schémas/modèles/fréquences via
`07_cross_scheme_utils.R`. Produit :
- l'importance comparée naive vs bloqué (une variable dont l'importance
  s'effondre en blocage suggère qu'elle servait de proxy de proximité
  spatio-temporelle plutôt qu'un vrai signal) ;
- une courbe de calibration par modèle (facettée par schéma) ;
- RMSE vs distance géographique et RMSE vs distance covariables, par
  fold -- si le nuage de points montre une pente positive nette, une
  partie du score naive est un artefact de fuite plutôt qu'un vrai
  pouvoir prédictif.

**`17_run_map_comparison.R`** : compare quantitativement (pas juste
visuellement) toutes les cartes produites par `12_` et `14_` pour
`TARGET_DATE_SINGLE` -- matrice de corrélation entre TOUTES les paires
de cartes, table de RMSE par paire (sur les pixels communs), et cartes
de différence pour 4 paires illustratives (RF vs XGB même schéma, RF
naive vs RF bloqué 20x20km, XGB avec/sans NA, RF vs randomForestSRC).
Ne réentraîne et ne reprédit rien -- lit les `.rds` numériques que `12_`
et `14_` sauvegardent maintenant en plus de leurs PNG
(`predictions_all_combined.rds`, `reconstruction_grid.rds`).

**`18_run_noise_robustness_test.R`** : pour chaque modèle x schéma,
entraîne UNE FOIS par fold (hyperparamètres déjà tunés par `10_`), puis
réévalue le même modèle sur des copies de plus en plus bruitées du même
jeu de test -- bruit gaussien indépendant sur chaque covariable
numérique, d'amplitude exprimée en fraction de l'écart-type de cette
covariable (calculé sur le train, jamais sur le test). Produit une
courbe RMSE-vs-bruit et une courbe de dégradation relative (%) par
modèle et par schéma, plus une comparaison directe CART/RF/XGB sur le
schéma naive. Répond à "quel modèle se dégrade le moins vite sous
perturbation" -- et permet de vérifier si un modèle robuste au bruit
l'est aussi au changement de domaine spatial/temporel (comparaison avec
les résultats de `16_`).

**Non implémenté (occlusion réelle, taux de NA artificiel)** : le script
`18_` couvre le bruit gaussien continu (perturbation de texture), mais
pas un test de type "masquer X% des pixels/covariables et mesurer la
chute" (occlusion dure plutôt que bruit continu). Dis-le moi si tu
veux ce test en complément -- il se brancherait sur le même principe
(`08_robustness.R`), en remplaçant l'ajout de bruit par une mise à `NA`
aléatoire d'une fraction des covariables du test.

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

## Importance des variables XGBoost : Gain vs SHAP

`importance_xgb()` (`03_models.R`, branchée dans `make_backend("xgb")`,
utilisée partout dans le pipeline) utilise **Gain** : la réduction
moyenne de la fonction de perte apportée par les splits sur chaque
variable. C'est le choix standard, mais il a deux limites : il ne dit
rien sur la **direction** de l'effet, et il répartit arbitrairement le
crédit entre variables corrélées (assez probable ici entre `Chla_total`
et les `*_totpig`).

**`importance_xgb_shap()`** (même fichier) calcule une alternative basée
sur les **valeurs SHAP** (`predict(model, newdata, predcontrib = TRUE)`) :
décomposition additive de chaque prédiction individuelle en contribution
par variable, avec signe. L'importance globale retournée est la moyenne
de la valeur absolue de ces contributions par variable.

Différence de signature avec `importance_xgb()` : SHAP a besoin de
données sur lesquelles calculer les contributions (`newdata`), pas
seulement du modèle -- elle n'est donc **pas branchée automatiquement**
dans `make_backend("xgb")` (l'interface générique `importance(model)`
est commune aux 4 backends). Pour l'utiliser à travers tous les folds
d'un schéma (agrégée, compatible avec `plot_importance_mean_sd()`
directement) :

```r
shap_imp <- compute_shap_importance_across_folds(scheme_xgb, best_params, fod_levels)
plot_importance_mean_sd(shap_imp, subtitle = "XGB - SHAP")
```

(`compute_shap_importance_across_folds()`, dans `04_diagnostics.R`,
calcule les SHAP sur le **test** de chaque fold, pas le train -- on veut
savoir sur quoi le modèle s'appuie pour ses prédictions hors échantillon.)

Pas encore branchée automatiquement dans `11_run_training.R` ni
`16_run_cross_scheme_analysis.R` (uniquement les fonctions de base pour
l'instant, comme demandé) -- dis-moi si tu veux que je l'y intègre en
complément de Gain (par exemple un fichier `importance_shap_all.csv` en
plus de `importance_all.csv`, et un plot comparatif Gain vs SHAP).

## Naive RS 80/20 : K-fold répété au lieu de Monte-Carlo (mise à jour)

Suite à la discussion sur Monte-Carlo CV vs K-fold (partition), le
naive RS 80/20 utilise maintenant par défaut du **K-fold répété**
(`build_naive_folds_kfold_repeated()`, `02_folds.R`) plutôt que des
tirages Monte-Carlo :

- `K` est dérivé de `frac_train` (`K = round(1 / (1 - frac_train))` --
  `K = 5` pour `frac_train = 0.8`, chaque fold test = 1/5 = 20%).
- Avec `NAIVE_N_FOLDS = 10`, ça donne **2 répétitions de 5-fold** :
  chaque répétition est une vraie partition disjointe (chaque
  observation testée une fois par répétition), avec une nouvelle
  permutation aléatoire à chaque répétition.
- Avantages sur Monte-Carlo : couverture complète et disjointe (pas
  d'observation jamais testée, pas de redondance), pas de
  sous-estimation de la variance inter-fold (les tests Monte-Carlo se
  chevauchaient), et une logique de partition cohérente avec les
  schémas bloqués.

**Basculer entre les deux méthodes** : `NAIVE_CV_METHOD` dans
`00_config.R` (`"kfold_repeated"` par défaut, `"monte_carlo"` pour
retrouver l'ancien comportement -- utile si tu veux comparer les deux
directement, ou reproduire d'anciens résultats).

`build_naive_folds()` reste l'interface publique inchangée (utilisée
partout ailleurs dans le pipeline) -- elle route juste vers l'une ou
l'autre des deux implémentations selon `NAIVE_CV_METHOD`. Aucun autre
fichier n'a besoin d'être modifié.

## Gradient latitudinal : 3 ajouts (diagnostic, blocage, buffers)

Suite à la discussion sur la stationnarité de 1er/2nd ordre :

**1. Plot résidu vs latitude** (`plot_residual_vs_latitude()`,
`05_plots.R`, branché dans `16_run_cross_scheme_analysis.R`) : nuage de
points + lissage LOESS du résidu (observé - prédit) en fonction de la
latitude, coloré par schéma, facetté par modèle. Utilise
`obs_pred_all.csv` déjà produit par `11_` -- aucun réentraînement.
Lecture : une courbe lissée plate autour de 0 -> le gradient latitudinal
est bien absorbé par les covariables. Une pente ou une forme
systématique -> biais résiduel lié à la latitude non modélisé,
généralement plus visible en blocage qu'en naive (la CV naive le masque
partiellement par fuite spatiale). Sortie : `residual_vs_latitude.png`.

**2. Blocage spatial stratifié par latitude** (`stratify_blocks_by_latitude()`,
`02_folds.R`, appelée automatiquement dans `build_spatial_folds()` dès
qu'il y a plus de blocs disponibles que de folds cibles) : au lieu d'un
tirage uniforme parmi tous les blocs éligibles, les blocs disponibles
sont d'abord répartis en 3 tertiles de latitude (bas/milieu/haut,
calculés sur la latitude moyenne de chaque bloc), puis le tirage se
fait à quotas égaux dans chacun des 3 tertiles (avec repli sur le pool
restant si un tertile manque de blocs). Le message console de
`build_spatial_folds()` indique maintenant explicitement "stratifiés
par latitude". Objectif : que la comparaison inter-folds (RMSE, etc.)
ne soit pas biaisée par le fait qu'un fold test se trouve, par hasard,
au centre ou au bord de la plage de latitude couverte.

**3. Variogramme des résidus détendancés** (`09_variogram.R` +
`19_run_variogram_analysis.R`) : `detrend_residuals()` retire une
tendance simple (`NASC ~ lat` par défaut -- juste une régression
linéaire, pas un modèle complet, pour isoler la structure spatiale
résiduelle propre sans mélanger "le modèle a mal appris" et "il reste
de la vraie autocorrélation locale") ; `compute_empirical_variogram()`
calcule ensuite la semi-variance par classe de distance sur un
sous-échantillon (le calcul est en O(n²), donc sous-échantillonné à
`n_sample = 3000` points par défaut -- ajustable). Le script produit un
variogramme spatial (km) et un variogramme temporel (jours), avec les
buffers actuels (`SPATIAL_RESOLUTIONS`, `TEMPORAL_RESOLUTIONS` dans
`00_config.R`) superposés en pointillés pour comparaison visuelle
directe : si le buffer tombe bien avant que la courbe plafonne, il est
probablement trop court (fuite résiduelle au-delà du buffer) ; s'il
tombe bien après, il est probablement trop large. Ne réentraîne aucun
modèle -- calcul purement descriptif sur les données d'entraînement.

## Échelles partagées entre plots (axes Y / colorbars identiques)

Pour comparer facilement modèles/schémas/fréquences, `SHARED_SCALE_SCOPE`
dans `00_config.R` (`"global"` par défaut, `"per_freq"` en alternative --
voir le commentaire dans le fichier pour le compromis) contrôle
l'étendue sur laquelle les échelles sont calculées.

**Ce qui est maintenant partagé** (même axe/colorbar sur tous les PNG
concernés) :
- **`11_run_training.R`** (restructuré en 2 phases : entraînement +
  collecte légère, puis calcul des échelles communes et génération de
  tous les plots) : RMSE (learning curve + barres par fold), R²,
  variance intra-fold, axes obs/pred, couleur des cartes de résidus
  (signée et absolue).
- **`12_/13_/14_`** (cartes de prédiction) : même échelle de couleur
  viridis = union de la plage de NASC **observée à l'entraînement** et
  de la plage **prédite**, partagée entre RF/XGB/randomForestSRC, tous
  schémas, toutes dates (`13_`).
- **`17_run_map_comparison.R`** : les 4 cartes de différence partagent
  une échelle symétrique commune (`max(|diff|)` calculé sur les 4 en
  même temps, avant de tracer).
- **`18_run_noise_robustness_test.R`** (restructuré en 2 phases) : RMSE
  absolu et dégradation relative (%) partagés entre CART/RF/XGB et tous
  les schémas -- légitime ici car c'est la même unité pour les 3 modèles.

**Ce qui reste volontairement NON partagé** (exceptions documentées, pas
des oublis) :
- **Importance des variables** (`07_importance_variables.png`) : partagée
  **au sein d'un même modèle** (entre schémas/fréquences), mais **jamais
  entre modèles différents** -- le Gain XGBoost et l'importance par
  permutation RF/CART ne sont pas la même unité ; forcer le même axe
  serait trompeur, pas juste une question d'esthétique.
- **Plots de métadonnées** (`11_distances_test_train.png`,
  `12_covariable_stats.png`, `13_covariable_variance.png`,
  `14_fod_distribution.png`, `15_distributions_par_fold.png`) : chaque
  config s'inspecte individuellement pour du contrôle qualité, pas en
  comparaison directe -- laissés auto-scalés. Dis-le moi si tu veux
  qu'ils soient unifiés aussi.
- **`timeseries_mean_prediction.png`** (`13_`, évolution temporelle de la
  moyenne spatiale prédite) : axe Y auto-scalé par fréquence/schéma --
  c'est une moyenne spatiale, pas une valeur pixel, donc moins
  directement comparable à la même échelle que les cartes ; RF vs XGB
  reste comparable *au sein* de chaque figure (même axe, deux courbes).
