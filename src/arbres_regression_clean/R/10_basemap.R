# =====================================================================
# 10_basemap.R -- fond de carte (continents + îles) pour les cartes de
# prédiction, avec cache local
# =====================================================================
# Nécessite les packages `sf` et `rnaturalearth` (+ `rnaturalearthdata`,
# installé automatiquement comme dépendance de rnaturalearth pour la
# couche "pays" à résolution moyenne). install.packages(c("sf",
# "rnaturalearth", "rnaturalearthdata")) si absent.
#
# Les Kerguelen ne sont PAS incluses dans la couche "pays" standard
# (trop petites) -- récupérées séparément via la couche Natural Earth
# "minor_islands" (ne_download(), nécessite un accès internet au
# PREMIER appel seulement ; mis en cache ensuite dans
# outputs_pipeline/basemap/ pour ne plus jamais retélécharger).
#
# En cas d'échec (pas d'accès internet, packages absents), la fonction
# retourne NULL et un warning explicite -- les cartes de prédiction
# continuent de fonctionner SANS fond de carte plutôt que de planter.
#
# IMPORTANT -- ce qui est LENT ici n'est PAS le téléchargement (quelques
# secondes) mais la fusion (st_union) des polygones de pays/îles à
# l'échelle MONDIALE, une opération géométrique coûteuse sans aucun
# indicateur de progression -- facile à confondre avec un blocage. On
# recadre donc les couches sur la zone d'étude (BASEMAP_BBOX) AVANT de
# les fusionner : la fusion ne porte plus que sur une poignée de
# polygones au lieu de ~250 pays + ~15000 îles, ce qui la rend quasi
# instantanée. Adapte BASEMAP_BBOX si ta zone d'étude change.

BASEMAP_CACHE_PATH <- file.path("outputs_pipeline", "basemap", "coastline_with_islands.rds")
BASEMAP_DOWNLOAD_TIMEOUT_SEC <- 15

# Boîte englobante (lon_min, lon_max, lat_min, lat_max) utilisée pour
# recadrer le fond de carte -- généreuse par rapport à la zone
# d'étude (~45-80°E, -55 à -30°N) pour ne rien couper aux bords des cartes.
BASEMAP_BBOX <- c(xmin = 30, xmax = 100, ymin = -60, ymax = -20)

# URL directe du shapefile "minor_islands" (contient les Kerguelen) --
# téléchargée manuellement avec download.file() de base R (mode "wininet"
# sous Windows) plutôt que via rnaturalearth::ne_download() (mecanisme
# "curl" different, qui peut rester bloque derriere certains proxy
# d'entreprise alors que download.file() de base passe sans probleme --
# constate en pratique sur ce projet).
BASEMAP_ISLANDS_URL <- "https://naciscdn.org/naturalearth/10m/physical/ne_10m_minor_islands.zip"

load_basemap_sf <- function(force_refresh = FALSE) {
  if (exists("BASEMAP_ENABLED") && isFALSE(BASEMAP_ENABLED)) {
    return(NULL)  # desactive explicitement -- aucune tentative de reseau
  }
  if (!force_refresh && file.exists(BASEMAP_CACHE_PATH)) {
    return(readRDS(BASEMAP_CACHE_PATH))
  }

  if (!requireNamespace("sf", quietly = TRUE) || !requireNamespace("rnaturalearth", quietly = TRUE)) {
    warning(
      "Packages 'sf' et/ou 'rnaturalearth' non installes -- pas de fond de ",
      "carte (continents/iles) sur les cartes de prediction. Installe avec : ",
      "install.packages(c('sf', 'rnaturalearth', 'rnaturalearthdata'))"
    )
    return(NULL)
  }
  suppressPackageStartupMessages(library(sf))  # necessaire pour que geom_sf() dispatche correctement

  land <- tryCatch({
    # bbox construite avec [[ ]] (valeur scalaire "nue") plutot que [ ]
    # (qui garde le nom "xmin" et provoquait un nom compose "xmin.xmin"
    # une fois passe dans c(xmin = ...) -- st_bbox() ne reconnaissait
    # alors plus les noms attendus, d'ou l'erreur "!anyNA(x) n'est pas
    # TRUE" en aval).
    bbox_sf <- sf::st_as_sfc(sf::st_bbox(c(
      xmin = BASEMAP_BBOX[["xmin"]], ymin = BASEMAP_BBOX[["ymin"]],
      xmax = BASEMAP_BBOX[["xmax"]], ymax = BASEMAP_BBOX[["ymax"]]
    ), crs = sf::st_crs(4326)))

    cat("Fond de carte : chargement des polygones de pays (donnees locales, pas de reseau)...\n")
    countries <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

    cat("Fond de carte : telechargement manuel des petites iles (Kerguelen incluses)...\n")
    zip_path <- tempfile(fileext = ".zip")
    dl_ok <- tryCatch({
      utils::download.file(BASEMAP_ISLANDS_URL, destfile = zip_path,
                            mode = "wb", quiet = TRUE, method = "auto",
                            timeout = BASEMAP_DOWNLOAD_TIMEOUT_SEC)
      TRUE
    }, error = function(e) FALSE)

    if (!dl_ok || !file.exists(zip_path) || file.size(zip_path) == 0) {
      stop("Telechargement manuel des iles echoue (voir BASEMAP_ISLANDS_URL).")
    }

    unzip_dir <- tempfile()
    dir.create(unzip_dir)
    utils::unzip(zip_path, exdir = unzip_dir)
    shp_file <- list.files(unzip_dir, pattern = "\\.shp$", full.names = TRUE)[1]
    if (is.na(shp_file)) stop("Aucun fichier .shp trouve dans l'archive des iles telechargee.")
    islands <- sf::st_read(shp_file, quiet = TRUE)

    cat("Fond de carte : recadrage sur la zone d'etude (evite une fusion mondiale, lente)...\n")
    countries_crop <- suppressWarnings(sf::st_crop(sf::st_make_valid(countries), bbox_sf))
    islands_crop   <- suppressWarnings(sf::st_crop(sf::st_make_valid(islands),   bbox_sf))

    cat("Fond de carte : fusion des polygones...\n")
    combined <- sf::st_union(sf::st_geometry(countries_crop), sf::st_geometry(islands_crop))
    sf::st_as_sf(combined)
  }, error = function(e) {
    warning(
      "Echec du chargement du fond de carte : ", conditionMessage(e),
      ". Les cartes seront generees SANS fond de carte."
    )
    NULL
  })

  if (!is.null(land)) {
    dir.create(dirname(BASEMAP_CACHE_PATH), showWarnings = FALSE, recursive = TRUE)
    saveRDS(land, BASEMAP_CACHE_PATH)
    cat("Fond de carte pret et mis en cache dans :", BASEMAP_CACHE_PATH, "\n")
  }
  land
}

