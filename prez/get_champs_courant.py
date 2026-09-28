"""
Telechargement de champs de courants CMEMS (uo, vo) pour le 2023-01-26
sur le secteur Sud-Ouest de l'ocean Indien : lon [45, 90] E, lat [-60, -30] N
puis trace d'une carte des courants (fond = vitesse, fleches = direction/intensite)

Prerequis :
    pip install copernicusmarine matplotlib cartopy numpy xarray

Identifiants CMEMS (compte gratuit sur https://data.marine.copernicus.eu/) :
    - soit en lancant une fois "copernicusmarine login" dans un terminal
      (les identifiants sont ensuite stockes localement)
    - soit via variables d'environnement :
        COPERNICUSMARINE_SERVICE_USERNAME
        COPERNICUSMARINE_SERVICE_PASSWORD
"""

import copernicusmarine
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# ---- 1. Parametres ----
DATASET_ID = "cmems_mod_glo_phy_my_0.083deg_P1D-m"  # reanalyse GLORYS12, couvre 2021
VARIABLES = ["uo", "vo"]                             # composantes du courant (est / nord)
DATE = "2023-01-26"

LON_MIN, LON_MAX = 45, 90
LAT_MIN, LAT_MAX = -60, -30

OUTPUT_DIR = "./cmems_data"
OUTPUT_FILE = "currents_2023-01-26_SO_indien.nc"

# ---- 2. Telechargement (subset direct depuis le catalogue CMEMS) ----
copernicusmarine.subset(
    dataset_id=DATASET_ID,
    variables=VARIABLES,
    minimum_longitude=LON_MIN,
    maximum_longitude=LON_MAX,
    minimum_latitude=LAT_MIN,
    maximum_latitude=LAT_MAX,
    start_datetime=f"{DATE}T00:00:00",
    end_datetime=f"{DATE}T00:00:00",
    minimum_depth=0,
    maximum_depth=1,          # surface uniquement
    output_directory=OUTPUT_DIR,
    output_filename=OUTPUT_FILE,
    force_download=True,
)

# ---- 3. Lecture du fichier telecharge ----
ds = xr.open_dataset(f"{OUTPUT_DIR}/{OUTPUT_FILE}")

uo = ds["uo"].squeeze()   # enleve les dimensions de taille 1 (time, depth)
vo = ds["vo"].squeeze()
lon = ds["longitude"].values
lat = ds["latitude"].values

speed = np.sqrt(uo**2 + vo**2)

# ---- 4. Carte des courants ----
fig = plt.figure(figsize=(10, 8))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([LON_MIN, LON_MAX, LAT_MIN, LAT_MAX], crs=ccrs.PlateCarree())

# Fond de carte : vitesse du courant (m/s)
mesh = ax.pcolormesh(
    lon, lat, speed,
    transform=ccrs.PlateCarree(), cmap="viridis", shading="auto"
)
cbar = plt.colorbar(mesh, ax=ax, orientation="vertical", pad=0.05, shrink=0.8)
cbar.set_label("Vitesse du courant (m/s)")

# Vecteurs de courant (sous-echantillonnes pour la lisibilite)
skip = 4
ax.quiver(
    lon[::skip], lat[::skip],
    uo.values[::skip, ::skip], vo.values[::skip, ::skip],
    transform=ccrs.PlateCarree(), color="white", scale=15, width=0.0025
)

ax.add_feature(cfeature.LAND, facecolor="lightgray", zorder=2)
ax.add_feature(cfeature.COASTLINE, zorder=3)
ax.gridlines(draw_labels=True, linewidth=0.3, color="gray", alpha=0.5)

ax.set_title(f"Champ de courants de surface - {DATE}\nSecteur Sud-Ouest de l'ocean Indien")

plt.tight_layout()
plt.savefig("carte_courants_2023-01-26.png", dpi=150)
plt.show()

print("Carte enregistree : carte_courants_2023-01-26.png")