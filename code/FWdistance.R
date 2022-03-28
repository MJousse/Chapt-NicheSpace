library(dplyr)
library(sf)
library(raster)

r <- getData("worldclim",var="bio",res=10, path = "data/raw/worldclim")
r <- r[[c(1, 5, 6, 12, 13, 14)]]
names(r) <- c("Tmean", "Tmax", "Tmin", "PrecMean", "PrecMax", "PrecMin")

# Europe
Europe <- st_read("data/raw/polygons/BiogeoRegions2016_shapefile/BiogeoRegions2016.shp") %>%
  filter(short_name != "outside")
Europe_pt <- do.call(rbind, st_sample(Europe, 100)) %>%
  as_tibble() %>% setNames(c("lon","lat")) %>%
  SpatialPoints(proj4string = CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"))

Europe_clims <- raster::extract(r, Europe_pt)

# Pyrenees
Pyrenees <- st_read("data/raw/polygons/EuropeanMountainAreas/m_massifs_v1.shp") %>%
  filter(name_mm == "Pyrenees")
Pyrenees_pt <- do.call(rbind, st_sample(Pyrenees, 100)) %>%
  as_tibble() %>% setNames(c("lon","lat")) %>%
  SpatialPoints(proj4string = CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"))

Pyrennees_clims <- raster::extract(r, Pyrenees_pt)

# Serengeti
Serengeti <- st_read("data/raw/polygons/Serengeti_Ecosystem/v3_serengeti_ecosystem.shp")
Serengeti_pt <- do.call(rbind, st_sample(Serengeti, 100)) %>%
  as_tibble() %>% setNames(c("lon","lat")) %>%
  SpatialPoints(proj4string = CRS("+proj=utm +zone=36 +south +ellps=clrk80 +units=m +no_defs"))

Serengeti_clims <- raster::extract(r, Serengeti_pt)

# High Arctic
HighArctic <- st_read("data/raw/polygons/CPCAD-BDCAPC_Dec2021.gdb/", layer = "CPCAD_BDCAPC_Dec2021") %>%
  filter(NAME_E == "Bylot Island Bird Sanctuary", BIOME == "T")
HighArctic_pt <- do.call(rbind, st_sample(HighArctic, 100)) %>%
  as_tibble() %>% setNames(c("lon","lat")) %>%
  SpatialPoints(proj4string = CRS("+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"))

HighArctic_clims <- raster::extract(r, HighArctic_pt)

# Environmental distance

