library(dplyr)
library(sf)
library(raster)
library(ggbiplot)

r <- getData("worldclim",var="bio",res=10, path = "data/raw/worldclim")
r <- r[[c(1, 5, 6, 12, 13, 14)]]
names(r) <- c("Tmean", "Tmax", "Tmin", "PrecMean", "PrecMax", "PrecMin")

# Europe
Europe <- st_read("data/raw/polygons/BiogeoRegions2016_shapefile/BiogeoRegions2016.shp") %>%
  filter(short_name != "outside")
Europe_pt <- do.call(rbind, st_sample(Europe, 100)) %>%
  as_tibble() %>% setNames(c("lon","lat")) %>%
  SpatialPoints(proj4string = CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs")) %>%
  spTransform(CRS("+init=epsg:3857"))

Europe_clims <- raster::extract(r, Europe_pt)

# Pyrenees
Pyrenees <- st_read("data/raw/polygons/EuropeanMountainAreas/m_massifs_v1.shp") %>%
  filter(name_mm == "Pyrenees")
Pyrenees_pt <- do.call(rbind, st_sample(Pyrenees, 100)) %>%
  as_tibble() %>% setNames(c("lon","lat")) %>%
  SpatialPoints(proj4string = CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs")) %>%
  spTransform(CRS("+init=epsg:3857"))

Pyrenees_clims <- raster::extract(r, Pyrenees_pt)

# Serengeti
Serengeti <- st_read("data/raw/polygons/Serengeti_Ecosystem/v3_serengeti_ecosystem.shp")
Serengeti_pt <- do.call(rbind, st_sample(Serengeti, 100)) %>%
  as_tibble() %>% setNames(c("lon","lat")) %>%
  SpatialPoints(proj4string = CRS("+proj=utm +zone=36 +south +ellps=clrk80 +units=m +no_defs")) %>%
  spTransform(CRS("+init=epsg:3857"))

Serengeti_clims <- raster::extract(r, Serengeti_pt)

# High Arctic
HighArctic <- st_read("data/raw/polygons/CPCAD-BDCAPC_Dec2021.gdb/", layer = "CPCAD_BDCAPC_Dec2021") %>%
  filter(NAME_E == "Bylot Island Bird Sanctuary", BIOME == "T")
HighArctic_pt <- do.call(rbind, st_sample(HighArctic, 100)) %>%
  as_tibble() %>% setNames(c("lon","lat")) %>%
  SpatialPoints(proj4string = CRS("+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")) %>%
  spTransform(CRS("+init=epsg:3857"))

HighArctic_clims <- raster::extract(r, HighArctic_pt)

# Environmental distance
env.pca <- prcomp(rbind(Europe_clims, Pyrenees_clims, Serengeti_clims, HighArctic_clims),
                  center = T, scale.= T)
summary(env.pca)
fws <- c(rep("Europe", 100), rep("Pyrenees", 100), rep("Serengeti", 100), rep("High Arctic", 100))
ggbiplot(env.pca, ellipse = T, groups = fws, obs.scale = 1, var.scale = 1) +
  scale_colour_manual(name = "Food web", values = c("royalblue4", "deepskyblue", "red2", "yellowgreen")) +
  theme_minimal() +
  theme(legend.position = "bottom")

env.centroid <- as.data.frame(env.pca$x) %>% cbind(fws) %>%
  group_by(fws) %>%
  summarise_all(mean)

env.dist <- dist(env.centroid[,c(2:7)])

# Geographic distance
coords <- rbind(Europe_pt@coords, Pyrenees_pt@coords, Serengeti_pt@coords, HighArctic_pt@coords)
geo.centroid <- as.data.frame(coords) %>% cbind(fws) %>%
  group_by(fws) %>%
  summarise_all(mean)

geo.dist <- dist(geo.centroid[,c(2,3)]/1000)

# Compositional distance
JaccardDissimilarity <- function(x, y){
  a = sum(x %in% y)
  b = sum(!(x %in% y))
  c = sum(!(y %in% x))
  return(1 - (a / (a+b+c)))
}

Europe.species <- read.csv("data/cleaned/EuroMWTaxo.csv")$Species
Pyrenees.species <- read.csv("data/cleaned/pyrenneesFWTaxo.csv")$Species
Serengeti.species <- read.csv("data/cleaned/SerengetiFWTaxo.csv")$Species
HighArctic.species <- read.csv("data/cleaned/HighArcticFWTaxo.csv")$Species

comp.dist <- c(JaccardDissimilarity(Europe.species, HighArctic.species),
               JaccardDissimilarity(Europe.species, Pyrenees.species),
               JaccardDissimilarity(Europe.species, Serengeti.species),
               JaccardDissimilarity(HighArctic.species, Pyrenees.species),
               JaccardDissimilarity(HighArctic.species, Serengeti.species),
               JaccardDissimilarity(Pyrenees.species, Serengeti.species))

# Everything together:
FWdist <- data.frame(FWs = c("Europe - High Arctic", "Europe - Pyrenees", "Europe - Serengeti",
                             "High Arctic - Pyrenees", "High Arctic - Serengeti",
                             "Pyrenees - Serengeti"),
                     Geo.Dist = c(geo.dist),
                     Env.Dist = c(env.dist),
                     Comp.Dist = comp.dist)
