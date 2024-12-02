# Step 04: Calculate dissimilarity between food webs
# 1. Calculate environmental Dissimilarity
# 2. Calculate geographic Distance
# 3. Calculate compositional dissimilarity
# 4. Calculate phylogenetic dissimilarity
# 5. Make plots of distance and map of food web

library(dplyr)
library(tidyr)
library(sf)
library(terra)
library(patchwork)
library(purrr)
library(geodata)
library(geos)
source("code/functions_distance.R")

# Environmental Dissimilarity ---------------------------------------------
# 1.In each food web, get worldclim data for N random points in the polygon
# 2.Do a PCA with the environmental variables
# 3.Calculated the distance of centroids in the PCA
N = 500 # number of points sampled per polygon
r <- worldclim_global(var="bio",res=10, path = "data/raw/worldclim")
names(r) <- c("Tmean", "DiuRange", "Isothermality", "Seasonality",
              "Tmax", "Tmin", "Trange", "Twetmonth", "Tdrymonth", 
              "Twarmmonth", "Tcoldmonth", "PrecMean", "PrecWetmonth", "PrecDrymonth",
              "PrecSeasonality", "PrecWetquat", "PrecDryquat", "PrecWarmquat",
              "PrecColdquat")

# europe
Europe y<- st_read("data/raw/polygons/BiogeoRegions2016_shapefile/BiogeoRegions2016.shp") %>%
  filter(short_name != "outside") %>% 
  st_union() ## st_union: all geometries are unioned together and an sfg or single-geometry sfc object is returned.
              # Unioning a set of overlapping polygons has the effect of merging the areas (i.e. the same effect as iteratively unioning all individual polygons together)


Europe_v <- vect(Europe) %>% project(crs(r))
Europe_clims <- terra::extract(r, Europe_v)

## with terra
library(tidyterra)
Europe x<- vect("data/raw/polygons/BiogeoRegions2016_shapefile/BiogeoRegions2016.shp") %>% 
  filter(short_name != "outside") %>%
  union() ## union: 


# pyrenees
Pyrenees <- st_read("data/raw/polygons/EuropeanMountainAreas/m_massifs_v1.shp") %>%
  filter(name_mm == "Pyrenees") %>% st_union()
Pyrenees_v <- vect(Pyrenees) %>% project(crs(r))
Pyrenees_clims <- terra::extract(r, Pyrenees_v)

# serengeti
Serengeti <- st_read("data/raw/polygons/Serengeti_Ecosystem/v3_serengeti_ecosystem.shp") %>% st_union()
Serengeti_v <- vect(Serengeti) %>% project(crs(r))
Serengeti_clims <- terra::extract(r, Serengeti_v)

# high arctic
HighArctic <- st_read("data/raw/polygons/HighArctic/QuebecLabrador50.shp") %>% st_union()
HighArctic_v <- vect(HighArctic) %>% project(crs(r))
HighArctic_clims <- terra::extract(r, HighArctic_v)

fws <- c(rep("Europe", nrow(Europe_clims)), rep("Pyrenees", nrow(Pyrenees_clims)), rep("Serengeti", nrow(Serengeti_clims)), rep("Nunavik", nrow(HighArctic_clims)))

# pca of climating variables
# env.pca <- prcomp(rbind(Europe_clims, Pyrenees_clims, Serengeti_clims, HighArctic_clims),
#                   center = T, scale.= T)
# summary(env.pca)
# p <- ggbiplot(env.pca, ellipse = T, groups = fws, obs.scale = 1, var.scale = 1) +
#   scale_colour_manual(name = "Food web", values = c("royalblue4", "deepskyblue", "red2", "yellowgreen")) +
#   theme_minimal() +
#   theme(legend.position = "bottom")
# ggsave("figures/SI/PCAenv.pdf", p)

# calculate environmental distances
env.centroid <- data.frame(rbind(Europe_clims, Pyrenees_clims, Serengeti_clims, HighArctic_clims)) %>% 
  scale() %>%
  as.data.frame() %>%
  cbind(fws) %>%
  group_by(fws) %>%
  summarise_all(mean, na.rm = T)
env.dist <- dist(env.centroid[,-c(1,2)])

# Geographic Distance -----------------------------------------------------
# 1. Get the centroid of each polygon
# 2. Transform to geodetic coordinates
# 3. Create a multipoint df object
# 4. Use the function st_distance to get the great-circle distance
sf_use_s2(FALSE)
europe_centroid <- Europe %>% st_transform(crs(r)) %>% st_centroid()
pyrenees_centroid <- Pyrenees %>% st_transform(crs(r)) %>% st_centroid()
serengeti_centroid <- Serengeti %>% st_transform(crs(r)) %>% st_centroid()
nunavik_centroid <- HighArctic %>% st_transform(crs(r)) %>% st_centroid()

geo.dist <- st_distance(c(europe_centroid, nunavik_centroid, pyrenees_centroid, serengeti_centroid))/1000

# Compositional Dissimilarity ---------------------------------------------
# jaccard dissimilarity
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

# Phylogenetic Dissimilarity ----------------------------------------------
# the mean phylogenetic distance to the nearest taxon
phydist <- as.matrix(read.csv("data/checkpoints/phylodist_mean.csv", row.names = 1))
# square root the phylogenetic distance matrix (Letten & Cornwell, 2014)
phydist <- sqrt(phydist)
colnames(phydist) <- rownames(phydist)
community <- data.frame(rbind(Europe = as.numeric(rownames(phydist) %in% Europe.species),
                              HighArctic = as.numeric(rownames(phydist) %in% HighArctic.species),
                              Pyrenees = as.numeric(rownames(phydist) %in% Pyrenees.species),
                              Serengeti = as.numeric(rownames(phydist) %in% Serengeti.species)))
colnames(community) <- rownames(phydist)
fw.comdistnt <- phylobetadiv(community, phydist)

# Visualization -----------------------------------------------------------
# put everything together
foodwebs <- c("Europe", "Nunavik", "Pyrenees", "Serengeti")
FWdist <- expand_grid(FW1 = foodwebs, FW2 = foodwebs)
FWdist$geo.dist <- c(as.matrix(geo.dist))
FWdist$env.dist <- c(as.matrix(env.dist))
FWdist$phylo.dist <- c(fw.comdistnt)
write.csv(FWdist, "data/checkpoints/FWdist.csv")

# geographic distance plot
p1 <- ggplot(FWdist, aes(x = FW1, y = FW2)) +
  geom_tile(aes(fill = geo.dist/1000), color = "white") +
  scale_fill_gradient2(low = "white", high = "red", name="Geographic Distance (10^3km)", na.value="transparent",
                       guide = guide_colourbar(title.position = "top", 
                                               barwidth = 8,
                                               barheight = 0.75)) +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   hjust = 1))+
  coord_fixed() +
  theme(axis.title = element_blank(),
        legend.justification = "center",
        legend.position = "top",
        legend.direction = "horizontal",
        text = element_text(size = 8))

# environmental dissimilarity plot
p2 <- ggplot(FWdist, aes(x = FW1, y = FW2)) +
  geom_tile(aes(fill = env.dist), color = "white") +
  scale_fill_gradient2(low = "white", high = "red", name="Environmental Dissimilarity", na.value="transparent",
                       guide = guide_colourbar(title.position = "top", 
                                               barwidth = 8,
                                               barheight = 0.75)) +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                   hjust = 1))+
  coord_fixed() +
  theme(axis.title = element_blank(),
        legend.justification = c(1, 0),
        legend.position = "top",
        legend.direction = "horizontal",
        text = element_text(size = 8))

# phylogenetic dissimilarity plot
p3 <- ggplot(FWdist, aes(x = FW1, y = FW2)) +
  geom_tile(aes(fill = phylo.dist), color = "white") +
  scale_fill_gradient2(low = "white", high = "red", name="Phylogenetic Dissimilarity", na.value="transparent",
                       guide = guide_colourbar(title.position = "top", 
                                               barwidth = 8,
                                               barheight = 0.75)) +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                   hjust = 1))+
  coord_fixed() +
  theme(axis.title = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal",
        text = element_text(size = 8))

# patchwork and save
p<-p1+p2+p3 +plot_layout(nrow =2)
ggsave("figures/SI/FWdist.png", p, device = "png")
