# Step xx: Measure dissimilarity between food webs
# 1. Calculate environmental Dissimilarity
# 2. Calculate geographic Distance
# 3. Calculate compositional dissimilarity
# 4. Calculate phylogenetic dissimilarity
# 5. Make plots of distance and map of food web

library(dplyr)
library(sf)
library(raster)
library(ggbiplot)
library(patchwork)
library(purrr)
source("code/DistanceFunctions.R")

# Environmental Dissimilarity ---------------------------------------------
# 1.In each food web, get worldclim data for N random points in the polygon
# 2.Do a PCA with the environmental variables
# 3.Calculated the distance of centroids in the PCA
N = 500 # number of points sampled per polygon
r <- getData("worldclim",var="bio",res=10, path = "data/raw/worldclim")
names(r) <- c("Tmean", "DiuRange", "Isothermality", "Seasonality",
              "Tmax", "Tmin", "Trange", "Twetmonth", "Tdrymonth", 
              "Twarmmonth", "Tcoldmonth", "PrecMean", "PrecWetmonth", "PrecDrymonth",
              "PrecSeasonality", "PrecWetquat", "PrecDryquat", "PrecWarmquat",
              "PrecColdquat")

# europe
Europe <- st_read("data/raw/polygons/BiogeoRegions2016_shapefile/BiogeoRegions2016.shp") %>%
  filter(short_name != "outside")
Europe_pt <- do.call(rbind, st_sample(Europe, N)) %>%
  as_tibble() %>% setNames(c("lon","lat")) %>%
  SpatialPoints(proj4string = CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs")) %>%
  spTransform(CRS("+init=epsg:3857"))
Europe_clims <- raster::extract(r, Europe_pt)

# pyrenees
Pyrenees <- st_read("data/raw/polygons/EuropeanMountainAreas/m_massifs_v1.shp") %>%
  filter(name_mm == "Pyrenees")
Pyrenees_pt <- do.call(rbind, st_sample(Pyrenees, N)) %>%
  as_tibble() %>% setNames(c("lon","lat")) %>%
  SpatialPoints(proj4string = CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs")) %>%
  spTransform(CRS("+init=epsg:3857"))
Pyrenees_clims <- raster::extract(r, Pyrenees_pt)

# serengeti
Serengeti <- st_read("data/raw/polygons/Serengeti_Ecosystem/v3_serengeti_ecosystem.shp") %>% st_union()
Serengeti_pt <- do.call(rbind, st_sample(Serengeti, N)) %>%
  as_tibble() %>% setNames(c("lon","lat")) %>%
  SpatialPoints(proj4string = CRS("+proj=utm +zone=36 +south +ellps=clrk80 +units=m +no_defs")) %>%
  spTransform(CRS("+init=epsg:3857"))
Serengeti_clims <- raster::extract(r, Serengeti_pt)

# high arctic
#HighArctic <- st_read("data/raw/polygons/CPCAD-BDCAPC_Dec2021.gdb/", layer = "CPCAD_BDCAPC_Dec2021") %>%
  #filter(NAME_E == "Bylot Island Bird Sanctuary", BIOME == "T")
#st_write(HighArctic, "data/raw/polygons/HighArctic/Bylot.shp")
HighArctic <- st_read("data/raw/polygons/HighArctic/Bylot.shp")
HighArctic_pt <- do.call(rbind, st_sample(HighArctic, N)) %>%
  as_tibble() %>% setNames(c("lon","lat")) %>%
  SpatialPoints(proj4string = CRS("+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")) %>%
  spTransform(CRS("+init=epsg:3857"))
HighArctic_clims <- raster::extract(r, HighArctic_pt)

# pca of climating variables
env.pca <- prcomp(rbind(Europe_clims, Pyrenees_clims, Serengeti_clims, HighArctic_clims),
                  center = T, scale.= T)
summary(env.pca)
fws <- c(rep("Europe", N), rep("Pyrenees", N), rep("Serengeti", N), rep("High Arctic", N))
p <- ggbiplot(env.pca, ellipse = T, groups = fws, obs.scale = 1, var.scale = 1) +
  scale_colour_manual(name = "Food web", values = c("royalblue4", "deepskyblue", "red2", "yellowgreen")) +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave("figures/SI/PCAenv.pdf", p)

# distance of the climatic centroid of each food web
env.centroid <- as.data.frame(env.pca$x) %>% cbind(fws) %>%
  group_by(fws) %>%
  summarise_all(mean)
env.dist <- dist(env.centroid[,c(2:7)])

# Geographic Distance -----------------------------------------------------
# 1. Find the centroid of the N points extracted earlier
# 2. Take the distance in km between all centroids
coords <- rbind(Europe_pt@coords, Pyrenees_pt@coords, Serengeti_pt@coords, HighArctic_pt@coords)
geo.centroid <- as.data.frame(coords) %>% cbind(fws) %>%
  group_by(fws) %>%
  summarise_all(mean)
geo.dist <- dist(geo.centroid[,c(2,3)]/1000) # distance in km


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
phydist <- as.matrix(read.csv("data/checkpoints/phylodist.csv", row.names = 1))
colnames(phydist) <- rownames(phydist)
community <- data.frame(rbind(Europe = as.numeric(rownames(phydist) %in% Europe.species),
                              HighArctic = as.numeric(rownames(phydist) %in% HighArctic.species),
                              Pyrenees = as.numeric(rownames(phydist) %in% Pyrenees.species),
                              Serengeti = as.numeric(rownames(phydist) %in% Serengeti.species)))
colnames(community) <- rownames(phydist)
fw.comdistnt <- phylobetadiv(community, phydist)


# Visualization -----------------------------------------------------------
# put everything together
FWdist <- data.frame(FW1 = c("Europe", "Europe", "Europe",
                             "High Arctic", "High Arctic",
                             "Pyrenees", "Europe", "Serengeti"),
                     FW2 = c("High Arctic", "Pyrenees", "Serengeti", "Pyrenees", "Serengeti", "Serengeti", "Europe", "Serengeti"),
                     Geo.Dist = c(geo.dist,NA,NA),
                     Env.Dist = c(env.dist,NA,NA),
                     Comp.Dist = c(comp.dist,NA,NA))
write.csv(FWdist, "data/checkpoints/FWdist.csv")

# geographic distance plot
p1 <- ggplot(FWdist, aes(x = FW1, y = FW2)) +
  geom_tile(aes(fill = Geo.Dist/1000), color = "white") +
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
  geom_tile(aes(fill = Env.Dist), color = "white") +
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

# compositional dissimilarity plot
p3 <- ggplot(FWdist, aes(x = FW1, y = FW2)) +
  geom_tile(aes(fill = Comp.Dist), color = "white") +
  scale_fill_gradient(low = "white", high = "red", name="Compositional Dissimilarity", na.value="transparent",
                       limit = c(0.8,1),
                      guide = guide_colourbar(title.position = "top", 
                                              barwidth = 8,
                                              barheight = 0.75)) +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                   hjust = 1))+
  coord_fixed() +
  theme(axis.title = element_blank(),text = element_text(size =8),
        legend.justification = c(1, 0),
        legend.position = "bottom",
        legend.direction = "horizontal")

# phylogenetic dissimilarity plot
df <- reshape2::melt(fw.comdistnt, value.name = "Phylo.Dist")
p4 <- ggplot(df, aes(x = Var1, y = Var2)) +
  geom_tile(aes(fill = Phylo.Dist), color = "white") +
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
p<-p1+p2+p3 + p4 + plot_layout(nrow =2)
ggsave("figures/SI/FWdist.png", p, device = "png")
