library(dplyr)
library(sf)
library(raster)
library(ggbiplot)

N = 500 # number of points sampled per polygon
r <- getData("worldclim",var="bio",res=10, path = "data/raw/worldclim")
names(r) <- c("Tmean", "DiuRange", "Isothermality", "Seasonality",
              "Tmax", "Tmin", "Trange", "Twetmonth", "Tdrymonth", 
              "Twarmmonth", "Tcoldmonth", "PrecMean", "PrecWetmonth", "PrecDrymonth",
              "PrecSeasonality", "PrecWetquat", "PrecDryquat", "PrecWarmquat",
              "PrecColdquat")

# Europe
Europe <- st_read("data/raw/polygons/BiogeoRegions2016_shapefile/BiogeoRegions2016.shp") %>%
  filter(short_name != "outside")
Europe_pt <- do.call(rbind, st_sample(Europe, N)) %>%
  as_tibble() %>% setNames(c("lon","lat")) %>%
  SpatialPoints(proj4string = CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs")) %>%
  spTransform(CRS("+init=epsg:3857"))

Europe_clims <- raster::extract(r, Europe_pt)

# Pyrenees
Pyrenees <- st_read("data/raw/polygons/EuropeanMountainAreas/m_massifs_v1.shp") %>%
  filter(name_mm == "Pyrenees")
Pyrenees_pt <- do.call(rbind, st_sample(Pyrenees, N)) %>%
  as_tibble() %>% setNames(c("lon","lat")) %>%
  SpatialPoints(proj4string = CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs")) %>%
  spTransform(CRS("+init=epsg:3857"))

Pyrenees_clims <- raster::extract(r, Pyrenees_pt)

# Serengeti
Serengeti <- st_read("data/raw/polygons/Serengeti_Ecosystem/v3_serengeti_ecosystem.shp") %>% st_union()
Serengeti_pt <- do.call(rbind, st_sample(Serengeti, N)) %>%
  as_tibble() %>% setNames(c("lon","lat")) %>%
  SpatialPoints(proj4string = CRS("+proj=utm +zone=36 +south +ellps=clrk80 +units=m +no_defs")) %>%
  spTransform(CRS("+init=epsg:3857"))

Serengeti_clims <- raster::extract(r, Serengeti_pt)

# High Arctic
#HighArctic <- st_read("data/raw/polygons/CPCAD-BDCAPC_Dec2021.gdb/", layer = "CPCAD_BDCAPC_Dec2021") %>%
  #filter(NAME_E == "Bylot Island Bird Sanctuary", BIOME == "T")
#st_write(HighArctic, "data/raw/polygons/HighArctic/Bylot.shp")
HighArctic <- st_read("data/raw/polygons/HighArctic/Bylot.shp")
HighArctic_pt <- do.call(rbind, st_sample(HighArctic, N)) %>%
  as_tibble() %>% setNames(c("lon","lat")) %>%
  SpatialPoints(proj4string = CRS("+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")) %>%
  spTransform(CRS("+init=epsg:3857"))

HighArctic_clims <- raster::extract(r, HighArctic_pt)

# Environmental distance
env.pca <- prcomp(rbind(Europe_clims, Pyrenees_clims, Serengeti_clims, HighArctic_clims),
                  center = T, scale.= T)
summary(env.pca)
fws <- c(rep("Europe", N), rep("Pyrenees", N), rep("Serengeti", N), rep("High Arctic", N))
p <- ggbiplot(env.pca, ellipse = T, groups = fws, obs.scale = 1, var.scale = 1) +
  scale_colour_manual(name = "Food web", values = c("royalblue4", "deepskyblue", "red2", "yellowgreen")) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("figures/SI/PCAenv.pdf", p)

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

# Functional distance (mean pairwise distance)
library(mFD)

traits <- read.csv("data/cleaned/SpeciesTraitsFull.csv", row.names = 1, stringsAsFactors = T)
rownames(traits) <- traits$Species
traits <- traits %>%
  dplyr::select(-Species, -Class, -Order, -Family, -Genus)

traits_cat <- data.frame(trait_name = colnames(traits),
                         trait_type = c("N", "Q", rep("F", 12), "Q", "Q", "Q", "F", "F", "F"),
                         fuzzy_name = c(NA, NA, rep("Habitat", 12), NA, NA, NA, rep("TrophicLevel", 3)))
community <- data.frame(rbind(Europe = as.numeric(rownames(traits) %in% Europe.species),
                   Pyrenees = as.numeric(rownames(traits) %in% Pyrenees.species),
                   Serengeti = as.numeric(rownames(traits) %in% Serengeti.species),
                   HighArctic = as.numeric(rownames(traits) %in% HighArctic.species)
                   ))
colnames(community) <- rownames(traits)

funct_dist <- funct.dist(traits, traits_cat, "gower", scale_euclid = "scale_center")
fspace <- quality.fspaces(funct_dist, maxdim_pcoa = 10)
sp_faxes_coord <- fspace$details_fspaces$sp_pc_coord
quality.fspaces.plot(fspace, quality_metric = "mad", 
                     fspaces_plot = c("pcoa_2d", "pcoa_3d", "pcoa_4d", "pcoa_5d", "pcoa_6d"),
                     gradient_deviation = c(neg = "darkblue", nul = "grey80", pos = "darkred"),
                     gradient_deviation_quality = c(low = "yellow", high = "red"))


fw_func_beta <- beta.fd.hill(as.matrix(community), funct_dist)

# Everything together:
FWdist <- data.frame(FW1 = c("Europe", "Europe", "Europe",
                             "High Arctic", "High Arctic",
                             "Pyrenees", "Europe", "Serengeti"),
                     FW2 = c("High Arctic", "Pyrenees", "Serengeti", "Pyrenees", "Serengeti", "Serengeti", "Europe", "Serengeti"),
                     Geo.Dist = c(geo.dist,NA,NA),
                     Env.Dist = c(env.dist,NA,NA),
                     Comp.Dist = c(comp.dist,NA,NA))

# Visualie
p1 <- ggplot(FWdist, aes(x = FW1, y = FW2)) +
  geom_tile(aes(fill = Geo.Dist/1000), color = "white") +
  scale_fill_gradient2(low = "white", high = "red", name="Geographic\nDistance\n(10^3km)", na.value="transparent") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   hjust = 1))+
  coord_fixed() +
  theme(axis.title = element_blank(),
        legend.justification = c(1, 0),
        legend.position = "bottom",
        legend.direction = "horizontal")

p2 <- ggplot(FWdist, aes(x = FW1, y = FW2)) +
  geom_tile(aes(fill = Env.Dist), color = "white") +
  scale_fill_gradient2(low = "white", high = "red", name="Environmental\nDissimilarity", na.value="transparent") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                   hjust = 1))+
  coord_fixed() +
  theme(axis.title = element_blank(),
        legend.justification = c(1, 0),
        legend.position = "bottom",
        legend.direction = "horizontal")

p3 <- ggplot(FWdist, aes(x = FW1, y = FW2)) +
  geom_tile(aes(fill = Comp.Dist), color = "white") +
  scale_fill_gradient(low = "white", high = "red", name="Compositional\nDissimilarity", na.value="transparent",
                       limit = c(0.8,1)) +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                   hjust = 1))+
  coord_fixed() +
  theme(axis.title = element_blank(),
        legend.justification = c(1, 0),
        legend.position = "bottom",
        legend.direction = "horizontal")

library(patchwork)
p<-p1+p2+p3
ggsave("figures/SI/FWdist.png", p, device = "png")

# Map
library(rnaturalearth)
library(rnaturalearthdata)
library(ggrepel)
world <- ne_countries(scale = "medium", returnclass = "sf")
sf_use_s2(FALSE)
Europe <- st_transform(Europe, crs(world)) %>% st_union()
Europe_centroid <- st_coordinates(st_centroid(Europe))
Pyrenees <- st_transform(Pyrenees, crs(world)) %>% st_union() 
Pyrenees_centroid <- st_coordinates(st_centroid(Pyrenees))
HighArctic <- st_transform(HighArctic, crs(world)) %>% st_union()
HighArctic_centroid <- cbind(st_coordinates(st_centroid(HighArctic)))
Serengeti <- st_transform(Serengeti, crs(world))
Serengeti_centroid <- st_coordinates(st_centroid(Serengeti))

p <- ggplot(size = 0.1) +
  geom_sf(data = world) +
  geom_sf(data = Europe, fill = "royalblue4", color = "transparent", alpha = 0.8) +
  geom_sf(data = Pyrenees, fill = "red3", alpha = 0.8) +
  geom_sf(data = HighArctic, fill = "deepskyblue", alpha = 0.8) +
  geom_sf(data = Serengeti, fill = "yellowgreen", alpha = 0.8) +
  geom_label_repel(data = as.data.frame(Pyrenees_centroid), aes(x = X, y = Y, label = "Pyrenees"), 
                  fontface = "bold", nudge_x = -15, nudge_y = 5, colour = "red3") +
  geom_label_repel(data = as.data.frame(Europe_centroid), aes(x = X, y = Y, label = "Europe"), 
                  fontface = "bold", nudge_x = -25, nudge_y = 15, colour = "royalblue4") +
  geom_label_repel(data = as.data.frame(HighArctic_centroid), aes(x = X, y = Y, label = "High Arctic"), 
                  fontface = "bold", nudge_x = 15, nudge_y = -5, colour = "deepskyblue") +
  geom_label_repel(data = as.data.frame(Serengeti_centroid), aes(x = X, y = Y, label = "Serengeti"), 
                  fontface = "bold", nudge_x = 10, nudge_y = 5, colour = "yellowgreen") +
  coord_sf(xlim = c(-85, 70), ylim = c(0, 80)) +
  theme_bw()

ggsave("figures/SI/FWmap.pdf", p)
