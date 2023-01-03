library(dplyr)
library(tidyr)
library(ggplot2)

# Appendix other performance metrics --------------------------------------
foodwebs <- c("Euro", "Pyrenees", "Arctic", "Serengeti")
foodwebs_labs <- c("Europe", "Pyrenees", "NorthQC", "Serengeti")
overall_performance <- read.csv("data/checkpoints/overall_performance_draws.csv", row.names = 1) %>%
  pivot_longer(cols = tpr:npv, names_to = "metric") %>%
  mutate(Source = factor(Source, levels = foodwebs, labels = foodwebs_labs), Target = factor(Target, levels = foodwebs, labels = foodwebs_labs), metric = factor(metric, levels = c("tpr", "tnr", "ppv", "npv")))

ggplot(overall_performance) +
  geom_density(aes(value, after_stat(scaled) ,fill = Source), alpha = .5) +
  facet_grid(Target~metric, scales = "free") +
  scale_fill_manual(values = c("royalblue4", "red3", "deepskyblue", "chartreuse4")) +
  labs(y = "Scaled density") +
  theme_classic()

ggsave("figures/SI/performance_metrics.png", width = 9, dpi = 600)

# how correlated is auc to aucpr and other metrics
library(corrplot)
overall_performance_othermetrics <- read.csv("data/checkpoints/overall_performance_draws.csv", row.names = 1) %>%
  group_by(Source, Target) %>%
  summarise_all(median)

overall_performance_auc <- read.csv("data/checkpoints/overall_performance.csv", row.names=1)

overall_performance <- full_join(overall_performance_othermetrics, overall_performance_auc)

Movr <- cor(overall_performance[,c("auc", "aucpr", "tpr", "tnr", "ppv", "npv")])

png("figures/SI/metricCorr_overall.png", width = 300, height = 300)
corrplot(Movr, method = "number", type = "upper", diag = F)
dev.off()

species_performance_othermetrics <- read.csv("data/checkpoints/species_performance_draws.csv", row.names = 1) %>%
  group_by(Source, Target, species) %>%
  summarise_all(median, na.rm = T)

species_performance_auc <- read.csv("data/checkpoints/species_performance.csv", row.names=1)

species_performance <- full_join(species_performance_othermetrics, species_performance_auc)
Msp <- cor(species_performance[,c("auc", "aucpr", "tpr", "tnr", "ppv", "npv")])

png("figures/SI/metricCorr_species.png", width = 300, height = 300)
corrplot(Msp, method = "number", type = "upper", diag = F)
dev.off()

# Sensitivity for sampling 500 pts
library(sf)
library(raster)
library(ggbiplot)
library(patchwork)
library(purrr)

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
Europe <- st_read("data/raw/polygons/BiogeoRegions2016_shapefile/BiogeoRegions2016.shp") %>%
  dplyr::filter(short_name != "outside")
Pyrenees <- st_read("data/raw/polygons/EuropeanMountainAreas/m_massifs_v1.shp") %>%
  dplyr::filter(name_mm == "Pyrenees")
Serengeti <- st_read("data/raw/polygons/Serengeti_Ecosystem/v3_serengeti_ecosystem.shp") %>% st_union()
HighArctic <- st_read("data/raw/polygons/HighArctic/QuebecLabrador50.shp")

fws <- c(rep("Europe", N), rep("Pyrenees", N), rep("Serengeti", N), rep("Nunavik", N))

envdist <- as.data.frame(caTools::combs(c("Europe", "Pyrenees", "Serengeti", "Nunavik"), k =2))

for (i in c(1:10)){
  # europe
  Europe_pt <- do.call(rbind, st_sample(Europe, N)) %>%
    as_tibble() %>% setNames(c("lon","lat")) %>%
    SpatialPoints(proj4string = CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs")) %>%
    spTransform(CRS("+init=epsg:3857"))
  Europe_clims <- raster::extract(r, Europe_pt)
  
  # pyrenees
  Pyrenees_pt <- do.call(rbind, st_sample(Pyrenees, N)) %>%
    as_tibble() %>% setNames(c("lon","lat")) %>%
    SpatialPoints(proj4string = CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs")) %>%
    spTransform(CRS("+init=epsg:3857"))
  Pyrenees_clims <- raster::extract(r, Pyrenees_pt)
  
  # serengeti
  Serengeti_pt <- do.call(rbind, st_sample(Serengeti, N)) %>%
    as_tibble() %>% setNames(c("lon","lat")) %>%
    SpatialPoints(proj4string = CRS("+proj=utm +zone=36 +south +ellps=clrk80 +units=m +no_defs")) %>%
    spTransform(CRS("+init=epsg:3857"))
  Serengeti_clims <- raster::extract(r, Serengeti_pt)
  
  # high arctic
  HighArctic_pt <- do.call(rbind, st_sample(HighArctic, N)) %>%
    as_tibble() %>% setNames(c("lon","lat")) %>%
    SpatialPoints(proj4string = CRS(" +proj=lcc +lat_1=49 +lat_2=77 +lat_0=63.390675 +lon_0=-91.86666666666666 +x_0=6200000 +y_0=3000000 +datum=NAD83 +units=m +no_defs")) %>%
    spTransform(CRS("+init=epsg:3857"))
  HighArctic_clims <- raster::extract(r, HighArctic_pt)
  
  # calculate environmental distances
  env.centroid <- data.frame(rbind(Europe_clims, Pyrenees_clims, Serengeti_clims, HighArctic_clims)) %>% 
    scale() %>%
    as.data.frame() %>%
    cbind(fws) %>%
    group_by(fws) %>%
    summarise_all(mean)
  env.dist <- dist(env.centroid[,-1])
  
  envdist <- cbind(envdist, c(env.dist))
}

plot_df <- data.frame(pair = rep(paste(envdist$V1, envdist$V2, sep = "-")),
                      meandist = apply(envdist[,c(3:12)], MARGIN = 1, mean, na.rm = T),
                      mindist = apply(envdist[,c(3:12)], MARGIN = 1, min, na.rm = T),
                      maxdist = apply(envdist[,c(3:12)], MARGIN = 1, max, na.rm = T))

ggplot(plot_df) +
  geom_pointrange(aes(y = pair, x = meandist, xmin = mindist, xmax = maxdist)) +
  labs(x = "Environmental distance", y = "Food web pair") +
  theme_classic()
ggsave("figures/SI/environmentaldistance_sensitivity.png", width = 6)
