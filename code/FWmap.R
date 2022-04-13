# Map of the Food Webs ----------------------------------------------------
library(dplyr)
library(sf)
library(tidyr)
library(ggrepel)
library(scatterpie)
library(maps)
sf_use_s2(FALSE)

Europe <- st_read("data/raw/polygons/BiogeoRegions2016_shapefile/BiogeoRegions2016.shp") %>%
  filter(short_name != "outside")
Pyrenees <- st_read("data/raw/polygons/EuropeanMountainAreas/m_massifs_v1.shp") %>%
  filter(name_mm == "Pyrenees")
Serengeti <- st_read("data/raw/polygons/Serengeti_Ecosystem/v3_serengeti_ecosystem.shp") %>% st_union()
HighArctic <- st_read("data/raw/polygons/HighArctic/Bylot.shp")

# Polygons
world <- map_data('world')
Europe <- st_transform(Europe, st_crs("EPSG:4326")) %>% st_union()
Europe_centroid <- st_coordinates(st_centroid(Europe))
Pyrenees <- st_transform(Pyrenees, st_crs("EPSG:4326")) %>% st_union() 
Pyrenees_centroid <- st_coordinates(st_centroid(Pyrenees))
HighArctic <- st_transform(HighArctic, st_crs("EPSG:4326")) %>% st_union()
HighArctic_centroid <- cbind(st_coordinates(st_centroid(HighArctic)))
Serengeti <- st_transform(Serengeti, st_crs("EPSG:4326"))
Serengeti_centroid <- st_coordinates(st_centroid(Serengeti))

# Curves
df <- expand_grid(from = c("Europe", "Pyrenees", "HighArctic", "Serengeti"), to = c("Europe", "Pyrenees", "HighArctic", "Serengeti")) %>%
  filter(from != to)
df[which(df$from == "Europe"), c("fromlong", "fromlat")] <- Europe_centroid
df[which(df$from == "Pyrenees"), c("fromlong", "fromlat")] <- Pyrenees_centroid
df[which(df$from == "HighArctic"), c("fromlong", "fromlat")] <- HighArctic_centroid
df[which(df$from == "Serengeti"), c("fromlong", "fromlat")] <- Serengeti_centroid
df[which(df$to == "Europe"), c("tolong", "tolat")] <- Europe_centroid
df[which(df$to == "Pyrenees"), c("tolong", "tolat")] <- Pyrenees_centroid
df[which(df$to == "HighArctic"), c("tolong", "tolat")] <- HighArctic_centroid
df[which(df$to == "Serengeti"), c("tolong", "tolat")]  <- Serengeti_centroid
df$col <- NA
df$col[which(df$from == "Europe")] <- "royalblue4"
df$col[which(df$from == "Pyrenees")] <- "red3"
df$col[which(df$from == "HighArctic")] <- "deepskyblue"
df$col[which(df$from == "Serengeti")] <- "yellowgreen"

# Species proportions
Europe_sp <- read.csv("data/cleaned/EuroMWTaxo.csv", row.names = 1)
Pyrenees_sp <- read.csv("data/cleaned/pyrenneesFWTaxo.csv", row.names = 1)
HighArctic_sp <- read.csv("data/cleaned/HighArcticFWTaxo.csv", row.names = 1)
Serengeti_sp <- read.csv("data/cleaned/SerengetiFWTaxo.csv", row.names = 1)
comp <- bind_rows(c(table(Europe_sp$Class) / nrow(Europe_sp), size = nrow(Europe_sp)),
          c(table(Pyrenees_sp$Class) / nrow(Pyrenees_sp), size = nrow(Pyrenees_sp)),
          c(table(HighArctic_sp$Class) / nrow(HighArctic_sp), size = nrow(HighArctic_sp)),
          c(table(Serengeti_sp$Class) / nrow(Serengeti_sp), size = nrow(Serengeti_sp))) %>%
  mutate(FW = c("Europe", "Pyrenees", "High Arctic", "Serengeti")) %>%
  pivot_longer(cols = -c(FW, size), values_to = "Proportion", names_to = "Class")
comp[is.na(comp)] <- 0
comp[which(comp$FW == "Europe"), c("long", "lat")] <- Europe_centroid
comp[which(comp$FW == "Pyrenees"), c("long", "lat")] <- Pyrenees_centroid
comp[which(comp$FW == "High Arctic"), c("long", "lat")] <- HighArctic_centroid
comp[which(comp$FW == "Serengeti"), c("long", "lat")] <- Serengeti_centroid
euro_barchart <- ggplotGrob(
  ggplot(comp[comp$FW == "Europe",])+
    geom_bar(aes(x=Class, fill = Class, y = Proportion),
             position='dodge',stat='identity', color = "black") +
    labs(x = NULL, y = NULL) + 
    theme(legend.position = "none", rect = element_blank(),
          line = element_blank(), text = element_blank()) 
)
pyrenees_barchart <- ggplotGrob(
  ggplot(comp[comp$FW == "Pyrenees",])+
    geom_bar(aes(x=Class, fill = Class, y = Proportion),
             position='dodge',stat='identity', color = "black") +
    labs(x = NULL, y = NULL) + 
    theme(legend.position = "none", rect = element_blank(),
          line = element_blank(), text = element_blank()) 
)
arctic_barchart <- ggplotGrob(
  ggplot(comp[comp$FW == "High Arctic",])+
    geom_bar(aes(x=Class, fill = Class, y = Proportion),
             position='dodge',stat='identity', color = "black") +
    labs(x = NULL, y = NULL) + 
    theme(legend.position = "none", rect = element_blank(),
          line = element_blank(), text = element_blank()) 
)
serengeti_barchart <- ggplotGrob(
  ggplot(comp[comp$FW == "Serengeti",])+
    geom_bar(aes(x=Class, fill = Class, y = Proportion),
             position='dodge',stat='identity', color = "black") +
    labs(x = NULL, y = NULL) + 
    theme(legend.position = "none", rect = element_blank(),
          line = element_blank(), text = element_blank()) 
)

p <- ggplot() +
  geom_polygon(data =world, aes(long, lat, group = group), color = "grey90", fill = "grey90") +
  geom_sf(data = Europe, fill = alpha("royalblue4", 0.5), color = "royalblue4") +
  geom_sf(data = Pyrenees, fill = alpha("red3", 0.5), color = "red3") +
  geom_sf(data = HighArctic, fill = alpha("deepskyblue", 0.5), color = "deepskyblue") +
  geom_sf(data = Serengeti, fill = alpha("yellowgreen", 0.5), color =  "yellowgreen") +
  geom_curve(data=df,
             aes(x=fromlong, y=fromlat, xend=tolong, yend=tolat, color = col),
             size=0.75,
             curvature=0.25, alpha = 0.6) + 
  scale_colour_identity() +
  geom_label_repel(data = as.data.frame(Pyrenees_centroid), aes(x = X, y = Y, label = "Pyrenees"), 
                   fontface = "bold", nudge_x = -15, nudge_y = -10, colour = "red3") +
  geom_label_repel(data = as.data.frame(Europe_centroid), aes(x = X, y = Y, label = "Europe"), 
                   fontface = "bold", nudge_x = 25, nudge_y = 20, colour = "royalblue4") +
  geom_label_repel(data = as.data.frame(HighArctic_centroid), aes(x = X, y = Y, label = "High Arctic"), 
                   fontface = "bold", nudge_x = -5, nudge_y = -10, colour = "deepskyblue") +
  geom_label_repel(data = as.data.frame(Serengeti_centroid), aes(x = X, y = Y, label = "Serengeti"), 
                   fontface = "bold", nudge_x = 18, nudge_y = 13, colour = "yellowgreen") +
  annotation_custom(euro_barchart, xmin = Europe_centroid[,1]+18, xmax = Europe_centroid[,1]+32, 
                    ymin = Europe_centroid[,2]+6, ymax = Europe_centroid[,2]+18) + 
  annotation_custom(pyrenees_barchart, xmin = Pyrenees_centroid[,1]-22, xmax = Pyrenees_centroid[,1]-8, 
                    ymin = Pyrenees_centroid[,2]-24, ymax = Pyrenees_centroid[,2]-12) + 
  annotation_custom(arctic_barchart, xmin = HighArctic_centroid[,1]-12, xmax = HighArctic_centroid[,1]+2, 
                    ymin = HighArctic_centroid[,2]-24, ymax = HighArctic_centroid[,2]-12) + 
  annotation_custom(serengeti_barchart, xmin = Serengeti_centroid[,1]+11, xmax = Serengeti_centroid[,1]+25, 
                    ymin = Serengeti_centroid[,2]-3, ymax = Serengeti_centroid[,2]+11) + 
  coord_sf(xlim = c(-85, 70), ylim = c(0, 80)) +
  theme_void()

ggsave("figures/SI/FWmap.pdf", p)
