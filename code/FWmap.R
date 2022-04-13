# Map of the Food Webs ----------------------------------------------------
library(dplyr)
library(sf)
library(tidyr)
library(ggrepel)
library(maps)
sf_use_s2(FALSE)

Europe <- st_read("data/raw/polygons/BiogeoRegions2016_shapefile/BiogeoRegions2016.shp") %>%
  filter(short_name != "outside")
Pyrenees <- st_read("data/raw/polygons/EuropeanMountainAreas/m_massifs_v1.shp") %>%
  filter(name_mm == "Pyrenees")
Serengeti <- st_read("data/raw/polygons/Serengeti_Ecosystem/v3_serengeti_ecosystem.shp") %>% st_union()
HighArctic <- st_read("data/raw/polygons/HighArctic/Bylot.shp")

world <- map_data('world')
Europe <- st_transform(Europe, st_crs("EPSG:4326")) %>% st_union()
Europe_centroid <- st_coordinates(st_centroid(Europe))
Pyrenees <- st_transform(Pyrenees, st_crs("EPSG:4326")) %>% st_union() 
Pyrenees_centroid <- st_coordinates(st_centroid(Pyrenees))
HighArctic <- st_transform(HighArctic, st_crs("EPSG:4326")) %>% st_union()
HighArctic_centroid <- cbind(st_coordinates(st_centroid(HighArctic)))
Serengeti <- st_transform(Serengeti, st_crs("EPSG:4326"))
Serengeti_centroid <- st_coordinates(st_centroid(Serengeti))

df <- expand_grid(from = c("Europe", "Pyrenees", "HighArctic", "Serengeti"), to = c("Europe", "Pyrenees", "HighArctic", "Serengeti")) %>%
  filter(from != to)
df[,c("tolong", "tolat")] <- NA
df$col <- NA
df[which(df$from == "Europe"), c("fromlong", "fromlat")] <- Europe_centroid
df[which(df$from == "Pyrenees"), c("fromlong", "fromlat")] <- Pyrenees_centroid
df[which(df$from == "HighArctic"), c("fromlong", "fromlat")] <- HighArctic_centroid
df[which(df$from == "Serengeti"), c("fromlong", "fromlat")] <- Serengeti_centroid
df[which(df$to == "Europe"), c("tolong", "tolat")] <- Europe_centroid
df[which(df$to == "Pyrenees"), c("tolong", "tolat")] <- Pyrenees_centroid
df[which(df$to == "HighArctic"), c("tolong", "tolat")] <- HighArctic_centroid
df[which(df$to == "Serengeti"), c("tolong", "tolat")]  <- Serengeti_centroid
df$col[which(df$from == "Europe")] <- "royalblue4"
df$col[which(df$from == "Pyrenees")] <- "red3"
df$col[which(df$from == "HighArctic")] <- "deepskyblue"
df$col[which(df$from == "Serengeti")] <- "yellowgreen"

ggplot() +
  geom_polygon(data =world, aes(long, lat, group = group), color = "grey90", fill = "grey90") +
  geom_sf(data = Europe, fill = alpha("royalblue4", 0.5), color = "royalblue4") +
  geom_sf(data = Pyrenees, fill = alpha("red3", 0.5), color = "red3") +
  geom_sf(data = HighArctic, fill = alpha("deepskyblue", 0.5), color = "deepskyblue") +
  geom_sf(data = Serengeti, fill = alpha("yellowgreen", 0.5), color =  "yellowgreen") +
  geom_curve(data=df,
             aes(x=fromlong, y=fromlat, xend=tolong, yend=tolat, color = col),
             size=0.75,
             curvature=0.3, alpha = 0.6) + 
  scale_colour_identity() +
  geom_label_repel(data = as.data.frame(Pyrenees_centroid), aes(x = X, y = Y, label = "Pyrenees"), 
                   fontface = "bold", nudge_x = -15, nudge_y = 5, colour = "red3") +
  geom_label_repel(data = as.data.frame(Europe_centroid), aes(x = X, y = Y, label = "Europe"), 
                   fontface = "bold", nudge_x = -25, nudge_y = 15, colour = "royalblue4") +
  geom_label_repel(data = as.data.frame(HighArctic_centroid), aes(x = X, y = Y, label = "High Arctic"), 
                   fontface = "bold", nudge_x = 15, nudge_y = -5, colour = "deepskyblue") +
  geom_label_repel(data = as.data.frame(Serengeti_centroid), aes(x = X, y = Y, label = "Serengeti"), 
                   fontface = "bold", nudge_x = 10, nudge_y = 5, colour = "yellowgreen") +
  coord_sf(xlim = c(-85, 70), ylim = c(0, 80)) +
  theme_void()

ggsave("figures/SI/FWmap.pdf", p)
