# Icons of each food web --------------------------------------------------
library(igraph)
library(dplyr)
library(NetIndices)
library(ggplot2)
set.seed(1234)
par(bg=NA)

# nunavik
nunavikFW <- read.csv("data/cleaned/HighArcticFW.csv", row.names = 1)
nunavikFW <- nunavikFW %>%
  transmute(resource = Prey, consumer = Predator)
nunavik.igraph <- simplify(graph_from_edgelist(as.matrix(nunavikFW)))
nunavik.matrix <- Matrix::as.matrix(as_adjacency_matrix(nunavik.igraph))

lay<-matrix(nrow = nrow(nunavik.matrix), ncol=2) # create a matrix with one column as runif, the other as trophic level
lay[,1]<-sample(c(1:nrow(nunavik.matrix)), nrow(nunavik.matrix))
lay[,2]<-TrophInd(nunavik.matrix)$TL-1 + runif(nrow(nunavik.matrix), min =0, max = 1)

png("figures/FWgraphs/nunavik.png")
par(bg=NA)
plot(nunavik.igraph, layout = lay, vertex.label=NA, vertex.size=10, edge.arrow.size=0, edge.width=4,
     vertex.color = "deepskyblue", edge.color = alpha("deepskyblue", 0.2))
dev.off()

# simplified
vid <- sample(c(1:length(V(nunavik.igraph))), 60)
nunavik.simplified <- induced_subgraph(nunavik.igraph, vid)
Isolated = which(degree(nunavik.simplified)==0)
lay <- lay[vid,]
if (length(Isolated !=0)){
  lay <- lay[-Isolated,]
}
nunavik.simplified <- delete.vertices(nunavik.simplified, Isolated)
png("figures/conceptual/simplifiedNunavik.png")
par(bg=NA)
plot(nunavik.simplified, layout = lay, vertex.label=NA, vertex.size=10, edge.arrow.size=0, edge.width=5,
     vertex.color = "deepskyblue", edge.color = alpha("deepskyblue", 0.2))
dev.off()

# mock prediction
png("figures/conceptual/predictedNunavik.png")
par(bg=NA)
plot(nunavik.simplified, layout = lay, vertex.label=NA, vertex.size=10, edge.arrow.size=0, edge.width=runif(length(E(nunavik.simplified)), min = 0, max = 10),
     vertex.color = "deepskyblue", edge.color = alpha("deepskyblue", 0.2))
dev.off()

# europe
EuropeMW <- read.csv("data/cleaned/EuroFWadults.csv", row.names = 1)
sp <- unique(c(EuropeMW$Predator, EuropeMW$Prey))
EuropeMW <- EuropeMW %>%
  transmute(resource = Prey, consumer = Predator) %>%
  filter(resource %in% sp, consumer %in% sp)

europe.igraph <- simplify(graph_from_edgelist(as.matrix(EuropeMW)))
europe.matrix <- Matrix::as.matrix(as_adjacency_matrix(europe.igraph))

lay<-matrix(nrow = nrow(europe.matrix), ncol=2) # create a matrix with one column as runif, the other as trophic level
lay[,1]<-sample(c(1:nrow(europe.matrix)), nrow(europe.matrix))
lay[,2]<-TrophInd(europe.matrix)$TL-1 + runif(nrow(europe.matrix), min =0, max = 1)

png("figures/FWgraphs/europe.png")
par(bg=NA)
plot(europe.igraph, layout = lay, vertex.label=NA, vertex.size=10, edge.arrow.size=0, edge.width=4,
     vertex.color = "royalblue4", edge.color = alpha("royalblue4", 0.2))
dev.off()

# simplified
vid <- sample(c(1:length(V(europe.igraph))), 75)
europe.simplified <- induced_subgraph(europe.igraph, vid)
Isolated = which(degree(europe.simplified)==0)
lay <- lay[vid,]
if (length(Isolated !=0)){
  lay <- lay[-Isolated,]
}
europe.simplified <- delete.vertices(europe.simplified, Isolated)
png("figures/conceptual/simplifiedEurope.png")
par(bg=NA)
plot(europe.simplified, layout = lay, vertex.label=NA, vertex.size=10, edge.arrow.size=0, edge.width=5,
     vertex.color = "royalblue4", edge.color = alpha("royalblue4", 0.2))
dev.off()

# mock prediction
png("figures/conceptual/predictedEurope.png")
par(bg=NA)
plot(europe.simplified, layout = lay, vertex.label=NA, vertex.size=10, edge.arrow.size=0, edge.width=runif(length(E(europe.simplified)), min = 0, max = 10),
     vertex.color = "royalblue4", edge.color = alpha("royalblue4", 0.2))
dev.off()

# serengeti
SerengetiFW <- read.csv("data/cleaned/SerengetiFW.csv", row.names = 1)
sp <- unique(c(SerengetiFW$Consumer_Species, SerengetiFW$Resource_Species))
SerengetiFW <- SerengetiFW %>%
  transmute(resource = Resource_Species, consumer = Consumer_Species) %>%
  filter(resource %in% sp, consumer %in% sp, !is.na(resource), !is.na(consumer))

serengeti.igraph <- simplify(graph_from_edgelist(as.matrix(SerengetiFW)))
serengeti.matrix <- Matrix::as.matrix(as_adjacency_matrix(serengeti.igraph))

lay<-matrix(nrow = nrow(serengeti.matrix), ncol=2) # create a matrix with one column as runif, the other as trophic level
lay[,1]<-sample(c(1:nrow(serengeti.matrix)), nrow(serengeti.matrix))
lay[,2]<-TrophInd(serengeti.matrix)$TL
lay[lay[,2]<1,2] <- 1
lay[,2] <- lay[,2]- 1 + runif(nrow(serengeti.matrix), min =0, max = 1)

png("figures/FWgraphs/serengeti.png")
par(bg=NA)
serengeti_graph <- plot(serengeti.igraph, layout = lay, vertex.label=NA, vertex.size=10, edge.arrow.size=0, edge.width=4,
                        vertex.color = "chartreuse4", edge.color = alpha("chartreuse4", 0.2))
dev.off()

# simplified
vid <- sample(c(1:length(V(serengeti.igraph))), 50)
serengeti.simplified <- induced_subgraph(serengeti.igraph, vid)
Isolated = which(degree(serengeti.simplified)==0)
lay <- lay[vid,]
if (length(Isolated !=0)){
  lay <- lay[-Isolated,]
}
serengeti.simplified <- delete.vertices(serengeti.simplified, Isolated)
png("figures/FWgraphs/simplifiedSerengeti.png")
par(bg=NA)
plot(serengeti.simplified, layout = lay, vertex.label=NA, vertex.size=10, edge.arrow.size=0, edge.width=5,
     vertex.color = "chartreuse4", edge.color = alpha("chartreuse4", 0.2))
dev.off()

# mock prediction
png("figures/conceptual/predictedSerengeti.png")
par(bg=NA)
plot(serengeti.simplified, layout = lay, vertex.label=NA, vertex.size=10, edge.arrow.size=0, edge.width=runif(length(E(serengeti.simplified)), min = 0, max = 10),
     vertex.color = "chartreuse4", edge.color = alpha("chartreuse4", 0.2))
dev.off()

# pyrenees
PyreneesFW <- read.csv("data/cleaned/pyrenneesFW.csv", row.names = 1)
sp <- unique(c(PyreneesFW$Predator, PyreneesFW$Prey))
PyreneesFW <- PyreneesFW %>%
  transmute(resource = Prey, consumer = Predator) %>%
  filter(resource %in% sp, consumer %in% sp)

pyrenees.igraph <- simplify(graph_from_edgelist(as.matrix(PyreneesFW)))
pyrenees.matrix <- Matrix::as.matrix(as_adjacency_matrix(pyrenees.igraph))

lay<-matrix(nrow = nrow(pyrenees.matrix), ncol=2) # create a matrix with one column as runif, the other as trophic level
lay[,1]<-sample(c(1:nrow(pyrenees.matrix)), nrow(pyrenees.matrix))
lay[,2]<-TrophInd(pyrenees.matrix)$TL - 1 + runif(nrow(pyrenees.matrix), min = 0, max = 1)

png("figures/FWgraphs/pyrenees.png")
par(bg=NA)
pyrenees_graph <- plot(pyrenees.igraph, layout = lay, vertex.label=NA, vertex.size=10, edge.arrow.size=0, edge.width=4,
                       vertex.color = "red3", edge.color = alpha("red3", 0.2))
dev.off()

# simplified
vid <- sample(c(1:length(V(pyrenees.igraph))), 40)
pyrenees.simplified <- induced_subgraph(pyrenees.igraph, vid)
Isolated = which(degree(pyrenees.simplified)==0)
lay <- lay[vid,]
if (length(Isolated !=0)){
  lay <- lay[-Isolated,]
}
pyrenees.simplified <- delete.vertices(pyrenees.simplified, Isolated)
png("figures/FWgraphs/simplifiedPyrenees.png")
par(bg=NA)
plot(pyrenees.simplified, layout = lay, vertex.label=NA, vertex.size=10, edge.arrow.size=0, edge.width=5,
     vertex.color = "red3", edge.color = alpha("red3", 0.2))
dev.off()

# mock prediction
png("figures/conceptual/predictedPyrenees.png")
par(bg=NA)
plot(pyrenees.simplified, layout = lay, vertex.label=NA, vertex.size=10, edge.arrow.size=0, edge.width=runif(length(E(pyrenees.simplified)), min = 0, max = 10),
     vertex.color = "red3", edge.color = alpha("red3", 0.2))
dev.off()

# Logistic regression icon ------------------------------------------------
library(tidybayes)
EuropeModel <- readRDS("~/OneDrive/Chapt-NicheSpace/models/EuroModel_brms.rds")
fe_only <- tibble(BM.prey = seq(-2, 2, length.out=100))
fe_only[,c("Omnivore.predator", "Carnivore.predator", "Habitat_breadth.predator", "Longevity.predator", "ClutchSize.predator", "Omnivore.prey", "Carnivore.prey",
           "Habitat_breadth.prey", "BM.predator", "Longevity.prey", "ClutchSize.prey", "ActivityTime.match",
           "Habitat.match", "BM.match")] <- 0
fe_only <- fe_only %>%
  add_epred_draws(EuropeModel,
                   re_formula = NA,
                   scale = "response", ndraws = 1e3)

fe_only_mean <- fe_only %>% 
  group_by(BM.prey) %>%
  summarize(.epred = mean(.epred))

p<-ggplot(fe_only,
       aes(x = BM.prey, y = .epred)) +
  stat_interval(aes(alpha = stat(level)), color = "royalblue4", .width = c(0.8, 0.95)) +
  geom_line(data = fe_only_mean, color = "royalblue4", size = 1) +
  scale_alpha_discrete(range = c(0.05, 0.1))+
  labs(x = "Prey body mass", y = "Interaction probability") +
  theme_minimal() + 
  theme(panel.grid = element_blank(), axis.line = element_line(size = 1),
        axis.text = element_blank(), legend.position = "none", axis.title = element_text(size = 12))
ggsave("figures/conceptual/logisiticreg.png", p, width = 2, height = 2)

# Map of the Food Webs ----------------------------------------------------
library(sf)
library(tidyr)
library(ggrepel)
library(maps)
library(rphylopic)
sf_use_s2(FALSE)

Europe <- st_read("data/raw/polygons/BiogeoRegions2016_shapefile/BiogeoRegions2016.shp") %>%
  filter(short_name != "outside")
Pyrenees <- st_read("data/raw/polygons/EuropeanMountainAreas/m_massifs_v1.shp") %>%
  filter(name_mm == "Pyrenees")
Serengeti <- st_read("data/raw/polygons/Serengeti_Ecosystem/v3_serengeti_ecosystem.shp") %>% st_union()
Nunavik <- st_read("data/raw/polygons/HighArctic/QuebecLabrador50.shp")

# polygons
world <- map_data('world')
#Europe <- st_transform(Europe, st_crs("EPSG:4326")) %>% st_union()
#st_write(Europe, "data/raw/polygons/Europe/Europe.shp")
Europe <- st_read("data/raw/polygons/Europe/Europe.shp")
Europe_centroid <- st_coordinates(st_centroid(Europe))
Pyrenees <- st_transform(Pyrenees, st_crs("EPSG:4326")) %>% st_union() 
Pyrenees_centroid <- st_coordinates(st_centroid(Pyrenees))
Nunavik <- st_transform(Nunavik, st_crs("EPSG:4326")) %>% st_union()
Nunavik_centroid <- cbind(st_coordinates(st_centroid(Nunavik)))
Serengeti <- st_transform(Serengeti, st_crs("EPSG:4326"))
Serengeti_centroid <- st_coordinates(st_centroid(Serengeti))

# curves
df <- expand_grid(from = c("Europe", "Pyrenees", "Nunavik", "Serengeti"), to = c("Europe", "Pyrenees", "Nunavik", "Serengeti")) %>%
  filter(from != to)
df[which(df$from == "Europe"), c("fromlong", "fromlat")] <- Europe_centroid
df[which(df$from == "Pyrenees"), c("fromlong", "fromlat")] <- Pyrenees_centroid
df[which(df$from == "Nunavik"), c("fromlong", "fromlat")] <- Nunavik_centroid
df[which(df$from == "Serengeti"), c("fromlong", "fromlat")] <- Serengeti_centroid
df[which(df$to == "Europe"), c("tolong", "tolat")] <- Europe_centroid
df[which(df$to == "Pyrenees"), c("tolong", "tolat")] <- Pyrenees_centroid
df[which(df$to == "Nunavik"), c("tolong", "tolat")] <- Nunavik_centroid
df[which(df$to == "Serengeti"), c("tolong", "tolat")]  <- Serengeti_centroid
df$size <- 0.5
df$size[which(df$from == "Europe" & df$to == "Nunavik")] <- 1


# species proportions
Europe_sp <- read.csv("data/cleaned/EuroMWTaxo.csv", row.names = 1)
Pyrenees_sp <- read.csv("data/cleaned/pyrenneesFWTaxo.csv", row.names = 1)
Nunavik_sp <- read.csv("data/cleaned/HighArcticFWTaxo.csv", row.names = 1)
Serengeti_sp <- read.csv("data/cleaned/SerengetiFWTaxo.csv", row.names = 1)
comp <- bind_rows(c(table(Europe_sp$Class) / nrow(Europe_sp), size = nrow(Europe_sp)),
          c(table(Pyrenees_sp$Class) / nrow(Pyrenees_sp), size = nrow(Pyrenees_sp)),
          c(table(Nunavik_sp$Class) / nrow(Nunavik_sp), size = nrow(Nunavik_sp)),
          c(table(Serengeti_sp$Class) / nrow(Serengeti_sp), size = nrow(Serengeti_sp))) %>%
  mutate(FW = c("Europe", "Pyrenees", "Nunavik", "Serengeti")) %>%
  pivot_longer(cols = -c(FW, size), values_to = "Proportion", names_to = "Class")
comp[is.na(comp)] <- 0
comp[which(comp$FW == "Europe"), c("long", "lat")] <- Europe_centroid
comp[which(comp$FW == "Pyrenees"), c("long", "lat")] <- Pyrenees_centroid
comp[which(comp$FW == "Nunavik"), c("long", "lat")] <- Nunavik_centroid
comp[which(comp$FW == "Serengeti"), c("long", "lat")] <- Serengeti_centroid
colorbar <- c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072")
euro_barchart <- ggplotGrob(
  ggplot(comp[comp$FW == "Europe",])+
    geom_bar(aes(x=Class, fill = Class, y = Proportion),
             position='dodge',stat='identity', color = "black") +
    labs(x = NULL, y = NULL) + 
    scale_fill_manual(values =colorbar) +
    theme(legend.position = "none", rect = element_blank(),
          line = element_blank(), text = element_blank()) 
)
pyrenees_barchart <- ggplotGrob(
  ggplot(comp[comp$FW == "Pyrenees",])+
    geom_bar(aes(x=Class, fill = Class, y = Proportion),
             position='dodge',stat='identity', color = "black") +
    labs(x = NULL, y = NULL) + 
    scale_fill_manual(values =colorbar) +
    theme(legend.position = "none", rect = element_blank(),
          line = element_blank(), text = element_blank()) 
)
nunavik_barchart <- ggplotGrob(
  ggplot(comp[comp$FW == "Nunavik",])+
    geom_bar(aes(x=Class, fill = Class, y = Proportion),
             position='dodge',stat='identity', color = "black") +
    labs(x = NULL, y = NULL) + 
    scale_fill_manual(values =colorbar) +
    theme(legend.position = "none", rect = element_blank(),
          line = element_blank(), text = element_blank()) 
)
serengeti_barchart <- ggplotGrob(
  ggplot(comp[comp$FW == "Serengeti",])+
    geom_bar(aes(x=Class, fill = Class, y = Proportion),
             position='dodge',stat='identity', color = "black") +
    labs(x = NULL, y = NULL) + 
    scale_fill_manual(values =colorbar) +
    theme(legend.position = "none", rect = element_blank(),
          line = element_blank(), text = element_blank()) 
)

legend <- data.frame(Class = c("Amphibia", "Aves", "Mammalia", "Reptilia"),
                     Proportion = c(0.15, 0.35, 0.25, 0.20))
frog <- image_data("c07ce7b7-5fb5-484f-83a0-567bb0795e18", size = "64")[[1]]
lizard <- image_data("83053aee-0f56-4cf3-bbfa-9207e6f13f46", size = "64")[[1]]
eagle <- image_data("92589388-08e3-422f-b452-aa7454411a9c", size = "64")[[1]]
lynx <- image_data("24f763a3-accf-44c9-9a08-71e9834047b7", size = "64")[[1]]

legend_barchart <- ggplotGrob(
  ggplot(legend)+
    geom_bar(aes(x=Class, fill = Class, y = Proportion),
             position='dodge', stat='identity', color = "black") +
    labs(x = NULL, y = NULL, title = "Community composition") + 
    scale_fill_manual(values = colorbar) +
    add_phylopic(frog, x = 1, y = legend$Proportion[1]+0.06, alpha = 1, ysize = 0.7, color = "black") +
    add_phylopic(eagle, x = 2, y = legend$Proportion[2]+0.02, alpha = 1, ysize = 0.8, color = "black") +
    add_phylopic(lynx, x = 3, y = legend$Proportion[3]+0.06, alpha = 1, ysize = 0.65, color = "black") +
    add_phylopic(lizard, x = 4, y = legend$Proportion[4]+0.03, alpha = 1, ysize = 0.55, color = "black") +
    lims(y = c(0,0.4))+
    theme(legend.position = "none", rect = element_blank(),
          line = element_blank(), text = element_blank(), title = element_text(size = 10, face = "bold"), plot.margin = unit(c(0, 0, 0, 0), "null")) 
)

p <- ggplot() +
  geom_polygon(data =world, aes(long, lat, group = group), color = "grey90", fill = "grey90") +
  geom_sf(data = Europe, fill = alpha("royalblue4", 0.5), color = "royalblue4") +
  geom_sf(data = Pyrenees, fill = alpha("red3", 0.5), color = "red3") +
  geom_sf(data = Nunavik, fill = alpha("deepskyblue", 0.5), color = "deepskyblue") +
  geom_sf(data = Serengeti, fill = alpha("chartreuse4", 0.5), color =  "chartreuse4") +
  geom_curve(data=df,
             aes(x=fromlong, y=fromlat, xend=tolong, yend=tolat, size = size), color = "black",
             curvature=0.25, alpha = 0.6) + 
  scale_colour_identity() +
  scale_size(range = c(0.5,2))+
  geom_label_repel(data = as.data.frame(Pyrenees_centroid), aes(x = X, y = Y, label = "Pyrenees"), 
                   fontface = "bold", nudge_x = -15, nudge_y = -10, colour = "red3") +
  geom_label_repel(data = as.data.frame(Europe_centroid), aes(x = X, y = Y, label = "Europe"), 
                   fontface = "bold", nudge_x = 13, nudge_y = 10, colour = "royalblue4") +
  geom_label_repel(data = as.data.frame(Nunavik_centroid), aes(x = X, y = Y, label = "Nunavik"), 
                   fontface = "bold", nudge_x = -5, nudge_y = -10, colour = "deepskyblue") +
  geom_label_repel(data = as.data.frame(Serengeti_centroid), aes(x = X, y = Y, label = "Serengeti"), 
                   fontface = "bold", nudge_x = 23, nudge_y = 20, colour = "chartreuse4") +
  annotation_custom(euro_barchart, xmin = Europe_centroid[,1]+3, xmax = Europe_centroid[,1]+23, 
                    ymin = Europe_centroid[,2]-13, ymax = Europe_centroid[,2]+8) + 
  annotation_custom(pyrenees_barchart, xmin = Pyrenees_centroid[,1]-25, xmax = Pyrenees_centroid[,1]-5, 
                    ymin = Pyrenees_centroid[,2]-27, ymax = Pyrenees_centroid[,2]-12) + 
  annotation_custom(nunavik_barchart, xmin = Nunavik_centroid[,1]-15, xmax = Nunavik_centroid[,1]+5, 
                    ymin = Nunavik_centroid[,2]-27, ymax = Nunavik_centroid[,2]-12) + 
  annotation_custom(serengeti_barchart, xmin = Serengeti_centroid[,1]+13, xmax = Serengeti_centroid[,1]+33, 
                    ymin = Serengeti_centroid[,2]-3, ymax = Serengeti_centroid[,2]+18) + 
  annotation_custom(legend_barchart, xmin = -90, xmax = -40, 
                    ymin = -5, ymax = 25) + 
  coord_sf(xlim = c(-85, 70), ylim = c(0, 80)) +
  theme_void() +
  theme(legend.position = "none")

ggsave("figures/conceptual/FWmap.png", p, scale = 2)

# # Performance ~ distances mini-maps ---------------------------------------
# overall_performance <- read.csv("data/checkpoints/overall_performance.csv", row.names = 1)
# overall_performance[overall_performance == "Arctic"] <- "HighArctic"
# overall_performance[overall_performance == "Euro"] <- "Europe"
# FWdist <- read.csv("data/checkpoints/FWdist.csv", row.names = 1)
# FWdist[FWdist == "High Arctic"] = "HighArctic"
# 
# df <- df %>%
#   left_join(overall_performance, by = c("from" = "Source", "to" = "Target")) %>%
#   left_join(FWdist, by = c("from" = "FW2", "to" = "FW1"))
# 
# centroids <- rbind(Europe_centroid, HighArctic_centroid, Pyrenees_centroid, Serengeti_centroid) %>%
#   as_tibble() %>%
#   mutate(col = c("royalblue4", "deepskyblue", "red3", "chartreuse4"))
# perform <- ggplot(df) +
#   geom_curve(data=df,
#              aes(x=fromlong, y=fromlat, xend=tolong, yend=tolat, color = col, size = auc),
#              curvature=0.25, alpha = 0.6) +
#   scale_size_continuous(range = c(0.5,12), limits = c(0.6,1))+
#   xlim(c(min(df$fromlong), max(df$fromlong)+5))+
#   geom_point(data = centroids, aes(x = X, y = Y, colour = col), size = 12, shape = 21, stroke = 5, fill = "white")+
#   scale_color_identity() +
#   theme_void() +
#   theme(legend.position = "none")
# 
# geo_dist <- ggplot(df) +
#   geom_curve(data=df,
#              aes(x=fromlong, y=fromlat, xend=tolong, yend=tolat, color = col, size = geo.dist),
#              curvature=0.25, alpha = 0.6) +
#   scale_size_continuous(range = c(0.5,12), limits = c(3500, 18000))+
#   xlim(c(min(df$fromlong), max(df$fromlong)+5))+
#   geom_point(data = centroids, aes(x = X, y = Y, colour = col), size = 12, shape = 21, stroke = 5, fill = "white")+
#   scale_color_identity() +
#   theme_void() +
#   theme(legend.position = "none")
#   
# env_dist <- ggplot(df) +
#   geom_curve(data=df,
#              aes(x=fromlong, y=fromlat, xend=tolong, yend=tolat, color = col, size = env.dist),
#              curvature=0.25, alpha = 0.6) +
#   scale_size_continuous(range = c(0.5,12), limits = c(1,10))+
#   xlim(c(min(df$fromlong), max(df$fromlong)+5))+
#   geom_point(data = centroids, aes(x = X, y = Y, colour = col), size = 12, shape = 21, stroke = 5, fill = "white")+
#   scale_color_identity() +
#   theme_void() +
#   theme(legend.position = "none")
# 
# phylo_dist <- ggplot(df) +
#   geom_curve(data=df,
#              aes(x=fromlong, y=fromlat, xend=tolong, yend=tolat, color = col, size = phylo.dist),
#              curvature=0.25, alpha = 0.6) +
#   scale_size_continuous(range = c(0.5,12), limits= c(0, 200))+
#   xlim(c(min(df$fromlong), max(df$fromlong)+5))+
#   geom_point(data = centroids, aes(x = X, y = Y, colour = col), size = 12, shape = 21, stroke = 5, fill = "white")+
#   scale_color_identity() +
#   theme_void() +
#   theme(legend.position = "none")
# 
# ggsave("figures/conceptual/performance_map.png", perform)
# ggsave("figures/conceptual/geodist_map.png", geo_dist)
# ggsave("figures/conceptual/envdist_map.png", env_dist)
# ggsave("figures/conceptual/phydist_map.png", phylo_dist)
