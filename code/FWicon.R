library(igraph)
library(dplyr)
library(NetIndices)
library(ggplot2)
set.seed(1234)
par(bg=NA)
# Actual food webs --------------------------------------------------------
# arctic
arcticFW <- read.csv("data/cleaned/HighArcticFW.csv", row.names = 1)
arcticFW <- arcticFW %>%
  transmute(resource = Prey, consumer = Predator)
arctic.igraph <- simplify(graph_from_edgelist(as.matrix(arcticFW)))
arctic.matrix <- Matrix::as.matrix(as_adjacency_matrix(arctic.igraph))

lay<-matrix(nrow = nrow(arctic.matrix), ncol=2) # create a matrix with one column as runif, the other as trophic level
lay[,1]<-sample(c(1:nrow(arctic.matrix)), nrow(arctic.matrix))
lay[,2]<-TrophInd(arctic.matrix)$TL-1 + runif(nrow(arctic.matrix), min =0, max = 1)

png("figures/FWgraphs/arctic.png")
par(bg=NA)
plot(arctic.igraph, layout = lay, vertex.label=NA, vertex.size=10, edge.arrow.size=0, edge.width=4,
     vertex.color = "deepskyblue", edge.color = alpha("deepskyblue", 0.2))
dev.off()

# mock prediction
png("figures/FWgraphs/predictedArctic.png")
par(bg=NA)
plot(arctic.igraph, layout = lay, vertex.label=NA, vertex.size=10, edge.arrow.size=0, edge.width=runif(length(E(arctic.igraph)), min = 0, max = 10),
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
png("figures/FWgraphs/simplifiedEurope.png")
par(bg=NA)
plot(europe.simplified, layout = lay, vertex.label=NA, vertex.size=10, edge.arrow.size=0, edge.width=5,
     vertex.color = "royalblue4", edge.color = alpha("royalblue4", 0.2))
dev.off()

# mock prediction
png("figures/FWgraphs/predictedEurope.png")
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
png("figures/FWgraphs/predictedSerengeti.png")
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
png("figures/FWgraphs/predictedPyrenees.png")
par(bg=NA)
plot(pyrenees.simplified, layout = lay, vertex.label=NA, vertex.size=10, edge.arrow.size=0, edge.width=runif(length(E(pyrenees.simplified)), min = 0, max = 10),
     vertex.color = "red3", edge.color = alpha("red3", 0.2))
dev.off()
