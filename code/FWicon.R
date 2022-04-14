library(igraph)
library(NetIndices)

arcticFW <- read.csv("data/cleaned/HighArcticFW.csv", row.names = 1)
arcticFW <- arcticFW %>%
  transmute(resource = Prey, consumer = Predator)
arctic.igraph <- graph_from_edgelist(as.matrix(arcticFW))
arctic.matrix <- Matrix::as.matrix(as_adjacency_matrix(arctic.igraph))

lay<-matrix(nrow = nrow(arctic.matrix), ncol=2) # create a matrix with one column as runif, the other as trophic level
lay[,1]<-sample(c(1:nrow(arctic.matrix)), nrow(arctic.matrix))
lay[,2]<-scale(TrophInd(arctic.matrix)$TL-1) + rnorm(nrow(arctic.matrix), mean =0,sd = 0.2)

plot(arctic.igraph, layout = lay, vertex.label=NA, vertex.size=10, edge.arrow.size=0, edge.width=.5,
     vertex.color = "deepskyblue", edge.color = alpha("deepskyblue", 0.3))

EuropeMW <- read.csv("data/cleaned/EuroFWadults.csv", row.names = 1)
sp <- unique(c(EuropeMW$Predator, EuropeMW$Prey))[sample(c(1:200))]
EuropeMW <- EuropeMW %>%
  transmute(resource = Prey, consumer = Predator) %>%
  filter(resource %in% sp, consumer %in% sp)

europe.igraph <- graph_from_edgelist(as.matrix(EuropeMW))
europe.matrix <- Matrix::as.matrix(as_adjacency_matrix(europe.igraph))

lay<-matrix(nrow = nrow(europe.matrix), ncol=2) # create a matrix with one column as runif, the other as trophic level
lay[,1]<-sample(c(1:nrow(europe.matrix)), nrow(europe.matrix))
lay[,2]<-scale(TrophInd(europe.matrix)$TL-1) + rnorm(nrow(europe.matrix), mean =0,sd = 0.2)

plot(europe.igraph, layout = lay, vertex.label=NA, vertex.size=10, edge.arrow.size=0, edge.width=.5,
     vertex.color = "royalblue4", edge.color = alpha("royalblue4", 0.3))

SerengetiFW <- read.csv("data/cleaned/SerengetiFW.csv", row.names = 1)
sp <- unique(c(SerengetiFW$Consumer_Species, SerengetiFW$Resource_Species))[sample(c(1:100))]
SerengetiFW <- SerengetiFW %>%
  transmute(resource = Resource_Species, consumer = Consumer_Species) %>%
  filter(resource %in% sp, consumer %in% sp)

serengeti.igraph <- graph_from_edgelist(as.matrix(SerengetiFW))
serengeti.matrix <- Matrix::as.matrix(as_adjacency_matrix(serengeti.igraph))

lay<-matrix(nrow = nrow(serengeti.matrix), ncol=2) # create a matrix with one column as runif, the other as trophic level
lay[,1]<-sample(c(1:nrow(serengeti.matrix)), nrow(serengeti.matrix))
lay[,2]<-scale(TrophInd(serengeti.matrix)$TL-1) + rnorm(nrow(serengeti.matrix), mean =0,sd = 0.2)

plot(serengeti.igraph, layout = lay, vertex.label=NA, vertex.size=10, edge.arrow.size=0, edge.width=.5,
     vertex.color = "chartreuse4", edge.color = alpha("chartreuse4", 0.3))

PyreneesFW <- read.csv("data/cleaned/pyrenneesFW.csv", row.names = 1)
sp <- unique(c(PyreneesFW$Predator, PyreneesFW$Prey))[sample(c(1:100))]
PyreneesFW <- PyreneesFW %>%
  transmute(resource = Prey, consumer = Predator) %>%
  filter(resource %in% sp, consumer %in% sp)

pyrenees.igraph <- graph_from_edgelist(as.matrix(PyreneesFW))
pyrenees.matrix <- Matrix::as.matrix(as_adjacency_matrix(pyrenees.igraph))

lay<-matrix(nrow = nrow(pyrenees.matrix), ncol=2) # create a matrix with one column as runif, the other as trophic level
lay[,1]<-sample(c(1:nrow(pyrenees.matrix)), nrow(pyrenees.matrix))
lay[,2]<-scale(TrophInd(pyrenees.matrix)$TL-1) + rnorm(nrow(pyrenees.matrix), mean =0,sd = 0.3)

plot(pyrenees.igraph, layout = lay, vertex.label=NA, vertex.size=10, edge.arrow.size=0, edge.width=.5,
     vertex.color = "red3", edge.color = alpha("red3", 0.2))
