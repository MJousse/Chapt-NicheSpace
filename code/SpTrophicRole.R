library(igraph)
library(cheddar)
library(bmotif)
library(multiweb)

arcticFW <- read.csv("data/cleaned/HighArcticFW.csv", row.names = 1)
pyreneesFW <- read.csv("data/cleaned/pyrenneesFW.csv", row.names = 1)


FW <- arcticFW

# centrality role
graph <- graph_from_edgelist(as.matrix(FW[,c(2,1)]))

nodes <- vertex.attributes(graph)$name
indegree <- centr_degree(graph, mode = "in")$res
outdegree <- centr_degree(graph, mode = "out")$res
betweeness <- centr_betw(graph)$res
closeness <- centr_clo(graph)$res
eigen <- centr_eigen(graph)$vector

# trophic level
ched_community <- Community(nodes = data.frame(node = nodes), 
                            trophic.links = data.frame(resource = FW[,2],
                                                       consumer = FW[,1]),
                            properties = list(title = "FW"))
tl <- TrophicLevels(ched_community)

# motifs role TODO
m <- as.matrix(as_adjacency_matrix(graph))
motif_role <- node_positions(m, weights_method = "none")
motif_role <- sweep(motif_role, MARGIN = 1, apply(attitude, 1, sum), FUN = "/")

# module-based role
modulerole <- calc_topological_roles(graph, ncores = 4)
