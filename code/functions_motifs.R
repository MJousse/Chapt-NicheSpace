motif_role <- function(m){
  sp <- colnames(m)
  n <- nrow(m)
  out <- matrix(0, nrow = n, ncol = 30)
  rownames(out) <- sp
  colnames(out) <- paste0("position", c(1:30))
  motifs <- rep(0,16)
  for (i in c(1:(n-2))){
    for (j in c((i+1):(n-1))){
      for (k in c((j+1):n)){
        subgraph <- m[c(i,j,k), c(i,j,k)]
        diag(subgraph) <- 0
        subigraph <- graph_from_adjacency_matrix(subgraph)
        if (is_connected(subigraph)){
          motif_id <- which(motifs(subigraph) == 1)
          motifs[motif_id] <- motifs[motif_id] + 1
          pos <- positions(subgraph, motif_id)
          out[i, pos[1]] = out[i, pos[1]] + 1
          out[j, pos[2]] = out[j, pos[2]] + 1
          out[k, pos[3]] = out[k, pos[3]] + 1
        }
      }
    }
  }
  return(list(motif_count = motifs, position_count = out))
}

positions <- function(subgraph, motif_id){
  pos <- rep(NA, 3)
  if (motif_id == 5){
    pos[which(rowSums(subgraph) == 0)] <- 1
    pos[which(rowSums(subgraph) == 1 && colSums(subgraph) == 1)] <- 2
    pos[which(colSums(subgraph) == 0)] <- 3
  } else if (motif_id == 8){
    pos[which(colSums(subgraph) == 2)] <- 4
    pos[which(rowSums(subgraph) == 2)] <- 5
    pos[which(rowSums(subgraph) == 1)] <- 6
  } else if (motif_id == 12){
    pos <- rep(7, 3)
  } else if (motif_id == 7){
    pos <- rep(8,3)
    pos[which(rowSums(subgraph) == 2)] <- 9
  } else if (motif_id == 3){
    pos <- rep(10,3)
    pos[which(colSums(subgraph) == 2)] <- 11
  } else if (motif_id == 9){
    pos <- rep(12,3)
    pos[which(colSums(subgraph) == 0)] <- 13
  } else if (motif_id == 14){
    pos <- rep(14,3)
    pos[which(colSums(subgraph) == 2)] <-15
  } else if (motif_id == 6){
    pos[which(colSums(subgraph) == 1)] <- 16
    pos[which(colSums(subgraph) == 0)] <- 17
    pos[which(colSums(subgraph) == 2)] <- 18
  } else if (motif_id == 10){
    pos[which(rowSums(subgraph) == 2)] <- 19
    pos[which(rowSums(subgraph) == 0)] <- 20
    pos[which(rowSums(subgraph) == 1)] <- 21
  } else if (motif_id == 13){
    pos[which(colSums(subgraph) == 2)] <- 22
    pos[which(rowSums(subgraph) == 1 && colSums(subgraph) == 1)] <- 23
    pos[which(rowSums(subgraph) == 2)] <- 24
  } else if (moti_id == 16){
    pos <- rep(25,3)
  } else if (motif_id == 15){
    pos[which(rowSums(subgraph) == 1 && colSums(subgraph) == 2)] <- 26
    pos[which(rowSums(subgraph) == 2 && colSums(subgraph) == 2)] <- 27
    pos[which(rowSums(subgraph) == 2 && colSums(subgraph) == 1)] <- 28
  } else if (motif_id == 11){
    pos <- rep(29,3)
    pos[which(rowSums(subgraph) == 2 && colSums(subgraph) == 2)] <- 30
  }
  return(pos)
}


A = matrix(c(0,1,1,0,0,0,0,1,0), nrow = 3)
A = graph_from_adjacency_matrix(A)
m <- which(motifs(A) == 1)
plot(A, vertex.size = 50, layout=layout_as_tree(A, flip.y = T), vertex.label = "")
plot(A, vertex.size = 50, layout=layout_nicely(A), vertex.label = "")
title(m,cex.main=3)
