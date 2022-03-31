motif_role(m){
  sp <- colnames(m)
  n <- nrow(m)
  out <- matrix(0, nrow = n, ncol = 30)
  rownames(out) <- sp
  colnames(out) <- paste0("position", c(1:30))
  for (i in c(1:n)){
    for (j in c(1:n)){
      if (i != j){
        for (k in c(1:n)){
          if (k != i & k != j){
            subgraph <- m[c(i,j,k), c(i,j,k)]
            subigraph <- graph_from_adjacency_matrix(subgraph)
            if (is_connected(subigraph)){
              motif_id <- which(motifs(subigraph) == 1)
              sp_positions <- positions(motif_id, subgraph)
            }
          }
        }
      }
    }
  }
}

positions(motif_id, subgraph){
  diag(subgraph) <- 0
  if(motif_id == 3 ){
    pred <- which.max(colSums(subgraph))
    positions <- rep(10,3)
    positions[pred] <- 11
    return(postions)
  } else if (motif_id == 5){
    positions <- rep(1,3)
    prey <- which(colSums(subgraph) == 0)
    int <- which(subgraph[prey,] == 1)
    positions[int] <- 2
    positions[prey] <- 3
  } else if (motif_id == 6)
}