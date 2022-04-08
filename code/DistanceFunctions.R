# Compositional distance
JaccardDissimilarity <- function(x, y){
  a = sum(x %in% y)
  b = sum(!(x %in% y))
  c = sum(!(y %in% x))
  return(1 - (a / (a+b+c)))
}

# Mean Nearest Taxon among communities
phylobetadiv <- function(community, phydist){
  communities <- rownames(community)
  ncommunity <- nrow(community)
  df <- matrix(0L, nrow = ncommunity, ncol = ncommunity)
  colnames(df) <- communities
  rownames(df) <- communities
  for (i in c(1:ncommunity)){
    for (j in c(1:ncommunity)){
      if (i != j){
        ispecies <- colnames(community)[community[i,] == 1]
        jspecies <- colnames(community)[community[j,] == 1]
        sp_mntd <- map_dbl(jspecies, mntd, ispecies, phydist)
        df[i,j] <- mean(sp_mntd, na.rm =  T)
      }
    }
  }
  return(df)
}

fnnd <- function(target_species, species_pool, funcdist){
  funcdist <- as.matrix(funcdist)
  if (target_species %in% colnames(funcdist)){
    min(funcdist[target_species, colnames(funcdist) %in% species_pool])
  } else {NA}
}

fmpd <- function(target_species, species_pool, funcdist){
  funcdist <- as.matrix(funcdist)
  if (target_species %in% colnames(funcdist)){
    mean(funcdist[target_species, colnames(funcdist) %in% species_pool])
  } else {NA}
}


mntd <- function(target_species, species_pool, phylodist){ # Shortest phylo distance
  if (target_species %in% colnames(phylodist)){
    min(phylodist[target_species, colnames(phylodist) %in% species_pool])
  } else {NA}
}