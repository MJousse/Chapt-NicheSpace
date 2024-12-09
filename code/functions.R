# Functions used throughout the different scripts

#' Jaccard similarity
#' 
#' @param pred_df dataframe. Traits of the predator.
#' @param prey_df dataframe. Traits of the prey.
#' @return a vector indicating the similarity between the traits of the predator and the prey.
Jaccard <- function(pred_df, prey_df){
  jaccard_s <- c()
  a <- rowSums(pred_df == 1 & prey_df == 1)
  b <- rowSums(pred_df == 1 & prey_df == 0)
  c <- rowSums(pred_df == 0 & prey_df == 1)
  jaccard_s <- a/(a+b+c)
  jaccard_s[is.na(jaccard_s)] <- 0.5
  return(jaccard_s)
}

#' Get the mode of a vector
#' 
#' @param v a vector
#' @return the value with the most occurrence
#' match: match returns a vector of the positions of (first) matches of its first argument in its second.
#'        %in% is a more intuitive interface as a binary operator, which returns a logical vector indicating if there is a match or not for its left operand.
#' which.max: Determines the location, i.e., index of the (first) minimum or maximum of a numeric (or logical) vector.
#' tabulate: tabulate takes the integer-valued vector bin and counts the number of times each integer occurs in it.
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

#' Transform traits into predictors
#' 
#' @param FW dataframe. Traits of the the predators and the prey.
#' @param prey_suffix character. Suffix after all column names for the traits of the prey.
#' @param predator_suffix character. Suffix after all column names for the traits of the predator.
#' @return a dataframe with the predictors of species interactions.
traits2predictors <- function(FW, prey_suffix = ".x", predator_suffix = ".y"){
  FW <- FW %>%
    rename_with(~ gsub(paste0("\\", predator_suffix, "$"), ".predator", .x)) %>%
    rename_with(~ gsub(paste0("\\", prey_suffix, "$"), ".prey", .x))
  foraging_predictor <- FW %>%
    select(Herbivore.predator, Omnivore.predator, Carnivore.predator,
           Habitat_breadth.predator = Habitat_breadth_IUCN.predator, 
           Order.predator, BM.predator = logBM.predator, Longevity.predator = logLongevity.predator,
           ClutchSize.predator = logClutchSize.predator)
  
  vulnerability_predictor <- FW %>%
    select(Herbivore.prey, Omnivore.prey, Carnivore.prey,
           Habitat_breadth.prey = Habitat_breadth_IUCN.prey, 
           Order.prey, BM.prey = logBM.prey, Longevity.prey = logLongevity.prey, 
           ClutchSize.prey = logClutchSize.prey)
  
  match_predictor <- data.frame(
    ActivityTime.match = FW$Diel_activity.prey == FW$Diel_activity.predator,
    Habitat.match = Jaccard(
      select_at(FW, vars(Forest.predator:Introduced.vegetation.predator)), 
      select_at(FW, vars(Forest.prey:Introduced.vegetation.prey))),
    BM.match = (FW$logBM.predator - FW$logBM.prey)^2
  )
  predictors <- data.frame(Predator = FW$Predator,
                           Prey = FW$Prey) %>%
    cbind(foraging_predictor) %>%
    cbind(vulnerability_predictor) %>%
    cbind(match_predictor) %>%
    mutate(ActivityTime.match = as.integer(ActivityTime.match))
  return(predictors)
}

scale2 <- function(x){(x - mean(x)) / (2*sd(x))}

get_predictors <- function(Species_List, FuncTraits){
  FW <- expand.grid(Species_List, Species_List)
  colnames(FW) <- c("Predator", "Prey")
  FW <- left_join(FW, FuncTraits, by = c("Prey" = "Species")) %>%
    left_join(FuncTraits, by = c("Predator" = "Species"))
  FW <- traits2predictors(FW)
  return(FW)
}

species_role <- function(FW, ncores = 4, nsim =10){
  # remove self-loop
  FW <- FW[FW$resource != FW$consumer,]
  # centrality role
  graph <- graph_from_edgelist(as.matrix(FW[,c("resource","consumer")])) # igraph
  nodes <- vertex.attributes(graph)$name # save species name
  indegree <- centr_degree(graph, mode = "in")$res
  outdegree <- centr_degree(graph, mode = "out")$res
  betweeness <- centr_betw(graph)$res
  closeness <- centr_clo(graph)$res
  eigen <- centr_eigen(graph)$vector
  
  # trophic level and omnivory
  tl <- trophiclevels(graph)
  
  # motifs role
  m <- get.adjacency(graph,sparse=FALSE)
  motif_role <-motif_role(m)
  # normalized
  #motif_role$position_count <- sweep(motif_role$position_count, MARGIN = 1, FUN = "/", rowSums(motif_role$position_count))
  
  # module-based role
  modulerole <- calc_topological_roles(graph, ncores = ncores, nsim = nsim) %>%
    group_by(node) %>%
    summarise(within_module_degree = median(within_module_degree, na.rm = T),
              among_module_conn = median(among_module_conn, na.rm = T))
  
  return(data.frame(species = nodes,
                   indegree,
                   outdegree,
                   betweeness,
                   closeness,
                   eigen,
                   tl,
                   modulerole[,c(2,3)],
                   motif_role$position_count,
                   row.names = c(1:length(nodes)))
  )
}

fw_properties <- function(FW, nsim){
  # centrality role
  graph <- graph_from_edgelist(as.matrix(FW[,c("resource","consumer")])) # igraph
  connectance <- ecount(graph)/vcount(graph)^2
  TLs <- trophiclevels(graph)
  meanTL <- mean(TLs$TL)
  maxTL <- max(TLs$TL)
  motifs <- motifs(graph, size = 3) 
  motifs <- motifs[c(5, 8, 12, 7, 3, 9, 14, 6, 10, 13, 16, 15, 11) ]# reorganize as in Stouffer et al 2012
  names(motifs) <- paste0("motif", c(1:13))
  diameter <- diameter(graph)
  n_clusters <- c()
  modularity <- c()
  for (i in c(1:nsim)){
    clusters <- cluster_spinglass(graph)
    n_clusters <- c(n_clusters, length(unique(clusters$membership)))
    modularity <- c(modularity, clusters$modularity)
  }
  return(c(connectance = connectance,
           meanTL = meanTL,
           maxTL = maxTL,
           motifs,
           diameter = diameter,
           n_clusters = n_clusters[which.max(modularity)],
           modularity = modularity[which.max(modularity)])
  )
}

make_predictions <- function(Model, newdata, ndraws = 100, allow_new_levels = TRUE){
  predictions <- predict(Model, newdata = newdata, allow_new_levels = allow_new_levels, ndraws = ndraws, summary = F)
  rownames(predictions) <- paste0("draws", c(1:ndraws))
  predictions <- as.data.frame(t(predictions))
  predictions$Estimate <- apply(predictions, MARGIN = 1, mean)
  predictions$Est.Error <- apply(predictions[,-ncol(predictions)], MARGIN = 1, sd)
  predictions <- select(newdata, Predator, Prey, interaction) %>%
    bind_cols(predictions)
  return(predictions)
}

trophiclevels <- function(graph){
  m <- fix_basal(graph) # the algorithm fails when the only prey of a group of species is only each other
  TrophInd(m)
}

fix_basal <- function(graph){
  # There is a problem with calculating trophic level when the prey of a group of species is only each other
  # for example: if C eats A and B, and A and B eat each other
  # m <- matrix(c(1,1,1,0,1,1,0,1,1), nrow = 3); TrophInd(m) returns negative trophic levels
  # solution is to remove the interactions between 'problematic' species (A and B in example)
  m <- get.adjacency(graph,sparse=FALSE)
  # find non-problematic basal species
  basal <- which(colSums(m)==0)
  # shortest path to basal
  sptb <- shortest.paths(graph,to = V(graph)[basal], mode = "in")
  sptb[is.infinite(sptb)] <- NA
  sptb <- apply(sptb, MARGIN = 1, min, na.rm = T)
  if (any(is.infinite(sptb))){
    # problematic species do not have paths to non-problematic basal species
    probl_sp <- which(is.infinite(sptb))
    if (length(probl_sp) > 1){
      # problematic *basal* species are the one that have the lowest trophic level
      tls <- round(TrophInd(m[probl_sp, probl_sp]),5)
      probl_basal <- rownames(tls)[tls$TL == min(tls$TL)]
      # remove interactions between problematic basal species
      m[probl_basal, probl_basal] <- 0
    } else {
      m[probl_sp, probl_sp] <- 0
    }
  }
  return(m)
}

performance_draws <- function(draws, interactions){
  out <- data.frame()
  for (i in c(1:ncol(draws))){
    tp <- sum(draws[,i] == 1 & interactions == 1)
    tn <- sum(draws[,i] == 0 & interactions == 0)
    fp <- sum(draws[,i] == 1 & interactions == 0)
    fn <- sum(draws[,i] == 0 & interactions == 1)
    out <- rbind(out, data.frame(draw = i, tpr = tp/(tp+fn), tnr = tn/(tn+fp), ppv = tp/(tp+fp), npv = tn/(fn+tn)))
  }
  return(out)
}
