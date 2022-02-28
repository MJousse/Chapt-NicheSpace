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

#' Match of Predator's Diet with the Prey
#' 
#' @param Small_Mam logical vector. Predators eats small Mammals?
#' @param Large_Mam logical vector. Predators eats large Mammals?
#' @param Herptile logical vector. Predators eats herptiles?
#' @param Bird_eggs logical vector. Predators eats bird's eggs?
#' @param Small_bird logical vector. Predators eats small birds?
#' @param Large_Bird logical vector. Predators eats large birds?
#' @param Class_prey character vector. Class of the prey.
#' @return A Logical vector. Does the diet of the predator match the class of the prey?
match_diet <- function(Small_Mam, Large_Mam, Herptile, Bird_eggs, Small_bird, Large_Bird, Class_prey){
  match <- ifelse(Class_prey == "Amphibia" & Herptile %in% c(1,2), 1,
                  ifelse(Class_prey == "Aves" & (Bird_eggs %in% c(1,2) | Small_bird %in% c(1,2) | Large_Bird %in% c(1,2)), 1,
                         ifelse(Class_prey == "Mammalia" & (Small_Mam %in% c(1,2) | Large_Mam %in% c(1,2)), 1,
                                ifelse(Class_prey == "Reptilia" & Herptile %in% c(1,2), 1, 0))))
  return(match)
}

#' Get the mode of a vector
#' 
#' @param v a vector
#' @return the value with the most occurrence
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

#' Transform traits into predictors
#' 
#' TODO
traits2predictors <- function(FW, prey_suffix = ".x", predator_suffix = ".y"){
  FW <- FW %>%
    rename_with(~ gsub(paste0("\\", predator_suffix, "$"), ".predator", .x)) %>%
    rename_with(~ gsub(paste0("\\", prey_suffix, "$"), ".prey", .x))
  foraging_predictor <- FW %>%
    select(Trophic_level.predator, Habitat_breadth.predator = Habitat_breadth_IUCN.predator, 
           Order.predator, BM.predator = logBM.predator, Longevity.predator = Max_longevity_d.predator,
           ClutchSize.predator = Litter_clutch_size.predator)
  
  vulnerability_predictor <- FW %>%
    select(Trophic_level.prey, Habitat_breadth.prey = Habitat_breadth_IUCN.prey, 
           Order.prey, BM.prey = logBM.prey, Longevity.prey = Max_longevity_d.prey, 
           ClutchSize.prey = Litter_clutch_size.prey)
  
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
    mutate_at(vars(Trophic_level.predator, Trophic_level.prey, ActivityTime.match), as.factor)
  return(predictors)
}

scale2 <- function(x){(x - mean(x)) / (2*sd(x))}

