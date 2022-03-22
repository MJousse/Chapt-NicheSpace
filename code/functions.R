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

make_predictions <- function(predictors, order_level, coef, model){
  x <- select(predictors, 
              -Predator, -Prey, -Order.predator, -Order.prey,
              -Herbivore.predator, -Herbivore.prey, -ClutchSize.predator,
              -ClutchSize.prey)
  
  x <- cbind(rep(1, nrow(predictors)), x) %>% as_data()
  predator_order <- as.numeric(factor(predictors$Order.pred, levels = order_level))
  linear_predictor <- rowSums(x * t(coef[,predator_order]))
  
  p <- ilogit(linear_predictor)
  
  predictions <- calculate(p, values = model, nsim = 100)
  predictions <- colMeans(predictions[[1]])
  predictions <- data.frame(Predator = as.factor(predictors$Predator), Prey = as.factor(predictors$Prey), prediction = predictions)
  return(predictions)
}
