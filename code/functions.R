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
