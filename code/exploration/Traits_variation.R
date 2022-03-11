rm(list = ls())
library(lme4)
library(nnet)
library(dplyr)
library(performance)

EuroTraits <- read.csv("data/cleaned/SpeciesTraitsFull.csv", row.names = 1, stringsAsFactors = T) %>%
  mutate(TrophicLevel = as.difelse(Herbivore == 1, "Herbivore", 
                               ifelse(Carnivore == 1, "Carnivore", 
                                      ifelse(Omnivore ==1, "Omnivore", NA))), .keep = "unused")

traits_col <- c(2:15, 20:23)

df <- data.frame()
for (itrait in traits_col){
  trait <- colnames(EuroTraits)[itrait]
  formula <- as.formula(paste0(trait, " ~ (1 | Class/Order/Family/Genus)"))
  if (length(unique(EuroTraits[,itrait])) == 2){
    model <- brm(formula, 
                 data = EuroTraits, family = bernoulli())
    var_partition <- icc(model, by_group = T)
  } else if (trait == "TrophicLevel"){
    model <- brm(formula, data = EuroTraits, family = categorical())
    var_partition <- icc(model, by_group = T)
  } else {
    model <- brm(formula, data = EuroTraits, family = gaussian())
    var_partition <- icc(model, by_group = T)
  }
  df <- rbind(df,
              c(trait, var_partion$ICC, (1 - sum(ICC))))
}
colnames(df) <- c("Trait", "Genus", "Family", "Order", "Class", "Residuals")
