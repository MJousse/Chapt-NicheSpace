# set up
rm(list = ls())
library(ROCR)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
load("data/checkpoints/BRT_predictions.RData")

foodwebs <- c("arctic", "europe", "pyrenees", "serengeti")
source("code/functions.R")

# Calculate performance ---------------------------------------------------
overall_performance <- expand_grid(Source = foodwebs, Target = foodwebs) %>%
  mutate(auc = NA, aucpr = NA, prevalence = NA)
overall_performance_draws <- data.frame()
species_performance <- data.frame()

# for each combinaison calculate the overall auc and aucpr
# in each target food web, calculate the auc and aucpr of each species
for (combination in c(1:nrow(overall_performance))){
  sourceFW <- overall_performance[combination, "Source"]
  targetFW <- overall_performance[combination, "Target"]
  predictions <- get(paste0(sourceFW, "BRT_", targetFW, "FW"))
  if (sourceFW == targetFW){
    predictions_test <- filter(predictions, testing == 1)
  } else {
    predictions_test <- predictions
  }
  overall_performance[combination, "auc"] <- performance(prediction(predictions_test$prediction, predictions_test$interaction), "auc")@y.values[[1]]
  overall_performance[combination, "aucpr"] <- performance(prediction(predictions_test$prediction, predictions_test$interaction), "aucpr")@y.values[[1]]
  overall_performance[combination, "prevalence"] <- sum(predictions_test$interaction)/nrow(predictions_test)
  fw_sp <- unique(predictions$Predator)
  for (species in fw_sp){
    sp_predictions <- filter(predictions, Predator == species | Prey == species)
    if(sum(sp_predictions$interaction != 0)){
      sp_auc <- performance(prediction(sp_predictions$prediction, sp_predictions$interaction), "auc")@y.values[[1]]
      sp_aucpr <- performance(prediction(sp_predictions$prediction, sp_predictions$interaction), "aucpr")@y.values[[1]]
      species_performance <- rbind(species_performance,
                                   data.frame(sourceFW, targetFW, species, auc = sp_auc, aucpr = sp_aucpr, prevalence = sum(sp_predictions$interaction)/nrow(sp_predictions)))
    }
  }
}

write.csv(species_performance, "data/checkpoints/species_brt_performance.csv")
write.csv(overall_performance, "data/checkpoints/overall_brt_performance.csv")
