library(ROCR)
library(dplyr)
library(tidyr)
rm(list = ls())
load("data/checkpoints/predictions.RData")

foodwebs <- c("Arctic", "Euro", "Pyrenees", "Serengeti")
overall_performance <- expand_grid(Source = foodwebs, Target = foodwebs) %>%
  mutate(auc = NA, aucpr = NA)
species_performance <- data.frame()

for (combination in c(1:nrow(overall_performance))){
  sourceFW <- overall_performance[combination, "Source"]
  targetFW <- overall_performance[combination, "Target"]
  predictions <- get(paste0(sourceFW, "_", targetFW, "_predictions"))
  overall_performance[combination, "auc"] <- performance(prediction(predictions$Estimate, predictions$interaction), "auc")@y.values[[1]]
  overall_performance[combination, "aucpr"] <- performance(prediction(predictions$Estimate, predictions$interaction), "aucpr")@y.values[[1]]
  fw_sp <- unique(predictions$Predator)
  for (species in fw_sp){
   sp_predictions <- filter(predictions, Predator == species | Prey == species)
   if(sum(sp_predictions$interaction != 0)){
     sp_auc <- performance(prediction(sp_predictions$Estimate, sp_predictions$interaction), "auc")@y.values[[1]]
     sp_aucpr <- performance(prediction(sp_predictions$Estimate, sp_predictions$interaction), "aucpr")@y.values[[1]]
     species_performance <- rbind(species_performance,
                                  data.frame(sourceFW, targetFW, species, auc = sp_auc, aucpr = sp_aucpr))
   }
  }
}

ggplot(species_performance) +
  geom_jitter(aes(x = Source, y = auc, colour = Source), alpha = 0.25) +
  geom_point(data = overall_performance, aes(x = Source, y = auc, fill = Source), shape = 21, size = 3) +
  facet_wrap(~Target, nrow = 1)


