# Step 11: Analyze the biases in food web properties predictions
# 1. Calculate the trophic roles of all species in all empirical webs
# Roles are related to centrality, trophic level, module-based and motifs-based.
# 2. Calculate the trophic roles of all species in predicted webs
# For each webs, use 100 posterior samples and use the mean.
# 3. Calulate the errors for each role metrics of all species
set.seed(16)
library(igraph)
library(dplyr)
library(tidyr)
library(NetIndices)
library(ggplot2)
library(patchwork)
source("code/functions.R")


# FW properties in predicted webs -----------------------------------------
library(foreach)
library(doParallel)
load("data/checkpoints/BRT_predictions.RData")
foodwebs <- c("arctic", "pyrenees", "serengeti", "europe")
combinations <- expand_grid(Source = foodwebs, Target = foodwebs)
predicted_properties <-c()
# for each combination of source and target webs, use 100 posterior sample
# calculate the roles of all species in these sample, and extract the mean
for (combination in c(1:nrow(combinations))){
  sourceFW <- as.character(combinations[combination, "Source"])
  targetFW <- as.character(combinations[combination, "Target"])
  print(Sys.time())
  print(paste0(sourceFW, " model predicting the ", targetFW, " food web..."))
  predictions <- get(paste0(sourceFW, "BRT_", targetFW, "FW"))
  cl <- makeCluster(16) 
  registerDoParallel(cl)
  properties <- foreach(i=c(1:100), .combine = rbind, 
                        .packages = c("igraph", "NetIndices", "dplyr", "tidyr")) %dopar% {
                          prediction <- data.frame(resource = predictions$Prey, 
                                                   consumer = predictions$Predator, 
                                                   interaction = rbinom(nrow(predictions),2,predictions$prediction))
                          prediction <- prediction[prediction$interaction == 1,]
                          fw_properties(prediction, nsim = 1)
                        }
  stopCluster(cl)
  fw_properties_mean <- apply(properties, MARGIN = 2, mean)
  fw_properties_sd <- apply(properties, MARGIN = 2, sd)
  predicted_properties <- rbind(predicted_properties, 
                                data.frame(t(fw_properties_mean)) %>% pivot_longer(everything(), names_to = "metric", values_to = "predicted") %>% mutate(targetFW = targetFW, sourceFW = sourceFW))
  write.csv(predicted_properties, file = "data/checkpoints/PredictedProperties_BRT.csv")
}
