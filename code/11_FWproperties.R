# Step 11: Analyze the biases in food web properties predictions
# 1. Calculate the trophic roles of all species in all empirical webs
# Roles are related to centrality, trophic level, module-based and motifs-based.
# 2. Calculate the trophic roles of all species in predicted webs
# For each webs, use 100 posterior samples and use the mean.
# 3. Calulate the errors for each role metrics of all species

library(igraph)
library(dplyr)
library(tidyr)
library(NetIndices)
library(ggplot2)
source("code/functions.R")

# Roles of all species in all empirical webs ------------------------------
# arctic
arcticFW <- read.csv("data/cleaned/HighArcticFW.csv", row.names = 1)
arcticFW <- arcticFW %>%
  transmute(resource = Prey, consumer = Predator)
arcticProperties <- fw_properties(arcticFW, nsim = 10)

# pyrenees
pyreneesFW <- read.csv("data/cleaned/pyrenneesFW.csv", row.names = 1)
pyreneesFW <- pyreneesFW %>%
  transmute(resource = Prey, consumer = Predator)
pyreneesProperties <- fw_properties(pyreneesFW, nsim = 10)

# serengeti
serengetiFW <- read.csv("data/cleaned/SerengetiFW.csv", row.names = 1)
serengetiFW <- serengetiFW %>%
  transmute(resource = Resource_Species, consumer = Consumer_Species) %>%
  drop_na()
serengetiProperties <- fw_properties(serengetiFW, nsim = 10)

# europe
europeMW <- read.csv("data/cleaned/EuroFWadults.csv", row.names = 1)
europeMW <- europeMW %>%
  transmute(resource = Prey, consumer = Predator)
europeProperties <- fw_properties(europeMW, nsim = 10)

# bind and save
empirical_properties <- bind_rows(
  data.frame(t(arcticProperties)) %>% pivot_longer(everything(), names_to = "metric", values_to = "empirical") %>% mutate(FW = "Arctic"),
  data.frame(t(pyreneesProperties)) %>% pivot_longer(everything(), names_to = "metric", values_to = "empirical") %>% mutate(FW = "Pyrenees"),
  data.frame(t(serengetiProperties)) %>% pivot_longer(everything(), names_to = "metric", values_to = "empirical") %>% mutate(FW = "Serengeti"),
  data.frame(t(europeProperties)) %>% pivot_longer(everything(), names_to = "metric", values_to = "empirical") %>% mutate(FW = "Europe")
)
write.csv(empirical_properties, file = "data/checkpoints/EmpiricalProperties.csv")

# Trophic roles in predicted webs -----------------------------------------
library(foreach)
library(doParallel)
load("data/checkpoints/predictions.RData")
foodwebs <- c("Arctic", "Euro", "Pyrenees", "Serengeti")
combinations <- expand_grid(Source = foodwebs, Target = foodwebs)
predicted_properties <-c()
# for each combination of source and target webs, use 100 posterior sample
# calculate the roles of all species in these sample, and extract the mean
for (combination in c(1:nrow(combinations))){
  sourceFW <- as.character(combinations[combination, "Source"])
  targetFW <- as.character(combinations[combination, "Target"])
  print(Sys.time())
  print(paste0(sourceFW, " model predicting the ", targetFW, " food web..."))
  predictions <- get(paste0(sourceFW, "_", targetFW, "_predictions"))
  cl <- makeCluster(8) 
  registerDoParallel(cl)
  properties <- foreach(i=c(1:100), .combine = rbind, 
                  .packages = c("igraph", "NetIndices", "dplyr", "tidyr")) %dopar% {
                    prediction <- data.frame(resource = predictions$Prey, 
                                             consumer = predictions$Predator, 
                                             interaction = predictions[,paste0("draws",i)])
                    prediction <- prediction[prediction$interaction == 1,]
                    fw_properties(prediction, nsim = 10)
                  }
  stopCluster(cl)
  fw_properties_mean <- apply(properties, MARGIN = 2, mean)
  fw_properties_sd <- apply(properties, MARGIN = 2, sd)
  predicted_properties <- rbind(predicted_properties, 
                                data.frame(t(fw_properties_mean)) %>% pivot_longer(everything(), names_to = "metric", values_to = "predicted") %>% mutate(targetFW = targetFW, sourceFW = sourceFW))
  write.csv(predicted_roles, file = "data/checkpoints/predicted_roles.csv")
}

# Compare the empirical roles and predicted roles -------------------------
empirical_roles <- read.csv("data/checkpoints/SpeciesRole.csv", row.names = 1)
predicted_roles <- read.csv("data/checkpoints/predicted_roles.csv", row.names = 1)

species_roles <- left_join(predicted_roles, empirical_roles,
                           by = c("species", "role", "targetFW" = "FW")) %>%
  drop_na() %>%
  group_by(role, targetFW, sourceFW) %>%
  mutate(predicted_scaled = (predicted - mean(empirical, na.rm = T))/sd(empirical, na.rm = T),
         empirical_scaled = (empirical - mean(empirical, na.rm = T))/sd(empirical, na.rm = T))
