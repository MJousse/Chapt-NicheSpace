# Step 10: Analyze the biases in trophic role predictions
# 1. Calculate the trophic roles of all species in all empirical webs
# Roles are related to centrality, trophic level, module-based and motifs-based.
# 2. Calculate the trophic roles of all species in predicted webs
# For each webs, use 100 posterior samples and use the mean.
# 3. Calulate the errors for each role metrics of all species

library(igraph)
library(dplyr)
library(tidyr)
library(NetIndices)
library(multiweb)
library(ggplot2)
source("code/functions.R")
source("code/functions_motifs.R")

# Roles of all species in all empirical webs ------------------------------
# arctic
arcticFW <- read.csv("data/cleaned/HighArcticFW.csv", row.names = 1)
arcticFW <- arcticFW %>%
  transmute(resource = Prey, consumer = Predator)
arcticRoles <- species_role(arcticFW, ncores = 16)

# pyrenees
pyreneesFW <- read.csv("data/cleaned/pyrenneesFW.csv", row.names = 1)
pyreneesFW <- pyreneesFW %>%
  transmute(resource = Prey, consumer = Predator)
pyreneesRoles <- species_role(pyreneesFW, ncores = 16)

# serengeti
serengetiFW <- read.csv("data/cleaned/SerengetiFW.csv", row.names = 1)
serengetiFW <- serengetiFW %>%
  transmute(resource = Resource_Species, consumer = Consumer_Species) %>%
  drop_na()
serengetiRoles <- species_role(serengetiFW, ncores = 16)

# europe
europeMW <- read.csv("data/cleaned/EuroFWadults.csv", row.names = 1)
europeMW <- europeMW %>%
  transmute(resource = Prey, consumer = Predator)
europeRoles <- species_role(europeMW, ncores = 16)

# bind and save
empirical_roles <- bind_rows(
  arcticRoles %>% pivot_longer(-species, names_to = "role", values_to = "empirical") %>% mutate(FW = "Arctic"),
  pyreneesRoles %>% pivot_longer(-species, names_to = "role", values_to = "empirical") %>% mutate(FW = "Pyrenees"),
  serengetiRoles %>% pivot_longer(-species, names_to = "role", values_to = "empirical") %>% mutate(FW = "Serengeti"),
  europeRoles %>% pivot_longer(-species, names_to = "role", values_to = "empirical") %>% mutate(FW = "Europe")
)
write.csv(empirical_roles, file = "data/checkpoints/SpeciesRole.csv")

# Trophic roles in predicted webs -----------------------------------------
library(foreach)
library(doParallel)
load("data/checkpoints/predictions.RData")
foodwebs <- c("Arctic", "Euro", "Pyrenees", "Serengeti")
combinations <- expand_grid(Source = foodwebs, Target = foodwebs) %>%
  filter(Target != "Euro")
predicted_roles <-c()
# for each combination of source and target webs, use 100 posterior sample
# calculate the roles of all species in these sample, and extract the mean
for (combination in c(1:nrow(combinations))){
  sourceFW <- as.character(combinations[combination, "Source"])
  targetFW <- as.character(combinations[combination, "Target"])
  print(Sys.time())
  print(paste0(sourceFW, " model predicting the ", targetFW, " food web..."))
  predictions <- get(paste0(sourceFW, "_", targetFW, "_predictions"))
  cl <- makeCluster(10) 
  registerDoParallel(cl)
  role <- foreach(i=c(1:100), .combine = rbind, 
                  .packages = c("igraph", "NetIndices", "multiweb", "dplyr", "tidyr"),
                  .export = c("motif_role", "positions")) %dopar% {
    prediction <- data.frame(resource = predictions$Prey, 
                              consumer = predictions$Predator, 
                              interaction = predictions[,paste0("draws",i)])
    prediction <- prediction[prediction$interaction == 1,]
    species_role(prediction, ncores = 0)
                  }
  stopCluster(cl)
  role_mean <- group_by(role, species) %>% summarise_all(mean, na.rm = T)
  role_sd <- group_by(role, species) %>% summarise_all(sd, na.rm = T)
  predicted_roles <- rbind(predicted_roles, 
                           pivot_longer(role_mean, -species, names_to = "role", values_to = "predicted") %>% mutate(targetFW = targetFW, sourceFW = sourceFW))
  write.csv(predicted_roles, file = "data/checkpoints/predicted_roles.csv")
}

# Compare the empirical roles and predicted roles -------------------------
species_roles <- left_join(predicted_roles, empirical_roles,
                           by = c("species", "role", "targetFW" = "FW"))

# calculate the error as the standardized mean difference
species_roles$error = (species_roles$predicted - species_roles$empirical)
sd_role_error <- species_roles %>% group_by(role) %>% summarise(sd = sd(error))
species_roles$error_std = species_roles$error / sd_role_error$sd[match(species_roles$role, sd_role_error$role)]
species_roles$insample <- species_roles$targetFW == species_roles$sourceFW

# plot
species_roles$role <- factor(species_roles$role, levels = unique(species_roles$role))
ggplot(species_roles, aes(x = role, y = error_std, fill = insample)) +
  geom_point() +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
