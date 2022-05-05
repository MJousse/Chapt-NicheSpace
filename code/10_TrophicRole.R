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
source("code/functions.R")
source("code/functions_motifs.R")

# Roles of all species in all empirical webs ------------------------------
# arctic
arcticFW <- read.csv("data/cleaned/HighArcticFW.csv", row.names = 1)
arcticFW <- arcticFW %>%
  transmute(resource = Prey, consumer = Predator)
arcticRoles <- species_role(arcticFW, ncores = 8)

# pyrenees
pyreneesFW <- read.csv("data/cleaned/pyrenneesFW.csv", row.names = 1)
pyreneesFW <- pyreneesFW %>%
  transmute(resource = Prey, consumer = Predator)
pyreneesRoles <- species_role(pyreneesFW, ncores = 8)

# serengeti
serengetiFW <- read.csv("data/cleaned/SerengetiFW.csv", row.names = 1)
serengetiFW <- serengetiFW %>%
  transmute(resource = Resource_Species, consumer = Consumer_Species) %>%
  drop_na()
serengetiRoles <- species_role(serengetiFW, ncores = 8)

# europe
europeMW <- read.csv("data/cleaned/EuroFWadults.csv", row.names = 1)
europeMW <- europeMW %>%
  transmute(resource = Prey, consumer = Predator)
europeRoles <- species_role(europeMW, ncores = 8)

# bind and save
empirical_roles <- bind_rows(
  arcticRoles %>% pivot_longer(-species, names_to = "role", values_to = "empirical") %>% mutate(FW = "Arctic"),
  pyreneesRoles %>% pivot_longer(-species, names_to = "role", values_to = "empirical") %>% mutate(FW = "Pyrenees"),
  serengetiRoles %>% pivot_longer(-species, names_to = "role", values_to = "empirical") %>% mutate(FW = "Serengeti"),
  europeRoles %>% pivot_longer(-species, names_to = "role", values_to = "empirical") %>% mutate(FW = "Europe")
)
write.csv(empirical_roles, file = "data/checkpoints/SpeciesRole.csv")


# Trophic roles in predicted webs -----------------------------------------
load("data/checkpoints/predictions.RData")
foodwebs <- c("Arctic", "Euro", "Pyrenees", "Serengeti")
combinations <- expand_grid(Source = foodwebs, Target = foodwebs)
predicted_roles <-c()

# for each combination of source and target webs, use 100 posterior sample
# calculate the roles of all species in these sample, and extract the mean
for (combination in c(1:nrow(combinations))){
  sourceFW <- combinations[combination, "Source"]
  targetFW <- combinations[combination, "Target"]
  predictions <- get(paste0(sourceFW, "_", targetFW, "_predictions")) %>%
    transmute(resource = Prey, consumer = Predator, interaction = Estimate)
  role <- c()
  for (i in c(1:100)){
    predictions$interaction <- rbinom(nrow(predictions), 1, predictions$interaction)
    predictions <- filter(predictions, interaction == 1)
    role <- rbind(role, species_role(predictions, ncores = 4))
  }
  role_mean <- group_by(role, species) %>% summarise_all(mean)
  role_sd <- group_by(role, species) %>% summarise_all(sd)
  predicted_roles <- rbind(predicted_roles, 
                           pivot_longer(role_mean, -species, names_to = "role", values_to = "predicted") %>% mutate(targetFW = targetFW, sourceFW = sourceFW))
}


# Compare the empirical roles and predicted roles -------------------------
species_roles <- left_join(predicted_roles, empirical_roles,
                           by = c("species", "role", "targetFW" = "FW"))

# calculate the error as the standardized mean difference
role_sd <- empirical_roles %>% group_by(role) %>% summarise(sd = sd(empirical))
species_roles$error = (species_roles$predicted - species_roles$empirical) / (role_sd$sd[match(species_roles$role, role_sd$role)])
species_roles$insample <- species_roles$targetFW == species_roles$sourceFW

# plot
ggplot(species_roles, aes(x = role, y = error, fill = insample)) +
  geom_point() +
  geom_boxplot()
