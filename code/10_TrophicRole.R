library(igraph)
library(dplyr)
library(tidyr)
library(NetIndices)
library(multiweb)
source("code/functions.R")
source("code/functions_motifs.R")

# roles in empirical food web
arcticFW <- read.csv("data/cleaned/HighArcticFW.csv", row.names = 1)
arcticFW <- arcticFW %>%
  transmute(resource = Prey, consumer = Predator)
arcticRoles <- species_role(arcticFW, ncores = 8)

pyreneesFW <- read.csv("data/cleaned/pyrenneesFW.csv", row.names = 1)
pyreneesFW <- pyreneesFW %>%
  transmute(resource = Prey, consumer = Predator)
pyreneesRoles <- species_role(pyreneesFW, ncores = 8)

serengetiFW <- read.csv("data/cleaned/SerengetiFW.csv", row.names = 1)
serengetiFW <- serengetiFW %>%
  transmute(resource = Resource_Species, consumer = Consumer_Species) %>%
  drop_na()
serengetiRoles <- species_role(serengetiFW, ncores = 8)

EuropeMW <- read.csv("data/cleaned/EuroFWadults.csv", row.names = 1)
EuropeMW <- EuropeMW %>%
  transmute(resource = Prey, consumer = Predator)
europeRoles <- species_role(EuropeMW, ncores = 8)

empirical_roles <- bind_rows(
  arcticRoles %>% pivot_longer(-species, names_to = "role", values_to = "empirical") %>% mutate(FW = "Arctic"),
  pyreneesRoles %>% pivot_longer(-species, names_to = "role", values_to = "empirical") %>% mutate(FW = "Pyrenees"),
  serengetiRoles %>% pivot_longer(-species, names_to = "role", values_to = "empirical") %>% mutate(FW = "Serengeti"),
  europeRoles %>% pivot_longer(-species, names_to = "role", values_to = "empirical") %>% mutate(FW = "Europe")
)

write.csv(empirical_roles, file = "data/checkpoints/SpeciesRole.csv")

# roles predicted
load("data/checkpoints/predictions.RData")
foodwebs <- c("Arctic", "Euro", "Pyrenees", "Serengeti")
combinations <- expand_grid(Source = foodwebs, Target = foodwebs)
predicted_roles <-c()

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

species_roles <- left_join(predicted_roles, empirical_roles,
                           by = c("species", "role", "targetFW" = "FW"))
role_sd <- empirical_roles %>% group_by(role) %>% summarise(sd = sd(empirical))

species_roles$error = (species_roles$predicted - species_roles$empirical) / (role_sd$sd[match(species_roles$role, role_sd$role)])
species_roles$insample <- species_roles$targetFW == species_roles$sourceFW

ggplot(species_roles, aes(x = role, y = error, fill = insample)) +
  geom_point() +
  geom_boxplot()
