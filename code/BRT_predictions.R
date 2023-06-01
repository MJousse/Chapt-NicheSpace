
# Predictions
rm(list = ls())
library(dismo)
library(gbm)
library(tidyr)
library(dplyr)
source("code/functions.R")
load("../../../OneDrive/Chapt-NicheSpace/models/brt_models.RData")
load("data/checkpoints/train_test_splits.RData")

# functional traits
FuncTraits <- read.csv("data/cleaned/SpeciesTraitsFull.csv", row.names = 1)

# european metaweb
EuroInteractions <- read.csv("data/cleaned/EuroFWadults.csv", row.names = 1)
EuroSpecies <- read.csv("data/cleaned/EuroMWTaxo.csv", row.names = 1) %>%
  distinct() %>%
  filter(!is.na(Species))
EuroMW <- get_predictors(EuroSpecies$Species, FuncTraits)
EuroInteractions$interaction <- 1
EuroMW <- left_join(EuroMW, EuroInteractions)
EuroMW$interaction[is.na(EuroMW$interaction)] <- 0
EuroMW$TL.predator <- ifelse(EuroMW$Herbivore.predator == 1, "Herbivore", ifelse(EuroMW$Omnivore.predator == 1, "Omnivore", "Carnivore"))
EuroMW$TL.prey <- ifelse(EuroMW$Herbivore.prey == 1, "Herbivore", ifelse(EuroMW$Omnivore.prey == 1, "Omnivore", "Carnivore"))

EuroMW <- select(EuroMW, -Order.prey, -Herbivore.predator, 
                 -Omnivore.predator, -Carnivore.predator, -Herbivore.prey, 
                 -Omnivore.prey, -Carnivore.prey) %>%
  mutate(Order.predator = factor(Order.predator, 
                                 levels = unique(FuncTraits$Order)),
         TL.predator = factor(TL.predator, levels = c("Herbivore", "Omnivore", "Carnivore")),
         TL.prey = factor(TL.prey, levels = c("Herbivore", "Omnivore", "Carnivore")),
         ActivityTime.match = factor(ActivityTime.match))


# high arctic
HighArcticInteractions <- read.csv("data/cleaned/HighArcticFW.csv", row.names = 1)
HighArcticSpecies <- read.csv("data/cleaned/HighArcticFWTaxo.csv", row.names = 1) %>%
  distinct() %>%
  filter(!is.na(Species))
HighArcticFW <- get_predictors(HighArcticSpecies$Species, FuncTraits)
HighArcticInteractions$interaction <- 1
HighArcticFW <- left_join(HighArcticFW, HighArcticInteractions)
HighArcticFW$interaction[is.na(HighArcticFW$interaction)] <- 0
HighArcticFW$Order.predator <- factor(HighArcticFW$Order.predator, 
                                      levels = unique(FuncTraits$Order))
HighArcticFW$TL.predator <- ifelse(HighArcticFW$Herbivore.predator == 1, "Herbivore", ifelse(HighArcticFW$Omnivore.predator == 1, "Omnivore", "Carnivore"))
HighArcticFW$TL.prey <- ifelse(HighArcticFW$Herbivore.prey == 1, "Herbivore", ifelse(HighArcticFW$Omnivore.prey == 1, "Omnivore", "Carnivore"))

HighArcticFW <- select(HighArcticFW, -Order.prey, -Herbivore.predator, 
                       -Omnivore.predator, -Carnivore.predator, -Herbivore.prey, 
                       -Omnivore.prey, -Carnivore.prey) %>%
  mutate(Order.predator = factor(Order.predator, 
                                 levels = unique(FuncTraits$Order)),
         TL.predator = factor(TL.predator, levels = c("Herbivore", "Omnivore", "Carnivore")),
         TL.prey = factor(TL.prey, levels = c("Herbivore", "Omnivore", "Carnivore")),
         ActivityTime.match = factor(ActivityTime.match))

# pyrenees
PyreneesInteractions <- read.csv("data/cleaned/pyrenneesFW.csv", row.names = 1)
PyreneesSpecies <- read.csv("data/cleaned/pyrenneesFWTaxo.csv", row.names = 1) %>%
  distinct() %>%
  filter(!is.na(Species))
PyreneesFW <- get_predictors(PyreneesSpecies$Species, FuncTraits)
PyreneesInteractions$interaction <- 1
PyreneesFW <- left_join(PyreneesFW, PyreneesInteractions)
PyreneesFW$interaction[is.na(PyreneesFW$interaction)] <- 0
PyreneesFW$Order.predator <- factor(PyreneesFW$Order.predator, 
                                    levels = unique(FuncTraits$Order))
PyreneesFW$TL.predator <- ifelse(PyreneesFW$Herbivore.predator == 1, "Herbivore", ifelse(PyreneesFW$Omnivore.predator == 1, "Omnivore", "Carnivore"))
PyreneesFW$TL.prey <- ifelse(PyreneesFW$Herbivore.prey == 1, "Herbivore", ifelse(PyreneesFW$Omnivore.prey == 1, "Omnivore", "Carnivore"))

PyreneesFW <- select(PyreneesFW, -Order.prey, -Herbivore.predator, 
                     -Omnivore.predator, -Carnivore.predator, -Herbivore.prey, 
                     -Omnivore.prey, -Carnivore.prey) %>%
  mutate(Order.predator = factor(Order.predator, 
                                 levels = unique(FuncTraits$Order)),
         TL.predator = factor(TL.predator, levels = c("Herbivore", "Omnivore", "Carnivore")),
         TL.prey = factor(TL.prey, levels = c("Herbivore", "Omnivore", "Carnivore")),
         ActivityTime.match = factor(ActivityTime.match))

# serengeti
SerengetiInteractions <- read.csv("data/cleaned/SerengetiFW.csv", row.names = 1) %>%
  select(Predator = Consumer_Species, Prey = Resource_Species)
SerengetiSpecies <- read.csv("data/cleaned/SerengetiFWTaxo.csv", row.names = 1) %>%
  distinct() %>%
  filter(!is.na(Species))
SerengetiFW <- get_predictors(SerengetiSpecies$Species, FuncTraits)
SerengetiInteractions$interaction <- 1
SerengetiFW <- left_join(SerengetiFW, SerengetiInteractions)
SerengetiFW$interaction[is.na(SerengetiFW$interaction)] <- 0
SerengetiFW$Order.predator <- factor(SerengetiFW$Order.predator, 
                                     levels = unique(FuncTraits$Order))
# prepare training set
SerengetiFW$TL.predator <- ifelse(SerengetiFW$Herbivore.predator == 1, "Herbivore", ifelse(SerengetiFW$Omnivore.predator == 1, "Omnivore", "Carnivore"))
SerengetiFW$TL.prey <- ifelse(SerengetiFW$Herbivore.prey == 1, "Herbivore", ifelse(SerengetiFW$Omnivore.prey == 1, "Omnivore", "Carnivore"))

SerengetiFW <- select(SerengetiFW, -Order.prey, -Herbivore.predator, 
                      -Omnivore.predator, -Carnivore.predator, -Herbivore.prey, 
                      -Omnivore.prey, -Carnivore.prey) %>%
  mutate(Order.predator = factor(Order.predator, 
                                 levels = unique(FuncTraits$Order)),
         TL.predator = factor(TL.predator, levels = c("Herbivore", "Omnivore", "Carnivore")),
         TL.prey = factor(TL.prey, levels = c("Herbivore", "Omnivore", "Carnivore")),
         ActivityTime.match = factor(ActivityTime.match))

# Predictions using the Arctic model --------------------------------------
# find the mean and sd of the predictors to scale new data
var2scale <- c("Habitat_breadth.predator", "BM.predator", "Longevity.predator", 
               "ClutchSize.predator", "Habitat_breadth.prey", "BM.prey",
               "Longevity.prey", "ClutchSize.prey", "Habitat.match", "BM.match")
Arctic_mean <- apply(HighArcticFW[,var2scale], MARGIN = 2, mean)
Arctic_sd <- apply(HighArcticFW[,var2scale], MARGIN = 2, sd)

# scale predictors
HighArcticFWscaled <- HighArcticFW
HighArcticFWscaled[,var2scale] <- sweep(HighArcticFWscaled[,var2scale], MARGIN = 2, Arctic_mean)
HighArcticFWscaled[,var2scale] <- sweep(HighArcticFWscaled[,var2scale], MARGIN = 2, FUN = "/", 2*Arctic_sd)
EuroMWscaled <- EuroMW
EuroMWscaled[,var2scale] <- sweep(EuroMWscaled[,var2scale], MARGIN = 2, Arctic_mean)
EuroMWscaled[,var2scale] <- sweep(EuroMWscaled[,var2scale], MARGIN = 2, FUN = "/", 2*Arctic_sd)
PyreneesFWscaled <- PyreneesFW
PyreneesFWscaled[,var2scale] <- sweep(PyreneesFWscaled[,var2scale], MARGIN = 2, Arctic_mean)
PyreneesFWscaled[,var2scale] <- sweep(PyreneesFWscaled[,var2scale], MARGIN = 2, FUN = "/", 2*Arctic_sd)
SerengetiFWscaled <- SerengetiFW
SerengetiFWscaled[,var2scale] <- sweep(SerengetiFWscaled[,var2scale], MARGIN = 2, Arctic_mean)
SerengetiFWscaled[,var2scale] <- sweep(SerengetiFWscaled[,var2scale], MARGIN = 2, FUN = "/", 2*Arctic_sd)

# predict arctic food webs
arcticBRT_arcticFW <- data.frame(Predator= HighArcticFWscaled$Predator, 
                                        Prey = HighArcticFWscaled$Prey,
                                        interaction = HighArcticFWscaled$interaction)
  
arcticBRT_arcticFW$prediction <- predict.gbm(brt_arctic, HighArcticFWscaled,
                                         n.trees=brt_arctic$gbm.call$best.trees, type="response")
arcticBRT_arcticFW$testing <- 0
arcticBRT_arcticFW$testing[testing_id_higharctic] <- 1

# predict european metaweb
arcticBRT_europeFW <- data.frame(Predator= EuroMWscaled$Predator, 
                                        Prey = EuroMWscaled$Prey,
                                        interaction = EuroMWscaled$interaction)

arcticBRT_europeFW$prediction <- predict.gbm(brt_arctic, EuroMWscaled,
                                                    n.trees=brt_arctic$gbm.call$best.trees, type="response")

# predict pyrenees food web
arcticBRT_pyreneesFW <- data.frame(Predator= PyreneesFWscaled$Predator, 
                                      Prey = PyreneesFWscaled$Prey,
                                      interaction = PyreneesFWscaled$interaction)

arcticBRT_pyreneesFW$prediction <- predict.gbm(brt_arctic, PyreneesFWscaled,
                                                  n.trees=brt_arctic$gbm.call$best.trees, type="response")

# predict serengeti food web
arcticBRT_serengetiFW <- data.frame(Predator= SerengetiFWscaled$Predator, 
                                          Prey = SerengetiFWscaled$Prey,
                                          interaction = SerengetiFWscaled$interaction)

arcticBRT_serengetiFW$prediction <- predict.gbm(brt_arctic, SerengetiFWscaled,
                                                      n.trees=brt_arctic$gbm.call$best.trees, type="response")

# Predictions using the European model --------------------------------------
# find the mean and sd of the predictors to scale new data
Europe_mean <- apply(EuroMW[,var2scale], MARGIN = 2, mean)
Europe_sd <- apply(EuroMW[,var2scale], MARGIN = 2, sd)

# scale predictors
HighArcticFWscaled <- HighArcticFW
HighArcticFWscaled[,var2scale] <- sweep(HighArcticFWscaled[,var2scale], MARGIN = 2, Europe_mean)
HighArcticFWscaled[,var2scale] <- sweep(HighArcticFWscaled[,var2scale], MARGIN = 2, FUN = "/", 2*Europe_sd)
EuroMWscaled <- EuroMW
EuroMWscaled[,var2scale] <- sweep(EuroMWscaled[,var2scale], MARGIN = 2, Europe_mean)
EuroMWscaled[,var2scale] <- sweep(EuroMWscaled[,var2scale], MARGIN = 2, FUN = "/", 2*Europe_sd)
PyreneesFWscaled <- PyreneesFW
PyreneesFWscaled[,var2scale] <- sweep(PyreneesFWscaled[,var2scale], MARGIN = 2, Europe_mean)
PyreneesFWscaled[,var2scale] <- sweep(PyreneesFWscaled[,var2scale], MARGIN = 2, FUN = "/", 2*Europe_sd)
SerengetiFWscaled <- SerengetiFW
SerengetiFWscaled[,var2scale] <- sweep(SerengetiFWscaled[,var2scale], MARGIN = 2, Europe_mean)
SerengetiFWscaled[,var2scale] <- sweep(SerengetiFWscaled[,var2scale], MARGIN = 2, FUN = "/", 2*Europe_sd)

# predict arctic food webs
europeBRT_arcticFW <- data.frame(Predator= HighArcticFWscaled$Predator, 
                                          Prey = HighArcticFWscaled$Prey,
                                          interaction = HighArcticFWscaled$interaction)

europeBRT_arcticFW$prediction <- predict.gbm(brt_europe, HighArcticFWscaled,
                                                      n.trees=brt_europe$gbm.call$best.trees, type="response")

# predict european metaweb
europeBRT_europeFW <- data.frame(Predator= EuroMWscaled$Predator, 
                                        Prey = EuroMWscaled$Prey,
                                        interaction = EuroMWscaled$interaction)

europeBRT_europeFW$prediction <- predict.gbm(brt_europe, EuroMWscaled,
                                                    n.trees=brt_europe$gbm.call$best.trees, type="response")
europeBRT_europeFW$testing <- 0
europeBRT_europeFW$testing[testing_id_euro] <- 1

# predict pyrenees food web
europeBRT_pyreneesFW <- data.frame(Predator= PyreneesFWscaled$Predator, 
                                            Prey = PyreneesFWscaled$Prey,
                                            interaction = PyreneesFWscaled$interaction)

europeBRT_pyreneesFW$prediction <- predict.gbm(brt_europe, PyreneesFWscaled,
                                                        n.trees=brt_europe$gbm.call$best.trees, type="response")

# predict serengeti food web
europeBRT_serengetiFW <- data.frame(Predator= SerengetiFWscaled$Predator, 
                                             Prey = SerengetiFWscaled$Prey,
                                             interaction = SerengetiFWscaled$interaction)

europeBRT_serengetiFW$prediction <- predict.gbm(brt_europe, SerengetiFWscaled,
                                                         n.trees=brt_europe$gbm.call$best.trees, type="response")

# Predictions using the Pyrenees model --------------------------------------
# find the mean and sd of the predictors to scale new data
Pyrenees_mean <- apply(PyreneesFW[,var2scale], MARGIN = 2, mean)
Pyrenees_sd <- apply(PyreneesFW[,var2scale], MARGIN = 2, sd)

# scale predictors
HighArcticFWscaled <- HighArcticFW
HighArcticFWscaled[,var2scale] <- sweep(HighArcticFWscaled[,var2scale], MARGIN = 2, Pyrenees_mean)
HighArcticFWscaled[,var2scale] <- sweep(HighArcticFWscaled[,var2scale], MARGIN = 2, FUN = "/", 2*Pyrenees_sd)
EuroMWscaled <- EuroMW
EuroMWscaled[,var2scale] <- sweep(EuroMWscaled[,var2scale], MARGIN = 2, Pyrenees_mean)
EuroMWscaled[,var2scale] <- sweep(EuroMWscaled[,var2scale], MARGIN = 2, FUN = "/", 2*Pyrenees_sd)
PyreneesFWscaled <- PyreneesFW
PyreneesFWscaled[,var2scale] <- sweep(PyreneesFWscaled[,var2scale], MARGIN = 2, Pyrenees_mean)
PyreneesFWscaled[,var2scale] <- sweep(PyreneesFWscaled[,var2scale], MARGIN = 2, FUN = "/", 2*Pyrenees_sd)
SerengetiFWscaled <- SerengetiFW
SerengetiFWscaled[,var2scale] <- sweep(SerengetiFWscaled[,var2scale], MARGIN = 2, Pyrenees_mean)
SerengetiFWscaled[,var2scale] <- sweep(SerengetiFWscaled[,var2scale], MARGIN = 2, FUN = "/", 2*Pyrenees_sd)

# predict arctic food webs
pyreneesBRT_arcticFW <- data.frame(Predator= HighArcticFWscaled$Predator, 
                                        Prey = HighArcticFWscaled$Prey,
                                        interaction = HighArcticFWscaled$interaction)

pyreneesBRT_arcticFW$prediction <- predict.gbm(brt_pyrenees, HighArcticFWscaled,
                                                    n.trees=brt_pyrenees$gbm.call$best.trees, type="response")

# predict european metaweb
pyreneesBRT_europeFW <- data.frame(Predator= EuroMWscaled$Predator, 
                                      Prey = EuroMWscaled$Prey,
                                      interaction = EuroMWscaled$interaction)

pyreneesBRT_europeFW$prediction <- predict.gbm(brt_pyrenees, EuroMWscaled,
                                                  n.trees=brt_pyrenees$gbm.call$best.trees, type="response")

# predict pyrenees food web
pyreneesBRT_pyreneesFW <- data.frame(Predator= PyreneesFWscaled$Predator, 
                                          Prey = PyreneesFWscaled$Prey,
                                          interaction = PyreneesFWscaled$interaction)

pyreneesBRT_pyreneesFW$prediction <- predict.gbm(brt_pyrenees, PyreneesFWscaled,
                                                      n.trees=brt_pyrenees$gbm.call$best.trees, type="response")
pyreneesBRT_pyreneesFW$testing <- 0
pyreneesBRT_pyreneesFW$testing[testing_id_pyrenees] <- 1

# predict serengeti food web
pyreneesBRT_serengetiFW <- data.frame(Predator= SerengetiFWscaled$Predator, 
                                           Prey = SerengetiFWscaled$Prey,
                                           interaction = SerengetiFWscaled$interaction)

pyreneesBRT_serengetiFW$prediction <- predict.gbm(brt_pyrenees, SerengetiFWscaled,
                                                       n.trees=brt_pyrenees$gbm.call$best.trees, type="response")

# Predictions using the Serengeti model --------------------------------------
# find the mean and sd of the predictors to scale new data
Serengeti_mean <- apply(EuroMW[,var2scale], MARGIN = 2, mean)
Serengeti_sd <- apply(EuroMW[,var2scale], MARGIN = 2, sd)

# scale predictors
HighArcticFWscaled <- HighArcticFW
HighArcticFWscaled[,var2scale] <- sweep(HighArcticFWscaled[,var2scale], MARGIN = 2, Serengeti_mean)
HighArcticFWscaled[,var2scale] <- sweep(HighArcticFWscaled[,var2scale], MARGIN = 2, FUN = "/", 2*Serengeti_sd)
EuroMWscaled <- EuroMW
EuroMWscaled[,var2scale] <- sweep(EuroMWscaled[,var2scale], MARGIN = 2, Serengeti_mean)
EuroMWscaled[,var2scale] <- sweep(EuroMWscaled[,var2scale], MARGIN = 2, FUN = "/", 2*Serengeti_sd)
PyreneesFWscaled <- PyreneesFW
PyreneesFWscaled[,var2scale] <- sweep(PyreneesFWscaled[,var2scale], MARGIN = 2, Serengeti_mean)
PyreneesFWscaled[,var2scale] <- sweep(PyreneesFWscaled[,var2scale], MARGIN = 2, FUN = "/", 2*Serengeti_sd)
SerengetiFWscaled <- SerengetiFW
SerengetiFWscaled[,var2scale] <- sweep(SerengetiFWscaled[,var2scale], MARGIN = 2, Serengeti_mean)
SerengetiFWscaled[,var2scale] <- sweep(SerengetiFWscaled[,var2scale], MARGIN = 2, FUN = "/", 2*Serengeti_sd)

# predict arctic food webs
serengetiBRT_arcticFW <- data.frame(Predator= HighArcticFWscaled$Predator, 
                                            Prey = HighArcticFWscaled$Prey,
                                            interaction = HighArcticFWscaled$interaction)

serengetiBRT_arcticFW$prediction <- predict.gbm(brt_serengeti, HighArcticFWscaled,
                                                        n.trees=brt_serengeti$gbm.call$best.trees, type="response")

# predict european metaweb
serengetiBRT_europeFW <- data.frame(Predator= EuroMWscaled$Predator, 
                                          Prey = EuroMWscaled$Prey,
                                          interaction = EuroMWscaled$interaction)

serengetiBRT_europeFW$prediction <- predict.gbm(brt_serengeti, EuroMWscaled,
                                                      n.trees=brt_serengeti$gbm.call$best.trees, type="response")

# predict pyrenees food web
serengetiBRT_pyreneesFW <- data.frame(Predator= PyreneesFWscaled$Predator, 
                                              Prey = PyreneesFWscaled$Prey,
                                              interaction = PyreneesFWscaled$interaction)

serengetiBRT_pyreneesFW$prediction <- predict.gbm(brt_serengeti, PyreneesFWscaled,
                                                          n.trees=brt_serengeti$gbm.call$best.trees, type="response")

# predict serengeti food web
serengetiBRT_serengetiFW <- data.frame(Predator= SerengetiFWscaled$Predator, 
                                               Prey = SerengetiFWscaled$Prey,
                                               interaction = SerengetiFWscaled$interaction)

serengetiBRT_serengetiFW$prediction <- predict.gbm(brt_serengeti, SerengetiFWscaled,
                                                           n.trees=brt_serengeti$gbm.call$best.trees, type="response")
serengetiBRT_serengetiFW$testing <- 0
serengetiBRT_serengetiFW$testing[testing_id_serengeti] <- 1

# Save everything ---------------------------------------------------------
save(arcticBRT_arcticFW, arcticBRT_europeFW, arcticBRT_pyreneesFW, arcticBRT_serengetiFW,
     europeBRT_arcticFW, europeBRT_europeFW, europeBRT_pyreneesFW, europeBRT_serengetiFW,
     pyreneesBRT_arcticFW, pyreneesBRT_europeFW, pyreneesBRT_pyreneesFW, pyreneesBRT_serengetiFW,
     serengetiBRT_arcticFW, serengetiBRT_europeFW, serengetiBRT_pyreneesFW, serengetiBRT_pyreneesFW,
     file = "data/checkpoints/BRT_predictions.RData")

