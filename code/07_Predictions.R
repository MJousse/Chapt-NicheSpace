# Step 07: Predict each food web (target food web) using each model (calibrated on source food web)
# 1. Load and standardize data
# For each model (one per food web):
# 2. Standardize predictors of target food web with the mean and sd of source food web
# 3. Make predictions for each target food web with 100 posterior samples
# Keep the mean, sd, 95% Credible interval of each prediction +
# the expected value, and whether the data was used for training
# 4. Save everything

# Predictions
rm(list = ls())
library(dplyr)
library(brms)
library(tidyr)
source("code/functions.R")

# Load model and food webs ------------------------------------------------
# models
ArcticModel <- readRDS("~/OneDrive/Chapt-NicheSpace/models/ArcticModel.rds")
EuropeModel <- readRDS("~/OneDrive/Chapt-NicheSpace/models/EuroModel.rds")
PyreneesModel <- readRDS("~/OneDrive/Chapt-NicheSpace/models/PyreneesModel.rds")
SerengetiModel <- readRDS("~/OneDrive/Chapt-NicheSpace/models/SerengetiModel.rds")

load("data/checkpoints/train_test_splits.RData")

# functional traits
FuncTraits <- read.csv("data/cleaned/SpeciesTraitsFull.csv", row.names = 1)

# Transform into predictors and scale
predictors <- get_predictors(FuncTraits$Species, FuncTraits)

predictors <- mutate_at(predictors, vars(Habitat_breadth.predator, BM.predator:ClutchSize.predator,
                                         Habitat_breadth.prey, BM.prey:ClutchSize.prey, 
                                         Habitat.match:BM.match), scale2)

# european metaweb
EuroInteractions <- read.csv("data/cleaned/EuroFWadults.csv", row.names = 1)
EuroSpecies <- read.csv("data/cleaned/EuroMWTaxo.csv", row.names = 1) %>%
  distinct() %>%
  filter(!is.na(Species))

# keep predictors of species within Europe
EuroMW <- filter(predictors, Predator %in% EuroSpecies$Species, Prey %in% EuroSpecies$Species)

# add response
EuroInteractions$interaction <- 1
EuroMW <- left_join(EuroMW, EuroInteractions)
EuroMW$interaction[is.na(EuroMW$interaction)] <- 0

# high arctic
HighArcticInteractions <- read.csv("data/cleaned/HighArcticFW.csv", row.names = 1)
HighArcticSpecies <- read.csv("data/cleaned/HighArcticFWTaxo.csv", row.names = 1) %>%
  distinct() %>%
  filter(!is.na(Species))

# keep predictors of species within Europe
HighArcticFW <- filter(predictors, Predator %in% HighArcticSpecies$Species, Prey %in% HighArcticSpecies$Species)

# add response
HighArcticInteractions$interaction <- 1
HighArcticFW <- left_join(HighArcticFW, HighArcticInteractions)
HighArcticFW$interaction[is.na(HighArcticFW$interaction)] <- 0

# pyrenees
PyreneesInteractions <- read.csv("data/cleaned/pyrenneesFW.csv", row.names = 1)
PyreneesSpecies <- read.csv("data/cleaned/pyrenneesFWTaxo.csv", row.names = 1) %>%
  distinct() %>%
  filter(!is.na(Species))

# keep predictors of species within Europe
PyreneesFW <- filter(predictors, Predator %in% PyreneesSpecies$Species, Prey %in% PyreneesSpecies$Species)

# add response
PyreneesInteractions$interaction <- 1
PyreneesFW <- left_join(PyreneesFW, PyreneesInteractions)
PyreneesFW$interaction[is.na(PyreneesFW$interaction)] <- 0

# serengeti
SerengetiInteractions <- read.csv("data/cleaned/SerengetiFW.csv", row.names = 1) %>%
  select(Predator = Consumer_Species, Prey = Resource_Species)
SerengetiSpecies <- read.csv("data/cleaned/SerengetiFWTaxo.csv", row.names = 1) %>%
  distinct() %>%
  filter(!is.na(Species))

# keep predictors of species within Europe
SerengetiFW <- filter(predictors, Predator %in% SerengetiSpecies$Species, Prey %in% SerengetiSpecies$Species)

# add response
SerengetiInteractions$interaction <- 1
SerengetiFW <- left_join(SerengetiFW, SerengetiInteractions)
SerengetiFW$interaction[is.na(SerengetiFW$interaction)] <- 0

# Predictions using the Arctic model --------------------------------------
# predict arctic food webs
Arctic_Arctic_predictions <- make_predictions(ArcticModel, newdata = HighArcticFW, 
                                              allow_new_levels = TRUE, ndraws = 100)
Arctic_Arctic_predictions$testing <- 0
Arctic_Arctic_predictions$testing[testing_id_higharctic] <- 1

# predict european metaweb
Arctic_Euro_predictions <- make_predictions(ArcticModel, newdata = EuroMW, 
                                            allow_new_levels = TRUE, ndraws = 100)

# predict pyrenees food web
Arctic_Pyrenees_predictions <- make_predictions(ArcticModel, newdata = PyreneesFW, 
                                                allow_new_levels = TRUE, ndraws = 100)

# predict serengeti food web
Arctic_Serengeti_predictions <- make_predictions(ArcticModel, newdata = SerengetiFW, 
                                                 allow_new_levels = TRUE, ndraws = 100)

# Predictions using the European model --------------------------------------
# predict arctic food webs
Euro_Arctic_predictions <- make_predictions(EuropeModel, newdata = HighArcticFW, 
                                            allow_new_levels = TRUE, ndraws = 100)

# predict european metaweb
Euro_Euro_predictions <- make_predictions(EuropeModel, newdata = EuroMW, 
                                          allow_new_levels = TRUE, ndraws = 100)
Euro_Euro_predictions$testing <- 0
Euro_Euro_predictions$testing[testing_id_euro] <- 1

# predict pyrenees food web
Euro_Pyrenees_predictions <- make_predictions(EuropeModel, newdata = PyreneesFW, 
                                              allow_new_levels = TRUE, ndraws = 100)

# predict serengeti food web
Euro_Serengeti_predictions <- make_predictions(EuropeModel, newdata = SerengetiFW, 
                                               allow_new_levels = TRUE, ndraws = 100)

# Predictions using the Pyrenees model --------------------------------------
# predict arctic food webs
Pyrenees_Arctic_predictions <- make_predictions(PyreneesModel, newdata = HighArcticFW, 
                                                allow_new_levels = TRUE, ndraws = 100)

# predict european metaweb
Pyrenees_Euro_predictions <- make_predictions(PyreneesModel, newdata = EuroMW, 
                                              allow_new_levels = TRUE, ndraws = 100)

# predict pyrenees food web
Pyrenees_Pyrenees_predictions <- make_predictions(PyreneesModel, newdata = PyreneesFW, 
                                                  allow_new_levels = TRUE, ndraws = 100)
Pyrenees_Pyrenees_predictions$testing <- 0
Pyrenees_Pyrenees_predictions$testing[testing_id_pyrenees] <- 1

# predict serengeti food web
Pyrenees_Serengeti_predictions <- make_predictions(PyreneesModel, newdata = SerengetiFW, 
                                                   allow_new_levels = TRUE, ndraws = 100)

# Predictions using the Serengeti model --------------------------------------
# predict arctic food webs
Serengeti_Arctic_predictions <- make_predictions(SerengetiModel, newdata = HighArcticFW, 
                                                 allow_new_levels = TRUE, ndraws = 100)

# predict european metaweb
Serengeti_Euro_predictions <- make_predictions(SerengetiModel, newdata = EuroMW, 
                                               allow_new_levels = TRUE, ndraws = 100)

# predict pyrenees food web
Serengeti_Pyrenees_predictions <- make_predictions(SerengetiModel, newdata = PyreneesFW, 
                                                   allow_new_levels = TRUE, ndraws = 100)

# predict serengeti food web
Serengeti_Serengeti_predictions <- make_predictions(SerengetiModel, newdata = SerengetiFW, 
                                                    allow_new_levels = TRUE, ndraws = 100)
Serengeti_Serengeti_predictions$testing <- 0
Serengeti_Serengeti_predictions$testing[testing_id_serengeti] <- 1

# Save everything ---------------------------------------------------------
save(Arctic_Arctic_predictions, Arctic_Euro_predictions, Arctic_Pyrenees_predictions, Arctic_Serengeti_predictions,
     Euro_Arctic_predictions, Euro_Euro_predictions, Euro_Pyrenees_predictions, Euro_Serengeti_predictions,
     Pyrenees_Arctic_predictions, Pyrenees_Euro_predictions, Pyrenees_Pyrenees_predictions, Pyrenees_Serengeti_predictions,
     Serengeti_Arctic_predictions, Serengeti_Euro_predictions, Serengeti_Pyrenees_predictions, Serengeti_Serengeti_predictions, 
     file = "~/OneDrive/Chapt-NicheSpace/predictions.RData")

# Predictions for and with the alternative Serengeti food web-------------
load("~/OneDrive/Chapt-NicheSpace/models/SerengetiModel2.RData")
Serengeti2_Serengeti_predictions <- make_predictions(SerengetiModel2, newdata = SerengetiFW2, 
                                                     allow_new_levels = TRUE, ndraws = 100)
Serengeti2_Serengeti_predictions$testing <- 0
Serengeti2_Serengeti_predictions$testing[testing_id_ser] <- 1

# predict arctic food webs
Serengeti2_Arctic_predictions <- make_predictions(SerengetiModel2, newdata = HighArcticFW, 
                                                 allow_new_levels = TRUE, ndraws = 100)

# predict european metaweb
Serengeti2_Euro_predictions <- make_predictions(SerengetiModel2, newdata = EuroMW, 
                                               allow_new_levels = TRUE, ndraws = 100)

# predict pyrenees food web
Serengeti2_Pyrenees_predictions <- make_predictions(SerengetiModel2, newdata = PyreneesFW, 
                                                   allow_new_levels = TRUE, ndraws = 100)

# predict arctic food webs
Arctic_Serengeti2_predictions <- make_predictions(ArcticModel, newdata = SerengetiFW2, 
                                                  allow_new_levels = TRUE, ndraws = 100)

# predict european metaweb
Euro_Serengeti2_predictions <- make_predictions(EuropeModel, newdata = SerengetiFW2, 
                                                allow_new_levels = TRUE, ndraws = 100)

# predict pyrenees food web
Pyrenees_Serengeti2_predictions <- make_predictions(PyreneesModel, newdata = SerengetiFW2, 
                                                    allow_new_levels = TRUE, ndraws = 100)

save(Serengeti2_Arctic_predictions, Serengeti2_Serengeti_predictions, Serengeti2_Euro_predictions,
     Serengeti2_Pyrenees_predictions, Arctic_Serengeti2_predictions, Euro_Serengeti2_predictions,
     Pyrenees_Serengeti2_predictions, file = "~/OneDrive/Chapt-NicheSpace/predictions_serengeti_alt.RData")