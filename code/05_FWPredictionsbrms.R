# Predictions
library(dplyr)
library(brms)
library(tidyr)
source("code/functions.R")

# Load model and food webs ------------------------------------------------
# models
ArcticModel <- readRDS("~/OneDrive/Chapt-NicheSpace/models/ArcticModel_brms.rds")
EuropeModel <- readRDS("~/OneDrive/Chapt-NicheSpace/models/EuroModel_brms.rds")
PyreneesModel <- readRDS("~/OneDrive/Chapt-NicheSpace/models/PyreneesModel_brms.rds")
SerengetiModel <- readRDS("~/OneDrive/Chapt-NicheSpace/models/SerengetiModel_brms.rds")

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
EuroMW$Order.predator <- factor(EuroMW$Order.predator, 
                                levels = unique(FuncTraits$Order))

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

# Predictions using the Arctic model --------------------------------------
# find the mean and sd of the predictors to scale new data
var2scale <- c("Habitat_breadth.predator", "BM.predator", "Longevity.predator", 
               "ClutchSize.predator", "Habitat_breadth.prey", "BM.prey",
               "Longevity.prey", "ClutchSize.prey", "Habitat.match", "BM.match")
Arctic_mean <- apply(HighArcticFW[,var2scale], MARGIN = 2, mean)
Arctic_sd <- apply(HighArcticFW[,var2scale], MARGIN = 2, sd)

# scale predictors
HighArcticFWscaled <- HighArcticFW
HighArcticFWscaled[,var2scale] <- sweep(HighArcticFW[,var2scale], MARGIN = 2, Arctic_mean)
HighArcticFWscaled[,var2scale] <- sweep(HighArcticFW[,var2scale], MARGIN = 2, FUN = "/", 2*Arctic_sd)
EuroMWscaled <- EuroMW
EuroMWscaled[,var2scale] <- sweep(EuroMW[,var2scale], MARGIN = 2, Arctic_mean)
EuroMWscaled[,var2scale] <- sweep(EuroMW[,var2scale], MARGIN = 2, FUN = "/", 2*Arctic_sd)
PyreneesFWscaled <- PyreneesFW
PyreneesFWscaled[,var2scale] <- sweep(PyreneesFW[,var2scale], MARGIN = 2, Arctic_mean)
PyreneesFWscaled[,var2scale] <- sweep(PyreneesFW[,var2scale], MARGIN = 2, FUN = "/", 2*Arctic_sd)
SerengetiFWscaled <- SerengetiFW
SerengetiFWscaled[,var2scale] <- sweep(SerengetiFW[,var2scale], MARGIN = 2, Arctic_mean)
SerengetiFWscaled[,var2scale] <- sweep(SerengetiFW[,var2scale], MARGIN = 2, FUN = "/", 2*Arctic_sd)

# predict arctic food webs
Arctic_Arctic_predictions <- predict(ArcticModel, newdata = HighArcticFWscaled, allow_new_levels = TRUE)
Arctic_Arctic_predictions <- select(HighArcticFW, Predator, Prey, interaction) %>%
  bind_cols(as.data.frame(Arctic_Arctic_predictions))
Arctic_Arctic_predictions$training <- ifelse(c(1:nrow(HighArcticFW) %in% as.numeric(rownames(ArcticModel$data))), 1, 0)

# predict european metaweb
Arctic_Euro_predictions <- predict(ArcticModel, newdata = EuroMWscaled, allow_new_levels = TRUE)
Arctic_Euro_predictions <- select(EuroMW, Predator, Prey, interaction) %>%
  bind_cols(as.data.frame(Arctic_Euro_predictions))

# predict pyrenees food web
Arctic_Pyrenees_predictions <- predict(ArcticModel, newdata = PyreneesFWscaled, allow_new_levels = TRUE)
Arctic_Pyrenees_predictions <- select(PyreneesFW, Predator, Prey, interaction) %>%
  bind_cols(as.data.frame(Arctic_Pyrenees_predictions))

# predict serengeti food web
Arctic_Serengeti_predictions <- predict(ArcticModel, newdata = SerengetiFWscaled, allow_new_levels = TRUE)
Arctic_Serengeti_predictions <- select(SerengetiFW, Predator, Prey, interaction) %>%
  bind_cols(as.data.frame(Arctic_Serengeti_predictions))

# Predictions using the European model --------------------------------------
# find the mean and sd of the predictors to scale new data
Europe_mean <- apply(EuroMW[,var2scale], MARGIN = 2, mean)
Europe_sd <- apply(EuroMW[,var2scale], MARGIN = 2, sd)

# scale predictors
HighArcticFWscaled <- HighArcticFW
HighArcticFWscaled[,var2scale] <- sweep(HighArcticFW[,var2scale], MARGIN = 2, Europe_mean)
HighArcticFWscaled[,var2scale] <- sweep(HighArcticFW[,var2scale], MARGIN = 2, FUN = "/", 2*Europe_sd)
EuroMWscaled <- EuroMW
EuroMWscaled[,var2scale] <- sweep(EuroMW[,var2scale], MARGIN = 2, Europe_mean)
EuroMWscaled[,var2scale] <- sweep(EuroMW[,var2scale], MARGIN = 2, FUN = "/", 2*Europe_sd)
PyreneesFWscaled <- PyreneesFW
PyreneesFWscaled[,var2scale] <- sweep(PyreneesFW[,var2scale], MARGIN = 2, Europe_mean)
PyreneesFWscaled[,var2scale] <- sweep(PyreneesFW[,var2scale], MARGIN = 2, FUN = "/", 2*Europe_sd)
SerengetiFWscaled <- SerengetiFW
SerengetiFWscaled[,var2scale] <- sweep(SerengetiFW[,var2scale], MARGIN = 2, Europe_mean)
SerengetiFWscaled[,var2scale] <- sweep(SerengetiFW[,var2scale], MARGIN = 2, FUN = "/", 2*Europe_sd)

# predict arctic food webs
Euro_Arctic_predictions <- predict(EuropeModel, newdata = HighArcticFWscaled, allow_new_levels = TRUE)
Euro_Arctic_predictions <- select(HighArcticFW, Predator, Prey, interaction) %>%
  bind_cols(as.data.frame(Euro_Arctic_predictions))

# predict european metaweb
Euro_Euro_predictions <- predict(EuropeModel, newdata = EuroMWscaled, allow_new_levels = TRUE)
Euro_Euro_predictions <- select(EuroMW, Predator, Prey, interaction) %>%
  bind_cols(as.data.frame(Euro_Euro_predictions))
Euro_Euro_predictions$training <- ifelse(c(1:nrow(EuroMW) %in% as.numeric(rownames(EuropeModel$data))), 1, 0)

# predict pyrenees food web
Euro_Pyrenees_predictions <- predict(EuropeModel, newdata = PyreneesFWscaled, allow_new_levels = TRUE)
Euro_Pyrenees_predictions <- select(PyreneesFW, Predator, Prey, interaction) %>%
  bind_cols(as.data.frame(Euro_Pyrenees_predictions))

# predict serengeti food web
Euro_Serengeti_predictions <- predict(EuropeModel, newdata = SerengetiFWscaled, allow_new_levels = TRUE)
Euro_Serengeti_predictions <- select(SerengetiFW, Predator, Prey, interaction) %>%
  bind_cols(as.data.frame(Euro_Serengeti_predictions))

# Predictions using the Pyrenees model --------------------------------------
# find the mean and sd of the predictors to scale new data
Pyrenees_mean <- apply(PyreneesFW[,var2scale], MARGIN = 2, mean)
Pyrenees_sd <- apply(PyreneesFW[,var2scale], MARGIN = 2, sd)

# scale predictors
HighArcticFWscaled <- HighArcticFW
HighArcticFWscaled[,var2scale] <- sweep(HighArcticFW[,var2scale], MARGIN = 2, Pyrenees_mean)
HighArcticFWscaled[,var2scale] <- sweep(HighArcticFW[,var2scale], MARGIN = 2, FUN = "/", 2*Pyrenees_sd)
EuroMWscaled <- EuroMW
EuroMWscaled[,var2scale] <- sweep(EuroMW[,var2scale], MARGIN = 2, Pyrenees_mean)
EuroMWscaled[,var2scale] <- sweep(EuroMW[,var2scale], MARGIN = 2, FUN = "/", 2*Pyrenees_sd)
PyreneesFWscaled <- PyreneesFW
PyreneesFWscaled[,var2scale] <- sweep(PyreneesFW[,var2scale], MARGIN = 2, Pyrenees_mean)
PyreneesFWscaled[,var2scale] <- sweep(PyreneesFW[,var2scale], MARGIN = 2, FUN = "/", 2*Pyrenees_sd)
SerengetiFWscaled <- SerengetiFW
SerengetiFWscaled[,var2scale] <- sweep(SerengetiFW[,var2scale], MARGIN = 2, Pyrenees_mean)
SerengetiFWscaled[,var2scale] <- sweep(SerengetiFW[,var2scale], MARGIN = 2, FUN = "/", 2*Pyrenees_sd)

# predict arctic food webs
Pyrenees_Arctic_predictions <- predict(PyreneesModel, newdata = HighArcticFWscaled, allow_new_levels = TRUE)
Pyrenees_Arctic_predictions <- select(HighArcticFW, Predator, Prey, interaction) %>%
  bind_cols(as.data.frame(Pyrenees_Arctic_predictions))

# predict european metaweb
Pyrenees_Euro_predictions <- predict(PyreneesModel, newdata = EuroMWscaled, allow_new_levels = TRUE)
Pyrenees_Euro_predictions <- select(EuroMW, Predator, Prey, interaction) %>%
  bind_cols(as.data.frame(Pyrenees_Euro_predictions))

# predict pyrenees food web
Pyrenees_Pyrenees_predictions <- predict(PyreneesModel, newdata = PyreneesFWscaled, allow_new_levels = TRUE)
Pyrenees_Pyrenees_predictions <- select(PyreneesFW, Predator, Prey, interaction) %>%
  bind_cols(as.data.frame(Pyrenees_Pyrenees_predictions))
Pyrenees_Pyrenees_predictions$training <- ifelse(c(1:nrow(PyreneesFW) %in% as.numeric(rownames(PyreneesModel$data))), 1, 0)

# predict serengeti food web
Pyrenees_Serengeti_predictions <- predict(PyreneesModel, newdata = SerengetiFWscaled, allow_new_levels = TRUE)
Pyrenees_Serengeti_predictions <- select(SerengetiFW, Predator, Prey, interaction) %>%
  bind_cols(as.data.frame(Pyrenees_Serengeti_predictions))

# Predictions using the Serengeti model --------------------------------------
# find the mean and sd of the predictors to scale new data
Serengeti_mean <- apply(EuroMW[,var2scale], MARGIN = 2, mean)
Serengeti_sd <- apply(EuroMW[,var2scale], MARGIN = 2, sd)

# scale predictors
HighArcticFWscaled <- HighArcticFW
HighArcticFWscaled[,var2scale] <- sweep(HighArcticFW[,var2scale], MARGIN = 2, Serengeti_mean)
HighArcticFWscaled[,var2scale] <- sweep(HighArcticFW[,var2scale], MARGIN = 2, FUN = "/", 2*Serengeti_sd)
EuroMWscaled <- EuroMW
EuroMWscaled[,var2scale] <- sweep(EuroMW[,var2scale], MARGIN = 2, Serengeti_mean)
EuroMWscaled[,var2scale] <- sweep(EuroMW[,var2scale], MARGIN = 2, FUN = "/", 2*Serengeti_sd)
PyreneesFWscaled <- PyreneesFW
PyreneesFWscaled[,var2scale] <- sweep(PyreneesFW[,var2scale], MARGIN = 2, Serengeti_mean)
PyreneesFWscaled[,var2scale] <- sweep(PyreneesFW[,var2scale], MARGIN = 2, FUN = "/", 2*Serengeti_sd)
SerengetiFWscaled <- SerengetiFW
SerengetiFWscaled[,var2scale] <- sweep(SerengetiFW[,var2scale], MARGIN = 2, Serengeti_mean)
SerengetiFWscaled[,var2scale] <- sweep(SerengetiFW[,var2scale], MARGIN = 2, FUN = "/", 2*Serengeti_sd)

# predict arctic food webs
Serengeti_Arctic_predictions <- predict(SerengetiModel, newdata = HighArcticFWscaled, allow_new_levels = TRUE)
Serengeti_Arctic_predictions <- select(HighArcticFW, Predator, Prey, interaction) %>%
  bind_cols(as.data.frame(Serengeti_Arctic_predictions))

# predict european metaweb
Serengeti_Euro_predictions <- predict(SerengetiModel, newdata = EuroMWscaled, allow_new_levels = TRUE)
Serengeti_Euro_predictions <- select(EuroMW, Predator, Prey, interaction) %>%
  bind_cols(as.data.frame(Serengeti_Euro_predictions))

# predict pyrenees food web
Serengeti_Pyrenees_predictions <- predict(SerengetiModel, newdata = PyreneesFWscaled, allow_new_levels = TRUE)
Serengeti_Pyrenees_predictions <- select(PyreneesFW, Predator, Prey, interaction) %>%
  bind_cols(as.data.frame(Serengeti_Pyrenees_predictions))

# predict serengeti food web
Serengeti_Serengeti_predictions <- predict(SerengetiModel, newdata = SerengetiFWscaled, allow_new_levels = TRUE)
Serengeti_Serengeti_predictions <- select(SerengetiFW, Predator, Prey, interaction) %>%
  bind_cols(as.data.frame(Serengeti_Serengeti_predictions))
Serengeti_Serengeti_predictions$training <- ifelse(c(1:nrow(SerengetiFW) %in% as.numeric(rownames(SerengetiModel$data))), 1, 0)
