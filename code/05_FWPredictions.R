# Step 5: Predict the three food webs with the model
# For each food web:
# 1. Prepare data
# 2. Make predictions
# 3. Plot predictons
# 4. Calulate performance

rm(list = ls())
set.seed(16)
library(bayesplot)
library(ROCR)
library(greta)
library(dplyr)
library(ggplot2)
library(tidyr)
source("code/functions.R")
load("data/models/GLMMadults_21032022.RData")

# Get the mean and sd of predictors in Europe for scaling -----------------
# load data
FuncTraits <- read.csv("data/cleaned/SpeciesTraitsFull.csv", row.names = 1)
EuroSpecies <- read.csv("data/cleaned/EuroMWTaxo.csv", row.names = 1) %>%
  distinct() %>%
  filter(!is.na(Species))

# transform traits into predictors
EuroMW <- get_predictors(EuroSpecies$Species, FuncTraits)

# get mean and sd of the predictors to scale
columns_to_scale <- c("Habitat_breadth.predator", "BM.predator", 
                     "ClutchSize.predator", "Longevity.predator", 
                     "Habitat_breadth.prey", "BM.prey", "ClutchSize.prey",
                     "Longevity.prey", "Habitat.match", "BM.match")
Predictor_means <- apply(select(EuroMW, all_of(columns_to_scale)),
                         MARGIN = 2, mean)
Predictor_sd <- apply(select(EuroMW, all_of(columns_to_scale)),
                         MARGIN = 2, sd)

# High Arctic food web ----------------------------------------------------
# load data
HighArcticInteractions <- read.csv("data/cleaned/HighArcticFW.csv", row.names = 1)
HighArcticInteractions$interaction <- 1
HighArcticSpecies <- read.csv("data/cleaned/HighArcticFWTaxo.csv", row.names = 1) %>%
  distinct() %>%
  filter(!is.na(Species))

# prepare predictors
HighArcticFW <- get_predictors(HighArcticSpecies$Species, FuncTraits)

# scale predictors
HighArcticFW[, columns_to_scale] <- sweep(HighArcticFW[, columns_to_scale], MARGIN = 2, Predictor_means)
HighArcticFW[, columns_to_scale] <- sweep(HighArcticFW[, columns_to_scale], MARGIN = 2, 2*Predictor_sd, FUN = "/")

# predictions
predictions <- make_predictions(HighArcticFW, unique(FuncTraits$Order), coef, GLMM) %>%
  left_join(HighArcticInteractions) %>% replace_na(list(interaction = 0))

# calculate performance
auc <- performance(prediction(predictions$predictions, predictions$interaction), "auc")@y.values[[1]]
aucpr <- performance(prediction(predictions$predictions, predictions$interaction), "aucpr")@y.values[[1]]

# plot predictions
ggplot(predictions) +
  geom_tile(aes(Predator, Prey, fill = predictions)) +
  geom_point(data = filter(predictions, interaction == 1), aes(Predator, Prey), size = 5) +
  scale_fill_distiller(palette = "YlGn") +
  theme(axis.text.x = element_text(angle = 90))

# Serengeti food web ------------------------------------------------------
# load data
SerengetiInteractions <- read.csv("data/cleaned/SerengetiFW.csv", row.names = 1) %>%
  select(Predator = Consumer_Species, Prey = Resource_Species)
SerengetiInteractions$interaction <- 1
SerengetiSpecies <- read.csv("data/cleaned/SerengetiFWTaxo.csv", row.names = 1) %>%
  distinct() %>%
  filter(!is.na(Species))

# prepare predictors
SerengetiFW <- get_predictors(SerengetiSpecies$Species, FuncTraits)

# scale predictors
SerengetiFW[, columns_to_scale] <- sweep(SerengetiFW[, columns_to_scale], MARGIN = 2, Predictor_means)
SerengetiFW[, columns_to_scale] <- sweep(SerengetiFW[, columns_to_scale], MARGIN = 2, 2*Predictor_sd, FUN = "/")

# make predictions
predictions <- make_predictions(SerengetiFW, unique(FuncTraits$Order), coef, GLMM) %>%
  left_join(SerengetiInteractions) %>% replace_na(list(interaction = 0))

# calculate performance
auc <- performance(prediction(predictions$predictions, predictions$interaction), "auc")@y.values[[1]]
aucpr <- performance(prediction(predictions$predictions, predictions$interaction), "aucpr")@y.values[[1]]

# plot predictions
ggplot(predictions) +
  geom_tile(aes(Predator, Prey, fill = predictions)) +
  geom_point(data = filter(predictions, interaction == 1), aes(Predator, Prey), size = 1) +
  scale_fill_distiller(palette = "YlGn") +
  theme(axis.text.x = element_text(angle = 90))

# Pyrenees food web -------------------------------------------------------
# load data
PyreneesInteractions <- read.csv("data/cleaned/pyrenneesFW.csv", row.names = 1)
PyreneesInteractions$interaction <- 1
PyreneesSpecies <- read.csv("data/cleaned/pyrenneesFWTaxo.csv", row.names = 1) %>%
  distinct() %>%
  filter(!is.na(Species))

# prepare predictors
PyreneesFW <- get_predictors(PyreneesSpecies$Species, FuncTraits) %>%
  na.omit()

# scale predictors
PyreneesFW[, columns_to_scale] <- sweep(PyreneesFW[, columns_to_scale], MARGIN = 2, Predictor_means)
PyreneesFW[, columns_to_scale] <- sweep(PyreneesFW[, columns_to_scale], MARGIN = 2, 2*Predictor_sd, FUN = "/")

# make predictions
predictions <- make_predictions(PyreneesFW, unique(FuncTraits$Order), coef, GLMM) %>%
  left_join(PyreneesInteractions) %>% replace_na(list(interaction = 0))

# measure performance
auc <- performance(prediction(predictions$predictions, predictions$interaction), "auc")@y.values[[1]]
aucpr <- performance(prediction(predictions$predictions, predictions$interaction), "aucpr")@y.values[[1]]

# plot predictions
ggplot(predictions) +
  geom_tile(aes(Predator, Prey, fill = predictions)) +
  geom_point(data = filter(predictions, interaction == 1), aes(Predator, Prey), size = 2) +
  scale_fill_distiller(palette = "YlGn") +
  theme(axis.text.x = element_text(angle = 90))
