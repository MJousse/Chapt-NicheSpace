# 1. Train BRTs using the same data as the HGLMs
# 2. Make predictions
# 3. Measure performance
# 4. Correlate to distance
# 5. Measure how well food web properties are predicted

# Set up
library(dismo)
library(gbm)
library(dplyr)
load("data/checkpoints/train_test_splits.RData")
source("code/functions.R")
FuncTraits <- read.csv("data/cleaned/SpeciesTraitsFull.csv", row.names = 1)

# Europe ------------------------------------------------------------------
# load data and standardize
EuroInteractions <- read.csv("data/cleaned/EuroFWadults.csv", row.names = 1)
EuroSpecies <- read.csv("data/cleaned/EuroMWTaxo.csv", row.names = 1) %>%
  distinct() %>%
  filter(!is.na(Species))

# transform trait into predictors for every species pair
EuroMW <- get_predictors(EuroSpecies$Species, FuncTraits)

# scale predictors
EuroMW <- mutate_at(EuroMW, vars(Habitat_breadth.predator, BM.predator:ClutchSize.predator,
                                 Habitat_breadth.prey, BM.prey:ClutchSize.prey, 
                                 Habitat.match:BM.match), scale2)

# add response
EuroInteractions$interaction <- 1
EuroMW <- left_join(EuroMW, EuroInteractions)
EuroMW$interaction[is.na(EuroMW$interaction)] <- 0

# prepare training set
EuroMW$TL.predator <- ifelse(EuroMW$Herbivore.predator == 1, "Herbivore", ifelse(EuroMW$Omnivore.predator == 1, "Omnivore", "Carnivore"))
EuroMW$TL.prey <- ifelse(EuroMW$Herbivore.prey == 1, "Herbivore", ifelse(EuroMW$Omnivore.prey == 1, "Omnivore", "Carnivore"))

EuroMW <- select(EuroMW,  -Predator, -Prey, -Order.prey, -Herbivore.predator, 
                       -Omnivore.predator, -Carnivore.predator, -Herbivore.prey, 
                       -Omnivore.prey, -Carnivore.prey) %>%
  mutate(Order.predator = factor(Order.predator, 
                                 levels = unique(FuncTraits$Order)),
         TL.predator = factor(TL.predator, levels = c("Herbivore", "Omnivore", "Carnivore")),
         TL.prey = factor(TL.prey, levels = c("Herbivore", "Omnivore", "Carnivore")),
         ActivityTime.match = factor(ActivityTime.match))

training <- EuroMW[training_id_euro, ]

brt_europe <- gbm.step(data=training, gbm.x = c(1,3:12,14,15), gbm.y = 13,
                       family = "bernoulli", tree.complexity = 5,
                       learning.rate = 0.01, bag.fraction = 0.5)

# Nunavik ------------------------------------------------------------------
# load data and standardize
HighArcticInteractions <- read.csv("data/cleaned/HighArcticFW.csv", row.names = 1)
HighArcticSpecies <- read.csv("data/cleaned/HighArcticFWTaxo.csv", row.names = 1) %>%
  distinct() %>%
  filter(!is.na(Species))

# transform trait into predictors for every species pair
HighArcticFW <- get_predictors(HighArcticSpecies$Species, FuncTraits)

# scale predictors
HighArcticFW <- mutate_at(HighArcticFW, vars(Habitat_breadth.predator, BM.predator:ClutchSize.predator,
                                             Habitat_breadth.prey, BM.prey:ClutchSize.prey, 
                                             Habitat.match:BM.match), scale2)

# add response
HighArcticInteractions$interaction <- 1
HighArcticFW <- left_join(HighArcticFW, HighArcticInteractions)
HighArcticFW$interaction[is.na(HighArcticFW$interaction)] <- 0

# prepare training set
HighArcticFW$TL.predator <- ifelse(HighArcticFW$Herbivore.predator == 1, "Herbivore", ifelse(HighArcticFW$Omnivore.predator == 1, "Omnivore", "Carnivore"))
HighArcticFW$TL.prey <- ifelse(HighArcticFW$Herbivore.prey == 1, "Herbivore", ifelse(HighArcticFW$Omnivore.prey == 1, "Omnivore", "Carnivore"))

HighArcticFW <- select(HighArcticFW,  -Predator, -Prey, -Order.prey, -Herbivore.predator, 
                   -Omnivore.predator, -Carnivore.predator, -Herbivore.prey, 
                   -Omnivore.prey, -Carnivore.prey) %>%
  mutate(Order.predator = factor(Order.predator, 
                                 levels = unique(FuncTraits$Order)),
         TL.predator = factor(TL.predator, levels = c("Herbivore", "Omnivore", "Carnivore")),
         TL.prey = factor(TL.prey, levels = c("Herbivore", "Omnivore", "Carnivore")),
         ActivityTime.match = factor(ActivityTime.match))

training <- HighArcticFW[training_id_higharctic, ]

brt_arctic <- gbm.step(data=training, gbm.x = c(1,3:12,14,15), gbm.y = 13,
                       family = "bernoulli", tree.complexity = 5,
                       learning.rate = 0.01, bag.fraction = 0.5)


# Pyrenees ----------------------------------------------------------------
# load data and standardize
PyreneesInteractions <- read.csv("data/cleaned/pyrenneesFW.csv", row.names = 1)
PyreneesSpecies <- read.csv("data/cleaned/pyrenneesFWTaxo.csv", row.names = 1) %>%
  distinct() %>%
  filter(!is.na(Species))

# transform trait into predictors for every species pair
PyreneesFW <- get_predictors(PyreneesSpecies$Species, FuncTraits)

# scale predictors
PyreneesFW <- mutate_at(PyreneesFW, vars(Habitat_breadth.predator, BM.predator:ClutchSize.predator,
                                         Habitat_breadth.prey, BM.prey:ClutchSize.prey, 
                                         Habitat.match:BM.match), scale2)

# add response
PyreneesInteractions$interaction <- 1
PyreneesFW <- left_join(PyreneesFW, PyreneesInteractions)
PyreneesFW$interaction[is.na(PyreneesFW$interaction)] <- 0

# prepare training set
PyreneesFW$TL.predator <- ifelse(PyreneesFW$Herbivore.predator == 1, "Herbivore", ifelse(PyreneesFW$Omnivore.predator == 1, "Omnivore", "Carnivore"))
PyreneesFW$TL.prey <- ifelse(PyreneesFW$Herbivore.prey == 1, "Herbivore", ifelse(PyreneesFW$Omnivore.prey == 1, "Omnivore", "Carnivore"))

PyreneesFW <- select(PyreneesFW,  -Predator, -Prey, -Order.prey, -Herbivore.predator, 
                       -Omnivore.predator, -Carnivore.predator, -Herbivore.prey, 
                       -Omnivore.prey, -Carnivore.prey) %>%
  mutate(Order.predator = factor(Order.predator, 
                                 levels = unique(FuncTraits$Order)),
         TL.predator = factor(TL.predator, levels = c("Herbivore", "Omnivore", "Carnivore")),
         TL.prey = factor(TL.prey, levels = c("Herbivore", "Omnivore", "Carnivore")),
         ActivityTime.match = factor(ActivityTime.match))

training <- PyreneesFW[training_id_pyrenees, ]

brt_pyrenees <- gbm.step(data=training, gbm.x = c(1,3:12,14,15), gbm.y = 13,
                       family = "bernoulli", tree.complexity = 5,
                       learning.rate = 0.01, bag.fraction = 0.5)

# Serengeti ---------------------------------------------------------------
# load data and standardize
SerengetiInteractions <- read.csv("data/cleaned/SerengetiFW.csv", row.names = 1) %>%
  select(Predator = Consumer_Species, Prey = Resource_Species)
SerengetiSpecies <- read.csv("data/cleaned/SerengetiFWTaxo.csv", row.names = 1) %>%
  distinct() %>%
  filter(!is.na(Species))

# transform trait into predictors for every species pairbf(interaction ~ 1 + . + (1 + .||Order.predator), family = bernoulli())

SerengetiFW <- get_predictors(SerengetiSpecies$Species, FuncTraits)

# scale predictors
SerengetiFW <- mutate_at(SerengetiFW, vars(Habitat_breadth.predator, BM.predator:ClutchSize.predator,
                                           Habitat_breadth.prey, BM.prey:ClutchSize.prey, 
                                           Habitat.match:BM.match), scale2)

# add response
SerengetiInteractions$interaction <- 1
SerengetiFW <- left_join(SerengetiFW, SerengetiInteractions)
SerengetiFW$interaction[is.na(SerengetiFW$interaction)] <- 0

# prepare training set
SerengetiFW$TL.predator <- ifelse(SerengetiFW$Herbivore.predator == 1, "Herbivore", ifelse(SerengetiFW$Omnivore.predator == 1, "Omnivore", "Carnivore"))
SerengetiFW$TL.prey <- ifelse(SerengetiFW$Herbivore.prey == 1, "Herbivore", ifelse(SerengetiFW$Omnivore.prey == 1, "Omnivore", "Carnivore"))

SerengetiFW <- select(SerengetiFW,  -Predator, -Prey, -Order.prey, -Herbivore.predator, 
                     -Omnivore.predator, -Carnivore.predator, -Herbivore.prey, 
                     -Omnivore.prey, -Carnivore.prey) %>%
  mutate(Order.predator = factor(Order.predator, 
                                 levels = unique(FuncTraits$Order)),
         TL.predator = factor(TL.predator, levels = c("Herbivore", "Omnivore", "Carnivore")),
         TL.prey = factor(TL.prey, levels = c("Herbivore", "Omnivore", "Carnivore")),
         ActivityTime.match = factor(ActivityTime.match))

training <- SerengetiFW[training_id_serengeti, ]

brt_serengeti <- gbm.step(data=training, gbm.x = c(1,3:12,14,15), gbm.y = 13,
                         family = "bernoulli", tree.complexity = 5,
                         learning.rate = 0.01, bag.fraction = 0.5)

save(brt_serengeti, brt_europe, brt_arctic, brt_pyrenees,
     file = "data/checkpoints/brt_models.RData")
