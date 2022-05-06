# Step 03: Train the predictive model on the European metaweb
# 1. Prepare the training dataset
# 2. Set up the GLMM with greta package
# 3. Train the model

rm(list = ls())
set.seed(16)
library(dplyr)
library(brms)
source("code/functions.R")
FuncTraits <- read.csv("data/cleaned/SpeciesTraitsFull.csv", row.names = 1)
brms_form <- bf(interaction ~ 1 + 
                  (Omnivore.predator + Carnivore.predator + Habitat_breadth.predator + BM.predator + Longevity.predator + ClutchSize.predator +
                     Omnivore.prey + Carnivore.prey + Habitat_breadth.prey + BM.prey + Longevity.prey + ClutchSize.prey + 
                     ActivityTime.match + Habitat.match + BM.match) + 
                  (1 + (Omnivore.predator + Carnivore.predator + Habitat_breadth.predator + BM.predator + Longevity.predator + ClutchSize.predator +
                          Omnivore.prey + Carnivore.prey + Habitat_breadth.prey + BM.prey + Longevity.prey + ClutchSize.prey + 
                          ActivityTime.match + Habitat.match + BM.match) || Order.predator), 
                family = bernoulli())


# Prepare training dataset ------------------------------------------------
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

# take 70% of positives and as many negatives for training
N = round(0.7 * sum(EuroMW$interaction))
training_positives <- sample(seq_len(sum(EuroMW$interaction)), size = N)
training_negatives <- sample(seq_len(sum(EuroMW$interaction==0)), size = N)
training_id <- as.numeric(c(
  rownames(EuroMW[EuroMW$interaction ==1,][training_positives,]),
  rownames(EuroMW[EuroMW$interaction ==0,][training_negatives,])
))
training <- EuroMW[training_id, ]

# Set up the bayesian glmm ------------------------------------------------
# extract predictors
training <- select(training,  -Predator, -Prey, -Order.prey, 
                   -Herbivore.predator,  -Herbivore.prey) %>%
  mutate(Order.predator = factor(Order.predator, 
                                 levels = unique(FuncTraits$Order)))

get_prior(brms_form, data = training)

model_priors <- c(
  prior(normal(0, 1), class = "Intercept"),
  prior(normal(0, 1), class = "b"),
  prior(cauchy(0, 5), class = "sd")
)

prior_predictions <- brm(formula = brms_form,
                         data = training,
                         prior = model_priors,
                         sample_prior = "only", inits = "0")
EuroModel <- brm(formula = brms_form,
                 data = training,
                 prior = model_priors, sample_prior = "no", 
                    cores = 4, backend = "cmdstan", threads = 4,
                    inits = "0", iter = 2000)

# save the model on OneDrive (too big for Github...)
saveRDS(EuroModel,
        file = paste0("~/OneDrive/Chapt-NicheSpace/models/EuroModel_brms.rds"))

# -------------------------------------------------------------------------

# Prepare training dataset ------------------------------------------------
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

# take 70% of positives and as many negatives for training
N = round(0.7 * sum(HighArcticFW$interaction))
training_positives <- sample(seq_len(sum(HighArcticFW$interaction)), size = N)
training_negatives <- sample(seq_len(sum(HighArcticFW$interaction==0)), size = N)
training_id <- as.numeric(c(
  rownames(HighArcticFW[HighArcticFW$interaction ==1,][training_positives,]),
  rownames(HighArcticFW[HighArcticFW$interaction ==0,][training_negatives,])
))
training <- HighArcticFW[training_id, ]

# Set up the bayesian glmm ------------------------------------------------
# extract predictors
training <- select(training,  -Predator, -Prey, -Order.prey, 
                   -Herbivore.predator,  -Herbivore.prey) %>%
  mutate(Order.predator = factor(Order.predator, 
                                 levels = unique(FuncTraits$Order)))

get_prior(brms_form, data = training)

model_priors <- c(
  prior(normal(0, 1), class = "Intercept"),
  prior(normal(0, 1), class = "b"),
  prior(cauchy(0, 5), class = "sd")
)

prior_predictions <- brm(formula = brms_form,
                         data = training,
                         prior = model_priors,
                         sample_prior = "only", inits = "0")

ArcticModel <- brm(formula = brms_form,
                   data = training,
                   prior = model_priors, sample_prior = "no", 
                   cores = 4, backend = "cmdstan", threads = 4,
                   inits = "0", iter = 2000)

# save the model
saveRDS(ArcticModel, 
        file = paste0("~/OneDrive/Chapt-NicheSpace/models/ArcticModel_brms.rds"))

# -------------------------------------------------------------------------

# Prepare training dataset ------------------------------------------------
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

# take 70% of positives and as many negatives for training
N = round(0.7 * sum(PyreneesFW$interaction))
training_positives <- sample(seq_len(sum(PyreneesFW$interaction)), size = N)
training_negatives <- sample(seq_len(sum(PyreneesFW$interaction==0)), size = N)
training_id <- as.numeric(c(
  rownames(PyreneesFW[PyreneesFW$interaction ==1,][training_positives,]),
  rownames(PyreneesFW[PyreneesFW$interaction ==0,][training_negatives,])
))
training <- PyreneesFW[training_id, ]

# Set up the bayesian glmm ------------------------------------------------
# extract predictors
training <- select(training,  -Predator, -Prey, -Order.prey, 
                   -Herbivore.predator,  -Herbivore.prey) %>%
  mutate(Order.predator = factor(Order.predator, 
                                 levels = unique(FuncTraits$Order)))

get_prior(brms_form, data = training)

model_priors <- c(
  prior(normal(0, 1), class = "Intercept"),
  prior(normal(0, 1), class = "b"),
  prior(cauchy(0, 5), class = "sd")
)

prior_predictions <- brm(formula = brms_form,
                         data = training,
                         prior = model_priors,
                         sample_prior = "only", inits = "0")

PyreneesModel <- brm(formula = brms_form,
                     data = training,
                     prior = model_priors, sample_prior = "no", 
                     cores = 4, backend = "cmdstan", threads = 4,
                     inits = "0", iter = 2000)


# save the model
saveRDS(PyreneesModel, 
        file = paste0("~/OneDrive/Chapt-NicheSpace/models/PyreneesModel_brms.rds"))

# Prepare training dataset ------------------------------------------------
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

# take 70% of positives and as many negatives for training
N = round(0.7 * sum(SerengetiFW$interaction))
training_positives <- sample(seq_len(sum(SerengetiFW$interaction)), size = N)
training_negatives <- sample(seq_len(sum(SerengetiFW$interaction==0)), size = N)
training_id <- as.numeric(c(
  rownames(SerengetiFW[SerengetiFW$interaction ==1,][training_positives,]),
  rownames(SerengetiFW[SerengetiFW$interaction ==0,][training_negatives,])
))
training <- SerengetiFW[training_id, ]

# Set up the bayesian glmm ------------------------------------------------
# extract predictors
training <- select(training,  -Predator, -Prey, -Order.prey, 
                   -Herbivore.predator,  -Herbivore.prey) %>%
  mutate(Order.predator = factor(Order.predator, 
                                 levels = unique(FuncTraits$Order)))

get_prior(brms_form, data = training)

model_priors <- c(
  prior(normal(0, 1), class = "Intercept"),
  prior(normal(0, 1), class = "b"),
  prior(cauchy(0, 5), class = "sd")
)

prior_predictions <- brm(formula = brms_form,
                         data = training,
                         prior = model_priors,
                         sample_prior = "only", inits = "0")

SerengetiModel <- brm(formula = brms_form,
                      data = training,
                      prior = model_priors, sample_prior = "no", 
                      cores = 4, backend = "cmdstan", threads = 4,
                      inits = "0", iter = 2000)


# save the model
saveRDS(SerengetiModel, 
     file = paste0("~/OneDrive/Chapt-NicheSpace/models/SerengetiModel_brms.rds"))
