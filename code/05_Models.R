# Step 05: Train a predictive model on each food web
# For each model (one per food web):
# 1. Load and standardize training data
# 1.1 load the data
# 1.2 standardize continuous predictors
# 1.3 take 70% of all interactions + equal number of non-interactions
# 2. Set up bayesian GLMM with brms
# 3. Train the model and save it on OneDrive (too big for GitHub)
# 4. Train a Boosted Regression Tree

rm(list = ls())
set.seed(16)
library(dplyr)
library(brms)
library(gbm)
source("code/functions.R")
FuncTraits <- read.csv("data/cleaned/SpeciesTraitsFull.csv", row.names = 1)
brms_form <- bf(interaction ~ 1 + 
                  (Omnivore.predator + Carnivore.predator + Habitat_breadth.predator + BM.predator + Longevity.predator + ClutchSize.predator +
                     Omnivore.prey + Carnivore.prey + Habitat_breadth.prey + BM.prey + Longevity.prey + ClutchSize.prey + 
                     ActivityTime.match + Habitat.match + BM.match) + 
                  (1 + (Omnivore.predator + Carnivore.predator + Habitat_breadth.predator + BM.predator + Longevity.predator + ClutchSize.predator 
                        + Omnivore.prey + Carnivore.prey + Habitat_breadth.prey + BM.prey + Longevity.prey + ClutchSize.prey + 
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

# keep 30% of the food web for validation controlling for prevalence
Npos = round(0.3 * sum(EuroMW$interaction == 1))
Nneg = round(0.3 * sum(EuroMW$interaction == 0))
testing_positives <- sample(seq_len(sum(EuroMW$interaction)), size = Npos)
testing_negatives <- sample(seq_len(sum(EuroMW$interaction==0)), size = Nneg)
testing_id_euro <- as.numeric(c(
  rownames(EuroMW[EuroMW$interaction ==1,][testing_positives,]),
  rownames(EuroMW[EuroMW$interaction ==0,][testing_negatives,])
))

# take the remaining positives and as many negatives for training
training <- EuroMW[-testing_id_euro, ]
N = sum(training$interaction)
training_negatives <- sample(seq_len(sum(training$interaction==0)), size = N)
training_id_euro <- as.numeric(c(
  rownames(training[training$interaction ==1,]),
  rownames(training[training$interaction ==0,][training_negatives,])
))
training <- EuroMW[training_id_euro, ]

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
                         sample_prior = "only", init = "0")
EuroModel <- brm(formula = brms_form,
                 data = training,
                 prior = model_priors, sample_prior = "no", 
                    cores = 4, backend = "cmdstan", threads = 4,
                    init = "0", iter = 2000)

# save the model on OneDrive (too big for Github...)
saveRDS(EuroModel,
        file = paste0("~/OneDrive/Chapt-NicheSpace/models/EuroModel_brms.rds"))

# train a BRT



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

# keep 30% of the food web for validation controlling for prevalence
Npos = round(0.3 * sum(HighArcticFW$interaction == 1))
Nneg = round(0.3 * sum(HighArcticFW$interaction == 0))
testing_positives <- sample(seq_len(sum(HighArcticFW$interaction)), size = Npos)
testing_negatives <- sample(seq_len(sum(HighArcticFW$interaction==0)), size = Nneg)
testing_id_higharctic <- as.numeric(c(
  rownames(HighArcticFW[HighArcticFW$interaction ==1,][testing_positives,]),
  rownames(HighArcticFW[HighArcticFW$interaction ==0,][testing_negatives,])
))

# take the remaining positives and as many negatives for training
training <- HighArcticFW[-testing_id_higharctic, ]
N = sum(training$interaction)
training_negatives <- sample(seq_len(sum(training$interaction==0)), size = N)
training_id_higharctic <- as.numeric(c(
  rownames(training[training$interaction ==1,]),
  rownames(training[training$interaction ==0,][training_negatives,])
))
training <- HighArcticFW[training_id_higharctic, ]

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
                         sample_prior = "only", init = "0")

ArcticModel <- brm(formula = brms_form,
                   data = training,
                   prior = model_priors, sample_prior = "no", 
                   cores = 4, backend = "cmdstan", threads = 4,
                   init = "0", iter = 2000)

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

# keep 30% of the food web for validation controlling for prevalence
Npos = round(0.3 * sum(PyreneesFW$interaction == 1))
Nneg = round(0.3 * sum(PyreneesFW$interaction == 0))
testing_positives <- sample(seq_len(sum(PyreneesFW$interaction)), size = Npos)
testing_negatives <- sample(seq_len(sum(PyreneesFW$interaction==0)), size = Nneg)
testing_id_pyrenees <- as.numeric(c(
  rownames(PyreneesFW[PyreneesFW$interaction ==1,][testing_positives,]),
  rownames(PyreneesFW[PyreneesFW$interaction ==0,][testing_negatives,])
))

# take the remaining positives and as many negatives for training
training <- PyreneesFW[-testing_id_pyrenees, ]
N = sum(training$interaction)
training_negatives <- sample(seq_len(sum(training$interaction==0)), size = N)
training_id_pyrenees <- as.numeric(c(
  rownames(training[training$interaction ==1,]),
  rownames(training[training$interaction ==0,][training_negatives,])
))
training <- PyreneesFW[training_id_pyrenees, ]

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
                         sample_prior = "only", init = "0")

PyreneesModel <- brm(formula = brms_form,
                     data = training,
                     prior = model_priors, sample_prior = "no", 
                     cores = 4, backend = "cmdstan", threads = 4,
                     init = "0", iter = 2000)

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

# keep 30% of the food web for validation controlling for prevalence
Npos = round(0.3 * sum(SerengetiFW$interaction == 1))
Nneg = round(0.3 * sum(SerengetiFW$interaction == 0))
testing_positives <- sample(seq_len(sum(SerengetiFW$interaction)), size = Npos)
testing_negatives <- sample(seq_len(sum(SerengetiFW$interaction==0)), size = Nneg)
testing_id_serengeti <- as.numeric(c(
  rownames(SerengetiFW[SerengetiFW$interaction ==1,][testing_positives,]),
  rownames(SerengetiFW[SerengetiFW$interaction ==0,][testing_negatives,])
))

# take the remaining positives and as many negatives for training
training <- SerengetiFW[-testing_id_serengeti, ]
N = sum(training$interaction)
training_negatives <- sample(seq_len(sum(training$interaction==0)), size = N)
training_id_serengeti <- as.numeric(c(
  rownames(training[training$interaction ==1,]),
  rownames(training[training$interaction ==0,][training_negatives,])
))
training <- SerengetiFW[training_id_serengeti, ]

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
                         sample_prior = "only", init = "0")

SerengetiModel <- brm(formula = brms_form,
                      data = training,
                      prior = model_priors, sample_prior = "no", 
                      cores = 4, backend = "cmdstan", threads = 4,
                      init = "0", iter = 2000)

# save the model
saveRDS(SerengetiModel, 
     file = paste0("~/OneDrive/Chapt-NicheSpace/models/SerengetiModel_brms.rds"))

# save the splits
save(testing_id_euro, training_id_euro,
     testing_id_pyrenees, training_id_pyrenees,
     testing_id_higharctic, training_id_higharctic,
     testing_id_serengeti, training_id_serengeti, 
     file = "data/checkpoints/train_test_splits.RData")
