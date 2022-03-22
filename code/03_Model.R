# Step 03: Train the predictive model on the European metaweb
# 1. Prepare the training dataset
# 2. Set up the GLMM with greta package
# 3. Train the model

rm(list = ls())
set.seed(16)
library(dplyr)
library(greta)
source("code/functions.R")

# Prepare training dataset ------------------------------------------------
# load data and standardize
EuroInteractions <- read.csv("data/cleaned/EuroFWadults.csv", row.names = 1)
FuncTraits <- read.csv("data/cleaned/SpeciesTraitsFull.csv", row.names = 1)
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
predictors <- select(training, 
                     -Predator, -Prey, -Order.predator, -Order.prey,
                     -Herbivore.predator, -Herbivore.prey, -interaction)

# add a column for the intercept
predictors <- cbind(rep(1, nrow(training)), predictors) %>% as_data()
y <- as_data(training$interaction)

# get levels for the random effect: predator order
predator_order <- as.numeric(factor(training$Order.pred, levels = unique(FuncTraits$Order)))

# number of predictors and orders
norder <- length(unique(FuncTraits$Order))
ntraits <- ncol(predictors)

# common effect
global_coef_mean <- normal(0, 1, ntraits)
global_coef_sd <- cauchy(0, 3, truncation = c(0.001, Inf), dim = ntraits)

# coef for each order is drawn from global mean and sd
coef <- normal(global_coef_mean, global_coef_sd, dim = c(ntraits))
for (i in c(2:norder)){
  coef <- cbind(coef, normal(global_coef_mean, global_coef_sd, dim = c(ntraits)))
}

# linear prediction
linear_predicton <- rowSums(predictors * t(coef[,predator_order]))

# ilogit of linear prediction
p <- ilogit(linear_predicton)

# interaction as Bernoulli trial
distribution(y) <- bernoulli(p)

# Train the model with greta ----------------------------------------------
m <- model(global_coef_mean, global_coef_sd)
EuroModel <- mcmc(m, n_samples = 20000, warmup = 10000, chains = 4)

# save the model on OneDrive (too big for Github...)
save(EuroModel, global_coef_mean, global_coef_sd, training_id, training, coef,  
     file =  paste0("~/OneDrive/Chapt-NicheSpace/models/EuroModel_", 
                    format(Sys.Date(), "%d%m%Y"), ".RData"))

# -------------------------------------------------------------------------

# Prepare training dataset ------------------------------------------------
# load data and standardize
HighArcticInteractions <- read.csv("data/cleaned/HighArcticFW.csv", row.names = 1)
HighArcticSpecies <- read.csv("data/cleaned/HighArcticFW.csv", row.names = 1) %>%
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
predictors <- select(training, 
                     -Predator, -Prey, -Order.predator, -Order.prey,
                     -Herbivore.predator, -Herbivore.prey, -interaction)

# add a column for the intercept
predictors <- cbind(rep(1, nrow(training)), predictors) %>% as_data()
y <- as_data(training$interaction)

# get levels for the random effect: predator order
predator_order <- as.numeric(factor(training$Order.pred, levels = unique(FuncTraits$Order)))

# number of predictors and orders
norder <- length(unique(FuncTraits$Order))
ntraits <- ncol(predictors)

# common effect
global_coef_mean <- normal(0, 1, ntraits)
global_coef_sd <- cauchy(0, 3, truncation = c(0.001, Inf), dim = ntraits)

# coef for each order is drawn from global mean and sd
coef <- normal(global_coef_mean, global_coef_sd, dim = c(ntraits))
for (i in c(2:norder)){
  coef <- cbind(coef, normal(global_coef_mean, global_coef_sd, dim = c(ntraits)))
}

# linear prediction
linear_predicton <- rowSums(predictors * t(coef[,predator_order]))

# ilogit of linear prediction
p <- ilogit(linear_predicton)

# interaction as Bernoulli trial
distribution(y) <- bernoulli(p)

# Train the model with greta ----------------------------------------------
m <- model(global_coef_mean, global_coef_sd)
ArcticModel <- mcmc(m, n_samples = 15000, warmup = 10000, chains = 4)

# save the model
save(ArcticModel, global_coef_mean, global_coef_sd, training_id, training, coef,  
     file =  paste0("~/OneDrive/Chapt-NicheSpace/models/ArcticModel_", 
                    format(Sys.Date(), "%d%m%Y"), ".RData"))

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
predictors <- select(training, 
                     -Predator, -Prey, -Order.predator, -Order.prey,
                     -Herbivore.predator, -Herbivore.prey, -interaction)

# add a column for the intercept
predictors <- cbind(rep(1, nrow(training)), predictors) %>% as_data()
y <- as_data(training$interaction)

# get levels for the random effect: predator order
predator_order <- as.numeric(factor(training$Order.pred, levels = unique(FuncTraits$Order)))

# number of predictors and orders
norder <- length(unique(FuncTraits$Order))
ntraits <- ncol(predictors)

# common effect
global_coef_mean <- normal(0, 1, ntraits)
global_coef_sd <- cauchy(0, 3, truncation = c(0.001, Inf), dim = ntraits)

# coef for each order is drawn from global mean and sd
coef <- normal(global_coef_mean, global_coef_sd, dim = c(ntraits))
for (i in c(2:norder)){
  coef <- cbind(coef, normal(global_coef_mean, global_coef_sd, dim = c(ntraits)))
}

# linear prediction
linear_predicton <- rowSums(predictors * t(coef[,predator_order]))

# ilogit of linear prediction
p <- ilogit(linear_predicton)

# interaction as Bernoulli trial
distribution(y) <- bernoulli(p)

# Train the model with greta ----------------------------------------------
m <- model(global_coef_mean, global_coef_sd)
PyreneesModel <- mcmc(m, n_samples = 20000, warmup = 10000, chains = 4)

# save the model
save(PyreneesModel, global_coef_mean, global_coef_sd, training_id, training, coef,  
     file =  paste0("~/OneDrive/Chapt-NicheSpace/models/PyreneesModel_", 
                    format(Sys.Date(), "%d%m%Y"), ".RData"))

# Prepare training dataset ------------------------------------------------
# load data and standardize
SerengetiInteractions <- read.csv("data/cleaned/SerengetiFW.csv", row.names = 1)
SerengetiSpecies <- read.csv("data/cleaned/SerengetiFWTaxo.csv", row.names = 1) %>%
  distinct() %>%
  filter(!is.na(Species))

# transform trait into predictors for every species pair
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
N = round(0.7 * sum(EuroMW$interaction))
training_positives <- sample(seq_len(sum(SerengetiFW$interaction)), size = N)
training_negatives <- sample(seq_len(sum(SerengetiFW$interaction==0)), size = N)
training_id <- as.numeric(c(
  rownames(SerengetiFW[SerengetiFW$interaction ==1,][training_positives,]),
  rownames(SerengetiFW[SerengetiFW$interaction ==0,][training_negatives,])
))
training <- SerengetiFW[training_id, ]

# Set up the bayesian glmm ------------------------------------------------
# extract predictors
predictors <- select(training, 
                     -Predator, -Prey, -Order.predator, -Order.prey,
                     -Herbivore.predator, -Herbivore.prey, -interaction)

# add a column for the intercept
predictors <- cbind(rep(1, nrow(training)), predictors) %>% as_data()
y <- as_data(training$interaction)

# get levels for the random effect: predator order
predator_order <- as.numeric(factor(training$Order.pred, levels = unique(FuncTraits$Order)))

# number of predictors and orders
norder <- length(unique(FuncTraits$Order))
ntraits <- ncol(predictors)

# common effect
global_coef_mean <- normal(0, 1, ntraits)
global_coef_sd <- cauchy(0, 3, truncation = c(0.001, Inf), dim = ntraits)

# coef for each order is drawn from global mean and sd
coef <- normal(global_coef_mean, global_coef_sd, dim = c(ntraits))
for (i in c(2:norder)){
  coef <- cbind(coef, normal(global_coef_mean, global_coef_sd, dim = c(ntraits)))
}

# linear prediction
linear_predicton <- rowSums(predictors * t(coef[,predator_order]))

# ilogit of linear prediction
p <- ilogit(linear_predicton)

# interaction as Bernoulli trial
distribution(y) <- bernoulli(p)

# Train the model with greta ----------------------------------------------
m <- model(global_coef_mean, global_coef_sd)
SerengetiModel <- mcmc(m, n_samples = 15000, warmup = 10000, chains = 4)

# save the model
save(SerengetiModel, global_coef_mean, global_coef_sd, training_id, training, coef,  
     file =  paste0("~/OneDrive/Chapt-NicheSpace/models/SerengetiModel_", 
                    format(Sys.Date(), "%d%m%Y"), ".RData"))
