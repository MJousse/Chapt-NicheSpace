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
                     -Herbivore.predator, -Herbivore.prey, -interaction,
                     -ClutchSize.prey, -ClutchSize.predator)

# add a column for the intercept
predictors <- cbind(rep(1, nrow(training)), predictors) %>% as_data()
y <- as_data(training$interaction)

# get levels for the random effect: predator order
predator_order <- as.numeric(factor(training$Order.pred, levels = unique(FuncTraits$Order)))

# number of predictors and orders
norder <- length(unique(FuncTraits$Order))
ntraits <- ncol(predictors)

# common effect
global_coef_mean <- normal(0,1, ntraits)
global_coef_sd <- cauchy(0, 5, truncation = c(0.001, Inf), dim = ntraits)

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
GLMM <- mcmc(m, n_samples = 15000, warmup = 10000, chains = 4)

# save the model
save(GLMM, global_coef_mean, global_coef_sd, training_id, training, coef,  
     file =  paste0("data/models/GLMM_", format(Sys.Date(), "%d%m%Y"), ".RData"))
