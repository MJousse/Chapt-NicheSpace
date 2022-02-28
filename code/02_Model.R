rm(list = ls())

library(dplyr)
library(greta)
source("code/functions.R")


# Load data and standardize -----------------------------------------------
EuroInteractions <- read.csv("data/cleaned/EuroFW.csv", row.names = 1)
FuncTraits <- read.csv("data/cleaned/SpeciesTraitsFull.csv", row.names = 1)
EuroSpecies <- read.csv("data/cleaned/EuroMWTaxo.csv", row.names = 1) %>%
  distinct() %>%
  filter(!is.na(Species))

# Transform trait into predictors for every species pair ------------------
EuroMW <- expand.grid(EuroSpecies$Species, EuroSpecies$Species) 
colnames(EuroMW) <- c("Predator", "Prey")

EuroMW <- left_join(EuroMW, FuncTraits, by = c("Prey" = "Species")) %>%
  left_join(FuncTraits, by = c("Predator" = "Species"))

EuroMW <- traits2predictors(EuroMW)

EuroMW <- mutate_at(EuroMW, vars(Habitat_breadth.predator, BM.predator:ClutchSize.predator,
                                 Habitat_breadth.prey, BM.prey:ClutchSize.prey, 
                                 Habitat.match:BM.match), scale2)

EuroInteractions$interaction <- 1

EuroMW <- left_join(EuroMW, EuroInteractions)
EuroMW$interaction[is.na(EuroMW$interaction)] <- 0

# Divide into training and testing data -----------------------------------
training_size = 0.7 * nrow(EuroMW)
train_ind <- sample(seq_len(nrow(EuroMW)), size = training_size)
training <- EuroMW[train_ind, ]
validation <- EuroMW[-train_ind, ]

# GLM ---------------------------------------------------------------------
predictors <- select(training, Trophic_level.predator, Habitat_breadth.predator, BM.predator, Longevity.predator,
                     ClutchSize.predator, Trophic_level.prey, Habitat_breadth.prey, BM.prey,
                     Longevity.prey, ClutchSize.prey, ActivityTime.match, Habitat.match, BM.match) %>% as_data()
y <- as_data(training$interaction)
predator_id <- as.numeric(factor(training$order.pred, levels = unique(FuncTraits$Order)))
norder <- length(unique(FuncTraits$Order))

global_alpha <- normal(0, 10, dim = 1)
global_alpha_sd <- uniform(0, 10, dim = 1) 
alpha <- normal(global_alpha, global_alpha_sd, dim = norders)

global_betas <- normal(0, 10, dim = 18)
global_betas_sd <- uniform(0, 10, dim = 18)
beta <- normal(global_betas, global_betas_sd, dim = c(18, norders))

linear_predictor <- alpha[predator_id] + predictors %*% beta[,predator_id]

# ilogit of linear predictor
p <- ilogit(linear_predictor)

# a single observation means our data are bernoulli distributed
distribution(y) <- bernoulli(p)

