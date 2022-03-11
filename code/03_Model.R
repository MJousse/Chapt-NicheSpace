rm(list = ls())
set.seed(16)
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

# GLMM --------------------------------------------------------------------
predictors <- select(training, 
                     -Predator, -Prey, -Order.predator, -Order.prey,
                     -Herbivore.predator, -Herbivore.prey, -interaction)

predictors <- cbind(rep(1, nrow(training)), predictors) %>% as_data()
y <- as_data(training$interaction)
predator_order <- as.numeric(factor(training$Order.pred, levels = unique(FuncTraits$Order)))

norder <- length(unique(FuncTraits$Order))
ntraits <- ncol(predictors)

# common effect
global_coef_mean <- normal(0,1, ntraits)
global_coef_sd <- cauchy(0, 5, truncation = c(0.001, Inf), dim = ntraits)
coef <- normal(global_coef_mean, global_coef_sd, dim = c(ntraits))
for (i in c(2:norder)){
  coef <- cbind(coef, normal(global_coef_mean, global_coef_sd, dim = c(ntraits)))
}

# parameters of the base model are a function of the order of the predator
linear_predictor <- rowSums(predictors * t(coef[,predator_order]))

# ilogit of linear predictor
p <- ilogit(linear_predictor)

# a single observation means our data are bernoulli distributed
distribution(y) <- bernoulli(p)

m <- model(global_coef_mean, global_coef_sd)

GLMM <- mcmc(m, n_samples = 25000, warmup = 10000, chains = 4)

save(GLMM, global_coef_mean, global_coef_sd, train_ind, training, coef,  
     file =  paste0("data/models/GLMM_", format(Sys.Date(), "%d%m%Y"), ".RData"))

GLMMextra <- extra_samples(GLMM, n_samples = 2000)
