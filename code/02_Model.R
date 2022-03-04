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

# # Divide into training and testing data -----------------------------------
# training_size = 0.01 * nrow(EuroMW)
# train_ind <- sample(seq_len(nrow(EuroMW)), size = training_size)
# training <- EuroMW[train_ind, ]
# validation <- EuroMW[-train_ind, ]
# 
# # GLMM --------------------------------------------------------------------
# TL <- select(training2, Omnivore.predator, Carnivore.predator, 
#              Omnivore.prey, Carnivore.prey)
# predictors2 <- select(training2, Habitat_breadth.predator, BM.predator, 
#                       Longevity.predator, ClutchSize.predator, 
#                       Habitat_breadth.prey, BM.prey, Longevity.prey, ClutchSize.prey,
#                       ActivityTime.match, Habitat.match, BM.match)
# 
# predictors <- select(training, 
#                      -Predator, -Prey, -Order.predator, -Order.prey,
#                      -Herbivore.predator, -Herbivore.prey, -interaction)
# predictors <- cbind(rep(1, nrow(training)), predictors) %>% as_data()
# y <- as_data(training$interaction)
# predator_order <- as.numeric(factor(training$Order.pred, levels = unique(FuncTraits$Order)))
# 
# norder <- length(unique(FuncTraits$Order))
# ntraits <- ncol(predictors)
# 
# # common effect
# common_coef1 <- normal(0, 1, dim = ntraits)
# 
# #random effect
# 
# order_coef_sd <- lognormal(0, 1)
# order_coef <- normal(0, order_coef_sd, dim = c(ntraits, norder-1))
# order_coef <- cbind(rep(0,ntraits), order_coef)
# 
# # parameters of the base model are a function of the order of the predator
# linear_predictor <- predictors %*% common_coef + rowSums(predictors * t(order_coef[,predator_order]))
# 
# # ilogit of linear predictor
# p <- ilogit(linear_predictor)
# 
# # a single observation means our data are bernoulli distributed
# distribution(y) <- bernoulli(p)
# 
# m <- model(common_coef, order_coef_sd, order_coef)
# 
# GLMM <- mcmc(m, n_samples = 1000, warmup = 1000, chains = 4)
# 
# library(bayesplot)
# mcmc_trace(GLMM, pars = paste0("order_coef[3,", c(1:10), "]"))
# mcmc_trace(GLMM, regex_pars = "common_coef")
# coef <- calculate(coef, values = GLMM2, nsim = 100)
# mcmc_trace(GLMM2, pars = paste0("order_coef[",c(1:16),",5]"))
# 

# GLMM test 2 -------------------------------------------------------------
training_size = 0.005 * nrow(EuroMW)
train_ind <- sample(seq_len(nrow(EuroMW)), size = training_size)
training2 <- EuroMW[train_ind, ]
validation2 <- EuroMW[-train_ind, ]

predictors2 <- select(training2, 
                     -Predator, -Prey, -Order.predator, -Order.prey,
                     -Herbivore.predator, -Herbivore.prey, -interaction)

predictors2 <- cbind(rep(1, nrow(training2)), predictors2) %>% as_data()
y2 <- as_data(training2$interaction)
predator_order2 <- as.numeric(factor(training2$Order.pred, levels = unique(FuncTraits$Order)))

norder <- length(unique(FuncTraits$Order))
ntraits <- ncol(predictors2)

# common effect
global_coef_mean <- normal(0,1, ntraits)
global_coef_sd <- uniform(0, 10, dim = ntraits)
coef <- normal(global_coef_mean, global_coef_sd, dim = c(ntraits))
for (i in c(2:norder)){
  coef <- cbind(coef, normal(global_coef_mean, global_coef_sd, dim = c(ntraits)))
}

# parameters of the base model are a function of the order of the predator
linear_predictor2 <- rowSums(predictors2 * t(coef[,predator_order2]))

# ilogit of linear predictor
p2 <- ilogit(linear_predictor2)

# a single observation means our data are bernoulli distributed
distribution(y2) <- bernoulli(p2)

m2 <- model(global_coef_mean, global_coef_sd)

GLMM2 <- mcmc(m2, n_samples = 5000, warmup = 5000, chains = 4)

library(bayesplot)
mcmc_trace(GLMM2, regex_pars = "global_coef_mean")
library(ROCR)
predict2 <- calculate(p2, values = GLMM2, nsim = 100)
prediction_m2 <- colMeans(predict2[[1]])
m2_auc <- performance(prediction(prediction_m2, training2$interaction), "auc")@y.values[[1]]
m2_aucpr <- performance(prediction(prediction_m2, training2$interaction), "aucpr")@y.values[[1]]

predict1 <- calculate(p, values = GLMM, nsim = 100)
prediction_m1 <- colMeans(predict1[[1]])
m1_auc <- performance(prediction(prediction_m1, training$interaction), "auc")@y.values[[1]]
m1_aucpr <- performance(prediction(prediction_m1, training$interaction), "aucpr")@y.values[[1]]
coda::gelman.diag(GLMM2)$psrf




mcmc_trace(GLMM2, pars = paste0("coef[3,", c(1:10), "]"))
mcmc_trace(GLMM2, pars = paste0("coef[",c(1:16),",5]"))

mcmc_intervals(GLMM, pars = "coef")


