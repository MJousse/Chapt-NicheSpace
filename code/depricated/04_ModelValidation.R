# Step 4: Model validation
# 1. Look at model convergence
# 2. Plot model coefficients
# 3. Measure performance on a European Metaweb validation dataset

rm(list = ls())
set.seed(16)
library(bayesplot)
library(ROCR)
library(greta)
library(dplyr)
library(ggplot2)
library(tidyr)
source("code/functions.R")
load("data/models/GLMM_21032022.RData")

# Check model convergence -------------------------------------------------
mcmc_trace(GLMM, regex_pars = "global_coef_mean")
mcmc_trace(GLMM, regex_pars = "global_coef_sd")
coda::gelman.diag(GLMM)$psrf

# Prepare validation dataset ----------------------------------------------
# load data
EuroInteractions <- read.csv("data/cleaned/EuroFWadults.csv", row.names = 1)
FuncTraits <- read.csv("data/cleaned/SpeciesTraitsFull.csv", row.names = 1)
EuroSpecies <- read.csv("data/cleaned/EuroMWTaxo.csv", row.names = 1) %>%
  distinct() %>%
  filter(!is.na(Species))

# transform traits into predictors
EuroMW <- get_predictors(EuroSpecies$Species, FuncTraits)
EuroMW <- mutate_at(EuroMW, vars(Habitat_breadth.predator, BM.predator:ClutchSize.predator,
                                 Habitat_breadth.prey, BM.prey:ClutchSize.prey, 
                                 Habitat.match:BM.match), scale2)

# add response
EuroInteractions$interaction <- 1
EuroMW <- left_join(EuroMW, EuroInteractions)
EuroMW$interaction[is.na(EuroMW$interaction)] <- 0

# extract predictor name
predictor_names <- c("Intercept", 
                     colnames(select(training, 
                                     -Predator, -Prey, -Order.predator, -Order.prey,
                                     -Herbivore.predator, -Herbivore.prey, -interaction)))

# prepare a validation dataset with the same prevalence than the metaweb
prev <- mean(EuroMW$interaction)
validation <- EuroMW[-training_id, ]
validation_interactions <- validation[validation$interaction == 1,] # use all interactions
N_noninteractions <- round((sum(validation$interaction) * (1-prev))/prev) # how many non-interactions
validation_noninteractions <- slice_sample(validation[validation$interaction == 0,],
                                           n = N_noninteractions) # add non-interactions
validation <- rbind(validation_interactions, validation_noninteractions)

# Plot model coefficients -------------------------------------------------
# global coefficient mean and sd
global_coef_mean <- data.frame(mean = summary(GLMM)$statistics[c(1:14),1])
global_coef_mean$sd <- summary(GLMM)$statistics[c(15:28),1]
global_coef_mean$predictor <- factor(predictor_names, 
                                     levels = rev(predictor_names))

# calculate order-specific coefficients
model_coef <- calculate(coef, values = GLMM)
model_coef <- model_coef$`11`
coef_order <- data.frame(t(matrix(apply(model_coef, MARGIN = 2, mean), ncol = 48, nrow = 14)))
colnames(coef_order) <- predictor_names
coef_order <- pivot_longer(coef_order, cols = everything(), names_to = "predictor", values_to = "mean")
coef_order$predictor <- factor(coef_order$predictor, 
                               levels = rev(predictor_names))
coef_order_sd <- data.frame(t(matrix(apply(model_coef, MARGIN = 2, sd), ncol = 48, nrow = 14)))
colnames(coef_order_sd) <- predictor_names
coef_order_sd <- pivot_longer(coef_order_sd, cols = everything(), names_to = "predictor", values_to = "sd")
coef_order_sd$predictor <- factor(coef_order_sd$predictor, 
                               levels = rev(predictor_names))

# plot the coefficients
ggplot()+
  geom_point(data = coef_order, aes(x = mean, y = predictor), alpha = 0.5, position = position_jitter(height = 0.4)) +
  geom_pointrange(data = global_coef_mean, aes(xmin = mean-1.96*sd, x = mean, xmax = mean+1.96*sd, y = predictor), colour = "red")

# Predict validation dataset ----------------------------------------------
# make predictions
predictions <- make_predictions(validation[,-ncol(validation)], unique(FuncTraits$Order), coef, GLMM)

# calculate roc-auc and pr-auc
auc <- performance(prediction(predictions$prediction, validation$interaction), "auc")@y.values[[1]]
aucpr <- performance(prediction(predictions$prediction, validation$interaction), "aucpr")@y.values[[1]]
