rm(list = ls())

set.seed(16)
library(bayesplot)
library(ROCR)
library(greta)
library(dplyr)
library(ggplot2)
library(tidyr)
source("code/functions.R")
load("data/models/GLMM_10032022.RData")

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

predictor_names <- c("Intercept", 
                     colnames(select(validation, 
                                     -Predator, -Prey, -Order.predator, -Order.prey,
                                     -Herbivore.predator, -Herbivore.prey, -interaction)))

# Check convergence -------------------------------------------------------
mcmc_trace(GLMM, regex_pars = "global_coef_mean")
mcmc_trace(GLMM, regex_pars = "global_coef_sd")
rhat <- coda::gelman.diag(GLMM)$psrf

# Plot coefficients -------------------------------------------------------
global_coef_mean <- data.frame(mean = summary(GLMM)$statistics[c(1:16),1])
global_coef_mean$sd <- summary(GLMM)$statistics[c(17:32),1]
global_coef_mean$predictor <- factor(predictor_names, 
                                     levels = rev(predictor_names))
model_coef <- calculate(coef, values = GLMM)
model_coef <- model_coef$`11`
coef_order <- data.frame(t(matrix(apply(model_coef, MARGIN = 2, mean), ncol = 48, nrow = 16)))
colnames(coef_order) <- predictor_names
coef_order <- pivot_longer(coef_order, cols = everything(), names_to = "predictor", values_to = "mean")
coef_order$predictor <- factor(coef_order$predictor, 
                               levels = rev(predictor_names))
coef_order_sd <- data.frame(t(matrix(apply(model_coef, MARGIN = 2, sd), ncol = 48, nrow = 16)))
colnames(coef_order_sd) <- predictor_names
coef_order_sd <- pivot_longer(coef_order_sd, cols = everything(), names_to = "predictor", values_to = "sd")
coef_order_sd$predictor <- factor(coef_order_sd$predictor, 
                               levels = rev(predictor_names))

ggplot()+
  geom_point(data = coef_order, aes(x = mean, y = predictor), alpha = 0.5, position = position_jitter(height = 0.4)) +
  geom_pointrange(data = global_coef_mean, aes(xmin = mean-1.96*sd, x = mean, xmax = mean+1.96*sd, y = predictor), colour = "red")

# Predict Validation data -------------------------------------------------
# Get a validation dataset with the same prevalence than the MW
prev <- mean(EuroMW$interaction)
validation <- EuroMW[-training_id, ]
validation_interactions <- validation[validation$interaction == 1,]
N_noninteractions <- round((sum(validation$interaction) * (1-prev))/prev) # how many non-interactions
validation_noninteractions <- slice_sample(validation[validation$interaction == 0,],
                                           n = N_noninteractions)
validation <- rbind(validation_interactions, validation_noninteractions)
predictors <- select(validation, 
                     -Predator, -Prey, -Order.predator, -Order.prey,
                     -Herbivore.predator, -Herbivore.prey, -interaction)

predictors <- cbind(rep(1, nrow(validation)), predictors) %>% as_data()
predator_order <- as.numeric(factor(validation$Order.pred, levels = unique(FuncTraits$Order)))

linear_predictor <- rowSums(predictors * t(coef[,predator_order]))

p <- ilogit(linear_predictor)

predictions <- calculate(p, values = GLMM, nsim = 100)
predictions <- colMeans(predictions[[1]])
auc <- performance(prediction(predictions, validation$interaction), "auc")@y.values[[1]]
aucpr <- performance(prediction(predictions, validation$interaction), "aucpr")@y.values[[1]]
