# Step 06: Check model behavior and compare coefficients
# For each model (one per food web):
# 1. Load the models
# 2. Check convergence with Rhat and trace plots
# 3. Posterior predictive check (predict realistic number of interactions)
# 4. Check R2
# 5. Compare coefficiens: *creates figure XX of manuscript*

rm(list = ls())
library(tidybayes)
library(bayesplot)
library(brms)
library(dplyr)
library(ggplot2)
library(tibble)

# Load the models ---------------------------------------------------------
ArcticModel <- readRDS("~/OneDrive/Chapt-NicheSpace/models/ArcticModel_brms.rds")
EuropeModel <- readRDS("~/OneDrive/Chapt-NicheSpace/models/EuroModel_brms.rds")
PyreneesModel <- readRDS("~/OneDrive/Chapt-NicheSpace/models/PyreneesModel_brms.rds")
SerengetiModel <- readRDS("~/OneDrive/Chapt-NicheSpace/models/SerengetiModel_brms.rds")

load("data/checkpoints/train_test_splits.RData")

# Check convergence -------------------------------------------------------
color_scheme_set("viridis")
summary(ArcticModel)
hist(rhat(ArcticModel), xlab = "Rhat", main = "Distribution of rank-normalized potential\nscale reduction factors for the Arctic model")
mcmc_plot(ArcticModel, variable = c('b_'), regex = T, type = 'rank_overlay') + theme(text = element_text(size = 8))
ggsave("figures/SI/trank_arctic.png", width = 8, height = 8)
summary(EuropeModel)
hist(rhat(EuropeModel), xlab = "Rhat", main = "Distribution of rank-normalized potential\nscale reduction factors for the Europe model")
mcmc_plot(EuropeModel, variable = c('b_'), regex = T, type = 'rank_overlay') + theme(text = element_text(size = 8))
ggsave("figures/SI/trank_europe.png", width = 8, height = 8)
summary(PyreneesModel)
hist(rhat(PyreneesModel), xlab = "Rhat", main = "Distribution of rank-normalized potential\nscale reduction factors for the Pyrenees model")
mcmc_plot(PyreneesModel, pars = c('b_'), type = 'rank_overlay') + theme(text = element_text(size = 8))
ggsave("figures/SI/trank_pyrenees.png", width = 8, height = 8)
summary(SerengetiModel)
hist(rhat(SerengetiModel), xlab = "Rhat", main = "Distribution of rank-normalized potential\nscale reduction factors for the Serengeti model")
mcmc_plot(SerengetiModel, pars = c('b_'), type = 'rank_overlay') + theme(text = element_text(size = 8))
ggsave("figures/SI/trank_serengeti.png", width = 8, height = 8)

# Posterior predictive checks ---------------------------------------------
color_scheme_set("purple")
pp_check(ArcticModel, ndraws = 100, type = "bars_grouped", group = "Order.predator") + theme(text = element_text(size = 8))
ggsave("figures/SI/ppcheck_arctic.png", width = 8, height = 8)
pp_check(EuropeModel, ndraws = 100, type = "bars_grouped", group = "Order.predator") + theme(text = element_text(size = 8))
ggsave("figures/SI/ppcheck_europe.png", width = 8, height = 8)
pp_check(PyreneesModel, ndraws = 100, type = "bars_grouped", group = "Order.predator") + theme(text = element_text(size = 8))
ggsave("figures/SI/ppcheck_pyrenees.png", width = 8, height = 8)
pp_check(SerengetiModel, ndraws = 100, type = "bars_grouped", group = "Order.predator") + theme(text = element_text(size = 7))
ggsave("figures/SI/ppcheck_serengeti.png", width = 8, height = 8)

# Bayes R2 ----------------------------------------------------------------
R2_arctic <- bayes_R2(ArcticModel)
R2_europe <- bayes_R2(EuropeModel)
R2_pyrenees <- bayes_R2(PyreneesModel)
R2_serengeti <- bayes_R2(SerengetiModel)

# Compare coefficients ----------------------------------------------------
# fixed effects
predictors <- c("Intercept", "Omnivore.predator", "Carnivore.predator", "Habitat_breadth.predator",
       "BM.predator", "Longevity.predator", "ClutchSize.predator", "Omnivore.prey", "Carnivore.prey",
       "Habitat_breadth.prey", "BM.prey", "Longevity.prey", "ClutchSize.prey", "ActivityTime.match",
       "Habitat.match", "BM.match")
foodwebs <- c("Arctic", "Europe", "Pyrenees", "Serengeti")
fixed_effects <- rbind(
  fixef(ArcticModel, pars = predictors),
  fixef(EuropeModel, pars = predictors),
  fixef(PyreneesModel, pars = predictors),
  fixef(SerengetiModel, pars = predictors)) %>%
  as.data.frame() %>%
  cbind(model = rep(foodwebs, each = length(predictors))) %>%
  cbind(coef = factor(rep(predictors, times = length(foodwebs)), levels = predictors))
  
ggplot(fixed_effects) +
  geom_pointrange(aes(x = coef, y = Estimate, ymin = Q2.5, ymax = Q97.5, color = model), position = position_dodge(width=0.9))

# random effects
ArcticRandomEffects <- ranef(ArcticModel, pars = predictors)
EuropeRandomEffects <- ranef(EuropeModel, pars = predictors)
PyreneesRandomEffects <- ranef(PyreneesModel, pars = predictors)
SerengetiRandomEffects <- ranef(SerengetiModel, pars = predictors)

RandomEffects <- data.frame(model = c(), order = c(), coef = factor(c(), levels = predictors), est = c(), Q2.5 = c(), Q97.5 = c())
for (ipredictor in c(1:length(predictors))){
  RandomEffecti <-
    bind_rows(
      ArcticRandomEffects$Order.predator[,,ipredictor] %>%
        as.data.frame() %>%
        rownames_to_column(var = "order") %>%
        mutate(model = "Arctic", coef = predictors[ipredictor]),
      EuropeRandomEffects$Order.predator[,,ipredictor] %>%
        as.data.frame() %>%
        rownames_to_column(var = "order") %>%
        mutate(model = "Europe", coef = predictors[ipredictor]),
      PyreneesRandomEffects$Order.predator[,,ipredictor] %>%
        as.data.frame() %>%
        rownames_to_column(var = "order") %>%
        mutate(model = "Pyrenees", coef = predictors[ipredictor]),
      SerengetiRandomEffects$Order.predator[,,ipredictor] %>%
        as.data.frame() %>%
        rownames_to_column(var = "order") %>%
        mutate(model = "Serengeti", coef = predictors[ipredictor]),
      ) %>% 
    mutate(coef = factor(coef, levels = predictors))
  # add fixed effect
  fe <- filter(fixed_effects, coef == predictors[ipredictor])
  RandomEffecti[,c(2,4,5)] <- RandomEffecti[,c(2,4,5)] + fe$Estimate[match(RandomEffecti$model, fe$model)]
  RandomEffects <- rbind(RandomEffects, RandomEffecti)
}

# Plot coefficients -------------------------------------------------------
coef_plot <- ggplot(aes(x = 1, y = exp(Estimate)), data = fixed_effects) +
  geom_point(data = fixed_effects, aes(fill = model), position = position_dodge(width=0.9), shape = 21, size = 2) +
  geom_linerange(data = fixed_effects, aes(ymin = exp(Q2.5), ymax = exp(Q97.5), color = model), position = position_dodge(width = 0.9)) +
  geom_hline(yintercept = 1)+
  scale_y_continuous(limits = c(0.002,30), trans = "log10")+
  scale_fill_manual(values = c("deepskyblue","royalblue4", "red3", "chartreuse4")) +
  scale_color_manual(values = c("deepskyblue","royalblue4", "red3", "chartreuse4")) +
  facet_wrap(vars(coef), scales = "free_x", nrow = 2)+
  labs(color = "Model", fill = "Model", y = "Odds ratio") + 
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())

ggsave("figures/SI/coef_plot.png", coef_plot, width = 7)
