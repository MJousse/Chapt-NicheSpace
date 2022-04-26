library(tidybayes)
library(bayesplot)

# Load model
ArcticModel <- readRDS("~/OneDrive/Chapt-NicheSpace/models/ArcticModel_brms.rds")
EuropeModel <- readRDS("~/OneDrive/Chapt-NicheSpace/models/EuroModel_brms.rds")
PyreneesModel <- readRDS("~/OneDrive/Chapt-NicheSpace/models/PyreneesModel_brms.rds")
SerengetiModel <- readRDS("~/OneDrive/Chapt-NicheSpace/models/SerengetiModel_brms.rds")

# Check convergence with Rhat and trace plots
summary(ArcticModel)
mcmc_trace(ArcticModel, regex_pars = "b_")
summary(EuropeModel)
mcmc_trace(EuropeModel, regex_pars = "b_")
summary(PyreneesModel)
mcmc_trace(PyreneesModel, regex_pars = "b_")
summary(SerengetiModel)
mcmc_trace(SerengetiModel, regex_pars = "b_")

# Posterior predictive checks
pp_check(ArcticModel, ndraws = 100, type = "bars_grouped", group = "Order.predator")
pp_check(EuropeModel, ndraws = 100, type = "bars_grouped", group = "Order.predator")
pp_check(PyreneesModel, ndraws = 100, type = "bars_grouped", group = "Order.predator")
pp_check(SerengetiModel, ndraws = 100, type = "bars_grouped", group = "Order.predator")

# Bayes R2
R2_arctic <- bayes_R2(ArcticModel)
R2_europe <- bayes_R2(EuropeModel)
R2_pyrenees <- bayes_R2(PyreneesModel)
R2_serengeti <- bayes_R2(SerengetiModel)

# Compare Coefficients
fixed_effects_names <- c("Intercept", "Omnivore.predator", "Carnivore.predator", "Habitat_breadth.predator",
       "BM.predator", "Longevity.predator", "ClutchSize.predator", "Omnivore.prey",
       "Habitat_breadth.prey", "BM.prey", "Longevity.prey", "ClutchSize.prey", "ActivityTime.match",
       "Habitat.match", "BM.match")
foodwebs <- c("Arctic", "Europe", "Pyrenees", "Serengeti")
fixed_effects <- rbind(
  fixef(ArcticModel, pars = fixed_effects_names),
  fixef(EuropeModel, pars = fixed_effects_names),
  fixef(PyreneesModel, pars = fixed_effects_names),
  fixef(SerengetiModel, pars = fixed_effects_names)) %>%
  as.data.frame() %>%
  cbind(model = rep(foodwebs, each = length(fixed_effects_names))) %>%
  cbind(coef = factor(rep(fixed_effects_names, times = length(foodwebs)), levels = fixed_effects_names))
  
ggplot(fixed_effects) +
  geom_pointrange(aes(x = coef, y = Estimate, ymin = Q2.5, ymax = Q97.5, color = FW), position = position_dodge(width=0.9))

