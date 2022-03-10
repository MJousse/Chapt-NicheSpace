rm(list = ls())
library(lme4)
library(MuMIn)

EuroTraits <- read.csv("data/cleaned/SpeciesTraitsFull.csv",, row.names = 1, stringsAsFactors = T)

# Habitat_breadth
habitat_breadth<- lmer(Habitat_breadth_IUCN ~ (1 | Class/Order/Family/Genus), data = EuroTraits, REML = F)

AIC.table  <- model.sel(habitat_breadth_0, habitat_breadth_1, habitat_breadth_2, habitat_breadth_3, habitat_breadth_4)

# Diel_activity
diel_activity <- glmer(Diel_activity ~ (1 | Class/Order/Family/Genus), data = EuroTraits, family = "binomial")
