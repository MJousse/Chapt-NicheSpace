# Step 09: Correlate species-specific performance to how distant it is from the species pool
# of the food webs on which the model has been trained on
# 1. Calculate the functional distance between all species
# 2. Calculate mean functional distance and phylogenetic distance to nearest taxon
# 3. Correlated species-specific performance to these distances.

rm(list =ls())
library(ape)
library(dplyr)
library(tidyr)
library(purrr)
library(mFD)
library(patchwork)
source("code/DistanceFunctions.R")

# species list
europe <- read.csv("data/cleaned/EuroMWTaxo.csv", row.names = 1)
pyrenees <- read.csv("data/cleaned/pyrenneesFWTaxo.csv", row.names = 1)
serengeti <- read.csv("data/cleaned/SerengetiFWTaxo.csv", row.names = 1) %>% drop_na()
arctic <- read.csv("data/cleaned/HighArcticFWTaxo.csv", row.names = 1)

# phylogenetic distance matrix
phydist <- as.matrix(read.csv("data/checkpoints/phylodist_mean.csv", row.names = 1))
colnames(phydist) <- gsub("\\.", " ", colnames(phydist))

# Functional distance matrix ----------------------------------------------
traits <- read.csv("data/cleaned/SpeciesTraitsFull.csv", row.names = 1, stringsAsFactors = T)

# create objects for mFD
rownames(traits) <- traits$Species
traits <- traits %>%
  dplyr::select(-Species, -Class, -Order, -Family, -Genus) %>%
  drop_na()
traits_cat <- data.frame(trait_name = colnames(traits),
                         trait_type = c("N", "Q", rep("F", 12), "Q", "Q", "Q", "F", "F", "F"),
                         fuzzy_name = c(NA, NA, rep("Habitat", 12), NA, NA, NA, rep("TrophicLevel", 3)))

# calculate pairwise distance
sp_funct.dist <- funct.dist(traits, traits_cat, metric = "gower", scale_euclid = "scale_center")

# Add distances to species performance ------------------------------------
species_performance <- read.csv("data/checkpoints/species_performance.csv", row.names = 1) %>%
  filter(Source != Target, auc != 1)

# mean nearest taxon distance
species_performance$mntd <- NA
species_performance$mntd[species_performance$Source == "Euro"] <- map_dbl(species_performance$species[species_performance$Source == "Euro"], mntd, europe$Species, phydist)
species_performance$mntd[species_performance$Source == "Pyrenees"] <- map_dbl(species_performance$species[species_performance$Source == "Pyrenees"], mntd, pyrenees$Species, phydist)
species_performance$mntd[species_performance$Source == "Serengeti"] <- map_dbl(species_performance$species[species_performance$Source == "Serengeti"], mntd, serengeti$Species, phydist)
species_performance$mntd[species_performance$Source == "Arctic"] <- map_dbl(species_performance$species[species_performance$Source == "Arctic"], mntd, arctic$Species, phydist)

# functional mean pairwise distance
species_performance$fmpd <- NA
species_performance$fmpd[species_performance$Source == "Euro"] <- map_dbl(species_performance$species[species_performance$Source == "Euro"], fmpd, europe$Species, sp_funct.dist)
species_performance$fmpd[species_performance$Source == "Pyrenees"] <- map_dbl(species_performance$species[species_performance$Source == "Pyrenees"], fmpd, pyrenees$Species, sp_funct.dist)
species_performance$fmpd[species_performance$Source == "Serengeti"] <- map_dbl(species_performance$species[species_performance$Source == "Serengeti"], fmpd, serengeti$Species, sp_funct.dist)
species_performance$fmpd[species_performance$Source == "Arctic"] <- map_dbl(species_performance$species[species_performance$Source == "Arctic"], fmpd, arctic$Species, sp_funct.dist)

# Correlate transferability to distance metrics ---------------------------
# glm with species-specific logit-auc and log(aucpr/prevalence) as response
# phylogenetic and functional distance as predictors
# source and target region as crossed random effect + prevalence as random effect

# transform responses
species_performance$logitauc <- log(species_performance$auc / (1-species_performance$auc))
species_performance$logaucpr <- log(species_performance$aucpr / species_performance$prevalence)

# scale predictors
species_performance$mntd <- as.vector(scale(species_performance$mntd))
species_performance$fmpd <- as.vector(scale(species_performance$fmpd))
species_performance$prevalence <- as.vector(scale(log(species_performance$prevalence)))

# model with logit auc as response using brms
auc_model <- brm(logitauc ~ mntd + fmpd + prevalence + (1|Source) + (1|Target),
                 data = species_performance,
                 prior = c(
                   prior(normal(0, 1), class = "Intercept"),
                   prior(normal(0, 1), class = "b"),
                   prior(cauchy(0, 5), class = "sd")
                 ), 
                 sample_prior = "no",
                 iter = 2000)

# model with log aucpr/prevalence as response using brms
aucpr_model <- brm(logaucpr ~ mntd + fmpd + prevalence + (1|Source) + (1|Target),
                   data = species_performance,
                   prior = c(
                     prior(normal(0, 1), class = "Intercept"),
                     prior(normal(0, 1), class = "b"),
                     prior(cauchy(0, 5), class = "sd")
                   ), 
                   sample_prior = "no",
                   iter = 2000)
