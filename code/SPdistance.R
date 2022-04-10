# Step xx: Measure distances between species in target food web to species in Europe
# 1. Calculate nearest taxon phylogenetic distance
# 2. Calculate mean functional distance
# 3. Make plots

library(ape)
library(dplyr)
library(tidyr)
library(purrr)
library(mFD)
library(patchwork)
load("code/DistanceFunctions.R")

# Phylogenetic distance to nearest taxon in Europe ------------------------
# species list
europe <- read.csv("data/cleaned/EuroMWTaxo.csv", row.names = 1)
pyrenees <- read.csv("data/cleaned/pyrenneesFWTaxo.csv", row.names = 1)
serengeti <- read.csv("data/cleaned/SerengetiFWTaxo.csv", row.names = 1) %>% drop_na()
arctic <- read.csv("data/cleaned/HighArcticFWTaxo.csv", row.names = 1)

# load matrix of phylogenetic distance
phydist <- as.matrix(read.csv("data/checkpoints/phylodist.csv", row.names = 1))

# calculate distance
arctic_mntd <- map_dbl(arctic$Species, mntd, europe$Species, phydist)
pyrenees_mntd <- map_dbl(pyrenees$Species, mntd, europe$Species, phydist)
serengeti_mntd <- map_dbl(serengeti$Species, mntd, europe$Species, phydist)

# Mean functional distance to species in Europe ---------------------------
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

# calculate each species mean functional distance
arctic_fmpd <- map_dbl(arctic$Species, fmpd, europe$Species, sp_funct.dist)
pyrenees_fmpd <- map_dbl(pyrenees$Species, fmpd, europe$Species, sp_funct.dist)
serengeti_fmpd <- map_dbl(serengeti$Species, fmpd, europe$Species, sp_funct.dist)


# Make plot ---------------------------------------------------------------
fills <- c("Serengeti" = "yellowgreen", "Pyrenees" = "red3", "Arctic" = "deepskyblue")

# functional distance
p1 <- ggplot() +
  geom_histogram(data = data.frame(serengeti_fmpd), aes(serengeti_fmpd, fill = "Serengeti"), alpha= 0.3) +
  geom_histogram(data = data.frame(pyrenees_fmpd), aes(pyrenees_fmpd, fill = "Pyrenees"), alpha= 0.3) +
  geom_histogram(data = data.frame(arctic_fmpd), aes(arctic_fmpd, fill = "Arctic"), alpha= 0.3) +
  theme_minimal() +
  labs(fill = "Food web", x = "Mean functional distance") +
  scale_fill_manual(values = fills) +
  theme(legend.position = "bottom",
        legend.direction = "horizontal")

# phylogenetic distance
p2 <- ggplot() +
  geom_histogram(data = data.frame(serengeti_mntd), aes(serengeti_mntd, fill = "Serengeti"), alpha= 0.3) +
  geom_histogram(data = data.frame(pyrenees_mntd), aes(pyrenees_mntd, fill = "Pyrenees"), alpha= 0.3) +
  geom_histogram(data = data.frame(arctic_mntd), aes(arctic_mntd, fill = "Arctic"), alpha= 0.3) +
  theme_minimal() +
  labs(fill = "Food web", x = "Phylogenetic distance to neareast taxon") +
  scale_fill_manual(values = fills) +
  theme(legend.position = "bottom",
        legend.direction = "horizontal")

# patchwork and save
p<-p1+p2
ggsave("figures/SI/SPdist.png", p, device = "png")
