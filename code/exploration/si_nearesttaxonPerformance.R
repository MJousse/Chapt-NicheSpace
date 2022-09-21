# SI - Phylogenetic distance categorical
library(tidyr)
library(dplyr)
library(purrr)
library(brms)
library(ggplot2)
library(tidybayes)
source("code/functions_distance.R")

# set taxonomic distance as 0 = same species, 1 = same genus, 2 = same family, 3 = same order, 4 = same class, 5 = different class
taxonomy <- read.csv("data/cleaned/SpeciesTraitsFull.csv", row.names = 1) %>%
  select(Species, Genus, Family, Order, Class)

taxonomy <- expand_grid(taxonomy, taxonomy, .name_repair = "minimal")
colnames(taxonomy) <- c("Species1", "Genus1", "Family1", "Order1", "Class1", "Species2", "Genus2", "Family2", "Order2", "Class2")

taxonomy$dist <- 5
taxonomy$dist[taxonomy$Class1 == taxonomy$Class2] <- 4
taxonomy$dist[taxonomy$Order1 == taxonomy$Order2] <- 3
taxonomy$dist[taxonomy$Family1 == taxonomy$Family2] <- 2
taxonomy$dist[taxonomy$Genus1 == taxonomy$Genus2] <- 1
taxonomy$dist[taxonomy$Species1 == taxonomy$Species2] <- 0

taxonomyDist <- taxonomy %>%
  select(Species1, Species2, dist) %>%
  pivot_wider(names_from = Species2, values_from = dist)
taxonomyDist <- select(taxonomyDist, -Species1) %>% as.matrix()
rownames(taxonomyDist) <- colnames(taxonomyDist)

# Add distances to species performance ------------------------------------
species_performance <- read.csv("data/checkpoints/species_performance.csv", row.names = 1) %>%
  filter(Source != Target)

# species list
europe <- read.csv("data/cleaned/EuroMWTaxo.csv", row.names = 1)
pyrenees <- read.csv("data/cleaned/pyrenneesFWTaxo.csv", row.names = 1)
serengeti <- read.csv("data/cleaned/SerengetiFWTaxo.csv", row.names = 1) %>% drop_na()
arctic <- read.csv("data/cleaned/HighArcticFWTaxo.csv", row.names = 1)

# mean nearest taxon distance
species_performance$mntd <- NA
species_performance$mntd[species_performance$Source == "Euro"] <- map_dbl(species_performance$species[species_performance$Source == "Euro"], mntd, europe$Species, taxonomyDist)
species_performance$mntd[species_performance$Source == "Pyrenees"] <- map_dbl(species_performance$species[species_performance$Source == "Pyrenees"], mntd, pyrenees$Species, taxonomyDist)
species_performance$mntd[species_performance$Source == "Serengeti"] <- map_dbl(species_performance$species[species_performance$Source == "Serengeti"], mntd, serengeti$Species, taxonomyDist)
species_performance$mntd[species_performance$Source == "Arctic"] <- map_dbl(species_performance$species[species_performance$Source == "Arctic"], mntd, arctic$Species, taxonomyDist)

# transform responses
species_performance$logitauc <- log(species_performance$auc / (1.001-species_performance$auc))
species_performance$logaucpr <- log(species_performance$aucpr / species_performance$prevalence)

species_performance$mntd <- factor(species_performance$mntd, ordered = T, labels = c("Shared species", "Same genus", "Same family", "Same order", "Same class", "Different class"))

auc_taxo <- brm(logitauc ~ mo(mntd) + (1|Source) + (1|Target),
                data = species_performance,
                prior = c(
                  prior(normal(0, 1), class = "Intercept"),
                  prior(normal(0, 1), class = "b"),
                  prior(cauchy(0, 5), class = "sd")
                ), 
                sample_prior = "no",
                iter = 2000)

taxo_fe <- tibble(mntd = factor(c(0:5))) %>%
  add_epred_draws(auc_taxo,
                  re_formula = NA, ndraws = 1e3) %>%
  mutate(x = mntd,
         y = exp(.epred)/(1+exp(.epred)))

taxo_fe <- taxo_fe %>% 
  group_by(x)  %>%
  summarize(median_qi(y, width = 0.95)) %>%
  select(x, y, ymin, ymax) %>%
  mutate(x = factor(x = x, labels = c("Shared species", "Same genus", "Same family", "Same order", "Same class", "Different class")))

ggplot(taxo_fe,
       aes(x = x, y = y)) +
  geom_jitter(aes(x = mntd, y = auc), data = species_performance, size = 0.2, alpha = 0.3, width = 0.2) +
  geom_pointrange(aes(ymin = ymin, ymax = ymax), shape = 21, fill = "white") +
  labs(y = "AUC", x = "Nearest taxon")+
  theme_minimal() +
  lims(y = c(0.25,1))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

ggsave("figures/SI/SpeciesPerformanceTaxo.png")
