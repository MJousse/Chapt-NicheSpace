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
library(tidybayes)
library(mFD)
library(patchwork)
library(brms)
library(ggplot2)
source("code/functions_distance.R")

# species list
europe <- read.csv("data/cleaned/EuroMWTaxo.csv", row.names = 1)
pyrenees <- read.csv("data/cleaned/pyrenneesFWTaxo.csv", row.names = 1)
serengeti <- read.csv("data/cleaned/SerengetiFWTaxo.csv", row.names = 1) %>% drop_na()
arctic <- read.csv("data/cleaned/HighArcticFWTaxo.csv", row.names = 1)

# phylogenetic distance matrix
phydist <- as.matrix(read.csv("data/checkpoints/phylodist_mean.csv", row.names = 1))
phydist <- sqrt(phydist)
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
  filter(Source != Target)

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

# remove species that were not in phylogeny
species_performance <- filter(species_performance, species %in% colnames(phydist))

write.csv(species_performance, "data/checkpoints/species_performance_trans.csv")

# Correlate transferability to distance metrics ---------------------------
# glm with species-specific logit-auc and log(aucpr/prevalence) as response
# phylogenetic and functional distance as predictors
# source and target region as crossed random effect + prevalence as random effect

# transform responses
species_performance$logitauc <- log(species_performance$auc / (1.001-species_performance$auc))
species_performance$logaucpr <- log(species_performance$aucpr / species_performance$prevalence)

# scale predictors
species_performance$mntd_sc <- as.vector(scale(species_performance$mntd))
species_performance$fmpd_sc <- as.vector(scale(species_performance$fmpd))
species_performance$prevalence_sc <- as.vector(scale(species_performance$prevalence))

# model with logit auc as response using brms
auc_mntd <- brm(logitauc ~ s(mntd_sc, k = 5) + (1|Source) + (1|Target),
                 data = species_performance,
                 prior = c(
                   prior(normal(0, 1), class = "Intercept"),
                   prior(normal(0, 1), class = "b"),
                   prior(cauchy(0, 5), class = "sd")
                 ), 
                 sample_prior = "no",
                 iter = 2000)

auc_fmpd <- brm(logitauc ~ s(fmpd_sc, k = 5) + (1|Source) + (1|Target),
                data = species_performance,
                prior = c(
                  prior(normal(0, 1), class = "Intercept"),
                  prior(normal(0, 1), class = "b"),
                  prior(cauchy(0, 5), class = "sd")
                ), 
                sample_prior = "no",
                iter = 2000)

auc_prevalence <- brm(logitauc ~ s(prevalence_sc, k = 5) + (1|Source) + (1|Target),
                data = species_performance,
                prior = c(
                  prior(normal(0, 1), class = "Intercept"),
                  prior(normal(0, 1), class = "b"),
                  prior(cauchy(0, 5), class = "sd")
                ), 
                sample_prior = "no",
                iter = 2000)

auc_all <- brm(logitauc ~ mntd_sc + fmpd_sc + prevalence_sc + (1|Source) + (1|Target),
                      data = species_performance,
                      prior = c(
                        prior(normal(0, 1), class = "Intercept"),
                        prior(normal(0, 1), class = "b"),
                        prior(cauchy(0, 5), class = "sd")
                      ), 
                      sample_prior = "no",
                      iter = 2000)

# model with log aucpr/prevalence as response using brms
aucpr_model <- brm(logaucpr ~ mntd_sc + fmpd_sc + prevalence_sc + (1|Source) + (1|Target),
                   data = species_performance,
                   prior = c(
                     prior(normal(0, 1), class = "Intercept"),
                     prior(normal(0, 1), class = "b"),
                     prior(cauchy(0, 5), class = "sd")
                   ), 
                   sample_prior = "no",
                   iter = 2000)

# plot performance to phylogenetic distance relation
phylo_fe <- tibble(mntd_sc = seq(min(species_performance$mntd_sc), 
                                   max(species_performance$mntd_sc), length.out=100), Source = NA, Target = NA) %>%
  add_epred_draws(auc_mntd,
                  re_formula = NA, ndraws = 1e3) %>%
  mutate(x = (mntd_sc * sd(species_performance$mntd) + mean(species_performance$mntd)),
         y = exp(.epred)/(1+exp(.epred)))

phylo_fe <- phylo_fe %>% 
  group_by(x)  %>%
  summarize(median_qi(y, width = 0.95)) %>%
  select(x, y, ymin, ymax)

p1 <- ggplot(phylo_fe,
             aes(x = x, y = y)) +
  geom_point(aes(x = mntd, y = auc), data = species_performance, size = 0.2, alpha = 0.3) +
  geom_lineribbon(aes(ymin = ymin, ymax = ymax), alpha= 0.5) +
  labs(y = "AUC", x = "Distance to nearest taxon")+
  theme_minimal() +
  lims(y = c(0.25,1))+
  theme(axis.line = element_line(size = 0.5), strip.text = element_text(size = 12),
        strip.background = element_rect(colour = "black"), panel.grid = element_blank())

# plot performance to functional distance relation
func_fe <- tibble(fmpd_sc = seq(min(species_performance$fmpd_sc), 
                                 max(species_performance$fmpd_sc), length.out=100), Source = NA, Target = NA) %>%
  add_epred_draws(auc_fmpd,
                  re_formula = NA, ndraws = 1e3) %>%
  mutate(x = (fmpd_sc * sd(species_performance$fmpd) + mean(species_performance$fmpd)),
         y = exp(.epred)/(1+exp(.epred)))

func_fe <- func_fe %>% 
  group_by(x)  %>%
  summarize(median_qi(y, width = 0.95)) %>%
  select(x, y, ymin, ymax)

p2 <- ggplot(func_fe,
             aes(x = x, y = y)) +
  geom_point(aes(x = fmpd, y = auc), data = species_performance, size = 0.2, alpha = 0.3) +
  geom_lineribbon(aes(ymin = ymin, ymax = ymax), alpha= 0.5) +
  labs(y = "AUC", x = "Mean functional distance")+
  theme_minimal() +
  lims(y = c(0.25,1)) +
  theme(axis.line = element_line(size = 0.5), panel.grid = element_blank(),
        axis.title.y = element_blank(), axis.text.y = element_blank())

# plot performance to prevalence
prev_fe <- tibble(prevalence_sc = seq(min(species_performance$prevalence_sc), 
                                max(species_performance$prevalence_sc), length.out=100), Source = NA, Target = NA) %>%
  add_epred_draws(auc_prevalence,
                  re_formula = NA, ndraws = 1e3) %>%
  mutate(x = prevalence_sc * sd(species_performance$prevalence) + mean(species_performance$prevalence),
         y = exp(.epred)/(1+exp(.epred)))

prev_fe <- prev_fe %>% 
  group_by(x)  %>%
  summarize(median_qi(y, width = 0.95)) %>%
  select(x, y, ymin, ymax)

p3 <- ggplot(prev_fe,
             aes(x = x, y = y)) +
  geom_point(aes(x = prevalence, y = auc), data = species_performance, size = 0.2, alpha = 0.3) +
  geom_lineribbon(aes(ymin = ymin, ymax = ymax), alpha= 0.5) +
  labs(y = "AUC", x = "Normalized degree")+
  theme_minimal() +
  lims(y = c(0.25,1)) +
  theme(axis.line = element_line(size = 0.5), panel.grid = element_blank(),
        axis.title.y = element_blank(), axis.text.y = element_blank())

p1 + p2 + p3

ggsave("figures/SpeciesPerformance.png", width = 9)

# mydatab <- data.frame(
#   mntd_sc = species_performance$mntd_sc,
#   Source = species_performance$Source,
#   Target = species_performance$Target,
#   prevalence_sc = 0,
#   fmpd_sc = 0
# )
# 
# species_performance_fitmntd <- cbind(species_performance, fitted(auc_model, mydatab, re_formula=NULL))
# species_performance_fitmntd$Estimate <- exp(species_performance_fitmntd$Estimate)/(1+exp(species_performance_fitmntd$Estimate))
# species_performance_fitmntd$Q2.5 <- exp(species_performance_fitmntd$Q2.5)/(1+exp(species_performance_fitmntd$Q2.5))
# species_performance_fitmntd$Q97.5 <- exp(species_performance_fitmntd$Q97.5)/(1+exp(species_performance_fitmntd$Q97.5))
# 
# p1 <- ggplot(species_performance_fitmntd) +
#   geom_point(aes(y = auc, x = mntd, colour = Source), size = 0.5, alpha = .3) +
#   geom_ribbon(aes(x = mntd, y = Estimate, ymin = Q2.5, ymax = Q97.5, fill = Source),
#               alpha = .3) +
#   geom_line(aes(x = mntd, y = Estimate, color = Source),
#             size = 1) +
#   scale_fill_manual(values = c("deepskyblue","royalblue4", "red3", "chartreuse4")) +
#   scale_color_manual(values = c("deepskyblue","royalblue4", "red3", "chartreuse4")) +
#   labs(x = "Phylogenetic distance", y = "", colour = "Model", fill = "Model") + 
#   facet_wrap(~Target, nrow = 1) +
#   theme_minimal() +
#   theme(axis.line = element_line(size = 0.5), strip.text = element_text(size = 10),
#         strip.background = element_rect(colour = "black"), axis.text = element_text(size = 8), axis.title = element_text(size = 10), legend.text = element_text(size = 8), legend.key.size = unit(0.5,"cm"))
# 
# mydatab <- data.frame(
#   mntd_sc = 0,
#   Source = species_performance$Source,
#   Target = species_performance$Target,
#   prevalence_sc = 0,
#   fmpd_sc = species_performance$fmpd_sc
# )
# 
# species_performance_fitfmpd <- cbind(species_performance, fitted(auc_model, mydatab, re_formula=NULL))
# species_performance_fitfmpd$Estimate <- exp(species_performance_fitfmpd$Estimate)/(1+exp(species_performance_fitfmpd$Estimate))
# species_performance_fitfmpd$Q2.5 <- exp(species_performance_fitfmpd$Q2.5)/(1+exp(species_performance_fitfmpd$Q2.5))
# species_performance_fitfmpd$Q97.5 <- exp(species_performance_fitfmpd$Q97.5)/(1+exp(species_performance_fitfmpd$Q97.5))
# 
# p2 <- ggplot(species_performance_fitfmpd) +
#   geom_point(aes(y = auc, x = fmpd, colour = Source), size = 0.5, alpha = .3) +
#   geom_ribbon(aes(x = fmpd, y = Estimate, ymin = Q2.5, ymax = Q97.5, fill = Source),
#               alpha = .3) +
#   geom_line(aes(x = fmpd, y = Estimate, color = Source),
#             size = 1) +
#   scale_fill_manual(values = c("deepskyblue","royalblue4", "red3", "chartreuse4")) +
#   scale_color_manual(values = c("deepskyblue","royalblue4", "red3", "chartreuse4")) +
#   labs(x = "Functional distance", y = "roc-auc", colour = "Model", fill = "Model") + 
#   facet_wrap(~Target, nrow = 1) +
#   theme_minimal() +
#   theme(axis.line = element_line(size = 0.5), strip.text = element_blank(), axis.text = element_text(size = 8), axis.title = element_text(size = 10), legend.text = element_text(size = 8), legend.key.size = unit(0.5,"cm"))
# 
# mydatab <- data.frame(
#   mntd_sc = 0,
#   Source = species_performance$Source,
#   Target = species_performance$Target,
#   prevalence_sc = species_performance$prevalence_sc,
#   fmpd_sc = 0
# )
# 
# species_performance_fitprev <- cbind(species_performance, fitted(auc_model, mydatab, re_formula=NULL))
# species_performance_fitprev$Estimate <- exp(species_performance_fitprev$Estimate)/(1+exp(species_performance_fitprev$Estimate))
# species_performance_fitprev$Q2.5 <- exp(species_performance_fitprev$Q2.5)/(1+exp(species_performance_fitprev$Q2.5))
# species_performance_fitprev$Q97.5 <- exp(species_performance_fitprev$Q97.5)/(1+exp(species_performance_fitprev$Q97.5))
# 
# p3 <- ggplot(species_performance_fitprev) +
#   geom_point(aes(y = auc, x = prevalence, colour = Source), size = 0.5, alpha = .3) +
#   geom_ribbon(aes(x = prevalence, y = Estimate, ymin = Q2.5, ymax = Q97.5, fill = Source),
#               alpha = .3) +
#   geom_line(aes(x = prevalence, y = Estimate, color = Source),
#             size = 1) +
#   scale_fill_manual(values = c("deepskyblue","royalblue4", "red3", "chartreuse4")) +
#   scale_color_manual(values = c("deepskyblue","royalblue4", "red3", "chartreuse4")) +
#   scale_x_continuous(breaks = c(0, 0.1, 0.2, 0.3))+
#   labs(x = "Prevalence", y = "", colour = "Model", fill = "Model") + 
#   facet_wrap(~Target, nrow = 1) +
#   theme_minimal() +
#   theme(axis.line = element_line(size = 0.5), strip.text = element_blank(), axis.text = element_text(size = 8), axis.title = element_text(size = 10), legend.text = element_text(size = 8), legend.key.size = unit(0.5,"cm"))
# 
# p<-p1/p2/p3 + plot_layout(guides = "collect") +
#   plot_annotation(title = 'Predicted food web', theme = theme(plot.title = element_text(hjust = 0.45, size = 10)))
# 
# ggsave("figures/SpeciesPerformance.png", p, width = 18, height = 12, units = "cm", dpi = 600)
