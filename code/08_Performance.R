# Step 08: Calculate the performance of models for each target food web
# For combination of mode-target food web
# 1. Calculate the overall roc-auc and pr-auc
# 2. Calculate the species specific roc-auc and pr-auc
# 3. Create plot of transferability
# 4. Correlate transferability to dissimilarity metrics

# set up
rm(list = ls())
library(ROCR)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(brms)
library(tidybayes)
load("data/checkpoints/predictions.RData")
foodwebs <- c("Arctic", "Euro", "Pyrenees", "Serengeti")
source("code/functions.R")

# Calculate performance ---------------------------------------------------
overall_performance <- expand_grid(Source = foodwebs, Target = foodwebs) %>%
  mutate(auc = NA, aucpr = NA, prevalence = NA)
overall_performance_draws <- data.frame()
species_performance <- data.frame()
species_performance_draws <- data.frame()

# for each combinaison calculate the overall auc and aucpr
# in each target food web, calculate the auc and aucpr of each species
for (combination in c(1:nrow(overall_performance))){
  sourceFW <- overall_performance[combination, "Source"]
  targetFW <- overall_performance[combination, "Target"]
  predictions <- get(paste0(sourceFW, "_", targetFW, "_predictions"))
  if (sourceFW == targetFW){
    predictions_test <- filter(predictions, training == 0)
  } else {
    predictions_test <- predictions
  }
  overall_performance[combination, "auc"] <- performance(prediction(predictions_test$Estimate, predictions_test$interaction), "auc")@y.values[[1]]
  overall_performance[combination, "aucpr"] <- performance(prediction(predictions_test$Estimate, predictions_test$interaction), "aucpr")@y.values[[1]]
  overall_performance[combination, "prevalence"] <- sum(predictions_test$interaction)/nrow(predictions_test)
  overall_performance_draws <- rbind(overall_performance_draws,
                                     data.frame(sourceFW = sourceFW, targetFW= targetFW,
                                           performance_draws(select(predictions_test, starts_with("draws")), predictions_test$interaction)))
  fw_sp <- unique(predictions$Predator)
  for (species in fw_sp){
   sp_predictions <- filter(predictions, Predator == species | Prey == species)
   if(sum(sp_predictions$interaction != 0)){
     sp_auc <- performance(prediction(sp_predictions$Estimate, sp_predictions$interaction), "auc")@y.values[[1]]
     sp_aucpr <- performance(prediction(sp_predictions$Estimate, sp_predictions$interaction), "aucpr")@y.values[[1]]
     species_performance <- rbind(species_performance,
                                  data.frame(sourceFW, targetFW, species, auc = sp_auc, aucpr = sp_aucpr, prevalence = sum(sp_predictions$interaction)/nrow(sp_predictions)))
     species_performance_draws <- rbind(species_performance_draws,
                                        data.frame(sourceFW = sourceFW, targetFW= targetFW, species,
                                                   performance_draws(select(sp_predictions, starts_with("draws")), sp_predictions$interaction)))
   }
  }
}
write.csv(species_performance, "data/checkpoints/species_performance.csv")
write.csv(overall_performance, "data/checkpoints/overall_performance.csv")
write.csv(species_performance_draws, "data/checkpoints/species_performance_draws.csv")
write.csv(overall_performance_draws, "data/checkpoints/overall_performance_draws.csv")

# Visualize transferability -----------------------------------------------
# add dissimilarity
species_performance <- read.csv("data/checkpoints/species_performance.csv", row.names = 1)
overall_performance <- read.csv("data/checkpoints/overall_performance.csv", row.names = 1)
FWdist <- read.csv("data/checkpoints/FWdist.csv", row.names = 1)
FWdist[!is.na(FWdist) & FWdist == "Nunavik"] <- "Arctic"
overall_performance[overall_performance == "Euro"] <- "Europe"
species_performance[species_performance == "Euro"] <- "Europe"
overall_performance <- left_join(overall_performance, FWdist, by = c("Source" = "FW2", "Target" = "FW1"))
species_performance <- left_join(species_performance, FWdist, by = c("Source" = "FW2", "Target" = "FW1"))

# geographic distance
p1 <- ggplot(species_performance) +
  geom_jitter(aes(x = geo.dist/1000, y = auc, colour = Source), alpha = 0.2, height = 0, width = 0.3) +
  geom_point(data = overall_performance, aes(x = geo.dist/1000, y = auc, fill = Source), shape = 21, size = 3) +
  scale_fill_manual(values = c("deepskyblue","royalblue4", "red3", "chartreuse4")) +
  scale_color_manual(values = c("deepskyblue","royalblue4", "red3", "chartreuse4")) +
  labs(y = "", x = "Geographic distance (10³km)", fill = "Source food web", colour = "Source food web")+
  facet_wrap(~Target, nrow = 1) +
  theme_minimal() +
  theme(axis.line = element_line(size = 0.5), strip.text = element_text(size = 12),
        strip.background = element_rect(colour = "black"))

# environmental distance
p2 <- ggplot(species_performance) +
  geom_jitter(aes(x = env.dist, y = auc, colour = Source), alpha = 0.2, height = 0, width = 0.2) +
  geom_point(data = overall_performance, aes(x = env.dist, y = auc, fill = Source), shape = 21, size = 3) +
  scale_fill_manual(values = c("deepskyblue","royalblue4", "red3", "chartreuse4")) +
  scale_color_manual(values = c("deepskyblue","royalblue4", "red3", "chartreuse4")) +
  labs(y = "roc-auc", x = "Environmental distance", fill = "Source food web", colour = "Source food web") +
  facet_wrap(~Target, nrow = 1) +
  theme_minimal() +
  theme(axis.line = element_line(size = 0.5), strip.text = element_blank())

# phylogenetic distance
p3 <- ggplot(species_performance) +
  geom_jitter(aes(x = phylo.dist, y = auc, colour = Source), alpha = 0.2, height = 0, width = 0.3) +
  geom_point(data = overall_performance, aes(x = phylo.dist, y = auc, fill = Source), shape = 21, size = 3) +
  scale_fill_manual(values = c("deepskyblue","royalblue4", "red3", "chartreuse4")) +
  scale_color_manual(values = c("deepskyblue","royalblue4", "red3", "chartreuse4")) +
  labs(y = "", x = "Phylogenetic distance", fill = "Source food web", colour = "Source food web") +
  facet_wrap(~Target, nrow = 1) +
  theme_minimal() +
  theme(axis.line = element_line(size = 0.5), strip.text = element_blank())

# everything together and save
p <- (p1 / p2 / p3)+ plot_layout(guides = "collect") +
  plot_annotation(title = 'Target food web', theme = theme(plot.title = element_text(hjust = 0.5)))
ggsave("figures/SI/ModelTransferabilityFull.png", p)

# Correlate transferability to distance metrics ---------------------------
# glm with overall logit-auc and log(aucpr/prevalence) as response
# distances as fixed effect
# source and target region as crossed random effect

# transform responses
overall_performance$logitauc <- log(overall_performance$auc / (1-overall_performance$auc))
overall_performance$logaucpr <- log(overall_performance$aucpr / overall_performance$prevalence)

# scale predictors
overall_performance$geo.dist_sc <- as.vector(scale(overall_performance$geo.dist))
overall_performance$env.dist_sc <- as.vector(scale(overall_performance$env.dist))
overall_performance$phylo.dist_sc <- as.vector(scale(overall_performance$phylo.dist))

# model with logit auc as response and only geographic distance using brms (total effect of geographic distances)
auc_geo_total <- brm(logitauc ~ geo.dist_sc + (1|Source) + (1|Target),
                 data = overall_performance,
                 prior = c(
                   prior(normal(0, 1), class = "Intercept"),
                   prior(normal(0, 1), class = "b"),
                   prior(cauchy(0, 5), class = "sd")
                 ), 
                 sample_prior = "no",
                 iter = 5000)

# model with logit auc as response and geographic distance, controlling for environmental and phylogenetic distances (direct effect of geographic distance)
auc_geo_direct <- brm(logitauc ~ geo.dist_sc + phylo.dist_sc + env.dist_sc + (1|Source) + (1|Target),
                     data = overall_performance,
                     prior = c(
                       prior(normal(0, 1), class = "Intercept"),
                       prior(normal(0, 1), class = "b"),
                       prior(cauchy(0, 5), class = "sd")
                     ), 
                     sample_prior = "no",
                     iter = 5000)

# model with logit auc as response and phylogenetic distance, controlling for geographic distances (effect of phylogenetic distance)
auc_phylo <- brm(logitauc ~ phylo.dist_sc + geo.dist_sc + (1|Source) + (1|Target),
               data = overall_performance,
               prior = c(
                 prior(normal(0, 1), class = "Intercept"),
                 prior(normal(0, 1), class = "b"),
                 prior(cauchy(0, 5), class = "sd")
               ), 
               sample_prior = "no",
               iter = 5000)

# model with logit auc as response and environmental distance, controlling for geographic distances (effect of environmental distance)
auc_env <- brm(logitauc ~ env.dist_sc + geo.dist_sc + (1|Source) + (1|Target),
                 data = overall_performance,
                 prior = c(
                   prior(normal(0, 1), class = "Intercept"),
                   prior(normal(0, 1), class = "b"),
                   prior(cauchy(0, 5), class = "sd")
                 ), 
                 sample_prior = "no",
                 iter = 5000)

# plot geographic distance
geo_fe_direct <- tibble(geo.dist_sc = seq(min(overall_performance$geo.dist_sc), 
                                max(overall_performance$geo.dist_sc), length.out=100), phylo.dist_sc = 0, env.dist_sc = 0, Source = NA, Target = NA) %>%
  add_epred_draws(auc_geo_direct,
                   re_formula = NA, ndraws = 1e4) %>%
  mutate(x = (geo.dist_sc * sd(overall_performance$geo.dist) + mean(overall_performance$geo.dist)),
         y = exp(.epred)/(1+exp(.epred)))

geo_fe_direct <- geo_fe_direct %>% 
  group_by(x)  %>%
  summarize(median_qi(y, width = 0.95)) %>%
  select(x, y, ymin, ymax)

geo_fe_total <- tibble(geo.dist_sc = seq(min(overall_performance$geo.dist_sc), 
                                          max(overall_performance$geo.dist_sc), length.out=100), Source = NA, Target = NA) %>%
  add_epred_draws(auc_geo_total,
                  re_formula = NA, ndraws = 1e4) %>%
  mutate(x = (geo.dist_sc * sd(overall_performance$geo.dist) + mean(overall_performance$geo.dist)),
         y = exp(.epred)/(1+exp(.epred)))

geo_fe_total <- geo_fe_total %>% 
  group_by(x)  %>%
  summarize(median_qi(y, width = 0.95)) %>%
  select(x, y, ymin, ymax)

p1 <- ggplot(geo_fe_direct,
       aes(x = x/1000, y = y)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha= 0.2, data = geo_fe_total, linetype = "dashed", color = "#46BAAA", fill = "#46BAAA") +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), color = "#183446", fill = "#183446", alpha =  0.2) +
  geom_line(data = geo_fe_total, size = 1, linetype = "dashed", color = "#46BAAA") +
  geom_line(size = 1, color = "#183446") +
  geom_point(aes(x = geo.dist/1000, y = auc), data = overall_performance) +
  labs(y = "AUC", x = "Geographic distance (10³km)")+
  theme_minimal() +
  lims(y = c(0.25,1))+
  theme(axis.line = element_line(size = 0.5), strip.text = element_text(size = 12),
        strip.background = element_rect(colour = "black"), panel.grid = element_blank())

# plot environmental distance
env_fe <- tibble(env.dist_sc = seq(min(overall_performance$env.dist_sc), 
                                   max(overall_performance$env.dist_sc), length.out=100), geo.dist_sc = 0, Source = NA, Target = NA) %>%
  add_epred_draws(auc_env,
                  re_formula = NA, ndraws = 1e4) %>%
  mutate(x = (env.dist_sc * sd(overall_performance$env.dist) + mean(overall_performance$env.dist)),
         y = exp(.epred)/(1+exp(.epred)))

env_fe <- env_fe %>% 
  group_by(x)  %>%
  summarize(median_qi(y, width = 0.95)) %>%
  select(x, y, ymin, ymax)

p2 <- ggplot(env_fe,
             aes(x = x, y = y)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), color = "#183446", fill = "#183446", alpha= 0.2) +
  geom_line(size = 1, color = "#183446") +
  geom_point(aes(x = env.dist, y = auc), data = overall_performance) +
  labs(y = "AUC", x = "Environmental distance")+
  lims(y = c(0.25,1))+
  theme_minimal() +
  theme(axis.line = element_line(size = 0.5), panel.grid = element_blank(),
        axis.title.y = element_blank(), axis.text.y = element_blank())

# plot phylogenetic distance
phylo_fe <- tibble(phylo.dist_sc = seq(min(overall_performance$phylo.dist_sc), 
                                   max(overall_performance$phylo.dist_sc), length.out=100), geo.dist_sc = 0, Source = NA, Target = NA) %>%
  add_epred_draws(auc_phylo,
                  re_formula = NA, ndraws = 1e4) %>%
  mutate(x = (phylo.dist_sc * sd(overall_performance$phylo.dist) + mean(overall_performance$phylo.dist)),
         y = exp(.epred)/(1+exp(.epred)))

phylo_fe <- phylo_fe %>% 
  group_by(x) %>%
  summarize(median_qi(y, width = 0.95)) %>%
  select(x, y, ymin, ymax)

p3 <- ggplot(phylo_fe,
             aes(x = x, y = y)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), color = "#183446", fill = "#183446", alpha= 0.2) +
  geom_line(size = 1, color = "#183446") +
  geom_point(aes(x = phylo.dist, y = auc), data = overall_performance) +
  labs(y = "AUC", x = "Phylogenetic distance")+
  lims(y = c(0.25,1))+
  theme_minimal() +
  theme(axis.line = element_line(size = 0.5), panel.grid = element_blank(),
        axis.title.y = element_blank(), axis.text.y = element_blank())

p1 + p2 + p3
ggsave("figures/ModelTransferability.png", width = 8, height = 4)

       