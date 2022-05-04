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
load("data/checkpoints/predictions.RData")
foodwebs <- c("Arctic", "Euro", "Pyrenees", "Serengeti")

# Calculate performance ---------------------------------------------------
overall_performance <- expand_grid(Source = foodwebs, Target = foodwebs) %>%
  mutate(auc = NA, aucpr = NA, prevalence = NA)
species_performance <- data.frame()

# for each combinaison calculate the overall auc and aucpr
# in each target food web, calculate the auc and aucpr of each species
for (combination in c(1:nrow(overall_performance))){
  sourceFW <- overall_performance[combination, "Source"]
  targetFW <- overall_performance[combination, "Target"]
  predictions <- get(paste0(sourceFW, "_", targetFW, "_predictions"))
  overall_performance[combination, "auc"] <- performance(prediction(predictions$Estimate, predictions$interaction), "auc")@y.values[[1]]
  overall_performance[combination, "aucpr"] <- performance(prediction(predictions$Estimate, predictions$interaction), "aucpr")@y.values[[1]]
  overall_performance[combination, "prevalence"] <- sum(predictions$interaction)/nrow(predictions)
  fw_sp <- unique(predictions$Predator)
  for (species in fw_sp){
   sp_predictions <- filter(predictions, Predator == species | Prey == species)
   if(sum(sp_predictions$interaction != 0)){
     sp_auc <- performance(prediction(sp_predictions$Estimate, sp_predictions$interaction), "auc")@y.values[[1]]
     sp_aucpr <- performance(prediction(sp_predictions$Estimate, sp_predictions$interaction), "aucpr")@y.values[[1]]
     species_performance <- rbind(species_performance,
                                  data.frame(sourceFW, targetFW, species, auc = sp_auc, aucpr = sp_aucpr, prevalence = sum(sp_predictions$interaction)/nrow(sp_predictions)))
   }
  }
}
write.csv(species_performance, "data/checkpoints/species_performance.csv")
write.csv(overall_performance, "data/checkpoints/overall_performance.csv")

# Visualize transferability -----------------------------------------------
# add dissimilarity
species_performance <- read.csv("data/checkpoints/species_performance.csv", row.names = 1)
FWdist <- read.csv("data/checkpoints/FWdist.csv", row.names = 1)
FWdist[!is.na(FWdist) & FWdist == "High Arctic"] <- "Arctic"
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
  geom_jitter(aes(x = phylo.dist, y = auc, colour = Source), alpha = 0.2, height = 0, width = 3) +
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
ggsave("figures/ModelTransferability.png", p)

# Correlate transferability to distance metrics ---------------------------
# glm with overall logit-auc and log(aucpr/prevalence) as response
# distances as fixed effect
# source and target region as crossed random effect

# transform responses
overall_performance$logitauc <- log(overall_performance$auc / (1-overall_performance$auc))
overall_performance$logaucpr <- log(overall_performance$aucpr / overall_performance$prevalence)

# scale predictors
overall_performance$geo.dist <- as.vector(scale(overall_performance$geo.dist))
overall_performance$env.dist <- as.vector(scale(overall_performance$env.dist))
overall_performance$phylo.dist <- as.vector(scale(overall_performance$phylo.dist))

# model with logit auc as response using brms
auc_model <- brm(logitauc ~ geo.dist + env.dist + phylo.dist + (1|Source) + (1|Target),
    data = overall_performance,
    prior = c(
      prior(normal(0, 1), class = "Intercept"),
      prior(normal(0, 1), class = "b"),
      prior(cauchy(0, 5), class = "sd")
    ), 
    sample_prior = "no",
    iter = 1000)

# model with log aucpr/prevalence as response using brms
aucpr_model <- brm(logaucpr ~ geo.dist + env.dist + phylo.dist + (1|Source) + (1|Target),
                 data = overall_performance,
                 prior = c(
                   prior(normal(0, 1), class = "Intercept"),
                   prior(normal(0, 1), class = "b"),
                   prior(cauchy(0, 5), class = "sd")
                 ), 
                 sample_prior = "no",
                 iter = 1000)