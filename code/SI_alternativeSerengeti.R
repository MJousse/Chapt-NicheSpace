rm(list = ls())
library(ROCR)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(purrr)
library(brms)
library(broom)
library(tidybayes)
library(prg)
foodwebs <- c("Arctic", "Euro", "Pyrenees", "Serengeti")
source("code/functions.R")

######### No Serengeti results ################
species_performance <- read.csv("data/checkpoints/species_performance.csv", row.names = 1) %>%
  filter(Source != "Serengeti", Target != "Serengeti")
overall_performance <- read.csv("data/checkpoints/overall_performance.csv", row.names = 1) %>%
  filter(Source != "Serengeti", Target != "Serengeti")
FWdist <- read.csv("data/checkpoints/FWdist.csv", row.names = 1)
FWdist[!is.na(FWdist) & FWdist == "Nunavik"] <- "Arctic"
overall_performance[overall_performance == "Euro"] <- "Europe"
species_performance[species_performance == "Euro"] <- "Europe"
overall_performance <- left_join(overall_performance, FWdist, by = c("Source" = "FW2", "Target" = "FW1"))
species_performance <- left_join(species_performance, FWdist, by = c("Source" = "FW2", "Target" = "FW1"))

### Figure 2 ###
# transform responses
overall_performance$logitauc <- log(overall_performance$auc / (1-overall_performance$auc))

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
  dplyr::select(x, y, ymin, ymax)

geo_fe_total <- tibble(geo.dist_sc = seq(min(overall_performance$geo.dist_sc), 
                                         max(overall_performance$geo.dist_sc), length.out=100), Source = NA, Target = NA) %>%
  add_epred_draws(auc_geo_total,
                  re_formula = NA, ndraws = 1e4) %>%
  mutate(x = (geo.dist_sc * sd(overall_performance$geo.dist) + mean(overall_performance$geo.dist)),
         y = exp(.epred)/(1+exp(.epred)))

geo_fe_total <- geo_fe_total %>% 
  group_by(x)  %>%
  summarize(median_qi(y, width = 0.95)) %>%
  dplyr::select(x, y, ymin, ymax)

p1 <- ggplot(geo_fe_direct,
             aes(x = x/1000, y = y)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha= 0.2, data = geo_fe_total, linetype = "dashed", color = "#46BAAA", fill = "#46BAAA") +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), color = "#183446", fill = "#183446", alpha =  0.2) +
  geom_line(data = geo_fe_total, size = 1, linetype = "dashed", color = "#46BAAA") +
  geom_line(size = 1, color = "#183446") +
  geom_point(aes(x = geo.dist/1000, y = auc), data = overall_performance) +
  labs(y = "AUC", x = "Geographic distance (10³km)")+
  theme_minimal() +
  lims(y = c(0,1))+
  theme(axis.line = element_line(linewidth = 0.5), strip.text = element_text(size = 12),
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
  dplyr::select(x, y, ymin, ymax)

p2 <- ggplot(env_fe,
             aes(x = x, y = y)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), color = "#183446", fill = "#183446", alpha= 0.2) +
  geom_line(size = 1, color = "#183446") +
  geom_point(aes(x = env.dist, y = auc), data = overall_performance) +
  labs(y = "AUC", x = "Environmental distance")+
  lims(y = c(0,1))+
  theme_minimal() +
  theme(axis.line = element_line(linewidth = 0.5), panel.grid = element_blank(),
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
  dplyr::select(x, y, ymin, ymax)

p3 <- ggplot(phylo_fe,
             aes(x = x, y = y)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), color = "#183446", fill = "#183446", alpha= 0.2) +
  geom_line(size = 1, color = "#183446") +
  geom_point(aes(x = phylo.dist, y = auc), data = overall_performance) +
  labs(y = "AUC", x = "Phylogenetic distance")+
  lims(y = c(0,1))+
  theme_minimal() +
  theme(axis.line = element_line(linewidth = 0.5), panel.grid = element_blank(),
        axis.title.y = element_blank(), axis.text.y = element_blank())

p1 + p2 + p3
ggsave("figures/SI/ModelTransferability_noser.png", width = 8, height = 4)

### figure 3 ###
species_performance <- fread("data/checkpoints/species_performance_trans.csv") %>%
  filter(Source != "Serengeti", Target != "Serengeti")

# transform responses
species_performance$logitauc <- log(species_performance$auc / (1.001-species_performance$auc))

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
  dplyr::select(x, y, ymin, ymax)

p1 <- ggplot(phylo_fe,
             aes(x = x, y = y)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), color = "#183446", fill = "#183446", alpha= 0.2) +
  geom_line(size = 1, color = "#183446") +
  geom_point(aes(x = mntd, y = auc), data = species_performance, size = 0.2, alpha = 0.3) +
  labs(y = "AUC", x = "Distance to nearest taxon")+
  theme_minimal() +
  lims(y = c(0,1))+
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
  dplyr::select(x, y, ymin, ymax)

p2 <- ggplot(func_fe,
             aes(x = x, y = y)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), color = "#183446", fill = "#183446", alpha= 0.2) +
  geom_line(size = 1, color = "#183446") +
  geom_point(aes(x = fmpd, y = auc), data = species_performance, size = 0.2, alpha = 0.3) +
  labs(y = "AUC", x = "Mean functional distance")+
  theme_minimal() +
  lims(y = c(0,1)) +
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
  dplyr::select(x, y, ymin, ymax)

p3 <- ggplot(prev_fe,
             aes(x = x, y = y)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), color = "#183446", fill = "#183446", alpha= 0.2) +
  geom_line(size = 1, color = "#183446") +
  geom_point(aes(x = prevalence, y = auc), data = species_performance, size = 0.2, alpha = 0.3) +
  labs(y = "AUC", x = "Normalized degree")+
  theme_minimal() +
  lims(y = c(0,1)) +
  theme(axis.line = element_line(size = 0.5), panel.grid = element_blank(),
        axis.title.y = element_blank(), axis.text.y = element_blank())

p1 + p2 + p3

ggsave("figures/SI/SpeciesPerformance_noser.png", width = 8, height = 4)

#### figure 4 ###
empirical_roles <- read.csv("data/checkpoints/SpeciesRole.csv", row.names = 1)
predicted_roles <- read.csv("data/checkpoints/predicted_roles.csv", row.names = 1)
empirical_roles$FW[empirical_roles$FW == "Europe"] <- "Euro"

species_roles <- left_join(predicted_roles, empirical_roles,
                           by = c("species", "role", "targetFW" = "FW")) %>%
  filter(targetFW != "Serengeti", sourceFW != "Serengeti")
roles <-  c("indegree", "outdegree", "betweeness", "closeness", "eigen", "within_module_degree", "among_module_conn", "position1", "position2", "position3", "position4", "position5", "position6", "position8", "position9", "position10", "position11")

fitted_models = species_roles %>% 
  drop_na() %>%
  nest(data = -c(role, sourceFW, targetFW)) %>%
  mutate(
    model = map(data, ~ lm(predicted ~ 1 + empirical, data = .x)),
    tidied = map(model, tidy),
    glanced = map(model, glance),
    augmented = map(model, augment, interval = "confidence")
  )
goodness_of_fit <- fitted_models %>% 
  unnest(glanced) %>%
  dplyr::select(targetFW, role, sourceFW, r.squared)

goodness_of_fit$insample <- factor(goodness_of_fit$targetFW == goodness_of_fit$sourceFW,
                                   levels = c(T,F), labels = c("within food web", "between food webs"))

goodness_of_fit2 <- goodness_of_fit %>%
  filter(role %in% roles) %>%
  mutate(role = factor(role, levels = roles))

r2_summary <- goodness_of_fit %>%
  filter(role %in% roles) %>%
  mutate(role = factor(role, levels = roles)) %>%
  group_by(role, insample) %>%
  summarise(r2_mean = mean(r.squared, na.rm = T), r2_min = min(r.squared, na.rm = T), r2_max = max(r.squared, na.rm = T))

ggplot(r2_summary, aes(x = role, y = r2_mean, colour = insample, fill = insample)) +
  geom_point(data = goodness_of_fit2, aes(y = r.squared, colour = insample, fill = insample), position=position_dodge(width=0.75), shape= 45, size = 6) + 
  geom_pointrange(aes(ymin = r2_min, ymax = r2_max, group = insample), position=position_dodge(width=0.75), shape= 21, size = 0.5) +
  scale_color_manual(values =  c("grey50","black")) +
  scale_fill_manual(values = c("white","black")) +
  geom_hline(yintercept = 0)+
  labs(y = "R²", x = "Species role", color = "Prediction", fill = "Prediction") +
  ylim(c(0,1))+
  theme_bw() +
  theme(strip.background = element_rect(fill = "transparent"), axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
        legend.title = element_blank(), axis.text.x = element_blank())

ggsave("figures/SI/SpeciesRolePerformance_noser.png", dpi = 600, width = 18, height = 9, units = "cm")

### Figure 5 ###
empirical_properties <- read.csv("data/checkpoints/EmpiricalProperties.csv", row.names = 1)
predicted_properties <- read.csv("data/checkpoints/PredictedProperties.csv", row.names = 1)

fw_properties <- left_join(predicted_properties, empirical_properties,
                           by = c("metric", "targetFW" = "FW")) %>%
  filter(targetFW != "Serengeti", sourceFW != "Serengeti") %>%
  drop_na() %>%
  mutate(error = (predicted-empirical)/empirical,
         insample = (targetFW == sourceFW))

fw_properties$metric <- factor(fw_properties$metric, levels = c("connectance", "maxTL", "meanTL", "n_clusters", "modularity", "diameter", paste0("motif", c(1:13))))
fw_properties$insample <- factor(fw_properties$insample, levels = c(T,F), labels = c("within food web", "between food webs"))

fw_properties_summary <- fw_properties %>%
  group_by(metric, insample) %>%
  summarise(error_mean = mean(error, na.rm = T), error_min = min(error, na.rm = T), error_max = max(error, na.rm = T))


p1 <- ggplot(filter(fw_properties, metric %in% c("connectance", "maxTL", "meanTL", "n_clusters", "modularity", "diameter")), aes(x = metric, y = error, colour = insample, fill = insample)) +
  geom_point(position=position_dodge(width=0.75), shape= 45, size = 6) +
  geom_pointrange(data = filter(fw_properties_summary, metric %in% c("connectance", "maxTL", "meanTL", "n_clusters", "modularity", "diameter")),
                  aes(y = error_mean, ymin = error_min, ymax = error_max, group = insample), position=position_dodge(width=0.75), shape= 21, size = 0.5) +
  scale_color_manual(values =  c("grey50","black")) +
  scale_fill_manual(values = c("white","black")) +
  geom_hline(yintercept = 0)+
  labs(y = "Relative error", x = "", color = "Prediction", fill = "Prediction") +
  scale_y_continuous(breaks = seq(0,15,5), limits = c(-1,15))+
  theme_bw() +
  theme(strip.background = element_rect(fill = "transparent"), axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)), legend.title = element_blank())

p2 <- ggplot(filter(fw_properties, metric %in% c("motif1", "motif2", "motif4", "motif5")), aes(x = metric, y = error, colour = insample, fill = insample)) +
  geom_point(position=position_dodge(width=0.75), shape= 45, size = 6) +
  geom_pointrange(data = filter(fw_properties_summary, metric %in% c("motif1", "motif2", "motif4", "motif5")),
                  aes(y = error_mean, ymin = error_min, ymax = error_max, group = insample), position=position_dodge(width=0.75), shape= 21, size = 0.5) +
  scale_color_manual(values =  c("grey50","black")) +
  scale_fill_manual(values = c("white","black")) +
  geom_hline(yintercept = 0)+
  labs(y = "", x = "", color = "Prediction", fill = "Prediction") +
  scale_y_continuous(breaks = seq(0,300,100), limits = c(-20,300))+
  theme_bw() +
  theme(strip.background = element_rect(fill = "transparent"), axis.title = element_blank(), legend.title = element_blank())

p1 + p2 + plot_layout(guides = "collect")

ggsave("figures/SI/FWproperties_noser.png", width = 18, height = 9, units = "cm", dpi = 600)

############ Alternative Serengeti results #####################
rm(list = ls())
load("~/OneDrive/Chapt-NicheSpace/predictions.RData")
load("~/OneDrive/Chapt-NicheSpace/predictions_serengeti_alt.RData")
foodwebs <- c("Arctic", "Euro", "Pyrenees", "Serengeti2")
source("code/functions.R")
Serengeti2_Serengeti2_predictions <- Serengeti2_Serengeti_predictions

overall_performance <- expand_grid(Source = foodwebs, Target = foodwebs) %>%
  mutate(auc = NA, aucpr = NA, prevalence = NA, aucprg = NA)

# for each combinaison calculate the overall auc and aucpr
# in each target food web, calculate the auc and aucpr of each species
for (combination in c(1:nrow(overall_performance))){
  sourceFW <- overall_performance[combination, "Source"]
  targetFW <- overall_performance[combination, "Target"]
  predictions <- get(paste0(sourceFW, "_", targetFW, "_predictions"))
  if (sourceFW == targetFW){
    predictions_test <- filter(predictions, testing == 1)
  } else {
    predictions_test <- predictions
  }
  overall_performance[combination, "auc"] <- performance(prediction(predictions_test$Estimate, predictions_test$interaction), "auc")@y.values[[1]]
  overall_performance[combination, "aucpr"] <- performance(prediction(predictions_test$Estimate, predictions_test$interaction), "aucpr")@y.values[[1]]
  overall_performance[combination, "aucprg"] <- calc_auprg(create_prg_curve(predictions_test$interaction, predictions_test$Estimate))
  overall_performance[combination, "prevalence"] <- sum(predictions_test$interaction)/nrow(predictions_test)
}
write.csv(overall_performance, "data/checkpoints/SI_performance_alternativeSerengeti.csv")

FWdist <- read.csv("data/checkpoints/FWdist.csv", row.names = 1)
overall_performance[overall_performance == "Serengeti2"] <- "Serengeti"
FWdist[!is.na(FWdist) & FWdist == "Nunavik"] <- "Arctic"
overall_performance[overall_performance == "Euro"] <- "Europe"
overall_performance <- left_join(overall_performance, FWdist, by = c("Source" = "FW2", "Target" = "FW1"))

### Figure 2 ###
# transform responses
overall_performance$logitauc <- log(overall_performance$auc / (1-overall_performance$auc))

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
  dplyr::select(x, y, ymin, ymax)

geo_fe_total <- tibble(geo.dist_sc = seq(min(overall_performance$geo.dist_sc), 
                                         max(overall_performance$geo.dist_sc), length.out=100), Source = NA, Target = NA) %>%
  add_epred_draws(auc_geo_total,
                  re_formula = NA, ndraws = 1e4) %>%
  mutate(x = (geo.dist_sc * sd(overall_performance$geo.dist) + mean(overall_performance$geo.dist)),
         y = exp(.epred)/(1+exp(.epred)))

geo_fe_total <- geo_fe_total %>% 
  group_by(x)  %>%
  summarize(median_qi(y, width = 0.95)) %>%
  dplyr::select(x, y, ymin, ymax)

p1 <- ggplot(geo_fe_direct,
             aes(x = x/1000, y = y)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha= 0.2, data = geo_fe_total, linetype = "dashed", color = "#46BAAA", fill = "#46BAAA") +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), color = "#183446", fill = "#183446", alpha =  0.2) +
  geom_line(data = geo_fe_total, size = 1, linetype = "dashed", color = "#46BAAA") +
  geom_line(size = 1, color = "#183446") +
  geom_point(aes(x = geo.dist/1000, y = auc), data = overall_performance) +
  labs(y = "AUC", x = "Geographic distance (10³km)")+
  theme_minimal() +
  lims(y = c(0,1))+
  theme(axis.line = element_line(linewidth = 0.5), strip.text = element_text(size = 12),
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
  dplyr::select(x, y, ymin, ymax)

p2 <- ggplot(env_fe,
             aes(x = x, y = y)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), color = "#183446", fill = "#183446", alpha= 0.2) +
  geom_line(size = 1, color = "#183446") +
  geom_point(aes(x = env.dist, y = auc), data = overall_performance) +
  labs(y = "AUC", x = "Environmental distance")+
  lims(y = c(0,1))+
  theme_minimal() +
  theme(axis.line = element_line(linewidth = 0.5), panel.grid = element_blank(),
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
  dplyr::select(x, y, ymin, ymax)

p3 <- ggplot(phylo_fe,
             aes(x = x, y = y)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), color = "#183446", fill = "#183446", alpha= 0.2) +
  geom_line(size = 1, color = "#183446") +
  geom_point(aes(x = phylo.dist, y = auc), data = overall_performance) +
  labs(y = "AUC", x = "Phylogenetic distance")+
  lims(y = c(0,1))+
  theme_minimal() +
  theme(axis.line = element_line(linewidth = 0.5), panel.grid = element_blank(),
        axis.title.y = element_blank(), axis.text.y = element_blank())

p1 + p2 + p3
ggsave("figures/SI/ModelTransferability_altser.png", width = 8, height = 4)

#### figure 4 ###
source("code/functions_motifs.R")
SerengetiFW2 <- SerengetiFW2 %>%
  filter(interaction == 1) %>%
  transmute(resource = Prey, consumer = Predator) %>%
  drop_na()
serengetiRoles <- species_role(SerengetiFW2, ncores = 4) %>% 
  pivot_longer(-species, names_to = "role", values_to = "empirical") %>% mutate(FW = "Serengeti")

empirical_roles <- read.csv("data/checkpoints/SpeciesRole.csv", row.names = 1) %>%
  filter(FW != "Serengeti") %>%
  rbind(serengetiRoles)
predicted_roles_alt <- read.csv("data/checkpoints/predicted_roles_seralt.csv", row.names = 1)
predicted_roles_alt[predicted_roles_alt == "Serengeti2"] <- "Serengeti"
predicted_roles <- read.csv("data/checkpoints/predicted_roles.csv", row.names = 1) %>%
  filter(targetFW != "Serengeti", sourceFW != "Serengeti") %>%
  rbind(predicted_roles_alt)

empirical_roles$FW[empirical_roles$FW == "Europe"] <- "Euro"

species_roles <- left_join(predicted_roles, empirical_roles,
                           by = c("species", "role", "targetFW" = "FW"))
roles <-  c("indegree", "outdegree", "betweeness", "closeness", "eigen", "within_module_degree", "among_module_conn", "position1", "position2", "position3", "position4", "position5", "position6", "position8", "position9", "position10", "position11")

fitted_models = species_roles %>% 
  drop_na() %>%
  nest(data = -c(role, sourceFW, targetFW)) %>%
  mutate(
    model = map(data, ~ lm(predicted ~ 1 + empirical, data = .x)),
    tidied = map(model, tidy),
    glanced = map(model, glance)
  )
goodness_of_fit <- fitted_models %>% 
  unnest(glanced) %>%
  dplyr::select(targetFW, role, sourceFW, r.squared)

goodness_of_fit$insample <- factor(goodness_of_fit$targetFW == goodness_of_fit$sourceFW,
                                   levels = c(T,F), labels = c("within food web", "between food webs"))

goodness_of_fit2 <- goodness_of_fit %>%
  filter(role %in% roles) %>%
  mutate(role = factor(role, levels = roles))

r2_summary <- goodness_of_fit %>%
  filter(role %in% roles) %>%
  mutate(role = factor(role, levels = roles)) %>%
  group_by(role, insample) %>%
  summarise(r2_mean = mean(r.squared, na.rm = T), r2_min = min(r.squared, na.rm = T), r2_max = max(r.squared, na.rm = T))

ggplot(r2_summary, aes(x = role, y = r2_mean, colour = insample, fill = insample)) +
  geom_point(data = goodness_of_fit2, aes(y = r.squared, colour = insample, fill = insample), position=position_dodge(width=0.75), shape= 45, size = 6) + 
  geom_pointrange(aes(ymin = r2_min, ymax = r2_max, group = insample), position=position_dodge(width=0.75), shape= 21, size = 0.5) +
  scale_color_manual(values =  c("grey50","black")) +
  scale_fill_manual(values = c("white","black")) +
  geom_hline(yintercept = 0)+
  labs(y = "R²", x = "Species role", color = "Prediction", fill = "Prediction") +
  ylim(c(0,1))+
  theme_bw() +
  theme(strip.background = element_rect(fill = "transparent"), axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
        legend.title = element_blank(), axis.text.x = element_blank())

ggsave("figures/SI/SpeciesRolePerformance_altser.png", dpi = 600, width = 18, height = 9, units = "cm")

### Figure 5 ###
serengetiProperties <- fw_properties(SerengetiFW2, nsim = 10) %>%
serengetiProperties <- data.frame(t(serengetiProperties)) %>% pivot_longer(everything(), names_to = "metric", values_to = "empirical") %>% mutate(FW = "Serengeti")
empirical_properties <- read.csv("data/checkpoints/EmpiricalProperties.csv", row.names = 1) %>%
  filter(FW != "Serengeti") %>%
  rbind(serengetiProperties)
predicted_properties_alt <- read.csv("data/checkpoints/PredictedProperties_seralt.csv", row.names = 1)
predicted_properties_alt[predicted_properties_alt == "Serengeti2"] <- "Serengeti"
predicted_properties <- read.csv("data/checkpoints/PredictedProperties.csv", row.names = 1) %>%
  filter(targetFW != "Serengeti", sourceFW != "Serengeti") %>%
  rbind(predicted_properties_alt)

fw_properties <- left_join(predicted_properties, empirical_properties,
                           by = c("metric", "targetFW" = "FW")) %>%
  drop_na() %>%
  mutate(error = (predicted-empirical)/empirical,
         insample = (targetFW == sourceFW))

fw_properties$metric <- factor(fw_properties$metric, levels = c("connectance", "maxTL", "meanTL", "n_clusters", "modularity", "diameter", paste0("motif", c(1:13))))
fw_properties$insample <- factor(fw_properties$insample, levels = c(T,F), labels = c("within food web", "between food webs"))

fw_properties_summary <- fw_properties %>%
  group_by(metric, insample) %>%
  summarise(error_mean = mean(error, na.rm = T), error_min = min(error, na.rm = T), error_max = max(error, na.rm = T))


p1 <- ggplot(filter(fw_properties, metric %in% c("connectance", "n_clusters", "modularity", "diameter")), aes(x = metric, y = error, colour = insample, fill = insample)) +
  geom_point(position=position_dodge(width=0.75), shape= 45, size = 6) +
  geom_pointrange(data = filter(fw_properties_summary, metric %in% c("connectance", "n_clusters", "modularity", "diameter")),
                  aes(y = error_mean, ymin = error_min, ymax = error_max, group = insample), position=position_dodge(width=0.75), shape= 21, size = 0.5) +
  scale_color_manual(values =  c("grey50","black")) +
  scale_fill_manual(values = c("white","black")) +
  geom_hline(yintercept = 0)+
  labs(y = "Relative error", x = "", color = "Prediction", fill = "Prediction") +
  scale_y_continuous(breaks = seq(0,20,5), limits = c(-1,20))+
  theme_bw() +
  theme(strip.background = element_rect(fill = "transparent"), axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)), legend.title = element_blank())

p2 <- ggplot(filter(fw_properties, metric %in% c("motif1", "motif2", "motif4", "motif5")), aes(x = metric, y = error, colour = insample, fill = insample)) +
  geom_point(position=position_dodge(width=0.75), shape= 45, size = 6) +
  geom_pointrange(data = filter(fw_properties_summary, metric %in% c("motif1", "motif2", "motif4", "motif5")),
                  aes(y = error_mean, ymin = error_min, ymax = error_max, group = insample), position=position_dodge(width=0.75), shape= 21, size = 0.5) +
  scale_color_manual(values =  c("grey50","black")) +
  scale_fill_manual(values = c("white","black")) +
  geom_hline(yintercept = 0)+
  labs(y = "", x = "", color = "Prediction", fill = "Prediction") +
  scale_y_continuous(breaks = seq(0,300,100), limits = c(-20,300))+
  theme_bw() +
  theme(strip.background = element_rect(fill = "transparent"), axis.title = element_blank(), legend.title = element_blank())

p1 + p2 + plot_layout(guides = "collect")

ggsave("figures/SI/FWproperties_altser.png", width = 18, height = 9, units = "cm", dpi = 600)
