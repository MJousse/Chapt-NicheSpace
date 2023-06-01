library(dplyr)
library(brms)
library(ggplot2)

rm(list = ls())
species_perf_brt <- read.csv("data/checkpoints/species_brt_performance.csv", row.names = 1) %>% mutate(model = "brt")
overall_perf_brt <- read.csv("data/checkpoints/overall_brt_performance.csv", row.names = 1) %>% mutate(model = "brt")

FWdist <- read.csv("data/checkpoints/FWdist.csv", row.names = 1)
FWdist[!is.na(FWdist) & FWdist == "Nunavik"] <- "Arctic"
FWdist$FW1 <- tolower(FWdist$FW1)
FWdist$FW2 <- tolower(FWdist$FW2)
overall_perf_brt <- left_join(overall_perf_brt, FWdist, by = c("Source" = "FW2", "Target" = "FW1"))
species_perf_brt <- left_join(species_perf_brt, FWdist, by = c("Source" = "FW2", "Target" = "FW1"))

# transform responses
overall_perf_brt$logitauc <- log(overall_perf_brt$auc / (1-overall_perf_brt$auc))
overall_perf_brt$logaucpr <- log(overall_perf_brt$aucpr / overall_perf_brt$prevalence)

# scale predictors
overall_perf_brt$geo.dist_sc <- as.vector(scale(overall_perf_brt$geo.dist))
overall_perf_brt$env.dist_sc <- as.vector(scale(overall_perf_brt$env.dist))
overall_perf_brt$phylo.dist_sc <- as.vector(scale(overall_perf_brt$phylo.dist))

# model with logit auc as response and only geographic distance using brms (total effect of geographic distances)
auc_geo_total <- brm(logitauc ~ geo.dist_sc + (1|Source) + (1|Target),
                     data = overall_perf_brt,
                     prior = c(
                       prior(normal(0, 1), class = "Intercept"),
                       prior(normal(0, 1), class = "b"),
                       prior(cauchy(0, 5), class = "sd")
                     ), 
                     sample_prior = "no",
                     iter = 5000)

# model with logit auc as response and geographic distance, controlling for environmental and phylogenetic distances (direct effect of geographic distance)
auc_geo_direct <- brm(logitauc ~ geo.dist_sc + phylo.dist_sc + env.dist_sc + (1|Source) + (1|Target),
                      data = overall_perf_brt,
                      prior = c(
                        prior(normal(0, 1), class = "Intercept"),
                        prior(normal(0, 1), class = "b"),
                        prior(cauchy(0, 5), class = "sd")
                      ), 
                      sample_prior = "no",
                      iter = 5000)

# model with logit auc as response and phylogenetic distance, controlling for geographic distances (effect of phylogenetic distance)
auc_phylo <- brm(logitauc ~ phylo.dist_sc + geo.dist_sc + (1|Source) + (1|Target),
                 data = overall_perf_brt,
                 prior = c(
                   prior(normal(0, 1), class = "Intercept"),
                   prior(normal(0, 1), class = "b"),
                   prior(cauchy(0, 5), class = "sd")
                 ), 
                 sample_prior = "no",
                 iter = 5000)

# model with logit auc as response and environmental distance, controlling for geographic distances (effect of environmental distance)
auc_env <- brm(logitauc ~ env.dist_sc + geo.dist_sc + (1|Source) + (1|Target),
               data = overall_perf_brt,
               prior = c(
                 prior(normal(0, 1), class = "Intercept"),
                 prior(normal(0, 1), class = "b"),
                 prior(cauchy(0, 5), class = "sd")
               ), 
               sample_prior = "no",
               iter = 5000)

load("data/checkpoints/glm_fixed_effects.RData")

# plot geographic distance
geo_fe_direct_brt <- tibble(geo.dist_sc = seq(min(overall_perf_brt$geo.dist_sc), 
                                          max(overall_perf_brt$geo.dist_sc), length.out=100), phylo.dist_sc = 0, env.dist_sc = 0, Source = NA, Target = NA) %>%
  add_epred_draws(auc_geo_direct,
                  re_formula = NA, ndraws = 1e4) %>%
  mutate(x = (geo.dist_sc * sd(overall_perf_brt$geo.dist) + mean(overall_perf_brt$geo.dist)),
         y = exp(.epred)/(1+exp(.epred)))

geo_fe_direct_brt <- geo_fe_direct_brt %>% 
  group_by(x)  %>%
  summarize(median_qi(y, width = 0.95)) %>%
  select(x, y, ymin, ymax)

geo_fe_total_brt <- tibble(geo.dist_sc = seq(min(overall_perf_brt$geo.dist_sc), 
                                         max(overall_perf_brt$geo.dist_sc), length.out=100), Source = NA, Target = NA) %>%
  add_epred_draws(auc_geo_total,
                  re_formula = NA, ndraws = 1e4) %>%
  mutate(x = (geo.dist_sc * sd(overall_perf_brt$geo.dist) + mean(overall_perf_brt$geo.dist)),
         y = exp(.epred)/(1+exp(.epred)))

geo_fe_total_brt <- geo_fe_total_brt %>% 
  group_by(x)  %>%
  summarize(median_qi(y, width = 0.95)) %>%
  select(x, y, ymin, ymax)

# direct effect
p1 <- ggplot(geo_fe_direct,
             aes(x = x/1000, y = y)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha= 0.2, data = geo_fe_direct_brt, linetype = "dashed", color = "#723D46", fill = "#723D46") +
  geom_line(data = geo_fe_direct_brt, size = 1, linetype = "dashed", color = "#723D46") +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), color = "#183446", fill = "#183446", alpha =  0.2) +
  geom_line(size = 1, color = "#183446") +
  labs(y = "AUC", x = "Geographic distance (10³km)")+
  theme_minimal() +
  lims(y = c(0.25,1))+
  theme(axis.line = element_line(linewidth = 0.5), strip.text = element_text(size = 12),
        strip.background = element_rect(colour = "black"), panel.grid = element_blank())

# plot environmental distance
env_fe_brt <- tibble(env.dist_sc = seq(min(overall_perf_brt$env.dist_sc), 
                                   max(overall_perf_brt$env.dist_sc), length.out=100), geo.dist_sc = 0, Source = NA, Target = NA) %>%
  add_epred_draws(auc_env,
                  re_formula = NA, ndraws = 1e4) %>%
  mutate(x = (env.dist_sc * sd(overall_perf_brt$env.dist) + mean(overall_perf_brt$env.dist)),
         y = exp(.epred)/(1+exp(.epred)))

env_fe_brt <- env_fe_brt %>% 
  group_by(x)  %>%
  summarize(median_qi(y, width = 0.95)) %>%
  select(x, y, ymin, ymax)

p2 <- ggplot(env_fe,
             aes(x = x, y = y)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha= 0.2, data = env_fe_brt, linetype = "dashed", color = "#723D46", fill = "#723D46") +
  geom_line(data = env_fe_brt, size = 1, linetype = "dashed", color = "#723D46") +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), color = "#183446", fill = "#183446", alpha =  0.2) +
  geom_line(size = 1, color = "#183446") +
  labs(y = "AUC", x = "Environmental distance")+
  theme_minimal() +
  lims(y = c(0.25,1))+
  theme(axis.line = element_line(linewidth = 0.5), strip.text = element_text(size = 12),
        strip.background = element_rect(colour = "black"), panel.grid = element_blank())

# plot phylogenetic distance
phylo_fe_brt <- tibble(phylo.dist_sc = seq(min(overall_perf_brt$phylo.dist_sc), 
                                       max(overall_perf_brt$phylo.dist_sc), length.out=100), geo.dist_sc = 0, Source = NA, Target = NA) %>%
  add_epred_draws(auc_phylo,
                  re_formula = NA, ndraws = 1e4) %>%
  mutate(x = (phylo.dist_sc * sd(overall_perf_brt$phylo.dist) + mean(overall_perf_brt$phylo.dist)),
         y = exp(.epred)/(1+exp(.epred)))

phylo_fe_brt <- phylo_fe_brt %>% 
  group_by(x) %>%
  summarize(median_qi(y, width = 0.95)) %>%
  select(x, y, ymin, ymax)

p3 <- ggplot(phylo_fe,
             aes(x = x, y = y)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), color = "#183446", fill = "#183446", alpha= 0.2) +
  geom_line(size = 1, color = "#183446")  +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha= 0.2, data = phylo_fe_brt, linetype = "dashed", color = "#723D46", fill = "#723D46") +
  geom_line(data = phylo_fe_brt, size = 1, linetype = "dashed", color = "#723D46") +
  labs(y = "AUC", x = "Phylogenetic distance")+
  lims(y = c(0.25,1))+
  theme_minimal() +
  theme(axis.line = element_line(linewidth = 0.5), panel.grid = element_blank(),
        axis.title.y = element_blank(), axis.text.y = element_blank())
