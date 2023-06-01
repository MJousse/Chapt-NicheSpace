library(dplyr)
library(brms)
library(ggplot2)

rm(list = ls())
species_perf_glm <- read.csv("data/checkpoints/species_performance.csv", row.names = 1) %>% mutate(model = "glmm")
overall_perf_glm <- read.csv("data/checkpoints/overall_performance.csv", row.names = 1) %>% mutate(model = "glmm")
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
