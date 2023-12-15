# Step 10: Analyze the biases in trophic role predictions
# 1. Calculate the trophic roles of all species in all empirical webs
# Roles are related to centrality, trophic level, module-based and motifs-based.
# 2. Calculate the trophic roles of all species in predicted webs
# For each webs, use 100 posterior samples and use the mean.
# 3. Calulate the errors for each role metrics of all species

library(igraph)
library(dplyr)
library(tidyr)
library(NetIndices)
library(multiweb)
library(ggplot2)
source("code/functions.R")
source("code/functions_motifs.R")

# Roles of all species in all empirical webs ------------------------------
# arctic
arcticFW <- read.csv("data/cleaned/HighArcticFW.csv", row.names = 1)
arcticFW <- arcticFW %>%
  transmute(resource = Prey, consumer = Predator)
arcticRoles <- species_role(arcticFW, ncores = 16)

# pyrenees
pyreneesFW <- read.csv("data/cleaned/pyrenneesFW.csv", row.names = 1)
pyreneesFW <- pyreneesFW %>%
  transmute(resource = Prey, consumer = Predator)
pyreneesRoles <- species_role(pyreneesFW, ncores = 16)

# serengeti
serengetiFW <- read.csv("data/cleaned/SerengetiFW.csv", row.names = 1)
serengetiFW <- serengetiFW %>%
  transmute(resource = Resource_Species, consumer = Consumer_Species) %>%
  drop_na()
serengetiRoles <- species_role(serengetiFW, ncores = 16)

# europe
europeMW <- read.csv("data/cleaned/EuroFWadults.csv", row.names = 1)
europeMW <- europeMW %>%
  transmute(resource = Prey, consumer = Predator)
europeRoles <- species_role(europeMW, ncores = 16)

# bind and save
empirical_roles <- bind_rows(
  arcticRoles %>% pivot_longer(-species, names_to = "role", values_to = "empirical") %>% mutate(FW = "Arctic"),
  pyreneesRoles %>% pivot_longer(-species, names_to = "role", values_to = "empirical") %>% mutate(FW = "Pyrenees"),
  serengetiRoles %>% pivot_longer(-species, names_to = "role", values_to = "empirical") %>% mutate(FW = "Serengeti"),
  europeRoles %>% pivot_longer(-species, names_to = "role", values_to = "empirical") %>% mutate(FW = "Europe")
)
write.csv(empirical_roles, file = "data/checkpoints/SpeciesRole.csv")

# Trophic roles in predicted webs -----------------------------------------
library(foreach)
library(doParallel)
load("data/checkpoints/predictions.RData")
foodwebs <- c("Arctic", "Pyrenees", "Serengeti", "Euro")
combinations <- expand_grid(Source = foodwebs, Target = foodwebs)
predicted_roles <-c()
# for each combination of source and target webs, use 100 posterior sample
# calculate the roles of all species in these sample, and extract the mean
for (combination in c(1:nrow(combinations))){
  sourceFW <- as.character(combinations[combination, "Source"])
  targetFW <- as.character(combinations[combination, "Target"])
  print(Sys.time())
  print(paste0(sourceFW, " model predicting the ", targetFW, " food web..."))
  predictions <- get(paste0(sourceFW, "_", targetFW, "_predictions"))
  cl <- makeCluster(8) 
  registerDoParallel(cl)
  role <- foreach(i=c(1:100), .combine = rbind, 
                  .packages = c("igraph", "NetIndices", "multiweb", "dplyr", "tidyr"),
                  .export = c("motif_role", "positions")) %dopar% {
    prediction <- data.frame(resource = predictions$Prey, 
                              consumer = predictions$Predator, 
                              interaction = predictions[,paste0("draws",i)])
    prediction <- prediction[prediction$interaction == 1,]
    species_role(prediction, ncores = 0, nsim=1)
                  }
  stopCluster(cl)
  role_mean <- group_by(role, species) %>% summarise_all(mean, na.rm = T)
  role_sd <- group_by(role, species) %>% summarise_all(sd, na.rm = T)
  predicted_roles <- rbind(predicted_roles, 
                           pivot_longer(role_mean, -species, names_to = "role", values_to = "predicted") %>% mutate(targetFW = targetFW, sourceFW = sourceFW))
  write.csv(predicted_roles, file = "data/checkpoints/predicted_roles.csv")
}

# Compare the empirical roles and predicted roles -------------------------
empirical_roles <- read.csv("data/checkpoints/SpeciesRole.csv", row.names = 1)
predicted_roles <- read.csv("data/checkpoints/predicted_roles.csv", row.names = 1)
empirical_roles$FW[empirical_roles$FW == "Europe"] <- "Euro"

species_roles <- left_join(predicted_roles, empirical_roles,
                           by = c("species", "role", "targetFW" = "FW")) 

species_roles_stats <- species_roles %>%
  group_by(role, targetFW, sourceFW) %>%
  reframe(predicted_mean = mean(predicted),
          empirical_mean = mean(empirical, na.rm = T),
          predicted_cv = sd(predicted)/mean(predicted),
          empirical_cv = sd(empirical, na.rm = T)/mean(empirical, na.rm = T)) %>%
  mutate(mean_ratio = predicted_mean/empirical_mean,
         cv_ratio = predicted_cv/empirical_cv) %>%
  filter(empirical_mean != 0)

# calculate the slope, intercept and r^2 for all role, targetFW and sourceFW
library(purrr)
library(broom)
library(brms)
fitted_models = species_roles %>% 
  drop_na() %>%
  nest(data = -c(role, sourceFW, targetFW)) %>%
  mutate(
    model = map(data, ~ lm(predicted ~ 1 + empirical, data = .x)),
    tidied = map(model, tidy),
    glanced = map(model, glance),
    augmented = map(model, augment, interval = "confidence")
  )

lm_coef <- fitted_models %>% 
  unnest(tidied) %>%
  dplyr::select(targetFW, sourceFW, role, term, estimate, std.error)
  
goodness_of_fit <- fitted_models %>% 
  unnest(glanced) %>%
  dplyr::select(targetFW, role, sourceFW, r.squared)

fitted <- fitted_models %>% 
  unnest(augmented) %>%
  dplyr::select(targetFW, sourceFW, role, empirical, .fitted, .upper, .lower)

# plot individual regressions
for (irole in unique(species_roles$role)){
  d <- filter(fitted, role == irole) %>%
    mutate(sourceFW = factor(sourceFW, levels = c("Arctic", "Euro", "Pyrenees", "Serengeti")),
           targetFW = factor(targetFW, levels = c("Arctic", "Euro", "Pyrenees", "Serengeti")))
  ggplot(d, aes(x = empirical, y = .fitted)) +
    geom_ribbon(aes(ymin = .lower, ymax = .upper, fill = sourceFW), alpha = 0.5) +
    geom_line(aes(colour = sourceFW)) +
    scale_color_manual(values = c("#595365", "#8ACB88", "#E4572E", "#FFBF46")) +
    scale_fill_manual(values = c("#595365", "#8ACB88", "#E4572E", "#FFBF46")) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    labs(x = "empirical", y = "predicted", colour = "Model", fill = "Model") + 
    facet_wrap(.~targetFW, nrow = 2, scales = "free") +
    theme_classic()
  ggsave(paste0("figures/SI/species_role/", irole, ".png"), scale = 3)
}

roles <-  c("indegree", "outdegree", "betweeness", "closeness", "eigen", "within_module_degree", "among_module_conn", "position1", "position2", "position3", "position4", "position5", "position6", "position8", "position9", "position10", "position11")

# plot R2
goodness_of_fit %>%
  mutate(role = factor(role, levels = rev(unique(goodness_of_fit$role))),
         sourceFW = factor(sourceFW, levels = c("Arctic", "Euro", "Pyrenees", "Serengeti")),
         targetFW = factor(targetFW, levels = c("Arctic", "Euro", "Pyrenees", "Serengeti"))) %>%
  filter(role %in% roles) %>%
  ggplot() +
  geom_point(aes(y = role, x = r.squared, color = sourceFW), alpha = 0.75, size = 2) +
  scale_color_manual(values = c("#595365", "#8ACB88", "#E4572E", "#FFBF46")) +
  labs(x = "R²", y = "Role", colour = "Model") +
  facet_wrap(.~targetFW, nrow = 1, scales = "free_x") +
  theme_classic() +
  theme(panel.grid.major.y = element_line())

ggsave(paste0("figures/SI/species_role/R2.png"), scale = 3)

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

ggsave("figures/SpeciesRolePerformance.png", dpi = 600, width = 18, height = 9, units = "cm")

# save coefficients and R²
lm_coef <- pivot_wider(lm_coef, names_from = term, values_from = c(estimate, std.error))
lm_results <- left_join(lm_coef, goodness_of_fit) %>%
  pivot_wider(names_from = targetFW, values_from = c("estimate_(Intercept)", "estimate_empirical", "std.error_(Intercept)", "std.error_empirical", "r.squared")) %>%
  arrange(role, sourceFW)

lm_results <- lm_results[,c(2,1,6,14,5,13,3,11,4,12,10,18,9,17,7,15,8,16,22,21,19,20)]
write.csv(lm_results, file = "data/checkpoints/lm_role_results.csv")


# using the alternative serengeti food web
load("~/../OneDrive - McGill University/Chapt-NicheSpace/predictions_serengeti_alt.RData")
combinations <- filter(combinations, Source == "Serengeti" | Target == "Serengeti")
combinations[combinations == "Serengeti"] <- 'Serengeti2'
# calculate the roles of all species in these sample, and extract the mean
for (combination in c(1:nrow(combinations))){
  sourceFW <- as.character(combinations[combination, "Source"])
  targetFW <- as.character(combinations[combination, "Target"])
  print(Sys.time())
  print(paste0(sourceFW, " model predicting the ", targetFW, " food web..."))
  predictions <- get(paste0(sourceFW, "_", targetFW, "_predictions"))
  cl <- makeCluster(8) 
  registerDoParallel(cl)
  role <- foreach(i=c(1:100), .combine = rbind, 
                  .packages = c("igraph", "NetIndices", "multiweb", "dplyr", "tidyr"),
                  .export = c("motif_role", "positions")) %dopar% {
                    prediction <- data.frame(resource = predictions$Prey, 
                                             consumer = predictions$Predator, 
                                             interaction = predictions[,paste0("draws",i)])
                    prediction <- prediction[prediction$interaction == 1,]
                    species_role(prediction, ncores = 0, nsim=1)
                  }
  stopCluster(cl)
  role_mean <- group_by(role, species) %>% summarise_all(mean, na.rm = T)
  role_sd <- group_by(role, species) %>% summarise_all(sd, na.rm = T)
  predicted_roles <- rbind(predicted_roles, 
                           pivot_longer(role_mean, -species, names_to = "role", values_to = "predicted") %>% mutate(targetFW = targetFW, sourceFW = sourceFW))
  write.csv(predicted_roles, file = "data/checkpoints/predicted_roles_seralt.csv")
}
