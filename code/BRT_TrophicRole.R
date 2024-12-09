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

# Trophic roles in predicted webs -----------------------------------------
library(foreach)
library(doParallel)
load("data/checkpoints/BRT_predictions.RData")
foodwebs <- c("arctic", "pyrenees", "serengeti", "europe")
combinations <- expand_grid(Source = foodwebs, Target = foodwebs)
predicted_roles <-c()
# for each combination of source and target webs, use 100 posterior sample
# calculate the roles of all species in these sample, and extract the mean
for (combination in c(1:nrow(combinations))){
  sourceFW <- as.character(combinations[combination, "Source"])
  targetFW <- as.character(combinations[combination, "Target"])
  print(Sys.time())
  print(paste0(sourceFW, " model predicting the ", targetFW, " food web..."))
  predictions <- get(paste0(sourceFW, "BRT_", targetFW, "FW"))
  cl <- makeCluster(12) 
  registerDoParallel(cl)
  role <- foreach(i=c(1:100), .combine = rbind, 
                  .packages = c("igraph", "NetIndices", "multiweb", "dplyr", "tidyr"),
                  .export = c("motif_role", "positions")) %dopar% {
                    prediction <- data.frame(resource = predictions$Prey, 
                                             consumer = predictions$Predator, 
                                             interaction = rbinom(nrow(predictions),2,predictions$prediction))
                    prediction <- prediction[prediction$interaction == 1,]
                    species_role(prediction, ncores = 0, nsim=1)
                  }
  stopCluster(cl)
  role_mean <- group_by(role, species) %>% summarise_all(mean, na.rm = T)
  role_sd <- group_by(role, species) %>% summarise_all(sd, na.rm = T)
  predicted_roles <- rbind(predicted_roles, 
                           pivot_longer(role_mean, -species, names_to = "role", values_to = "predicted") %>% mutate(targetFW = targetFW, sourceFW = sourceFW))
  write.csv(predicted_roles, file = "data/checkpoints/predicted_roles_BRT.csv")
}

# Compare the empirical roles and predicted roles -------------------------
empirical_roles <- read.csv("data/checkpoints/SpeciesRole.csv", row.names = 1)
predicted_roles <- read.csv("data/checkpoints/predicted_roles_BRT.csv", row.names = 1)
empirical_roles$FW <- tolower(empirical_roles$FW)

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

roles <-  c("indegree", "outdegree", "betweeness", "closeness", "eigen", "within_module_degree", "among_module_conn", "position1", "position2", "position3", "position4", "position5", "position6", "position8", "position9", "position10", "position11")

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

ggsave("figures/SI/SpeciesRolePerformance_BRT.png", dpi = 600, width = 18, height = 9, units = "cm")

# save coefficients and R²
lm_coef <- pivot_wider(lm_coef, names_from = term, values_from = c(estimate, std.error))
lm_results <- left_join(lm_coef, goodness_of_fit) %>%
  pivot_wider(names_from = targetFW, values_from = c("estimate_(Intercept)", "estimate_empirical", "std.error_(Intercept)", "std.error_empirical", "r.squared")) %>%
  arrange(role, sourceFW)

lm_results <- lm_results[,c(2,1,6,14,5,13,3,11,4,12,10,18,9,17,7,15,8,16,22,21,19,20)]
write.csv(lm_results, file = "tables/SI/BRT_lm_role_results.csv")
