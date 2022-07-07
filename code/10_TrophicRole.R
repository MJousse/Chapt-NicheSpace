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
  cl <- makeCluster(2) 
  registerDoParallel(cl)
  role <- foreach(i=c(1:100), .combine = rbind, 
                  .packages = c("igraph", "NetIndices", "multiweb", "dplyr", "tidyr"),
                  .export = c("motif_role", "positions")) %dopar% {
    prediction <- data.frame(resource = predictions$Prey, 
                              consumer = predictions$Predator, 
                              interaction = predictions[,paste0("draws",i)])
    prediction <- prediction[prediction$interaction == 1,]
    species_role(prediction, ncores = 0)
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

species_roles <- left_join(predicted_roles, empirical_roles,
                           by = c("species", "role", "targetFW" = "FW")) %>%
  drop_na() %>%
  group_by(role, targetFW, sourceFW) %>%
  mutate(predicted_scaled = (predicted - mean(empirical, na.rm = T))/sd(empirical, na.rm = T),
         empirical_scaled = (empirical - mean(empirical, na.rm = T))/sd(empirical, na.rm = T))

# calculate the slope, intercept and r^2 for all role, targetFW and sourceFW
library(purrr)
library(broom)
fitted_models = species_roles %>% 
  drop_na() %>%
  nest(data = -c(role, targetFW, sourceFW)) %>%
  mutate(model = map(data, ~ lm(predicted_scaled ~ empirical_scaled, data = .x)),
         tidied = map(model, tidy)
  ) %>%
  unnest(tidied) %>%
  dplyr::select(role, targetFW, sourceFW, term, estimate, std.error)

correlations <- species_roles %>% 
  drop_na() %>%
  group_by(role, targetFW, sourceFW) %>%
  summarise(correlation = cor(predicted, empirical))

# plot
correlations$role <- factor(correlations$role, levels = rev(unique(empirical_roles$role)))
p<-ggplot(subset(correlations, role %in% c("indegree", "outdegree", "betweeness", "closeness", "eigen", "TL", "OI", "within_module_degree", "among_module_conn", "position1", "position2", "position3", "position4", "position5", "position6", "position8", "position9", "position10", "position11")), aes(x = role, y = correlation, fill = sourceFW)) +
  geom_point(shape = 21, size = 3, alpha=0.8) +
  scale_fill_manual(values =  c("deepskyblue","royalblue4", "red3", "chartreuse4")) +
  coord_flip() +
  geom_hline(yintercept = 0)+
  facet_grid(~targetFW) +
  labs(y = "Correlation", x = "Species role", fill = "Model", title = "Predicted food web") +
  ylim(c(-0.5,1))+
  theme_bw() +
  theme(strip.background = element_rect(fill = "transparent"), plot.title = element_text(hjust = 0.5))

ggsave("figures/SpeciesRoleCorrelation.png", dpi = 600)

intercepts <- filter(fitted_models, term == "(Intercept)")
intercepts$role <- factor(intercepts$role, levels = rev(unique(empirical_roles$role)))

ggplot(subset(intercepts, role %in% c("indegree", "outdegree", "betweeness", "closeness", "eigen", "TL", "OI", "within_module_degree", "among_module_conn", "position1", "position2", "position3", "position4", "position5", "position6", "position8", "position9", "position10", "position11")), aes(x = role, y = estimate, fill = sourceFW)) +
  geom_point(shape = 21, size = 4, alpha=0.8) +
  scale_fill_manual(values =  c("deepskyblue","royalblue4", "red3", "chartreuse4")) +
  coord_flip() +
  scale_y_continuous(limits= c(-10, 20)) +
  geom_hline(yintercept = 0)+
  facet_grid(~targetFW) +
  labs(y = "Intercept", x = "Species role") +
  theme_bw()

slopes <- filter(fitted_models, term == "empirical_scaled")
slopes$role <- factor(slopes$role, levels = levels(intercepts$role))
ggplot(subset(slopes, role %in% c("indegree", "outdegree", "betweeness", "closeness", "eigen", "TL", "OI", "within_module_degree", "among_module_conn", "position1", "position2", "position3", "position4", "position5", "position6", "position8", "position9", "position10", "position11")), aes(x = role, y = estimate, fill = sourceFW)) +
  geom_point(shape = 21, size = 4, alpha=0.8) +
  scale_fill_manual(values =  c("deepskyblue","royalblue4", "red3", "chartreuse4")) +
  scale_y_continuous(limits = c(-3,35))+
  coord_flip() +
  geom_hline(yintercept = 0)+
  facet_grid(~targetFW) +
  labs(y = "Slope", x = "Species role") +
  theme_bw()
