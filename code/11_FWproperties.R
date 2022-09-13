# Step 11: Analyze the biases in food web properties predictions
# 1. Calculate the trophic roles of all species in all empirical webs
# Roles are related to centrality, trophic level, module-based and motifs-based.
# 2. Calculate the trophic roles of all species in predicted webs
# For each webs, use 100 posterior samples and use the mean.
# 3. Calulate the errors for each role metrics of all species
set.seed(16)
library(igraph)
library(dplyr)
library(tidyr)
library(NetIndices)
library(ggplot2)
source("code/functions.R")

# Roles of all species in all empirical webs ------------------------------
# arctic
arcticFW <- read.csv("data/cleaned/HighArcticFW.csv", row.names = 1)
arcticFW <- arcticFW %>%
  transmute(resource = Prey, consumer = Predator)
arcticProperties <- fw_properties(arcticFW, nsim = 10)

# pyrenees
pyreneesFW <- read.csv("data/cleaned/pyrenneesFW.csv", row.names = 1)
pyreneesFW <- pyreneesFW %>%
  transmute(resource = Prey, consumer = Predator)
pyreneesProperties <- fw_properties(pyreneesFW, nsim = 10)

# serengeti
serengetiFW <- read.csv("data/cleaned/SerengetiFW.csv", row.names = 1)
serengetiFW <- serengetiFW %>%
  transmute(resource = Resource_Species, consumer = Consumer_Species) %>%
  drop_na()
serengetiProperties <- fw_properties(serengetiFW, nsim = 10)

# europe
europeMW <- read.csv("data/cleaned/EuroFWadults.csv", row.names = 1)
europeMW <- europeMW %>%
  transmute(resource = Prey, consumer = Predator)
europeProperties <- fw_properties(europeMW, nsim = 10)

# bind and save
empirical_properties <- bind_rows(
  data.frame(t(arcticProperties)) %>% pivot_longer(everything(), names_to = "metric", values_to = "empirical") %>% mutate(FW = "Arctic"),
  data.frame(t(pyreneesProperties)) %>% pivot_longer(everything(), names_to = "metric", values_to = "empirical") %>% mutate(FW = "Pyrenees"),
  data.frame(t(serengetiProperties)) %>% pivot_longer(everything(), names_to = "metric", values_to = "empirical") %>% mutate(FW = "Serengeti"),
  data.frame(t(europeProperties)) %>% pivot_longer(everything(), names_to = "metric", values_to = "empirical") %>% mutate(FW = "Euro")
)
write.csv(empirical_properties, file = "data/checkpoints/EmpiricalProperties.csv")

# FW properties in predicted webs -----------------------------------------
library(foreach)
library(doParallel)
load("data/checkpoints/predictions.RData")
foodwebs <- c("Arctic", "Pyrenees", "Serengeti", "Euro")
combinations <- expand_grid(Source = foodwebs, Target = foodwebs)
predicted_properties <-c()
# for each combination of source and target webs, use 100 posterior sample
# calculate the roles of all species in these sample, and extract the mean
for (combination in c(1:nrow(combinations))){
  sourceFW <- as.character(combinations[combination, "Source"])
  targetFW <- as.character(combinations[combination, "Target"])
  print(Sys.time())
  print(paste0(sourceFW, " model predicting the ", targetFW, " food web..."))
  predictions <- get(paste0(sourceFW, "_", targetFW, "_predictions"))
  cl <- makeCluster(16) 
  registerDoParallel(cl)
  properties <- foreach(i=c(1:100), .combine = rbind, 
                  .packages = c("igraph", "NetIndices", "dplyr", "tidyr")) %dopar% {
                    prediction <- data.frame(resource = predictions$Prey, 
                                             consumer = predictions$Predator, 
                                             interaction = predictions[,paste0("draws",i)])
                    prediction <- prediction[prediction$interaction == 1,]
                    fw_properties(prediction, nsim = 1)
                  }
  stopCluster(cl)
  fw_properties_mean <- apply(properties, MARGIN = 2, mean)
  fw_properties_sd <- apply(properties, MARGIN = 2, sd)
  predicted_properties <- rbind(predicted_properties, 
                                data.frame(t(fw_properties_mean)) %>% pivot_longer(everything(), names_to = "metric", values_to = "predicted") %>% mutate(targetFW = targetFW, sourceFW = sourceFW))
  write.csv(predicted_properties, file = "data/checkpoints/PredictedProperties.csv")
}

# Compare the empirical roles and predicted roles -------------------------
empirical_properties <- read.csv("data/checkpoints/EmpiricalProperties.csv", row.names = 1)
predicted_properties <- read.csv("data/checkpoints/PredictedProperties.csv", row.names = 1)

fw_properties <- left_join(predicted_properties, empirical_properties,
                           by = c("metric", "targetFW" = "FW")) %>%
  drop_na() %>%
  mutate(error = (predicted-empirical)/empirical,
         insample = (targetFW == sourceFW))

fw_properties$metric <- factor(fw_properties$metric, levels = c("connectance", "maxTL", "meanTL", "n_clusters", "modularity", "diameter", paste0("motif", c(1:13))))
fw_properties$insample <- factor(fw_properties$insample, levels = c(T,F), labels = c("within food web", "between food webs"))

p1 <- ggplot(filter(fw_properties, metric %in% c("connectance", "maxTL", "meanTL", "n_clusters", "modularity", "diameter")), aes(x = metric, y = error, colour = insample, fill = insample)) +
  geom_point(position=position_dodge(width=0.75), shape= 21, size = 2, alpha = 0.7) +
  scale_color_manual(values =  c("grey50","black")) +
  scale_fill_manual(values = c("white","black")) +
  geom_hline(yintercept = 0)+
  labs(y = "Relative error", x = "", color = "Prediction", fill = "Prediction") +
  scale_y_continuous(breaks = seq(0,15,5), limits = c(-1,15))+
  theme_bw() +
  theme(strip.background = element_rect(fill = "transparent"), axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)), legend.title = element_blank())

p2 <- ggplot(filter(fw_properties, metric %in% c("motif1", "motif2", "motif4", "motif5")), aes(x = metric, y = error, colour = insample, fill = insample)) +
  geom_point(position=position_dodge(width=0.75), shape= 21, size = 2, alpha = 0.7) +
  scale_color_manual(values =  c("grey50","black")) +
  scale_fill_manual(values = c("white","black")) +
  geom_hline(yintercept = 0)+
  labs(y = "", x = "", color = "Prediction", fill = "Prediction") +
  scale_y_continuous(breaks = seq(0,300,100), limits = c(-20,300))+
  theme_bw() +
  theme(strip.background = element_rect(fill = "transparent"), axis.title = element_blank(), legend.title = element_blank())

p1 + p2 + plot_layout(guides = "collect")

ggsave("figures/FWproperties.png", width = 18, height = 9, units = "cm")

ggplot(filter(fw_properties, !is.infinite(error)), aes(x = metric, y = error, fill = sourceFW)) +
  geom_pointrange(aes(ymin = error, ymax = error, group = insample), position=position_dodge(width=0.75), shape= 21, size = 0.5) +
  scale_fill_manual(values =  c("deepskyblue","royalblue4", "red3", "chartreuse4")) +
  facet_grid(~targetFW, scales = "free") +
  coord_flip() +
  geom_hline(yintercept = 0)+
  labs(y = "Relative error", x = "Property", color = "Prediction", fill = "Prediction") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "transparent"), axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)))
