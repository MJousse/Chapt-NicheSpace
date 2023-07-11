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
library(patchwork)
source("code/functions.R")


# FW properties in predicted webs -----------------------------------------
library(foreach)
library(doParallel)
load("data/checkpoints/BRT_predictions.RData")
foodwebs <- c("arctic", "pyrenees", "serengeti", "europe")
combinations <- expand_grid(Source = foodwebs, Target = foodwebs)
predicted_properties <-c()
# for each combination of source and target webs, use 100 posterior sample
# calculate the roles of all species in these sample, and extract the mean
for (combination in c(1:nrow(combinations))){
  sourceFW <- as.character(combinations[combination, "Source"])
  targetFW <- as.character(combinations[combination, "Target"])
  print(Sys.time())
  print(paste0(sourceFW, " model predicting the ", targetFW, " food web..."))
  predictions <- get(paste0(sourceFW, "BRT_", targetFW, "FW"))
  cl <- makeCluster(16) 
  registerDoParallel(cl)
  properties <- foreach(i=c(1:100), .combine = rbind, 
                        .packages = c("igraph", "NetIndices", "dplyr", "tidyr")) %dopar% {
                          prediction <- data.frame(resource = predictions$Prey, 
                                                   consumer = predictions$Predator, 
                                                   interaction = rbinom(nrow(predictions),2,predictions$prediction))
                          prediction <- prediction[prediction$interaction == 1,]
                          fw_properties(prediction, nsim = 1)
                        }
  stopCluster(cl)
  fw_properties_mean <- apply(properties, MARGIN = 2, mean)
  fw_properties_sd <- apply(properties, MARGIN = 2, sd)
  predicted_properties <- rbind(predicted_properties, 
                                data.frame(t(fw_properties_mean)) %>% pivot_longer(everything(), names_to = "metric", values_to = "predicted") %>% mutate(targetFW = targetFW, sourceFW = sourceFW))
  write.csv(predicted_properties, file = "data/checkpoints/PredictedProperties_BRT.csv")
}

# Compare the empirical roles and predicted roles -------------------------
empirical_properties <- read.csv("data/checkpoints/EmpiricalProperties.csv", row.names = 1)
predicted_properties <- read.csv("data/checkpoints/PredictedProperties_BRT.csv", row.names = 1)

empirical_properties$FW <- tolower(empirical_properties$FW)
empirical_properties$FW[empirical_properties$FW == "euro"] <- "europe"

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
  geom_point(position=position_dodge(width=0.75), shape= 45, size = 3) +
  geom_pointrange(data = filter(fw_properties_summary, metric %in% c("connectance", "n_clusters", "modularity", "diameter")),
                  aes(y = error_mean, ymin = error_min, ymax = error_max, group = insample), position=position_dodge(width=0.75), shape= 21, size = 0.5) +
  scale_color_manual(values =  c("grey50","black")) +
  scale_fill_manual(values = c("white","black")) +
  geom_hline(yintercept = 0)+
  labs(y = "Relative error", x = "", color = "Prediction", fill = "Prediction") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "transparent"), axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)), legend.title = element_blank())

p2 <- ggplot(filter(fw_properties, metric %in% c("motif1", "motif2", "motif4", "motif5")), aes(x = metric, y = error, colour = insample, fill = insample)) +
  geom_point(position=position_dodge(width=0.75), shape= 45, size = 3) +
  geom_pointrange(data = filter(fw_properties_summary, metric %in% c("motif1", "motif2", "motif4", "motif5")),
                  aes(y = error_mean, ymin = error_min, ymax = error_max, group = insample), position=position_dodge(width=0.75), shape= 21, size = 0.5) +
  scale_color_manual(values =  c("grey50","black")) +
  scale_fill_manual(values = c("white","black")) +
  geom_hline(yintercept = 0)+
  labs(y = "", x = "", color = "Prediction", fill = "Prediction") +
  scale_y_continuous(breaks = seq(0,100,50), limits = c(-20,100))+
  theme_bw() +
  theme(strip.background = element_rect(fill = "transparent"), axis.title = element_blank(), legend.title = element_blank())

p1 + p2 + plot_layout(guides = "collect")

ggsave("figures/SI/FWproperties_BRT.png", width = 18, height = 9, units = "cm")
