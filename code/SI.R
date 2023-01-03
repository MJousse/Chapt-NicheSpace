library(dplyr)
library(tidyr)
library(ggplot2)

# Appendix other performance metrics --------------------------------------
foodwebs <- c("Euro", "Pyrenees", "Arctic", "Serengeti")
foodwebs_labs <- c("Europe", "Pyrenees", "NorthQC", "Serengeti")
overall_performance <- read.csv("data/checkpoints/overall_performance_draws.csv", row.names = 1) %>%
  pivot_longer(cols = tpr:npv, names_to = "metric") %>%
  mutate(Source = factor(Source, levels = foodwebs, labels = foodwebs_labs), Target = factor(Target, levels = foodwebs, labels = foodwebs_labs), metric = factor(metric, levels = c("tpr", "tnr", "ppv", "npv")))

ggplot(overall_performance) +
  geom_density(aes(value, after_stat(scaled) ,fill = Source), alpha = .5) +
  facet_grid(Target~metric, scales = "free") +
  scale_fill_manual(values = c("royalblue4", "red3", "deepskyblue", "chartreuse4")) +
  labs(y = "Scaled density") +
  theme_classic()

ggsave("figures/SI/performance_metrics.png", width = 9, dpi = 600)

# how correlated is auc to aucpr and other metrics
library(corrplot)
overall_performance_othermetrics <- read.csv("data/checkpoints/overall_performance_draws.csv", row.names = 1) %>%
  group_by(Source, Target) %>%
  summarise_all(median)

overall_performance_auc <- read.csv("data/checkpoints/overall_performance.csv", row.names=1)

overall_performance <- full_join(overall_performance_othermetrics, overall_performance_auc)

Movr <- cor(overall_performance[,c("auc", "aucpr", "tpr", "tnr", "ppv", "npv")])

png("figures/SI/metricCorr_overall.png", width = 300, height = 300)
corrplot(Movr, method = "number", type = "upper", diag = F)
dev.off()

species_performance_othermetrics <- read.csv("data/checkpoints/species_performance_draws.csv", row.names = 1) %>%
  group_by(Source, Target, species) %>%
  summarise_all(median, na.rm = T)

species_performance_auc <- read.csv("data/checkpoints/species_performance.csv", row.names=1)

species_performance <- full_join(species_performance_othermetrics, species_performance_auc)
Msp <- cor(species_performance[,c("auc", "aucpr", "tpr", "tnr", "ppv", "npv")])

png("figures/SI/metricCorr_species.png", width = 300, height = 300)
corrplot(Msp, method = "number", type = "upper", diag = F)
dev.off()
