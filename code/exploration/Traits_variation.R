rm(list = ls())
library(lme4)
library(nnet)
library(dplyr)
library(performance)
library(tidyr)
library(ggplot2)

EuroTraits <- read.csv("data/cleaned/SpeciesTraitsFull.csv", row.names = 1, stringsAsFactors = T)

traits_col <- c(2:15, 20:24)

df <- data.frame()
for (itrait in traits_col){
  trait <- colnames(EuroTraits)[itrait]
  formula <- as.formula(paste0(trait, " ~ (1 | Class/Order/Family/Genus)"))
  if (length(unique(EuroTraits[,itrait])) == 2){
    model <- brm(formula, 
                 data = EuroTraits, family = bernoulli(), cores = 4)
    var_partition <- icc(model, by_group = T)
  }
  else {
    model <- brm(formula, data = EuroTraits, family = gaussian(), cores = 4)
    var_partition <- icc(model, by_group = T)
  }
  df <- rbind(df,
              c(trait, var_partition$ICC, (1 - sum(var_partition$ICC))))
}
colnames(df) <- c("Trait", "Class", "Order", "Family", "Genus", "Residuals")
df[,c(2:6)] <- as.numeric(df[,c(2:6)])

df %>% pivot_longer(-Trait, names_to = "Level", values_to = "ICC") %>%
  mutate(ICC = as.numeric(ICC),
         Trait = factor(Trait, levels = unique(df$Trait)),
         Level = factor(Level, levels = c("Residuals", "Genus", "Family", "Order", "Class"))) %>%
  ggplot()+
  geom_col(aes(fill = Level, x = Trait, y = ICC)) +
  scale_fill_brewer(type = "qual", palette = 6) +
  theme_default() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank())

# How is the food web structured?
EuroMW <- read.csv("data/cleaned/EuroFW.csv", row.names = 1)
EuroTaxo <- read.csv("data/cleaned/EuroMWTaxo.csv", row.names = 1)

as_adjancy <- function(interaction_list){
  Slist <- unique(c(interaction_list[,1], interaction_list[,2]))
  A <- matrix(0, nrow = length(Slist), ncol = length(Slist))
  colnames(A) <- Slist
  rownames(A) <- Slist
  for (i in c(1:nrow(interaction_list))){
    A[interaction_list[i,1], interaction_list[i,2]] <- 1
  }
  return(A)
}

AICgroup <- function(adjancy_matrix, group){
  S <- nrow(adjancy_matrix)
  k <- length(unique(group))
  logLik <- 0
  for (i in unique(group)){
    Si <- sum(group == i)
    for (j in unique(group)){
      Sj <- sum(group == j)
      Lij <- sum(adjancy_matrix[which(group == i), which(group == j)])
      pij <- Lij / (Si*Sj)
      if (pij < 1 & pij >0){
        logLik <- logLik + (Lij * log10(pij) + (Si*Sj-Lij) * log10(1-pij))
      }
    }
  }
  AIC <- 2 * k^2 + 2*S -2 * logLik
  return(c(logLik,AIC))
}

class <- as.numeric(as.factor(EuroTaxo$Class[match(rownames(A), EuroTaxo$Species)]))
order <- as.numeric(as.factor(EuroTaxo$Order[match(rownames(A), EuroTaxo$Species)]))
family <- as.numeric(as.factor(EuroTaxo$Family[match(rownames(A), EuroTaxo$Species)]))
genus <- as.numeric(as.factor(EuroTaxo$Genus[match(rownames(A), EuroTaxo$Species)]))
