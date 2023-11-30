# SI threat Serengeti node's as the mean across species:
library(dplyr)
library(rgbif)
# Serengeti food web ------------------------------------------------------
# raw data
SerengetiNodes<- read.csv("data/raw/FW/Serengeti_nodes.csv")
# Remove useless columns
SerengetiNodes <- SerengetiNodes %>% select(Taxa = Taxa..species..family.or.order., Node = Node.no.)
SerengetiNodes_clean <- data.frame()
# get gbif and node id for each species
for (i in c(1:nrow(SerengetiNodes))){
  if (SerengetiNodes$Taxa[i] != "-"){
    taxa <- name_backbone(SerengetiNodes$Taxa[i])
    if (!is.null(taxa$class) && taxa$class %in% c("Mammalia", "Reptilia", "Amphibia", "Aves")) {
      Class <- ifelse(is.null(taxa$class), NA, taxa$class)
      Order <- ifelse(is.null(taxa$order), NA, taxa$order)
      Family <- ifelse(is.null(taxa$family), NA, taxa$family)
      Genus <- ifelse(is.null(taxa$genus), NA, taxa$genus)
      Species <- ifelse(is.null(taxa$species), NA, taxa$species)
      if(SerengetiNodes$Node[i] != ""){
        current_node = SerengetiNodes$Node[i]
      }
      SerengetiNodes_clean <- rbind(
        SerengetiNodes_clean,
        data.frame(Original = SerengetiNodes$Taxa[i], Class, Order, Family, Genus, Species, Node = current_node)
      )
    }
  }
}

# Some manual correction
SerengetiNodes_clean[SerengetiNodes_clean$Original == "Tadarida pumila", c("Genus", "Species")] <- c("Chaerephon", "Chaerephon pumilus")
SerengetiNodes_clean[SerengetiNodes_clean$Original == "Nycticeius hirundo", c("Genus", "Species")] <- c("Scotoecus", "Scotoecus hirundo")
SerengetiNodes_clean[SerengetiNodes_clean$Original == "Tadarida africana", c("Genus", "Species")] <- c("Tadarida", "Tadarida fulminans")
SerengetiNodes_clean[SerengetiNodes_clean$Original == "Pachydactylus turneri", c("Genus", "Species")] <- c("Chondrodactylus", "Chondrodactylus turneri")
SerengetiNodes_clean[SerengetiNodes_clean$Original == "Lygosoma afrum", c("Genus", "Species")] <- c("Mochlus", "Mochlus sundevallii")
SerengetiNodes_clean[SerengetiNodes_clean$Original == "Mabuya striata", c("Genus", "Species")] <- c("Trachylepis", "Trachylepis striata")
SerengetiNodes_clean[SerengetiNodes_clean$Original == "Ma buya maculilabris", c("Genus", "Species")] <- c("Trachylepis", "Trachylepis maculilabris")
SerengetiNodes_clean[SerengetiNodes_clean$Original == "Procavia johnstoni", c("Genus", "Species")] <- c("Procavia", "Procavia capensis")

# remove duplicates
SerengetiNodes_clean <- SerengetiNodes_clean[!duplicated(SerengetiNodes_clean$Species),]

# clean food web
SerengetiInteractions <- read.csv("data/raw/FW/Serengeti_interactions.csv")
SerengetiInteractions <- SerengetiInteractions %>% transmute(Predator = as.character(Consumer), Prey = as.character(Resource)) %>%
  filter(Predator %in% SerengetiNodes_clean$Node & Prey %in% SerengetiNodes_clean$Node)

# get predictors
library(brms)
source("code/functions.R")
FuncTraits <- read.csv("data/cleaned/SpeciesTraitsFull.csv", row.names = 1)
brms_form <- bf(interaction ~ 1 + 
                  (Omnivore.predator + Carnivore.predator + Habitat_breadth.predator + BM.predator + Longevity.predator + ClutchSize.predator +
                     Omnivore.prey + Carnivore.prey + Habitat_breadth.prey + BM.prey + Longevity.prey + ClutchSize.prey + 
                     ActivityTime.match + Habitat.match + BM.match) + 
                  (1 + (Omnivore.predator + Carnivore.predator + Habitat_breadth.predator + BM.predator + Longevity.predator + ClutchSize.predator 
                        + Omnivore.prey + Carnivore.prey + Habitat_breadth.prey + BM.prey + Longevity.prey + ClutchSize.prey + 
                          ActivityTime.match + Habitat.match + BM.match) || Order.predator), 
                family = bernoulli())


# Transform into predictors and scale
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

SerengetiNodes_traits <- left_join(SerengetiNodes_clean, FuncTraits, by = "Species")

SerengetiNodes_traits <- SerengetiNodes_traits %>%
  group_by(Node) %>%
  summarise(Order = getmode(Order.x),
         Diel_activity = getmode(Diel_activity),
         Habitat_breadth_IUCN = mean(Habitat_breadth_IUCN, na.rm = T),
         Forest = getmode(Forest), Savanna = getmode(Savanna), Shrubland = getmode(Shrubland), Grassland = getmode(Grassland), Wetland = getmode(Wetland), Rocky.areas = getmode(Rocky.areas), Caves.and.subterranean = getmode(Caves.and.subterranean), Desert = getmode(Desert), Marine = getmode(Marine),
         Marine.intertidal.or.coastal.supratidal = getmode(Marine.intertidal.or.coastal.supratidal), Artificial = getmode(Artificial), Introduced.vegetation = getmode(Introduced.vegetation),
         logBM = mean(logBM, na.rm = T), logLongevity = mean(logLongevity, na.rm = T), logClutchSize = mean(logClutchSize, na.rm = T), Herbivore = getmode(Herbivore), Omnivore = getmode(Omnivore), Carnivore = getmode(Carnivore)) %>%
  filter(!is.na(Diel_activity))
SerengetiFW2 <- expand.grid(SerengetiNodes_traits$Node, SerengetiNodes_traits$Node)
colnames(SerengetiFW2) <- c("Predator", "Prey")
SerengetiFW2 <- left_join(SerengetiFW2, SerengetiNodes_traits, by = c("Prey" = "Node")) %>%
  left_join(SerengetiNodes_traits, by = c("Predator" = "Node"))
SerengetiFW2 <- traits2predictors(SerengetiFW2)

predictors <- get_predictors(FuncTraits$Species, FuncTraits)
var2scale <- c("Habitat_breadth.predator", "BM.predator", "Longevity.predator", 
               "ClutchSize.predator", "Habitat_breadth.prey", "BM.prey",
               "Longevity.prey", "ClutchSize.prey", "Habitat.match", "BM.match")
pool_mean <- apply(predictors[,var2scale], MARGIN = 2, mean)
pool_sd <- apply(predictors[,var2scale], MARGIN = 2, sd)

# scale predictors
SerengetiFW2[,var2scale] <- sweep(SerengetiFW2[,var2scale], MARGIN = 2, pool_mean)
SerengetiFW2[,var2scale] <- sweep(SerengetiFW2[,var2scale], MARGIN = 2, FUN = "/", 2*pool_sd)

# training model
# add response
SerengetiInteractions$interaction <- 1
SerengetiFW2 <- left_join(SerengetiFW2, SerengetiInteractions)
SerengetiFW2$interaction[is.na(SerengetiFW2$interaction)] <- 0

# keep 30% of the food web for validation controlling for prevalence
Npos = round(0.3 * sum(SerengetiFW2$interaction == 1))
Nneg = round(0.3 * sum(SerengetiFW2$interaction == 0))
testing_positives <- sample(seq_len(sum(SerengetiFW2$interaction)), size = Npos)
testing_negatives <- sample(seq_len(sum(SerengetiFW2$interaction==0)), size = Nneg)
testing_id_ser <- as.numeric(c(
  rownames(SerengetiFW2[SerengetiFW2$interaction ==1,][testing_positives,]),
  rownames(SerengetiFW2[SerengetiFW2$interaction ==0,][testing_negatives,])
))

# take the remaining positives and as many negatives for training
training <- SerengetiFW2[-testing_id_ser, ]
N = sum(training$interaction)
training_negatives <- sample(seq_len(sum(training$interaction==0)), size = N)
training_id_ser <- as.numeric(c(
  rownames(training[training$interaction ==1,]),
  rownames(training[training$interaction ==0,][training_negatives,])
))
training <- SerengetiFW2[training_id_ser, ]

training <- select(training,  -Predator, -Prey, -Order.prey, 
                   -Herbivore.predator,  -Herbivore.prey) %>%
  mutate(Order.predator = factor(Order.predator, 
                                 levels = unique(FuncTraits$Order)))

get_prior(brms_form, data = training)

model_priors <- c(
  prior(normal(0, 1), class = "Intercept"),
  prior(normal(0, 1), class = "b"),
  prior(cauchy(0, 5), class = "sd")
)

prior_predictions <- brm(formula = brms_form,
                         data = training,
                         prior = model_priors,
                         sample_prior = "only", init = "0")
SerengetiModel2 <- brm(formula = brms_form,
                 data = training,
                 prior = model_priors, sample_prior = "no", 
                 cores = 4, backend = "cmdstan", threads = 4,
                 init = "0", iter = 2000, refresh = 50)

# save the model on OneDrive (too big for Github...)
save(SerengetiModel2, SerengetiFW2, testing_id_ser, file = "~/OneDrive/Chapt-NicheSpace/models/SerengetiModel2.RData")
