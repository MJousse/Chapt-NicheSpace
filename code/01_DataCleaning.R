rm(list = ls())

library(rgbif)
library(purrr)
library(dplyr)
library(tidyr)
library(igraph)
source("code/functions.R")

# Standardize Food Web ----------------------------------------------------

### European Metaweb
EuroMW <- read.csv("data/raw/FW/TetraEU_pairwise_interactions.csv", sep = ";")

gbif_names <- map_df(
  unique(c(EuroMW$sourceTaxonName, EuroMW$targetTaxonName)),
  name_backbone, phylum = "Chordata") # Clean species names

EuroMW <- data.frame(Predator = gbif_names$species[match(EuroMW$sourceTaxonName, gbif_names$canonicalName)], 
                           Prey = gbif_names$species[match(EuroMW$targetTaxonName, gbif_names$canonicalName)]) %>%
  distinct() %>%
  filter(!is.na(Predator), !is.na(Prey))

write.csv(EuroMW, "data/cleaned/EuroFW.csv")

# Save clean species name with taxonomy
EuroMw_species <- gbif_names %>%
  select(Species = species, Class = class, Order = order, Family = family, Genus = genus)

write.csv(EuroMw_species, "data/cleaned/EuroMWTaxo.csv")

### Serengeti Food Web
SerengetiNodes<- read.csv("data/raw/FW/Serengeti_nodes.csv")
# Remove useless columns
SerengetiNodes <- SerengetiNodes %>% select(Taxa = Taxa..species..family.or.order., Node = Node.no.)
SerengetiNodes_clean <- data.frame()

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

# Remove duplicates
SerengetiNodes_clean <- SerengetiNodes_clean[!duplicated(SerengetiNodes_clean$Species),]

# clean food web
SerengetiInteractions <- read.csv("data/raw/FW/Serengeti_interactions.csv")
SerengetiInteractions <- SerengetiInteractions %>% filter(Consumer %in% SerengetiNodes_clean$Node & Resource %in% SerengetiNodes_clean$Node)
SerengetiInteractions_clean <- data.frame()
for (i in c(1:nrow(SerengetiInteractions))){
  if(all(SerengetiInteractions[i,] %in% SerengetiNodes_clean$Node)){
    Consumer <- SerengetiNodes_clean[SerengetiNodes_clean$Node == SerengetiInteractions$Consumer[i], c("Class", "Order", "Family", "Genus", "Species")] %>% 
      rename_all(function(x) paste0("Consumer_",x))
    Resource <- SerengetiNodes_clean[SerengetiNodes_clean$Node == SerengetiInteractions$Resource[i], c("Class", "Order", "Family", "Genus", "Species")] %>% 
      rename_all(function(x) paste0("Resource_",x))
    SerengetiInteractions_clean <- rbind(SerengetiInteractions_clean,
                                tidyr::crossing(Consumer, Resource))
  }
}
write.csv(SerengetiInteractions_clean, "data/cleaned/SerengetiFW.csv")

# Save clean species name with taxonomy
SerengetiFW_species <- SerengetiNodes_clean %>%
  select(Species, Class, Order, Family, Genus)

write.csv(SerengetiFW_species, "data/cleaned/SerengetiFWTaxo.csv")

### Pyrenees Food Web
pyrenneesFW <- read.graph("data/raw/FW/pyrenees-network.graphml", format = "graphml")
pyrenneesFW <- as.data.frame(as_edgelist(pyrenneesFW))
colnames(pyrenneesFW) <- c("Prey", "Predator")

gbif_names <- map_df(
  unique(c(pyrenneesFW$Prey, pyrenneesFW$Predator)),
  name_backbone, phylum = "Chordata") # Clean species names

pyrenneesFW <- data.frame(Predator = gbif_names$species[match(pyrenneesFW$Predator, gbif_names$canonicalName)], 
                     Prey = gbif_names$species[match(pyrenneesFW$Prey, gbif_names$canonicalName)]) %>%
  distinct() %>%
  filter(!is.na(Predator), !is.na(Prey))

write.csv(pyrenneesFW, "data/cleaned/pyrenneesFW.csv")

# Save clean species name with taxonomy
PyrenneesFW_species <- gbif_names %>%
  select(Species = species, Class = class, Order = order, Family = family, Genus = genus) %>%
  filter(!is.na(Class))

write.csv(PyrenneesFW_species, "data/cleaned/pyrenneesFWTaxo.csv")

### High Arctic Food Web
arcticFW_peak <- read.csv("data/raw/FW/HighArctic_peakyear.csv")
arcticFW_crash <- read.csv("data/raw/FW/HighArctic_crashyear.csv")
threshold <- 0

# Combine crash and peak years
arcticFW <- rbind(arcticFW_crash, arcticFW_peak) %>%
  group_by(Consumer, Resource) %>%
  summarise(Diet = max(Diet)) %>%
  filter(Diet > threshold) %>%
  select(-Diet) %>% ungroup()

# Clean nodes
arctic_nodes <- read.csv("data/raw/FW/HighArctic_nodes.csv")
gbif_names <- map_df(arctic_nodes$Species,
  name_backbone, phylum = "Chordata") # Clean species names
arctic_nodes$Species <- gbif_names$species

arcticFW_clean <- left_join(arcticFW, arctic_nodes, by = c("Consumer" = "Functiongroup")) %>%
  left_join(arctic_nodes, by = c("Resource" = "Functiongroup"), suffix = c("Consumer", "Resource")) %>%
  select(Predator = SpeciesConsumer, Prey = SpeciesResource)

write.csv(arcticFW_clean, "data/cleaned/HighArcticFW.csv")

# Save clean species name with taxonomy
arcticFW_species <- gbif_names %>%
  select(Species = species, Class = class, Order = order, Family = family, Genus = genus) %>%
  filter(!is.na(Class))

write.csv(arcticFW_species, "data/cleaned/HighArcticFWTaxo.csv")
