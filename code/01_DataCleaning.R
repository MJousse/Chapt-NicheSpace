# Step 01: Clean the food web data
# For each of the three food webs: 
# 1. extract the raw data
# 2. get gbif name of each species
# 3. standardize format into a dataframe with a column for predator and prey

rm(list = ls())
library(rgbif)
library(purrr)
library(dplyr)
library(tidyr)
library(igraph)
source("./functions.R")

# European metaweb --------------------------------------------------------
# raw data
EuroMW <- read.csv("../data/raw/FW/TetraEU_pairwise_interactions.csv", sep = ";") #tetraeu web

# # get gbif names ##long on computer
# gbif_names <- map_df(
#   unique(c(EuroMW$sourceTaxonName, EuroMW$targetTaxonName)),
#   name_backbone, phylum = "Chordata") # Clean species names ##note map_df is supserseded

gbif_names <- list_rbind(lapply(unique(c(EuroMW$sourceTaxonName, EuroMW$targetTaxonName)), name_backbone, phylum = "Chordata"))


# create different version of the fw
EuroMWadults <- filter(EuroMW, targetLifestageName == "adults")
EuroMWjuveniles <- filter(EuroMW, targetLifestageName == "larvae or young")
EuroMWeggs <- filter(EuroMW, targetLifestageName == "eggs")

EuroMWadults <- data.frame(Predator = gbif_names$species[match(EuroMWadults$sourceTaxonName, gbif_names$canonicalName)], 
                           Prey = gbif_names$species[match(EuroMWadults$targetTaxonName, gbif_names$canonicalName)]) %>%
  distinct() %>%
  filter(!is.na(Predator), !is.na(Prey))

EuroMWjuveniles <- data.frame(Predator = gbif_names$species[match(EuroMWjuveniles$sourceTaxonName, gbif_names$canonicalName)], 
                           Prey = gbif_names$species[match(EuroMWjuveniles$targetTaxonName, gbif_names$canonicalName)]) %>%
  distinct() %>%
  filter(!is.na(Predator), !is.na(Prey))

EuroMWeggs <- data.frame(Predator = gbif_names$species[match(EuroMWeggs$sourceTaxonName, gbif_names$canonicalName)], 
                           Prey = gbif_names$species[match(EuroMWeggs$targetTaxonName, gbif_names$canonicalName)]) %>%
  distinct() %>%
  filter(!is.na(Predator), !is.na(Prey))

# save cleaned fw
write.csv(EuroMWadults, "data/cleaned/EuroFWadults.csv")
write.csv(EuroMWjuveniles, "data/cleaned/EuroFWjuveniles.csv")
write.csv(EuroMWeggs, "data/cleaned/EuroFWeggs.csv")

# save species name with taxonomy
EuroMw_species <- gbif_names %>%
  select(Species = species, Class = class, Order = order, Family = family, Genus = genus)
write.csv(EuroMw_species, "data/cleaned/EuroMWTaxo.csv")


# Serengeti food web ------------------------------------------------------
# raw data
SerengetiNodes<- read.csv("../data/raw/FW/Serengeti_nodes.csv")
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
SerengetiInteractions <- read.csv("../data/raw/FW/Serengeti_interactions.csv")
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

# save cleaned food web
write.csv(SerengetiInteractions_clean, "data/cleaned/SerengetiFW.csv")

# save species name with taxonomy
SerengetiFW_species <- SerengetiNodes_clean %>%
  select(Species, Class, Order, Family, Genus)
write.csv(SerengetiFW_species, "data/cleaned/SerengetiFWTaxo.csv")

# Pyrenees food web -------------------------------------------------------
# raw data
#pyrenneesFW <- read.graph("../data/raw/FW/pyrenees-network.graphml", format = "graphml") #read.graph deprecated
pyrenneesFW <- read_graph("../data/raw/FW/pyrenees-network.graphml", format = "graphml") 

pyrenneesFW <- as.data.frame(as_edgelist(pyrenneesFW))
colnames(pyrenneesFW) <- c("Prey", "Predator")

# # get gbif species name
# gbif_names <- map_df(
#   unique(c(pyrenneesFW$Prey, pyrenneesFW$Predator)),
#   name_backbone, phylum = "Chordata") %>%
#   filter(class %in% c("Amphibia", "Reptilia", "Mammalia", "Aves")) # map_df deprecated

gbif_names <- list_rbind(lapply(unique(c(pyrenneesFW$Prey, pyrenneesFW$Predator)), name_backbone, phylum = "Chordata")) %>%
  filter(class %in% c("Amphibia", "Reptilia", "Mammalia", "Aves"))


# standardize food web
pyrenneesFW <- data.frame(Predator = gbif_names$species[match(pyrenneesFW$Predator, gbif_names$canonicalName)], 
                     Prey = gbif_names$species[match(pyrenneesFW$Prey, gbif_names$canonicalName)]) %>%
  distinct() %>%
  filter(!is.na(Predator), !is.na(Prey))

# save cleaned food web
write.csv(pyrenneesFW, "data/cleaned/pyrenneesFW.csv")

# save species name with taxonomy
PyrenneesFW_species <- gbif_names %>%
  select(Species = species, Class = class, Order = order, Family = family, Genus = genus) %>%
  filter(!is.na(Class) & Class %in% c("Amphibia", "Reptilia", "Aves", "Mammalia"))
write.csv(PyrenneesFW_species, "data/cleaned/pyrenneesFWTaxo.csv")

# High Arctic food web ----------------------------------------------------
# raw data
arcticFW <- read.csv("../data/raw/FW/Tundra_Nunavik_Trophic_relationships.txt", sep = "\t", row.names = 1)

# get gbif species name
arctic_species <- rownames(arcticFW)
# gbif_names <- map_df(arctic_species,
#   name_backbone, phylum = "Chordata") # Clean species names #map_df deprecated
gbif_names <- list_rbind(lapply(arctic_species, name_backbone, phylum = "Chordata"))

arctic_vertebrates_i <- gbif_names$class %in% c("Aves", "Mammalia", "Reptilia", "Amphibia") & !is.na(gbif_names$species)
arctic_vertebrates <- gbif_names[arctic_vertebrates_i,]

# keep only vertebrates
arcticFW <- arcticFW[arctic_vertebrates_i, arctic_vertebrates_i]

# standardize food web format
arctic_interactions <- which(arcticFW == 1, arr.ind = T, useNames = F)
arcticFW_clean <- data.frame(Predator = arctic_vertebrates$species[arctic_interactions[,2]],
                             Prey = arctic_vertebrates$species[arctic_interactions[,1]])

# save cleaned food web
write.csv(arcticFW_clean, "data/cleaned/HighArcticFW.csv")

# save clean species name with taxonomy
arcticFW_species <- arctic_vertebrates %>%
  dplyr::select(Species = species, Class = class, Order = order, Family = family, Genus = genus) %>%
  filter(!is.na(Class))
write.csv(arcticFW_species, "data/cleaned/HighArcticFWTaxo.csv")
