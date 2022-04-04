# I want to calculate the phylo distance between a given species and a species pool (distance to nearest taxon)

library(ape)
library(dplyr)
library(tidyr)
library(purrr)
library(rgbif)
sample = 1

# load phylogenies
mammaltrees <- read.nexus("data/raw/phylogeny/mammals/output.nex")
squamatetrees <- read.nexus("data/raw/phylogeny/squamates/output.nex")
birdtrees <- read.nexus("data/raw/phylogeny/birds/output.nex")
amphibiantrees <- read.nexus("data/raw/phylogeny/amphibians/output.nex")
turtletrees <- read.nexus("data/raw/phylogeny/turtles/sampleFromPosterior.100.tre")

# Correct tip labels
label_key <- read.csv("../phylogeny/vertlife_gbif_key.csv", row.names = 1)
mam_lab <- mammaltrees[[sample]]$tip.label
mam_lab <- gsub("_", " ", mam_lab)
mammaltrees[[sample]]$tip.label <- label_key$gbif[match(mam_lab, label_key$vertlife)]
squamate_lab <- squamatetrees[[sample]]$tip.label
squamate_lab <- gsub("_", " ", squamate_lab)
squamatetrees[[sample]]$tip.label <- label_key$gbif[match(squamate_lab, label_key$vertlife)]
bird_lab <- birdtrees[[sample]]$tip.label
bird_lab <- gsub("_", " ", bird_lab)
birdtrees[[sample]]$tip.label <- label_key$gbif[match(bird_lab, label_key$vertlife)]
amphibian_lab <- amphibiantrees[[sample]]$tip.label
amphibian_lab <- gsub("_", " ", amphibian_lab)
amphibiantrees[[sample]]$tip.label <- label_key$gbif[match(amphibian_lab, label_key$vertlife)]
turtle_lab <- turtletrees[[sample]]$tip.label
turtle_lab <- sub("_[^_]+$", "", turtle_lab)
turtle_lab <- gsub("_", " ", turtle_lab)
turtle_lab <- map_df(turtle_lab, name_backbone, order = "Testudines")$species
turtletree <- turtletrees[[sample]]
turtletree <- keep.tip(turtletree, which(turtle_lab %in% species &!is.na(turtle_lab)))
turtletree$tip.label <- turtle_lab[turtle_lab %in% species &!is.na(turtle_lab)]

# create big phylogeny
tetrapodtrees <- mammaltrees[[sample]] + squamatetrees[[sample]] + birdtrees[[sample]] + amphibiantrees[[sample]] + turtletree

# create community matrix
europe <- read.csv("data/cleaned/EuroMWTaxo.csv", row.names = 1)
pyrenees <- read.csv("data/cleaned/pyrenneesFWTaxo.csv", row.names = 1)
serengeti <- read.csv("data/cleaned/SerengetiFWTaxo.csv", row.names = 1) %>% drop_na()
arctic <- read.csv("data/cleaned/HighArcticFWTaxo.csv", row.names = 1)
species <- unique(c(europe$Species, pyrenees$Species, serengeti$Species, arctic$Species)) 
species <- species[species %in% tetrapodtrees$tip.label]

comm <- c()%>%
  rbind(species %in% europe$Species) %>%
  rbind(species %in% pyrenees$Species) %>%
  rbind(species %in% serengeti$Species) %>%
  rbind(species %in% arctic$Species) %>%
  as.data.frame(row.names = c("Europe", "Pyrenees", "Serengeti", "Arctic"))

colnames(comm) <- species

# Phylogenetic distance matrix
phydist <- cophenetic(tetrapodtrees)

minphylodist <- function(target_species, species_pool, phylodist){
  if (target_species %in% colnames(phylodist)){
    min(phylodist[target_species, colnames(phylodist) %in% species_pool])
  } else {NA}
}

arcticdist <- map_dbl(arctic$Species, minphylodist, europe$Species, phydist)
pyreneesdist <- map_dbl(pyrenees$Species, minphylodist, europe$Species, phydist)
serengetidist <- map_dbl(serengeti$Species, minphylodist, europe$Species, phydist)
