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
traits <- read.csv("data/cleaned/SpeciesTraitsFull.csv", row.names = 1, stringsAsFactors = T)

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
turtletree <- keep.tip(turtletree, which(turtle_lab %in% traits$Species &!is.na(turtle_lab)))
turtletree$tip.label <- turtle_lab[turtle_lab %in% traits$Species &!is.na(turtle_lab)]

# create big phylogeny
tetrapodtrees <- mammaltrees[[sample]] + squamatetrees[[sample]] + birdtrees[[sample]] + amphibiantrees[[sample]] + turtletree

# species list in FW
europe <- read.csv("data/cleaned/EuroMWTaxo.csv", row.names = 1)
pyrenees <- read.csv("data/cleaned/pyrenneesFWTaxo.csv", row.names = 1)
serengeti <- read.csv("data/cleaned/SerengetiFWTaxo.csv", row.names = 1) %>% drop_na()
arctic <- read.csv("data/cleaned/HighArcticFWTaxo.csv", row.names = 1)

# Phylogenetic distance matrix
phydist <- cophenetic(tetrapodtrees)
write.csv(phydist, "data/checkpoints/phylodist.csv")

mntd <- function(target_species, species_pool, phylodist){ # Shortest phylo distance
  if (target_species %in% colnames(phylodist)){
    min(phylodist[target_species, colnames(phylodist) %in% species_pool])
  } else {NA}
}

arctic_mntd <- map_dbl(arctic$Species, mntd, europe$Species, phydist)
pyrenees_mntd <- map_dbl(pyrenees$Species, mntd, europe$Species, phydist)
serengeti_mntd <- map_dbl(serengeti$Species, mntd, europe$Species, phydist)

# functional distance 
library(mFD)
rownames(traits) <- traits$Species
traits <- traits %>%
  dplyr::select(-Species, -Class, -Order, -Family, -Genus) %>%
  drop_na()

traits_cat <- data.frame(trait_name = colnames(traits),
                         trait_type = c("N", "Q", rep("F", 12), "Q", "Q", "Q", "F", "F", "F"),
                         fuzzy_name = c(NA, NA, rep("Habitat", 12), NA, NA, NA, rep("TrophicLevel", 3)))

sp_funct.dist <- funct.dist(traits, traits_cat, metric = "gower", scale_euclid = "scale_center")

fnnd <- function(target_species, species_pool, funcdist){
  funcdist <- as.matrix(funcdist)
  if (target_species %in% colnames(funcdist)){
    min(funcdist[target_species, colnames(funcdist) %in% species_pool])
  } else {NA}
}

fmpd <- function(target_species, species_pool, funcdist){
  funcdist <- as.matrix(funcdist)
  if (target_species %in% colnames(funcdist)){
    mean(funcdist[target_species, colnames(funcdist) %in% species_pool])
  } else {NA}
}

arctic_fnnd <- map_dbl(arctic$Species, fnnd, europe$Species, sp_funct.dist)
pyrenees_fnnd <- map_dbl(pyrenees$Species, fnnd, europe$Species, sp_funct.dist)
serengeti_fnnd <- map_dbl(serengeti$Species, fnnd, europe$Species, sp_funct.dist)

arctic_fmpd <- map_dbl(arctic$Species, fmpd, europe$Species, sp_funct.dist)
pyrenees_fmpd <- map_dbl(pyrenees$Species, fmpd, europe$Species, sp_funct.dist)
serengeti_fmpd <- map_dbl(serengeti$Species, fmpd, europe$Species, sp_funct.dist)

df <- data.frame(Arctic = arctic_fnnd, Pyrenees = pyrenees_fnnd, Serengeti = serengeti_fnnd)
fills <- c("Serengeti" = "yellowgreen", "Pyrenees" = "red3", "Arctic" = "deepskyblue")
p1 <- ggplot() +
  geom_histogram(data = data.frame(serengeti_fmpd), aes(serengeti_fmpd, fill = "Serengeti"), alpha= 0.3) +
  geom_histogram(data = data.frame(pyrenees_fmpd), aes(pyrenees_fmpd, fill = "Pyrenees"), alpha= 0.3) +
  geom_histogram(data = data.frame(arctic_fmpd), aes(arctic_fmpd, fill = "Arctic"), alpha= 0.3) +
  theme_minimal() +
  labs(fill = "Food web", x = "Mean functional distance") +
  scale_fill_manual(values = fills) +
  theme(legend.position = "bottom",
        legend.direction = "horizontal")
  
p2 <- ggplot() +
  geom_histogram(data = data.frame(serengeti_mntd), aes(serengeti_mntd, fill = "Serengeti"), alpha= 0.3) +
  geom_histogram(data = data.frame(pyrenees_mntd), aes(pyrenees_mntd, fill = "Pyrenees"), alpha= 0.3) +
  geom_histogram(data = data.frame(arctic_mntd), aes(arctic_mntd, fill = "Arctic"), alpha= 0.3) +
  theme_minimal() +
  labs(fill = "Food web", x = "Phylogenetic distance to neareast taxon") +
  scale_fill_manual(values = fills) +
  theme(legend.position = "bottom",
        legend.direction = "horizontal")

library(patchwork)
p<-p1+p2

ggsave("figures/SI/SPdist.png", p, device = "png")
  
