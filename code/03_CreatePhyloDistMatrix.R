# Step 03: Create matrix of phylogenetic distance
# 1. Clean turtle labels
# 2. Calculate distance from 100 posterior samples
# 3. Take the mean and sd

library(purrr)
library(ape)
library(rgbif)

# Clean tips label --------------------------------------------------------
# load phylogenies
mammaltrees <- read.nexus("data/raw/phylogeny/mammals/output.nex")
squamatetrees <- read.nexus("data/raw/phylogeny/squamates/output.nex")
birdtrees <- read.nexus("data/raw/phylogeny/birds/output.nex")
amphibiantrees <- read.nexus("data/raw/phylogeny/amphibians/output.nex")
turtletrees <- read.nexus("data/raw/phylogeny/turtles/sampleFromPosterior.100.tre")
# get traits
traits <- read.csv("data/cleaned/SpeciesTraitsFull.csv", row.names = 1, stringsAsFactors = T)

# Correct tip labels
label_key <- read.csv("../phylogeny/vertlife_gbif_key.csv", row.names = 1) ##Cannot find in repo
turtle_lab <- turtletrees[[1]]$tip.label
turtle_lab <- sub("_[^_]+$", "", turtle_lab)
turtle_lab <- gsub("_", " ", turtle_lab)
turtle_labkey <- data.frame(gbif = map_df(turtle_lab, name_backbone, order = "Testudines")$species,
                            phylo = turtle_lab)

# Calculate distance for each of the 100 posterior tree -------------------
phydist <- list()

for (tree in c(1:100)){ #note there are 100 elements in each
  # mammals
  mam_tree <- mammaltrees[[tree]]
  mam_lab <- mam_tree$tip.label
  mam_lab <- gsub("_", " ", mam_lab)
  mam_tree$tip.label <- label_key$gbif[match(mam_lab, label_key$vertlife)]

  # squamates
  squam_tree <- squamatetrees[[tree]]
  squamate_lab <- squam_tree$tip.label
  squamate_lab <- gsub("_", " ", squamate_lab)
  squam_tree$tip.label <- label_key$gbif[match(squamate_lab, label_key$vertlife)]

  # birds
  bird_tree <- birdtrees[[tree]]
  bird_lab <- bird_tree$tip.label
  bird_lab <- gsub("_", " ", bird_lab)
  bird_tree$tip.label <- label_key$gbif[match(bird_lab, label_key$vertlife)]

  # amphibians
  amphi_tree <- amphibiantrees[[tree]]
  amphibian_lab <- amphi_tree$tip.label
  amphibian_lab <- gsub("_", " ", amphibian_lab)
  amphi_tree$tip.label <- label_key$gbif[match(amphibian_lab, label_key$vertlife)]

  # turtles
  turtle_tree <- turtletrees[[tree]]
  turtle_lab <- turtle_tree$tip.label
  turtle_lab <- sub("_[^_]+$", "", turtle_lab)
  turtle_lab <- gsub("_", " ", turtle_lab)
  turtle_tree <- keep.tip(turtle_tree, which(turtle_lab %in% traits$Species &!is.na(turtle_lab)))
  turtle_lab <- turtle_lab[which(turtle_lab %in% traits$Species &!is.na(turtle_lab))]
  turtle_tree$tip.label <- turtle_lab

  # combine everything
  tetrapodtree <- mam_tree + squam_tree + bird_tree + amphi_tree + turtle_tree
  
  sp <- as.character(traits$Species[traits$Species %in% tetrapodtree$tip.label])
  
  # calculate distance matrix
  dist_matrix <- cophenetic(tetrapodtree)
  phydist[[tree]] <- dist_matrix[sp, sp]
}

# Calculate mean and sd of distances over all 100 samples -----------------
phydist_mean <- apply(simplify2array(phydist), 1:2, mean)
phydist_sd <- apply(simplify2array(phydist), 1:2, sd)

write.csv(phydist_mean, "data/checkpoints/phylodist_mean.csv")
write.csv(phydist_sd, "data/checkpoints/phylodist_sd.csv")

