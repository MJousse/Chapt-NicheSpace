library(igraph)
library(dplyr)
library(tidyr)
library(NetIndices)
library(multiweb)
source("code/functions.R")
source("code/functions_motifs.R")

arcticFW <- read.csv("data/cleaned/HighArcticFW.csv", row.names = 1)
arcticFW <- arcticFW %>%
  transmute(resource = Prey, consumer = Predator)
arcticRoles <- species_role(arcticFW, ncores = 8)

pyreneesFW <- read.csv("data/cleaned/pyrenneesFW.csv", row.names = 1)
pyreneesFW <- pyreneesFW %>%
  transmute(resource = Prey, consumer = Predator)
pyreneesRoles <- species_role(pyreneesFW, ncores = 8)

serengetiFW <- read.csv("data/cleaned/SerengetiFW.csv", row.names = 1)
serengetiFW <- serengetiFW %>%
  transmute(resource = Resource_Species, consumer = Consumer_Species) %>%
  drop_na()
serengetiRoles <- species_role(serengetiFW, ncores = 8)

EuropeMW <- read.csv("data/cleaned/EuroFWadults.csv", row.names = 1)
EuropeMW <- EuropeMW %>%
  transmute(resource = Prey, consumer = Predator)
EuropeRoles <- species_role(EuropeMW, ncores = 8)

save(arcticRoles, pyreneesRoles, serengetiRoles, EuropeRoles, file = "data/checkpoints/SpeciesRole.RData")
