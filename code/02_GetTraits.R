rm(list = ls())

library(dplyr)

# Add functional traits for the species in each FW ------------------------

### Get global trait database 

# Etard et al. 2020 (https://doi.org/10.1111/geb.13184)
Amphi_traits <- read.csv("data/raw/traits/traits_etard2020/Amphibians.csv")
Amphi_traits$Species <- map_df(Amphi_traits$Best_guess_binomial, name_backbone, class = "Amphibia")$species
Bird_traits <- read.csv("data/raw/traits/traits_etard2020/Birds.csv")
Bird_traits$Species <- map_df(Bird_traits$Best_guess_binomial, name_backbone, class = "Aves")$species
Mam_traits <- read.csv("data/raw/traits/traits_etard2020/Mammals.csv")
Mam_traits$Species <- map_df(Mam_traits$Best_guess_binomial, name_backbone, class = "Mammalia")$species
Rept_traits <- read.csv("data/raw/traits/traits_etard2020/Reptiles.csv")
Rept_traits$Species <- map_df(Rept_traits$Best_guess_binomial, name_backbone, class = "Reptilia")$species

Tetrapod_traits <- bind_rows(Amphi_traits, Bird_traits, Mam_traits, Rept_traits)

write.csv(Tetrapod_traits, "data/checkpoints/TetrapodTraits.csv")

### European metaweb
Euro_species <- read.csv("data/cleaned/EuroMWTaxo.csv", row.names = 1)
Euro_traits <- Tetrapod_traits %>%
  select(-Order, - Family, -Genus) %>%
  right_join(Euro_species, by = c("Best_guess_binomial" = "Species"))

### Serengeti FW
Serengeti_species <- read.csv("data/cleaned/SerengetiFWTaxo.csv", row.names = 1)
Serengeti_traits <- Tetrapod_traits %>%
  select(-Order, - Family, -Genus) %>%
  right_join(Serengeti_species, by = c("Best_guess_binomial" = "Species"))

### Pyrennees FW
Pyrennees_species <- read.csv("data/cleaned/pyrenneesFWTaxo.csv", row.names = 1)
Pyrennees_traits <- Tetrapod_traits %>%
  select(-Order, - Family, -Genus) %>%
  right_join(Pyrennees_species, by = c("Best_guess_binomial" = "Species"))

### High Arctic FW
HighArctic_species <- read.csv("data/cleaned/HighArcticFWTaxo.csv", row.names = 1)
HighArctic_traits <- Tetrapod_traits %>%
  select(-Order, - Family, -Genus) %>%
  right_join(HighArctic_species, by = c("Best_guess_binomial" = "Species"))

