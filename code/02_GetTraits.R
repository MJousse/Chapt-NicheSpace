rm(list = ls())

library(dplyr)
library(rgbif)
library(purrr)
library(missForest)
source("code/functions.R")

# Add functional traits for the species in each FW ------------------------

### Get global trait database - Pretty long

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
  select(-Order, -Family, -Genus, -Best_guess_binomial) %>%
  right_join(Euro_species, by = "Species") %>%
  filter(!is.na(Species))

# Take the mode of non-numeric columns
Euro_traits_nonnum <- Euro_traits %>%
  select_at(vars(Species, Trophic_level:Other.Unknown, Class:Genus)) %>%
  group_by(Species) %>%
  summarise_all(getmode) %>%
  mutate_at(vars(-Species), as.factor)

# Take the mean of numeric columns
Euro_traits_num <- Euro_traits %>%
  select_at(vars(Species, Body_length_mm:Litter_clutch_size, Adult_svl_cm:Generation_length_d)) %>%
  group_by(Species) %>%
  summarise_all(mean, na.rm = T) %>%
  mutate_at(vars(-Species), as.numeric)

Euro_traits <- full_join(Euro_traits_nonnum, Euro_traits_num)

### Serengeti FW
Serengeti_species <- read.csv("data/cleaned/SerengetiFWTaxo.csv", row.names = 1)
Serengeti_traits <- Tetrapod_traits %>%
  select(-Order, -Family, -Genus, -Best_guess_binomial) %>%
  right_join(Serengeti_species, by = "Species") %>%
  filter(!is.na(Species))

# Take the mode of non-numeric columns
Serengeti_traits_nonnum <- Serengeti_traits %>%
  select_at(vars(Species, Trophic_level:Other.Unknown, Class:Genus)) %>%
  group_by(Species) %>%
  summarise_all(getmode) %>%
  mutate_at(vars(-Species), as.factor)

# Take the mean of numeric columns
Serengeti_traits_num <- Serengeti_traits %>%
  select_at(vars(Species, Body_length_mm:Litter_clutch_size, Adult_svl_cm:Generation_length_d)) %>%
  group_by(Species) %>%
  summarise_all(mean, na.rm = T) %>%
  mutate_at(vars(-Species), as.numeric)

Serengeti_traits <- full_join(Serengeti_traits_nonnum, Serengeti_traits_num)

### Pyrennees FW
Pyrennees_species <- read.csv("data/cleaned/pyrenneesFWTaxo.csv", row.names = 1)
Pyrennees_traits <- Tetrapod_traits %>%
  select(-Order, -Family, -Genus, -Best_guess_binomial) %>%
  right_join(Pyrennees_species, by = "Species") %>%
  filter(!is.na(Species))

# Take the mode of non-numeric columns
Pyrennees_traits_nonnum <- Pyrennees_traits %>%
  select_at(vars(Species, Trophic_level:Other.Unknown, Class:Genus)) %>%
  group_by(Species) %>%
  summarise_all(getmode) %>%
  mutate_at(vars(-Species), as.factor)

# Take the mean of numeric columns
Pyrennees_traits_num <- Pyrennees_traits %>%
  select_at(vars(Species, Body_length_mm:Litter_clutch_size, Adult_svl_cm:Generation_length_d)) %>%
  group_by(Species) %>%
  summarise_all(mean, na.rm = T) %>%
  mutate_at(vars(-Species), as.numeric)

Pyrennees_traits <- full_join(Pyrennees_traits_nonnum, Pyrennees_traits_num)

### High Arctic FW
HighArctic_species <- read.csv("data/cleaned/HighArcticFWTaxo.csv", row.names = 1)
HighArctic_traits <- Tetrapod_traits %>%
  select(-Order, -Family, -Genus, -Best_guess_binomial) %>%
  right_join(HighArctic_species, by = "Species") %>%
  filter(!is.na(Species))

# Take the mode of non-numeric columns
HighArctic_traits_nonnum <- HighArctic_traits %>%
  select_at(vars(Species, Trophic_level:Other.Unknown, Class:Genus)) %>%
  group_by(Species) %>%
  summarise_all(getmode) %>%
  mutate_at(vars(-Species), as.factor)

# Take the mean of numeric columns
HighArctic_traits_num <- HighArctic_traits %>%
  select_at(vars(Species, Body_length_mm:Litter_clutch_size, Adult_svl_cm:Generation_length_d)) %>%
  group_by(Species) %>%
  summarise_all(mean, na.rm = T) %>%
  mutate_at(vars(-Species), as.numeric)

HighArctic_traits <- full_join(HighArctic_traits_nonnum, HighArctic_traits_num)

# Checkpoint
write.csv(Euro_traits, "data/checkpoints/EuroTraits.csv")
write.csv(Serengeti_traits, "data/checkpoints/SerengetiTraits.csv")
write.csv(Pyrennees_traits, "data/checkpoints/PyrenneesTraits.csv")
write.csv(HighArctic_traits, "data/checkpoints/HighArcticTraits.csv")

Euro_traits <- read.csv("data/checkpoints/EuroTraits.csv", row.names = 1)
Serengeti_traits <- read.csv("data/checkpoints/SerengetiTraits.csv", row.names = 1)
Pyrennees_traits <- read.csv("data/checkpoints/PyrenneesTraits.csv", row.names = 1)
HighArctic_traits <- read.csv("data/checkpoints/HighArcticTraits.csv", row.names = 1)


### Impute missing data for each class separated
Tetrapods_traits <- bind_rows(Euro_traits, Serengeti_traits, Pyrennees_traits, HighArctic_traits) %>%
  distinct() %>%
  mutate_at(vars(Trophic_level:Diel_activity, Artificial_habitat_use:Genus), as.factor)

# amphibians
Amphi_traits <- filter(Tetrapods_traits, Class == "Amphibia") %>%
  droplevels()
sp <- Amphi_traits$Species
traits <- Amphi_traits %>% select(-Species)

Amphi_traits_full <- missForest(traits)
Amphi_traits_full <- data.frame(Species=sp, Amphi_traits_full$ximp) %>%
  left_join(select(Amphi_traits, Species, Genus, Family))

# birds
Bird_traits <- filter(Tetrapods_traits, Class == "Aves") %>%
  droplevels()
sp <- Bird_traits$Species
traits <- Bird_traits %>% select(-Species, -Genus, -Family,- Habitat_breadth_IUCN)

Bird_traits_full <- missForest(traits)
Bird_traits_full <- data.frame(Species=sp, Bird_traits_full$ximp) %>%
  left_join(select(Bird_traits, Species, Genus, Family))

# mammals
Mam_traits <- filter(Tetrapods_traits, Class == "Mammalia") %>%
  droplevels()
sp <- Mam_traits$Species
traits <- Mam_traits %>% select(-Species, -Genus)

Mam_traits_full <- missForest(traits)
Mam_traits_full <- data.frame(Species=sp, Mam_traits_full$ximp) %>%
  left_join(select(Mam_traits, Species, Genus, Family))

# reptiles
Rept_traits <- filter(Tetrapods_traits, Class == "Reptilia") %>%
  droplevels()
sp <- Rept_traits$Species
traits <- Rept_traits %>% select(-Species, -Genus)

Rept_traits_full <- missForest(traits)
Rept_traits_full <- data.frame(Species=sp, Rept_traits_full$ximp) %>%
  left_join(select(Rept_traits, Species, Genus, Family))

# put everything together
Tetrapods_traits <- bind_rows(Amphi_traits_full, Bird_traits_full, Mam_traits_full, Rept_traits_full) %>%
  select(-Adult_svl_cm, -Longevity_d, -Generation_length_d, -Body_length_mm)

write.csv(Tetrapods_traits, "data/cleaned/SpeciesTraitsFull.csv")
