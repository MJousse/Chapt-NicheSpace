rm(list = ls())

library(dplyr)
library(rgbif)
library(purrr)
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

write.csv(Euro_traits, "data/cleaned/EuroTraits.csv")
write.csv(Serengeti_traits, "data/cleaned/SerengetiTraits.csv")
write.csv(Pyrennees_traits, "data/cleaned/PyrenneesTraits.csv")
write.csv(HighArctic_traits, "data/cleaned/HighArcticTraits.csv")
