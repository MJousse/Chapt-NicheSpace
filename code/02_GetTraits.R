# Step 02: Clean trait data
# 1. Load raw data from Etard et al. 2020 (https://doi.org/10.1111/geb.13184)
# 2. Get GBIF species name
# 3. For each food web, extract traits
# 4. Impute missing data for each class separately
# 5. Keep traits we want and log-transform some of them

rm(list = ls())
library(dplyr)
library(rgbif)
library(purrr)
library(missForest)
source("code/functions.R")

# Clean trait data --------------------------------------------------------
# load data
Amphi_traits <- read.csv("data/raw/traits/traits_etard2020/Amphibians.csv")
Bird_traits <- read.csv("data/raw/traits/traits_etard2020/Birds.csv")
Mam_traits <- read.csv("data/raw/traits/traits_etard2020/Mammals.csv")
Rept_traits <- read.csv("data/raw/traits/traits_etard2020/Reptiles.csv")

# get GBIF species name - Pretty long
Amphi_traits$Species <- map_df(Amphi_traits$Best_guess_binomial, name_backbone, class = "Amphibia")$species
Bird_traits$Species <- map_df(Bird_traits$Best_guess_binomial, name_backbone, class = "Aves")$species
Mam_traits$Species <- map_df(Mam_traits$Best_guess_binomial, name_backbone, class = "Mammalia")$species
Rept_traits$Species <- map_df(Rept_traits$Best_guess_binomial, name_backbone, class = "Reptilia")$species

# since that was long, save the results
Tetrapod_traits <- bind_rows(Amphi_traits, Bird_traits, Mam_traits, Rept_traits)
write.csv(Tetrapod_traits, "data/checkpoints/TetrapodTraits.csv")


# European metaweb --------------------------------------------------------
# load species
Euro_species <- read.csv("data/cleaned/EuroMWTaxo.csv", row.names = 1)
Euro_traits <- Tetrapod_traits %>%
  select(-Order, -Family, -Genus, -Best_guess_binomial) %>%
  right_join(Euro_species, by = "Species") %>%
  filter(!is.na(Species))

# take the mode of non-numeric columns
Euro_traits_nonnum <- Euro_traits %>%
  select_at(vars(Species, Trophic_level:Other.Unknown, Class:Genus)) %>%
  group_by(Species) %>%
  summarise_all(getmode) %>%
  mutate_at(vars(-Species), as.factor)

# take the mean of numeric columns
Euro_traits_num <- Euro_traits %>%
  select_at(vars(Species, Body_length_mm:Litter_clutch_size, Adult_svl_cm:Generation_length_d)) %>%
  group_by(Species) %>%
  summarise_all(mean, na.rm = T) %>%
  mutate_at(vars(-Species), as.numeric)

# save the results
Euro_traits <- full_join(Euro_traits_nonnum, Euro_traits_num)
write.csv(Euro_traits, "data/checkpoints/EuroTraits.csv")

# Serengeti food web ------------------------------------------------------
# load species
Serengeti_species <- read.csv("data/cleaned/SerengetiFWTaxo.csv", row.names = 1)
Serengeti_traits <- Tetrapod_traits %>%
  select(-Order, -Family, -Genus, -Best_guess_binomial) %>%
  right_join(Serengeti_species, by = "Species") %>%
  filter(!is.na(Species))

# take the mode of non-numeric columns
Serengeti_traits_nonnum <- Serengeti_traits %>%
  select_at(vars(Species, Trophic_level:Other.Unknown, Class:Genus)) %>%
  group_by(Species) %>%
  summarise_all(getmode) %>%
  mutate_at(vars(-Species), as.factor)

# take the mean of numeric columns
Serengeti_traits_num <- Serengeti_traits %>%
  select_at(vars(Species, Body_length_mm:Litter_clutch_size, Adult_svl_cm:Generation_length_d)) %>%
  group_by(Species) %>%
  summarise_all(mean, na.rm = T) %>%
  mutate_at(vars(-Species), as.numeric)

# save the results
Serengeti_traits <- full_join(Serengeti_traits_nonnum, Serengeti_traits_num)
write.csv(Serengeti_traits, "data/checkpoints/SerengetiTraits.csv")

# Pyrenees food web -------------------------------------------------------
# load species
Pyrennees_species <- read.csv("data/cleaned/pyrenneesFWTaxo.csv", row.names = 1)
Pyrennees_traits <- Tetrapod_traits %>%
  select(-Order, -Family, -Genus, -Best_guess_binomial) %>%
  right_join(Pyrennees_species, by = "Species") %>%
  filter(!is.na(Species))

# take the mode of non-numeric columns
Pyrennees_traits_nonnum <- Pyrennees_traits %>%
  select_at(vars(Species, Trophic_level:Other.Unknown, Class:Genus)) %>%
  group_by(Species) %>%
  summarise_all(getmode) %>%
  mutate_at(vars(-Species), as.factor)

# take the mean of numeric columns
Pyrennees_traits_num <- Pyrennees_traits %>%
  select_at(vars(Species, Body_length_mm:Litter_clutch_size, Adult_svl_cm:Generation_length_d)) %>%
  group_by(Species) %>%
  summarise_all(mean, na.rm = T) %>%
  mutate_at(vars(-Species), as.numeric)

# save the results
Pyrennees_traits <- full_join(Pyrennees_traits_nonnum, Pyrennees_traits_num)
write.csv(Pyrennees_traits, "data/checkpoints/PyrenneesTraits.csv")

# High Arctic food web ----------------------------------------------------
# load species
HighArctic_species <- read.csv("data/cleaned/HighArcticFWTaxo.csv", row.names = 1)
HighArctic_traits <- Tetrapod_traits %>%
  select(-Order, -Family, -Genus, -Best_guess_binomial) %>%
  right_join(HighArctic_species, by = "Species") %>%
  filter(!is.na(Species))

# take the mode of non-numeric columns
HighArctic_traits_nonnum <- HighArctic_traits %>%
  select_at(vars(Species, Trophic_level:Other.Unknown, Class:Genus)) %>%
  group_by(Species) %>%
  summarise_all(getmode) %>%
  mutate_at(vars(-Species), as.factor)

# take the mean of numeric columns
HighArctic_traits_num <- HighArctic_traits %>%
  select_at(vars(Species, Body_length_mm:Litter_clutch_size, Adult_svl_cm:Generation_length_d)) %>%
  group_by(Species) %>%
  summarise_all(mean, na.rm = T) %>%
  mutate_at(vars(-Species), as.numeric)

# save the results
HighArctic_traits <- full_join(HighArctic_traits_nonnum, HighArctic_traits_num)
write.csv(HighArctic_traits, "data/checkpoints/HighArcticTraits.csv")

# Checkpoint
Euro_traits <- read.csv("data/checkpoints/EuroTraits.csv", row.names = 1)
Serengeti_traits <- read.csv("data/checkpoints/SerengetiTraits.csv", row.names = 1)
Pyrennees_traits <- read.csv("data/checkpoints/PyrenneesTraits.csv", row.names = 1)
HighArctic_traits <- read.csv("data/checkpoints/HighArcticTraits.csv", row.names = 1)


# Check trait coverage ----------------------------------------------------
# europe
Euro_coverage <- select(Euro_traits, -Species, -Genus, -Class, -Order, -Family) %>% #remove taxonomy
  mutate(habitat = Forest) %>% # let's threat habitat as 1 trait
  select(Trophic_level, Diel_activity, habitat, Body_mass_g, Longevity_d, Litter_clutch_size)
Euro_coverage <- sum(!is.na(Euro_coverage)) / (ncol(Euro_coverage) * nrow(Euro_coverage))

# arctic
Arctic_coverage <- select(HighArctic_traits, -Species, -Genus, -Class, -Order, -Family) %>% #remove taxonomy
  mutate(habitat = Forest) %>% # let's threat habitat as 1 trait
  select(Trophic_level, Diel_activity, habitat, Body_mass_g, Longevity_d, Litter_clutch_size)
Arctic_coverage <- sum(!is.na(Arctic_coverage)) / (ncol(Arctic_coverage) * nrow(Arctic_coverage))

# pyrenees
Pyrenees_coverage <- select(Pyrennees_traits, -Species, -Genus, -Class, -Order, -Family) %>% #remove taxonomy
  mutate(habitat = Forest) %>% # let's threat habitat as 1 trait
  select(Trophic_level, Diel_activity, habitat, Body_mass_g, Longevity_d, Litter_clutch_size)
Pyrenees_coverage <- sum(!is.na(Pyrenees_coverage)) / (ncol(Pyrenees_coverage) * nrow(Pyrenees_coverage))

# serengeti
Serengeti_coverage <- select(Serengeti_traits, -Species, -Genus, -Class, -Order, -Family) %>% #remove taxonomy
  mutate(habitat = Forest) %>% # let's threat habitat as 1 trait
  select(Trophic_level, Diel_activity, habitat, Body_mass_g, Longevity_d, Litter_clutch_size)
Serengeti_coverage <- sum(!is.na(Serengeti_coverage)) / (ncol(Serengeti_coverage) * nrow(Serengeti_coverage))

# Impute missing data for each class separately ----------------------------
# Join data from all food webs
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
traits <- Bird_traits %>% select(-Species, -Genus, -Family)
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

# Put everything together -------------------------------------------------
# log transform habitat breadth, body mass, longevity, clutch size
# separate trophic levels
# keep traits of interest
Tetrapods_traits <- bind_rows(Amphi_traits_full, Bird_traits_full, Mam_traits_full, Rept_traits_full) %>%
  mutate(Habitat_breadth_IUCN = log(Habitat_breadth_IUCN),
         logBM = log(Body_mass_g),
         logLongevity = log(Max_longevity_d),
         logClutchSize = log(Litter_clutch_size),
         Herbivore = ifelse(Trophic_level == "Herbivore", 1,0),
         Omnivore = ifelse(Trophic_level == "Omnivore", 1,0),
         Carnivore = ifelse(Trophic_level == "Carnivore", 1, 0)) %>%
  select(-Adult_svl_cm, -Longevity_d, -Generation_length_d, -Body_length_mm, -Maturity_d, -Artificial_habitat_use, -Other.Unknown, -Body_mass_g,
         -Trophic_level, -Max_longevity_d, -Litter_clutch_size)

# save final 
write.csv(Tetrapods_traits, "data/cleaned/SpeciesTraitsFull.csv")