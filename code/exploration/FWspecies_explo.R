rm(list = ls())
library(dplyr)
library(tidyr)
library(ggplot2)


EuroTraits <- read.csv("data/cleaned/EuroTraits.csv", row.names = 1)
SerengetiTraits <- read.csv("data/cleaned/SerengetiTraits.csv", row.names = 1)
PyrenneesTraits <- read.csv("data/cleaned/PyrenneesTraits.csv", row.names = 1)
HighArcticTraits <- read.csv("data/cleaned/HighArcticTraits.csv", row.names = 1)

##### Taxonomic overlap
# Pyrennees:
pyr <- data.frame(FW = rep("Pyrennees", nrow(PyrenneesTraits)))
pyr$species = PyrenneesTraits$Species %in% EuroTraits$Species
pyr$genus <- PyrenneesTraits$Genus %in% EuroTraits$Genus & !pyr$species
pyr$family <- PyrenneesTraits$Family %in% EuroTraits$Family & !pyr$species & !pyr$genus
pyr$order <- PyrenneesTraits$Order %in% EuroTraits$Order & !pyr$species & !pyr$genus  & !pyr$family
pyr$class <- !pyr$species & !pyr$genus  & !pyr$family & !pyr$order

# Serengeti:
ser <- data.frame(FW = rep("Serengeti", nrow(SerengetiTraits)))
ser$species = SerengetiTraits$Species %in% EuroTraits$Species
ser$genus <- SerengetiTraits$Genus %in% EuroTraits$Genus & !ser$species
ser$family <- SerengetiTraits$Family %in% EuroTraits$Family & !ser$species & !ser$genus
ser$order <- SerengetiTraits$Order %in% EuroTraits$Order & !ser$species & !ser$genus  & !ser$family
ser$class <- !ser$species & !ser$genus  & !ser$family & !ser$order

# High Arctic:
HA <- data.frame(FW = rep("HighArctic", nrow(HighArcticTraits)))
HA$species = HighArcticTraits$Species %in% EuroTraits$Species
HA$genus <- HighArcticTraits$Genus %in% EuroTraits$Genus & !HA$species
HA$family <- HighArcticTraits$Family %in% EuroTraits$Family & !HA$species & !HA$genus
HA$order <- HighArcticTraits$Order %in% EuroTraits$Order & !HA$species & !HA$genus  & !HA$family
HA$class <- !HA$species & !HA$genus  & !HA$family & !HA$order

# Put everything together
df <- rbind(pyr, ser, HA) %>%
  group_by(FW) %>%
  summarise(species = sum(species), genus = sum(genus), family = sum(family), order = sum(order), class = sum(class)) %>%
  pivot_longer(species:class, names_to = "Level") %>%
  mutate(Level = factor(Level, levels = c("class", "order", "family", "genus", "species")),
         FW = factor(FW, levels = c("Pyrennees", "HighArctic", "Serengeti")))
p1 <- ggplot(df, aes( y=value, x=FW)) + 
  geom_bar(position="fill", stat="identity", aes(fill  = Level)) +
  scale_fill_brewer(type = "seq") +
  geom_text(data = data.frame(FW = factor(c("Pyrennees", "HighArctic", "Serengeti")),
                              N = c(nrow(PyrenneesTraits), nrow(HighArcticTraits), nrow(SerengetiTraits))),
            aes(x = FW, y = 1.02, label = paste(N, "species"))) +
  labs(title = "Taxonomic overlap with European metaweb")

ggsave("figures/exploration/TaxoOverlap.png", p1, device = "png")

##### Trait coverage
rm(df, HA, pyr, ser)

# EuroMW
euro <- data.frame(FW = rep("Europe", nrow(EuroTraits)))
euro$TL = !is.na(EuroTraits$Trophic_level)
euro$ActivityTime = !is.na(EuroTraits$Diel_activity)
euro$Habitat = !is.na(EuroTraits$Forest)
euro$Size = !is.na(EuroTraits$Body_mass_g) | !is.na(EuroTraits$Body_length_mm) | !is.na(EuroTraits$Adult_svl_cm)
euro$ClutchSize = !is.na(EuroTraits$Litter_clutch_size)
euro$LifeHistory = !is.na(EuroTraits$Max_longevity_d) | !is.na(EuroTraits$Maturity_d) | !is.na(EuroTraits$Longevity_d) | !is.na(EuroTraits$Generation_length_d)

# Pyrennees
pyr <- data.frame(FW = rep("Pyrennees", nrow(PyrenneesTraits)))
pyr$TL = !is.na(PyrenneesTraits$Trophic_level)
pyr$ActivityTime = !is.na(PyrenneesTraits$Diel_activity)
pyr$Habitat = !is.na(PyrenneesTraits$Forest)
pyr$Size = !is.na(PyrenneesTraits$Body_mass_g) | !is.na(PyrenneesTraits$Body_length_mm) | !is.na(PyrenneesTraits$Adult_svl_cm)
pyr$ClutchSize = !is.na(PyrenneesTraits$Litter_clutch_size)
pyr$LifeHistory = !is.na(PyrenneesTraits$Max_longevity_d) | !is.na(PyrenneesTraits$Maturity_d) | !is.na(PyrenneesTraits$Longevity_d) | !is.na(PyrenneesTraits$Generation_length_d)

# Serengeti
ser <- data.frame(FW = rep("Serengeti", nrow(SerengetiTraits)))
ser$TL = !is.na(SerengetiTraits$Trophic_level)
ser$ActivityTime = !is.na(SerengetiTraits$Diel_activity)
ser$Habitat = !is.na(SerengetiTraits$Forest)
ser$Size = !is.na(SerengetiTraits$Body_mass_g) | !is.na(SerengetiTraits$Body_length_mm) | !is.na(SerengetiTraits$Adult_svl_cm)
ser$ClutchSize = !is.na(SerengetiTraits$Litter_clutch_size)
ser$LifeHistory = !is.na(SerengetiTraits$Max_longevity_d) | !is.na(SerengetiTraits$Maturity_d) | !is.na(SerengetiTraits$Longevity_d) | !is.na(SerengetiTraits$Generation_length_d)

# HighArctic
HA <- data.frame(FW = rep("HighArctic", nrow(HighArcticTraits)))
HA$TL = !is.na(HighArcticTraits$Trophic_level)
HA$ActivityTime = !is.na(HighArcticTraits$Diel_activity)
HA$Habitat = !is.na(HighArcticTraits$Forest)
HA$Size = !is.na(HighArcticTraits$Body_mass_g) | !is.na(HighArcticTraits$Body_length_mm) | !is.na(HighArcticTraits$Adult_svl_cm)
HA$ClutchSize = !is.na(HighArcticTraits$Litter_clutch_size)
HA$LifeHistory = !is.na(HighArcticTraits$Max_longevity_d) | !is.na(HighArcticTraits$Maturity_d) | !is.na(HighArcticTraits$Longevity_d) | !is.na(HighArcticTraits$Generation_length_d)

# Put everything together
df <- rbind(euro, pyr, ser, HA) %>%
  group_by(FW) %>%
  summarise(TL = mean(TL), ActivityTime = mean(ActivityTime), Habitat = mean(Habitat), 
            Size = mean(Size), ClutchSize = mean(ClutchSize), LifeHistory = mean(LifeHistory)) %>%
  pivot_longer(TL:LifeHistory, names_to = "Trait", values_to = "Coverage") %>%
  mutate(Trait = factor(Trait, levels = c("TL", "ActivityTime", "Habitat", "Size", "ClutchSize", "LifeHistory")),
         FW = factor(FW, levels = c("Europe", "Pyrennees", "HighArctic", "Serengeti")))
p2 <- ggplot(df, aes(y=Coverage, x=FW)) + 
  geom_bar(position=position_dodge(), stat="identity", aes(fill  = Trait)) +
  scale_fill_brewer(type = "qual") +
  coord_flip() +
  labs(title = "Taxonomic overlap with European metaweb")

ggsave("figures/exploration/TraitCoverage.png", p2, device = "png")
