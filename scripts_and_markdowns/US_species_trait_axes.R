## SCRIPT

## Development: Patrick Alexander Walkden

## Description: US species functional traits --> for the identified US species I need to get their traits in the 
##              format ready to generate the functional trait spaces.

rm(list = ls())

require(tidyverse)
require(magrittr)
require(TPD)


### US species 

US_species <- readRDS("outputs/continential_US_bird_species.rds")

## load in the trait data 

AVONET_full <-  read.csv("../../../(2020-_Natural_History_Museum/Datasets/GBD/GBD_BiometricsRaw_combined_11_Aug_2021_MASTER.csv")

species_specimens <- AVONET_full %>% dplyr::filter(Species1_BirdLife %in% US_species)


## are all species within AVONET


## Define teh traits that make up teh seperate partitions

foraging_traits = c("Beak.Length_Culmen", "Beak.Length_Nares" , "Beak.Width" ,"Beak.Depth")
locomotory_traits = c("Tarsus.Length", "Wing.Length","Secondary1","Tail.Length" )


## perform the two-step pca - some specimens will be dropped if after some imputation there are still incomplete cases

species_specimens <- species_specimens[,c("Species1_BirdLife",foraging_traits,locomotory_traits)]


source("functions/two_step_pca.R")

two_step_dataframe <- two_step_PCA(dataframe = species_specimens, means = FALSE,
                                   foraging_traits = foraging_traits,
                                   locomotory_traits = locomotory_traits)



#### for TPD computation each species needs greater than four specimens so lets check which have sufficient


n_specimens <- two_step_dataframe %>% dplyr::group_by(Birdlife_Name) %>% dplyr::mutate(n_specimen = n())


full_specimens <- n_specimens %>% dplyr::filter(n_specimen >= 4) 

full_species <- full_specimens %>% dplyr::distinct(Birdlife_Name) %>% pull()
## so we have 3532 species that have sufficient specimens to create full TPDs



### how many with fewer that 4 but greater than one specimen
fewer_specimens <- n_specimens %>% dplyr::filter(n_specimen < 4 & n_specimen > 1)

fewer_species <- fewer_specimens %>% dplyr::distinct(Birdlife_Name) %>% pull()
### 32 species with fewer

### now only those with a single specimen
single_specimen <- n_specimens %>% dplyr::filter(n_specimen == 1)

single_species <- single_specimen %>% dplyr::distinct(Birdlife_Name) %>% pull()
## nine species with a single specimen

791+ 29 + 26 



### this doesn't add up because some species don't have any trait data so data is taken from a very close relative

BL_traits <- read.csv("../../../(2020-_Natural_History_Museum/Datasets/GBD/AVONET1_BirdLife.csv") %>%
  dplyr::mutate(Species1 = ifelse(Species1 == "Gorsachius magnificus", "Oroanassa magnifica", paste(Species1)),
                Species1 = ifelse(Species1 == "Psittacula krameri", "Alexandrinus krameri", paste(Species1)),
                Species1 = ifelse(Species1 == "Arachnothera hypogrammica", "Kurochkinegramma hypogrammica", paste(Species1)),
                Species1 = ifelse(Species1 == "Psittacula eupatria", "Palaeornis eupatria", paste(Species1)),
                Species1 = ifelse(Species1 == "Psittacula finschii", "Himalayapsitta finschii", paste(Species1)),
                Species1 = ifelse(Species1 == "Psittacula roseata", "Himalayapsitta roseata", paste(Species1)),
                Species1 = ifelse(Species1 == "Psittacula longicauda", "Belocercus longicaudus", paste(Species1)),
                Species1 = ifelse(Species1 == "Hapalocrex flaviventer", "Laterallus flaviventer", paste(Species1)),
                Species1 = ifelse(Species1 == "Hylocharis chrysura", "Amazilia chrysura", paste(Species1)),
                Species1 = ifelse(Species1 == "Porzana spiloptera", "Laterallus spilopterus", paste(Species1)),
                Species1 = ifelse(Species1 == "Hylocharis cyanus", "Amazilia cyanus", paste(Species1)),
                Species1 = ifelse(Species1 == "Chrysuronia oenone", "Amazilia oenone", paste(Species1)),
                Species1 = ifelse(Species1 == "Hylocharis eliciae", "Amazilia eliciae", paste(Species1)),
                Species1 = ifelse(Species1 == "Juliamyia julie", "Amazilia julie", paste(Species1)),
                Species1 = ifelse(Species1 == "Lepidopyga goudoti", "Amazilia goudoti", paste(Species1)),
                Species1 = ifelse(Species1 == "Lepidopyga coeruleogularis", "Amazilia coeruleogularis", paste(Species1)),
                Species1 = ifelse(Species1 == "Psittacula himalayana", "Himalayapsitta himalayana", paste(Species1)),
                Species1 = ifelse(Species1 == "Psittacula columboides", "Nicopsitta columboides", paste(Species1)),
                Species1 = ifelse(Species1 == "Psittacula cyanocephala", "Himalayapsitta cyanocephala", paste(Species1)),
                Species1 = ifelse(Species1 == "Psittacula calthrapae", "Nicopsitta calthrapae", paste(Species1)))



no_specimen_species <- US_species[!(US_species %in% c(full_species,fewer_species,single_species))]



imputed_traits <- BL_traits %>%
  dplyr::filter(Species1 %in% no_specimen_species)



### okay will need to add this into the species specimen data frame and perform the two step again.

colnames(imputed_traits)[2] <- "Species1_BirdLife"
imputed_traits <- imputed_traits[,colnames(species_specimens)]


### if there are multiple specimens for some species but they all are missing values for a single trait we can
### possibly impute some more values for the multiple specimens


species_specimens <- species_specimens %>% dplyr::filter(!(Species1_BirdLife %in% no_specimen_species)) %>%
  rbind(imputed_traits)


two_step_dataframe <- two_step_PCA(species_specimens,foraging_traits = foraging_traits,locomotory_traits = locomotory_traits,
                                   means = FALSE)


### are all species now considered?

all(US_species %in% two_step_dataframe$Birdlife_Name)


single_species <- c(single_species,no_specimen_species)

######
### Now to calculate the TPDs to see whether all can be assessed

## this contains a lot of useful functions for computing TPDs

source("functions/TPD_computation_functions.R")

### trait ranges define the size of the morphospace that species are projected into
### I use 2.5% buffer on all species traits 

trait_ranges <- trait_range_calculation(traits = full_specimens, range = 0.025)


## alpha dictates the amount of cut off point for the TPD so we are takking 99% percentile of the derived probability density functions. 

full_tpds <-
  TPDs(
    species = full_specimens$Birdlife_Name,
    traits = full_specimens[, c(2:4)],
    trait_ranges = trait_ranges, alpha = 0.99
  )


check_species <- TPD_check(full_tpds)

## good all full species can be projected



fewer_specimens <- fewer_specimens %>% dplyr::group_by(Birdlife_Name) %>% dplyr::summarise(meanfor = mean(foraging),
                                                                                           meanloco = mean(locomotory),
                                                                                           meanbody = mean(body_size),
                                                                                           sdfor = sd(foraging),
                                                                                           sdloco = sd(locomotory),
                                                                                           sdbody = sd(body_size))


fewer_tpds <- TPDsMean(species = fewer_specimens$Birdlife_Name, means = fewer_specimens[,c(2:4)],
                       sds = fewer_specimens[,c(5:7)], alpha = 0.99, trait_ranges = trait_ranges)


mean_check_species <- TPD_check(fewer_tpds)

## mean check species -- identified no problematic species


### unconstrained bandwidth estimator as used by Carmona et al.based on all specimens
single_bandwidth <- sqrt(diag(Hpi.diag(two_step_dataframe[,c(2:4)])))


## get the data in the format necessary taking the mean values of the specimens available
single_traits <- two_step_dataframe %>% dplyr::filter(Birdlife_Name %in% single_species) %>% dplyr::group_by(Birdlife_Name) %>%
  dplyr::summarise(meanfor = mean(foraging),
                   meanloco = mean(locomotory),
                   meanbody = mean(body_size)) %>% dplyr::ungroup() %>%
  dplyr::mutate(sdfor = single_bandwidth[1],
                sdloco = single_bandwidth[2],
                sdbody = single_bandwidth[3])


## calculate TPDS
single_tpds <- TPDsMean(species = single_traits$Birdlife_Name, means = single_traits[,c(2:4)], sds = single_traits[,c(5:7)],
                        trait_ranges = trait_ranges, alpha = 0.99)


## great it's got nothing in it!

final_tpd_check <- TPD_check(single_tpds)

### none -  brilliant!



US_trait_list <- list()

US_trait_list$complete_traits <- data.frame(full_specimens)
US_trait_list$partial_traits <- data.frame(fewer_specimens)
US_trait_list$single_traits <- data.frame(single_traits)
US_trait_list$full_specimen_dataframe <- two_step_dataframe


### final checks to see whether all species are represented and only represented once.


#### are all species considered here 

full_sp <- unique(US_trait_list$complete$Birdlife_Name)
partial_sp <- unique(US_trait_list$partial$Birdlife_Name)
single_sp <- unique(US_trait_list$single$Birdlife_Name)

all_species <- c(full_sp,partial_sp,single_sp)

## great

all(all_species %in% US_species)


## all species just appear once?

any(full_sp %in% partial_sp)
any(full_sp %in% single_sp)
any(partial_sp %in% single_sp)


### great, now save

write_rds(file = "outputs/US_species_full_traits.rds", US_trait_list)


