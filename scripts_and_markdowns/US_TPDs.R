## SCRIPT

## Development: Patrick Alexander Walkden

## Description:
## Script to generate the full species hypervolume based on the US species list and the trait values gathered in the
## previous script

rm(list = ls())

require(tidyverse)
require(TPD)
source("functions/TPD_computation_functions.R")
require(doParallel)


### US_species & traits

US_species <- readRDS("outputs/continential_US_bird_species.rds")

US_traits <- readRDS("outputs/US_species_full_traits.rds")


### let's make the full hypervolume for the US species 

trait_ranges <- get_species_trait_ranges(species = US_species,traits = US_traits, range = 0.055)


full_TPD <- species_TPD(species = US_species,trait_ranges = trait_ranges,traits = US_traits)


### we also want to create a species x cell matrix of this hypervolume 

## get the cells of the morphological trait space to which each species was projected onto 

eval_grid <- full_TPD$data$evaluation_grid %>% set_colnames(c("locomotion","foraging","body"))


## using 8 cores combine all the species TPDs into a single matrix 

registerDoParallel(cores = 8)


sp_trait_density <- foreach(sp = names(full_TPD$TPDs),
                            .combine = "cbind",
                            .inorder = FALSE,
                            .packages = c("tidyverse")) %dopar%{
                              
                              mat <- matrix(full_TPD[["TPDs"]][[sp]],
                                            ncol = 1)
                              
                              
                              colnames(mat) <- sp
                              
                              
                              return(mat)
                              
                              
                            }

registerDoSEQ()
closeAllConnections()


comm <- matrix(rep(1,length(US_species)), nrow = 1, dimnames = list("comm",US_species))

TPD_full <- TPDc(full_TPD, comm)

full_occupancy <- matrix(TPD_full$TPDc$TPDc$comm, ncol = 1, dimnames = list(c(1:125000),"full_values"))


##combine grid and species occupancy probabilities 
species_trait_density <- cbind(eval_grid,full_occupancy, sp_trait_density)

### then filter out those cells which do no have any occupancy
species_trait_density <-  species_trait_density[which(rowSums(species_trait_density[,-c(1:3)]) >  0),]


write_rds(species_trait_density, file = "outputs/TPD_value_matrix.rds")

