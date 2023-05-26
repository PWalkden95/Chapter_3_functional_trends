## SCRIPT

## Development: Patrick Alexander Walkden

## Description:
## This script is going to create a map of functional endemism by taking each cell of occupied trait space and 
## assigning it a distictiveness score which is an inverse of the probability of occupancy, so that the more
## likley the cell is occupied the lower its score and vice versa. 

## This value is then evenly assigned to all raster cells defining the ranges of all species that could express
## that trait combination

rm(list = ls())


require(tidyverse)
require(terra)
require(sf)
require(raster)
source("functions/polygon_handling_functions.R")
require(doParallel)

## blank map of US to ptoject onto

continental_us <- readRDS("outputs/continental_US_polygon.rds")


blank_map <- terra::rast("../../../(2020-_Natural_History_Museum/Datasets/Environmental_Variables/wc2.1_30s_bio_6.tif")


blank_map <- terra::aggregate(x = blank_map, fact = 30)

crop_map <- blank_map %>% terra::mask(mask = continental_us) %>% terra::crop(y = c(-179.1435,0,25,83.116522)) %>%
  terra::classify(cbind(-Inf,Inf,0))

plot(crop_map)


## species polygon list

range_files <- list.files("../../../(2020-_Natural_History_Museum/Datasets/Birdlife_Maps/Shapefiles/PREDICTS_BL/",
                          full.names = TRUE)

species_from_path <- function(path){
  path_sp <- unlist(str_split(path, pattern = "/"))
  path_sp <- gsub(x = path_sp[length(path_sp)],pattern = "Birdlife_Range_data.rds",replacement = "")  
  
  return(path_sp)
  
}

range_file_species <- apply(data.frame(range_files),MARGIN = 1,FUN = species_from_path)



## load in the TPD occupancy matrix 

US_TPD <- readRDS("outputs/TPD_value_matrix.rds")




functional_endemism_map <- foreach(i = 1:nrow(US_TPD),
                                   .combine = "rbind",.inorder = FALSE,
                                   .packages = c("tidyverse","terra","sf")) %do% {

TPD_cell_weight <- 1 - US_TPD[i,"full_values"]

cell_species <- colnames(US_TPD)[which(as.numeric(US_TPD[i,c(5:ncol(US_TPD))]) > 0)]

## now get the polygons of the species that can express that trait combination

species_polygon <- c()



for(sp in cell_species){
  
  data <- readRDS(range_files[which(range_file_species == sp)])
  
  data <- data %>% dplyr::filter(presence %in% c(1,2,3),
                                 origin %in% c(1,2,3),
                                 seasonal %in% c(1,2,3))

  
  species_polygon <- rbind(species_polygon,data)
  
  }


shape <- terra::vect(spatial_combine_polygons(species_polygon))


sp_endemism_map <- terra::rasterize(x = shape,y = crop_map, fun = "sum", background = NA)


sp_endemism_map <- sp_endemism_map * (TPD_cell_weight/sum(terra::values(sp_endemism_map) > 0, na.rm = TRUE))



mat <- matrix(terra::values(sp_endemism_map),nrow = 1,dimnames = list(i,1:length(terra::values(crop_map))))



return(mat)

}




new_map <- crop_map

terra::values(new_map) <- colSums(functional_endemism_map, na.rm = TRUE)


