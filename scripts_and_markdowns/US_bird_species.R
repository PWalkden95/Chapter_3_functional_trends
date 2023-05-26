## SCRIPT 

## Development: Patrick Alexander Walkden

## This script is designed to identify priorty regions of the US for functional diversity. To generate this map
## I will create a map similar to the IUCNs range-rarity map that maps for each species the proportion of the 
## species's range contained within a given pixel.

## The functional diversity equivalent will be the proportion of a species' functional trait space contained 
## within a given pixel, then weighted by the proportion of occupancy of that cell. That's the theory anyway.

## First, because I'll be focusing on the continential US I want to just identify those species that overlap 
## with the US at somepoint of their range. So I need a polygon of the US and then overlap the species ranges
## maybe their are some ebird ranges too? but initially I will ocntinue just with Birdlife ranges. 

rm(list = ls())

require(terra)
require(tidyverse)
require(sf)
require(doParallel)


#### get the country shapefiles


countries <-
  sf::st_read("../../../(2020-_Natural_History_Museum/Datasets/Country_shapefiles/all_countries.shp")


continental_us <- sf::st_union(countries$geometry[c(234,40)]) %>% terra::vect() %>% 
  terra::crop(y = c(-179.1435,0,25,83.116522))


saveRDS(object =  continental_us, file = "outputs/continental_US_polygon.rds")


## there we have the continental US and canada. Now we want to go through the species ranges and idenitfy
## whether they overlap with this polygon

range_files <- list.files("../../../(2020-_Natural_History_Museum/Datasets/Birdlife_Maps/Shapefiles/PREDICTS_BL/",
                          full.names = TRUE)


### Function to combine polygons 

########## combine spatial polygons ready for analysis

spatial_combine_polygons <-
  function(geometry) {
    if (any(
      class(geometry$Shape)[1] == "sfc_MULTISURFACE",
      class(geometry$Shape)[1] == "sfc_GEOMETRY"
    )) {
      for (k in 1:NROW(geometry)) {
        geometry$Shape[[k]] <- st_cast(geometry$Shape[[k]], "MULTIPOLYGON")
      }
    }
    
    
    shape <-
      st_combine(geometry$Shape)
    
    
    
    st_crs(shape) <-
      "+proj=longlat +datum=WGS84 +no_defs"
    
    return(shape)
  }




US_species <- foreach(species = range_files,.combine = "c",.inorder = FALSE,.packages = c("tidyverse",
                                                                                          "sf",
                                                                                          "terra")) %do% {
                                                                                            


                                                                                            
## load in species range

data <- try(readRDS(species))


if(any(class(data) == "try-error")){
  path_chr <- unlist(str_split(species,pattern = "/"))
  path_chr <- path_chr[length(path_chr)]
  sp <- gsub(pattern = "Birdlife_Range_data.rds",replacement = "",x = path_chr)
  
  
  sp <- paste("error",sp,sep = "_")
  return(sp)
} else {
  


  
  
species_name <- data$binomial[1]





## filter polygons to those that have presence (1 Extant, 2 Probably extant, 3 possibly extant)
## filter polygons to those that have origin (1 Native, 2 Reintroduced, 3 introduced)
## filter polygons to those that have seasonality (1 Resident, 2 Breeding, 3 Non-breeding)


data <- data %>% dplyr::filter(presence %in% c(1,2,3),
                               origin %in% c(1,2,3),
                               seasonal %in% c(1,2,3))

if(nrow(data) == 0){
  return(NA)
} else {

shape <- terra::vect(spatial_combine_polygons(data))

# sf_use_s2(TRUE)
# if (!st_is_valid(shape)) {
#   sf_use_s2(FALSE)
# }


if(terra::is.related(shape,continental_us,"intersects")){
  sp <- species_name 
} else {
  sp <- NA
}

return(sp)

}

}
                                                                                          }


US_species <- na.omit(US_species)


error_species <- grep(pattern = "error",x = US_species,value = TRUE)


species_csv <- data.frame(species = gsub(x = error_species,pattern = "error_",replacement = ""))

write.csv(species_csv, file = "C:/Users/patri/Desktop/sp.csv")


US_species <- c(US_species,"Campylorhynchus brunneicapillus","Xanthocephalus xanthocephalus")

US_species <- grep(pattern = "error",x = US_species,value = TRUE, invert = TRUE)

## save species list 

write_rds(x = US_species, file = "outputs/continential_US_bird_species.rds")


