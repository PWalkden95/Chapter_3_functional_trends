### This is a script for a function that conducts the two-step PCA using the bird morphological
## traits. 

# INPUT: dataframe with species name and morphological traits.
#OUTPUT: dataframe with the three meaningful axes of variation: locomotory, forgaing and size niche


## the function can work with specimen data or just trait means. trait means just simple performs the 
## two-step PCA on the trait values while dealing with the specimen data the function attempts to fill
## in gaps in the specimen data (because it only takes full cases and we want to use as much of the information
## as possible). To do this the value is filled in as the mean of all other available specimens.


two_step_PCA <- function(dataframe, means = TRUE, 
                         foraging_traits,
                         locomotory_traits){

  colnames(dataframe)[1] <- "Birdlife_Name"
  ### if we are just evaluating trait means then it is easy to incorporate as we don't have to deal with
  ### missing data for some specimens.
  
  if(means) {
  
  for_data <- prcomp(scale(log(dataframe[,foraging_traits]))) ## foraging PCA -- values scaled and logged before
  
  
  loco_data <- prcomp(scale(log(dataframe[,locomotory_traits]))) ## locomotry PCA -- values scaled and logged before
  
  
  body_size_pc <- prcomp(data.frame(for_data$x[,1],loco_data$x[,1]))  ## PCA on PC1 of ecah previous PCA
    
  
  
  pca_data <- data.frame(Birdlife_Name = dataframe[,1], ### combine 
                         foraging = for_data$x[,2], 
                         locomotory = loco_data$x[,2],
                         body_size = body_size_pc$x[,1])
    
return(pca_data)
} else  {
  
  
  
  drop_rows <- c() #### rows where specimens have no trait information so drop 
  for(i in 1:nrow(dataframe)){
    if(all(is.na(dataframe[i,-1]))){
      drop_rows <- c(drop_rows,i)  
    }
  }
  
  
  
  dataframe <- dataframe[-drop_rows,]
  
  ## now we need to consider when we have multiple specimens how we deal with missing data. Especially,
  ## when specimens have the data for the majority of traits but missing one. The PCA requires complete
  ## cases so it would be a shame to lose data. As a fix those specimens that are missing data for some traits
  ## we will insert the mean value of the other specimens.
  
    for(trait in c(foraging_traits,locomotory_traits)){
        
    data <- dataframe %>% dplyr::select(Birdlife_Name,all_of(trait))
    
    na_sp <- unique(data[is.na(data[,2]),1])
  
    for(species in na_sp){
      
      species_trait_mean <- mean(((data %>% dplyr::filter(Birdlife_Name == species)))[,2],na.rm = TRUE)
      
 
    
      dataframe[dataframe$Birdlife_Name == species & is.na(dataframe[,trait]),trait] <- species_trait_mean
        
    }
  }
  
  ## okay now that we have the imputed data -- I don't know whether it's the best way but it's been done.
  ## then it is how we work with the differential number of specimens. This will actually be done later
  ## just do the two-step bit for now.
  
  
  ##  get rid of the specimens that still have the any NAs within their traits 
  
  dataframe <- dataframe %>% na.omit() %>% dplyr::group_by(Birdlife_Name) %>%
    dplyr::mutate(number_of_specimens = n())
  
  
  ### then do PCA as previous
  
  for_data <- prcomp(scale(log(dataframe[,foraging_traits])))
  
  
  loco_data <- prcomp(scale(log(dataframe[,locomotory_traits])))
  
  
  body_size_pc <- prcomp(data.frame(for_data$x[,1],loco_data$x[,1])) 
  
  
  
  pca_data <- data.frame(Birdlife_Name = dataframe[,1], 
                         foraging = for_data$x[,2], 
                         locomotory = loco_data$x[,2],
                         body_size = body_size_pc$x[,1])
  
  return(pca_data)
  
  
}
}
