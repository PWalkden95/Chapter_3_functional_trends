#### This is a script that will contain functions that are relevant to computing TPDs and alter the trait ranges,
#### test whether some species were able to be computed, and others



## INPUT: desired buffer to add to trait ranges for the TPD computation ; trait list for the species that TPDs are t be computed for

## OUTPUT: a list of length three containing the trait range values to input into the TPD function


trait_range_calculation <- function(range = 0.15, traits) {
  n_range_2 <- dist(c(min(traits[, 2]), max(traits[, 2])))[1]
  n_range_3 <- dist(c(min(traits[, 3]), max(traits[, 3])))[1]
  n_range_4 <- dist(c(min(traits[, 4]), max(traits[, 4])))[1]
  
  trait_ranges <-
    list(c(
      min(traits[, 2]) - (range * n_range_2),
      max(traits[, 2]) + (range * n_range_2)
    ),
    c(
      min(traits[, 3]) - (range * n_range_3),
      max(traits[, 3]) + (range * n_range_3)
    ),
    c(
      min(traits[, 4]) - (range * n_range_4),
      max(traits[, 4]) + (range * n_range_4)
    ))
  return(trait_ranges)
}


###################### Function to find out which species a TPD couldn't be calculated for with the range and division combination given


## INPUT: species TPDs that have been calculated usin gthe TPDs or TPDSmean functions in the tPD package
## OUTPUT: a character vector of the species that do not create viable TPDs for

## DESCRIPTION: a function to check whether the TPDs creates viable hypervolumes


TPD_check <- function(TPD) {
  Trait_density <- TPD
  
  fucked_species <- c()
  for (i in 1:length(Trait_density$TPDs)) {
    if (any(is.nan(Trait_density$TPDs[[i]]))) {
      fucked_species <- c(fucked_species, names(Trait_density$TPDs)[[i]])
    }
  }
  return(fucked_species)
}


##################################
##################################

### need study species trait ranges

### a broader function that gets a group of species and goes through the three different
### trait lists to get the trait ranges

## INPUT: species vecotr ; trait list; and desired trait ranges

## OUTPUT: list of length three with the desired trait ranges



get_species_trait_ranges <- function(species, traits, range = 0.15) {
  full_species <- unique(traits$complete$Birdlife_Name)
  partial_species <- unique(traits$partial$Birdlife_Name)
  single_species <- unique(traits$single$Birdlife_Name)
  
  species_traits <- c()
  
  for (sp in species) {
    if (sp %in% full_species) {
      species_traits <-
        rbind(species_traits,
              data.frame(Birdlife_Name = sp, traits$complete[traits$complete$Birdlife_Name == sp, c(2:4)]))
    }
    
    if (sp %in% partial_species) {
      species_traits <-
        rbind(
          species_traits,
          data.frame(
            Birdlife_Name = sp,
            foraging = traits$partial[traits$partial$Birdlife_Name == sp, 2],
            locomotory = traits$partial[traits$partial$Birdlife_Name == sp, 3],
            body_size = traits$partial[traits$partial$Birdlife_Name == sp, 4]
          )
        )
    }
    
    if (sp %in% single_species) {
      species_traits <-
        rbind(
          species_traits,
          data.frame(
            Birdlife_Name = sp,
            foraging = traits$single[traits$single$Birdlife_Name == sp, 2],
            locomotory = traits$single[traits$single$Birdlife_Name == sp, 3],
            body_size = traits$single[traits$single$Birdlife_Name == sp, 4]
          )
        )
    }
    
  }
  
  study_species_trait_ranges <-
    trait_range_calculation(range = range, traits = species_traits)
  
  
  return(study_species_trait_ranges)
  
}


##################################
#################################  Species TPD function


## function that allows for both TPDs and TPDsmean TPD computation and combines the
## resulting TPDs.

## INPUT: species vector; trait ranges; traits

## OUTPUT: TPDs class object with all the species TPDs


species_TPD <- function(species, trait_ranges, traits) {
  ## if sds and therefore using species information
  
  
  
  
  mtpd <- FALSE ### do we have to do meanTPDs
  
  ### are there any species that require mean tpd calculation
  if (any(
    species %in% c(
      traits$partial_traits$Birdlife_Name,
      traits$single_traits$Birdlife_Name
    )
  )) {
    # YES
    mtpd <- TRUE
    
    ### what species are they
    Mean_sp <-
      species[species %in% c(traits$partial_traits$Birdlife_Name,
                             traits$single_traits$Birdlife_Name)]
    
    
    ## combine all partial traits into the same df representing means and bandwidths to use filter for the species in the community and convert to df
    
    mean_TPD_dat <-
      rbind(traits$partial_traits, traits$single_traits) %>% dplyr::filter(Birdlife_Name %in% Mean_sp) %>% data.frame()
    
    
    ### then calculate mean TPDs
    
    mean_TPD <-
      TPDsMean(
        species = mean_TPD_dat[, 1],
        means = mean_TPD_dat[, c(2:4)],
        sds = mean_TPD_dat[, c(5:7)],
        trait_ranges = trait_ranges
      )
    
    check <- TPD_check(mean_TPD)
    
    
    if (!is_null(check)) {
      print("print CHECK ISSUE -- MEAN")
      
      stop()
    }
    
  }
  
  ## for those species in the community that have sufficient traits extract their traits
  
  comm_traits <-
    traits[["complete_traits"]] %>% dplyr::filter(Birdlife_Name %in% species)
  
  ## claculate TPD
  
  trait_density <-
    TPDs(species = comm_traits[, 1],
         traits = comm_traits[, c(2:4)],
         trait_ranges = trait_ranges)
  
  
  check <- TPD_check(trait_density)
  
  
  CHECK <- FALSE
  if (!is_null(check)) {
    print("print CHECK ISSUE -- COMPLETE")
    
    CHECK <- TRUE
    
    
    check_df <- comm_traits[comm_traits$Birdlife_Name %in% check,] %>% dplyr::group_by(Birdlife_Name) %>%
      dplyr::summarise(meanfor = mean(foraging), meanloco = mean(locomotory), meanbody = mean(body_size),
                       sdforag = sd(foraging), sdloco = sd(locomotory), sdbody = sd(body_size)) %>% data.frame()
    
    
    
    check_TPD <- TPDsMean(
      species = check_df[,1],
      means = check_df[,c(2:4)],
      sds = check_df[,c(5:7)],
      trait_ranges = trait_ranges
    )
    
    
  }
  
  
  
  
  
  
  ## if there were TPDs calculated using the mean and sd join up the data so that both are represented in the same list
  
  if (mtpd) {
    trait_density$data$species <-
      c(trait_density$data$species, mean_TPD$data$species)
    trait_density$TPDs <- c(trait_density$TPDs, mean_TPD$TPDs)
    trait_density$data$traits <-
      rbind(trait_density$data$traits, mean_TPD$data$means)
  }
  
  if(CHECK){
    
    
    check_rows <- which(trait_density$data$species %in% check)
    
    trait_density$data$species <- trait_density$data$species[-check_rows]
    
    trait_density$TPDs <-  trait_density$TPDs[-which(names(trait_density$TPDs) %in% check)]
    
    trait_density$data$traits <- trait_density$data$traits[-check_rows,]
    
    
    
    trait_density$data$species <-
      c(trait_density$data$species, check_TPD$data$species)
    trait_density$TPDs <- c(trait_density$TPDs, check_TPD$TPDs)
    trait_density$data$traits <-
      rbind(trait_density$data$traits, check_TPD$data$means)
    
    
    
  }
  
  
  
  return(trait_density)
}

##########################################
########################################

### getting the beta-diversity sites community data

## INPUT: the sites that are to be compared; predicts database; trait list 

## OUTPUT: a dataframe of three columns that contain information on the dissimilarity
## between communities and the proportion of that dissimilarity that is attributed to
## nestedness and turnover.




get_beta_diversity <- function(site_one,site_two, predicts, traits){
  
  data <- predicts %>% dplyr::filter(SSBS %in% c(site_one,site_two)) %>% dplyr::select(Birdlife_Name,RelativeAbundance,SSBS)
  
  site_species <- unique(data$Birdlife_Name)
  
  beta_trait_ranges <- get_species_trait_ranges(species = site_species,traits = traits,range = 0.025)
  
  
  comm_data <- data.frame(Birdlife_Name = site_species) %>% 
    dplyr::left_join(data[data$SSBS == site_one, c("Birdlife_Name","RelativeAbundance")]) %>% 
    dplyr::rename(site_one = RelativeAbundance) %>%
    dplyr::left_join(data[data$SSBS == site_two, c("Birdlife_Name","RelativeAbundance")]) %>%
    dplyr::rename(site_two = RelativeAbundance)
  
  comm_data$site_one <- ifelse(is.na(comm_data$site_one), 0, comm_data$site_one)
  comm_data$site_two <- ifelse(is.na(comm_data$site_two), 0, comm_data$site_two)
  
  rownames(comm_data) <- comm_data$Birdlife_Name
  comm_data <- comm_data %>% dplyr::select(-Birdlife_Name)
  
  
  beta_species <- data %>% dplyr::filter(SSBS %in% c(site_one,site_two)) %>% dplyr::distinct(Birdlife_Name) %>% pull()
  
  beta_species_tpd <- species_TPD(species = beta_species,trait_ranges = beta_trait_ranges, traits = traits)
  
  comm <- t(comm_data[beta_species,])
  
  if(any(as.numeric(rowSums(comm)) == 0)){
    beta_frame <- data.frame(dissimilarity = NA,
                             beta_shared = NA,
                             beta_not_shared = NA)
    
    
    return(beta_frame)
  }
  
  community_tpds <- TPDc(beta_species_tpd, t(comm_data[beta_species,]))
  
  overlap <- TPD::dissim(community_tpds)
  
  beta_frame <- data.frame(dissimilarity = overlap$communities$dissimilarity[2],
                           beta_shared = overlap$communities$P_shared[2],
                           beta_not_shared = overlap$communities$P_non_shared[2])
  
  
  return(beta_frame)
  
}
