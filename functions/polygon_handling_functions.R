## FUNCTIONS

## Development: Patrick Alexander Walkden

## Description: Functions that help the handling of polygons of the birdlife maps 


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

