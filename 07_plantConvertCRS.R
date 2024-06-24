library(sf)
library(tidyverse)
library(scales)
library(rnaturalearth)
library(patchwork)

### function to clean and intersect polygons (species occurrence) ---
covertCRS <- function(PathtoPlant) {
    
    # Read and clean plant shapefile
    plant <- st_read(PathtoPlant) %>%
        # presence: 1--Extant, 2--Probably Extant, 3--Possibly Extant, 
        # 4--Possibly Extinct, 5--Extinct, 6--Presence Uncertain
        # origin: 1--Native, 2--Reintroduced, 3--Introduced, 
        # 4--Vagrant, 5--Origin Uncertain, 6--Assisted Uncertain
        filter(
            # !terrestial == 'true', 
            presence %in% c(1, 2, 3),
            !origin %in% c(3, 4))
    print('Here1!')
    # Transform and validate plant shapefile
    plant_esri <- st_transform(plant, "ESRI:54017")
    plant_valid <- st_make_valid(plant_esri)
    
    part <- basename(PathtoPlant)
    filename <- file.path(paste0('result/spResult/07_ESRI54017_', part))
    # Write the grid as a shapefile
    st_write(plant_valid, filename)
}

# plant range data
# plant_shapefile_paths <- c('data/spData/PLANTS/PLANTS_PART1.shp', # 5395 unique sp; 4933 terrestial
#                            'data/spData/PLANTS/PLANTS_PART2.shp')

plant_shapefile_paths <- c('data/spData/PLANTS/PLANTS_PART2.shp')
walk(plant_shapefile_paths, covertCRS)

df <- read_csv('data/spData/redlist_species_data/points_data.csv') # 36803; 36306 not in polygons
df2 <- df %>% 
    dplyr::select(sci_name, presence, origin, longitude, latitude) %>% 
    filter(presence %in% c(1, 2, 3),
           !origin %in% c(3, 4))

plantPoints.sf <- st_as_sf(df2, coords=c("longitude", "latitude"))
st_crs(plantPoints.sf) <- 4326
plantPoints.sf <- st_transform(plantPoints.sf, "ESRI:54017")
filename <- 'result/spResult/07_ESRI54017_plantPoint.shp'
st_write(plantPoints.sf, filename)


