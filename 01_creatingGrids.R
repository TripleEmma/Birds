### created by Emma
### last modified: 2024-04-24
### create multiple grids with different resolutions


### load libraries ---
library(tidyverse)
library(units)
library(sf)

# rm(list = ls())

# Define the world boundary polygon
worldSingleCell <- st_polygon(list(
    rbind(c(-180, -90),
          c(180, -90),
          c(180, 90),
          c(-180, 90),
          c(-180, -90)
    )
))
# convert to sfc object with CRS as "EPSG:4326"
worldSingleCell_sfc <- st_sfc(worldSingleCell, crs = "EPSG:4326")

# Function to create grids
gridsList <- function(edgeLen, worldSingleCell_sfc){
    
    # set the length of edge
    edge <- units::as_units(edgeLen, "km")
    # Create grid
    grid <- worldSingleCell_sfc %>%  
        st_transform("ESRI:54017") %>% 
        st_make_grid(cellsize = c(edge, edge)) %>% 
        st_sf() %>% 
        st_cast()
    grid$gridNo <- paste0('grid', seq(nrow(grid)))
    return(grid)
}

# Define edge lengths
edges <- seq(10, 800, 15)

# Create list of grids with varying edge lengths
gridsLists <- map(edges, gridsList, worldSingleCell_sfc)

output_dir <- "result/gridsShapeResult"
# Create the directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE)

# Function to write grid as a shapefile
write_grid_shapefile <- function(grid, edge) {
    # Define the filename for the shapefile
    filename <- file.path(output_dir, paste0('edges', edge, ".shp"))
    # Write the grid as a shapefile
    write_sf(grid, filename)
}

# Use walk2() to iterate over the list of grids and write each grid as a shapefile
walk2(gridsLists, edges, write_grid_shapefile)


# grid %>%
#     st_coordinates() %>%
#     as_tibble() %>%
#     slice(seq(1, n(), 5)) %>%
#     mutate(area = as.numeric(st_area(grid))) %>%
#     ggplot(aes(y = area, x = Y)) + 
#     geom_point() +
#     labs(x = 'latitude', title = 'Area of grid cells as a function of latitude',
#          subtitle = 'Using an equal-area projection as basis for the grid.')
