### created by Emma
### last modified: 2024-04-25
### count species richness per grid using IUCN data 
## and also plot both diversity

### load libraries ----
library(sf)
library(tidyverse)
library(scales)
library(rnaturalearth)
library(patchwork)

### function to clean and intersect polygons (species occurrence) ---
gridSP <- function(PathtoPlant, grids) {
    
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
    print('Here2!')
    # Perform intersection
    overlaps <- st_intersection(grids, plant_valid) %>% 
        st_collection_extract(type = c("POLYGON"))
    print('Here3!')
    # Select relevant columns
    overlap_selected <- overlaps %>%
        st_drop_geometry() %>% 
        select(gridNo, sci_name, presence, origin, class, order_, family, genus, terrestial) %>%
        distinct(gridNo, sci_name, .keep_all = TRUE) 
    
    return(overlap_selected)
}


### load grids data ----
# shapeFile <- "result/gridsShapeResult/edges250.shp"
shapeFile <- "result/gridsShapeResult/edges280.shp"
grid <- read_sf(shapeFile)

# plant range data
plant_shapefile_paths <- c('data/spData/PLANTS/PLANTS_PART1.shp', # 5395 unique sp; 4933 terrestial
                           'data/spData/PLANTS/PLANTS_PART2.shp') # 5286 unique sp; 4829 terrestial
# plant1 <- st_read('data/spData/PLANTS/PLANTS_PART1.shp')
# length(unique(plant1$sci_name))
# (plant1_terrestial <- as_tibble(plant1) %>% 
#         filter(terrestial=='true') %>% 
#         select(sci_name) %>% pull() %>% 
#         unique() %>% 
#         length())
# 
# plant2 <- st_read('data/spData/PLANTS/PLANTS_PART2.shp')
# length(unique(plant2$sci_name))
# (plant2_terrestial <- as_tibble(plant2) %>% 
#         filter(terrestial=='true') %>% 
#         select(sci_name) %>% pull() %>% 
#         unique() %>% 
#         length())

# Apply function to each plant shapefile
overlap_selected_list <- map(plant_shapefile_paths, gridSP, grid)
rangeSP <- bind_rows(overlap_selected_list) %>% 
    distinct(gridNo, sci_name, .keep_all = TRUE)

# plant point data
df <- read_csv('data/spData/redlist_species_data/points_data.csv') # 36803; 36306 not in polygons
df2 <- df %>% 
    dplyr::select(sci_name, presence, origin, longitude, latitude) %>% 
    filter(presence %in% c(1, 2, 3),
           !origin %in% c(3, 4))

plantPoints.sf <- st_as_sf(df2, coords=c("longitude", "latitude"))
st_crs(plantPoints.sf) <- 4326
plantPoints.sf <- st_transform(plantPoints.sf, "ESRI:54017")
grid_plant <- st_join(plantPoints.sf, grid)
pointSP <- st_drop_geometry(grid_plant) %>% 
    filter(!is.na(gridNo)) %>% 
    distinct(gridNo, sci_name, .keep_all = TRUE)

gridSP <- full_join(rangeSP, pointSP, 
                    by = c("gridNo", "sci_name", "presence", "origin")) %>% 
    distinct(gridNo, sci_name, .keep_all = TRUE)
gridSP$geometry <- NULL
# outputFile = paste0('result/spResult/07_speciesPerGrid.csv')
outputFile = paste0('result/spResult/07_speciesPerGrid280.csv')
write_csv(gridSP, outputFile)

SRgrids <- gridSP %>% group_by(gridNo) %>% summarise(SR=n())
SRgrids_shp <- grid %>% left_join(SRgrids, by = 'gridNo')
# Load world map data
world_map <- ne_countries(scale = "medium", returnclass = "sf")
world_map <- st_transform(world_map, "ESRI:54017")
p1 <- ggplot() + 
    geom_sf(data = world_map, fill = "transparent", color = "black") +
    geom_sf(data = SRgrids_shp, aes(fill = SR)) +
    theme_minimal() +
    scale_fill_distiller(name = 'Species Richness',
                         palette = "YlGnBu",
                         breaks = pretty_breaks(),
                         direction = 1) +
    theme(legend.position = 'bottom',
          plot.title = element_text(size = 18)) +
    guides(fill = guide_colorbar(barwidth = 20, barheight = 1)) +
    ggtitle('Global plant species richness')
ggsave('result/07_SReveryGrid280.png', width = 11.9, height = 8.39)



## combine both SR and GD ----
# load genetic diversity
# GDgrids <- read_delim('result/genResult/05_piEveryGridSeedOnlyGrid325.txt', delim = '\t') %>% 
#     filter(GD < 0.4) 
# bothDiv_shp <- SRgrids_shp %>% 
#     left_join(GDgrids) %>%
#     filter(!is.na(GD))
# # Load world map data
# world_map <- ne_countries(scale = "medium", returnclass = "sf")
# world_map <- st_transform(world_map, "ESRI:54017")
# # plot
# 
# p2 <- ggplot() + 
#     geom_sf(data = world_map, fill = "transparent", color = "black") +
#     geom_sf(data = bothDiv_shp, aes(fill = GD)) +
#     theme_minimal() +
#     scale_fill_distiller(name = expression(paste("Genetic Diversity ", pi)),
#                          palette = "YlGnBu",
#                          breaks = pretty_breaks(),
#                          direction = 1) +
#     theme(legend.position = 'bottom',
#           plot.title = element_text(size = 18)) +
#     guides(fill = guide_colorbar(barwidth = 20, barheight = 1)) +
#     ggtitle('Global genetic diversity')
# 
# p3 <- ggplot() +
#     geom_sf(data = world_map, fill = "transparent", color = "black") +
#     geom_sf(data = bothDiv_shp, aes(fill = SR)) +
#     theme_minimal() +
#     scale_fill_distiller(name = 'Species Richness',
#                          palette = "YlGnBu",
#                          breaks = pretty_breaks(),
#                          direction = 1) +
#     theme(legend.position = 'bottom',
#           plot.title = element_text(size = 18)) +
#     guides(fill = guide_colorbar(barwidth = 20, barheight = 1)) +
#     ggtitle('Global species richness')
# p2/p3
# ggsave('result/bothDivSeedOnlyGrid325.png', width = 11.9, height = 8.39)
