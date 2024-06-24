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
    plant_valid <- read_sf(PathtoPlant) 
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
plant_shapefile_paths <- c('result/spResult/07_ESRI54017_PLANTS_PART1.shp', # 5395 unique sp; 4933 terrestial
                           'result/spResult/07_ESRI54017_PLANTS_PART2.shp') # 5286 unique sp; 4829 terrestial
overlap_selected_list <- map(plant_shapefile_paths, gridSP, grid)
rangeSP <- bind_rows(overlap_selected_list) %>% 
    distinct(gridNo, sci_name, .keep_all = TRUE)

# plant point data
plantPoints.sf <- read_sf('result/spResult/07_ESRI54017_plantPoint.shp')
grid_plant <- st_join(plantPoints.sf, grid)
pointSP <- st_drop_geometry(grid_plant) %>% 
    filter(!is.na(gridNo)) %>% 
    distinct(gridNo, sci_name, .keep_all = TRUE)

gridSP <- full_join(rangeSP, pointSP, 
                    by = c("gridNo", "sci_name", "presence", "origin")) %>% 
    distinct(gridNo, sci_name, .keep_all = TRUE)
gridSP$geometry <- NULL
# outputFile = paste0('result/spResult/07_speciesPerGrid.csv')
outputFile <- 'result/spResult/07_plantPerGrid280.csv'
write_csv(gridSP, outputFile)
gridSP <- read_csv('result/spResult/07_plantPerGrid280.csv')
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
                         palette = 'GnBu',
                         breaks = pretty_breaks(),
                         direction = 1) +
    theme(legend.position = 'bottom',
          plot.title = element_text(size = 18)) +
    guides(fill = guide_colorbar(barwidth = 20, barheight = 1)) +
    ggtitle('Global Plant Species Richness')
ggsave('result/07_plantSReveryGrid280.png', width = 11.9, height = 8.39)



## combine both SR and GD ----
# load genetic diversity
GDgrids <- read_delim('result/genResult/05_piEveryGrid', delim = '\t') %>%
    filter(GD < 0.4)
bothDiv_shp <- SRgrids_shp %>%
    left_join(GDgrids) %>%
    filter(!is.na(GD))
# Load world map data
world_map <- ne_countries(scale = "medium", returnclass = "sf")
world_map <- st_transform(world_map, "ESRI:54017")

p2 <- ggplot() +
    geom_sf(data = world_map, fill = "transparent", color = "black") +
    geom_sf(data = bothDiv_shp, aes(fill = SR)) +
    theme_minimal() +
    scale_fill_distiller(name = 'Species Richness',
                         palette = 'GnBu',
                         breaks = pretty_breaks(),
                         direction = 1) +
    theme(legend.position = 'bottom',
          plot.title = element_text(size = 18)) +
    guides(fill = guide_colorbar(title.position = "top",
                                 title.vjust = 0.5,
                                 title.hjust = 0.5,
                                 barwidth = 20, barheight = 1))  +
    ggtitle('Plant Species Richness')


filename = paste0('result/07_plantGD_Grid280.png')
p1/p2
ggsave(filename, width = 8.96, height = 9.72)


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
