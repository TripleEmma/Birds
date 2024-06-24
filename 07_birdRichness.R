library(sf)
library(terra)
library(tidyverse)
library(scales)
library(rnaturalearth)
library(patchwork)
library(geodata)


rm(list = ls())

birdsESRI54017 <- read_sf('result/spResult/07_ESRI54017_birds.gpkg')
birdsESRI54017Validated <- st_is_valid(birdsESRI54017)
birdsESRI54017Selected <- birdsESRI54017 %>%
    dplyr::select(sisid, sci_name, presence, origin, seasonal, Shape_Length, Shape_Area)
birdsESRI54017ValidatedRows <- birdsESRI54017Selected[birdsESRI54017Validated,]
# birdsESRI54017InvalidatedRows <- birdsESRI54017Selected[!birdsESRI54017Validated,]
# plot(birdsESRI54017InvalidatedRows$geom, axes=T)

shapeFile <- "result/gridsShapeResult/edges280.shp"
grid <- read_sf(shapeFile)
overlaps <- st_intersection(grid, birdsESRI54017ValidatedRows) %>% 
    st_collection_extract(type = c("POLYGON"))
# Warning message:
#     attribute variables are assumed to be spatially constant throughout all geometries 
overlapsSubset <- overlaps %>% 
    st_drop_geometry() %>% 
    dplyr::filter(
        presence %in% c(1, 2, 3),
        !origin %in% c(3, 4)) %>% 
    distinct(gridNo, sci_name, .keep_all = TRUE)

overlapsSubset$geometry <- NULL
outputFile = paste0('result/spResult/07_birdPerGrid280.csv')
write_csv(overlapsSubset, outputFile)
SRgrids <- overlapsSubset %>% group_by(gridNo) %>% summarise(SR=n())
SRgrids_shp <- grid %>% left_join(SRgrids, by = 'gridNo')

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
    ggtitle('Global Bird Species Richness')
ggsave('result/07_birdSReveryGrid280.png')

ggplot(SRgrids, aes(x = log10(SR))) + 
    geom_histogram(bins = 300, alpha = 0.7) + 
    xlab('log10(species richness)') +
    ylab('# grids') + 
    theme_light() + 
    theme(axis.title = element_text(size = 20))
ggsave('result/07_birdSReveryGrid280_hist.png')



GDfile <- 'result/genResult/05_piEveryGrid280.txt'
GDgrids <- read_delim(GDfile, delim = '\t') %>%
    filter(GD < 0.4)
env_GD_shp <- SRgrids_shp %>%
    left_join(GDgrids, by = 'gridNo') %>%
    filter(!is.na(GD))
# filename = paste0('result/envResult/shapeFiles/08_landscape_GD.shp')
filename = paste0('result/envResult/shapeFiles/07_birdGDGrid280.shp')
st_write(env_GD_shp, filename, append=FALSE)

p2 <- ggplot() +
    geom_sf(data = world_map, fill = "transparent", color = "black") +
    geom_sf(data = env_GD_shp, aes(fill = SR)) +
    theme_minimal() +
    scale_fill_distiller(name = 'Species Richness',
                         palette = "YlGnBu",
                         breaks = pretty_breaks(),
                         direction = 1) +
    theme(legend.position = 'bottom',
          plot.title = element_text(size = 18)) +
    guides(fill = guide_colorbar(title.position = "top",
                                 title.vjust = 0.5,
                                 title.hjust = 0.5,
                                 barwidth = 20, barheight = 1))  +
    ggtitle('Bird Species Richness')

filename = paste0('result/07_birdGDgrid280.png')
p1/p2
ggsave(filename, width = 8.96, height = 9.72)



# gdb_dir <- "data/spData/birds/BOTW_2023_1/BOTW.gdb/"
# gdb_layers <- st_layers(gdb_dir)
# gdb_layers
# layer_name <- gdb_layers$name[1]
# birdsRange <- st_read(dsn = gdb_dir, layer = layer_name)
# birdsESRI54017 <- st_transform(birdsRange, "ESRI:54017")
# filename <- 'result/spResult/07_ESRI54017_birds.gpkg'
# write_sf(birdsESRI54017, filename)
# 
# birdsESRI54017 <- read_sf('result/spResult/07_ESRI54017_birds.gpkg')
# nrow(birdsESRI54017)
# birdsESRI54017Cast <- st_cast(birdsESRI54017, "MULTIPOLYGON")# No error or message


