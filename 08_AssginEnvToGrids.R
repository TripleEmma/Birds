### created by Emma
### last modified: 2024-04-26
### Assign environment factor value to grids

### load libraries ----
library(sf)
library(terra)
library(tidyverse)
library(scales)
library(rnaturalearth)
library(patchwork)
library(geodata)
library(spatialEco)

rm(list = ls())
namesToEnv <- c('solar', 'preci', 'temp', 'DiurnalR',
                'soilPH', 'soilN', 'footprint', 'population',
                'simpsonInd', 'pet', 'soilP')

source('08_FunctionAssginEnvToGrids.R')

### load the grids ----
shapeFile <- "result/gridsShapeResult/edges280.shp"
grid <- read_sf(shapeFile)
# grid4326 <- st_transform(grid, 'EPSG:4326')
GDfile <- 'result/genResult/05_piEveryGrid280.txt'

### load solar data ----
solarFiles <- paste0("data/envData/climate/wc2.1_5m_srad/wc2.1_5m_srad_", sprintf("%02d", 1:12), ".tif")
solarRaster <- mean(rast(solarFiles), na.rm = TRUE)
assign_to_grid_cell(solarRaster, 'solar', grid, GDfile)

### load precipitation, temperature, DiurnalRange, data ----
climates <- vector(mode = 'list', length = 3)
climates[[1]] <- rast('data/envData/climate/worldclim/wc2.1_5m/wc2.1_5m_bio_1.tif') 
climates[[2]] <- rast('data/envData/climate/worldclim/wc2.1_5m/wc2.1_5m_bio_2.tif')
climates[[3]] <- rast('data/envData/climate/worldclim/wc2.1_5m/wc2.1_5m_bio_12.tif')
# 5 minutes equals to roughly 10km at the equator
climatesType <- c('temp', 'DiurnalR', 'preci')
walk2(climates, climatesType, assign_to_grid_cell, grid, GDfile)

### load soil data ----
# extend is different
# extendExt <- ext(-180, 180, -90, 90)
# envRaster <- extend(envRaster, extendExt)
soils <- vector(mode = 'list', length = 3)
soils[[1]] <- rast('data/envData/soil/soil_world/nitrogen_0-5cm_mean_30s.tif')
soils[[2]] <- rast('data/envData/soil/soil_world/phh2o_0-5cm_mean_30s.tif')
soils[[3]] <- rast('data/envData/soil/Global_distribution_of_soil_phosphorus_retention_potential/Global_distribution_of_soil_phosphorus_retention_potential.tif')
soilsType <- c('soilN', 'soilPH', 'soilP')
walk2(soils, soilsType, assign_to_grid_cell, grid, GDfile)

### load disturbance data ----
disturbance <- vector(mode = 'list', length = 2)
disturbance[[1]] <- rast('data/envData/disturbance/wildareas-v3-2009-human-footprint_geo.tif')
disturbance[[2]] <- rast('data/envData/disturbance/pop/gpw_v4_population_density_rev11_2010_10m.tif')
disturbanceType <- c('footprint', 'population')
walk2(disturbance, disturbanceType, assign_to_grid_cell, grid, GDfile)


### load potential evapotranspiration data ----
petRaster <- rast('data/envData/Global-AI_ET0_v3_annual/et0_v3_yr.tif')
assign_to_grid_cell(petRaster, 'pet', grid, GDfile)


### load landscape data ----
coverDFs <- function(cover, grid){
    print(cover)
    # cover <- "trees"
    envRaster <- landcover(var=cover, path='data/envData/landCover/')
    envRasterESRI54017 <- terra::project(envRaster, 'ESRI:54017')
    # for categorical data, use default method (nearest neighbor method) for the reprojection.
    grid_ras <- rasterize(vect(grid), envRasterESRI54017, field = 'gridNo')
    grid_env <- zonal(envRasterESRI54017, 
                      grid_ras, 
                      fun = "mean", 
                      na.rm = TRUE)
    grid_env
}

covers <- list("trees", "grassland", "shrubs", 
              "cropland", "built", "bare", "snow", 
              "water", "wetland", "mangroves", "moss")

coverDFs <- map(covers, coverDFs, grid)
names(coverDFs) <- unlist(covers)
coverMerged <- reduce(coverDFs, left_join)
simpson <- function(row) {
    if (all(is.na(row))) {
        return(NA)
    }
    1 - sum(row^2, na.rm = TRUE)
}

coverSimpsonDiv <- coverMerged %>%
    rowwise() %>%
    mutate(simpsonInd = simpson(c_across(where(is.numeric)))) %>%
    unnest(simpsonInd)
grid_env_shp <- grid %>%
    left_join(coverSimpsonDiv, by = 'gridNo') %>%
    filter(!is.na(simpsonInd))

filename = paste0('result/envResult/shapeFiles/08_landscape.shp')
write_sf(grid_env_shp, filename, append=FALSE)
grid_env_shp <- read_sf('result/envResult/shapeFiles/08_landscape.shp')
world_map <- ne_countries(scale = "medium", returnclass = "sf")
world_map <- st_transform(world_map, "ESRI:54017")

p1 <- ggplot() +
    geom_sf(data = world_map, fill = "transparent", color = "black") +
    geom_sf(data = grid_env_shp, aes(fill = simpsonInd)) +
    theme_minimal() +
    scale_fill_distiller(name = legends[["simpsonInd"]],
                         palette = 'GnBu',
                         breaks = pretty_breaks(),
                         direction = 1) +
    theme(legend.position = 'bottom',
          plot.title = element_text(size = 18)) +
    guides(fill = guide_colorbar(title.position = "top",
                                 title.vjust = 0.5,
                                 title.hjust = 0.5,
                                 barwidth = 20, barheight = 1)) +
    ggtitle(titles[["simpsonInd"]])

# only grids with GD
GDgrids <- read_delim(GDfile, delim = '\t') %>%
    filter(GD < 0.4)
env_GD_shp <- grid_env_shp %>%
    left_join(GDgrids, by = 'gridNo') %>%
    filter(!is.na(GD))
# filename = paste0('result/envResult/shapeFiles/08_landscape_GD.shp')
filename = paste0('result/envResult/shapeFiles/08_landscape_GDGrid280.shp')
st_write(env_GD_shp, filename, append=FALSE)

p2 <- ggplot() +
    geom_sf(data = world_map, fill = "transparent", color = "black") +
    geom_sf(data = env_GD_shp, aes(fill = simpsonInd)) +
    theme_minimal() +
    scale_fill_distiller(name = legends[["simpsonInd"]],
                         palette = 'GnBu',
                         breaks = pretty_breaks(),
                         direction = 1) +
    theme(legend.position = 'bottom',
          plot.title = element_text(size = 18)) +
    guides(fill = guide_colorbar(title.position = "top",
                                 title.vjust = 0.5,
                                 title.hjust = 0.5,
                                 barwidth = 20, barheight = 1))  +
    ggtitle(titles2[["simpsonInd"]])

filename = paste0('result/08_', namesToEnvFull[["simpsonInd"]], '_Grid280.png')
p1/p2
ggsave(filename, width = 8.96, height = 9.72)



elevation <- elevation_global(0.5, 'data/envData/elevation/')
# trying URL 'https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_30s_elev.zip'
# Content type 'application/zip' length 338188707 bytes (322.5 MB)
# r.tri <- spatialEco::tri(elevation)
envRasterESRI54017 <- terra::project(elevation, 'ESRI:54017', method = "mode")
names(envRasterESRI54017) <- 'ruggedness'
r.tri <- spatialEco::tri(envRasterESRI54017)
grid_ras <- rasterize(vect(grid), r.tri, field = 'gridNo')
grid_env <- terra::zonal(r.tri , 
                         grid_ras, 
                         fun = "mean", 
                         na.rm = TRUE)
grid_env_shp <- grid %>%
    left_join(grid_env, by = 'gridNo') %>%
    filter(!is.na(ruggedness))
filename = paste0('result/envResult/shapeFiles/08_ruggedness.shp')
write_sf(grid_env_shp, filename, append=FALSE)

p1 <- ggplot() +
    geom_sf(data = world_map, fill = "transparent", color = "black") +
    geom_sf(data = grid_env_shp, aes(fill = ruggedness)) +
    theme_minimal() +
    scale_fill_distiller(name = 'Terrain Ruggedness Index',
                         palette = 'GnBu',
                         breaks = pretty_breaks(),
                         direction = 1) +
    theme(legend.position = 'bottom',
          plot.title = element_text(size = 18)) +
    guides(fill = guide_colorbar(title.position = "top",
                                 title.vjust = 0.5,
                                 title.hjust = 0.5,
                                 barwidth = 20, barheight = 1)) +
    ggtitle('Global Terrain Ruggedness Index')

GDgrids <- read_delim(GDfile, delim = '\t') %>%
    filter(GD < 0.4)
env_GD_shp <- grid_env_shp %>%
    left_join(GDgrids, by = 'gridNo') %>%
    filter(!is.na(GD))
# filename = paste0('result/envResult/shapeFiles/08_landscape_GD.shp')
filename = paste0('result/envResult/shapeFiles/08_ruggedness_GDGrid280.shp')
st_write(env_GD_shp, filename, append=FALSE)

p2 <- ggplot() +
    geom_sf(data = world_map, fill = "transparent", color = "black") +
    geom_sf(data = env_GD_shp, aes(fill = simpsonInd)) +
    theme_minimal() +
    scale_fill_distiller(name = 'Terrain Ruggedness Index',
                         palette = 'GnBu',
                         breaks = pretty_breaks(),
                         direction = 1) +
    theme(legend.position = 'bottom',
          plot.title = element_text(size = 18)) +
    guides(fill = guide_colorbar(title.position = "top",
                                 title.vjust = 0.5,
                                 title.hjust = 0.5,
                                 barwidth = 20, barheight = 1))  +
    ggtitle('Terrain Ruggedness Index')

filename = paste0('result/08_ruggedness_Grid280.png')
p1/p2
ggsave(filename, width = 8.96, height = 9.72)
