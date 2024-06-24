### created by Emma
### last modified: 2024-04-25
### plot genetic diveristy pi 

### load libraries ----
library(sf)
library(tidyverse)
library(scales)
library(rnaturalearth)
library(patchwork)

### load grids data ----
shapeFile <- "result/gridsShapeResult/edges280.shp"
grid <- st_read(shapeFile)

### load genetic diversity ----
GDgrids <- read_delim('result/genResult/05_piEveryGrid10.txt', delim = '\t')

tooDiverse <- GDgrids %>% filter(GD >= 0.4) %>% pull(gridNo)
b <- grid %>% filter(gridNo %in% tooDiverse)
world_map <- ne_countries(scale = "medium", returnclass = "sf")
world_map <- st_transform(world_map, "ESRI:54017")
ggplot() + 
    geom_sf(data = world_map, fill = "transparent", color = "black") +
    geom_sf(data = b, fill = 'red') +
    theme_minimal() +
    ggtitle('pi larger than 0.4') 
ggsave('result/06_tooGDdiverseGrid10.png')

### plot genetic diversity ----
GDgrids <- GDgrids %>% filter(GD < 0.4) %>% arrange(GD) # 125 plots
GDgrids_shp <- grid %>% left_join(GDgrids, by = 'gridNo') %>% 
    # filter(gridNo %in% goodGrids) %>% 
    filter(!is.na(GD))

# Load world map data
world_map <- ne_countries(scale = "medium", returnclass = "sf")
world_map <- st_transform(world_map, "ESRI:54017")
# plot
plotTitle <- paste0('Global genetic diversity (', nrow(GDgrids_shp), ' grids)')
p1 <- ggplot() + 
    geom_sf(data = world_map, fill = "transparent", color = "black") +
    geom_sf(data = GDgrids_shp, aes(fill = GD)) +
    theme_minimal() +
    scale_fill_distiller(name = expression(paste("Genetic Diversity ", pi)),
                         palette = "YlGnBu",
                         breaks = pretty_breaks(),
                         direction = 1) +
    theme(legend.position = 'bottom',
          plot.title = element_text(size = 18)) +
    guides(fill = guide_colorbar(barwidth = 20, barheight = 1)) +
    ggtitle(plotTitle)
ggsave('result/06_GDeveryGrid280All.png')


inputFile <- 'result/genResult/03_selectedRecords_COI_GoodGrid10.csv'
goodGrids <- read_csv(inputFile) %>% 
    dplyr::select(gridNo) %>% 
    distinct() %>% 
    pull()

goodGridsGD_shp <- GDgrids_shp %>% filter(gridNo %in% goodGrids)

allGrids <- paste0('All Grids (', nrow(GDgrids_shp), ' grids)')
goodGrids <- paste0('Grids with at least 3 species (', nrow(goodGridsGD_shp), ' grids)')
GDgrids_shp$Type <- allGrids
goodGridsGD_shp$Type <- goodGrids

# Combine the data frames
allTogether <- rbind(GDgrids_shp, goodGridsGD_shp)

# Convert the type column to a factor with specified levels
allTogether$Type <- factor(allTogether$Type, levels = c(allGrids, goodGrids))

ggplot() + 
    geom_sf(data = world_map, fill = "transparent", color = "black") +
    geom_sf(data = allTogether, aes(fill = GD)) +
    theme_minimal() +
    scale_fill_distiller(name = expression(paste("Genetic Diversity ", pi)),
                         palette = "YlGnBu",
                         breaks = pretty_breaks(),
                         direction = 1) +
    theme(legend.position = 'bottom',
          plot.title = element_text(size = 18)) +
    guides(fill = guide_colorbar(barwidth = 20, barheight = 1)) +
    facet_wrap(~Type, dir = "v")

ggsave('result/06_GDeveryGrid280Both.png')


inputFile <- 'result/genResult/03_selectedRecords_COI_Grid10.csv'
df <- read_csv(inputFile) %>% 
    group_by(gridNo) %>% 
    summarise(spNums = n_distinct(species)) %>% 
    ungroup() %>% 
    left_join(GDgrids)

dfMean <- df %>% 
    group_by(spNums) %>% 
    summarise(GDavg = mean(GD)) %>% 
    ungroup()
    
ggplot_method <- expression(paste(italic('`geom_smooth()` using method = \'loess\' and formula = \'y ~ x\'')))  
ggplot(df, aes(x= sqrt(spNums), y = GD)) + 
    geom_point() + 
    geom_point(dfMean, mapping = aes(x=sqrt(spNums), y = GDavg), color = 'red') +
    geom_line(dfMean, mapping = aes(x=sqrt(spNums), y = GDavg), color = 'red') +
    geom_smooth() +
    ggtitle('Grid size 10km * 10km') +
    # labs(caption = ggplot_method) +
    ylab(expression(paste("Genetic Diversity ", pi))) + 
    xlab(expression(sqrt("# Species"))) # + 
    # theme(plot.caption = element_text(size = 13.6, colour = 'grey'))
ggsave('result/06_piConvergeGrid10.png')    
