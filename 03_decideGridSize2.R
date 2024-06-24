### created by Emma
### last modified: 2024-04-24
### decide the grid size according to:
## 1. the sequences per species in each grid; 
## 2. the standard deviation of seqeunces per species in each grid;
## 3. the proportion of species with 5 sequences in different grid;


### load libraries ----
library(tidyverse)
library(sf)
library(scales)
library(rnaturalearth)


### load gene data ----
geneMeta <- read_delim("result/genResult/02_geneMeta.txt", 
                       delim = "\t", col_names = TRUE) %>% 
    dplyr::select(ID, latitude, longitude, gene, species, class, order, family, 
                  accession, source_id, basis_of_record, 
                  coordinate_uncertainty_in_meters, issue, flag) # 1971 species, 28407 records

# filter some records
geneMetaSelected <- geneMeta %>% 
    dplyr::select(ID, gene, latitude, longitude, species, class, issue, flag) %>% 
    # filter(!class %in% c('Sphagnopsida', 'Charophyceae', 'Chlorophyceae', 'Ulvophyceae')) %>% 
    # filter(class %in% c('Cycadopsida', 'Liliopsida', 'Gnetopsida',
    #                     'Magnoliopsida', 'Pinopsida')) %>% 
    filter(!grepl('INVALID', issue)) %>% 
    filter(flag != 'g') # 1571 species, 10867 records

# get the gene with most abundant records
gene_counts <- geneMetaSelected %>% 
    group_by(gene) %>% 
    summarise(n = n()) %>% 
    arrange(desc(n))  # 10 different genes; COI 10650 records; 

gene_with_highest_count <- gene_counts$gene[1]
geneMetaUsed <- geneMetaSelected %>% 
    filter(gene == gene_with_highest_count) # 1561 species; 10650 records

# change to sf of the genes records
geneCoor <- geneMetaUsed %>% dplyr::select(ID, longitude, latitude)
geneMetaUsed.sf <- st_as_sf(geneMetaUsed, coords=c("longitude", "latitude"))
st_crs(geneMetaUsed.sf) <- 4326
geneMetaUsed.sf <- st_transform(geneMetaUsed.sf, "ESRI:54017")

### load grids data ----
input_dir <- "result/gridsShapeResult"
shapeFiles <- list.files(input_dir, pattern = "\\.shp$", full.names = TRUE)
gridsLists <- list()
for (shp in shapeFiles) {
    gridsLists[[basename(shp)]] <- st_read(shp)
}

### function to assign gene sequence to each grid ----
gridGenesDistribution <- function(grid, geneMetaUsed.sf){
    # grid <- gridsLists$edges675.sh
    # grid4326 <- st_transform(grid, "EPSG:4326")
    # gridGenes <- st_join(grid4326, geneMetaUsed.sf)
    # gridGenes <- st_join(grid, geneMetaUsed.sf)
    gridGenes <- st_join(geneMetaUsed.sf, grid)
    
    gridGenesDFsubset <- st_drop_geometry(gridGenes) %>% 
        group_by(gridNo, species) %>% 
        summarize(SeqNum = n()) %>% 
        filter(SeqNum > 4) %>% 
        ungroup() %>% 
        left_join(st_drop_geometry(gridGenes))
    
    grid_SeqPerSp <- gridGenesDFsubset %>% 
        group_by(gridNo) %>%
        summarize(unique_species = n_distinct(species),
                  seqNum = n_distinct(ID)) %>% 
        mutate(seqPerSp = seqNum / unique_species)
    
    grid_SeqNumSD <- gridGenesDFsubset %>% 
        group_by(gridNo, species) %>%
        summarize(seqNum = n_distinct(ID)) %>% 
        ungroup() %>% 
        group_by(gridNo) %>% 
        summarize(SeqNumSD = sd(seqNum)) 
    
    gridSeq <- grid_SeqPerSp %>% 
        left_join(grid_SeqNumSD) %>% 
        filter(!is.na(SeqNumSD)) # sd is NA when there is only one species in the grid
    
    return(list(gridSeq = gridSeq, 
                gridGenesDFsubset = gridGenesDFsubset))
}

### apply the function and plot results ----
gridGenesDistributionList <- map(gridsLists, gridGenesDistribution, geneMetaUsed.sf)

first_elements <- map(gridGenesDistributionList, pluck, 1)
merged_df <- bind_rows(first_elements, .id = "original_dataframe") %>% 
    mutate(edgeLen = as.numeric(gsub("\\D", "", original_dataframe))) %>% 
    arrange(edgeLen) %>% 
    mutate(edgeLen = factor(edgeLen)) %>%
    pivot_longer(cols = c(seqPerSp, SeqNumSD), 
                 names_to = "statType", 
                 values_to = "stat")
merged_df$statType <- factor(merged_df$statType, levels = c("seqPerSp", "SeqNumSD"))

ggplot(merged_df) +
    geom_boxplot(aes(x = edgeLen, y = stat), outlier.size = 0.5) + 
    xlab('Grid Size (km)') + ylab('') +
    theme(strip.text = element_text(size = 14, face = "bold")) +
    facet_wrap(~statType, dir = "v", 
               labeller = labeller(
                   statType = c("seqPerSp" = "Sequence per Species per Grid", 
                                "SeqNumSD" = "Sequence Number Standard Deviation per Grid")),
               scales = "free_y") + 
    theme_light() +
    theme(strip.text = element_text(size = 14, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1))  +
    geom_hline(data = subset(merged_df, statType == "seqPerSp"),
               aes(yintercept = 6), color = "blue") +
    geom_hline(data = subset(merged_df, statType == "seqPerSp"),
               aes(yintercept = 7), color = "red") +
    geom_hline(data = subset(merged_df, statType == "SeqNumSD"),
               aes(yintercept = 1), color = "red") +
    geom_hline(data = subset(merged_df, statType == "SeqNumSD"),
               aes(yintercept = 2), color = "blue")
    
ggsave('result/03_selectGridSize1_5sequences.png', 
       width = 10, height = 6)

second_elements <- map(gridGenesDistributionList, pluck, 2)

spPerGrids <- function(gridGenesDFsubset, spNum){
    # spNum <- 3

    # gridGenesDFsubset contains species with at least 5 sequences on grid-based
    spNum_goodGrids <- gridGenesDFsubset %>%
        group_by(gridNo) %>%
        summarise(sp_n = n_distinct(species)) %>%
        filter(sp_n >= spNum) %>% 
        ungroup() %>% 
        left_join(gridGenesDFsubset)
    
    uniqueGridNum <- length(unique(gridGenesDFsubset$gridNo))
    uniqueGridSps <- length(unique(gridGenesDFsubset$species))
    goodGridNum <- length(unique(spNum_goodGrids$gridNo))
    goodGridSps <- length(unique(spNum_goodGrids$species))
    
    
    return(c(uniqueGridNum, uniqueGridSps,
             goodGridNum, goodGridSps))
}

gridSpDistributionList <- map(second_elements, spPerGrids, 3)
merged_df2 <- bind_rows(gridSpDistributionList, .id = "original_dataframe") 
rownames(merged_df2) <- c('totalGrids', 'totalSpecies', 'goodGrids', 'goodSpecies')
merged_df3 <- as.data.frame(t(merged_df2)) %>% 
    rownames_to_column('edgeLen') %>% 
    mutate(edgeLen = as.numeric(gsub("\\D", "", edgeLen))) %>% 
    arrange(edgeLen) %>% 
    mutate(edgeLen = factor(edgeLen))
               
ggplot(merged_df3) +
    geom_bar(mapping = aes(x = edgeLen, y = totalGrids), 
             stat = "identity", fill = "grey", alpha = 0.5) + 
    geom_bar(mapping = aes(x = edgeLen, y = goodGrids), 
             stat = "identity", fill = "lightgreen") +
    geom_line(mapping = aes(x = edgeLen, y = totalSpecies *(1/2)), 
              linewidth = 2, color = "grey", group = 1) +
    geom_line(mapping = aes(x = edgeLen, y = goodSpecies *(1/2)), 
              linewidth = 2, color = "lightgreen", group = 1) +
    geom_hline(aes(yintercept = 30), color = "yellow", linewidth = 1.2) +
    scale_y_continuous(name = "Numbers of grids", 
                       breaks = c(0, 30, 50, 100, 150),
                       sec.axis = sec_axis(~. * 2, name = "Numbers of species")) + 
    xlab('Grid Size (km)') +
    theme_light() +
    theme(strip.text = element_text(size = 14, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title = element_text(size = 12))
ggsave('result/03_selectGridSize2_5sequences.png',
       width = 10, height = 6)


### select sequences ----
shapeFile <- "result/gridsShapeResult/edges10.shp"
grid <- st_read(shapeFile)
gridGenes <- st_join(geneMetaUsed.sf, grid)
geneCoor <- geneMetaUsed %>% dplyr::select(ID, longitude, latitude)
df_used_group <- st_drop_geometry(gridGenes) %>% 
    group_by(gridNo, species) %>%
    summarize(SeqNum = n()) %>% 
    filter(SeqNum > 4) %>% 
    ungroup() %>% 
    left_join(st_drop_geometry(gridGenes)) %>% 
    left_join(geneCoor) %>% 
    dplyr::select(ID, gene, species, issue, flag, gridNo, longitude, latitude)

outputFile = paste0('result/genResult/03_selectedRecords_', gene_with_highest_count, '_Grid10.csv')
write_csv(df_used_group, outputFile)

spNum <- 3
df_used_group2 <- st_drop_geometry(gridGenes) %>% 
    group_by(gridNo, species) %>%
    summarize(SeqNum = n()) %>% 
    filter(SeqNum > 4) %>% 
    ungroup() %>% 
    group_by(gridNo) %>% 
    summarise(speciesNum = n_distinct(species)) %>% 
    filter(speciesNum >= spNum) %>% 
    left_join(st_drop_geometry(gridGenes)) %>% 
    left_join(geneCoor) %>% 
    dplyr::select(ID, gene, species, issue, flag, gridNo, longitude, latitude)
outputFile = paste0('result/genResult/03_selectedRecords_', gene_with_highest_count, '_GoodGrid10.csv')
write_csv(df_used_group2, outputFile)

### plot species distribution in each plot ----
a <- df_used_group %>% # df_used_group each species has at least 4 seq on grid-based
    dplyr::select(gridNo, species) %>% 
    group_by(gridNo) %>% 
    summarise(n = n_distinct(species)) %>% 
    ungroup() %>% 
    arrange(desc(n))
frequency_table <- table(a$n)
png("result/03_SpPerPlotGrid10.png", width = 2000, height = 1250, res = 300)
plot(frequency_table, xlab = '# Species per grid', ylab = '# Grids')
title('Most grids have only one species')
dev.off()

b <- df_used_group %>%    # df_used_group each species has at least 4 seq on grid-based
    select(species, gridNo) %>% 
    distinct() %>% 
    group_by(species) %>% 
    summarise(occurrence = n()) %>%
    ungroup() %>% 
    arrange(desc(occurrence))
frequency_table <- table(b$occurrence)
png("result/03_PlotPerSpGrid10.png", width = 2000, height = 1250, res = 300)
plot(frequency_table, xlab = '# Grids', ylab = '# Species')
title('Most species occur only in one grid')
dev.off()

# highest: grid4763, grid6850, grid7111
top10 <- a[1:10,]
c <- grid %>% filter(gridNo %in% head(a$gridNo, 10))
order_df1 <- match(top10$gridNo, c$gridNo)
c_reordered <- c[order_df1, ]
c_reordered$top <- factor(paste0('top', seq(nrow(c))), 
                          levels = paste0('top', seq(nrow(c))))

world_map <- ne_countries(scale = "medium", returnclass = "sf")
world_map <- st_transform(world_map, "ESRI:54017")
ggplot() + 
    geom_sf(data = world_map, fill = "transparent", color = "black") +
    # geom_sf(data = c_reordered, aes(fill=top)) +
    geom_sf_text(data = c_reordered, aes(label=top)) +
    theme_minimal() +
    ggtitle('Species hotspots of from phylogatR') +
    labs(fill = "Hotspot",
         subtitle = "         Grid size 10 * 10 km^2") +
    theme(plot.subtitle = element_text(hjust = 0, size = 12),
          plot.title = element_text(size = 24))
ggsave("result/03_selectGridSize3_hotspotGrid10.png")


####### looked into plots
# z <- df_used_group %>% filter(gridNo=='grid2419') %>% 
#     dplyr::select(ID, gene, species, gridNo) %>% 
#     left_join(geneMeta)
# View(z)
