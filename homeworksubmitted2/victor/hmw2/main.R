wd <- '~/projects/ubc_related/courses/bayes2025homework/homeworksubmitted/victor/hmw2'

library(hypervolume)
library(doFuture)
library(dplyr)
library(gpkg)
library(terra)

run_hypervolume <- FALSE
process_gpkg <- FALSE

# load input data
bodymass <- read.csv(file.path(wd, 'input', 'carnivorebodymass.csv'))
teeth <- read.csv(file.path(wd, 'input', 'carnivoreteeth.csv'))

# load supplementary tables of the original paper
table1 <- read.csv(file.path(wd, 'input', 'table_1.csv'))
colnames(table1) <- c('species_A', 'species_B', 'family', 'time', 
                      'rangesize_A', 'rangesize_B', 'overlap', 'mass_A', 'mass_B')
table2 <- read.csv(file.path(wd, 'input', 'table_2.csv'))
colnames(table2) <- c('species_A', 'species_B', 'M1_A', 'M1_B', 
                      'PM4_A', 'PM4_B', 'CsupL_A', 'CsupL_B')
data_supp <- inner_join(table1, table2, by = c('species_A', 'species_B'))

# a first idea is to compute an index of functionnal diversity, based on the tooth length
# this can be done by computing hypervolume similarity
teeth_filtered <- teeth %>%
  dplyr::select('Species', 'PM4', 'CsupL') %>%
  na.omit() %>%
  group_by(Species) %>%
  mutate(nobs = n()) %>%
  filter(log(nobs) >= 2) # min number of obs. to compute a 2D hypervolume 

# compute a 2-dimensional functional hypervolume per species
if(run_hypervolume){
  tb <- Sys.time()
  plan(multisession, workers = 18)
  y <- foreach(s = unique(teeth_filtered$Species)) %dofuture% {
    fun_hpv <- hypervolume(teeth_filtered[teeth_filtered$Species == s, c('PM4', 'CsupL')], method='gaussian')
    saveRDS(fun_hpv, file = file.path(wd, 'output/trait_hypervolume', paste0(s, '.rds')))
  }
  plan(sequential);gc()
  te <- Sys.time()
  print(te-tb) # less than 2 minutes on my computer
}

# for each species pair in the original paper, compute the similarity between the hypervolumes
data_supp <- data_supp %>%
  filter(species_A %in% teeth_filtered$Species & species_B %in% teeth_filtered$Species)
data_supp$fun_dissimilarity <- NA
for(i in 1:nrow(data_supp)){
  hpv_A <- readRDS(file.path(wd, 'output/trait_hypervolume', paste0(data_supp[i, 'species_A'], '.rds')))
  hpv_B <- readRDS(file.path(wd, 'output/trait_hypervolume', paste0(data_supp[i, 'species_B'], '.rds')))
  
  hpv_set <- hypervolume_set(hpv_A, hpv_B, check.memory=FALSE)
  hpv_stats <- hypervolume_overlap_statistics(hpv_set)
  
  data_supp[i, 'fun_dissimilarity'] <- 1-hpv_stats['sorensen']
}
plot(data_supp$overlap ~ data_supp$fun_dissimilarity)
# unfortunately, we don't have a lot of data

# let's try with more species! not only the closest pairs of the paper
# we thus need to compute the range overlap for these species combinations
if(process_gpkg){
  # data come from the Mammal Diversity Databse
  g <-  geopackage('/home/victor/Downloads/MDD/MDD_Carnivora.gpkg', connect = TRUE)
  mdd <- gpkg_vect(g, 'MDD_Carnivora')
  mdd_select <- mdd[mdd$sciname %in% unique(teeth_filtered$Species)]
  writeVector(mdd_select, filename =file.path(wd, 'output', 'carnivor_ranges.shp'))
  RSQLite::dbDisconnect(g) 
  rm(g, mdd, mdd_select)
  gc()
}else{
  mdd_ranges <- vect(file.path(wd, 'output', 'carnivor_ranges.shp'))
}

# unique pairs of species
species_comb <- t(combn(unique(teeth_filtered$Species), 2))

# create our own dataframe!
newdata <- data.frame(species_A = NA, species_B = NA, 
                      range_A = NA, range_B = NA,
                      range_overlap_area = NA, fun_similarity = NA)
# .mdd_ranges <- wrap(mdd_ranges) # needed for parallel computation
# plan(multisession, workers = 2)
for(i in 1:nrow(species_comb)){
  cat(paste0(i, '\n'))
  #mdd_ranges <- rast(.mdd_ranges)
  # compute range metrics
  range1 <-  mdd_ranges[mdd_ranges$sciname %in% species_comb[i,1]]
  range2 <-  mdd_ranges[mdd_ranges$sciname %in% species_comb[i,2]]
  int12 <- terra::intersect(range1, range2)
  newdata[i, 'species_A'] <- species_comb[i,1]
  newdata[i, 'species_B'] <- species_comb[i,2]
  range_A <- expanse(range1, unit = 'km')
  newdata[i, 'range_A'] <- ifelse(length(range_A) == 0, NA, range_A)
  range_B <- expanse(range2, unit = 'km')
  newdata[i, 'range_B'] <- ifelse(length(range_B) == 0, NA, range_B)
  int_area <- expanse(int12, unit = 'km')
  newdata[i, 'range_overlap_area'] <- ifelse(length(int_area) == 0, 0, int_area)
  rm(range1, range2, int12);gc()
  
  # compute functional dissimilarity metrics
  hpv1 <- readRDS(file.path(wd, 'output/trait_hypervolume', paste0(species_comb[i,1], '.rds')))
  hpv2 <- readRDS(file.path(wd, 'output/trait_hypervolume', paste0(species_comb[i,2], '.rds')))
  hpv_set <- hypervolume_set(hpv1, hpv2, check.memory=FALSE)
  hpv_stats <- hypervolume_overlap_statistics(hpv_set)
  newdata[i, 'fun_similarity'] <- hpv_stats['sorensen']
}
saveRDS(newdata, file.path(wd, 'output', 'newdata.rds'))

newdata2 <- na.omit(newdata) %>% filter(fun_similarity != 0)
newdata2$spat_similarity <- 2*newdata2$range_overlap_area/(newdata2$range_A+newdata2$range_B)
newdata2 <- newdata2 %>% 
  rowwise() %>% 
  mutate(range_overlap_jdavies = range_overlap_area/max(c(range_A, range_B)),
         fun_dissimilarity = 1-fun_similarity)
plot(newdata2$fun_dissimilarity ~ (newdata2$spat_similarity))

plot(newdata2$spat_similarity ~ newdata2$fun_similarity)
plot(newdata2$range_overlap_jdavies ~ newdata2$fun_dissimilarity)

