## add dummy variable within/between basins

library(sp)
library(raster)
library(ggplot2)
library(dplyr)
library(tidyverse)

getwd()

## directory for figures
out.dir <- "/Users/katieirving/OneDrive - SCCWRP/Documents - Katie’s MacBook Pro/git/sYNGEO_2023/Figures/"

## upload fish abundance and site data
originaldata <- read.csv("input_data/Bio/fishdata_selection_basins_same_time_window_10262020.csv")
head(originaldata)

## take only sites and basins
sites <- originaldata %>%
  dplyr::select(SiteID, HydroBasin) %>%
  distinct()

## all site synchrony

sync <- read.csv("sync/02_between_all_sites_single_traits_CWM_CWV_CV.csv")
head(sync)

## join first site basin
site1 <- left_join(sync, sites, by = c("Site_ID1" = "SiteID"))
head(site1)

## change names so we know it's for site 1
site1 <- site1 %>%
  dplyr::select(-X) %>%
  rename(HydroBasin1 = HydroBasin)

## join site 2

site2 <- left_join(site1, sites, by = c("Site_ID2" = "SiteID"))
head(site2)

## change names so we know it's for site 2
site2 <- site2 %>%
  # select(-X) %>%
  rename(HydroBasin2 = HydroBasin)

head(site2)


## create dummy variable. If hydrobasin 1 matches hydrobasin 2 = within site (1), if not then between sites (0)

con <- site2 %>%
  mutate(Connectivity = ifelse(HydroBasin1 == HydroBasin2, 1, 0))

head(con)

write.csv(con, "sync/03_func_sync_connectivity.csv")

rm(con)
# calculate distance ------------------------------------------------------

sync <- read.csv("sync/03_func_sync_connectivity.csv")
head(sync)

## define pairs
syncsites <- sync %>%
  dplyr::select(Site_ID1, Site_ID2, Pair, Connectivity) %>%
  # filter(TraitGroup == "FoodAquisition") %>%
  # tidyr::pivot_wider(names_from = Trait, values_from = Correlation) %>%
  # mutate(Pair = paste(Site_ID1, ".", Site_ID2, sep="")) %>%
  mutate(Euclid_Dist_Meters = 0, Similarity = 0, MeanLat = 0, MeanLon = 0) %>%
  # dplyr::select(-FeedingGroup, ) %>%
  distinct()

head(syncsites)

pairs <- unique(syncsites$Pair)
tail(pairs)
length(pairs) ##  207299

## site coords

fish_ab <- read.csv("input_data/Bio/fishdata_selection_basins_same_time_window_10262020.csv")
head(fish_ab)

## get coords
SiteCoords <- fish_ab %>%
  dplyr::select(SiteID, Latitude, Longitude) %>%
  distinct()

## loop over pairs - takes a long ass time, go do something else for a bit...

for(p in 1:length(pairs)) {

  ## get pair from sync data

  pairx <- syncsites %>%
    filter(Pair == pairs[p])
  # pairx
  ## define sites
  S1 <- pairx$Site_ID1
  S2 <- pairx$Site_ID2

  # S1

  ## get coords for each site
  CoordsS1 <- SiteCoords %>%
    filter(SiteID == S1) %>%
    dplyr::select(Longitude, Latitude, SiteID)

  CoordsS2 <- SiteCoords %>%
    filter(SiteID == S2) %>%
    dplyr::select(Longitude, Latitude, SiteID)

  sp::coordinates(CoordsS1) <- c("Longitude", "Latitude")
  sp::coordinates(CoordsS2) <- c("Longitude", "Latitude")

  #Make a distance matrix
  dst <- pointDistance(CoordsS1,CoordsS2, lonlat=TRUE)
  # str(dst)
  # get mean latitude/longitude
  MeanLat <- (CoordsS1$Latitude+CoordsS2$Latitude)/2
  MeanLon <- (CoordsS1$Longitude+CoordsS2$Longitude)/2

  ## add to dataframe
  syncsites[p,5] <- dst
  syncsites[p,7] <- MeanLat
  syncsites[p,8] <- MeanLon
 # head(syncsites)

}

head(syncsites)

save(syncsites, file = "sync/03_all_pair_distances.RData")



# Combine with single trait DF --------------------------------------------

sync <- read.csv("sync/02_between_all_sites_single_traits_CWM_CWV_CV.csv")
head(sync)

load(file = "sync/03_all_pair_distances.RData") ## syncsites

## make Pair column
# sync <- sync %>%
#   unite(Pair, Site_ID1:Site_ID2,sep = ".", remove=F)

## take only distance columns
sync_sub <- syncsites %>%
  dplyr::select(Pair:MeanLon)

## join
all_sync <- left_join(sync, sync_sub, by = "Pair")
head(all_sync)
## convert to similarities

syncDF <- all_sync %>%
  group_by(Region, Trait) %>%
  mutate(MaxDist = max(Euclid_Dist_Meters)) %>%
  mutate(Similarity = 1-(Euclid_Dist_Meters/MaxDist))

head(syncDF)


save(syncDF, file = "sync/03_sync_traits_CWM_CWV_distances.RData")


# Join with temp and flow -------------------------------------------------

sync <- read.csv( "sync/02_between_all_sites_temp_synchrony.csv")
head(sync)
dim(sync)

load(file = "sync/03_all_pair_distances.RData") ## syncsites

## make Pair column
sync <- sync %>%
  dplyr::select(-X) %>%
  unite(Pair, Site_ID1:Site_ID2,sep = ".", remove=F)

## take only distance columns
sync_sub <- syncsites %>%
  dplyr::select(Pair:MeanLon)

## join
all_sync <- left_join(sync, sync_sub, by = "Pair")
head(all_sync)

## convert to similarities

syncDF <- all_sync %>%
  group_by(Region, Metric) %>%
  mutate(MaxDist = max(Euclid_Dist_Meters)) %>%
  mutate(Similarity = 1-(Euclid_Dist_Meters/MaxDist))

head(syncDF)


save(syncDF, file = "sync/03_sync_temp_distances.RData")

