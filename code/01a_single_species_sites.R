### removing sites with only 1 species

library(tidylog)
library(tidyverse)


## upload abundances
load(file="sync/01_fish_abundances.RData")

head(fish_ab_rel_int)
names(fish_ab_rel_int)


## tally number of species per site

fishTally <- fish_ab_rel_int %>%
  select(Species, SiteID, BioRealm, Country) %>%
  distinct() %>%
  group_by(Country, SiteID) %>% 
  tally() %>%
  arrange(Country, SiteID) 

fishTally

## list single species sites

singleSp <- fishTally %>%
  filter(n==1) #%>%
  distinct(SiteID)

singleSpl <- unique(singleSp$SiteID)
length(singleSpl) ## 49

## sites with 2 species
singleSp2 <- fishTally %>%
  filter(n<=2) 

singleSp2l <- unique(singleSp2$SiteID)
length(singleSp2l) ## 103


## sites with 2 species
singleSp3 <- fishTally %>%
  filter(n<=3) 

singleSp3l <- unique(singleSp3$SiteID)
length(singleSp3l) ## 151
singleSp3$SiteID

test <- fish_ab_rel_int %>%
  filter(SiteID %in% singleSp3$SiteID)

test2 <- fish_ab_rel_int %>%
  filter(!SiteID %in% singleSp3$SiteID)

## save out

save(singleSpl, file="output_data/01a_sites_1_species.RData")
save(singleSp2l, file="output_data/01a_sites_2_species.RData")
save(singleSp3l, file="output_data/01a_sites_3_species.RData")


# Get sites and coords ----------------------------------------------------

fish_ab <- read.csv("input_data/Bio/fishdata_selection_basins_same_time_window_10262020.csv")

head(fish_ab)

## remove basins - Sweden 2080030650 and 2080031160
## keep only origin Ohio and LTRM in mississippi
## change to relative abundance
## missing trait values - remove fish with less than 2 traits (check that it's less than 5%)
## remove basins - Sweden 2080030650 and 2080031160
## keep only origin Ohio and LTRM in mississippi
basins_remove <- c(2080031160, 2080030650)
origin_remove <- c("Tennessee", "ltr")

fish_ab <- fish_ab %>%
  filter(!HydroBasin %in% basins_remove, !ORIGIN %in% origin_remove) 


## sites with 3 species
load(file="output_data/01a_sites_3_species.RData")
head(singleSp3l)


head(sites)
length(unique(sites$SiteID)) # 616

## main df
load(file = "sync/04_temp_pref_env_dist_no_dupl_pairs_WC.RData")
head(allsyncx)

sites2 <- allsyncx %>%
  select(Site_ID1, Site_ID2) %>%
  distinct()
unique(sites2$Site_ID1)


test <- c(unique(sites2$Site_ID1), unique(sites2$Site_ID2))
length(test)

## select columns
sites <- fish_ab %>%
  select(SiteID, Latitude, Longitude, BioRealm,  HydroBasin, Country) %>%
  filter(!SiteID %in% singleSp3l) %>% ## remove 3 species sites
  filter(SiteID %in% test) %>% ## make sure it matches main data
  distinct()

sites

write.csv(sites, "output_data/01a_all_sites_616.csv")
