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

## save out

save(singleSpl, file="output_data/01a_sites_1_species.RData")
save(singleSp2l, file="output_data/01a_sites_2_species.RData")
save(singleSp3l, file="output_data/01a_sites_3_species.RData")


