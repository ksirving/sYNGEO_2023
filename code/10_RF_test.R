## trying a RF on data

library(tidyverse)
library(tidylog)

library(randomForest)

set.seed(234) ## reproducibility

# NUMBER OF BOOTSTRAP REPLICATES
b=10001

sp=0.6 ## for dependence plots

# Upload and Format Data --------------------------------------------------

load(file = "sync/04_temp_pref_env_dist_no_dupl_pairs_WC.RData")
head(allsyncx)
names(allsyncx)

## change names 

allsyncx <- allsyncx %>%
  rename(overlap = diversity2, overlapScaled = diversity3, overlapBayes = diversityBayes) 

## connectiovity as a factor and change name 
allsyncx$Connectivity <- as.factor(allsyncx$Connectivity)
allsyncx$Connectivity <- recode_factor(allsyncx$Connectivity,  "1" = "Within Basin", "0" = "Between Basin") 

# ## within basin - make longer to get site names for membership model
allsyncxWithin <- allsyncx %>%
  filter(Connectivity == "Within Basin") %>%
  pivot_longer(Site_ID1:Site_ID2, names_to = "SiteNumber", values_to = "SiteName") 

## orignial data
originaldata <- read.csv("input_data/Bio/fishdata_selection_basins_same_time_window_10262020.csv") %>%
  select(SiteID, HydroBasin, Waterbody) %>% distinct() %>%
  rename(SiteName = SiteID)
head(originaldata)

## add hydrobasin, waterbody to synchrony

allsyncxWithin <- inner_join(allsyncxWithin, originaldata, by = "SiteName") %>%
  drop_na()

names(allsyncxWithin)

## format for RF
rf.data <- allsyncxWithin %>%
  select(Sync, Region, HydroBasin, WCDistkm, annual_avgC,MeanLat,distance,overlap,  Waterbody) 


## run model
?randomForest
rf.final <- randomForest(Sync ~., data = rf.data, mtry = 2, nodesize=5,
                          importance=TRUE, norm.votes=TRUE, proximity=TRUE)

rf.final$importance


