## format for figures and stats

# packages

library(tidyverse)
library(tidylog)
library("easystats")
library(scales)
getwd()
## directory for figures
out.dir <- "/Users/katieirving/OneDrive - SCCWRP/Documents - Katie’s MacBook Pro/git/sYNGEO_2023/Figures/"


# Temp pref synchrony -----------------------------------------------------


load(file = "sync/03_sync_traits_CWM_CWV_distances.RData")

funcsync <- syncDF %>%
  rename(Sync = synchrony) %>%
  dplyr::select(-X) %>%
  filter( !Site_ID1 == Site_ID2) ## remove pairs comprised of same sites

head(funcsync)

# round(range(allsync$distance),digits = 2)
## checking missing pairs script 5a

# singSpeciesorig <- funcsync %>%
#   filter(Pair %in% otherSites) ## all sync NA values
# dim(singSpeciesorig) ##  1763

# Temp synchrony ----------------------------------------------------------

load(file="sync/03_sync_temp_distances.RData")

tempsyncS <- syncDF %>%
  dplyr::select(-corCoef) %>%
  distinct() %>%
  pivot_wider(names_from = Metric, values_from = synchrony) %>%
  filter( !Site_ID1 == Site_ID2) ## remove pairs comprised of same sites

tempsyncC <- syncDF %>%
  dplyr::select(-synchrony) %>%
  distinct() %>%
  pivot_wider(names_from = Metric, values_from = corCoef) %>%
  filter( !Site_ID1 == Site_ID2) %>% ## remove pairs comprised of same sites
  dplyr::select(Pair, Region, annual_avg) %>%
  rename(annual_avgC = annual_avg)
  
tempsync <- full_join(tempsyncS, tempsyncC, by = c("Pair", "Region"))
## checking missing pairs script 5a

# singSpeciesorig <- tempsync %>%
#   filter(Pair %in% otherSites) 
# dim(singSpeciesorig) ##  1763

# Join --------------------------------------------------------------------

allsync <- full_join(funcsync, tempsync, by = c("Pair","Site_ID1", "Site_ID2", "Region", "Connectivity", "Euclid_Dist_Meters",
                                                "Similarity", "MeanLat", "MeanLon", "MaxDist")) %>%
  mutate(DistKM = Euclid_Dist_Meters/1000)

save(allsync, file = "sync/04_temp_pref_env_dist_ALL.RData")

round(range(na.omit(allsyncx$distance)),digits = 2)

# Format sites ------------------------------------------------------------

## some pairs are the same sites but other way around
head(allsync)
length(unique(allsync$Pair))

sites <- allsync %>%
  ungroup() %>%
  # filter(Region == paste(region[r])) %>%
  # dplyr::select(Site_ID1, Site_ID2,Pair, Region) %>% 
  separate(Pair, into =c("Site_ID2a", "Site_ID1a"), remove = F) %>%
  mutate(rev =  paste0(Site_ID1a, ".", Site_ID2a, sep = "")) # #calculate reverse entries

sites$no <- seq_along(sites$Pair) #number the entries
sites$whererev <- match(sites$rev, sites$Pair) #identify where reversed entries occur
sites$whererev[sites$whererev>sites$no] <- NA #remove the first of each duplicated pair 
sites$Pair[!is.na(sites$whererev)] <- sites$rev[!is.na(sites$whererev)] #replace duplicates

all_sites <- unique(sites$Pair) ## pairs with no reverse duplicate
length(all_sites) ## 207299
dim(sites)

## subset original df 
names(sites)
allsyncx  <- sites %>%
  ungroup() %>%
  dplyr::select(-Site_ID2a, -Site_ID1a, -rev,-no,-whererev) %>%
  distinct(Pair, .keep_all = TRUE) %>%
  drop_na(Sync)

save(allsyncx, file = "sync/04_temp_pref_env_dist_no_dupl_pairs.RData")

# sum(is.na(allsyncx$distance))

# watercourse V euclid distance -------------------------------------------

library(scales)
## euclid distance
load(file = "sync/04_temp_pref_env_dist_no_dupl_pairs.RData") ## allsyncx
head(allsyncx)

## format euclid - remove sites not in same basin
syncsites <- allsyncx %>%
  mutate(DistMetersEuclid = Euclid_Dist_Meters) %>%
  filter(Connectivity == 1) %>%
  filter( !Site_ID1 == Site_ID2) %>% ## remove pairs comrised of same sites
  dplyr::select(Pair, DistMetersEuclid)

## water course distance
watersites <- read.csv2("input_data/sites/Sites_DistancesRiverATLAS_280621.csv")
head(watersites)

## format water course data

watersites <- watersites %>%
  dplyr::mutate(Pair = paste(SiteID_Orig, SiteID_Dest, sep =".")) %>%
  dplyr::mutate(DistMetersWater = as.numeric(TotLong_Meters)) %>%
  dplyr::select(Pair, DistMetersWater, BioRealm, HydroBasin, Country)

head(watersites)
head(syncsites)

length(unique(syncsites$Pair)) # 6076
length(unique(watersites$Pair)) # 17221
sum(watersites$Pair %in% syncsites$Pair) ## ~200 site pairs missing - could be related to filtering of species/traits/sites

all_sites <- left_join(syncsites, watersites, by = "Pair")
head(all_sites)

str(all_sites)

cor(all_sites$DistMetersEuclid, all_sites$DistMetersWater, use = "complete.obs") ## 0.93

## make data long for plot

# all_sites_long <- all_sites %>%
#   pivot_longer(DistMetersEuclid:DistMetersWater, names_to = "Type", values_to = "Meters")
# 
# head(all_sites_long)

t1 <- ggplot(all_sites, aes(x=DistMetersWater/1000, y = DistMetersEuclid/1000)) +
  geom_point() +
  scale_y_log10(name="Log Eucliean Distance (km)", labels = comma) +
  scale_x_log10(name="Log Water Course Distance (km)", labels = comma) 

file.name1 <- paste0(out.dir, "watercourse_v_euclid.jpg")
ggsave(t1, filename=file.name1, dpi=300, height=5, width=6)
t1

ggplot(all_sites, aes(x=DistMetersWater/1000, y = DistMetersEuclid/1000)) +
  geom_point() +
  scale_y_continuous(name="Eucliean Distance (km)", labels = comma, limits = c(0, 1000)) +
  scale_x_continuous(name="Water Course Distance (km)", labels = comma, limits = c(0, 3000))


# Means of each per basin -------------------------------------------------

head(all_sites)

meanDist <- all_sites %>%
  pivot_longer(DistMetersEuclid:DistMetersWater, names_to = "Type", values_to = "distance") %>%
  mutate(BiogeoRegion = case_when(Country %in% c("FIN", "SWE", "GBR", "ESP", "FRA") ~ "Europe",
                                  Country== "AUS" ~ "Oceania", Country == "USA" ~ "USA")) %>%
  group_by(BiogeoRegion, HydroBasin, Type) %>%
  summarise(MeanDist = mean(distance)) %>%
  mutate(MeanDist = MeanDist/1000, HydroBasin = as.character(HydroBasin)) %>%
  pivot_wider(names_from = Type, values_from = MeanDist) %>%
  rename(DistKMEuclid = DistMetersEuclid, DistKMWater = DistMetersWater) %>%
  mutate(DistDifference = DistKMWater-DistKMEuclid) %>%
  drop_na()


## get overall means
rowAdd <- c("All", "Mean", mean(meanDist$DistKMEuclid), 
            mean(meanDist$DistKMWater), mean(meanDist$DistDifference))

rowAdd
write.csv(meanDist, "output_data/03_mean_distances.csv")

# Checking weird sites ----------------------------------------------------

## problem is some euclid distance is longer than water course - water course should always be longer

## get ratio of watercourse/euclid distance

head(all_sites)

ratios <- all_sites %>%
  mutate(Ratio = DistMetersEuclid/DistMetersWater)

head(ratios)

ratios_big <- ratios %>%
  filter(Ratio > 1) %>%
  separate(Pair, into = c("Site_ID1", "Site_ID2"), remove = F)

dim(ratios_big) ## 242

write.csv(ratios_big, "output_data/03_watercourse_smaller_than_euclid_dist.csv")

head(ratios_big)
## add coords

pairs <- ratios_big$Pair
pairs



p=1

## get pair from sync data

pairx <- ratios_big %>%
  filter(Pair == pairs[p])
pairx
## define sites
S1 <- pairx$Site_ID1
S2 <- pairx$Site_ID2

S1

## get coords for each site
CoordsS1 <- SiteCoords %>%
  filter(SiteID == S1) %>%
  dplyr::select(Longitude, Latitude, SiteID)

CoordsS1

CoordsS2 <- SiteCoords %>%
  filter(SiteID == S2) %>%
  dplyr::select(Longitude, Latitude, SiteID)

CoordsS2

sp::coordinates(CoordsS1) <- c("Longitude", "Latitude")
sp::coordinates(CoordsS2) <- c("Longitude", "Latitude")

#Make a distance matrix
dst <- pointDistance(CoordsS1,CoordsS2, lonlat=TRUE)
dst
# str(dst)
# get mean latitude/longitude
MeanLat <- (CoordsS1$Latitude+CoordsS2$Latitude)/2
MeanLon <- (CoordsS1$Longitude+CoordsS2$Longitude)/2

## add to dataframe
syncsites[p,8] <- dst
syncsites[p,10] <- MeanLat
syncsites[p,11] <- MeanLon
head(syncsites)


# Add water course distance to synchrony ----------------------------------

head(watersites)

## convert to km

watersites <- watersites %>%
  mutate(WCDistkm = DistMetersWater/1000) %>%
  dplyr::select(-BioRealm)


head(allsyncx)


syncDF <- left_join(allsyncx, watersites, by = c("Pair") ) 
syncDF
syncDF <- syncDF %>%
  filter(Connectivity == 1) %>%
  mutate(DistKM = Euclid_Dist_Meters/1000 ) %>%
  dplyr::select(Pair, DistKM, WCDistkm)
  # pivot_longer(c(Euclid_Dist_KM, WCDistkm), names_to = "Dist_Type", values_to = "KM")


# 
# sm1a <- ggplot(syncDF, aes(x=KM, y=synchrony, color = Dist_Type)) +
#   geom_smooth(method = "lm") +
#   # facet_wrap(~Trait) +
#   scale_color_discrete(name = "Distance Type", labels = c("Euclidean", "Water Course")) +
#   scale_y_continuous(name="Synchrony") +
#   scale_x_log10(name="Log Distance (km)", labels = comma) 
# sm1a
# 
# 
# file.name1 <- paste0(out.dir, "Euclid_v_waterCourse_distance.jpg")
# ggsave(sm1a, filename=file.name1, dpi=300, height=5, width=6)

head(syncDF)


## upload 
load(file = "sync/04_temp_pref_env_dist_no_dupl_pairs.RData")
head(allsyncx)


## join water course distance

allsyncx <- full_join(allsyncx, syncDF, by = c("Pair", "DistKM"))


save(allsyncx, file="sync/04_temp_pref_env_dist_no_dupl_pairs_WC.RData")

round(range(allsyncx$distance),digits = 2)

