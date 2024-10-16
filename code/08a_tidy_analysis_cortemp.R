library(ggeffects)
library(tidyverse)
library(tidylog)
library(cowplot) #for manuscript ready figures
library(sjPlot) #for plotting lmer and glmer mods
library(sjmisc)
library(sjlabelled)
library(lmerTest) ## from Luke et al 2016 - evaluating significance in linear mixed-effects models in R
library("easystats") ## multicolinearality
library(lmerMultiMember)
library(effects)
library(sjstats) #use for r2 functions
library("scales")
library(performance)
library(lme4)

# library("devtools")
# install_github("jvparidon/lmerMultiMember")

load(file = "sync/04_temp_pref_env_dist_no_dupl_pairs_WC.RData")
head(allsyncx)
names(allsyncx)

## sites with 3 species
load(file="output_data/01a_sites_3_species.RData")
head(singleSp3l)
# length(singleSp3l) ## 151 - removed earlier in process

## change names and filter single species sites

allsyncx <- allsyncx %>%
  rename(overlap = diversity2, overlapScaled = diversity3, overlapBayes = diversityBayes) %>%
  filter(!Site_ID1 %in% singleSp3l, !Site_ID2 %in% singleSp3l) #%>%

## remove NAs from temp column
allsyncx <- allsyncx %>%
  drop_na(annual_avg)

## directory for figures
out.dir <- "/Users/katieirving/OneDrive - SCCWRP/Documents - Katie’s MacBook Pro/git/sYNGEO_2023/Figures/Tidy/"

cor(allsyncx$annual_avg, allsyncx$annual_avgC) ## 0.998


#  Make z scores and histograms ----------------------------------------------------------

allsyncx <- allsyncx %>%
  group_by(Connectivity) %>% ## not grouping by z scores to enable back transformation
  mutate(ZTemp = (annual_avg - mean(annual_avg))/sd(annual_avg), ## z scores
         ZTempCor = (annual_avgC - mean(annual_avgC))/sd(annual_avgC),
         ZDistance = (distance - mean(distance))/sd(distance),
         ZDistanceAbs = (distanceAbs - mean(distanceAbs))/sd(distanceAbs),
         ZOverlap = (overlap - mean(overlap))/sd(overlap),
         ZOverlapScaled = (overlapScaled - mean(overlapScaled))/sd(overlapScaled),
         ZOverlapBayes = (overlapBayes - mean(overlapBayes))/sd(overlapBayes),
         WCDistkmsqrt = sqrt(WCDistkm),
         DistKMLog = log(DistKM+1),
         DistKMsqrt = sqrt(DistKM)) ## square root distances

## number of sites
length(unique(allsyncx$Site_ID2)) ## 613
length(unique(allsyncx$Site_ID1)) ## 613

hist(allsyncx$Sync)

hist(log(allsyncx$distance))
hist(allsyncx$distance)

hist(allsyncx$overlap)
hist(log(allsyncx$overlap))

hist(allsyncx$overlapScaled)
hist(allsyncx$overlapBayes)
# hist(log(allsyncx$overlapScaled))

hist(allsyncx$overlap)
hist(log(allsyncx$overlap))

hist(log(allsyncx$annual_avg))
hist(allsyncx$annual_avg)

hist(allsyncx$DistKM)
hist(allsyncx$WCDistkm)
hist(allsyncx$WCDistkmsqrt)

# CORRELATIONS ------------------------------------------------------------

## all
allsyncxCor <- allsyncx %>%
  ungroup() %>%
  select(distance, overlap, overlapScaled, overlapBayes, annual_avg, DistKM)

syncCor <- cor(allsyncxCor)
syncCor
write.csv(syncCor, "output_data/Tidy/08_cor_dist_overlap_temp_n4_sites.csv")

# distance    overlap annual_avg     DistKM
# distance    1.0000000 -0.6033331 -0.6745187  0.6549164
# overlap    -0.6033331  1.0000000  0.3591466 -0.3927173
# annual_avg -0.6745187  0.3591466  1.0000000 -0.8928775
# DistKM      0.6549164 -0.3927173 -0.8928775  1.0000000

## within basin
allsyncxCorWithin <- allsyncx %>%
  ungroup() %>%
  filter(Connectivity == 1) %>%
  select(distance, overlap, annual_avg, DistKM, WCDistkm)

syncCorwithin <- cor(allsyncxCorWithin, use = "complete.obs")
syncCorwithin
#             distance     overlap  annual_avg      DistKM    WCDistkm
# distance    1.0000000 -0.74592662 -0.12016829  0.16799780  0.11205326
# overlap    -0.7459266  1.00000000  0.04334475 -0.09658227 -0.02503594
# annual_avg -0.1201683  0.04334475  1.00000000 -0.90538556 -0.90044643
# DistKM      0.1679978 -0.09658227 -0.90538556  1.00000000  0.93188165
# WCDistkm    0.1120533 -0.02503594 -0.90044643  0.93188165  1.00000000

write.csv(syncCorwithin, "output_data/Tidy/08_cor_dist_overlap_temp_within_Basin_n4_sites.csv")

## between basin
allsyncxCorBetween <- allsyncx %>%
  ungroup() %>%
  filter(Connectivity == 0) %>%
  select(distance, overlap, annual_avg, DistKM)

syncCorBetween <- cor(allsyncxCorBetween)
syncCorBetween

# distance    overlap annual_avg     DistKM
# distance    1.0000000 -0.5905045 -0.6631846  0.6417011
# overlap    -0.5905045  1.0000000  0.3285168 -0.3592465
# annual_avg -0.6631846  0.3285168  1.0000000 -0.8866043
# DistKM      0.6417011 -0.3592465 -0.8866043  1.0000000

write.csv(syncCorBetween, "output_data/Tidy/08_cor_dist_overlap_temp_between_Basin_n4_sites.csv")

syncCorwithin
syncCorBetween

### z scores
allsyncxCor <- allsyncx %>%
  ungroup() %>%
  select(ZDistance, ZOverlap, ZTemp, DistKM)

syncCor <- cor(na.omit(allsyncxCor))
syncCor
write.csv(syncCor, "output_data/Tidy/08_cor_dist_overlap_temp_n4_sites_zscores.csv")

# ZDistance   ZOverlap      ZTemp     DistKM
# ZDistance  1.0000000 -0.5999192 -0.6262589  0.5949952
# ZOverlap  -0.5999192  1.0000000  0.3093605 -0.3314823
# ZTemp     -0.6262589  0.3093605  1.0000000 -0.8244516
# DistKM     0.5949952 -0.3314823 -0.8244516  1.00000000

## within basin
allsyncxCorWithin <- allsyncx %>%
  ungroup() %>%
  filter(Connectivity == 1) %>%
  select(ZDistance, ZOverlap, ZTemp, DistKM)
allsyncxCorWithin
syncCorwithin <- cor(allsyncxCorWithin, use = "complete.obs")
syncCorwithin
# ?cor

#           ZDistance    ZOverlap       ZTemp     DistKM
# ZDistance  1.0000000 -0.75490908 -0.13668478  0.1418257
# ZOverlap  -0.7549091  1.00000000  0.08565081 -0.1275028
# ZTemp     -0.1366848  0.08565081  1.00000000 -0.4478271
# DistKM     0.1418257 -0.12750279 -0.44782708  1.0000000

write.csv(syncCorwithin, "output_data/Tidy/08_cor_dist_overlap_temp_within_Basin_n4_sites_zscores.csv")

## between basin
allsyncxCorBetween <- allsyncx %>%
  ungroup() %>%
  filter(Connectivity == "Between Basin") %>%
  select(ZDistance, ZOverlap, ZTemp, DistKM)

syncCorBetween <- cor(allsyncxCorBetween)
syncCorBetween

#             ZDistance   ZOverlap      ZTemp     DistKM
# ZDistance  1.000000 -0.5918300 -0.6518110  0.6306360
# ZOverlap  -0.591830  1.0000000  0.3210365 -0.3506969
# ZTemp     -0.651811  0.3210365  1.0000000 -0.8705121
# DistKM     0.630636 -0.3506969 -0.8705121  1.0000000

write.csv(syncCorBetween, "output_data/Tidy/08_cor_dist_overlap_temp_between_Basin_n4_sites_zscores.csv")

## correaltions per region

unique(allsyncxWithin$Region)

corsEU <- allsyncxWithin %>%
  ungroup() %>%
  filter(Region == "Europe") %>%
  select(distance, overlap, annual_avgC, Sync, WCDistkm)

cEU <- cor(corsEU, use = "complete.obs")

write.csv(cEU, "output_data/Tidy/08_cor_dist_abs_vals_EU.csv")

corsUS <- allsyncxWithin %>%
  ungroup() %>%
  filter(Region == "USA") %>%
  select(distance, overlap, annual_avgC, Sync, WCDistkm)

cUS <- cor(corsUS, use = "complete.obs")
cUS

write.csv(cUS, "output_data/Tidy/08_cor_dist_abs_vals_USA.csv")

corsAU <- allsyncxWithin %>%
  ungroup() %>%
  filter(Region == "Oceania") %>%
  select(distance, overlap, annual_avgC, Sync, WCDistkm)
corsAU
cAU <- cor(corsAU, use = "complete.obs")
cAU

write.csv(cUS, "output_data/Tidy/08_cor_dist_abs_vals_Oceania.csv")




# Prep Data for models ----------------------------------------------------

## connectiovity as a factor and change name 
allsyncx$Connectivity <- as.factor(allsyncx$Connectivity)
allsyncx$Connectivity <- recode_factor(allsyncx$Connectivity,  "1" = "Within Basin", "0" = "Between Basin") 

# ## within basin - make longer to get site names for membership model
allsyncxWithin <- allsyncx %>%
  filter(Connectivity == "Within Basin") %>%
  pivot_longer(Site_ID1:Site_ID2, names_to = "SiteNumber", values_to = "SiteName") 

allsyncxWithin
# ## between basin -  make longer to get site names for membership model
allsyncxBetween <- allsyncx %>%
  filter(Connectivity == "Between Basin") %>%
  pivot_longer(Site_ID1:Site_ID2, names_to = "SiteNumber", values_to = "SiteName") 

round(range(allsyncxBetween$distance),digits = 2) ## 0.00 0.31

# N basins etc ------------------------------------------------------------

## upload fish abundance and site data
originaldata <- read.csv("input_data/Bio/fishdata_selection_basins_same_time_window_10262020.csv") %>%
  dplyr::select(SiteID, BioRealm, HydroBasin, Species) 
head(originaldata)

sBa <- allsyncxWithin %>% 
  dplyr::select(SiteName, Pair, Region, WCDistkm) %>%
  rename(SiteID = SiteName)

length(unique(sBa$SiteID)) ## 616

## join
test <- right_join(originaldata, sBa, by = "SiteID", relationship = "many-to-many")

head(test)
length(unique(test$SiteID)) ## 616



test1 <- test %>% group_by(Region) %>% summarise(nSites = length(unique(SiteID)),
                                                nSpecies = length(unique(Species)),
                                                nPair = length(unique(Pair)),
                                                nBasins = length(unique(HydroBasin)),
                                                mWCDist = mean(na.omit(WCDistkm)))
test1

write.csv(test1, "output_data/Tidy/08_n_sites_species_per_region.csv")
unique(test$BioRealm)


## look at overlap percentiles, check wc/eu distances

round(quantile(allsyncxWithin$ZOverlap, probs = c(0.04, 0.25,0.5,0.75, 0.9)), digits = 3)
round(quantile(allsyncxWithin$Sync, probs = c(0.25,0.5,0.75, 0.9)), digits = 3)

OverT <- allsyncxWithin %>%
  filter(ZOverlap <= -1.450) %>%
  select(Region, Sync:overlap, annual_avgC:WCDistkm)

min(allsyncx$ZOverlap)
names(allsyncxWithin)

## sumarry stats per region

SumT <- OverT %>%
  pivot_longer(c(distance, overlap, annual_avgC, Sync, DistKM, WCDistkm), names_to = "Variable", values_to = "Values") %>%
  group_by(Region, Variable) %>%
  summarise(MinVals = min(na.omit(Values)), MaxVals = max(na.omit(Values)),
            MeanVals = mean(na.omit(Values)), SdtError = sder(na.omit(Values))) %>%
  pivot_longer(MinVals:SdtError, names_to = "stat", values_to = "statVal") %>%
  mutate(statVal = round(statVal, digits = 2)) %>%
  pivot_wider(names_from = Variable, values_from = statVal)

# Histograms of all vars

OverT_wide <- OverT %>% 
  pivot_longer(c(Sync:WCDistkm), names_to = "Var", values_to = "Vals")
  
OverT_wide

## plot all hists per region
ggplot(OverT_wide, aes(Vals, colour = Var), fill = NA) +
  geom_histogram() +
  facet_grid(rows = vars(Region), cols = vars(Var), scales = "free")

## plot relationships with sync
unique(relsT$name)

names(allsyncxWithin)
relsT <- allsyncxWithin %>%
  select(Region, Sync, distance, overlapScaled, annual_avgC, WCDistkm) %>%
  pivot_longer(distance:WCDistkm) %>%
  mutate(name = factor(name, levels = c("annual_avgC", "distance", "overlapScaled","WCDistkm"),
                       labels = c("Temp Synchrony", "Trait Distance", "Trait Overlap", "Water Course (km)" )))

## plot all relationships per region
p1 <- ggplot(relsT, aes(y=Sync, x=value, colour = name)) +
  geom_smooth(method = "lm") +
  facet_grid(rows = vars(Region), cols = vars(name), scales = "free") +
  labs(y = "Trait Synchrony") +
  theme(legend.position = "none",
        legend.title = element_blank(),
        legend.text=element_text(size=15),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        # axis.title.y = element_text("Trait Synchrony"),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15)) 

p1

file.name1 <- paste0(out.dir, "08_all_relationships_raw_data.jpg")
ggsave(p1, filename=file.name1, dpi=600, height=6, width=10)

## take out region
p2 <- ggplot(relsT, aes(y=Sync, x=value, colour = name)) +
  geom_smooth(method = "lm") +
  facet_wrap(~name, scales = "free") 

p2

file.name1 <- paste0(out.dir, "08_all_relationships_raw_data_no region.jpg")
ggsave(p2, filename=file.name1, dpi=300, height=5, width=6)
## overlap in eu is negative 

## plot relationships with overlap

relsT2 <- allsyncxWithin %>%
  select(Sync:overlapScaled, Region,annual_avgC:WCDistkm) %>%
  pivot_longer(c(Sync, distance, annual_avgC:WCDistkm))

## plot all relationships per region
ggplot(relsT2, aes(y=overlapScaled, x=value, colour = name)) +
  geom_smooth(method = "lm") +
  facet_grid(rows = vars(Region), cols = vars(name), scales = "free") +
  scale_y_continuous(limits=c(0,1))


## plot all relationships per region - points
ggplot(relsT2, aes(y=overlap, x=value, colour = name)) +
  geom_point() +
  facet_grid(rows = vars(Region), cols = vars(name), scales = "free") +
  scale_y_continuous(limits=c(0,1))


names(relsT2)

## look at values again
OverT2 <- allsyncxWithin %>%
  filter(ZOverlap >= 1.322, Sync >=0.833) 

## join with original data to get site info
head(originaldata)

test <- inner_join(OverT2, originaldata, by = c("SiteName" = "SiteID"))

length(unique(test$HydroBasin))

# Extract Values ----------------------------------------------------------
sder <- function(x) sd(x)/sqrt(length(x))

## values of drivers
range(allsyncxWithin$Sync) ## 0.06639871 0.99291205
round(range(allsyncxWithin$distance),digits = 2) ## 0.00 0.16
round(range(allsyncxWithin$overlap),digits = 2) ## 0.00 0.97
round(range(allsyncxWithin$annual_avgC),digits = 2) ## 0.82 1.00
round(range(allsyncxWithin$overlapScaled),digits = 2) ## 0.00 0.99
names(allsyncxWithin)

ranges <- allsyncxWithin %>%
  pivot_longer(c(distance, overlap,overlapScaled, annual_avgC, Sync, WCDistkm), names_to = "Variable", values_to = "Values") %>%
  group_by(Variable) %>%
  summarise(MinVals = min(na.omit(Values)), MaxVals = max(na.omit(Values)), 
            MeanVals = mean(na.omit(Values)), SdtError = sder(na.omit(Values)))
ranges
## save 
write.csv(ranges, "output_data/Tidy/08_values_overall.csv")

## per region
rangesRegion <- allsyncxWithin %>%
  pivot_longer(c(distance, overlap, overlapScaled, annual_avgC, Sync, WCDistkm), names_to = "Variable", values_to = "Values") %>%
  group_by(Region, Variable) %>%
  summarise(MinVals = min(Values), MaxVals = max(Values), MeanVals = mean(Values), SdtError = sder(Values))
rangesRegion

## save
write.csv(rangesRegion, "output_data/Tidy/08_values_per_region.csv")

## wc distances per region

rangesRegionDist <- allsyncxWithin %>%
  pivot_longer(c(DistKM, WCDistkm), names_to = "Variable", values_to = "Values") %>%
  group_by(Region, Variable) %>%
  summarise(MinVals = min(na.omit(Values)), MaxVals = max(na.omit(Values)),
                          MeanVals = mean(na.omit(Values)), SdtError = sder(na.omit(Values))) %>%
  pivot_longer(MinVals:SdtError, names_to = "stat", values_to = "statVal") %>%
  mutate(statVal = round(statVal, digits = 2)) %>%
  pivot_wider(names_from = stat, values_from = statVal)

rangesRegionDist

## save
write.csv(rangesRegionDist, "output_data/Tidy/08_geo_distances_per_region.csv")


# Models: Within Basin spatial decay ----------------------------------------------------

## within Basin - water course distance
## create memberships
Wa <- lmerMultiMember::weights_from_vector(allsyncxWithin$Region)
Wj <- Matrix::fac2sparse(allsyncxWithin$SiteName)  # convert single membership vars to an indicator matrix with fac2sparse()
Waj <- interaction_weights(Wa, Wj)

mem_mixedwc <- lmerMultiMember::lmer(Sync ~ WCDistkmsqrt
                                     + (1 | Region ) + ## add predictors here to get random effect per region
                                       + (1 | RegionXSiteName), 
                                     memberships = list(Region = Wa, RegionXSiteName = Waj), 
                                     REML = T,
                                     data = allsyncxWithin)



## coefs
summary(mem_mixedwc, ddf = "Satterthwaite")
anova(mem_mixedwc, ddf = "Satterthwaite")
r2_nakagawa(mem_mixedwc) ## 0.15
check_singularity(mem_mixedwc) ## False
performance::icc(mem_mixedwc, by_group = T) ## 0.12

## get and save coefs
modwc <- data.frame(coef(summary(mem_mixedwc, ddf = "Satterthwaite")))

df0 <- NULL
df0$Effects <- rownames(modwc)
df0$Estimates <- modwc$Estimate ## estimates
df0$DegreesOfFreedom <- modwc$df ## df
df0$TStat <- modwc$t.value ## t statistics
df0$Pvalues <- modwc$Pr...t.. ## p values
df0 <- as.data.frame(df0)
df0

## fitted values add to df
allsyncxWithinx <- allsyncxWithin %>%
  drop_na(WCDistkmsqrt) ## remove WC NAs

class(mem_mixedwc) <- "lmerModLmerTest"

allsyncxWithinx$mem_mixedwc <- fitted(mem_mixedwc)

set_theme(base = theme_classic(), #To remove the background color and the grids
          theme.font = 'serif',   #To change the font type
          axis.title.size = 1.3,  #To change axis title size
          axis.textsize.x = 1.2,    #To change x axis text size
          # axis.angle.x = 60,      #To change x axis text angle
          # axis.hjust.x = 1,
          # axis.ticksmar = 50,
          axis.textsize.y = 1)  #To change y axis text size


S2 <- ggplot(allsyncxWithinx, aes(x = WCDistkmsqrt, y = mem_mixedwc)) +
  geom_point(aes(y=Sync, col = Region), size = 0.01) +
  geom_smooth(method = "lm") +
  scale_y_continuous(name="Spatial Synchrony in Thermal Preference (CWM)") +
  scale_x_continuous(name="Water Course Distance (sqrt KM)")
S2

file.name1 <- paste0(out.dir, "08_Fitted_Sync_Over_WCDistance_Within_sqrt.jpg")
ggsave(S2, filename=file.name1, dpi=300, height=5, width=6)



# Biotic Variables: mixed membership model --------------------------------

## distance: Within Basin

Wa <- lmerMultiMember::weights_from_vector(allsyncxWithin$Region)
Wj <- Matrix::fac2sparse(allsyncxWithin$SiteName)  # convert single membership vars to an indicator matrix with fac2sparse()
Waj <- interaction_weights(Wa, Wj)

mem_mixed0 <- lmerMultiMember::lmer(Sync ~  ZDistance*ZTempCor
                                    + (1 | Region ) + ## add predictors here to get random effect per region
                                      (1 | RegionXSiteName), 
                                    memberships = list(Region = Wa, RegionXSiteName = Waj), 
                                    REML = T,
                                    data = allsyncxWithin)
# 
# mem_mixed0a <- lmerMultiMember::lmer(Sync ~  ZDistanceAbs*ZTempCor
#                                     + (1 | Region ) + ## add predictors here to get random effect per region
#                                       (1 | RegionXSiteName), 
#                                     memberships = list(Region = Wa, RegionXSiteName = Waj), 
#                                     REML = T,
#                                     data = allsyncxWithin)

summary(mem_mixed0, ddf = "Satterthwaite")
anova(mem_mixed0, ddf = "Satterthwaite")
r2_nakagawa(mem_mixed0) ## 0.151
r2 <- r2_nakagawa(mem_mixed0)[1]
check_singularity(mem_mixed0) ## False
performance::icc(mem_mixed0, by_group = T) ## 0.122

## get and save coefs
mod0 <- data.frame(coef(summary(mem_mixed0, ddf = "Satterthwaite")))

df0 <- NULL
df0$Effects <- rownames(mod0)
df0$Estimates <- mod0$Estimate ## estimates
df0$DegreesOfFreedom <- mod0$df ## df
df0$TStat <- mod0$t.value ## t statistics
df0$Pvalues <- mod0$Pr...t.. ## p values
df0 <- as.data.frame(df0)

df0$R2 <- r2$R2_conditional ## 0.153
df0$ICCRxSN <- performance::icc(mem_mixed0, by_group = T)[1,2] ## 0.124
df0$ICCR <- performance::icc(mem_mixed0, by_group = T)[2,2]

df0

write.csv(df0, "Tables/08_distance_coefs.csv")

### plots 
class(mem_mixed0) <- "lmerModLmerTest"
devtools::install_github("strengejacke/strengejacke")

data(iris)
m2 <- lm(Sepal.Length ~ Sepal.Width + Petal.Length + Species, data = iris)
iris
# variable names as labels, but made "human readable"
# separating dots are removed
plot_model(m2)


estsdw <- sjPlot::plot_model(mem_mixed0, 
                             type = "est",
                             axis.lim = c(-0.05, 0.05),
                             vline.color = "grey30",
                             show.values=TRUE, show.p=TRUE,
                             axis.labels = c("Interaction", "Temperature Synchrony", "Thermal Distance"),
                             title = "")

estsdw <- estsdw + scale_y_continuous(limits = c(-0.05, 0.05))

file.name1 <- paste0(out.dir, "08_withinBasin_effect_sizes_all_expsync_n4_sites_distance_zscores_cortemp.jpg")
ggsave(estsdw, filename=file.name1, dpi=300, height=8, width=10)

# estsTab <- sjPlot::tab_model(mem_mixed0,
#                              show.re.var= TRUE,
#                              pred.labels =c("(Intercept)", "ZDistance", "ZTemp"),
#                              dv.labels= "Drivers of Thermal Synchrony")
# 
# estsTab


theme_set(theme_sjplot())

### make individual plots of interactions
mean(allsyncxWithin$overlap)
quantile(allsyncxWithin$ZTempCor, probs = c(0.05,0.5,0.95))
round(quantile(allsyncxWithin$ZOverlap, probs = seq(0.05,1,0.05)), digits = 3)
round(quantile(allsyncxWithin$ZDistance, probs = seq(0.05,1,0.05)), digits = 3)

## overlap v temp
tempDis <- sjPlot::plot_model(mem_mixed0, type="pred", terms= c("ZDistance", "ZTempCor [-2.5088766, 0.4314952, 0.4841359]"),
                              axis.title = c( "Trait Distance (Z Score)","Thermal Trait Synchrony"), 
                              legend.title = "Environmental synchrony")
tempDis


DisTemp <- sjPlot::plot_model(mem_mixed0, type="pred", terms= c("ZTempCor", "ZDistance [-1.013,  -0.892,  -0.759,  -0.608,  -0.422,  0.489,   1.036,    2.034,  4.744 ]"),
                              axis.title = c("Environmental synchrony (Z Score)" ,"Thermal Trait Synchrony"), 
                              legend.title = "Trait Distance")
DisTemp

# allsyncxBetweenx$mem_mixed0 <- fitted(mem_mixed0)

## Overlap: Within Basin
# mem_mixed1 <- lmerMultiMember::lmer(Sync ~  ZOverlap*ZTempCor
#                                     + (1 | Region ) + ## add predictors here to get random effect per region
#                                       (1 | RegionXSiteName), 
#                                     memberships = list(Region = Wa, RegionXSiteName = Waj), 
#                                     REML = T,
#                                     data = allsyncxWithin)

mem_mixed1 <- lmerMultiMember::lmer(Sync ~  ZOverlapScaled*ZTempCor
                                    + (1 | Region ) + ## add predictors here to get random effect per region
                                      (1 | RegionXSiteName),
                                    memberships = list(Region = Wa, RegionXSiteName = Waj),
                                    REML = T,
                                    data = allsyncxWithin)

# mem_mixed1b <- lmerMultiMember::lmer(Sync ~  ZOverlapBayes*ZTempCor
#                                      + (1 | Region ) + ## add predictors here to get random effect per region
#                                        (1 | RegionXSiteName), 
#                                      memberships = list(Region = Wa, RegionXSiteName = Waj), 
#                                      REML = T,
#                                      data = allsyncxWithin)


summary(mem_mixed1, ddf = "Satterthwaite")

anova(mem_mixed1, ddf = "Satterthwaite")
r2_nakagawa(mem_mixed1)[1] ## 0.153
r2 <- r2_nakagawa(mem_mixed1)[1]
r2$R2_conditional
check_singularity(mem_mixed1) ## False
performance::icc(mem_mixed1, by_group = T)[1,2] ## 0.124

mem_mixed1

## get and save coefs
mod1 <- data.frame(coef(summary(mem_mixed1, ddf = "Satterthwaite")))

df1 <- NULL
df1$Effects <- rownames(mod1)
df1$Estimates <- mod1$Estimate ## estimates
df1$DegreesOfFreedom <- mod1$df ## df
df1$TStat <- mod1$t.value ## t statistics
df1$Pvalues <- mod1$Pr...t.. ## p values
df1 <- as.data.frame(df1)

df1$R2 <- r2$R2_conditional ## 0.153
df1$ICCRxSN <- performance::icc(mem_mixed1, by_group = T)[1,2] ## 0.124
df1$ICCR <- performance::icc(mem_mixed1, by_group = T)[2,2]

df1

write.csv(df1, "Tables/08_overlap_coefs.csv")


### plots 
class(mem_mixed1) <- "lmerModLmerTest"
# class(mem_mixed1)

estsow <- sjPlot::plot_model(mem_mixed1, 
                             type = "est",
                             axis.lim = c(-0.05, 0.05),
                             vline.color = "grey",
                             show.values=TRUE, show.p=TRUE,
                             axis.labels = c("Interaction", "Temperature Synchrony", "Thermal Overlap"),
                             title = "")

estsow  <- estsow + scale_y_continuous(limits = c(-0.05, 0.05))

file.name1 <- paste0(out.dir, "08_withinBasin_effect_sizes_all_expsync_n4_sites_overlap_zscores_corTemp.jpg")
ggsave(estsow, filename=file.name1, dpi=300, height=8, width=10)

# estsTab <- sjPlot::tab_model(mem_mixed1,
#                              show.re.var= TRUE,
#                              # pred.labels =c("(Intercept)", "ZOverlap", "ZTempCor"),
#                              dv.labels= "Drivers of Thermal Synchrony")

set_theme(base = theme_classic(), #To remove the background color and the grids
          theme.font = 'serif',   #To change the font type
          axis.title.size = 1.3,  #To change axis title size
          axis.textsize.x = 1.2,    #To change x axis text size
          # axis.angle.x = 60,      #To change x axis text angle
          # axis.hjust.x = 1,
          # axis.ticksmar = 50,
          axis.textsize.y = 1)  #To change y axis text size

theme_set(theme_sjplot())

### make individual plots of interactions
quantile(allsyncxWithin$ZTemp, probs = c(0.05,0.5,0.95))
round(quantile(allsyncxWithin$ZOverlapScaled, probs = c(0.05,0.5,0.95)), digits = 3)
round(quantile(allsyncxWithin$ZOverlap, probs = c(0.05,0.5,0.95)), digits = 3)
round(quantile(allsyncxWithin$ZOverlapBayes, probs = c(0.05,0.5,0.95)), digits = 3)
# round(quantile(allsyncxWithin$ZDistance, probs = c(0.05,0.5,0.95)), digits = 3)

## overlap v temp
tempDiv <- sjPlot::plot_model(mem_mixed1b, type="pred", terms= c("ZOverlapBayes", "ZTempCor [-2.513, 0.436, 0.49]"),
                              axis.title = c( "Trait Overlap (Z Score)","Thermal Trait Synchrony"), 
                              legend.title = "Environmental synchrony")
tempDiv


DivTemp <- sjPlot::plot_model(mem_mixed1b, type="pred", terms= c("ZTempCor", " ZOverlapBayes [-1.794, 0.092, 1.457 ]"),
                              axis.title = c("Environmental synchrony (Z Score)" ,"Thermal Trait Synchrony"), 
                              legend.title = "Trait Overlap")
DivTemp
## overlap v temp
tempDiv <- sjPlot::plot_model(mem_mixed1a, type="pred", terms= c("ZOverlapScaled", "ZTempCor [-1.780, 0.048, 1.560]"),
                              axis.title = c( "Trait Overlap (Z Score)","Thermal Trait Synchrony"), 
                              legend.title = "Environmental synchrony")
tempDiv


DivTemp <- sjPlot::plot_model(mem_mixed1a, type="pred", terms= c("ZTempCor", " ZOverlapScaled [-2.089, 0.261, 1.101]"),
                              axis.title = c("Environmental synchrony (Z Score)" ,"Thermal Trait Synchrony"), 
                              legend.title = "Trait Overlap")
DivTemp


# Figures -----------------------------------------------------------------

## estimates
eststplot <- plot_grid(list(estsdw, estsow)) 

file.name1 <- paste0(out.dir, "08_all_estimates.jpg")
ggsave(eststplot, filename=file.name1, dpi=600, height=5, width=7)

## Distance
distplot <- plot_grid(list(estsdw, DisTemp)) ## grid figures

file.name1 <- paste0(out.dir, "08_distance_zscores_between_within_ests_rels_corTemp.jpg")
ggsave(distplot, filename=file.name1, dpi=300, height=12, width=18)

## Overlap

divplot <- plot_grid(list(estsow, estsob, DivTemp,DivTempB)) ## grid figures

file.name1 <- paste0(out.dir, "08_overlap_zscores_between_within_ests_rels_corTemp.jpg")
ggsave(divplot, filename=file.name1, dpi=300, height=12, width=18)


# Within Basin figures: incomplete ----------------------------------------------------

## data 
allsyncxWithin 
## distance 
allsyncxWithin$DistanceWithinFit <- fitted(mem_mixed0)
range(allsyncxWithin$ZDistance) ## -1.348625  4.768400
round(quantile(allsyncxWithin$ZDistance, probs = c(0.05,0.5,0.95)), digits = 3)
## ZDistance [-1.094, -0.318, 2.004]

### add categories of zscore - quantiles
## 5% = -1.348625 to -1.094
## 50% = -0.8 to 0
## 95% = 2.004 to 4.768400
?plot_model

## overlap
allsyncxWithin$OverlapWithinFit <- fitted(mem_mixed1)
##  ZOverlap [-1.436, 0.053, 1.499]

range(allsyncxWithin$ZOverlap)
## -2.100778  1.928271

round(quantile(allsyncxWithin$ZOverlap, probs = c(0.05,0.5,0.95)), digits = 3)

## 5% = < -1.436
## 50th = >-0.8 < 0.2
## 95th = > 1.499

## back transform the values
mean.vari <- mean(allsyncxWithin$annual_avg)
sd.vari <- sd(allsyncxWithin$annual_avg)
var.back <- (allsyncxWithin$ZTemp * sd.vari) + mean.vari


names(allsyncxWithin)
### add categories of 
FigureData <- allsyncxWithin %>%
  select(Sync, Pair, Region, distance, overlap, annual_avg, ZDistance, ZOverlap, ZTemp, DistanceWithinFit, OverlapWithinFit) %>%
  mutate(DistanceThresh = ifelse(ZDistance < -1.094, "5thPercentile", NA),
         DistanceThresh = ifelse(ZDistance > -0.4 & ZDistance < -0.320, "50thPercentile", DistanceThresh),
         DistanceThresh = ifelse(ZDistance > 2.004 , "95thPercentile", DistanceThresh)) %>%
  mutate(OverlapThresh = ifelse(ZOverlap < -1.436, "5thPercentile", NA),
         OverlapThresh = ifelse(ZOverlap > -0.8 & ZOverlap < 0.2, "50thPercentile", OverlapThresh),
         OverlapThresh = ifelse(ZOverlap > 1.499 , "95thPercentile", OverlapThresh))

FigureData <- allsyncxWithin %>%
  select(Sync, Pair, Region, distance, overlap, annual_avg, ZDistance, ZOverlap, ZTemp, DistanceWithinFit, OverlapWithinFit) %>%
  mutate(DistancePercentiles = percent_rank(ZDistance)*100)

head(FigureData)  

FigureDataDist <- FigureData %>%
  drop_na(DistanceThresh)

min(allsyncxWithin$ZTemp)
max(allsyncxWithin$ZTemp)

DisTempx <- sjPlot::plot_model(mem_mixed0, type="pred", terms= c("ZTempCor", "ZDistance [-1.094, -0.318, 2.004]"),
                               axis.title = c("Environmental synchrony (Z Score)" ,"Thermal Trait Synchrony"), 
                               legend.title = "Trait Distance Percentile") + 
  scale_color_discrete(labels = c("5th", "50th", "95th")) +
  scale_x_continuous( sec.axis = sec_axis(~.*sd.vari+mean.vari, name="Environmental synchrony (Original)" ))

DisTempx

unique(FigureData$Region)
Waj
new_df       <- expand.grid(ZTemp = seq(-10,1.11,0.1), ZDistance = c(-1.094, -0.318, 2.004), 
                            Region = c("Europe",  "Oceania", "USA"), RegionXSiteName = Waj)
new_df
predictions
?fitted
predictions  <- fitted(mem_mixed0, se.fit = T)
new_df$mpg   <- predictions$fit
new_df$upper <- new_df$mpg + 1.96 * predictions$se.fit
new_df$lower <- new_df$mpg - 1.96 * predictions$se.fit
new_df$disp  <- factor(new_df$disp)

ggplot(new_df, aes(hp, mpg)) +
  geom_ribbon(aes(ymax = upper, ymin = lower, fill = disp), alpha = 0.3) +
  geom_line(aes(color = disp)) +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1")

ggplot(FigureDataDist, aes(x=TempBack, y=DistanceWithinFit, group = DistanceThresh, color = DistanceThresh)) +
  geom_smooth(method = lm)

?stat_smooth

### nothing works!!!!


#  Within Basin figures: complete ----------------------------------------------------------
library(wesanderson)
# install.packages("ggsci")
library("ggsci")

set_theme(base = theme_classic(), #To remove the background color and the grids
          theme.font = 'serif',   #To change the font type
          axis.title.size = 1.3,  #To change axis title size
          axis.textsize.x = 1.2,    #To change x axis text size
          # axis.angle.x = 60,      #To change x axis text angle
          # axis.hjust.x = 1,
          # axis.ticksmar = 50,
          axis.textsize.y = 1)  #To change y axis text size

theme_set(theme_sjplot())

## fitted values add to df
allsyncxWithinx <- allsyncxWithin %>%
  drop_na(WCDistkmsqrt) ## remove WC NAs

## checking distance data
testD <- allsyncxWithinx %>%
  # filter(ZDistance <= 0) %>% ## more similar
  filter(ZDistance >= 0) %>% ## more different
  select(Region, Sync, distance, annual_avgC, ZDistance) #%>%
  # pivot_longer(distance:annual_avgC)

ggplot(data = testD, aes(x = distance, y = Sync)) +
  stat_smooth(method = "lm") +
  geom_point() +
  facet_wrap(~Region)

ggplot(data = testD, aes(x = distance)) +
  geom_histogram()+
  facet_wrap(~Region, scales = "free")

## checking overlap data
testO <- allsyncxWithinx %>%
  filter(ZOverlapScaled >= 0) %>% ## more similar
  # filter(ZOverlap <= 0) %>% ## more different
  select(Region, Sync, overlapScaled, annual_avgC, ZDistance) #%>%
# pivot_longer(distance:annual_avgC)

ggplot(data = testO, aes(x = overlapScaled, y = Sync)) +
  stat_smooth(method = "lm") +
  facet_wrap(~Region)

ggplot(data = testD, aes(x = distance)) +
  geom_histogram()+
  facet_wrap(~Region, scales = "free")


# build formula to get  from allsyncxWithin$ZTempCor
orginal_values <- allsyncxWithinx$ZTempCor
original_mean <- mean(allsyncxWithin$ZTempCor) ## original mean
original_sd <- sd(allsyncxWithin$ZTempCor) ## original sd


target_values <- allsyncxWithin$annual_avgC
target_mean <- mean(allsyncxWithin$annual_avgC) ## target mean
target_sd <- sd(allsyncxWithin$annual_avgC) ## target sd
target_sd
target_mean
## find the scale factor 

scFact <- 1/original_sd ## 1.000206
scFact


## test transformation
test <- (scFact*(orginal_values-original_mean)) *target_sd + target_mean
test2 <- (target_sd*(orginal_values-original_mean)/original_sd+target_mean)
test
sd(test) ## correct!
mean(test) ## correct!

# cor(target_values, test) ## 0.633
## values don't exactly match, but mean and sd are correct

## try on figure values
df$x
test <- (scFact*(df$x-original_mean)) *target_sd + target_mean
test ## sync predctions going over 1, need to remove or change z scores. can only use integers but ztemp cor values should stop ~0.05

## prioginal values of overlap

round(quantile(allsyncxWithin$overlapScaled, probs = c(0.05,0.5,0.95)), digits = 3)
#    5%   50%   95% 
# 0.000 0.459 0.856 
hist(allsyncxWithin$overlapScaled)

test <- 1 + allsyncxWithin$ZOverlapScaled -min(allsyncxWithin$ZOverlapScaled)
hist(test)
min(allsyncxWithin$ZOverlapScaled)
min(test)

round(quantile(allsyncxWithin$distance, probs = c(0.05,0.5,0.95)), digits = 3)
#    5%   50%   95% 
# 0.001 0.018 0.074 
hist(allsyncxWithin$distance)
hist(allsyncxWithin$ZDistance)

## predict the values of sync with the interactions
quantile(allsyncxWithin$ZTempCor, probs = c(0.05,0.5,0.95))
round(quantile(allsyncxWithin$ZOverlapScaled, probs = c(0.05,0.5,1)), digits = 3)
round(quantile(allsyncxWithin$ZDistance, probs = c(0.05,0.5,0.95)), digits = 3)

## overlap interaction
df <- ggpredict(mem_mixed1, terms= c("ZTempCor [-10:1]", "ZOverlapScaled [-1.780, 0.048, 2.119]"))
?ggpredict
## transform the env sync values using above formula
df$x2 <- (scFact*(df$x-original_mean)) *target_sd + target_mean ## back transform values
df$group
cols <- c("#999999", "#E69F00","#56B4E9")

## plot using ggplot 
DivTemp  <- ggplot(df, aes(x2, predicted)) + 
  geom_line(aes(linetype=group, color=group, linewidth = 0)) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill=group), alpha=0.15) +
  scale_color_manual(values = cols, labels = c("Low", "Median", "High"), name = "Community Similarity") +
  #scale_color_manual(values= wes_palette("Zissou1", n = 3), labels = c("5th", "50th", "95th")) +
  scale_fill_manual(values = cols,labels = c("Low", "Median", "High"), name = "Community Similarity") +
  scale_linetype_manual(values = c("dotted", "dashed", "solid"), labels = c("Low", "Median", "High"), name = "Community Similarity") +
  scale_x_continuous(name = "Temperature Synchrony", limits = c(0.82,1)) +
  scale_y_continuous(name = "Thermal Synchrony") + theme(axis.text=element_text(size=20,face="bold"),
                                                               axis.title=element_text(size = 20,face="bold"), 
                                                               legend.text=element_text(size=20, face="bold"), 
                                                               legend.title=element_text(size = 20,face="bold")) +
  guides(linewidth = "none")
  
DivTemp 

## save out
file.name1 <- paste0(out.dir, "08_biotic_interactions_overlap_all_species_in_10_year_corTemp_Jan2024.jpg")
ggsave(DivTemp, filename=file.name1, dpi=300, height=12, width=18)

## overlap sungle effects

# build formula 
orginal_valuesO <- allsyncxWithinx$ZOverlapScaled
original_meanO <- mean(allsyncxWithin$ZOverlapScaled) ## original mean
original_sdO <- sd(allsyncxWithin$ZOverlapScaled) ## original sd

target_valuesO <- allsyncxWithin$overlapScaled
target_meanO <- mean(allsyncxWithin$overlapScaled) ## target mean
target_sdO <- sd(allsyncxWithin$overlapScaled) ## target sd


## scale factor
scFactO <- 1/original_sdO ## 1.000206
scFactO

dfT <- as.data.frame(ggpredict(mem_mixed1, terms= c("ZTempCor [-10:1]")))
dfO <- as.data.frame(ggpredict(mem_mixed1, terms= c("ZOverlapScaled")))

## transform the env sync values using above formula
dfT$x2 <- (scFact*(dfT$x-original_mean)) *target_sd + target_mean ## back transform values

## transform the overlap values using above formula
dfO$x2 <- (scFactO*(dfO$x-original_meanO)) *target_sdO + target_meanO ## back transform values

## combine data
dfT$Type <- "Temperature Synchrony"
dfO$Type <- "Thermal Overlap"

df1 <- bind_rows(dfT, dfO)
df1
palette.colors()
## plot using ggplot 
DivTempO  <- ggplot(df1, aes(x2, predicted)) + 
  geom_line(aes( color=Type, linewidth = 0)) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill=Type), alpha=0.15) +
  scale_color_manual(values = c("#009E73", "#CC79A7")) +
  # #scale_color_manual(values= wes_palette("Zissou1", n = 3), labels = c("5th", "50th", "95th")) +
  scale_fill_manual(values = c("#009E73", "#CC79A7")) +
  facet_wrap(~Type, scales = "free") +
  # scale_linetype_manual(values = c("solid", "dashed", "dotted"), labels = c("Low similarity", "Median Similarity", "High Similarity"), name = "") +
  scale_x_continuous(name = "") +
  scale_y_continuous(name = "Thermal Synchrony") + 
  theme(axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size = 20,face="bold"),
        legend.text=element_text(size=20, face="bold"),
        legend.title=element_text(size = 20,face="bold"),
        strip.text = element_text(size = 20,face="bold")) +
  guides(linewidth = "none") +
  theme(legend.position = "bottom", legend.title = element_blank())

DivTempO

## save out
file.name1 <- paste0(out.dir, "08_overlap_all_species_in_10_year_corTemp_single_effects.jpg")
ggsave(DivTempO, filename=file.name1, dpi=300, height=12, width=18)

## distance
df <- ggpredict(mem_mixed0, terms= c("ZTempCor [-10:1]", "ZDistance [-1.013,  -0.310, 2.034]"))
## transform the env sync values using above formula
df$x2 <- (scFact*(df$x-original_mean)) *target_sd + target_mean ## back transform values
cols

df$group
## plot using ggplot 
DistTemp  <- ggplot(df, aes(x2, predicted)) + 
  geom_line(aes(linetype=group, color=group, linewidth = 0)) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill=group), alpha=0.15) +
  # scale_fill_discrete(breaks=c(2.034, -0.31, -1.013)) +
  scale_color_manual(breaks=c(2.034, -0.31, -1.013), values = cols, labels = c("Low", "Median", "High"), name = "Community Similarity") +
  #scale_color_manual(values= wes_palette("Zissou1", n = 3), labels = c("5th", "50th", "95th")) +
  scale_fill_manual(breaks=c(2.034, -0.31, -1.013), values = cols, labels = c("Low", "Median", "High"), name = "Community Similarity") +
  scale_linetype_manual(breaks=c(2.034, -0.31, -1.013), values = c("dotted", "dashed", "solid"), labels = c("Low", "Median", "High"), name = "Community Similarity") +
  scale_x_continuous(name = "Temperature Synchrony", limits = c(0.82,1)) +
  scale_y_continuous(name = "Thermal Synchrony") + theme(axis.text=element_text(size=20,face="bold"),
                                                             axis.title=element_text(size = 20,face="bold"), 
                                                             legend.text=element_text(size=20, face="bold"), 
                                                             legend.title=element_text(size = 20,face="bold")) +
  guides(linewidth = "none")


DistTemp 

## save out
file.name1 <- paste0(out.dir, "08_biotic_interactions_dist_corTemp_Jan2024.jpg")
ggsave(DistTemp, filename=file.name1, dpi=300, height=12, width=18)

## distance sungle effects

# build formula 
orginal_valuesD <- allsyncxWithinx$ZDistance
original_meanD <- mean(allsyncxWithin$ZDistance) ## original mean
original_sdD <- sd(allsyncxWithin$ZDistance) ## original sd

target_valuesD <- allsyncxWithin$distance
target_meanD <- mean(allsyncxWithin$distance) ## target mean
target_sdD <- sd(allsyncxWithin$distance) ## target sd


## scale factor
scFactD <- 1/original_sdD ## 1.000206
scFactD

dfT <- as.data.frame(ggpredict(mem_mixed0, terms= c("ZTempCor [-10:1]")))
dfD <- as.data.frame(ggpredict(mem_mixed0, terms= c("ZDistance")))

## transform the env sync values using above formula
dfT$x2 <- (scFact*(dfT$x-original_mean)) *target_sd + target_mean ## back transform values

## transform the overlap values using above formula
dfD$x2 <- (scFactD*(dfD$x-original_meanD)) *target_sdD + target_meanD ## back transform values

## combine data
dfT$Type <- "Temperature Synchrony"
dfD$Type <- "Thermal Distance"

df2 <- bind_rows(dfT, dfD, dfO)
df2

## plot using ggplot 
DisTempAll  <- ggplot(df2, aes(x2, predicted)) + 
  geom_line(aes( color=Type, linewidth = 0)) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill=Type), alpha=0.15) +
  scale_color_manual(values = c("#009E73", "#CC79A7", "#F0E442")) +
  # #scale_color_manual(values= wes_palette("Zissou1", n = 3), labels = c("5th", "50th", "95th")) +
  scale_fill_manual(values = c("#009E73", "#CC79A7", "#F0E442")) +
  facet_wrap(~Type, scales = "free_x") +
  # scale_linetype_manual(values = c("solid", "dashed", "dotted"), labels = c("Low similarity", "Median Similarity", "High Similarity"), name = "") +
  scale_x_continuous(name = "") +
  scale_y_continuous(name = "Thermal Synchrony") + 
  theme(axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size = 20,face="bold"),
        legend.text=element_text(size=20, face="bold"),
        legend.title=element_text(size = 20,face="bold"),
        strip.text = element_text(size = 20,face="bold")) +
  guides(linewidth = "none") +
  theme(legend.position = "none", legend.title = element_blank())

DisTempAll

## save out
file.name1 <- paste0(out.dir, "08_distance_overlap_all_species_in_10_year_corTemp_single_effects.jpg")
ggsave(DisTempAll, filename=file.name1, dpi=300, height=12, width=18)

## join together
# WithinPlot <- cowplot::plot_grid(DistTemp, DivTemp, labels="auto", align = "hv") ## grid figures
# WithinPlot
# 
# 
# ## save out
# file.name1 <- paste0(out.dir, "08_biotic_interactions_dist_overlap_corTemp_Jan2024.jpg")
# ggsave(WithinPlot, filename=file.name1, dpi=300, height=12, width=18)

## estimates

## overlap
# estsow <- sjPlot::plot_model(mem_mixed1, 
#                              show.values=TRUE, show.p=TRUE,
#                              title="") +
#   scale_y_continuous(limits = c(-0.02, 0.02)) +
#   theme(text = element_text(size=25))
# 
# estsow
# 
# ## distance
# 
# estsdw <- sjPlot::plot_model(mem_mixed0, 
#                              show.values=TRUE, show.p=TRUE,
#                              title="") +
#   scale_y_continuous(limits = c(-0.02, 0.02))+
#   theme(text = element_text(size=25))
# 
# estsdw

### join together by biotic variable

## distance

DistPlot <- cowplot::plot_grid(estsdw, DistTemp, labels="auto") ## grid figures
DistPlot

## save out
file.name1 <- paste0(out.dir, "08_distance_ests_rels_corTemp.jpg")
ggsave(DistPlot, filename=file.name1, dpi=600, height=15, width=25)

## overlap
DivPlot <- cowplot::plot_grid(estsow, DivTemp, labels="auto") ## grid figures
DivPlot

## save out
file.name1 <- paste0(out.dir, "08_overlap_ests_rels_corTemp.jpg")
ggsave(DivPlot, filename=file.name1, dpi=600, height=15, width=20)


# Between basin plots: Appendix -------------------------------------------

## back transform zscores
## temp
mean.vari.temp <- mean(allsyncxBetween$annual_avg)
sd.vari.temp <- sd(allsyncxBetween$annual_avg)
var.back.temp <- (allsyncxBetween$ZTemp * sd.vari.temp) + mean.vari.temp

## distance
mean.vari.dist <- mean(allsyncxBetween$distance)
sd.vari.dist <- sd(allsyncxBetween$distance)
var.back.dist <- (allsyncxBetween$ZDistance * sd.vari.dist) + mean.vari.dist

## overlap/diversity
mean.vari.div <- mean(allsyncxBetween$overlap)
sd.vari.div <- sd(allsyncxBetween$overlap)
var.back.div <- (allsyncxBetween$ZOverlap * sd.vari.div) + mean.vari.div

## set background
set_theme(base = theme_classic(), #To remove the background color and the grids
          theme.font = 'serif',   #To change the font type
          axis.title.size = 1.3,  #To change axis title size
          axis.textsize.x = 1.2,    #To change x axis text size
          # axis.angle.x = 60,      #To change x axis text angle
          # axis.hjust.x = 1,
          # axis.ticksmar = 50,
          axis.textsize.y = 1)  #To change y axis text size

theme_set(theme_sjplot())

## Between basin distance
DisTempx <- sjPlot::plot_model(mem_mixed3, type="pred", terms= c("ZTemp", "ZDistance [-1.210, -0.199, 1.957]"),
                               axis.title = c("Environmental synchrony (Z Score)" ,"Thermal Trait Synchrony"), 
                               legend.title = "Trait Distance Percentile", title = "") + 
  scale_color_discrete(labels = c("5th", "50th", "95th")) + 
  # theme(legend.position = "none") +
  scale_x_continuous( sec.axis = sec_axis(~.*sd.vari+mean.vari, name="Environmental synchrony (Original)" ))+
  theme(text = element_text(size=25))

DisTempx

## Between basin overlap
DivTemp <- sjPlot::plot_model(mem_mixed2, type="pred", terms= c("ZTemp", "ZOverlap [-0.616, -0.601, 2.393]"),
                              axis.title = c("Environmental synchrony (Z Score)" ,"Thermal Trait Synchrony"), 
                              legend.title = "Trait Overlap Percentile", title = "") + 
  scale_color_discrete(labels = c("5th", "50th", "95th")) +
  scale_x_continuous( sec.axis = sec_axis(~.*sd.vari.temp+mean.vari.temp, name="Environmental synchrony (Original)" ))+
  theme(text = element_text(size=25))

DivTemp

## join together
BetweenPlot <- cowplot::plot_grid(DisTempx, DivTemp, labels="auto") ## grid figures
BetweenPlot

## save out
file.name1 <- paste0(out.dir, "08_biotic_interactions_dist_overlap_between.jpg")
ggsave(BetweenPlot, filename=file.name1, dpi=300, height=12, width=18)

## estimates

## overlap
estsow <- sjPlot::plot_model(mem_mixed2, 
                             show.values=TRUE, show.p=TRUE,
                             title="") +
  scale_y_continuous(limits = c(-0.02, 0.02))+
  theme(text = element_text(size=25))

estsow

## distance

estsdw <- sjPlot::plot_model(mem_mixed3, 
                             show.values=TRUE, show.p=TRUE,
                             title="") +
  scale_y_continuous(limits = c(-0.02, 0.02))+
  theme(text = element_text(size=25))

estsdw

### join together by biotic variable

## distance

DistPlot <- cowplot::plot_grid(estsdw, DisTempx, labels="auto") ## grid figures
DistPlot

## save out
file.name1 <- paste0(out.dir, "08_distance_ests_rels_between.jpg")
ggsave(DistPlot, filename=file.name1, dpi=300, height=15, width=20)

## overlap
DivPlot <- cowplot::plot_grid(estsow, DivTemp, labels="auto") ## grid figures
DivPlot

## save out
file.name1 <- paste0(out.dir, "08_overlap_ests_rels_between.jpg")
ggsave(DivPlot, filename=file.name1, dpi=300, height=15, width=20)


