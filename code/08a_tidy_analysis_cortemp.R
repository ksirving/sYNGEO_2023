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
library(ggpubr)


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

## summary stats per region

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
                       labels = c("Temp Synchrony", "Thermal Distance", "Thermal Overlap", "Water Course (km)" )))

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


## fitted values add to df
allsyncxWithinx <- allsyncxWithin %>%
  drop_na(WCDistkmsqrt) ## remove WC NAs

class(mem_mixedwc) <- "lmerModLmerTest"

allsyncxWithinx$mem_mixedwc <- fitted(mem_mixedwc)

set_theme(base = theme_classic(), #To remove the background color and the grids
          theme.font = 'sans',   #To change the font type
          axis.title.size = 1.3,  #To change axis title size
          axis.textsize.x = 1.2,    #To change x axis text size
          # axis.angle.x = 60,      #To change x axis text angle
          # axis.hjust.x = 1,
          # axis.ticksmar = 50,
          axis.textsize.y = 1)  #To change y axis text size


S2 <- ggplot(allsyncxWithinx, aes(x = WCDistkmsqrt, y = mem_mixedwc)) +
  geom_point(aes(y=Sync, col = Region), size = 0.01) +
  geom_smooth(method = "lm",color = "grey20", size=0.5) +
  scale_y_continuous(name="Thermal Synchrony") +
  scale_x_continuous(name="Water Course Distance (sqrt KM)") +
  labs(subtitle = "R² = 0.15, ICC = 0.12") +
  theme(text = element_text( size = 15)) + ## family = "Helvetica"
  theme(legend.key=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=4))) 
  
S2
file.name1 <- paste0(out.dir, "08_Fitted_Sync_Over_WCDistance_Within_sqrt.jpg")
ggsave(S2, filename=file.name1, dpi=300, height=5, width=8)

allsyncxWithinx
mem_mixedwc
summary(mem_mixedwc)
range(allsyncxWithinx$WCDistkm)
## get marginal effect of decrease in y with KMs

# Define the fixed effect coefficient from your model
beta1 <- -0.002105

# Choose distances (in km) where you want to evaluate marginal effects
km_values <- c(1, 5, 10, 20, 50, 100, 200, 500,1000,2000)

# Calculate marginal effect of a 1 km increase at each distance
marginal_effects <- sapply(km_values, function(x) {
  delta_y <- beta1 * (sqrt(x + 1) - sqrt(x))
  return(delta_y)
})

# Combine into a data frame for easy viewing
results <- data.frame(
  km_start = km_values,
  km_end = km_values + 1,
  marginal_effect = marginal_effects
)

# Show the results
print(results)

## all effect sizes
marginal_effects <- beta1 * (sqrt(allsyncxWithinx$WCDistkm + 1) - sqrt(allsyncxWithinx$WCDistkm))
marginal_effects

# Calculate the average marginal effect across all km values
average_effect <- mean(marginal_effects, na.rm = TRUE)


# Print result
cat("Average marginal effect per km:", round(average_effect, 6), "units of y\n")

# Biotic Variables: mixed membership model --------------------------------

## distance: Within Basin

Wa <- lmerMultiMember::weights_from_vector(allsyncxWithin$Region)
Wj <- Matrix::fac2sparse(allsyncxWithin$SiteName)  # convert single membership vars to an indicator matrix with fac2sparse()
Waj <- interaction_weights(Wa, Wj)

mem_mixed0 <- lmerMultiMember::lmer(Sync ~  ZTempCor*ZDistance
                                    + (1 | Region ) + ## add predictors here to get random effect per region
                                      (1 | RegionXSiteName), 
                                    memberships = list(Region = Wa, RegionXSiteName = Waj), 
                                    REML = T,
                                    data = allsyncxWithin)

summary(mem_mixed0, ddf = "Satterthwaite")
anova(mem_mixed0, ddf = "Satterthwaite")
r2_nakagawa(mem_mixed0) ## 0.145
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

## plot estimates
estsdw <- sjPlot::plot_model(mem_mixed0, 
                             type = "est",
                             axis.lim = c(-0.05, 0.05),
                             vline.color = "grey30",
                             show.values=T, show.p=TRUE,
                             axis.labels = c("Interaction", "Thermal Distance", "Temperature Synchrony"),
                             title = "")

estsdw <- estsdw + scale_y_continuous(limits = c(-0.025, 0.025))
estsdw

# Extract the plot data
plot_data <- ggplot_build(estsdw)$data[[4]]  # Typically the layer with text (values + stars)
plot_data

## Retain only stars, remove numbers, and set to NA if no stars
plot_data$label2 <- sub("^.*?(\\*+)$", "\\1", plot_data$label)  # Keep only the stars
plot_data$label2 <- ifelse(grepl("\\*", plot_data$label2), plot_data$label2, NA)  # Set to NA if no stars

set_theme(base = theme_classic(), #To remove the background color and the grids
          theme.font = 'sans',   #To change the font type
          axis.title.size = 1.3,  #To change axis title size
          axis.textsize.x = 1.2,    #To change x axis text size
          # axis.angle.x = 60,      #To change x axis text angle
          # axis.hjust.x = 1,
          # axis.ticksmar = 50,
          axis.textsize.y = 1)  #To change y axis text size


# Ensure the group variable is discrete
if ("group" %in% names(plot_data)) {
  plot_data$group <- as.factor(plot_data$group)  # Convert group to a factor if necessary
  plot_data$group <- droplevels(plot_data$group)
}

## plot again without values and ps
plot <- sjPlot::plot_model(mem_mixed0, 
                           type = "est",
                           axis.lim = c(-0.05, 0.05),
                           vline.color = "grey",
                           show.values=F, show.p=F,
                           # value.size = 0,
                           axis.labels = c("Interaction", "Thermal Distance", "Temperature Synchrony"),
                           title = "",
                           axis.title = c("Beta-coefficient", "")) 

plot <- plot + 
  scale_y_continuous(limits = c(-0.025, 0.025)) 
  
plot
# Re-plot with updated text labels
customized_distance <- plot +
  # labs(x = "", y = "Beta-coefficient") +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        legend.text=element_text(size=15), 
        legend.title=element_text(size=15))+
  geom_text(
    data = plot_data, aes(x = x, y = y, label = label2),
    size = 4, hjust = -0.5, inherit.aes = FALSE
  ) 

customized_distance


file.name1 <- paste0(out.dir, "08_withinBasin_effect_sizes_all_expsync_n4_sites_distance_zscores_cortemp.jpg")
ggsave(customized_distance, filename=file.name1, dpi=300, height=8, width=10)


## Overlap: Within Basin

mem_mixed1 <- lmerMultiMember::lmer(Sync ~  ZTempCor*ZOverlapScaled
                                    + (1 | Region ) + ## add predictors here to get random effect per region
                                      (1 | RegionXSiteName),
                                    memberships = list(Region = Wa, RegionXSiteName = Waj),
                                    REML = T,
                                    data = allsyncxWithin)


summary(mem_mixed1, ddf = "Satterthwaite")

anova(mem_mixed1, ddf = "Satterthwaite")
r2_nakagawa(mem_mixed1)[1] ## 0.15
r2 <- r2_nakagawa(mem_mixed1)[1]
r2$R2_conditional
check_singularity(mem_mixed1) ## False
performance::icc(mem_mixed1, by_group = T)[1,2] ## 0.13

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

## estimates
estsow <- sjPlot::plot_model(mem_mixed1, 
                             type = "est",
                             axis.lim = c(-0.05, 0.05),
                             vline.color = "grey",
                             show.values=T, show.p=T,
                             # value.size = 0,
                             axis.labels = c("Interaction", "Thermal Overlap", "Temperature Synchrony"),
                             title = "",
                             axis.title = c("Beta-coefficient", "")) 

estsow  <- estsow + scale_y_continuous(limits = c(-0.025, 0.025))
estsow

# Extract the plot data
plot_data <- ggplot_build(estsow)$data[[4]]  # Typically the layer with text (values + stars)
plot_data

## Retain only stars, remove numbers, and set to NA if no stars
plot_data$label2 <- sub("^.*?(\\*+)$", "\\1", plot_data$label)  # Keep only the stars
plot_data$label2 <- ifelse(grepl("\\*", plot_data$label2), plot_data$label2, NA)  # Set to NA if no stars


# Ensure the group variable is discrete
if ("group" %in% names(plot_data)) {
  plot_data$group <- as.factor(plot_data$group)  # Convert group to a factor if necessary
  plot_data$group <- droplevels(plot_data$group)
}

## plot again without values and ps
plotO <- sjPlot::plot_model(mem_mixed1, 
                             type = "est",
                             axis.lim = c(-0.05, 0.05),
                             vline.color = "grey",
                             show.values=F, show.p=F,
                             # value.size = 0,
                             axis.labels = c("Interaction", "Thermal Overlap", "Temperature Synchrony"),
                             title = "",
                            axis.title = c("Beta-coefficient", "")) 
plotO
plotO <- plotO + scale_y_continuous(limits = c(-0.025, 0.025))

# Re-plot with updated text labels
customized_overlap<- plotO +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        legend.text=element_text(size=15), 
        legend.title=element_text(size=15))+
  geom_text(
    data = plot_data, aes(x = x, y = y, label = label2),
    size = 4, hjust = -0.5, inherit.aes = FALSE
  ) 

customized_overlap
customized_distance

file.name1 <- paste0(out.dir, "08_withinBasin_effect_sizes_all_expsync_n4_sites_overlap_zscores_corTemp.jpg")
ggsave(customized_overlap, filename=file.name1, dpi=300, height=8, width=10)


### combine estimates plots
combined_plot <- ggarrange(customized_distance,
                           customized_overlap,
                           nrow = 2,
                           ncol = 1,
                           labels = "AUTO")
combined_plot

## save
file.name1 <- paste0(out.dir, "08_all_estimates.jpg")
ggsave(combined_plot, filename=file.name1, dpi=600, height=8, width=10)

#  Within Basin figures: complete ----------------------------------------------------------
library(wesanderson)
# install.packages("ggsci")
library("ggsci")

set_theme(base = theme_classic(), #To remove the background color and the grids
          theme.font = 'sans',   #To change the font type
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

## target
target_values <- allsyncxWithin$annual_avgC
target_mean <- mean(allsyncxWithin$annual_avgC) ## target mean
target_sd <- sd(allsyncxWithin$annual_avgC) ## target sd

## find the scale factor 

scFact <- 1/original_sd ## 1.000206
scFact


## test transformation
test <- (scFact*(orginal_values-original_mean)) *target_sd + target_mean
sd(test) ## correct!
mean(test) ## correct!

## predict the values of sync with the interactions
quantile(allsyncxWithin$ZTempCor, probs = c(0.05,0.5,0.95))
round(quantile(allsyncxWithin$ZOverlapScaled, probs = c(0.05,0.5,1)), digits = 3)
round(quantile(allsyncxWithin$ZDistance, probs = c(0.05,0.5,0.95)), digits = 3)

## overlap interaction
dfO <- ggpredict(mem_mixed1, terms= c("ZTempCor [-10:1]", "ZOverlapScaled [-1.780, 0.048, 2.119]"))
dfO
## transform the env sync values using above formula
dfO$x2 <- (scFact*(dfO$x-original_mean)) *target_sd + target_mean ## back transform values
dfO$group
cols <- c( "#56B4E9","#E69F00","#999999") 
dfO
## plot using ggplot 
DivTemp  <- ggplot(dfO, aes(x2, predicted)) + 
  geom_line(aes(linetype=group, color=group, linewidth = 0)) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill=group), alpha=0.15) +
  scale_color_manual(values = cols, labels = c(  "High","Median", "Low"), name = "Thermal Dissimilarity") +
  #scale_color_manual(values= wes_palette("Zissou1", n = 3), labels = c("5th", "50th", "95th")) +
  scale_fill_manual(values = cols,labels = c( "High","Median", "Low"), name = "Thermal Dissimilarity") +
  scale_linetype_manual(values = c("solid", "dashed", "dotted"), labels = c( "High","Median", "Low"), name = "Thermal Dissimilarity") +
  scale_x_continuous(name = "Temperature Synchrony", limits = c(0.82,1)) +
  scale_y_continuous(name = "Thermal Synchrony") + theme(axis.text=element_text(size=20,face="bold"),
                                                               axis.title=element_text(size = 20,face="bold"), 
                                                               legend.text=element_text(size=20, face="bold"), 
                                                               legend.title=element_text(size = 20,face="bold")) +
  guides(linewidth = "none") +
  guides(color = guide_legend(override.aes = list(size=10))) 
  
DivTemp  ## overlap - low values = high dissimilarity
## save out
file.name1 <- paste0(out.dir, "08_biotic_interactions_overlap_all_species_in_10_year_corTemp_Jan2024.jpg")
ggsave(DivTemp, filename=file.name1, dpi=300, height=12, width=18)

## overlap single effects

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
dfT$x2 <- (scFact*(dfT$x-original_mean)) *target_sd + target_mean ## back transform temp values

## transform the overlap values using above formula
dfO$x2 <- (scFactO*(dfO$x-original_meanO)) *target_sdO + target_meanO ## back transform overlap values

## combine data
dfT$Type <- "Temperature Synchrony"
dfO$Type <- "Thermal Overlap"

df1 <- bind_rows(dfT, dfO)

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
  theme(text = element_text( size = 15)) + ## family = "Helvetica"
  guides(linewidth = "none") +
  theme(legend.position = "bottom", legend.title = element_blank())

DivTempO

## save out
file.name1 <- paste0(out.dir, "08_overlap_all_species_in_10_year_corTemp_single_effects.jpg")
ggsave(DivTempO, filename=file.name1, dpi=300, height=12, width=18)

## distance interations
dfd <- ggpredict(mem_mixed0, terms= c("ZTempCor [-10:1]", "ZDistance [-1.013,  -0.310, 2.034]"))
dfd$group
## transform the env sync values using above formula
dfd$x2 <- (scFact*(dfd$x-original_mean)) *target_sd + target_mean ## back transform values

cols2 <- c("#999999",  "#E69F00", "#56B4E9") ## grey = high similarity/low dissimilarity, blue = low similarity/high dissimilarity

## plot using ggplot 
DistTemp  <- ggplot(dfd, aes(x2, predicted)) + 
  geom_line(aes(linetype=group, color=group, linewidth = 0)) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill=group), alpha=0.15) +
  # scale_fill_discrete(breaks=c(2.034, -0.31, -1.013)) +
  scale_color_manual(values = cols2, labels = c("Low", "Median", "High"), name = "Thermal Dissimilarity") +
  #scale_color_manual(values= wes_palette("Zissou1", n = 3), labels = c("5th", "50th", "95th")) +
  scale_fill_manual(values = cols2, labels = c("Low", "Median", "High"), name = "Thermal Dissimilarity") +
  scale_linetype_manual( values = c( "dotted","dashed", "solid"), labels = c("Low", "Median", "High"), name = "Thermal Dissimilarity") +
  scale_x_continuous(name = "Temperature Synchrony", limits = c(0.82,1)) +
  scale_y_continuous(name = "Thermal Synchrony") + theme(axis.text=element_text(size=20,face="bold"),
                                                             axis.title=element_text(size = 20,face="bold"), 
                                                             legend.text=element_text(size=20, face="bold"), 
                                                             legend.title=element_text(size = 20,face="bold")) +
  guides(linewidth = "none")

#breaks=c(-1.013, -0.31, 2.034), 
DistTemp ## distance - high values = high dissimilarity - figure correct 
DivTemp ## overlap - low values = high dissimilarity

## save out
file.name1 <- paste0(out.dir, "08_biotic_interactions_dist_corTemp_Jan2024.jpg")
ggsave(DistTemp, filename=file.name1, dpi=300, height=12, width=18)

## distance single effects

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
dfD
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
WithinPlot <- cowplot::plot_grid(DistTemp, DivTemp, labels="auto", align = "hv") ## grid figures
WithinPlot


# Plot interactions together -----------------------------------------------

## overlap interaction
dfO <- ggpredict(mem_mixed1, terms= c("ZTempCor [-10:1]", "ZOverlapScaled [-1.780, 0.048, 2.119]"))
## transform the env sync values using above formula
dfO$x2 <- (scFact*(dfO$x-original_mean)) *target_sd + target_mean ## back transform values

## distance interations
dfd <- ggpredict(mem_mixed0, terms= c("ZTempCor [-10:1]", "ZDistance [-1.013,  -0.310, 2.034]"))
dfd$group
## transform the env sync values using above formula
dfd$x2 <- (scFact*(dfd$x-original_mean)) *target_sd + target_mean ## back transform values

dfO


## make into dataframe 
dfO_df <- as.data.frame(dfO) %>%
  mutate(ThermalNiche = "Overlap")

dfd_df <- as.data.frame(dfd) %>%
  mutate(ThermalNiche = "Distance")

DistTemp ## distance - high values = high dissimilarity  
DivTemp ## overlap - low values = high dissimilarity

## combine

dfAll <- bind_rows(dfd_df, dfO_df) %>%
  mutate(Labels = case_when(group %in% c(-1.78, 2.034) ~ "High",
                            group %in% c(-0.31, 0.048) ~ "Median",
                            group %in% c(-1.013, 2.119)~ "Low")) %>%
  mutate(Labels = factor(Labels, levels = c("Low", "Median", "High"))) %>%
  mutate(FacetNames = ifelse(ThermalNiche == "Distance", paste("a) Distance"), paste("b) Overlap")))

head(dfAll)

unique(dfAll$Labels)

set_theme(base = theme_classic(), #To remove the background color and the grids
          theme.font = 'sans')   #To change the font type
          # axis.title.size = 1.3,  #To change axis title size
          # axis.textsize.x = 1.2,    #To change x axis text size
          # axis.angle.x = 60,      #To change x axis text angle
          # axis.hjust.x = 1,
          # axis.ticksmar = 50,
          # axis.textsize.y = 1)  #To change y axis text size



cols <- c("#56B4E9", "#E69F00", "#999999")

unique(dfAll$group)

## plot using ggplot 
InterPlot  <- ggplot(dfAll, aes(x2, predicted)) + 
  geom_line(aes(linetype=Labels, color=Labels, linewidth = 0)) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill=Labels), alpha=0.15) +
  labs(subtitle = "Distance: R² = 0.15, ICC = 0.12,  Overlap: R² = 0.15, ICC = 0.13") +
  # scale_fill_discrete(breaks=c(2.034, -0.31, -1.013)) +
  facet_wrap(~FacetNames) +
  scale_color_manual(values = cols, name = "Thermal Dissimilarity") +
  #scale_color_manual(values= wes_palette("Zissou1", n = 3), labels = c("5th", "50th", "95th")) +
  scale_fill_manual(values = cols, name = "Thermal Dissimilarity") +
  scale_linetype_manual(values = c("solid", "dashed", "dotted"), name = "Thermal Dissimilarity") +
  scale_x_continuous(name = "Temperature Synchrony", limits = c(0.82,1)) +
  scale_y_continuous(name = "Thermal Synchrony", limits = c(0.2,0.7), breaks = seq(0.2,0.7, by = 0.2)) + 
  theme(axis.text=element_text(size=30),
    axis.title=element_text(size = 30), 
    legend.text=element_text(size=30), 
    legend.title=element_text(size = 30),
    strip.text = element_text(size = 30)) +
  theme(text = element_text( size = 30)) + ## family = "Helvetica"
  guides(linewidth = "none") +
  guides(color = guide_legend(override.aes = list(size=20))) 



InterPlot

## save out
file.name1 <- paste0(out.dir, "08_biotic_interactions_both.jpg")
ggsave(InterPlot, filename=file.name1, dpi=600, height=12, width=18)



