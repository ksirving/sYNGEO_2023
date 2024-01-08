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

load(file = "sync/04_temp_pref_env_dist_no_dupl_pairs_WC.RData")
head(allsyncx)

## sites with 3 species
load(file="output_data/01a_sites_3_species.RData")
head(singleSp3l)
length(singleSp3l) ## 151

## change names and filter single species sites

allsyncx <- allsyncx %>%
  rename(overlap = diversity2) %>%
  filter(!Site_ID1 %in% singleSp3l, !Site_ID2 %in% singleSp3l) #%>%

## remove NAs from temp column
allsyncx <- allsyncx %>%
  drop_na(annual_avg)

## directory for figures
out.dir <- "/Users/katieirving/OneDrive - SCCWRP/Documents - Katie’s MacBook Pro/git/sYNGEO_2023/Figures/Tidy/"

#  Make z scores and histograms ----------------------------------------------------------

allsyncx <- allsyncx %>%
  group_by(Connectivity, Region) %>%
  mutate(ZTemp = (annual_avg - mean(annual_avg))/sd(annual_avg), ## z scores
         ZDistance = (distance - mean(distance))/sd(distance),
         ZOverlap = (overlap - mean(overlap))/sd(overlap),
         WCDistkmsqrt = sqrt(WCDistkm),
         DistKMLog = log(DistKM+1),
         DistKMsqrt = sqrt(DistKM)) ## square root distances

## number of sites
length(unique(allsyncx$Site_ID2)) ## 613

hist(allsyncx$Sync)

hist(log(allsyncx$distance))
hist(allsyncx$distance)

hist(allsyncx$overlap)
hist(log(allsyncx$overlap))

hist(log(allsyncx$annual_avg))
hist(allsyncx$annual_avg)

# CORRELATIONS ------------------------------------------------------------

## all
allsyncxCor <- allsyncx %>%
  ungroup() %>%
  select(distance, overlap, annual_avg, DistKM)

syncCor <- cor(allsyncxCor)

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

# distance     overlap  annual_avg      DistKM    WCDistkm
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
  filter(Connectivity == "Within Basin") %>%
  select(ZDistance, ZOverlap, ZTemp, DistKM)

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


# Prep Data for models ----------------------------------------------------

## connectiovity as a factor and change name 
allsyncx$Connectivity <- as.factor(allsyncx$Connectivity)
allsyncx$Connectivity <- recode_factor(allsyncx$Connectivity,  "1" = "Within Basin", "0" = "Between Basin") 

# ## within basin - make longer to get site names for membership model
allsyncxWithin <- allsyncx %>%
  filter(Connectivity == "Within Basin") %>%
  pivot_longer(Site_ID1:Site_ID2, names_to = "SiteNumber", values_to = "SiteName") 

# ## between basin -  make longer to get site names for membership model
allsyncxBetween <- allsyncx %>%
  filter(Connectivity == "Between Basin") %>%
  pivot_longer(Site_ID1:Site_ID2, names_to = "SiteNumber", values_to = "SiteName") 


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

## fitted values add to df
allsyncxWithinx <- allsyncxWithin %>%
  drop_na(WCDistkmsqrt) ## remove WC NAs

class(mem_mixedwc) <- "lmerModLmerTest"

allsyncxWithinx$mem_mixedwc <- fitted(mem_mixedwc)

S2 <- ggplot(allsyncxWithinx, aes(x = WCDistkmsqrt, y = mem_mixedwc)) +
  geom_point(aes(y=Sync, col = Region), size = 0.01) +
  geom_smooth(method = "loess") +
  scale_y_continuous(name="Trait Synchrony") +
  scale_x_continuous(name="Water Course Distance (sqrt KM)")
S2

file.name1 <- paste0(out.dir, "08_Fitted_Sync_Over_WCDistance_Within_sqrt.jpg")
ggsave(S2, filename=file.name1, dpi=300, height=5, width=6)

## between basin - euclidean distance

## create memberships
Wa <- lmerMultiMember::weights_from_vector(allsyncxBetween$Region)
Wj <- Matrix::fac2sparse(allsyncxBetween$SiteName)  # convert single membership vars to an indicator matrix with fac2sparse()
Waj <- interaction_weights(Wa, Wj)

mem_mixedb <- lmerMultiMember::lmer(Sync ~ DistKMsqrt
                                     + (1 | Region ) + ## add predictors here to get random effect per region
                                       + (1 | RegionXSiteName), 
                                     memberships = list(Region = Wa, RegionXSiteName = Waj), 
                                     REML = T,
                                     data = allsyncxBetween)



## coefs
# summary(mem_mixedb, ddf = "Satterthwaite")
anova(mem_mixedb, ddf = "Satterthwaite")
r2_nakagawa(mem_mixedb) ## 0.15
check_singularity(mem_mixedb) ## False
performance::icc(mem_mixedb, by_group = T) ## 0.12

## fitted values add to df
allsyncxBetweenx <- allsyncxBetween %>%
  drop_na(DistKMsqrt) ## remove WC NAs

class(mem_mixedb) <- "lmerModLmerTest"

allsyncxBetweenx$mem_mixedb <- fitted(mem_mixedb)

S3 <- ggplot(allsyncxBetweenx, aes(x = DistKMsqrt, y = mem_mixedb)) +
  geom_point(aes(y=Sync, col = Region), size = 0.01) +
  geom_smooth(method = "loess", se = F) +
  scale_y_continuous(name="Trait Synchrony") +
  scale_x_continuous(name="Euclidean Distance (sqrt KM)")
S3

file.name1 <- paste0(out.dir, "08_Fitted_Sync_Over_EUCDistance_Between_sqrt.jpg")
ggsave(S3, filename=file.name1, dpi=300, height=5, width=6)

## combine plots
# spdecay <- plot_grid(list(S2, S3)) ## grid figures
spdecay <- cowplot::plot_grid(S2,S3, align="h")

file.name1 <- paste0(out.dir, "08_spatial_decay_within_between.jpg")
ggsave(spdecay, filename=file.name1, dpi=300, height=5, width=12)


# Biotic Variables: mixed membership model --------------------------------

## distance: Within Basin

Wa <- lmerMultiMember::weights_from_vector(allsyncxWithin$Region)
Wj <- Matrix::fac2sparse(allsyncxWithin$SiteName)  # convert single membership vars to an indicator matrix with fac2sparse()
Waj <- interaction_weights(Wa, Wj)

mem_mixed0 <- lmerMultiMember::lmer(Sync ~  ZDistance*ZTemp
                                    + (1 | Region ) + ## add predictors here to get random effect per region
                                      (1 | RegionXSiteName), 
                                    memberships = list(Region = Wa, RegionXSiteName = Waj), 
                                    REML = T,
                                    data = allsyncxWithin)
allsyncxWithin$SiteName

summary(mem_mixed0, ddf = "Satterthwaite")
anova(mem_mixed0, ddf = "Satterthwaite")
r2_nakagawa(mem_mixed0) ## 0.151
check_singularity(mem_mixed0) ## False
performance::icc(mem_mixed0, by_group = T) ## 0.122


### plots 
class(mem_mixed0) <- "lmerModLmerTest"

estsdw <- sjPlot::plot_model(mem_mixed0, 
                           show.values=TRUE, show.p=TRUE,
                           title="Drivers of Thermal Synchrony")

estsdw
file.name1 <- paste0(out.dir, "08_withinBasin_effect_sizes_all_expsync_n4_sites_distance_zscores.jpg")
ggsave(estsdw, filename=file.name1, dpi=300, height=8, width=10)

# estsTab <- sjPlot::tab_model(mem_mixed0, 
#                              show.re.var= TRUE, 
#                              pred.labels =c("(Intercept)", "ZDistance", "ZTemp"),
#                              dv.labels= "Drivers of Thermal Synchrony")
# 
# estsTab

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
mean(allsyncxWithin$overlap)
quantile(allsyncxWithin$ZTemp, probs = c(0.05,0.5,0.95))
round(quantile(allsyncxWithin$ZOverlap, probs = c(0.05,0.5,0.95)), digits = 3)
round(quantile(allsyncxWithin$ZDistance, probs = c(0.05,0.5,0.95)), digits = 3)

## overlap v temp
tempDis <- sjPlot::plot_model(mem_mixed0, type="pred", terms= c("ZDistance", "ZTemp [-1.86, 0.3, 0.92]"),
                              axis.title = c( "Trait Distance (Z Score)","Thermal Trait Synchrony"), 
                              legend.title = "Environmental synchrony")
tempDis


DisTemp <- sjPlot::plot_model(mem_mixed0, type="pred", terms= c("ZTemp", "ZDistance [-1.094, -0.318, 2.004]"),
                              axis.title = c("Environmental synchrony (Z Score)" ,"Thermal Trait Synchrony"), 
                              legend.title = "Trait Distance")
DisTemp

allsyncxBetweenx$mem_mixedb <- fitted(mem_mixedb)

## Overlap: Within Basin
mem_mixed1 <- lmerMultiMember::lmer(Sync ~  ZOverlap*ZTemp
                                    + (1 | Region ) + ## add predictors here to get random effect per region
                                      (1 | RegionXSiteName), 
                                    memberships = list(Region = Wa, RegionXSiteName = Waj), 
                                    REML = T,
                                    data = allsyncxWithin)


summary(mem_mixed1, ddf = "Satterthwaite")
anova(mem_mixed1, ddf = "Satterthwaite")
r2_nakagawa(mem_mixed1) ## 0.153
check_singularity(mem_mixed1) ## False
performance::icc(mem_mixed1, by_group = T) ## 0.124


### plots 
class(mem_mixed1) <- "lmerModLmerTest"

estsow <- sjPlot::plot_model(mem_mixed1, 
                           show.values=TRUE, show.p=TRUE,
                           title="Drivers of Thermal Synchrony")

estsow
file.name1 <- paste0(out.dir, "08_withinBasin_effect_sizes_all_expsync_n4_sites_overlap_zscores.jpg")
ggsave(estsow, filename=file.name1, dpi=300, height=8, width=10)

# estsTab <- sjPlot::tab_model(mem_mixed1, 
#                              show.re.var= TRUE, 
#                              pred.labels =c("(Intercept)", "overlap", "annual avg"),
#                              dv.labels= "Drivers of Thermal Synchrony")
# 
# estsTab

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
round(quantile(allsyncxWithin$ZOverlap, probs = c(0.05,0.5,0.95)), digits = 3)
round(quantile(allsyncxWithin$ZDistance, probs = c(0.05,0.5,0.95)), digits = 3)

## overlap v temp
tempDiv <- sjPlot::plot_model(mem_mixed1, type="pred", terms= c("ZOverlap", "ZTemp [-1.86, 0.3, 0.92]"),
                              axis.title = c( "Trait Overlap (Z Score)","Thermal Trait Synchrony"), 
                              legend.title = "Environmental synchrony")
tempDiv


DivTemp <- sjPlot::plot_model(mem_mixed1, type="pred", terms= c("ZTemp", " ZOverlap [-1.436, 0.053, 1.499]"),
                              axis.title = c("Environmental synchrony (Z Score)" ,"Thermal Trait Synchrony"), 
                              legend.title = "Trait Overlap")
DivTemp



## create memberships
Wa <- lmerMultiMember::weights_from_vector(allsyncxBetween$Region)
Wj <- Matrix::fac2sparse(allsyncxBetween$SiteName)  # convert single membership vars to an indicator matrix with fac2sparse()
Waj <- interaction_weights(Wa, Wj)

## Overlap: Between basin
mem_mixed2 <- lmerMultiMember::lmer(Sync ~  ZOverlap*ZTemp
                                    + (1 | Region ) + ## add predictors here to get random effect per region
                                      (1 | RegionXSiteName), 
                                    memberships = list(Region = Wa, RegionXSiteName = Waj), 
                                    REML = T,
                                    data = allsyncxBetween)


summary(mem_mixed2, ddf = "Satterthwaite")
anova(mem_mixed2, ddf = "Satterthwaite")
r2_nakagawa(mem_mixed2) ## 0.148
check_singularity(mem_mixed2) ## False
performance::icc(mem_mixed2, by_group = T) ## 0.113

# lattice::qqmath(mem_mixed1)
# plot(mem_mixed1, type=c("p","smooth"), col.line=1)
# check_model(mem_mixed1)

### plots 
class(mem_mixed2) <- "lmerModLmerTest"
# sjPlot::plot_model(mem_mixed) 
estsob <- sjPlot::plot_model(mem_mixed2, 
                           show.values=TRUE, show.p=TRUE,
                           title="Drivers of Thermal Synchrony")

estsob
file.name1 <- paste0(out.dir, "08_betweenBasin_effect_sizes_all_expsync_n4_sites_overlap_zscores.jpg")
ggsave(estsob, filename=file.name1, dpi=300, height=8, width=10)

# estsTab <- sjPlot::tab_model(mem_mixed2, 
#                              show.re.var= TRUE, 
#                              pred.labels =c("(Intercept)", "overlap", "annual avg"),
#                              dv.labels= "Drivers of Thermal Synchrony")
# 
# estsTab

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
mean(allsyncxBetween$overlap)
quantile(allsyncxBetween$ZTemp, probs = c(0.05,0.5,0.95))
round(quantile(allsyncxBetween$ZOverlap, probs = c(0.05,0.5,0.95)), digits = 3)
round(quantile(allsyncxBetween$ZDistance, probs = c(0.05,0.5,0.95)), digits = 3)

## overlap v temp
tempDivB <- sjPlot::plot_model(mem_mixed2, type="pred", terms= c("ZOverlap", "ZTemp [-1.96, 0.16, 1.15]"),
                              axis.title = c( "Trait Overlap (Z Score)","Thermal Trait Synchrony"), 
                              legend.title = "Environmental synchrony")
tempDivB


DivTempB <- sjPlot::plot_model(mem_mixed2, type="pred", terms= c("ZTemp", " ZOverlap [-0.616, -0.601, 2.393]"),
                              axis.title = c("Environmental synchrony (Z Score)" ,"Thermal Trait Synchrony"), 
                              legend.title = "Trait Overlap")
DivTempB




## Distance: Between Basin

mem_mixed3 <- lmerMultiMember::lmer(Sync ~  ZDistance*ZTemp
                                    + (1 | Region ) + ## add predictors here to get random effect per region
                                      (1 | RegionXSiteName), 
                                    memberships = list(Region = Wa, RegionXSiteName = Waj), 
                                    REML = T,
                                    data = allsyncxBetween)


summary(mem_mixed3, ddf = "Satterthwaite")
anova(mem_mixed3, ddf = "Satterthwaite")
r2_nakagawa(mem_mixed3) ## 0.15
check_singularity(mem_mixed3) ## False
performance::icc(mem_mixed3, by_group = T) ## 0.12

# lattice::qqmath(mem_mixed1)
# plot(mem_mixed1, type=c("p","smooth"), col.line=1)
# check_model(mem_mixed1)

### plots 
class(mem_mixed3) <- "lmerModLmerTest"
# sjPlot::plot_model(mem_mixed) 
estsdb <- sjPlot::plot_model(mem_mixed3, 
                           show.values=TRUE, show.p=TRUE,
                           title="Drivers of Thermal Synchrony")

estsdb
file.name1 <- paste0(out.dir, "08_betweenBasin_effect_sizes_all_expsync_n4_sites_distance_zscores.jpg")
ggsave(estsdb, filename=file.name1, dpi=300, height=8, width=10)

# estsTab <- sjPlot::tab_model(mem_mixed3, 
#                              show.re.var= TRUE, 
#                              pred.labels =c("(Intercept)", "overlap", "annual avg"),
#                              dv.labels= "Drivers of Thermal Synchrony")
# 
# estsTab

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
mean(allsyncxWithin$overlap)
quantile(allsyncxBetween$ZTemp, probs = c(0.05,0.5,0.95))
round(quantile(allsyncxBetween$ZOverlap, probs = c(0.05,0.5,0.95)), digits = 3)
round(quantile(allsyncxBetween$ZDistance, probs = c(0.05,0.5,0.95)), digits = 3)

## overlap v temp
tempDisB <- sjPlot::plot_model(mem_mixed3, type="pred", terms= c("ZDistance", "ZTemp [-1.96, 0.16, 1.15]"),
                              axis.title = c( "Trait Distance (Z Score)","Thermal Trait Synchrony"), 
                              legend.title = "Environmental synchrony")
tempDisB


DisTempB <- sjPlot::plot_model(mem_mixed3, type="pred", terms= c("ZTemp", "ZDistance [-1.210, -0.199, 1.957]"),
                              axis.title = c("Environmental synchrony (Z Score)" ,"Thermal Trait Synchrony"), 
                              legend.title = "Trait Distance")
DisTempB


# Figures -----------------------------------------------------------------

## Distance
distplot <- plot_grid(list(estsdw, estsdb, DisTemp,DisTempB)) ## grid figures

file.name1 <- paste0(out.dir, "08_distance_zscores_between_within_ests_rels.jpg")
ggsave(distplot, filename=file.name1, dpi=300, height=12, width=18)

## Overlap

divplot <- plot_grid(list(estsow, estsob, DivTemp,DivTempB)) ## grid figures

file.name1 <- paste0(out.dir, "08_overlap_zscores_between_within_ests_rels.jpg")
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
  select(Sync, Pair, Region, distance, overlap, annual_avg, ZDistance, ZOverlap, ZTemp, TempBack,DistanceWithinFit:DistanceThresh) %>%
  mutate(DistanceThresh = ifelse(ZDistance < -1.094, "5thPercentile", NA)),
         DistanceThresh = ifelse(ZDistance > -0.4 & ZDistance < -0.320, "50thPercentile", DistanceThresh),
         DistanceThresh = ifelse(ZDistance > 2.004 , "95thPercentile", DistanceThresh)) %>%
  mutate(OverlapThresh = ifelse(ZOverlap < -1.436, "5thPercentile", NA),
         OverlapThresh = ifelse(ZOverlap > -0.8 & ZOverlap < 0.2, "50thPercentile", OverlapThresh),
         OverlapThresh = ifelse(ZOverlap > 1.499 , "95thPercentile", OverlapThresh))

FigureData <- allsyncxWithin %>%
  select(Sync, Pair, Region, distance, overlap, annual_avg, ZDistance, ZOverlap, ZTemp, TempBack,DistanceWithinFit:DistanceThresh) %>%
  mutate(DistancePercentiles = percent_rank(ZDistance)*100)

head(FigureData)  

FigureDataDist <- FigureData %>%
  drop_na(DistanceThresh)

min(allsyncxWithin$ZTemp)
max(allsyncxWithin$ZTemp)
?scale_x_continuous
DisTempx <- sjPlot::plot_model(mem_mixed0, type="pred", terms= c("ZTemp", "ZDistance [-1.094, -0.318, 2.004]"),
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

## only thing that works so far is a secondary axis

## back transform zscores
## temp
mean.vari.temp <- mean(allsyncxWithin$annual_avg)
sd.vari.temp <- sd(allsyncxWithin$annual_avg)
var.back.temp <- (allsyncxWithin$ZTemp * sd.vari.temp) + mean.vari.temp

## distance
mean.vari.dist <- mean(allsyncxWithin$distance)
sd.vari.dist <- sd(allsyncxWithin$distance)
var.back.dist <- (allsyncxWithin$ZDistance * sd.vari.dist) + mean.vari.dist

## overlap/diversity
mean.vari.div <- mean(allsyncxWithin$overlap)
sd.vari.div <- sd(allsyncxWithin$overlap)
var.back.div <- (allsyncxWithin$ZOverlap * sd.vari.div) + mean.vari.div

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

## within basin distance
DisTempx <- sjPlot::plot_model(mem_mixed0, type="pred", terms= c("ZTemp", "ZDistance [-1.094, -0.318, 2.004]"),
                               axis.title = c("Environmental synchrony (Z Score)" ,"Thermal Trait Synchrony"), 
                               legend.title = "Distance Percentile", title = "") + 
  scale_color_discrete(labels = c("5th", "50th", "95th")) + 
  # theme(legend.position = "none") +
  scale_x_continuous( sec.axis = sec_axis(~.*sd.vari+mean.vari, name="Environmental synchrony (Original)" )) +
  theme(text = element_text(size=25))

DisTempx

## within basin overlap
DivTemp <- sjPlot::plot_model(mem_mixed1, type="pred", terms= c("ZTemp", " ZOverlap [-1.436, 0.053, 1.499]"),
                              axis.title = c("Environmental synchrony (Z Score)" ,"Thermal Trait Synchrony"), 
                              legend.title = "Overlap Percentile", title = "") + 
  scale_color_discrete(labels = c("5th", "50th", "95th")) +
  scale_x_continuous( sec.axis = sec_axis(~.*sd.vari.temp+mean.vari.temp, name="Environmental synchrony (Original)" )) +
  theme(text = element_text(size=25))
DivTemp

## join together
WithinPlot <- cowplot::plot_grid(DisTempx, DivTemp, labels="auto") ## grid figures
WithinPlot
?plot_grid

## save out
file.name1 <- paste0(out.dir, "08_biotic_interactions_dist_overlap.jpg")
ggsave(WithinPlot, filename=file.name1, dpi=300, height=12, width=18)

## estimates

## overlap
estsow <- sjPlot::plot_model(mem_mixed1, 
                             show.values=TRUE, show.p=TRUE,
                             title="") +
  scale_y_continuous(limits = c(-0.02, 0.02)) +
  theme(text = element_text(size=25))

estsow

## distance

estsdw <- sjPlot::plot_model(mem_mixed0, 
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
file.name1 <- paste0(out.dir, "08_distance_ests_rels.jpg")
ggsave(DistPlot, filename=file.name1, dpi=600, height=15, width=20)

## overlap
DivPlot <- cowplot::plot_grid(estsow, DivTemp, labels="auto") ## grid figures
DivPlot

## save out
file.name1 <- paste0(out.dir, "08_overlap_ests_rels.jpg")
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


