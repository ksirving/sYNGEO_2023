### mixed effects models

# install.packages("haven")
# library(haven)

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
# library(ecodist)
# 

# 
# 
# 
# 
# library(lmerTest) ## from Luke et al 2016 - evaluating significance in linear mixed-effects models in R
# library("easystats") ## multicolinearality
# library(lmerMultiMember)
# # install_github("jvparidon/lmerMultiMember")
# library("scales")
# library(devtools)



# library(devtools)
# install_github("jvparidon/lmerMultiMember")
# install.packages("performance")

# install.packages("remotes")
# remotes::install_github("reumandc/mms")
getwd()
# install.packages("PopGenReport")
# library(PopGenReport)

load(file = "sync/04_temp_pref_env_dist_no_dupl_pairs_WC.RData")
head(allsyncx)

## sites with 1 species
load(file="output_data/01a_sites_3_species.RData")
head(singleSp3l)
length(singleSp3l) ## 151

## change names and filter single species sites

allsyncx <- allsyncx %>%
  rename(overlap = diversity2) %>%
  # group_by(Region) %>%
  # mutate(distance1 = rescale(distance))# %>%
  filter(!Site_ID1 %in% singleSp3l, !Site_ID2 %in% singleSp3l) #%>%
  # mutate(distance2 = (distance - mean(na.omit(distance)))/sd(na.omit(distance))) #%>% ## makes negative numbers
  # pivot_wider(names_from=Metric, values_from = Values)


# range(allsyncx$overlap)
# range(allsyncx$distance) ## might need to scale this to 1
# 
# ### check NAs
# sum(is.na(allsyncx)) ## 170581
# sum(is.na(allsyncx$WCDistkm)) ## 169863 - most are WC distance
# 170581-169863 ## 718
# 
# which(is.na(allsyncx$Sync))
# which(is.na(allsyncx$distance))
# which(is.na(allsyncx$overlap))
# which(is.na(allsyncx$Trait))
# which(is.na(allsyncx$annual_avg))
# sum(is.na(allsyncx$annual_avg)) ## here are the NAs 718
# 
# ### look to see why missing
# ind <- which(is.na(allsyncx$annual_avg)) 
# 
# mis <- allsyncx[ind,] ## all sites paired with same site, can remove

## remove NAs from temp column
allsyncx <- allsyncx %>%
  drop_na(annual_avg)

out.dir <- "/Users/katieirving/OneDrive - SCCWRP/Documents - Katie’s MacBook Pro/git/sYNGEO_2023/Figures/"


# CORRELATIONS ------------------------------------------------------------

## all
allsyncxCor <- allsyncx %>%
  select(distance, overlap, annual_avg, DistKM)

syncCor <- cor(allsyncxCor)
syncCor
write.csv(syncCor, "output_data/06_cor_dist_overlap_temp_n4_sites.csv")

# distance    overlap annual_avg     DistKM
# distance    1.0000000 -0.5534340 -0.6496278  0.6225063
# overlap    -0.5534340  1.0000000  0.3202372 -0.3441639
# annual_avg -0.6496278  0.3202372  1.0000000 -0.9031770
# DistKM      0.6225063 -0.3441639 -0.9031770  1.0000000

## within basin
allsyncxCorWithin <- allsyncx %>%
  filter(Connectivity == 1) %>%
  select(distance, overlap, annual_avg, DistKM, WCDistkm)

syncCorwithin <- cor(allsyncxCorWithin, use = "complete.obs")
syncCorwithin
# ?cor

#               distance     overlap annual_avg     DistKM    WCDistkm
# distance    1.00000000 -0.65278800 -0.1003695  0.1468929  0.09476159
# overlap    -0.65278800  1.00000000 -0.0278078 -0.0122894  0.04630234
# annual_avg -0.10036947 -0.02780780  1.0000000 -0.9058248 -0.90126546
# DistKM      0.14689286 -0.01228940 -0.9058248  1.0000000  0.93261998
# WCDistkm    0.09476159  0.04630234 -0.9012655  0.9326200  1.00000000

write.csv(syncCorwithin, "output_data/06_cor_dist_overlap_temp_within_Basin_n4_sites.csv")

## between basin
allsyncxCorBetween <- allsyncx %>%
  filter(Connectivity == 0) %>%
  select(distance, overlap, annual_avg, DistKM)

syncCorBetween <- cor(allsyncxCorBetween)
syncCorBetween

# distance    overlap annual_avg     DistKM
# distance    1.0000000 -0.5446567 -0.6415716  0.6127346
# overlap    -0.5446567  1.0000000  0.2976611 -0.3179755
# annual_avg -0.6415716  0.2976611  1.0000000 -0.8987739
# DistKM      0.6127346 -0.3179755 -0.8987739  1.0000000

write.csv(syncCorBetween, "output_data/06_cor_dist_overlap_temp_between_Basin_n4_sites.csv")

syncCorwithin
syncCorBetween

# Models: Within Basin spatial decay ----------------------------------------------------

names(allsyncx)
head(allsyncx)
dim(allsyncx)

hist(log(allsyncx$overlap+1))
hist(allsyncx$overlap)
hist(allsyncx$distance)

## connectiovity as a factor and change name 
allsyncx$Connectivity <- as.factor(allsyncx$Connectivity)
allsyncx$Connectivity <- recode_factor(allsyncx$Connectivity,  "1" = "Within Basin", "0" = "Between Basin") 

## distance as log and sqrt

allsyncx <- allsyncx %>%
  mutate(DistKMLog = log(DistKM+1),
         DistKMsqrt = sqrt(DistKM+1))

range(allsyncx$DistKM)

### check sync of 1
withinB <- allsyncx %>%
  filter(Connectivity == "Within Basin")

BetweenB <- allsyncx %>%
  filter(Connectivity == "Between Basin", Sync ==1) 

length(unique(BetweenB$Site_ID2)) ## 3 sites with sync of 1
length(unique(BetweenB$Pair)) ##  4 pairs with sync of 1
length(unique(allsyncx$Site_ID2)) ## 613


round(range(na.omit(withinB$DistKM), digits =2))

round(range(na.omit(BetweenB$DistKM), digits =2))


# ## make longer to get site names for membership model
allsyncxWithin <- allsyncx %>%
  filter(Connectivity == "Within Basin") %>%
  pivot_longer(Site_ID1:Site_ID2, names_to = "SiteNumber", values_to = "SiteName") 

Wa <- lmerMultiMember::weights_from_vector(allsyncxWithin$Region)
Wj <- Matrix::fac2sparse(allsyncxWithin$SiteName)  # convert single membership vars to an indicator matrix with fac2sparse()
Waj <- interaction_weights(Wa, Wj)

## mixed model with distance and connectivity
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8784019/#:~:text=As%20linear%20mixed%2Deffects%20models,associated%20with%20a%20random%20effect.
## paper suggesting that fewer than 5 levels in Random effect is ok when only looking at fixed effects, but not good when analysing random effects

mem_mixed0 <- lmerMultiMember::lmer(log(Sync) ~ DistKM + WCDistkm
                                      + (1 | Region ) + ## add predictors here to get random effect per region
                                      + (1 | RegionXSiteName), 
                                    memberships = list(Region = Wa, RegionXSiteName = Waj), 
                                    REML = T,
                                    data = allsyncxWithin)
## coefs
summary(mem_mixed0, ddf = "Satterthwaite")
anova(mem_mixed0, ddf = "Satterthwaite")
r2_nakagawa(mem_mixed0) ## 0.15
check_singularity(mem_mixed0) ## False

### plots 
class(mem_mixed0) <- "lmerModLmerTest"

# sjPlot::plot_model(mem_mixed) 
ests <- sjPlot::plot_model(mem_mixed0, 
                           show.values=TRUE, show.p=TRUE,
                           title="Drivers of Thermal Synchrony")

ests
file.name1 <- paste0(out.dir, "withinBasin_Spatial_decay_effect_sizes_n4_sites.jpg")
ggsave(ests, filename=file.name1, dpi=300, height=8, width=10)

spatPlotED <- sjPlot::plot_model(mem_mixed0, type="pred", terms=c("DistKM"),
                                 axis.title = c("Euclidean Distance KM", "Thermal Trait Synchrony"), pred.type="re", 
                                 ci.lvl=NA, show.data  = T, dot.size = 0.2)
spatPlotED

spatPlotWC <- sjPlot::plot_model(mem_mixed0, type="pred", terms=c("WCDistkm"),
                                 axis.title = c("Water Course Distance KM", "Thermal Trait Synchrony"), pred.type="re", 
                                 ci.lvl=NA, show.data  = T, dot.size = 0.2)
spatPlotWC

## non significant interactions
spatDec <- plot_grid(list(spatPlotED, spatPlotWC))## grid figures

#terms=c("Euclid_Dist_Meters", "Connectivity")
file.name1 <- paste0(out.dir, "withinBasin_Spatial_decay_plots_n4_sites.jpg")
ggsave(spatDec, filename=file.name1, dpi=300, height=8, width=10)




# lattice::qqmath(mem_mixed1)
# plot(mem_mixed1, type=c("p","smooth"), col.line=1)
# check_model(mem_mixed1)





# Membership Model: within Basin Region as fixed ---------------------------------------

# Wa <- lmerMultiMember::weights_from_vector(allsyncxWithin$Region)
# Wj <- Matrix::fac2sparse(allsyncxWithin$SiteName)  # convert single membership vars to an indicator matrix with fac2sparse()
# Waj <- interaction_weights(Wa, Wj)
round(range(allsyncxWithin$overlap),digits =2)
## model with new diversity measure

mem_mixed0 <- lmerMultiMember::lmer(Sync ~  (overlap*annual_avg + distance*annual_avg)
                                    + (1 | Region ) + ## add predictors here to get random effect per region
                                      (1 | RegionXSiteName), 
                                    memberships = list(Region = Wa, RegionXSiteName = Waj), 
                                    REML = T,
                                    data = allsyncxWithin)


summary(mem_mixed0, ddf = "Satterthwaite")
anova(mem_mixed0, ddf = "Satterthwaite")
r2_nakagawa(mem_mixed0) ## 0.16
check_singularity(mem_mixed0) ## False

# lattice::qqmath(mem_mixed1)
# plot(mem_mixed1, type=c("p","smooth"), col.line=1)
# check_model(mem_mixed1)

### plots 
class(mem_mixed0) <- "lmerModLmerTest"
# sjPlot::plot_model(mem_mixed) 
ests <- sjPlot::plot_model(mem_mixed0, 
                           show.values=TRUE, show.p=TRUE,
                           title="Drivers of Thermal Synchrony")

ests
file.name1 <- paste0(out.dir, "withinBasin_effect_sizes_all_expsync_n4_sites.jpg")
ggsave(ests, filename=file.name1, dpi=300, height=8, width=10)

estsTab <- sjPlot::tab_model(mem_mixed0, 
                             show.re.var= TRUE, 
                             pred.labels =c("(Intercept)", "DistKM", "annual avg", "diversity", "Connectivity"),
                             dv.labels= "Drivers of Thermal Synchrony")

estsTab

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
quantile(allsyncxWithin$annual_avg, probs = c(0.05,0.5,0.95))
round(quantile(allsyncxWithin$overlap, probs = c(0.05,0.5,0.95)), digits = 3)
round(quantile(allsyncxWithin$distance, probs = c(0.05,0.5,0.95)), digits = 3)

## overlap v temp
tempDiv <- sjPlot::plot_model(mem_mixed0, type="pred", terms= c("overlap", "annual_avg [0.96, 0.98, 0.999]"),
                              axis.title = c( "Trait Overlap","Thermal Trait Synchrony"), 
                              legend.title = "Environmental synchrony")
tempDiv


DivTemp <- sjPlot::plot_model(mem_mixed0, type="pred", terms= c("annual_avg", " overlap [0.0, 0.39, 0.85]"),
                              axis.title = c("Environmental synchrony" ,"Thermal Trait Synchrony"), 
                              legend.title = "Trait Overlap")
DivTemp

## significant interactions
CompPlotSigDiv <- plot_grid(list(tempDiv, DivTemp)) ## grid figures

file.name1 <- paste0(out.dir, "withinBasin_overlap_sig_interactions_n4_sites.jpg")
ggsave(CompPlotSigDiv, filename=file.name1, dpi=300, height=12, width=18)

## distance v temp
tempDist <- sjPlot::plot_model(mem_mixed0, type="pred", terms= c("distance", "annual_avg [0.96, 0.98, 0.999]"),
                               axis.title = c("Trait Distance", "Thermal Trait Synchrony"), 
                               legend.title = "Environmental synchrony")
tempDist

DistTemp <- sjPlot::plot_model(mem_mixed0, type="pred", terms= c("annual_avg", "distance [0.001, 0.018, 0.082]"),
                              axis.title = c("Environmental synchrony" ,"Thermal Trait Synchrony"), 
                              legend.title = "Trait Distance")
DistTemp

## significant interactions
CompPlotSigDist <- plot_grid(list( tempDist, DistTemp)) ## grid figures

file.name1 <- paste0(out.dir, "withinBasin_distance_sig_interactions_n4_sites.jpg")
ggsave(CompPlotSigDist, filename=file.name1, dpi=300, height=12, width=18)
# CompPlotInSig
# file.name1 <- paste0(out.dir, "comp_mod_insig_interactions_region_removed_v2.jpg")
# ggsave(CompPlotInSig, filename=file.name1, dpi=300, height=10, width=15)


# Model selection: within Basin -------------------------------------------

### temp
mem_mixed1 <- lmerMultiMember::lmer(Sync ~  annual_avg  #Connectivity
                                    + (1 | Region ) + ## add predictors here to get random effect per region
                                      (1 | RegionXSiteName), 
                                    memberships = list(Region = Wa, RegionXSiteName = Waj), 
                                    REML = T,
                                    data = allsyncxWithin)
mem_mixed1
## add in interactions between other variables
#(annual_avg*log(diversity2) + distance*log(diversity2)) * Connectivity
# (annual_avg +distance+diversity2))*Connectivity
summary(mem_mixed1, ddf = "Satterthwaite")
anova(mem_mixed1, ddf = "Satterthwaite")
r2_nakagawa(mem_mixed1) ## 0.19

### plots 
class(mem_mixed1) <- "lmerModLmerTest"
# sjPlot::plot_model(mem_mixed) 
ests <- sjPlot::plot_model(mem_mixed1, 
                           show.values=TRUE, show.p=TRUE,
                           title="Drivers of Thermal Synchrony")

ests
file.name1 <- paste0(out.dir, "withinBasin_effect_sizes_temp_n4_sites.jpg")
ggsave(ests, filename=file.name1, dpi=300, height=8, width=10)


tempMod <- sjPlot::plot_model(mem_mixed1, type="pred", terms=c("annual_avg"),
                                 axis.title = c("Environmental Synchrony", "Thermal Trait Synchrony"), pred.type="re", 
                                 ci.lvl=NA, show.data  = T, dot.size = 0.2)
tempMod

file.name1 <- paste0(out.dir, "withinBasin_relationship_temp.jpg")
ggsave(tempMod, filename=file.name1, dpi=300, height=8, width=10)

### distance
mem_mixed2 <- lmerMultiMember::lmer(Sync ~  distance  #Connectivity
                                    + (1 | Region ) + ## add predictors here to get random effect per region
                                      (1 | RegionXSiteName), 
                                    memberships = list(Region = Wa, RegionXSiteName = Waj), 
                                    REML = T,
                                    data = allsyncxWithin)
mem_mixed2
## add in interactions between other variables
#(annual_avg*log(diversity2) + distance*log(diversity2)) * Connectivity
# (annual_avg +distance+diversity2))*Connectivity
summary(mem_mixed2, ddf = "Satterthwaite")
anova(mem_mixed2, ddf = "Satterthwaite")
r2_nakagawa(mem_mixed2) ## 0.2

### plots 
class(mem_mixed2) <- "lmerModLmerTest"
# sjPlot::plot_model(mem_mixed) 
ests <- sjPlot::plot_model(mem_mixed2, 
                           show.values=TRUE, show.p=TRUE,
                           title="Drivers of Thermal Synchrony")

ests
file.name1 <- paste0(out.dir, "withinBasin_effect_sizes_distance_n4_sites.jpg")
ggsave(ests, filename=file.name1, dpi=300, height=8, width=10)

distMod <- sjPlot::plot_model(mem_mixed2, type="pred", terms=c("distance"),
                              axis.title = c("Trait Distance", "Thermal Trait Synchrony"), pred.type="re", 
                              ci.lvl=NA, show.data  = T, dot.size = 0.2)
distMod

file.name1 <- paste0(out.dir, "withinBasin_relationship_distance_n4_sites.jpg")
ggsave(distMod, filename=file.name1, dpi=300, height=8, width=10)


### overlap
mem_mixed3 <- lmerMultiMember::lmer(Sync ~  overlap  #Connectivity
                                    + (1 | Region ) + ## add predictors here to get random effect per region
                                      (1 | RegionXSiteName), 
                                    memberships = list(Region = Wa, RegionXSiteName = Waj), 
                                    REML = T,
                                    data = allsyncxWithin)
mem_mixed3
## add in interactions between other variables
#(annual_avg*log(diversity2) + distance*log(diversity2)) * Connectivity
# (annual_avg +distance+diversity2))*Connectivity
summary(mem_mixed3, ddf = "Satterthwaite")
anova(mem_mixed3, ddf = "Satterthwaite")
r2_nakagawa(mem_mixed3) ## 0.2

### plots 
class(mem_mixed3) <- "lmerModLmerTest"
# sjPlot::plot_model(mem_mixed) 
ests <- sjPlot::plot_model(mem_mixed3, 
                           show.values=TRUE, show.p=TRUE,
                           title="Drivers of Thermal Synchrony")

ests
file.name1 <- paste0(out.dir, "withinBasin_effect_sizes_overlap_n4_sites.jpg")
ggsave(ests, filename=file.name1, dpi=300, height=8, width=10)

ovMod <- sjPlot::plot_model(mem_mixed3, type="pred", terms=c("overlap"),
                              axis.title = c("Trait Overlap", "Thermal Trait Synchrony"), pred.type="re", 
                              ci.lvl=NA, show.data  = T, dot.size = 0.2)
ovMod

file.name1 <- paste0(out.dir, "withinBasin_relationship_overlap_n4_sites.jpg")
ggsave(ovMod, filename=file.name1, dpi=300, height=8, width=10)


### temp and distance

mem_mixed4 <- lmerMultiMember::lmer(Sync ~  annual_avg * distance  #Connectivity
                                    + (1 | Region ) + ## add predictors here to get random effect per region
                                      (1 | RegionXSiteName), 
                                    memberships = list(Region = Wa, RegionXSiteName = Waj), 
                                    REML = T,
                                    data = allsyncxWithin)
mem_mixed4
## add in interactions between other variables
#(annual_avg*log(diversity2) + distance*log(diversity2)) * Connectivity
# (annual_avg +distance+diversity2))*Connectivity
summary(mem_mixed4, ddf = "Satterthwaite")
anova(mem_mixed4, ddf = "Satterthwaite")
r2_nakagawa(mem_mixed4) ## 0.19

### plots 
class(mem_mixed4) <- "lmerModLmerTest"
# sjPlot::plot_model(mem_mixed) 
ests <- sjPlot::plot_model(mem_mixed4, 
                           show.values=TRUE, show.p=TRUE,
                           title="Drivers of Thermal Synchrony")

ests

file.name1 <- paste0(out.dir, "withinBasin_effect_sizes_temp_distance_n4_sites.jpg")
ggsave(ests, filename=file.name1, dpi=300, height=8, width=10)

Mod4 <- sjPlot::plot_model(mem_mixed4, type="pred", terms=c("distance", "annual_avg [0.96, 0.98, 0.999]"),
                            axis.title = c("Trait Distance", "Thermal Trait Synchrony"), pred.type="re", 
                            ci.lvl=NA, show.data  = T, dot.size = 0.2)
Mod4

file.name1 <- paste0(out.dir, "withinBasin_relationship_temp_distance_n4_sites.jpg")
ggsave(Mod4, filename=file.name1, dpi=300, height=8, width=10)

### temp and overlap

mem_mixed5 <- lmerMultiMember::lmer(Sync ~  annual_avg * overlap  #Connectivity
                                    + (1 | Region ) + ## add predictors here to get random effect per region
                                      (1 | RegionXSiteName), 
                                    memberships = list(Region = Wa, RegionXSiteName = Waj), 
                                    REML = T,
                                    data = allsyncxWithin)
mem_mixed5
## add in interactions between other variables
#(annual_avg*log(diversity2) + distance*log(diversity2)) * Connectivity
# (annual_avg +distance+diversity2))*Connectivity
summary(mem_mixed5, ddf = "Satterthwaite")
anova(mem_mixed5, ddf = "Satterthwaite")
r2_nakagawa(mem_mixed5) ## 0.19

### plots 
class(mem_mixed5) <- "lmerModLmerTest"
# sjPlot::plot_model(mem_mixed) 
ests <- sjPlot::plot_model(mem_mixed5, 
                           show.values=TRUE, show.p=TRUE,
                           title="Drivers of Thermal Synchrony")

ests

file.name1 <- paste0(out.dir, "withinBasin_effect_sizes_temp_distance_n4_sites.jpg")
ggsave(ests, filename=file.name1, dpi=300, height=8, width=10)

Mod5 <- sjPlot::plot_model(mem_mixed5, type="pred", terms=c("overlap", "annual_avg [0.96, 0.98, 0.999]"),
                           axis.title = c("Trait Overlap", "Thermal Trait Synchrony"), pred.type="re", 
                           ci.lvl=NA, show.data  = T, dot.size = 0.2)
Mod5

file.name1 <- paste0(out.dir, "withinBasin_relationship_temp_overlap_n4_sites.jpg")
ggsave(Mod5, filename=file.name1, dpi=300, height=8, width=10)


######

## plot
set_theme(base = theme_classic(), #To remove the background color and the grids
          theme.font = 'serif',   #To change the font type
          axis.title.size = 1.3,  #To change axis title size
          axis.textsize.x = 1.2,    #To change x axis text size
          # axis.angle.x = 60,      #To change x axis text angle
          # axis.hjust.x = 1,
          # axis.ticksmar = 50,
          axis.textsize.y = 1)  #To change y axis text size

theme_set(theme_sjplot())

quantile(allsyncx$annual_avg, probs = c(0.05,0.5,0.95))
round(quantile(allsyncx$overlap, probs = c(0.05,0.5,0.95)), digits = 3)
round(quantile(allsyncx$distance, probs = c(0.05,0.5,0.95)), digits = 3)


tempDiv <- sjPlot::plot_model(mem_mixed0, type="pred", terms= c("overlap",  "annual_avg [0.68, 0.92, 0.998]"),
                              axis.title = c("Trait Overlap", "Thermal Trait Synchrony"), 
                              legend.title = "Environmental synchrony")


tempDist <- sjPlot::plot_model(mem_mixed0, type="pred", terms= c("distance",  "annual_avg [0.68, 0.92, 0.998]"),
                              axis.title = c("Trait Distance", "Thermal Trait Synchrony"), 
                              legend.title = "Environmental synchrony")

Divtemp <- sjPlot::plot_model(mem_mixed0, type="pred", terms= c("overlap",  "annual_avg [0.68, 0.92, 0.998]"),
                              axis.title = c("Environmental Synchrony", "Thermal Trait Synchrony"), 
                              legend.title = "Trait Overlap")


Disttemp <- sjPlot::plot_model(mem_mixed0, type="pred", terms= c("annual_avg",  "distance [0.84, 0.98, 1.13]"),
                               axis.title = c("Environmental Synchrony", "Thermal Trait Synchrony"), 
                               legend.title = "Trait Distance")

tempDist
## significant interactions
CompPlotSig <- plot_grid(list(tempDiv, tempDist, Divtemp, Disttemp)) ## grid figures

##insignificant interactions
CompPlotInSig <- plot_grid(list(tempDiv, tempDist, distCon, tempDistCon)) ## grid figures

#terms=c("Euclid_Dist_Meters", "Connectivity")
file.name1 <- paste0(out.dir, "comp_mod_sig_interactions_region_connectivity_removed_v2.jpg")
ggsave(CompPlotSig, filename=file.name1, dpi=300, height=12, width=18)

file.name1 <- paste0(out.dir, "comp_mod_insig_interactions_region_removed_v2.jpg")
ggsave(CompPlotInSig, filename=file.name1, dpi=300, height=10, width=15)

?plot_model
tempDist <- compPlot1[[2]]
tempCon <- compPlot1[[3]]
tempReg <- compPlot1[[4]]
DivRegion <- compPlot1[[5]]
DistRegion<- compPlot1[[6]]
ConReg <- compPlot1[[7]]
tempDivReg <- compPlot1[[8]]
tempDistReg <- compPlot1[[9]]
tempConReg <- compPlot1[[10]]



## non significant interactions
CompPlotnon <- plot_grid(list(tempCon, tempReg, ConReg, tempConReg))## grid figures

#terms=c("Euclid_Dist_Meters", "Connectivity")
file.name1 <- paste0(out.dir, "comp_mod_interactions_non_significant_region_removed_v2.jpg")
ggsave(CompPlotnon, filename=file.name1, dpi=300, height=8, width=10)



# compPlot2 <- sjPlot::plot_model(mem_mixed0, type="int", terms=c("annual_avg","distance", "Connectivity","Region"),
#                                 pred.type="re", 
#                                 ci.lvl=NA)
# #axis.title = c("Temp sync", "Thermal Trait Synchrony"),
# compPlot2
# #terms=c("Euclid_Dist_Meters", "Connectivity")
# file.name1 <- paste0(out.dir, "distance_composition_model_by_region.jpg")
# ggsave(compPlot2, filename=file.name1, dpi=300, height=8, width=10)


# conPlot3 <- sjPlot::plot_model(mem_mixed0, type="pred", terms=c("annual_avg","overlap", "distance","Connectivity"),
#                                 pred.type="re", 
#                                 ci.lvl=NA)
# 
# conPlot3
# #terms=c("Euclid_Dist_Meters", "Connectivity")
# file.name1 <- paste0(out.dir, "composition_model_by_connectivity.jpg")
# ggsave(conPlot3, filename=file.name1, dpi=300, height=8, width=10)
# 
# 
# RegPlot3 <- sjPlot::plot_model(mem_mixed0, type="pred", terms=c("annual_avg","overlap", "Region","Connectivity"),
#                                 pred.type="re", 
#                                ci.lvl=NA)
# 
# RegPlot3
# #terms=c("Euclid_Dist_Meters", "Connectivity")
# file.name1 <- paste0(out.dir, "composition_model_by_connectivityRegion.jpg")
# ggsave(RegPlot3, filename=file.name1, dpi=300, height=8, width=10)



