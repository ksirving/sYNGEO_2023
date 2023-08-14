## sensitivty analysis of single species sites

## work plan
## remove single species sites
## compare with full data mixed model results
## remove 2 & 3 species sites
## compare with 1 species and full data results


## packages

library(tidyverse)
library(tidylog)
library(ecodist)

library(performance)
library(lme4)

library(cowplot) #for manuscript ready figures
library(sjPlot) #for plotting lmer and glmer mods
library(sjmisc) 
library(sjlabelled)

library(effects)
library(sjstats) #use for r2 functions
library(lmerTest) ## from Luke et al 2016 - evaluating significance in linear mixed-effects models in R
library("easystats") ## multicolinearality
library(lmerMultiMember)
library("scales")


# format sites for Stefano ------------------------------------------------

## read in Stefanos sites
stefSites <- read.csv("input_data/Env/TableSites.csv")
head(stefSites)


## upload fish abundance and site data
originaldata <- read.csv("input_data/Bio/fishdata_selection_basins_same_time_window_10262020.csv")
head(originaldata)


# Load and format data ----------------------------------------------------


## load data
load(file = "output_data/sync/04_temp_pref_env_dist_no_dupl_pairs.RData")
head(allsyncx)

allsyncx <- allsyncx %>%
  rename(overlap = diversity2)

## some NAs in sync. find and remove - from sites with 1 species, cannot compute synchrony
ind <- which(is.na(allsyncx$Sync))
allsyncx[ind,]
allsyncx <- allsyncx[-ind,]
sum(is.na(allsyncx))

allsyncx <- na.omit(allsyncx)

dim(allsyncx) ## 205592

out.dir <- "/Users/katieirving/OneDrive - SCCWRP/Documents - Katie’s MacBook Pro/git/sYNGEO_Func_Sync_V3/Figures/"


# Remove single species sites ---------------------------------------------

## upload list 

load(file="output_data/01a_sites_1_species.RData")
head(singleSpl)
length(singleSpl) ## 49

allsyncx <- allsyncx %>%
  filter(!Site_ID1 %in% singleSpl, !Site_ID2 %in% singleSpl) ## remove sites from both sites in pair

dim(allsyncx) ## 176731


# Correlation -------------------------------------------------------------

## all
allsyncxCor <- allsyncx %>%
  select(distance, overlap, annual_avg, DistKM)

syncCor <- cor(allsyncxCor)

write.csv(syncCor, "output_data/06_cor_dist_overlap_temp_single_sp.csv")

## within basin
allsyncxCorWithin <- allsyncx %>%
  filter(Connectivity == 1) %>%
  select(distance, overlap, annual_avg, DistKM)

syncCorwithin <- cor(allsyncxCorWithin)

write.csv(syncCorwithin, "output_data/06_cor_dist_overlap_temp_within_Basin_single_sp.csv")

## between basin
allsyncxCorBetween <- allsyncx %>%
  filter(Connectivity == 0) %>%
  select(distance, overlap, annual_avg, DistKM)

syncCorBetween <- cor(allsyncxCorBetween)

write.csv(syncCorBetween, "output_data/06_cor_dist_overlap_temp_between_Basin_single_sp.csv")

syncCor
syncCorwithin
syncCorBetween


# Modes: spatial decay ----------------------------------------------------

head(allsyncx)
## connectivity as a factor and change name 
allsyncx$Connectivity <- as.factor(allsyncx$Connectivity)
allsyncx$Connectivity <- recode_factor(allsyncx$Connectivity,  "1" = "Within Basin", "0" = "Between Basin") 


## distance as log and sqrt

allsyncx <- allsyncx %>%
  mutate(DistKMLog = log(DistKM+1),
         DistKMsqrt = sqrt(DistKM+1))

range(allsyncx$DistKM)

# ## make longer to get site names for membership model
allsyncx <- allsyncx %>%
  pivot_longer(Site_ID1:Site_ID2, names_to = "SiteNumber", values_to = "SiteName") 

Wa <- lmerMultiMember::weights_from_vector(allsyncx$Region)
Wj <- Matrix::fac2sparse(allsyncx$SiteName)  # convert single membership vars to an indicator matrix with fac2sparse()
Waj <- interaction_weights(Wa, Wj)

## recode region
allsyncx$Region <- as.factor(allsyncx$Region)
allsyncx$Region <- recode_factor(allsyncx$Region,  "USA" = "USA", "Oceania" = "Oceania", "Europe" = "Europe") 

## mixed model with distance and connectivity
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8784019/#:~:text=As%20linear%20mixed%2Deffects%20models,associated%20with%20a%20random%20effect.
## paper suggesting that fewer than 5 levels in Random effect is ok when only looking at fixed effects, but not good when analysing random effects

mem_mixed0 <- lmerMultiMember::lmer(log(Sync) ~ DistKM  + Connectivity
                                    + (1 | Region ) + ## add predictors here to get random effect per region
                                      + (1 | RegionXSiteName), 
                                    memberships = list(Region = Wa, RegionXSiteName = Waj), 
                                    REML = T,
                                    data = allsyncx)

spatPlot2a <- sjPlot::plot_model(mem_mixed0, type="pred", terms=c("DistKM", "Connectivity"),
                                 axis.title = c("Euclidean Distance KM", "Thermal Trait Synchrony"), pred.type="re", 
                                 ci.lvl=NA, show.data  = T, dot.size = 0.2)

spatPlot2a 
file.name1 <- paste0(out.dir, "spatial_decay_mem_mod_exp_single_sp.jpg")
ggsave(spatPlot2a, filename=file.name1, dpi=300, height=8, width=10)


summary(mem_mixed0, ddf = "Satterthwaite")
anova(mem_mixed0, ddf = "Satterthwaite")
r2_nakagawa(mem_mixed0) 
check_singularity(mem_mixed0) ## false


### plots 
class(mem_mixed0) <- "lmerModLmerTest"
# sjPlot::plot_model(mem_mixed) 
ests <- sjPlot::plot_model(mem_mixed0, 
                           show.values=TRUE, show.p=TRUE,
                           title="Drivers of Thermal Synchrony")

ests
file.name1 <- paste0(out.dir, "effect_sizes_spatial_decay_single_sp.jpg")
ggsave(ests, filename=file.name1, dpi=300, height=8, width=10)


# Membership Model: Drivers ---------------------------------------


## model with new diversity measure

mem_mixed0 <- lmerMultiMember::lmer(Sync ~  (overlap*annual_avg + distance*annual_avg)*Connectivity
                                    + (1 | Region ) + ## add predictors here to get random effect per region
                                      (1 | RegionXSiteName), 
                                    memberships = list(Region = Wa, RegionXSiteName = Waj), 
                                    REML = T,
                                    data = allsyncx)

summary(mem_mixed0, ddf = "Satterthwaite")
anova(mem_mixed0, ddf = "Satterthwaite")
r2_nakagawa(mem_mixed0) 
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
file.name1 <- paste0(out.dir, "effect_sizes_diversity2_overlap_region_removed_single_sp.jpg")
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

quantile(allsyncx$annual_avg, probs = c(0.05,0.5,0.95))
round(quantile(allsyncx$overlap, probs = c(0.05,0.5,0.95)), digits = 3)
round(quantile(allsyncx$distance, probs = c(0.05,0.5,0.95)), digits = 3)


tempDiv <- sjPlot::plot_model(mem_mixed0, type="pred", terms= c("overlap",  "annual_avg [0.68, 0.92, 0.998]"),
                              axis.title = c("Trait Overlap", "Thermal Trait Synchrony"), 
                              legend.title = "Environmental synchrony", pred.type="re", 
                              ci.lvl=NA, show.data  = T, dot.size = 0.2)

tempDiv
tempDist <- sjPlot::plot_model(mem_mixed0, type="pred", terms= c("distance",  "annual_avg [0.68, 0.92, 0.998]"),
                               axis.title = c("Trait Distance", "Thermal Trait Synchrony"), 
                               legend.title = "Environmental synchrony", pred.type="re", 
                               ci.lvl=NA, show.data  = T, dot.size = 0.2)
tempDist
Disttemp <- sjPlot::plot_model(mem_mixed0, type="pred", terms= c("annual_avg",  "distance [0.84, 0.98, 1.13]"),
                               axis.title = c("Environmental Synchrony", "Thermal Trait Synchrony"), 
                               legend.title = "Trait Distance", pred.type="re", 
                               ci.lvl=NA, show.data  = T, dot.size = 0.2)


divCon <- sjPlot::plot_model(mem_mixed0, type="pred", terms= c("overlap", "Connectivity"),
                             axis.title = c("Trait Overlap", "Thermal Trait Synchrony"), 
                             legend.title = "Connectivity", pred.type="re", 
                             ci.lvl=NA, show.data  = T, dot.size = 0.2)
divCon

distCon <- sjPlot::plot_model(mem_mixed0, type="pred", terms= c("distance", "Connectivity"),
                              axis.title = c("Trait Distance", "Thermal Trait Synchrony"), 
                              legend.title = "Connectivity", pred.type="re", 
                              ci.lvl=NA, show.data  = T, dot.size = 0.2)
distCon

tempDistCon <- sjPlot::plot_model(mem_mixed0, type="pred", terms= c("distance","Connectivity", "annual_avg [0.68, 0.92, 0.998]" ),
                                  axis.title = c("Trait Distance", "Thermal Trait Synchrony"), 
                                  legend.title = "Connectivity", pred.type="re")
tempDistCon

tempDistCon2 <- sjPlot::plot_model(mem_mixed0, type="pred", terms= c("annual_avg","Connectivity", "distance [0.02, 0.056, 0.187]" ),
                                   axis.title = c("Environmental synchrony", "Thermal Trait Synchrony"), 
                                   legend.title = "Connectivity", pred.type="re")

tempDivCon <- sjPlot::plot_model(mem_mixed0, type="pred", terms= c("overlap","Connectivity", "annual_avg [0.68, 0.92, 0.998]" ),
                                 axis.title = c("Trait Overlap", "Thermal Trait Synchrony"), 
                                 legend.title = "Connectivity")
tempDivCon
tempDivCon2 <- sjPlot::plot_model(mem_mixed0, type="pred", terms= c("annual_avg","Connectivity", "overlap [0, 0.035, 0.795]" ),
                                  axis.title = c("Environmental synchrony", "Thermal Trait Synchrony"), 
                                  legend.title = "Connectivity")


tempDivCon2
?plot_model
## significant interactions
CompPlotSig <- plot_grid(list(tempDist, distCon, tempDistCon, tempDistCon2))

##insignificant interactions
# CompPlotInSig <- plot_grid(list(tempDiv, tempDist, distCon, tempDistCon)) ## grid figures

#terms=c("Euclid_Dist_Meters", "Connectivity")
file.name1 <- paste0(out.dir, "comp_mod_sig_interactions_region_removed_single_sp.jpg")
ggsave(CompPlotSig, filename=file.name1, dpi=300, height=12, width=18)

# file.name1 <- paste0(out.dir, "comp_mod_insig_interactions_region_removed_single_sp.jpg")
# ggsave(CompPlotInSig, filename=file.name1, dpi=300, height=10, width=15)
# 
# ### model selection 
# mem_mixed1 <- lmerMultiMember::lmer(Sync ~  annual_avg*(overlap+distance)  #Connectivity
#                                     + (1 | Region ) + ## add predictors here to get random effect per region
#                                       (1 | RegionXSiteName), 
#                                     memberships = list(Region = Wa, RegionXSiteName = Waj), 
#                                     REML = T,
#                                     data = allsyncx)

## add in interactions between other variables
#(annual_avg*log(diversity2) + distance*log(diversity2)) * Connectivity
# (annual_avg +distance+diversity2))*Connectivity
summary(mem_mixed1, ddf = "Satterthwaite")
anova(mem_mixed1, ddf = "Satterthwaite")
r2_nakagawa(mem_mixed1) 

### plots 
class(mem_mixed1) <- "lmerModLmerTest"
# sjPlot::plot_model(mem_mixed) 
ests <- sjPlot::plot_model(mem_mixed1, 
                           show.values=TRUE, show.p=TRUE,
                           title="Drivers of Thermal Synchrony")

ests

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
                              legend.title = "Environmental synchrony", pred.type="re", 
                              ci.lvl=NA, show.data  = T, dot.size = 0.2)


tempDist <- sjPlot::plot_model(mem_mixed0, type="pred", terms= c("distance",  "annual_avg [0.68, 0.92, 0.998]"),
                               axis.title = c("Trait Distance", "Thermal Trait Synchrony"), 
                               legend.title = "Environmental synchrony", pred.type="re", 
                               ci.lvl=NA, show.data  = T, dot.size = 0.2)

Disttemp <- sjPlot::plot_model(mem_mixed0, type="pred", terms= c("annual_avg",  "distance [0.84, 0.98, 1.13]"),
                               axis.title = c("Environmental Synchrony", "Thermal Trait Synchrony"), 
                               legend.title = "Trait Distance", pred.type="re", 
                               ci.lvl=NA, show.data  = T, dot.size = 0.2)

Divtemp <- sjPlot::plot_model(mem_mixed0, type="pred", terms= c("overlap",  "annual_avg [0.68, 0.92, 0.998]"),
                              axis.title = c("Environmental Synchrony", "Thermal Trait Synchrony"), 
                              legend.title = "Trait Overlap", pred.type="re", 
                              ci.lvl=NA, show.data  = T, dot.size = 0.2)

distCon <- sjPlot::plot_model(mem_mixed0, type="pred", terms= c("distance", "Connectivity"),
                              axis.title = c("Trait Distance", "Thermal Trait Synchrony"), 
                              legend.title = "Connectivity", pred.type="re", 
                              ci.lvl=NA, show.data  = T, dot.size = 0.2)




tempDistCon <- sjPlot::plot_model(mem_mixed0, type="pred", terms= c("distance","Connectivity", "annual_avg [0.68, 0.92, 0.998]" ),
                                  axis.title = c("Trait Distance", "Thermal Trait Synchrony"), 
                                  legend.title = "Connectivity", pred.type="re", 
                                  ci.lvl=NA, show.data  = T, dot.size = 0.2)
tempDistCon

tempDistCon2 <- sjPlot::plot_model(mem_mixed0, type="pred", terms= c("annual_avg","Connectivity", "distance [0.02, 0.056, 0.187]" ),
                                   axis.title = c("Environmental synchrony", "Thermal Trait Synchrony"), 
                                   legend.title = "Connectivity", pred.type="re", 
                                   ci.lvl=NA, show.data  = T, dot.size = 0.2)

tempDist
## significant interactions
CompPlotSig <- plot_grid(list(tempDist, distCon, tempDistCon, tempDistCon2)) ## grid figures

##insignificant interactions
CompPlotInSig <- plot_grid(list(tempDiv, tempDist, tempDistCon)) ## grid figures

#terms=c("Euclid_Dist_Meters", "Connectivity")
file.name1 <- paste0(out.dir, "comp_mod_sig_interactions_region_connectivity_removed_single_sp.jpg")
ggsave(CompPlotSig, filename=file.name1, dpi=300, height=12, width=18)

# file.name1 <- paste0(out.dir, "comp_mod_insig_interactions_region_removed_single_sp.jpg")
# ggsave(CompPlotInSig, filename=file.name1, dpi=300, height=10, width=15)
# 
# ?plot_model
# tempDist <- compPlot1[[2]]
# tempCon <- compPlot1[[3]]
# tempReg <- compPlot1[[4]]
# DivRegion <- compPlot1[[5]]
# DistRegion<- compPlot1[[6]]
# ConReg <- compPlot1[[7]]
# tempDivReg <- compPlot1[[8]]
# tempDistReg <- compPlot1[[9]]
# tempConReg <- compPlot1[[10]]
# 
# 
# 
# ## non significant interactions
# CompPlotnon <- plot_grid(list(tempCon, tempReg, ConReg, tempConReg))## grid figures
# 
# #terms=c("Euclid_Dist_Meters", "Connectivity")
# file.name1 <- paste0(out.dir, "comp_mod_interactions_non_significant_region_removed.jpg")
# ggsave(CompPlotnon, filename=file.name1, dpi=300, height=8, width=10)

# Remove sites with 2 species ---------------------------------------------

load(file = "output_data/sync/04_temp_pref_env_dist_no_dupl_pairs.RData")
head(allsyncx)

allsyncx <- allsyncx %>%
  rename(overlap = diversity2)

## some NAs in sync. find and remove - from sites with 1 species, cannot compute synchrony
ind <- which(is.na(allsyncx$Sync))
allsyncx[ind,]
allsyncx <- allsyncx[-ind,]
sum(is.na(allsyncx))

allsyncx <- na.omit(allsyncx)
## upload list 

load(file="output_data/01a_sites_2_species.RData")
head(singleSp2l)
length(singleSp2l) ## 103

allsyncx <- allsyncx %>%
  filter(!Site_ID1 %in% singleSp2l, !Site_ID2 %in% singleSp2l) ## remove sites from both sites in pair

dim(allsyncx) ## 146410


# Correlation -------------------------------------------------------------

## all
allsyncxCor <- allsyncx %>%
  select(distance, overlap, annual_avg, DistKM)

syncCor <- cor(allsyncxCor)

write.csv(syncCor, "output_data/06_cor_dist_overlap_temp_2_sp.csv")

## within basin
allsyncxCorWithin <- allsyncx %>%
  filter(Connectivity == 1) %>%
  select(distance, overlap, annual_avg, DistKM)

syncCorwithin <- cor(allsyncxCorWithin)

write.csv(syncCorwithin, "output_data/06_cor_dist_overlap_temp_within_Basin_2_sp.csv")

## between basin
allsyncxCorBetween <- allsyncx %>%
  filter(Connectivity == 0) %>%
  select(distance, overlap, annual_avg, DistKM)

syncCorBetween <- cor(allsyncxCorBetween)

write.csv(syncCorBetween, "output_data/06_cor_dist_overlap_temp_between_Basin_2_sp.csv")

syncCor
syncCorwithin
syncCorBetween


# Modes: spatial decay ----------------------------------------------------

head(allsyncx)
## connectivity as a factor and change name 
allsyncx$Connectivity <- as.factor(allsyncx$Connectivity)
allsyncx$Connectivity <- recode_factor(allsyncx$Connectivity,  "1" = "Within Basin", "0" = "Between Basin") 


## distance as log and sqrt

allsyncx <- allsyncx %>%
  mutate(DistKMLog = log(DistKM+1),
         DistKMsqrt = sqrt(DistKM+1))

range(allsyncx$DistKM)

# ## make longer to get site names for membership model
allsyncx <- allsyncx %>%
  pivot_longer(Site_ID1:Site_ID2, names_to = "SiteNumber", values_to = "SiteName") 

Wa <- lmerMultiMember::weights_from_vector(allsyncx$Region)
Wj <- Matrix::fac2sparse(allsyncx$SiteName)  # convert single membership vars to an indicator matrix with fac2sparse()
Waj <- interaction_weights(Wa, Wj)

## recode region
allsyncx$Region <- as.factor(allsyncx$Region)
allsyncx$Region <- recode_factor(allsyncx$Region,  "USA" = "USA", "Oceania" = "Oceania", "Europe" = "Europe") 

## mixed model with distance and connectivity
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8784019/#:~:text=As%20linear%20mixed%2Deffects%20models,associated%20with%20a%20random%20effect.
## paper suggesting that fewer than 5 levels in Random effect is ok when only looking at fixed effects, but not good when analysing random effects

mem_mixed0 <- lmerMultiMember::lmer(log(Sync) ~ DistKM  + Connectivity
                                    + (1 | Region ) + ## add predictors here to get random effect per region
                                      + (1 | RegionXSiteName), 
                                    memberships = list(Region = Wa, RegionXSiteName = Waj), 
                                    REML = T,
                                    data = allsyncx)

spatPlot2a <- sjPlot::plot_model(mem_mixed0, type="pred", terms=c("DistKM", "Connectivity"),
                                 axis.title = c("Euclidean Distance KM", "Thermal Trait Synchrony"), pred.type="re", 
                                 ci.lvl=NA, show.data  = T, dot.size = 0.2)

spatPlot2a 
file.name1 <- paste0(out.dir, "spatial_decay_mem_mod_exp_2_sp.jpg")
ggsave(spatPlot2a, filename=file.name1, dpi=300, height=8, width=10)


summary(mem_mixed0, ddf = "Satterthwaite")
anova(mem_mixed0, ddf = "Satterthwaite")
r2_nakagawa(mem_mixed0) 
check_singularity(mem_mixed0) ## false


### plots 
class(mem_mixed0) <- "lmerModLmerTest"
# sjPlot::plot_model(mem_mixed) 
ests <- sjPlot::plot_model(mem_mixed0, 
                           show.values=TRUE, show.p=TRUE,
                           title="Drivers of Thermal Synchrony")

ests
file.name1 <- paste0(out.dir, "effect_sizes_spatial_decay_2_sp.jpg")
ggsave(ests, filename=file.name1, dpi=300, height=8, width=10)


# Membership Model: Drivers ---------------------------------------


## model with new diversity measure

mem_mixed0 <- lmerMultiMember::lmer(Sync ~  (overlap*annual_avg + distance*annual_avg)*Connectivity
                                    + (1 | Region ) + ## add predictors here to get random effect per region
                                      (1 | RegionXSiteName), 
                                    memberships = list(Region = Wa, RegionXSiteName = Waj), 
                                    REML = T,
                                    data = allsyncx)

summary(mem_mixed0, ddf = "Satterthwaite")
anova(mem_mixed0, ddf = "Satterthwaite")
r2_nakagawa(mem_mixed0) 
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
file.name1 <- paste0(out.dir, "effect_sizes_diversity2_overlap_region_removed_2_sp.jpg")
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

quantile(allsyncx$annual_avg, probs = c(0.05,0.5,0.95))
round(quantile(allsyncx$overlap, probs = c(0.05,0.5,0.95)), digits = 3)
round(quantile(allsyncx$distance, probs = c(0.05,0.5,0.95)), digits = 3)


tempDiv <- sjPlot::plot_model(mem_mixed0, type="pred", terms= c("overlap",  "annual_avg [0.68, 0.92, 0.998]"),
                              axis.title = c("Trait Overlap", "Thermal Trait Synchrony"), 
                              legend.title = "Environmental synchrony", pred.type="re", 
                              ci.lvl=NA, show.data  = T, dot.size = 0.2)


tempDist <- sjPlot::plot_model(mem_mixed0, type="pred", terms= c("distance",  "annual_avg [0.68, 0.92, 0.998]"),
                               axis.title = c("Trait Distance", "Thermal Trait Synchrony"), 
                               legend.title = "Environmental synchrony", pred.type="re", 
                               ci.lvl=NA, show.data  = T, dot.size = 0.2)

Disttemp <- sjPlot::plot_model(mem_mixed0, type="pred", terms= c("annual_avg",  "distance [0.84, 0.98, 1.13]"),
                               axis.title = c("Environmental Synchrony", "Thermal Trait Synchrony"), 
                               legend.title = "Trait Distance", pred.type="re", 
                               ci.lvl=NA, show.data  = T, dot.size = 0.2)

tempCon <- sjPlot::plot_model(mem_mixed0, type="pred", terms= c("annual_avg", "Connectivity"),
                              axis.title = c("Environmental synchrony", "Thermal Trait Synchrony"), 
                              legend.title = "Connectivity", pred.type="re", 
                              ci.lvl=NA, show.data  = T, dot.size = 0.2)

divCon <- sjPlot::plot_model(mem_mixed0, type="pred", terms= c("overlap", "Connectivity"),
                             axis.title = c("Trait Overlap", "Thermal Trait Synchrony"), 
                             legend.title = "Connectivity", pred.type="re", 
                             ci.lvl=NA, show.data  = T, dot.size = 0.2)
divCon

distCon <- sjPlot::plot_model(mem_mixed0, type="pred", terms= c("distance", "Connectivity"),
                              axis.title = c("Trait Distance", "Thermal Trait Synchrony"), 
                              legend.title = "Connectivity", pred.type="re", 
                              ci.lvl=NA, show.data  = T, dot.size = 0.2)
distCon

tempDistCon <- sjPlot::plot_model(mem_mixed0, type="pred", terms= c("distance","Connectivity", "annual_avg [0.68, 0.92, 0.998]" ),
                                  axis.title = c("Trait Distance", "Thermal Trait Synchrony"), 
                                  legend.title = "Connectivity", pred.type="re")
tempDistCon

tempDistCon2 <- sjPlot::plot_model(mem_mixed0, type="pred", terms= c("annual_avg","Connectivity", "distance [0.02, 0.056, 0.187]" ),
                                   axis.title = c("Environmental synchrony", "Thermal Trait Synchrony"), 
                                   legend.title = "Connectivity", pred.type="re")

tempDivCon <- sjPlot::plot_model(mem_mixed0, type="pred", terms= c("overlap","Connectivity", "annual_avg [0.68, 0.92, 0.998]" ),
                                 axis.title = c("Trait Overlap", "Thermal Trait Synchrony"), 
                                 legend.title = "Connectivity", pred.type="re")

tempDivCon2 <- sjPlot::plot_model(mem_mixed0, type="pred", terms= c("annual_avg","Connectivity", "overlap [0, 0.035, 0.795]" ),
                                  axis.title = c("Environmental synchrony", "Thermal Trait Synchrony"), 
                                  legend.title = "Connectivity", pred.type="re")


tempDivCon2
?plot_model
## significant interactions
CompPlotSig <- plot_grid(list(tempDist,tempCon, distCon,  tempDistCon, tempDistCon2))

##insignificant interactions
# CompPlotInSig <- plot_grid(list(tempDiv, tempDist, distCon, tempDistCon)) ## grid figures

#terms=c("Euclid_Dist_Meters", "Connectivity")
file.name1 <- paste0(out.dir, "comp_mod_sig_interactions_region_removed_2_sp.jpg")
ggsave(CompPlotSig, filename=file.name1, dpi=300, height=12, width=18)

# file.name1 <- paste0(out.dir, "comp_mod_insig_interactions_region_removed_single_sp.jpg")
# ggsave(CompPlotInSig, filename=file.name1, dpi=300, height=10, width=15)
# 
# ### model selection 
# mem_mixed1 <- lmerMultiMember::lmer(Sync ~  annual_avg*(overlap+distance)  #Connectivity
#                                     + (1 | Region ) + ## add predictors here to get random effect per region
#                                       (1 | RegionXSiteName), 
#                                     memberships = list(Region = Wa, RegionXSiteName = Waj), 
#                                     REML = T,
#                                     data = allsyncx)

## add in interactions between other variables
#(annual_avg*log(diversity2) + distance*log(diversity2)) * Connectivity
# (annual_avg +distance+diversity2))*Connectivity
summary(mem_mixed1, ddf = "Satterthwaite")
anova(mem_mixed1, ddf = "Satterthwaite")
r2_nakagawa(mem_mixed1) 

### plots 
class(mem_mixed1) <- "lmerModLmerTest"
# sjPlot::plot_model(mem_mixed) 
ests <- sjPlot::plot_model(mem_mixed1, 
                           show.values=TRUE, show.p=TRUE,
                           title="Drivers of Thermal Synchrony")

ests

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

Disttemp <- sjPlot::plot_model(mem_mixed0, type="pred", terms= c("annual_avg",  "distance [0.84, 0.98, 1.13]"),
                               axis.title = c("Environmental Synchrony", "Thermal Trait Synchrony"), 
                               legend.title = "Trait Distance")

Divtemp <- sjPlot::plot_model(mem_mixed0, type="pred", terms= c("overlap",  "annual_avg [0.68, 0.92, 0.998]"),
                              axis.title = c("Environmental Synchrony", "Thermal Trait Synchrony"), 
                              legend.title = "Trait Overlap")

distCon <- sjPlot::plot_model(mem_mixed0, type="pred", terms= c("distance", "Connectivity"),
                              axis.title = c("Trait Distance", "Thermal Trait Synchrony"), 
                              legend.title = "Connectivity")




tempDistCon <- sjPlot::plot_model(mem_mixed0, type="pred", terms= c("distance","Connectivity", "annual_avg [0.68, 0.92, 0.998]" ),
                                  axis.title = c("Trait Distance", "Thermal Trait Synchrony"), 
                                  legend.title = "Connectivity")
tempDistCon

tempDistCon2 <- sjPlot::plot_model(mem_mixed0, type="pred", terms= c("annual_avg","Connectivity", "distance [0.02, 0.056, 0.187]" ),
                                   axis.title = c("Environmental synchrony", "Thermal Trait Synchrony"), 
                                   legend.title = "Connectivity")

tempDist
## significant interactions
CompPlotSig <- plot_grid(list(tempDist, distCon, tempDistCon, tempDistCon2)) ## grid figures

##insignificant interactions
CompPlotInSig <- plot_grid(list(tempDiv, tempDist, tempDistCon)) ## grid figures

#terms=c("Euclid_Dist_Meters", "Connectivity")
file.name1 <- paste0(out.dir, "comp_mod_sig_interactions_region_connectivity_removed_2_sp.jpg")
ggsave(CompPlotSig, filename=file.name1, dpi=300, height=12, width=18)

# file.name1 <- paste0(out.dir, "comp_mod_insig_interactions_region_removed_single_sp.jpg")
# ggsave(CompPlotInSig, filename=file.name1, dpi=300, height=10, width=15)
# 
# ?plot_model
# tempDist <- compPlot1[[2]]
# tempCon <- compPlot1[[3]]
# tempReg <- compPlot1[[4]]
# DivRegion <- compPlot1[[5]]
# DistRegion<- compPlot1[[6]]
# ConReg <- compPlot1[[7]]
# tempDivReg <- compPlot1[[8]]
# tempDistReg <- compPlot1[[9]]
# tempConReg <- compPlot1[[10]]
# 
# 
# 
# ## non significant interactions
# CompPlotnon <- plot_grid(list(tempCon, tempReg, ConReg, tempConReg))## grid figures
# 
# #terms=c("Euclid_Dist_Meters", "Connectivity")
# file.name1 <- paste0(out.dir, "comp_mod_interactions_non_significant_region_removed.jpg")
# ggsave(CompPlotnon, filename=file.name1, dpi=300, height=8, width=10)

