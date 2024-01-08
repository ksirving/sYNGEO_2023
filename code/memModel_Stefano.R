## mulimembership model code 
# https://github.com/jvparidon/lmerMultiMember

## packages

library(lmerMultiMember) ## membership model package
library(lmerTest) ## from Luke et al 2016 - evaluating significance in linear mixed-effects models in R
library(lme4)
library(sjstats) #use for r2 functions
library(performance) ## for ICC
library(sjPlot) #for plotting lmer and glmer mods


## create memberships
Wa <- lmerMultiMember::weights_from_vector(allsyncxWithin$Region) ## memberhip of region
Wj <- Matrix::fac2sparse(allsyncxWithin$SiteName)  # convert single membership vars to an indicator matrix with fac2sparse() - SiteName is one of the site pairs
Waj <- interaction_weights(Wa, Wj)

mem_mixed <- lmerMultiMember::lmer(Sync ~ WCDistkmsqrt
                                     + (1 | Region ) + ## add predictors here to get random effect per region
                                       + (1 | RegionXSiteName), 
                                     memberships = list(Region = Wa, RegionXSiteName = Waj), 
                                     REML = T,
                                     data = allsyncxWithin)

## coefs
summary(mem_mixed, ddf = "Satterthwaite")
anova(mem_mixed, ddf = "Satterthwaite")
r2_nakagawa(mem_mixed) 
check_singularity(mem_mixed)
performance::icc(mem_mixed, by_group = T) 

## trick lmer into think it's an lme model to use sjplot
class(mem_mixedwc) <- "lmerModLmerTest"

## estimates
estsdw <- sjPlot::plot_model(mem_mixed, 
                             show.values=TRUE, show.p=TRUE,
                             title="Drivers of Thermal Synchrony")

## example relationships
tempDis <- sjPlot::plot_model(mem_mixed, type="pred", terms= c("ZDistance", "ZTempCor [-1.86, 0.3, 0.92]"),
                              axis.title = c( "Trait Distance (Z Score)","Thermal Trait Synchrony"), 
                              legend.title = "Environmental synchrony")
