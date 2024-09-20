### synchrony on single traits
# install.packages("overlapping")
## packages
library(tidyverse)
library(reshape2)
library(tidyr)
library(dplyr)
library(synchrony)
library(codyn)
# library(rfishbase)
library(munfold)
library(data.table)
library(gdata)
library(here)
library(tidylog)
library(overlapping)
# 
# update.packages("bayestestR")
# install.packages("bayestestR")

getwd()
## upload fish abundance and site data

originaldata <- read.csv("input_data/Bio/fishdata_selection_basins_same_time_window_10262020.csv")
head(originaldata)

# test <- originaldata %>% filter(SiteID == "S6654")
# test

## upload <=3 species sites to remove

load(file="output_data/01a_sites_3_species.RData")
head(singleSp3l)

## upload and format community weighted mean traits - all groups

trait_matrix <- read.csv("output_data/01_trt_single_traits_interpolated_cwm_cmv.csv")
head(trait_matrix)
## combine all groups

all_groups <- trait_matrix %>%
  mutate(BiogeoRegion = case_when(Country %in% c("FIN", "SWE", "GBR", "ESP", "FRA") ~ "Europe",
                                  Country== "AUS" ~ "Oceania", Country == "USA" ~ "USA")) %>%
  filter(!SiteID %in% singleSp3l)
head(all_groups)

## temp preference per region

all_groups %>% group_by(BiogeoRegion) %>% summarise(MeanTemp = mean(CWM))
names(all_groups)

## remove <3 species sites
all_groups

## save out sites
# finalSites <- all_groups %>%
#   select(HydroBasin:Longitude, BiogeoRegion) %>%
#   distinct()
# 
# write.csv(finalSites, "input_data/bio/02_Sites.csv")
  
source("https://raw.githubusercontent.com/alaindanet/fishcom/master/R/synchrony.R")


# between all sites - within biogeogrpahic region -----------------------

## define biogeogrpahic region
regionsID<-unique(all_groups$BiogeoRegion) # 3 regions
synchronyx = NULL
region = 1
ax=1


  ### loop over regions
  for (region in 2:length(regionsID)) {
    
    print(regionsID[region])
    
    basindata<-all_groups[all_groups$BiogeoRegion==regionsID[region],]
    # head(basindata)
    basindata <- basindata[order(basindata$SiteID),]
    
    ### loop over axis
    Ntraits<-unique(basindata$Trait)
 
    
    for (ax in 1: length(Ntraits)) {
      
      
      trait_data<-basindata[basindata$Trait==unique(basindata$Trait)[ax],]
      # sum(is.na(trait_data))
      # head(trait_data)
      years <- unique( trait_data$Year)
      # years
      
      # make df wide - mean
      trait_CWM  <- trait_data %>% 
        dplyr::select(-c(X, site_year, Country, HydroBasin, CWV,CV, Latitude, Longitude) ) %>%
        spread(SiteID, CWM) #%>% ## some NAs appear here, don't have all trait scores for all site years
      # head( trait_CWM)
      # make df wide - variance
      trait_CWV  <- trait_data %>% 
        dplyr::select(-c(X, site_year, Country, HydroBasin, CWM,CV,  Latitude, Longitude) ) %>%
        spread(SiteID, CWV) #%>% ## some NAs appear here, don't have all trait scores for all site years
      # head(trait_CWV)
      # # make df wide - CV
      # trait_CV  <- trait_data %>% 
      #   dplyr::select(-c(X, site_year, Country, HydroBasin, CWM,CWV,  Latitude, Longitude) ) %>%
      #   spread(SiteID, CV) #%>% ## some NAs appear here, don't have all trait scores for all site years
      # # head(trait_CV)
      # # remove non value columns
      trait_CWM <- (trait_CWM)[,-c(1:3)]
      trait_CWV <- (trait_CWV)[,-c(1:3)]
      # trait_CV <- (trait_CV)[,-c(1:3)]
   
      ### synchrony 
      cc <- expand.grid(colnames(trait_CWM), colnames(trait_CWM), KEEP.OUT.ATTRS = FALSE)
      cc
      synchrony <- sapply(seq_len(nrow(cc)), function(k) {
       
        print("synchrony")
        print(k)
        
        
        i <- cc[k,1]
        j <- cc[k,2]
        
        sync_mat <- matrix(
          c(trait_CWM[, i],trait_CWM[,j]),
          nrow = 10,
          byrow = F,
          dimnames = list(years,
            c("site1", "site2")
          )
        )
        # sync_mat
        compute_synchrony(cov(sync_mat))
      })
      
      # synchrony
      # head(synchrony)
      
      ## diversity 2: overlapping index on CWM
      # diversity2 <- sapply(seq_len(nrow(cc)), function(k) {
      #   
      #   print("diversity2")
      #   print(k)
      #   
      #   i <- cc[k,1]
      #   j <- cc[k,2]
      #   
      #   ## format data as list
      #   div_mat2 <- list(site1 = trait_CWM[, i], site2 = trait_CWM[,j])
      #   
      #   # sync_mat
      #   overlapping::overlap(div_mat2)$OV #, type = "2"
      #   # plot(bayestestR::overlap(div_mat2$site1, div_mat2$site2))
      # })
      # 
      # diversity2
      # diversity2 - v2
      
      ## overlap using type =2
      # diversity3 <- sapply(seq_len(nrow(cc)), function(k) {
      # 
      #   print("diversity3")
      #   print(k)
      #   
      #   i <- cc[k,1]
      #   j <- cc[k,2]
      # 
      #   ## format data as list
      #   div_mat2 <- list(site1 = scale(trait_CWM[, i]), site2 = scale(trait_CWM[,j]))
      # 
      #   # bayestestR::overlap(div_mat2$site1, div_mat2$site2) ## overalp index
      #   overlapping::overlap(div_mat2)$OV #, type = "2"
      #   # bayestestR::overlap(div_mat2$site1,div_mat2$site2, method_density = "KernSmooth")[[1]]
      #   # plot(bayestestR::overlap(div_mat2$site1, div_mat2$site2))
      # })
      
      # ggplot(test, aes(x = value,fill = name)) +
      #   geom_density( alpha = 0.7) #+
      #   labs(title = "Kernel Density Plot of Salary",
      #        x = "Salary",
      #        y = "Density")
      
      diversityBayes <- sapply(seq_len(nrow(cc)), function(k) {
        
        print("diversityBayes")
        print(k)
        
        i <- cc[k,1]
        j <- cc[k,2]
        
        ## format data as list
        div_mat2 <- list(site1 = trait_CWM[, i], site2 = trait_CWM[,j])
        # div_mat2
   
        bayestestR::overlap(scale(div_mat2$site1), scale(div_mat2$site2), method_density = "KernSmooth")[[1]] ## overalp index
        # plot(bayestestR::overlap(scale(div_mat2$site1), scale(div_mat2$site2), method_auc = "spline"))
        # bayestestR::overlap(scale(div_mat2$site1), scale(div_mat2$site2), method_density = "KernSmooth", method_auc = "spline", extend=F)
        # overlapping::overlap(div_mat2, type = "2")$OV #, type = "2"
        # plot(bayestestR::overlap(div_mat2$site1, div_mat2$site2))
      })

      # cor(diversity3, diversity2) ## negative
      ### distance: Difference in temporal average of Community Weighted Mean
      
      # cc <- expand.grid(colnames(trait_CWM), colnames(trait_CWM), KEEP.OUT.ATTRS = FALSE)

      distance <- sapply(seq_len(nrow(cc)), function(k) {
        
        print("distance")
        print(k)
        
        i <- cc[k,1]
        j <- cc[k,2]
        
        dist_mat <- matrix(
          c(trait_CWM[, i],trait_CWM[,j]),
          nrow = 10,
          byrow = F,
          dimnames = list(years,
                          c("site1", "site2")
          )
        )
        # dist_mat
        abs((mean(dist_mat[,1])-mean(dist_mat[,2]))/mean(dist_mat))
      
      })
      
  
      ## distance without dividing by mean
      distanceAbs <- sapply(seq_len(nrow(cc)), function(k) {
        
        print("distanceAbs")
        print(k)
        
        i <- cc[k,1]
        j <- cc[k,2]
        
        dist_mat <- matrix(
          c(trait_CWM[, i],trait_CWM[,j]),
          nrow = 10,
          byrow = F,
          dimnames = list(years,
                          c("site1", "site2")
          )
        )
        
        mean((dist_mat[,1])-mean(dist_mat[,2]))#/mean(dist_mat)
        
      })
     
      # round(cor(distance, distanceAbs, use = "pairwise"), digits = 2)

      ## combine all
      synchrony <- cbind(cc, synchrony, distance, distanceAbs, diversityBayes)
      # synchrony
      ## add traits and region
      synchrony <- synchrony %>%
        mutate(Trait = Ntraits[ax], Region = regionsID[region])
   
      ## add tio main DF
      synchronyx <- rbind(synchronyx, synchrony)
      
      
    }
    
    
  }
  
  
synchrony_axis <- synchronyx %>%
  rename(Site_ID1 = Var1, Site_ID2 = Var2) %>%
  mutate(Pair = paste0(Site_ID1, ".", Site_ID2))

## negative values in distance

head(synchrony_axis)
length(unique(synchrony_axis$Trait))
length(unique(synchrony_axis$Region))

## add to other overlap metrics
synchrony_axisX <- read.csv("sync/02_between_all_sites_single_traits_CWM_CWV_CV.csv")

synchrony_axisX <- synchrony_axisX %>%
  select(diversity2, diversity3)

## join together
synchrony_axis <- cbind(synchrony_axis, synchrony_axisX)


head(synchrony_axis)


###save results
write.csv(synchrony_axis, "sync/02_between_all_sites_single_traits_CWM_CWV_CV.csv")

## checking diverasity2 nas
# synchrony_axis <- read.csv( "output_data/sync/02_between_all_sites_single_traits_CWM_CWV_CV.csv")
# head(synchrony_axis)
# 
# ind <- which(is.na(synchrony_axis$diversity2))
# length(ind) ## 62697
# test <- synchrony_axis[ind,]
# 
# test$Pair[1:10]

# Climate synchrony-------------------------------------------------------------------------

melt_clim_raw <- read.csv(file="input_data/Env/air_annual_and_summer_avg.csv") %>%
  rename(SiteID = siteid)

head(melt_clim_raw)

#load(file="input_data/Env/flow_data_melt_raw_new_sites.RData")

# S10203.S10089 e.g. missing pair, both sites missing 

sites <- all_groups %>%
  select(SiteID, HydroBasin, Country) %>%
  distinct()

sites <- sites %>%
  mutate(BiogeoRegion = case_when(Country %in% c("FIN", "SWE", "GBR", "ESP", "FRA") ~ "Europe",
                                  Country== "AUS" ~ "Oceania", Country == "USA" ~ "USA"))

tmean <- read.csv("/Users/katieirving/Library/CloudStorage/OneDrive-SCCWRP/Documents - Katie’s MacBook Pro/git/sYNGEO_Env_Extraction/Output/Sites_tmean_av_new_sites_2004_2013.csv") %>%
  pivot_longer(X2004:X2013, names_to = "Year", values_to = "temp") %>% mutate(metric = "annual_avg") %>%
  mutate(Year = gsub ( "X", "", Year))

head(tmean)

# tmin <- read.csv("/Users/katieirving/Library/CloudStorage/OneDrive-SCCWRP/Documents - Katie’s MacBook Pro/git/sYNGEO_Env_Extraction/Output/Sites_tmin_av_new_sites_2004_2013.csv")
# tmax <- read.csv("/Users/katieirving/Library/CloudStorage/OneDrive-SCCWRP/Documents - Katie’s MacBook Pro/git/sYNGEO_Env_Extraction/Output/Sites_tmax_av_new_sites_2004_2013.csv")

head(sites)
dim(sites)

## join sites to climate data
clim_sites <- full_join(tmean, sites, by="SiteID")

unique(clim_sites$BiogeoRegion)

clim_sites <- na.omit(clim_sites)
dim(clim_sites)
length(unique(clim_sites$SiteID)) 


regionsID<-unique(clim_sites$BiogeoRegion) # 3 regions
synchronyx = NULL
region = 2
ax = 1

regionsID


### loop over regions
for (region in 1:length(regionsID)) {
  
  basindata<-clim_sites[clim_sites$BiogeoRegion==regionsID[region],]
  # head(basindata)
  basindata <- basindata[order(basindata$SiteID),]

  ### loop over axis
  Ntraits<-unique(basindata$metric)
  Ntraits
  
  for (ax in 1: length(Ntraits)) {

    trait_data<-basindata[basindata$metric==unique(basindata$metric)[ax],]
    years <- unique(trait_data$Year)
    
    # make df wide - mean
    trait_temp  <- trait_data %>% 
      dplyr::select(-c( Country, HydroBasin) ) %>%
      spread(SiteID, temp) #%>% ## some NAs appear here, don't have all trait scores for all site years

    # format data
    trait_temp <- as.data.frame(trait_temp[,-c(1:3)])
    trait_temp
    ### synchrony 

    cc <- expand.grid(colnames(trait_temp), colnames(trait_temp), KEEP.OUT.ATTRS = FALSE)

    synchrony <- sapply(seq_len(nrow(cc)), function(k) {
      i <- cc[k,1]
      j <- cc[k,2]
 
      sync_mat <- matrix(
        c(trait_temp[, i],trait_temp[,j]),
        nrow = 10,
        byrow = F,
        dimnames = list(years,
                        c("site1", "site2")
        )
      )
      sync_mat
      compute_synchrony(cov(sync_mat))
      

    })
    
    corCoef <- sapply(seq_len(nrow(cc)), function(k) {
      i <- cc[k,1]
      j <- cc[k,2]
      
      sync_mat <- matrix(
        c(trait_temp[, i],trait_temp[,j]),
        nrow = 10,
        byrow = F,
        dimnames = list(years,
                        c("site1", "site2")
        )
      )
  
      cor(sync_mat)[1,2]
      ?cor
    })
    
    synchrony<- cbind(cc,synchrony, corCoef) 
    synchrony
    
    ## add traits and region
    synchrony <- synchrony %>%
      mutate(Metric = Ntraits[ax], Region = regionsID[region])

    ## add to main DF
    synchronyx <- rbind(synchronyx, synchrony)
    
    
  }
  
  
}

synchrony_axis <- synchronyx %>%
  rename(Site_ID1 = Var1, Site_ID2 = Var2) %>%
  mutate(Pair = paste0(Site_ID1, ".", Site_ID2))

## negative values in distance

head(synchrony_axis)
length(unique(synchrony_axis$Metric))
length(unique(synchrony_axis$Region))


###save results
write.csv(synchrony_axis, "sync/02_between_all_sites_temp_synchrony.csv")

