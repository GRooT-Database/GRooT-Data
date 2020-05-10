### R code for obtaining error risk for species for traits in the GRooT Database & ###
### to calculate calculate mean, median, first and third quantiles per species ###

# load data: GRooTFullVersion #
 
GRooTFullVersion<- read.csv("C:/Users/GRooTFullVersion.csv", header=T, na.strings=c("", "NA")) ####Final
str(GRooTFullVersion)
names(GRooTFullVersion)

GRooTFullVersion<-GRooTFullVersion[,c(1:71)]

require(tidyverse)

##################
### error risk ###
##################

### data entries with information only at genus level ###
### these information is not included to calculate error risks at species level ###
genusGRooT<-dplyr::filter(GRooTFullVersion,is.na(speciesTNRS)) 
genusGRooT$errorRiskEntries<-""
genusGRooT$errorRisk<-""

### data entries with information at species level ###
### these data was used to calculate error risk at species level ###
speciesGRooT<-dplyr::filter(GRooTFullVersion,!is.na(speciesTNRS)) 

names(speciesGRooT)

#scale_this <- function(x) as.vector(scale(x))

### error risks calculated for trait in which logarithmic transformation is required ###

speciesGRooTlog<- speciesGRooT %>% 
  dplyr::select(GRooTID, genusTNRS, speciesTNRS, traitName, traitValue) %>%
  group_by(genusTNRS, speciesTNRS, traitName) %>% 
  filter(traitName %in% c("Root_cortex_thickness", "Root_stele_diameter", "Root_stele_fraction", "Root_vessel_diameter",
                          "Root_branching_density", "Root_branching_ratio", "Root_C_N_ratio",
                          "Root_Ca_concentration", "Root_K_concentration", "Root_Mg_concentration",
                          "Root_Mn_concentration", "Root_N_concentration", "Root_N_P_ratio", "Root_P_concentration",
                          "Root_lifespan_mean", "Root_lifespan_median", "Root_litter_mass_loss_rate", "Root_production",
                          "Root_turnover_rate", "Mean_Root_diameter", "Root_dry_matter_content", "Root_tissue_density",
                          "Specific_root_area", "Specific_root_length", "Specific_root_respiration",
                          "Coarse_root_fine_root_mass_ratio", "Fine_root_mass_leaf_mass_ratio", "Root_length_density_volume",
                          "Root_mass_density", "Rooting_depth")) %>% 
  mutate(errorRiskEntries = n()) %>% 
  mutate(traitValuelog2=log2(traitValue+0.0001)) %>% ###0.0001 was added to include values = 0
  mutate(meanSpp=mean(traitValuelog2), sdSpp=sd(traitValuelog2)) %>%
  group_by(traitName) %>% 
  mutate(SDSppAvg=mean(sdSpp, na.rm =T)) %>%
  mutate(errorRisk=((meanSpp-traitValuelog2)/SDSppAvg)) %>%
  dplyr::select(GRooTID, genusTNRS, speciesTNRS, traitName, traitValue, errorRiskEntries, errorRisk)
  
####Error risk of zero means that only one observations is present for the that specific trait and species combination##

####Error risk for trait which follow a normal distribution###

speciesGRooTother<- speciesGRooT %>% 
  dplyr::select(GRooTID, genusTNRS, speciesTNRS, traitName, traitValue) %>%
  group_by(genusTNRS, speciesTNRS, traitName) %>% 
  filter(traitName %in% c("Root_xylem_vessel_number", "Root_mass_fraction", "Root_C_concentration",
                          "Root_lignin_concentration", "Root_total_structural_carbohydrate_concentration",
                          "Lateral_spread", "Root_mycorrhizal colonization", "Net_nitrogen_uptake_rate")) %>% 
  mutate(errorRiskEntries = n()) %>% 
  mutate(meanSpp=mean(traitValue), sdSpp=sd(traitValue)) %>%
  group_by(traitName) %>% 
  mutate(SDSppAvg=mean(sdSpp, na.rm =T)) %>%
  mutate(errorRisk=((meanSpp-traitValue)/SDSppAvg)) %>%
  dplyr::select(GRooTID, genusTNRS, speciesTNRS, traitName, traitValue, errorRiskEntries, errorRisk)

speciesRisk<-rbind(speciesGRooTlog, speciesGRooTother)

### Zero values for error risk are produced when only 1 data entry is available or ###
### when all data entries have the same value for the species ###

### merge error risk with other information in the database ###
speciesTotal<-merge(speciesGRooT, speciesRisk, by=c("GRooTID", "genusTNRS", "speciesTNRS", "traitName", "traitValue"))

### join the data at species and genus level ###
GRooTFull<-rbind(speciesTotal, genusGRooT) 

####change the order of columns#####

GRooTFull<-GRooTFull %>% arrange(GRooTID)

names(GRooTFull)

GRooTFullVersion <- GRooTFull[, c(1, 6:17, 2:3, 18:71, 4:5, 72:73 )]

names(GRooTFullVersion)


setwd("C:/Users/kiran/Dropbox/sROOT_mine/sROOT_database/")
write.csv(GRooTFullVersion, file = "GRooTFullVersionNew.csv",row.names=FALSE, na="")



########################################################################
### Calculate mean, median, first and third percentiles per species ####
########################################################################

###data entries which contain info at species level ###
speciesGRooT<-dplyr::filter(GRooTFullVersion,!is.na(speciesTNRS)) 

### Note 1: We calculated species by site values first to account for potential pseudo-replication and 
### variability in data entries' resolutions in GRooT (i.e., mean values versus individual observations).
### However, the user can calculate mean values using all the data entries by removing "studySite". 
### Note 2: Activate filter (by uncommenting L84 and 85 or L 101 and 102) is you want to: select for 
### belowground entities, error risk, or specific traits.

speciesGRooT$errorRisk<-as.numeric(speciesGRooT$errorRisk) 

###For trait that are not normal distributed, mean values can be calculated by using log transform values and back transform or by using means in the original units###

GRooTAggregateSpeciesVersion<- speciesGRooT %>% 
  ##filter(belowgroundEntities == "FR") %>% #if you are interested only in particular entities (other option below)
  ##filter(between(errorRisk, -4, 4)) %>% #if you want to filter by error risk values
  mutate(studySite= paste(referencesAbbreviated, decimalLatitude, decimalLongitud, locationID, location)) %>%
  dplyr::select(studySite, genusTNRS, speciesTNRS, traitName, traitValue) %>% 
  group_by(studySite, genusTNRS, speciesTNRS, traitName) %>% 
  filter(traitName %in% c("Root_cortex_thickness", "Root_stele_diameter", "Root_stele_fraction", "Root_vessel_diameter",
                          "Root_branching_density", "Root_branching_ratio", "Root_C_N_ratio",
                          "Root_Ca_concentration", "Root_K_concentration", "Root_Mg_concentration",
                          "Root_Mn_concentration", "Root_N_concentration", "Root_N_P_ratio", "Root_P_concentration",
                          "Root_lifespan_mean", "Root_lifespan_median", "Root_litter_mass_loss_rate", "Root_production",
                          "Root_turnover_rate", "Mean_Root_diameter", "Root_dry_matter_content", "Root_tissue_density",
                          "Specific_root_area", "Specific_root_length", "Specific_root_respiration",
                          "Coarse_root_fine_root_mass_ratio", "Fine_root_mass_leaf_mass_ratio", "Root_length_density_volume",
                          "Root_mass_density", "Rooting_depth")) %>%
  mutate(valuesLog=log(traitValue+0.0001)) %>%
  summarise(meanStudySite = mean(traitValue), meanStudySiteLog = mean(valuesLog)) %>% 
  #summarise(meanStudySiteLog = mean(valuesLog)) %>%
  group_by(genusTNRS, speciesTNRS, traitName) %>% 
  summarise(entriesStudySite = n(), meanSpecies = mean(meanStudySite), meanSpeciesExp = exp(mean(meanStudySiteLog)), 
            medianSpecies = median(meanStudySite), firstQuantile = quantile(meanStudySite, probs = c(0.25)),
            thirdQuantile = quantile(meanStudySite, probs = c(0.75)))  

###For trait normal distributed###

GRooTAggregateSpeciesVersion2<- speciesGRooT %>% 
  ##filter(belowgroundEntities == "FR") %>% #if you are interested only in particular entities (other option below)
  ##filter(between(errorRisk, -4, 4)) %>% #if you want to filter by error risk values
  mutate(studySite= paste(referencesAbbreviated, decimalLatitude, decimalLongitud, locationID, location)) %>%
  dplyr::select(studySite, genusTNRS, speciesTNRS, traitName, traitValue) %>% 
  group_by(studySite, genusTNRS, speciesTNRS, traitName) %>% 
  filter(traitName %in% c("Root_xylem_vessel_number", "Root_mass_fraction", "Root_C_concentration",
                          "Root_lignin_concentration", "Root_total_structural_carbohydrate_concentration",
                          "Lateral_spread", "Root_mycorrhizal colonization", "Net_nitrogen_uptake_rate")) %>%
  summarise(meanStudySite = mean(traitValue)) %>% 
  group_by(genusTNRS, speciesTNRS, traitName) %>% 
  summarise(entriesStudySite = n(), meanSpecies = mean(meanStudySite), 
            medianSpecies = median(meanStudySite), firstQuantile = quantile(meanStudySite, probs = c(0.25)),
            thirdQuantile = quantile(meanStudySite, probs = c(0.75)))  


#####option b for belowground entities ####

GRooTAggregateSpeciesVersion1<- speciesGRooT %>% 
  ##filter(between(errorRisk, -4, 4)) %>% ## if you want to filter error risk
  mutate(studySite= paste(referencesAbbreviated, decimalLatitude, decimalLongitud, locationID, location)) %>%
  select(studySite, belowgroundEntities, genusTNRS, speciesTNRS, traitName, traitValue, errorRisk) %>% 
  group_by(studySite, belowgroundEntities, genusTNRS, speciesTNRS, traitName) %>% 
  #filter(traitName %in% c("Root_N_concentration", "Mean_Root_diameter", "Root_dry_matter_content",
  #                        "Root_tissue_density", "Specific_root_length")) %>%
  summarise(meanStudySite = mean(traitValue)) %>% 
  group_by(belowgroundEntities, genusTNRS, speciesTNRS, traitName) %>% 
  summarise(entriesStudySite = n(), meanSpecies = mean(meanStudySite), 
            medianSpecies = median(meanStudySite), firstQuantile = quantile(meanStudySite, probs = c(0.25)),
            thirdQuantile = quantile(meanStudySite, probs = c(0.75)))  







