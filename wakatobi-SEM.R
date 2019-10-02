# PATH ANALYSIS of Wakatobi SES
# Requires permission (shared Google drive folder) and access (internet) to Wakatobi Data

rm(list=ls())

#update.packages(ask = FALSE, checkBuilt = TRUE)

library(piecewiseSEM)
library(nlme) # Version 3.1.122
library(googledrive)
library(reshape2)
library(ggplot2)
library(corrplot)
library(car) #VIF calculations
library(vegan)
library(codyn) #Simpson's evenness calculation


# FIRST: setwd for where you want outputs saved: 
outdir<-"~/Analyses/_RESULTS/SES-SEM/"


# Data files are in GoogleDrive
drive_auth() # Will require you to sign into Google account and grant permission to tidyverse for access 


# input / munge fish data
# get file ID from Google Drive's "shareable link" for the file: https://drive.google.com/open?id=12-DNIlHdoVT2JiWpoF2fu-Drl28RzjVG
drive_download(as_id("12-DNIlHdoVT2JiWpoF2fu-Drl28RzjVG"), overwrite=TRUE) # Saves file to working directory 
fishdat<-read.csv("cleaned_wakatobi_fish_uvc.csv") 
file.remove("cleaned_wakatobi_fish_uvc.csv") # Now that it's loaded into R, can delete file that was just downloaded



# RESPONSE VARIABLES:

# Calculate total fish biomass at site level (site_id column), averaged across three transects
fish.mass1<-aggregate(biomass_g ~ site_id + transect, data=fishdat, FUN=sum)
fish.mass<-aggregate(biomass_g ~ site_id, data=fish.mass1, FUN=mean)

# Calculate fish biomass by functional group at site level, averaged across three transects
#fun.mass1<-aggregate(biomass_g ~ site_id + transect + trophic_group, data=fishdat, FUN=sum)
#fun.mass<-aggregate(biomass_g ~ site_id + trophic_group, data=fun.mass1, FUN=mean)


# Calculate species diversity and test all of them
# See Morris et al. (Ecology and Evolution) for discussion simultaneously considering analyses
# of multiple indices can provide greater insight

# 1 - richness 
spcount.transect<-aggregate(scientific_name ~ site_id + transect, data=fishdat, FUN=length)
names(spcount.transect)[3]<-"no_of_species" 
# Note: no_of_species is a bit of a misnomer; most likely an undercount. EXAMPLE: "Acanthurs spp" could be used for more than one species, but only counted as one species here
fish.rich<-aggregate(no_of_species ~ site_id, data=spcount.transect, FUN=mean)

# 2- shannon diversity (aka H') - using "TROPHIC_GROUP" or FUNCTIONAL DIVERSITY
### need to re-check below for diversity calculation
countPerSp.transect<-aggregate(number_of_fish ~ site_id + trophic_group + transect, data=fishdat, FUN=sum)
countPerSp.site<-aggregate(number_of_fish ~ site_id + trophic_group, data=countPerSp.transect, FUN=mean)
# Convert to matrix for calculating diversity:
countPerSp.mat<-acast(countPerSp.site, site_id~trophic_group, value.var = "number_of_fish", FUN=sum)
# Replace NAs with 0s
countPerSp.mat[is.na(countPerSp.mat)]<-0
fish.shan<-diversity(countPerSp.mat, index="shannon")
fish.shan<-as.data.frame(cbind(as.numeric(names(fish.shan)), fish.shan))
names(fish.shan)<-c("site_id", "shannon")


# 3 - use INVERSE of simpsons diversity (aka D2), more commonly used than original D1 index (see Morris et al. 2014)
fish.isim<-diversity(countPerSp.mat, index="invsimpson")
fish.isim<-as.data.frame(cbind(as.numeric(names(fish.isim)), fish.isim))
names(fish.isim)<-c("site_id", "invsimpson")

# 4 - calculate Simpson's Evenness (use codyn package) - requires dataframe of counts (not matrix)
fish.even<-community_structure(countPerSp.site, abundance.var="number_of_fish", replicate.var="site_id", metric="SimpsonEvenness")
fish.even<-fish.even[,-2] # Remove richness column

### OTHER POTENTIAL RESPONSE VARIABLES: functional trait diversity using "fun.mass" above?


# NEXT: plot histograms of all response variables, consider log transforming some variables

# Total biomass
setwd(outdir)
pdf(file="plot_histogram.totalbiomass.pdf")
hist(fish.mass[,"biomass_g"],xlab="Biomass", main="Histogram of Site-Level Biomass")
dev.off()

# Try log biomass
pdf(file="plot_histogram.LOGtotalbiomass.pdf")
hist(log10(fish.mass[,"biomass_g"]),xlab="log Biomass", main="Histogram of Site-Level log Biomass")
dev.off()

# Option to use log biomass as response variable:
fish.logmass<-fish.mass
fish.logmass$biomass_g<-log10(fish.logmass[,"biomass_g"])
names(fish.logmass)[2]<-"log_biomass_g"

# Richness
setwd(outdir)
pdf(file="plot_histogram.richness.pdf")
hist(fish.rich[,"no_of_species"],xlab="Richness", main="Histogram of Site-Level Species Richness")
dev.off()

# Shannon
setwd(outdir)
pdf(file="plot_histogram.shannon.pdf")
hist(fish.shan[,"shannon"],xlab="Shannon Diversity", main="Histogram of Site-Level Diversity")
dev.off()

# inverse Simpson
setwd(outdir)
pdf(file="plot_histogram.invsimpson.pdf")
hist(fish.isim[,"invsimpson"],xlab="Inverse Simpson's Diversity", main="Histogram of Site-Level Diversity")
dev.off()

# Simpson's Evenness
setwd(outdir)
pdf(file="plot_histogram.evenness.pdf")
hist(fish.even[,"SimpsonEvenness"],xlab="Simpson's Evenness", main="Histogram of Site-Level Diversity")
dev.off()

# input / munge benthic cover data: https://drive.google.com/open?id=1ba04__uY3alCvHNXI1CLmQmInstLwach
drive_download(as_id("1ba04__uY3alCvHNXI1CLmQmInstLwach"), overwrite=TRUE) # Saves file to working directory 
coraldat<-read.csv("raw_coral_data_wakatobi_may_2018.csv") 
file.remove("raw_coral_data_wakatobi_may_2018.csv")

table(coraldat$Dive.Site, coraldat$Transect) # 3 Transects per site (100 measurements per transect)

coral.tmp<-aggregate(Life.Form ~ Dive.Site, data=coraldat, FUN = table)
#write.csv(agg.tmp, file="benthicCountsTable.csv", row.names=FALSE, quote = FALSE)

# CLEANING:
# cm -> CM
# s -> S
coraldat$Life.Form[coraldat$Life.Form=="s"]<-"S"
coraldat$Life.Form[coraldat$Life.Form=="cm"]<-"CM"
coraldat$Life.Form<-factor(coraldat$Life.Form) # Reset levels


agg.clean<-aggregate(Life.Form ~ Dive.Site, data=coraldat, FUN = table)
setwd(outdir)
write.csv(agg.clean, file="data_wakatobi_benthicCountsTable-allcategories.csv", row.names=FALSE, quote = FALSE)

site.counts<-as.data.frame(table(coraldat$Dive.Site))
benth.site<-(agg.clean[,-1] / site.counts$Freq)
rownames(benth.site)<-as.character(agg.clean$Dive.Site)

# Percent Cov Key:
# ACB -   Acropora branching
# ACD - Acropora digitate
# ACE -    Acropora encrusting
# ACT -   Acropora tabulate
# CB -   Branching coral
# CE -   Coral encrusting
# CF -    Foliose coral
# CHE - Heliopora
# CM - Massive coral
# CME - Millephora
# CMR - Mushroom coral
# CS - Submassive coral
# DCA -   Dead coral with algae
# MA -   Macro Algae
# OT -    Others
# R -    Ruble
# RCK - Rock
# S - Sand
# SC - Soft Coral
# SP - Sponge



# AGGREGATE All Hard Corals and R+RCK+S (abiotic)
hardCor<-c("ACB", "ACD", "ACE", "ACT", "CB", "CE", "CF", "CHE", "CM", "CME", "CMR", "CS")
hard.set<-benth.site[,colnames(benth.site) %in% hardCor]
hard.sum<-rowSums(hard.set)

# For algae, keep "Dead coral with Algae" separate from "Macroalgae" - ie, don't sum
#allAlg<-c("DCA", "MA")
#alg.set<-percentcov[,colnames(percentcov) %in% allAlg]
#alg.sum<-rowSums(alg.set)

allAbio<-c("R", "RCK", "S")
abio.set<-benth.site[,colnames(benth.site) %in% allAbio]
abio.sum<-rowSums(abio.set)

##############################################################################
###### FOR NOW, remove:
# (1) all taxa and morphology-specific coral data (i.e., all columns from "ACB" (column 1) to "CS" (column 12))
# (2) OTHER (very little variation)
# (3) Sponge - not interested
benth.site<-subset(benth.site, select=c(DCA, MA, R, RCK, SC, S))


benthcov.site<-cbind(benth.site, hard.sum, abio.sum)
colnames(benthcov.site)[(ncol(benthcov.site)-1):(ncol(benthcov.site))]<-c("All_HardCoral", "All_Abiotic")
setwd(outdir)
write.csv(benthcov.site, "data_wakatobi_benthicPercentCover-allcategories.csv", quote = FALSE)


# input / munge rugosity data: https://drive.google.com/open?id=1bB-UTzGzF2CEJhrUsJH9xxZ067khCqgw
drive_download(as_id("1bB-UTzGzF2CEJhrUsJH9xxZ067khCqgw"), overwrite=TRUE)
rugdat<-read.csv("raw_rugosity_data_wakatobi_may_2018.csv")
file.remove("raw_rugosity_data_wakatobi_may_2018.csv")
rug.site<-aggregate(Rugosity ~ Site.Name, data=rugdat, FUN=mean)
setwd(outdir)
write.csv(rug.site, "data_wakatobi_benthicRugosity.csv", quote=FALSE, row.names=FALSE)


# input human population data: https://drive.google.com/open?id=1DcVqeVEx6yGksBqLGzbWhU8WxN6UcYBz

#distWeighted.dat<-read.csv("_humanPopData/data_wakatobiHumans_distanceWeighted.csv") # LEAVE THESE OUT FOR NOW
#distToLandings.dat<-read.csv("_humanPopData/data_wakatobiHumans_distanceToLandingsSite.csv") # LEAVE THESE OUT FOR NOW
drive_download(as_id("1DcVqeVEx6yGksBqLGzbWhU8WxN6UcYBz"), overwrite=TRUE) # Saves file to working directory 
humanDensity.dat<-read.csv("data_wakatobiHumans_areaWeightedDensityMetrics_5_km_buffer.csv") # weights each village's population density by its area to get "total population" within 5km buffer
file.remove("data_wakatobiHumans_areaWeightedDensityMetrics_5_km_buffer.csv")




# input oceanographic (and other) variables from MSEC: https://drive.google.com/open?id=12CErWykopoj2_gpQI47XYHbUEocOdOSr
drive_download(as_id("12CErWykopoj2_gpQI47XYHbUEocOdOSr"), overwrite=TRUE) # Saves file to working directory 
msec.dat<-read.csv("msec_out_5km.csv")
file.remove("msec_out_5km.csv")


# input SST data: https://drive.google.com/open?id=1ROPUFf6yi6r78vw9eTOa78WyqxFfYuBK
drive_download(as_id("1ROPUFf6yi6r78vw9eTOa78WyqxFfYuBK"), overwrite=TRUE) # Saves file to working directory 
sst.dat<-read.csv("Wakatobi_2018_SSTExtract.csv")
file.remove("Wakatobi_2018_SSTExtract.csv")


# input cleaned fish flow data: https://drive.google.com/open?id=1PRrdjBQ-aWsjKwO5gHV91VXKZ4cRxE_G
drive_download(as_id("1PRrdjBQ-aWsjKwO5gHV91VXKZ4cRxE_G"), overwrite=TRUE) # Saves file to working directory 
trip.dat<-read.csv("Wakatobi-landings_201909124_TRIP_cleanedFishFlows.csv")
file.remove("Wakatobi-landings_201909124_TRIP_cleanedFishFlows.csv")

# Trim down dataset
trip.dat<-subset(trip.dat, select=c(trip_id, fishing_grnd1,
                          fishflow_abund,
                          landings_sold_personally_abund,
                          landings_sold_Papalele_abund,
                          landings_sold_Pengumpul_abund,
                          landings_eaten_abund,
                          landings_given_abund))




## input aggregation file for landings trips: https://drive.google.com/open?id=1PkaXlA1r1RA6tUWX7Tm3sk3SPm7kJxMf
drive_download(as_id("1PkaXlA1r1RA6tUWX7Tm3sk3SPm7kJxMf"), overwrite=TRUE)
trip.agg<-read.csv("aggregationKey-FishingGround_PC.csv")
file.remove("aggregationKey-FishingGround_PC.csv")


# Aggregate (calculate mean) of groups of fishing grounds based on column: new_fg
# First remove trips with no fishing ground information (n=4)
trip.dat[!(trip.dat$fishing_grnd1 %in% trip.agg$original_fg),]
trip.dat<-trip.dat[(trip.dat$fishing_grnd1 %in% trip.agg$original_fg),]

# Change column name
names(trip.dat)[grep("fishing_grnd1", names(trip.dat))]<-"original_fg"
trip.dat<-merge(trip.dat, trip.agg, by="original_fg")

table(trip.dat$new_fg)

# How do total landings affect fish response
landings_sumtot<-aggregate(trip.dat$fishflow_abund ~ trip.dat$new_fg, FUN = sum )
names(landings_sumtot)<-c("location", "landings_sum_tot")
# Doesn't seem like an appropriate driver, needs to be normalized by area of fishing ground

landings_meantot<-aggregate(trip.dat$fishflow_abund ~ trip.dat$new_fg, FUN = mean )
names(landings_meantot)<-c("location", "landings_mean_tot")
table(trip.dat$new_fg)
# Better measure of fishing pressure than sum total?
# If area of fishing ground is proportional to number of trips, then yes
# Or if assumption is that fishing EFFORT is equal everywhere (Badjao fish everywhere equally), then yes 

# How do landings eaten/given vs sold (to anyone) affect fish response
# Note: ZERO landings eaten or given in this dataset



# How do "subsistence" landings (eaten, given, sold_personally - ie, no middle-person?) affect fish response
landings_personal<-aggregate(trip.dat$landings_sold_personally_abund +
                         trip.dat$landings_eaten_abund +
                         trip.dat$landings_given_abund ~ trip.dat$new_fg, FUN = sum )
names(landings_personal)<-c("location", "landings_sum_personal")
landings_personal<-merge(landings_personal, landings_sumtot, by="location")
landings_personal$landings_personal_prop<-landings_personal$landings_sum_personal / landings_personal$landings_sum_tot
landings_personal<-subset(landings_personal, select=-landings_sum_tot)

# How do landings sold to pengumpul (GO OFF ISLAND) affect fish response? - include sum, proportion, and mean off-island catch
landings_pengumpul<-aggregate(trip.dat$landings_sold_Pengumpul_abund ~ trip.dat$new_fg, FUN = sum )
names(landings_pengumpul)<-c("location", "landings_sum_pengumpul")
landings_pengumpul<-merge(landings_pengumpul, landings_sumtot, by="location")
landings_pengumpul$landings_pengumpul_prop<-landings_pengumpul$landings_sum_pengumpul / landings_pengumpul$landings_sum_tot

landings_pengumpul_mean<-aggregate(trip.dat$landings_sold_Pengumpul_abund ~ trip.dat$new_fg, FUN = mean )
names(landings_pengumpul_mean)<-c("location", "landings_mean_pengumpul")
landings_pengumpul<-merge(landings_pengumpul, landings_pengumpul_mean, by="location")

landings_pengumpul<-subset(landings_pengumpul, select=-landings_sum_tot)


# How do landings sold to papalele affect fish response?
landings_papalele<-aggregate(trip.dat$landings_sold_Papalele_abund ~ trip.dat$new_fg, FUN = sum )
names(landings_papalele)<-c("location", "landings_sum_papalele")
landings_papalele<-merge(landings_papalele, landings_sumtot, by="location")
landings_papalele$landings_papalele_prop<-landings_papalele$landings_sum_papalele / landings_papalele$landings_sum_tot
landings_papalele<-subset(landings_papalele, select=-landings_sum_tot)
# Papalele = sold ON-ISLAND

# How do landings sold to pengumpul or papalele affect fish response?
landings_market<-aggregate(trip.dat$landings_sold_Papalele_abund +
                             trip.dat$landings_sold_Pengumpul_abund ~ trip.dat$new_fg, FUN=sum)
names(landings_market)<-c("location", "landings_sum_market")
landings_market<-merge(landings_market, landings_sumtot, by="location")
landings_market$landings_market_prop<-landings_market$landings_sum_market / landings_market$landings_sum_tot 
landings_market<-subset(landings_market, select=-landings_sum_tot)

# How do landings sold to papalele (SOLD ON ISLAND) affect fish response? - include sum, proportion, and mean off-island catch
landings_onisland<-aggregate(trip.dat$landings_sold_personally_abund +
                               trip.dat$landings_eaten_abund +
                               trip.dat$landings_given_abund +
                               trip.dat$landings_sold_Papalele_abund ~ trip.dat$new_fg, FUN = sum )
names(landings_onisland)<-c("location", "landings_sum_onisland")
landings_onisland<-merge(landings_onisland, landings_sumtot, by="location")
landings_onisland$landings_onisland_prop<-landings_onisland$landings_sum_onisland / landings_onisland$landings_sum_tot

landings_onisland_mean<-aggregate(trip.dat$landings_sold_personally_abund +
                                    trip.dat$landings_eaten_abund +
                                    trip.dat$landings_given_abund +
                                    trip.dat$landings_sold_Papalele_abund ~ trip.dat$new_fg, FUN = mean )
names(landings_onisland_mean)<-c("location", "landings_mean_onisland")
landings_onisland<-merge(landings_onisland, landings_onisland_mean, by="location")


landings_onisland<-subset(landings_onisland, select=-landings_sum_tot)



landings.dat<-merge(landings_meantot, landings_sumtot, by="location")
landings.dat<-merge(landings.dat, landings_personal, by="location")
landings.dat<-merge(landings.dat, landings_pengumpul, by="location")
landings.dat<-merge(landings.dat, landings_papalele, by="location")
landings.dat<-merge(landings.dat, landings_market, by="location")
landings.dat<-merge(landings.dat, landings_onisland, by="location")


###### Merge fish, oceanographic (MSEC), human pop data, rugosity, benthic cover, SST AND catch data using "site journal.xlsx" as site key: https://drive.google.com/open?id=1SNHtCmszbl6SYMPng1RLCDQVmap3e27n
drive_download(as_id("1SNHtCmszbl6SYMPng1RLCDQVmap3e27n"), overwrite=TRUE) # Saves file to working directory 
site.key<-read.csv("site journal-CLEANED-siteNames-removedsite17-decimalDegrees-meanVisibility.csv")
file.remove("site journal-CLEANED-siteNames-removedsite17-decimalDegrees-meanVisibility.csv")
site.key<-subset(site.key, select=c("site_id", "Site.Name", "lat_dd", "long_dd", "exposed", "u_visibility", "type_reef", "location"))

################################################################################################################
################################################################################################################




responseDF<-as.data.frame(cbind(fish.response=c("fish.logmass", "fish.mass", "fish.rich", "fish.shan", "fish.isim", "fish.even"),
                                fish.col=c("log_biomass_g", "biomass_g", "no_of_species", "shannon", "invsimpson", "SimpsonEvenness"),
                                fish.title=c("log Total Biomass (g)", "Total Biomass (g)", "Richness", "Shannon Diversity (H')", "Inverse Simpson's Diversity (D2)", "Simpson's Evenness (E)" )))


#### Next identify fish response object name, column name, and title
## CHOICES:
## for biomass, set as fish.mass 
## for log biomass, set as fish.logmass
## for species richness, set as fish.rich 
## for shannon diversity, set as fish.shan
## for inverse simpson's, set as fish.isim
## for simpson's evenness, set as fish.even
fish.response<-"fish.mass" # Set response here
responseRow<-grep(fish.response, responseDF$fish.response)
fish.title<-as.character(responseDF[responseRow, "fish.title"])
fish.col<-as.character(responseDF[responseRow, "fish.col"])

# Setting "rish.response", "fish.title", and "fish.col" above allows for the remainder of code below to be flexible based on desired response variable



# Merge all data: 
##### Do this in the following order: fish, fishing grounds (catch), oceanographic (MSEC), human pop data, rugosity, benthic cover, SST

# Merge fish data
dat.tmp<-merge(site.key, get(fish.response), by="site_id")

# Merge catch data
dat.tmp<-merge(dat.tmp, landings.dat, by="location", all.x=TRUE)
# note: all.x=TRUE because one UVC site (Furake on Hoga Island) was on a research station where there is zero fishing (no landings data)
dat.tmp[is.na(dat.tmp)]<-0 # Replace NAs for Furake site with 0

# Merge MSEC-SESYNC (oceanographic) data
###### FOR NOW, remove all human population data and other unnecessary columns
msec.datOnly<-subset(msec.dat, select=-c(no, long, lat, 
                                         npp_flag, 
                                         land_area_5km,
                                         wave_ww3_res, 
                                         pop1990_5km, pop2010_5km, pop2015_5km, pop2000_5km, dist_market))
dat.tmp<-merge(dat.tmp, msec.datOnly, by="Site.Name")

# Merge human pop data:
dat.tmp<-merge(dat.tmp, humanDensity.dat, by="Site.Name", all.x = TRUE)

# Replace NA human pop data with 0s
dat.tmp[is.na(dat.tmp)]<-0  

# Continue merging human data ## LEAVE THESE OUT FOR NOW
#dat.tmp<-merge(dat.tmp, distWeighted.dat, by="Site.Name")
#dat.tmp<-merge(dat.tmp, distToLandings.dat, by="Site.Name")

# Merge rugosity data
dat.tmp<-merge(dat.tmp, rug.site, by="Site.Name")



# Merge benthic cover data
# BENTHCOV.SITE is class=matrix; need to wrangle this into class=data.frame
benth.tmp<-as.data.frame(cbind(row.names(benthcov.site), benthcov.site))
names(benth.tmp)[1]<-"Site.Name"
benth.cols<-names(benth.tmp)[-1]
benth.tmp[benth.cols]<-apply(benth.tmp[benth.cols], MARGIN=2, FUN=as.character)
benth.tmp[benth.cols]<-apply(benth.tmp[benth.cols], MARGIN=2, FUN=as.numeric)
dat.tmp<-merge(dat.tmp, benth.tmp, by="Site.Name")

# Merge SST data (FINAL MERGE): 
sst.datOnly<-subset(sst.dat, select=c(site_id, SST_stdev, SST_50Perc, SST_98perc, SST_2perc, SST_kurtosis, SST_skewness))
alldat.site<-merge(dat.tmp, sst.datOnly, by="site_id")


### DIVIDE human metrics data by reef area
alldat.site$Population_2017<-alldat.site$Population_2017/alldat.site$reef_area_5km
alldat.site$No_of_Fishermen<-alldat.site$No_of_Fishermen/alldat.site$reef_area_5km
alldat.site$Row_Boats<-alldat.site$Row_Boats/alldat.site$reef_area_5km
alldat.site$Total_Motorboats<-alldat.site$Total_Motorboats/alldat.site$reef_area_5km
tmp.col<-grep("reef_area_5km", names(alldat.site))
alldat.site<-alldat.site[,-tmp.col]

setwd(outdir)
write.csv(alldat.site, "data_wakatobi_allDataMerged.csv", quote=FALSE, row.names = FALSE)

# QUESTION: how sensitive is SEM to (response) variable(s) having normal distribution?

# NEXT: create scatterplots

# Select all response + predictor variables + hierarchical variables
sitecols<-names(alldat.site) %in% c("site_id", "Site.Name", "lat_dd", "long_dd")
scatter.final<-alldat.site[,!sitecols]
loc.col<-grep("location", names(scatter.final))
response.col<-grep(fish.col, names(scatter.final))


scatter.names<-names(scatter.final)[-c(loc.col, response.col)]


########################################################################################
########################################################################################
########################################################################################
# Set graph names here: should match object scatter.names
scatter.titles<-c( # site journal columns
                  "Exposure", "Visibility", "Reef Type",  
                  
                  # Landings dat
                  "Landings per Trip", "Total Landings", 
                  "Total Personal", "Proportion Personal",
                  "Total Pengumpul", "Proportion Pengumpul", "Mean Pengumpul",
                  "Total Papalele", "Proportion Papalele",
                  "Total Market", "Proportion Market",
                  "Total On-Island", "Proportion On-Island", "Mean On-Island",
                  
                  # MSEC columns
                  "Mean NPP", "Min NPP", "Max NPP", "NPP SD", "Interannual NPP SD",
                  "Mean Wave Energy", "SD of Wave Energy", "Interannual SD of Wave Energy", "Wind Fetch", 
                  
                  # getWakatobiPop columns
                  "Population within 5km", "Fishers within 5km", "Rowboats within 5km", "Motorboats within 5km",
                  
                  # benthic columns
                  "Rugosity",
                  "Dead Coral With Algae", "Macroalgae", "Rubble", "Rock", "Soft Coral", "Sand", 
                  "All Hard Coral", "All Abiotic Substrate",
                  
                  # SST columns
                  "SST SD", "SST 50th Percentile", "SST 98th Percentile", "SST 2nd Percentile", "SST Kurtosis", "SST Skewness"
                  )

# Create scatterplots for fish.col
for(i in 1:length(scatter.names))
{
  newfile=paste("plot_scatter_", fish.col, "_vs_", scatter.names[i], ".pdf", sep="")
  p<-ggplot(data=scatter.final, aes(x=get(scatter.names[i]), y=get(fish.col))) + 
    geom_point(aes(x=get(scatter.names[i]), y=get(fish.col), color=location, shape=type_reef), size=2) +
    labs(y=fish.title, x=scatter.titles[i]) +  
    theme_light()+
    theme(text = element_text(size=18), 
          legend.position="right")
  
  pdf(file=newfile)
  print(p)      
  dev.off()
}

##### CALCULATE CORRELATIONS and VIF
##### NOTE: Correlations and VIF are unaffected by scaling (which preserves the shape and spread of data)
##### i.e., scaling can wait until the last step

##### Remove discreet variables
discreetvar<-match(c("wave_wind_fetch", "type_reef", "exposed"), scatter.names)
corr.names<-scatter.names[-discreetvar]
corr.final<-scatter.final[corr.names]


cor.dat<-cor(corr.final)
setwd(outdir)
write.csv(cor.dat, file="_Table_CorrelationsPearson.csv")
cor.test<-abs(cor.dat)>0.5
write.csv(cor.test, file="_Table_CorrelationsTestPearson.csv")

# Spearman's for all ranked, ordered vars (e.g., DACOR data)
cor.spear<-cor(corr.final, method="spearman")
write.csv(cor.spear, file="_Table_CorrelationsSpearman.csv")
cor.spear.test<-abs(cor.spear)>0.5
write.csv(cor.spear.test, file="_Table_CorrelationsTestSpearman.csv")

# Generate p values and confidence intervals for each correlation pair
pvals<-cor.mtest(corr.final, conf.level=0.95)

pdf(file="_Figure_CorrelationVisualPearson.pdf")
corrplot.mixed(cor.dat, upper="circle", lower="number", tl.pos="lt", tl.col="black", tl.cex=0.7, lower.col="black", addCoefasPercent=TRUE, number.cex=0.7, p.mat=pvals$p, sig.level=0.05, insig="blank", diag="n")
dev.off()

pdf(file="_Figure_CorrelationVisualSpearman.pdf")
corrplot.mixed(cor.spear, upper="circle", lower="number", tl.pos="lt", tl.col="black", tl.cex=0.7, lower.col="black", addCoefasPercent=TRUE, number.cex=0.7, p.mat=pvals$p, sig.level=0.05, insig="blank", diag="n")
dev.off()

##### QUESTION: do we have to worry about collinearity / multicollinearity in SEM?
# Supplementary table 6 here only tests correlation of exogenous variables (only those at the very start of a causation path): https://royalsocietypublishing.org/doi/suppl/10.1098/rspb.2016.0561#secSuppl
# Other studies, only address multicollinearity (but not sure which variables were tested): https://www.pnas.org/content/pnas/113/22/6230.full.pdf

##### ANSWER: as per Jonathan Lefcheck's email, both collinearity and multicollinearity are potential issues (but will only affect tests of significance, not parameter estimates)
## i.e., Both Rugosity and Hard Corals should not be used together as predictors of fish biomass

##### List of correlated variables:

## Suites of variables: Choose one of each
## All NPP variables with each other
## All wave variables with each other
## All human population variables with each other
## Some SST variables with each other: 
# SD with 98th and 2nd percentile and kurtosis
# 50th Percentile with 2nd percentile and skewness
# 98th Percentile with kurtosis and skewness
# 2nd percentile with kurtosis and skewness

## Other correlations:
## Visibility vs. 2nd percentile SST
## All wave variables vs. Rugosity, Soft Coral, Hard Coral, and some SST variables (SD, kurtosis, 2nd percentile)
## Rugosity vs. Hard Corals and SST_skew: 
## Soft corals vs Hard Corals and some SST variables (SD, 2nd and 98th percentile, and kurtosis)
## Hard coral vs. All abiotic benthos and SST_kurtosis
## Rubble with all human population variables: Don't use rubble

## Landings: TBD
# Currently, personal use vs. sold to market (papalele + pengumpul); Consider changing to: on-island (papalele + personal) vs off-island (pengumpul)

## STRATEGY: FIRST, Create path diagram, then trim based on correlations

### NOW, test for MULTICOLLINEARITY (calculate VIF)
### Seems like test for multicollinearity should be constructed for EACH linear model found within the SEM

#### Define PSEM equations here
#### Then call them later again for actual PSEM function

### REMINDER: for biomass, response.col indicates columns for raw and logged response variable
analysis.col<-grep(fish.col, names(alldat.site))

form1<-as.formula(paste(names(alldat.site)[analysis.col], " ~ All_HardCoral + landings_mean_onisland + Population_2017", sep=""))
form2<-as.formula("landings_mean_onisland ~ Population_2017")
form3<-as.formula("All_HardCoral ~ Population_2017")




#form3<-as.formula("Rugosity ~ All_HardCoral + wave_mean + Population_2017 ")  
#Had to remove rugosity because it was correlated with hard coral

fit1 <- lm(form1, data=alldat.site)
fit2 <- lm(form2, data=alldat.site)
fit3 <- lm(form3, data=alldat.site)


# Note if some path equations only contain two variables, VIF test below is invalid
vif.test1<-vif(fit1)
vif.test2<-vif(fit2)
vif.test3<-vif(fit3)


### SUBSET only model variables from alldat
vars1<-all.vars(form1)
vars2<-all.vars(form2)
vars3<-all.vars(form3)
model.vars<-unique(c(vars1, vars2, vars3))


# Aggregate to fishing ground level, scale and center data:
sem.vars.ground<-aggregate(sem.vars.site[model.vars], FUN=mean, by=list(sem.vars.site$location))
sem.ground.scaled<-as.data.frame(apply(sem.vars.ground[,-1], 2, scale))
sem.ground.scaled<-cbind(sem.vars.ground[,1], sem.ground.scaled)
names(sem.ground.scaled)[1]<- "location"

waka.psem<-psem(lm(form1, data=sem.ground.scaled), 
                lm(form2, data=sem.ground.scaled),
                lm(form3, data=sem.ground.scaled))


basisSet(waka.psem)
# NOTE: So far, there are no independence claims so basis set = 0

setwd(outdir)
txtname<-paste("stats_wakatobiSEM_", fish.col, ".txt", sep="")
sink(txtname)
print(summary(waka.psem, .progressBar = F))
sink()


# coefficients should already be standardized since data was already scaled
#coefs(waka.psem, standardize="scale")


#### RANDOM EFFECTS:
### Use site-level (unaggregated) data and incorporate random effects by location
# i.e., start with sem.dat.site

sem.site.scaled<-as.data.frame(apply(sem.vars.site[model.vars], 2, scale))
sem.site.scaled<-cbind(sem.vars.site$location, sem.site.scaled)
names(sem.site.scaled)[1]<- "location"

## For examples on hierarchical model specification see:
## http://www.rensenieuwenhuis.nl/r-sessions-21-multilevel-model-specification-nlme/

# Try to get simple models to converge before adding more complexity:
# Only random intercepts by location
wakarandom.psem<-psem(lme(form1, random = ~ 1 | location, data=sem.site.scaled), 
                      lme(form2, random = ~ 1 | location, data=sem.site.scaled),
                      lme(form3, random = ~ 1 | location, data=sem.site.scaled))


setwd(outdir)
txtname<-paste("stats_wakatobiSEM_", fish.col, "_randomIntercepts.txt", sep="")
sink(txtname)
print(summary(wakarandom.psem, .progressBar = F))
sink()

# For non-convergence problems, https://stats.stackexchange.com/questions/40647/lme-error-iteration-limit-reached
# use lmeControl to help with convergence?
#ctrl<-lmeControl(opt = "optim", maxIter)
#wakarandom.psem<-psem(lme(form1, random = ~ 1 | location, data=sem.site.scaled, method="ML", control=ctrl), 
#                      lme(form2, random = ~ 1 | location, data=sem.site.scaled, method="ML", control=ctrl),
#                      lme(form3, random = ~ 1 | location, data=sem.site.scaled, method="ML", control=ctrl))

## Does removing repeated measures help? NO
## sem.site.scaled2<-sem.site.scaled[!sem.site.scaled$location %in% c("hoga island", "tomia island", "wanci island", "komponaone island"),]
#wakarandom.psem<-psem(lme(form1, random = ~ 1 | location, data=sem.site.scaled2), 
#                      lme(form2, random = ~ 1 | location, data=sem.site.scaled2),
#                      lme(form3, random = ~ 1 | location, data=sem.site.scaled2))


## Random intercepts for site-level + 
## predictor that is allowed to vary between groups (e.g., effect of coral) + 
## group-level predictor (e.g., effect of total landings)
wakarandom.psem<-psem(lme(form1, random = ~ 1 + All_HardCoral | location, data=sem.site.scaled, method="ML"), 
                      #lme(form2, random = ~ 1 + Population_2017 | location, data=sem.site.scaled),
                      lme(form3, random = ~ 1 + All_HardCoral | location, data=sem.site.scaled, method="ML"))
     
wakarandom.psem<-psem(lme(form1, random = ~ 1 + All_HardCoral | location, data=sem.site.scaled, method="ML", control=ctrl), 
                      #lme(form2, random = ~ 1 + Population_2017 | location, data=sem.site.scaled),
                      lme(form3, random = ~ 1 + All_HardCoral | location, data=sem.site.scaled, method="ML", control=ctrl))


setwd(outdir)
txtname<-paste("stats_wakatobiSEM_", fish.col, "_randomInterceptsAndSlopes.txt", sep="")
sink(txtname)
print(summary(wakarandom.psem, .progressBar = F))
sink()

