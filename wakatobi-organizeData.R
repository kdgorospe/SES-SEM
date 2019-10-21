# Requires permission (shared Google drive folder) and access (internet) to Wakatobi Data
rm(list=ls())

#update.packages(ask = FALSE, checkBuilt = TRUE)
library(googledrive)
library(tidyverse)
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
fish.mass<-fishdat %>%
  group_by(site_id, transect) %>%
  summarise(transect.mass=sum(biomass_g, na.rm=TRUE)) %>%
  summarise(biomass_g=mean(transect.mass))

# Old code using base R functions:
#fish.mass1<-aggregate(biomass_g ~ site_id + transect, data=fishdat, FUN=sum)
#ish.mass<-aggregate(biomass_g ~ site_id, data=fish.mass1, FUN=mean)

# Calculate species diversity and test all of them
# See Morris et al. (Ecology and Evolution) for discussion simultaneously considering analyses
# of multiple indices can provide greater insight

# 1 - richness 
fish.rich<-fishdat %>% 
  group_by(site_id, transect) %>%
  summarise(transect.count=n())  %>%
  summarise(no_of_species=mean(transect.count))

### LEFT OFF HERE: Current data frame has missing trophic group data
# 2- shannon diversity (aka H') - using "TROPHIC_GROUP" or FUNCTIONAL DIVERSITY
fishdat %>%
  filter(is.na(functional_group)) %>%
  distinct(scientific_name)

# The missing data is counted explicitly when using tidyverse functions:
fishdat %>%
  group_by(site_id, transect, trophic_group) %>%
  summarise(count.transect=sum(number_of_fish, na.rm=TRUE)) %>%
  group_by(site_id, trophic_group)
  

### OLD CODE using base R functions
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

# merge all fish data together
fish.dat<-merge(fish.mass, fish.logmass, by="site_id")
fish.dat<-merge(fish.dat, fish.rich, by="site_id")
fish.dat<-merge(fish.dat, fish.shan, by="site_id")
fish.dat<-merge(fish.dat, fish.isim, by="site_id")
fish.dat<-merge(fish.dat, fish.even, by="site_id")

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


# input 5km human population data: https://drive.google.com/open?id=1DcVqeVEx6yGksBqLGzbWhU8WxN6UcYBz

#drive_download(as_id("1DcVqeVEx6yGksBqLGzbWhU8WxN6UcYBz"), overwrite=TRUE) # Saves file to working directory 
#humanDensity.dat<-read.csv("data_wakatobiHumans_areaWeightedDensityMetrics_5_km_buffer.csv") # weights each village's population density by its area to get "total population" within 5km buffer
#file.remove("data_wakatobiHumans_areaWeightedDensityMetrics_5_km_buffer.csv")

# input 2.5km human population data: https://drive.google.com/open?id=1xf91oaXfqp-BDKIPcgK0bW-AfvJ3joHY
drive_download(as_id("1xf91oaXfqp-BDKIPcgK0bW-AfvJ3joHY"), overwrite=TRUE) # Saves file to working directory 
humanDensity.dat<-read.csv("data_wakatobiHumans_areaWeightedDensityMetrics_2.5_km_buffer.csv") # weights each village's population density by its area to get "total population" within 5km buffer
file.remove("data_wakatobiHumans_areaWeightedDensityMetrics_2.5_km_buffer.csv")

# input 10 km human population data: https://drive.google.com/open?id=1c05yh4thZTqhEPMA0SrdRCuM8n86b1G1
drive_download(as_id("1c05yh4thZTqhEPMA0SrdRCuM8n86b1G1"), overwrite=TRUE) # Saves file to working directory 
humanDensity.dat<-read.csv("data_wakatobiHumans_areaWeightedDensityMetrics_10_km_buffer.csv") # weights each village's population density by its area to get "total population" within 5km buffer
file.remove("data_wakatobiHumans_areaWeightedDensityMetrics_10_km_buffer.csv")


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

landings_personal_mean<-aggregate(trip.dat$landings_sold_personally_abund +
                               trip.dat$landings_eaten_abund +
                               trip.dat$landings_given_abund ~ trip.dat$new_fg, FUN = mean )
names(landings_personal_mean)<-c("location", "landings_mean_personal")
landings_personal<-merge(landings_personal, landings_personal_mean, by="location")


# How do landings sold to pengumpul (GO OFF ISLAND) affect fish response? - include sum, proportion, and mean off-island catch
landings_pengumpul<-aggregate(trip.dat$landings_sold_Pengumpul_abund ~ trip.dat$new_fg, FUN = sum )
names(landings_pengumpul)<-c("location", "landings_sum_pengumpul")
landings_pengumpul<-merge(landings_pengumpul, landings_sumtot, by="location")
landings_pengumpul$landings_pengumpul_prop<-landings_pengumpul$landings_sum_pengumpul / landings_pengumpul$landings_sum_tot
landings_pengumpul<-subset(landings_pengumpul, select=-landings_sum_tot)

landings_pengumpul_mean<-aggregate(trip.dat$landings_sold_Pengumpul_abund ~ trip.dat$new_fg, FUN = mean )
names(landings_pengumpul_mean)<-c("location", "landings_mean_pengumpul")
landings_pengumpul<-merge(landings_pengumpul, landings_pengumpul_mean, by="location")


# How do landings sold to papalele affect fish response? # Papalele = sold ON-ISLAND
landings_papalele<-aggregate(trip.dat$landings_sold_Papalele_abund ~ trip.dat$new_fg, FUN = sum )
names(landings_papalele)<-c("location", "landings_sum_papalele")
landings_papalele<-merge(landings_papalele, landings_sumtot, by="location")
landings_papalele$landings_papalele_prop<-landings_papalele$landings_sum_papalele / landings_papalele$landings_sum_tot
landings_papalele<-subset(landings_papalele, select=-landings_sum_tot)

landings_papalele_mean<-aggregate(trip.dat$landings_sold_Papalele_abund ~ trip.dat$new_fg, FUN = mean )
names(landings_papalele_mean)<-c("location", "landings_mean_papalele")
landings_papalele<-merge(landings_papalele, landings_papalele_mean, by="location")


# How do landings sold to ANY fish market (i.e., pengumpul + papalele) affect fish response?
landings_market<-aggregate(trip.dat$landings_sold_Papalele_abund +
                             trip.dat$landings_sold_Pengumpul_abund ~ trip.dat$new_fg, FUN=sum)
names(landings_market)<-c("location", "landings_sum_market")
landings_market<-merge(landings_market, landings_sumtot, by="location")
landings_market$landings_market_prop<-landings_market$landings_sum_market / landings_market$landings_sum_tot 
landings_market<-subset(landings_market, select=-landings_sum_tot)

landings_market_mean<-aggregate(trip.dat$landings_sold_Papalele_abund +
                             trip.dat$landings_sold_Pengumpul_abund ~ trip.dat$new_fg, FUN=mean)
names(landings_market_mean)<-c("location", "landings_mean_market")
landings_market<-merge(landings_market, landings_market_mean, by="location")


# How do landings kept ON ISLAND (i.e., sold to papalele, sold personally, eaten, or given) affect fish response? - include sum, proportion, and mean off-island catch
landings_onisland<-aggregate(trip.dat$landings_sold_personally_abund +
                               trip.dat$landings_eaten_abund +
                               trip.dat$landings_given_abund +
                               trip.dat$landings_sold_Papalele_abund ~ trip.dat$new_fg, FUN = sum )
names(landings_onisland)<-c("location", "landings_sum_onisland")
landings_onisland<-merge(landings_onisland, landings_sumtot, by="location")
landings_onisland$landings_onisland_prop<-landings_onisland$landings_sum_onisland / landings_onisland$landings_sum_tot
landings_onisland<-subset(landings_onisland, select=-landings_sum_tot)

landings_onisland_mean<-aggregate(trip.dat$landings_sold_personally_abund +
                                    trip.dat$landings_eaten_abund +
                                    trip.dat$landings_given_abund +
                                    trip.dat$landings_sold_Papalele_abund ~ trip.dat$new_fg, FUN = mean )
names(landings_onisland_mean)<-c("location", "landings_mean_onisland")
landings_onisland<-merge(landings_onisland, landings_onisland_mean, by="location")


landings.dat<-merge(landings_meantot, landings_sumtot, by="location")
landings.dat<-merge(landings.dat, landings_personal, by="location")
landings.dat<-merge(landings.dat, landings_pengumpul, by="location")
landings.dat<-merge(landings.dat, landings_papalele, by="location")
landings.dat<-merge(landings.dat, landings_market, by="location")
landings.dat<-merge(landings.dat, landings_onisland, by="location")

## Rearrange so that "sums", "means", and "proportions" are together:
location_cols<-grep("location", names(landings.dat))
sum_cols<-grep("sum", names(landings.dat))
mean_cols<-grep("mean", names(landings.dat))
prop_cols<-grep("prop", names(landings.dat))

landings.dat<-landings.dat[c(location_cols, mean_cols, sum_cols, prop_cols)]



###### Merge fish, oceanographic (MSEC), human pop data, rugosity, benthic cover, SST AND catch data using "site journal.xlsx" as site key: https://drive.google.com/open?id=1SNHtCmszbl6SYMPng1RLCDQVmap3e27n
drive_download(as_id("1SNHtCmszbl6SYMPng1RLCDQVmap3e27n"), overwrite=TRUE) # Saves file to working directory 
site.key<-read.csv("site journal-CLEANED-siteNames-removedsite17-decimalDegrees-meanVisibility.csv")
file.remove("site journal-CLEANED-siteNames-removedsite17-decimalDegrees-meanVisibility.csv")
site.key<-subset(site.key, select=c("site_id", "Site.Name", "lat_dd", "long_dd", "exposed", "u_visibility", "type_reef", "location"))


# Merge all data: 
##### Do this in the following order: fish, fishing grounds (catch), oceanographic (MSEC), human pop data, rugosity, benthic cover, SST

# Merge fish data
dat.tmp<-merge(site.key, fish.dat, by="site_id")

# Merge catch data
dat.tmp<-merge(dat.tmp, landings.dat, by="location", all.x=TRUE)
# note: all.x=TRUE because one UVC site (Furake on Hoga Island) was on a research station where there is zero fishing (no landings data)
dat.tmp[is.na(dat.tmp)]<-0 # Replace NAs for Furake site with 0

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

# Merge MSEC-SESYNC (oceanographic) data
###### FOR NOW, remove all human population data and other unnecessary columns
msec.datOnly<-subset(msec.dat, select=-c(no, long, lat, 
                                         npp_flag, 
                                         land_area_5km,
                                         wave_ww3_res, 
                                         pop1990_5km, pop2010_5km, pop2015_5km, pop2000_5km, dist_market))
dat.tmp<-merge(dat.tmp, msec.datOnly, by="Site.Name")


# Merge SST data (FINAL MERGE): 
sst.datOnly<-subset(sst.dat, select=c(site_id, SST_stdev, SST_50Perc, SST_98perc, SST_2perc, SST_kurtosis, SST_skewness))
dat.tmp<-merge(dat.tmp, sst.datOnly, by="site_id")

# Merge human pop data:
alldat.site<-merge(dat.tmp, humanDensity.dat, by="Site.Name", all.x = TRUE)
# Replace NA human pop data with 0s
alldat.site[is.na(alldat.site)]<-0  


### DIVIDE human metrics data by reef area
#alldat.site$Population_2017<-alldat.site$Population_2017/alldat.site$reef_area_5km
#alldat.site$No_of_Fishermen<-alldat.site$No_of_Fishermen/alldat.site$reef_area_5km
#alldat.site$Row_Boats<-alldat.site$Row_Boats/alldat.site$reef_area_5km
#alldat.site$Total_Motorboats<-alldat.site$Total_Motorboats/alldat.site$reef_area_5km
#tmp.col<-grep("reef_area_5km", names(alldat.site))
#alldat.site<-alldat.site[,-tmp.col]

setwd(outdir)
write.csv(alldat.site, "data_wakatobi_allDataMerged.csv", quote=FALSE, row.names = FALSE)
