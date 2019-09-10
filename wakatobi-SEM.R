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


# Calculate total fish biomass at site level (site_id column), averaged across three transects
fish.mass1<-aggregate(biomass_g ~ site_id + transect, data=fishdat, FUN=sum)
fish.mass<-aggregate(biomass_g ~ site_id, data=fish.mass1, FUN=mean)


# Calculate fish biomass by functional group at site level, averaged across three transects
fun.mass1<-aggregate(biomass_g ~ site_id + transect + trophic_group, data=fishdat, FUN=sum)
fun.mass<-aggregate(biomass_g ~ site_id + trophic_group, data=fun.mass1, FUN=mean)


# Calculate species diversity and test all of them
# See Morris et al. (Ecology and Evolution) for discussion simultaneously considering analyses
# of multiple indices can provide greater insight

# 1 - richness 
spcount.transect<-aggregate(scientific_name ~ site_id + transect, data=fishdat, FUN=length)
names(spcount.transect)[3]<-"no_of_species" 
# Note: no_of_species is a bit of a misnomer; most likely an undercount. EXAMPLE: "Acanthurs spp" could be used for more than one species, but only counted as one species here
fish.rich<-aggregate(no_of_species ~ site_id, data=spcount.transect, FUN=mean)

# 2- shannon diversity (aka H')
### need to re-check below for diversity calculation
countPerSp.transect<-aggregate(number_of_fish ~ site_id + scientific_name + transect, data=fishdat, FUN=sum)
countPerSp.site<-aggregate(number_of_fish ~ site_id + scientific_name, data=countPerSp.transect, FUN=mean)
# Convert to matrix for calculating diversity:
countPerSp.mat<-acast(countPerSp.site, site_id~scientific_name, value.var = "number_of_fish", FUN=sum)
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

### OTHER POTENTIAL RESPONSE VARIABLES: functional trait diversity


# For now, focus on total fish biomass for analysis: 

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

# input LANDINGS TRIP data: https://drive.google.com/open?id=10C0WT0Psz00bHkw8S-NGbPU9FR7Ufgzv
# Path file: /Users/KGthatsme/Projects/Google Drive/Wakatobi-SEMAnalysis/_landingsData/Wakatobi-landings_20190910_TRIP-cleanedFishingGrounds.csv
drive_download(as_id("10C0WT0Psz00bHkw8S-NGbPU9FR7Ufgzv"), overwrite=TRUE)
trip.dat<-read.csv("Wakatobi-landings_20190910_TRIP_cleanedFishingGrounds.csv")
file.remove("Wakatobi-landings_20190910_TRIP_cleanedFishingGrounds.csv")

# Select relevant columns:
trip.dat<-subset(trip.dat, select=c(trip_id, fishing_grnd1, landing_no, landing_unit,
                          landings_sold_personally_no, landings_sold_personally_unit,
                          landings_sold_Papalele_no, landings_sold_Papalele_unit, 
                          landings_sold_Pengumpul_no, landings_sold_Pengumpul_unit,
                          landings_eaten_no, landings_eaten_unit,
                          landings_given_no, landings_given_unit))

summary(trip.dat)
# Need to clean column: "landings_sold_personally_no"; inputted as "factor" and should be "numeric"

levels(trip.dat$landings_sold_personally_no)
trip.dat$landings_sold_personally_no<-as.character(trip.dat$landings_sold_personally_no)

fix_index<-grep("1 box kecil", trip.dat$landings_sold_personally_no)
trip.dat$landings_sold_personally_no[fix_index]<-"1"
trip.dat$landings_sold_personally_unit[fix_index]<-"box kecil"


fix_index<-grep("21 ekor jual sendiri", trip.dat$landings_sold_personally_no) # selling 21 fish personally
trip.dat$landings_sold_personally_no[fix_index]<-"21"
trip.dat$landings_sold_personally_unit[fix_index]<-"fish"

fix_index<-grep("jual ke pasar dan pengumpul", trip.dat$landings_sold_personally_no) # selling to market and collectors - ie, need to move entire landings of "1 box" to "Pengumpul" column
trip.dat$landings_sold_Pengumpul_no[fix_index]<-1
trip.dat$landings_sold_Pengumpul_unit[fix_index]<-"box"
trip.dat$landings_sold_personally_no[fix_index]<-"0"

trip.dat$landings_sold_personally_no<-as.numeric(trip.dat$landings_sold_personally_no)

# Convert all NAs in "no." columns to 0s
no_cols<-grep("_no", names(trip.dat))
for (i in no_cols)
{
  trip.dat[,i][is.na(trip.dat[,i])]<-0
}


# Standardize all landing_units: 
levels(trip.dat$landing_unit)
# First, set as character so that new factors can be added (e.g., can't add small bucket since it's not currently a factor level)
trip.dat$landing_unit<-as.character(trip.dat$landing_unit)
trip.dat$landing_unit[grep("box kecil", trip.dat$landing_unit)]<-"small box"
trip.dat$landing_unit[grep("smal box", trip.dat$landing_unit)]<-"small box"
trip.dat$landing_unit[grep("bucket kecil", trip.dat$landing_unit)]<-"small bucket"
trip.dat$landing_unit[grep("ekor", trip.dat$landing_unit)]<-"fish"
trip.dat$landing_unit[grep("fish ", trip.dat$landing_unit)]<-"fish"
# Final check
table(trip.dat$landing_unit)
# Now, reset as factor
trip.dat$landing_unit<-as.factor(trip.dat$landing_unit)
levels(trip.dat$landing_unit)


# Standardize all landings_sold_personally_units: 
levels(trip.dat$landings_sold_personally_unit)
trip.dat$landings_sold_personally_unit<-as.character(trip.dat$landings_sold_personally_unit)
trip.dat$landings_sold_personally_unit[grep("box kecil", trip.dat$landings_sold_personally_unit)]<-"small box"
trip.dat$landings_sold_personally_unit[grep("smal box", trip.dat$landings_sold_personally_unit)]<-"small box"
trip.dat$landings_sold_personally_unit[grep("ekor", trip.dat$landings_sold_personally_unit)]<-"fish"
trip.dat$landings_sold_personally_unit[!(trip.dat$landings_sold_personally_unit %in% c("basket", "box", "bucket", "small box", "small bucket", "fish"))]<-NA
# Final check
table(trip.dat$landings_sold_personally_unit)
# Now, reset as factor
trip.dat$landings_sold_personally_unit<-as.factor(trip.dat$landings_sold_personally_unit)
levels(trip.dat$landings_sold_personally_unit)


# Standardize all landings_sold_Papalele_units: 
levels(trip.dat$landings_sold_Papalele_unit)
trip.dat$landings_sold_Papalele_unit<-as.character(trip.dat$landings_sold_Papalele_unit)
trip.dat$landings_sold_Papalele_unit[grep("box kecil", trip.dat$landings_sold_Papalele_unit)]<-"small box"
trip.dat$landings_sold_Papalele_unit[grep("live fish", trip.dat$landings_sold_Papalele_unit)]<-"fish"
trip.dat$landings_sold_Papalele_unit[!(trip.dat$landings_sold_Papalele_unit %in% c("basket", "box", "bucket", "small box", "small bucket", "fish"))]<-NA
# Final check
table(trip.dat$landings_sold_Papalele_unit)
# Now, reset as factor
trip.dat$landings_sold_Papalele_unit<-as.factor(trip.dat$landings_sold_Papalele_unit)
levels(trip.dat$landings_sold_Papalele_unit)



# Standardize all landings_sold_Pengumpul_units: 
levels(trip.dat$landings_sold_Pengumpul_unit)
trip.dat$landings_sold_Pengumpul_unit<-as.character(trip.dat$landings_sold_Pengumpul_unit)
trip.dat$landings_sold_Pengumpul_unit[!(trip.dat$landings_sold_Pengumpul_unit %in% c("basket", "box", "bucket", "small box", "small bucket", "fish"))]<-NA
# Final check
table(trip.dat$landings_sold_Pengumpul_unit)
# Now, reset as factor
trip.dat$landings_sold_Pengumpul_unit<-as.factor(trip.dat$landings_sold_Pengumpul_unit)
levels(trip.dat$landings_sold_Pengumpul_unit)




# Standardize all landings_eaten_units: 
levels(trip.dat$landings_eaten_unit)
trip.dat$landings_eaten_unit<-as.character(trip.dat$landings_eaten_unit)
trip.dat$landings_eaten_unit[!(trip.dat$landings_eaten_unit %in% c("basket", "box", "bucket", "small box", "small bucket", "fish"))]<-NA
# Final check
table(trip.dat$landings_eaten_unit)
# Now, reset as factor
trip.dat$landings_eaten_unit<-as.factor(trip.dat$landings_eaten_unit)
levels(trip.dat$landings_eaten_unit)


# Standardize all landings_given_units: 
levels(trip.dat$landings_given_unit)
trip.dat$landings_given_unit<-as.character(trip.dat$landings_given_unit)
trip.dat$landings_given_unit[!(trip.dat$landings_given_unit %in% c("basket", "box", "bucket", "small box", "small bucket", "fish"))]<-NA
# Final check
table(trip.dat$landings_given_unit)
# Now, reset as factor
trip.dat$landings_given_unit<-as.factor(trip.dat$landings_given_unit)
levels(trip.dat$landings_given_unit)


# Conversion of units of volume to abundance (as per Melati)
# essentially converting all units to "fish units"
box          <- 58
basket       <- 16
small_box    <- 19
bucket       <- 14
small_bucket <- 7

unit_cols<-grep("unit", names(trip.dat))
for (i in unit_cols)
{
  # make new column for fish abundance units
  trip.dat$newcol<-0
  trip.dat$newcol[grep("box", trip.dat[,i])]<-box
  trip.dat$newcol[grep("basket", trip.dat[,i])]<-basket
  trip.dat$newcol[grep("small box", trip.dat[,i])]<-small_box
  trip.dat$newcol[grep("bucket", trip.dat[,i])]<-bucket
  trip.dat$newcol[grep("small bucket", trip.dat[,i])]<-small_bucket
  
  # for data already in units of fish abundance (i.e., unit column == "fish"), insert into "newcol" whatever number is in column "landing_no"
  trip.dat$newcol[grep("fish", trip.dat[,i])]<-trip.dat$landing_no[grep("fish", trip.dat[,i])]
  
  # rename newcol
  renamecol<-grep("newcol", names(trip.dat))
  names(trip.dat)[renamecol]<-paste(names(trip.dat)[i], "_abund", sep="")
}

# Calculate all fish flow abundances by multiplying "no." column with "unit_abund" column
trip.dat$landings_abund<-trip.dat$landing_no * trip.dat$landing_unit_abund
trip.dat$landings_sold_personally_abund<-trip.dat$landings_sold_personally_no * trip.dat$landings_sold_personally_unit_abund
trip.dat$landings_sold_Papalele_abund<-trip.dat$landings_sold_Papalele_no * trip.dat$landings_sold_Papalele_unit_abund
trip.dat$landings_sold_Pengumpul_abund<-trip.dat$landings_sold_Pengumpul_no * trip.dat$landings_sold_Pengumpul_unit_abund
trip.dat$landings_eaten_abund<-trip.dat$landings_eaten_no * trip.dat$landings_eaten_unit_abund
trip.dat$landings_given_abund<-trip.dat$landings_given_no * trip.dat$landings_given_unit_abund

### CHECK that landings_abund = sum(all other landings_abund)
fishflowQC<-trip.dat$landings_abund != trip.dat$landings_sold_personally_abund + 
  trip.dat$landings_sold_Papalele_abund + 
  trip.dat$landings_sold_Pengumpul_abund + 
  trip.dat$landings_eaten_abund + 
  trip.dat$landings_given_abund

### LEFT OFF HERE - use rownames output below to figure out which trips need to be re-entered from raw data
trip.dat_needsQC<-trip.dat[fishflowQC,]
rownames(trip.dat_needsQC)

## input aggregation file for landings trips: https://drive.google.com/open?id=1PkaXlA1r1RA6tUWX7Tm3sk3SPm7kJxMf
drive_download(as_id("1PkaXlA1r1RA6tUWX7Tm3sk3SPm7kJxMf"), overwrite=TRUE)
trip.agg<-read.csv("aggregationKey-FishingGround.csv")
file.remove("aggregationKey-FishingGround.csv")


# Aggregate (calculate mean) of groups of fishing grounds based on column: new_fg
trip.dat$landings_abund




# NEXT: Merge fish, oceanographic (MSEC), human pop data, rugosity, benthic cover, SST using "site journal.xlsx" as site key: https://drive.google.com/open?id=1SNHtCmszbl6SYMPng1RLCDQVmap3e27n
drive_download(as_id("1SNHtCmszbl6SYMPng1RLCDQVmap3e27n"), overwrite=TRUE) # Saves file to working directory 
site.key<-read.csv("site journal-CLEANED-siteNames-removedsite17-decimalDegrees-meanVisibility.csv")
file.remove("site journal-CLEANED-siteNames-removedsite17-decimalDegrees-meanVisibility.csv")
site.key<-subset(site.key, select=c("site_id", "Site.Name", "lat_dd", "long_dd", "exposed", "u_visibility", "type_reef", "location"))

# Merge all data: 
##### Do this in the following order: fish, oceanographic (MSEC), human pop data, rugosity, benthic cover, SST

################################################################################################################
################################################################################################################


responseDF<-as.data.frame(cbind(fish.response=c("fish.mass", "fish.rich", "fish.shan", "fish.isim", "fish.even"),
                                fish.col=c("biomass_g", "no_of_speces", "shannon", "invsimpson", "SimpsonEvenness"),
                                fish.title=c("Total Biomass (g)", "Richness", "Shannon Diversity (H')", "Inverse Simpson's Diversity (D2)", "Simpson's Evenness (E)" )))


#### First, identify fish response object name, column name, and title
## CHOICES:
## for biomass, set as fish.mass 
## for species richness, set as fish.rich 
## for shannon diversity, set as fish.shan
## for inverse simpson's, set as fish.isim
## for simpson's evenness, set as fish.even
fish.response<-"fish.mass" # Set response here
responseRow<-grep(fish.response, responseDF$fish.response)
fish.title<-as.character(responseDF[responseRow, "fish.title"])
fish.col<-as.character(responseDF[responseRow, "fish.col"])

# Setting "rish.response", "fish.title", and "fish.col" above allows for the remainder of code below to be flexible based on desired response variable



# Merge fish data
dat.tmp<-merge(site.key, get(fish.response), by="site_id")

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

# NEXT: plot histograms, consider log transforming some variables
setwd(outdir)
pdf(file="plot_histogram.totalbiomass.pdf")
hist(alldat.site[,fish.col],xlab=fish.title, main=paste("Histogram of Site-Level ", fish.title, sep=""))
dev.off()

# Try log biomass
pdf(file="plot_histogram.LOGtotalbiomass.pdf")
hist(log10(alldat.site[,fish.col]),xlab=paste("log ", fish.title, sep=""), main=paste("Histogram of Site-Level log ", fish.title, sep=""))
dev.off()

# Include logbiomass as a potential response variable in data.frame
alldat.site$log_biomass_g<-log10(alldat.site[,fish.col])

# QUESTION: how sensitive is SEM to (response) variable(s) having normal distribution?

# NEXT: create scatterplots

# Select all response + predictor variables + hierarchical variables
scatter.final<-alldat.site[,-c(1:4)]
loc.col<-grep("location", names(scatter.final))
bio.col<-grep(fish.col, names(scatter.final))

#### LEFT OFF HERE - replace bio.col with fish.col
scatter.names<-names(scatter.final)[-c(loc.col, bio.col)]


########################################################################################
########################################################################################
########################################################################################
# Set graph names here: should match object scatter.names
scatter.titles<-c( # site journal columns
                  "Exposure", "Visibility", "Reef Type",  
                  
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

# Create scatterplots for raw total biomass
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

# ReDO scatterplots for logbiomass:

for(i in 1:length(scatter.names))
{
  newfile=paste("plot_scatter_LOG_", fish.col, "_vs_", scatter.names[i], ".pdf", sep="")
  p<-ggplot(data=scatter.final, aes(x=get(scatter.names[i]), y=get(paste("log_", fish.col, sep="")))) + 
    geom_point(aes(x=get(scatter.names[i]), y=get(paste("log_", fish.col, sep="")), color=location, shape=type_reef), size=2) +
    labs(y=paste("log", fish.title, sep=""), x=scatter.titles[i]) +  
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


## STRATEGY: FIRST, Create path diagram, then trim based on correlations

## Notes: Remove rugosity (correlated with hard corals)

### NOW, test for MULTICOLLINEARITY (calculate VIF)
### Seems like test for multicollinearity should be constructed for EACH linear model found within the SEM

#### Define PSEM equations here
#### Then call them later again for actual PSEM function

#### LEFT OFF HERE

form1<-as.formula("log_biomass_g ~ Population_2017 + All_HardCoral + MA + wave_interann_sd + SST_98perc + npp_mean")
form2<-as.formula("All_HardCoral ~ Population_2017 + wave_interann_sd + SST_98perc + npp_mean")
form3<-as.formula("MA ~ Population_2017 + SST_98perc")

#form3<-as.formula("Rugosity ~ All_HardCoral + wave_mean + Population_2017 ")  
#Had to remove rugosity because it was correlated with hard coral

fit1 <- lm(form1, data=alldat.site)
fit2 <- lm(form2, data=alldat.site)
fit3 <- lm(form3, data=alldat.site)

vif.test1<-vif(fit1)
vif.test2<-vif(fit2)
vif.test3<-vif(fit3)


### MODIFY AND REPEAT:
### Remove waves (tried mean, sd, and interann_sd and ALL fail multicollinearity test)
form1<-as.formula("log_biomass_g ~ Population_2017 + All_HardCoral + MA + SST_98perc + npp_mean")
form2<-as.formula("All_HardCoral ~ Population_2017 + SST_98perc + npp_mean")
form3<-as.formula("MA ~ Population_2017 + SST_98perc")


fit1 <- lm(form1, data=alldat.site)
fit2 <- lm(form2, data=alldat.site)
fit3 <- lm(form3, data=alldat.site)

vif.test1<-vif(fit1)
vif.test2<-vif(fit2)
vif.test3<-vif(fit3)


### TWO POTENTIAL HIEARCHIES for analysis: location or reef_type (i.e., atoll vs finging reef)
### Modify "location" to be one of 6 islands/atolls: Wanci, Kaledupa, Tomia, Binongko, Kaledupa Atoll, Kapota Atoll
### See: map_wakatobi_FishSiteNames.pdf for REFERENCE
### Use this info for hierarchical analysis below:
alldat.site$location<-as.character(alldat.site$location)
alldat.site$location[alldat.site$location=="hoga"]<-"kaledupa"
alldat.site$location[alldat.site$location=="kapota island"]<-"wanci-wanci"
alldat.site$location[alldat.site$location=="wanci island"]<-"wanci-wanci"
alldat.site$location[alldat.site$location=="komponaone island"]<-"wanci-wanci"
alldat.site$location<-as.factor(alldat.site$location)

# Scale and center data:
sem.vars<-c(corr.names, "biomass_g", "log_biomass_g")
#sem.vars<-c( "biomass_g", "log_biomass_g", 
#             "u_visibility", "Rugosity", "All_HardCoral", "All_Abiotic", "MA",
#             "Population_2017", "No_of_Fishermen", "Row_Boats", "Total_Motorboats", 
#             "npp_mean", "npp_sd", "npp_interann_sd",
#             "wave_mean", "wave_sd", "wave_interann_sd", 
#             "SST_98perc")

# See SiteRankOrder-forVariousMSECOceanographMetrics.csv: 
# Choice of NPP makes a big difference in rank-order of sites
# Choice of wave is less important

sem.dat<-alldat.site[sem.vars]
sem.dat.scaled<-as.data.frame(apply(sem.dat, 2, scale))
sem.dat.scaled<-cbind(alldat.site$Site.Name, alldat.site$location, sem.dat.scaled)
names(sem.dat.scaled)[1:2]<-c("Site.Name", "location")






waka.psem<-psem(lm(form1, data=sem.dat.scaled), 
                lm(form2, data=sem.dat.scaled),
                lm(form3, data=sem.dat.scaled))
  
  
# FIT piecewise SEM (based on DAG-EcoOnlySEM.JPG): no random effects - create list of structured equations
#waka.psem<-psem(lm(log_biomass_g ~ Row_Boats + All_HardCoral + Rugosity, data = sem.dat.scaled),
#                lm(All_HardCoral ~ Row_Boats, data = sem.dat.scaled),
#                lm(Rugosity ~ All_HardCoral + Row_Boats, data = sem.dat.scaled)
#)


basisSet(waka.psem)
# NOTE: So far, there are no independence claims so basis set = 0

setwd(outdir)
txtname<-paste("stats_wakatobiSEM_ecology_and_Population2017.txt", sep="")
sink(txtname)
print(summary(waka.psem, .progressBar = F))
sink()

# coefficients should already be standardized since data was already scaled
#coefs(waka.psem, standardize="scale")


### Incorporate random effects by location
# With island-level random effects, here random intercepts are modeled for each island 


form_a<-as.formula(paste("log_biomass_g ~ ", as.symbol(humanmetrics[i]), "+ All_HardCoral + Rugosity"))
form_b<-as.formula(paste("All_HardCoral ~ ", as.symbol(humanmetrics[i])))  
form_c<-as.formula(paste("Rugosity ~ All_HardCoral + ", as.symbol(humanmetrics[i])))  

wakarandom.psem<-psem(lme(form_a, random = ~ 1 | location, data=sem.dat.scaled), 
                lme(form_b, random = ~ 1 | location, data=sem.dat.scaled),
                lme(form_c, random = ~ 1 | location, data=sem.dat.scaled))


#wakarandom.psem<-psem(lme(log_biomass_g ~ Population_2017 + All_HardCoral + Rugosity, random = ~ 1 | location, data = sem.dat.scaled),
#                lme(All_HardCoral ~ Population_2017, random = ~ 1 | location, data = sem.dat.scaled),
#                lme(Rugosity ~ All_HardCoral + Population_2017, random = ~ 1 | location, data = sem.dat.scaled)
#)


setwd(outdir)
txtname<-paste("stats_wakatobiSEM_ecologyOnly-randomIslandIntercepts_", humanmetrics[i], ".txt", sep="")
sink(txtname)
print(summary(wakarandom.psem, .progressBar = F))
sink()

### New Code: Still need to try: Random slopes and interecepts for all predictors
#wakarand_slope_intercept.psem<-psem(lme(form_a, random = ~ 1 + as.symbol(humanmetrics)[i] + All_HardCoral + Rugosity | location, data=sem.dat.scaled), 
#                      lme(form_b, random = ~ 1 + as.symbol(humanmetrics)[i] | location, data=sem.dat.scaled),
#                      lme(form_c, random = ~ 1 + All_HardCoral + as.symbol(humanmetrics)[i] | location, data=sem.dat.scaled))


#setwd(outdir)
#txtname<-paste("stats_wakatobiSEM_ecologyOnly-randomIslandSlopesAndIntercepts_", humanmetrics[i], ".txt", sep="")
#sink(txtname)
#print(summary(wakarand_slope_intercept.psem, .progressBar = F))
#sink()





# NEXT consider converting biomass (kg) to biomass/ha (or some other standardized unit)
# NEXT Look up lme syntax for how to model hierarchies - e.g., sites within islands
# NEXT Look up lme syntax for how to model random intercepts and slopes




### Below is script for processing IPB data:
# Import IPB fish data
#ecodat<-read.csv("IPB_Wakatobi_2016.csv")

# ADD ISLAND ID Column
#ecodat$ISLAND<-as.character(ecodat$site) # must be character vector for strsplit to work
#island.split<-unlist(strsplit(ecodat$ISLAND, split="_")) 
#ecodat$ISLAND<-island.split[seq(from=1, to=length(island.split), by=2)] #Only keep island names
#ecodat$ISLAND<-as.factor(ecodat$ISLAND)

#setwd("~/Analyses_notGit/fish-otakotak/indo-dat/Wakatobi/_coralCSVs")
#coral.files<-list.files()

#for (i in 1:length(coral.files))
#{
#  setwd("~/Analyses_notGit/fish-otakotak/indo-dat/Wakatobi/_coralCSVs")
#  corali<-read.csv(coral.files[i])
#  # identify "Transisi" columns
#  trans_cols<-grep("Transisi", names(corali))
#    for (j in trans_cols)
#    {
#      subtrahend<-corali[,j] # minuend - subtrahend = difference
#      minuend<-corali[2:length(corali[,j]), j]
#      minuend<-append(minuend, NA) # Add NA so vectors are the same length
#      difference<-minuend-subtrahend
#      difference<-difference[-(length(difference))] # Remove NA so vector has same length as replacement
#      corali[2:length(corali[,j]), j]<-difference
#    }
#  csvsplit<-unlist(strsplit(coral.files[i], split=".", fixed=TRUE))
#  csvname=paste(csvsplit[1], "-Lengths", ".csv", sep="")
#  setwd(outdir)
#  write.csv(corali, csvname, row.names=FALSE)
#}

# Read in processed coral data: 1 - calculate benthic cover and 2 - merge with fish data
#all.files<-list.files()
#length.files<-all.files[grep("Lengths", all.files)]

#for (i in 1:length(length.files))
#{
#  setwd(outdir)
#  lengthi<-read.csv(length.files[i])
#  lifeform_cols<-grep("Lifeform", names(lengthi))
# onecol_lifeforms<-NA
#  for (j in lifeform_cols)
#  {
#    lengthi[j]<-as.factor(trimws(lengthi[,j])) # Remove trailing/leading white space for life form columns (i.e., fix spreadsheet problems)
#    onecol_lifeforms<-append(onecol_lifeforms, as.character(lengthi[,j])) # LEFT OFF HERE
#  }
#}


