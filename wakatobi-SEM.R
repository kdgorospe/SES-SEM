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

# Use logbiomass as response variable in data.frame
fish.mass$biomass_g<-log10(fish.mass[,"biomass_g"])
names(fish.mass)[2]<-"log_biomass_g"

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



##### FISH FLOW CLEANING:
##### 1 - check for duplicate trips

# First remove all duplicate rows
dim(trip.dat)
trip.dat<-unique(trip.dat)
dim(trip.dat)

tripID.tab<-as.data.frame(table(trip.dat$trip_id))
tripID.check<-tripID.tab[(tripID.tab$Freq!=1),]$Var1

# "1 11 DAHLAN" and "1 11 ACER + DULETES" are both duplicated rows (only difference is that some columns are NAs and others are blank)
delete1<-grep("1 11 DAHLAN", trip.dat$trip_id)[1]
trip.dat<-trip.dat[-delete1,]

delete2<-grep("1 11 ACER", trip.dat$trip_id)[1]
trip.dat<-trip.dat[-delete2,]
# Now, only "11 20 BAGON" is a duplicated tripID

##### FISH FLOW CLEANING (contd):
##### 2 - Need to clean column: "landings_sold_personally_no"; inputted as "factor" and should be "numeric"
summary(trip.dat)
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



##### FISH FLOW CLEANING (contd):
##### 3 - Fix spelling errors and standardize all landing_units: 
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

##### FISH FLOW CLEANING (contd):
##### 4 - check which fish flow data are missing (is.empty tests for 0s and NAs)
# First convert all blank spaces to NAs
trip.dat[trip.dat==""]<-NA

# Convert all NAs in "no." columns to 0s
no_cols<-grep("_no", names(trip.dat))
for (i in no_cols)
{
  trip.dat[,i][is.na(trip.dat[,i])]<-0
}

# Which trips have "0" fish flow information? (columns: sold, eaten, or given)
# all columns except "landings_no"
flow.dat<-trip.dat[no_cols]
notflow<-grep("landing_no", names(flow.dat))
flow.dat<-flow.dat[,-notflow]
trip.dat$FlowSums<-rowSums(flow.dat)

## 10 trips have no fish flow data (throw these out)
trip.dat<-trip.dat[trip.dat$FlowSums!=0,]


##### FISH FLOW CLEANING (contd):
##### 5 - Quality Control: Check that sum of fish flows = total landings
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
  
  # for data already in units of fish abundance (i.e., unit column == "fish"), insert the number 1 (since units are already in "fish", number of units x 1 = number of fish)
  trip.dat$newcol[grep("fish", trip.dat[,i])]<-1
  
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
trip.dat$fishflow_abund <- trip.dat$landings_sold_personally_abund + 
  trip.dat$landings_sold_Papalele_abund + 
  trip.dat$landings_sold_Pengumpul_abund + 
  trip.dat$landings_eaten_abund + 
  trip.dat$landings_given_abund

fishflowQC<-trip.dat$fishflow_abund != trip.dat$landings_abund
trip.dat_needsQC<-trip.dat[fishflowQC,]

# These are the trip IDs that need fish flow QC
tripID_needsQC<-trip.dat_needsQC$trip_id

# How much do landings abundance data differ from sum of fish flow data?
summary(abs(trip.dat_needsQC$landings_abund - trip.dat_needsQC$fishflow_abund))
# Median: 19.5
# Median 72.6

# Without absolute value:
summary(trip.dat_needsQC$landings_abund - trip.dat_needsQC$fishflow_abund)
# On average, landings data were 11.477 greater than sum of fish flow data (Median: -11.750)

summary(trip.dat)
# If we keep these in the dataset, majority of fish flow data are zeroes

# For now, REMOVE THESE:
trip.dat<-trip.dat[!(trip.dat$trip_id %in% tripID_needsQC),]
summary(abs(trip.dat$landings_abund - trip.dat$fishflow_abund))
summary(trip.dat)

##### FISH FLOW CLEANING (contd):
##### 6 - Quality Control: Check for outliers
plot(trip.dat$landings_abund)
plot(trip.dat$landings_sold_personally_abund)
plot(trip.dat$landings_sold_Papalele_abund) # ON ISLAND
plot(trip.dat$landings_sold_Pengumpul_abund) # OFF ISLAND
# How many trips had catch sold off-island?
sum(trip.dat$landings_sold_Pengumpul_abund!=0) # 32 out of 211
plot(trip.dat$landings_eaten_abund)
plot(trip.dat$landings_given_abund)

# Trim down dataset
trip.dat<-subset(trip.dat, select=c(trip_id, fishing_grnd1,
                          landings_abund, 
                          landings_sold_personally_abund,
                          landings_sold_Papalele_abund,
                          landings_sold_Pengumpul_abund,
                          landings_eaten_abund,
                          landings_given_abund))

##### END FISH FLOW CLEANING
#### Export CSV?


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
landings_sumtot<-aggregate(trip.dat$landings_abund ~ trip.dat$new_fg, FUN = sum )
names(landings_sumtot)<-c("location", "landings_sum_tot")
# Doesn't seem like an appropriate driver, needs to be normalized by area of fishing ground

landings_meantot<-aggregate(trip.dat$landings_abund ~ trip.dat$new_fg, FUN = mean )
names(landings_meantot)<-c("location", "landings_mean_tot")
table(trip.dat$new_fg)
# Better measure of fishing pressure than sum total?
# If area of fishing ground is proportional to number of trips, then yes
# Or if assumption is that fishing EFFORT is equal everywhere (Badjao fish everywhere equally), then yes 

# How do landings eaten/given vs sold (to anyone) affect fish response
# Note: ZERO landings eaten or given in this dataset

# How do landings more-directly benefiting fishers (eaten, given, sold_personally - ie, no middle-person?) affect fish response
landings_personal<-aggregate(trip.dat$landings_sold_personally_abund +
                         trip.dat$landings_eaten_abund +
                         trip.dat$landings_given_abund ~ trip.dat$new_fg, FUN = sum )
names(landings_personal)<-c("location", "landings_sum_personal")
landings_personal<-merge(landings_personal, landings_sumtot, by="location")
landings_personal$landings_personal_prop<-landings_personal$landings_sum_personal / landings_personal$landings_sum_tot


# How do landings sold to pengumpul affect fish response?
landings_pengumpul<-aggregate(trip.dat$landings_sold_Pengumpul_abund ~ trip.dat$new_fg, FUN = sum )
names(landings_pengumpul)<-c("location", "landings_sum_pengumpul")
landings_pengumpul<-merge(landings_pengumpul, landings_sumtot, by="location")
landings_pengumpul$landings_pengumpul_prop<-landings_pengumpul$landings_sum_pengumpul / landings_pengumpul$landings_sum_tot
# Pengumpul = sold OFF-ISLAND

# How do landings sold to papalele affect fish response?
landings_papalele<-aggregate(trip.dat$landings_sold_Papalele_abund ~ trip.dat$new_fg, FUN = sum )
names(landings_papalele)<-c("location", "landings_sum_papalele")
landings_papalele<-merge(landings_papalele, landings_sumtot, by="location")
landings_papalele$landings_papalele_prop<-landings_papalele$landings_sum_papalele / landings_papalele$landings_sum_tot
# Papalele = sold ON-ISLAND

# How do landings sold to pengumpul or papalele affect fish response?
landings_market<-aggregate(trip.dat$landings_sold_Papalele_abund +
                             trip.dat$landings_sold_Pengumpul_abund ~ trip.dat$new_fg, FUN=sum)
names(landings_market)<-c("location", "landings_sum_market")
landings_market<-merge(landings_market, landings_sumtot, by="location")
landings_market$landings_market_prop<-landings_market$landings_sum_market / landings_market$landings_sum_tot 

landings.dat<-merge(landings_meantot, landings_personal, by="location")
landings.dat<-merge(landings.dat, landings_pengumpul, by="location")
landings.dat<-merge(landings.dat, landings_papalele, by="location")
landings.dat<-merge(landings.dat, landings_market, by="location")
landings.dat<-subset(landings.dat, select=c(location, landings_mean_tot, landings_personal_prop, landings_pengumpul_prop, landings_papalele_prop, landings_market_prop))

###### Merge fish, oceanographic (MSEC), human pop data, rugosity, benthic cover, SST AND catch data using "site journal.xlsx" as site key: https://drive.google.com/open?id=1SNHtCmszbl6SYMPng1RLCDQVmap3e27n
drive_download(as_id("1SNHtCmszbl6SYMPng1RLCDQVmap3e27n"), overwrite=TRUE) # Saves file to working directory 
site.key<-read.csv("site journal-CLEANED-siteNames-removedsite17-decimalDegrees-meanVisibility.csv")
file.remove("site journal-CLEANED-siteNames-removedsite17-decimalDegrees-meanVisibility.csv")
site.key<-subset(site.key, select=c("site_id", "Site.Name", "lat_dd", "long_dd", "exposed", "u_visibility", "type_reef", "location"))

################################################################################################################
################################################################################################################




responseDF<-as.data.frame(cbind(fish.response=c("fish.mass", "fish.rich", "fish.shan", "fish.isim", "fish.even"),
                                fish.col=c("log_biomass_g", "no_of_species", "shannon", "invsimpson", "SimpsonEvenness"),
                                fish.title=c("log Total Biomass (g)", "Richness", "Shannon Diversity (H')", "Inverse Simpson's Diversity (D2)", "Simpson's Evenness (E)" )))


#### Next identify fish response object name, column name, and title
## CHOICES:
## for biomass, set as fish.mass 
## for species richness, set as fish.rich 
## for shannon diversity, set as fish.shan
## for inverse simpson's, set as fish.isim
## for simpson's evenness, set as fish.even
fish.response<-"fish.even" # Set response here
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
                  "Landings per trip", "Personal", "Pengumpul", "Papalele", "Market",
                  
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

form1<-as.formula(paste(names(alldat.site)[analysis.col], " ~ All_HardCoral + Population_2017 + landings_mean_tot", sep=""))
form2<-as.formula("All_HardCoral ~ Population_2017")
form3<-as.formula("landings_mean_tot ~ landings_market_prop")




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
sem.vars<-c(fish.col, "location", "All_HardCoral", "landings_mean_tot", "Population_2017", "landings_market_prop")
sem.vars.site<-alldat.site[sem.vars]
sem.dat.site<-subset(sem.vars.site, select=-location)

# Aggregate to fishing ground level, scale and center data:
sem.vars.ground<-aggregate(sem.dat.site, FUN=mean, by=list(alldat.site$location))
names(sem.vars.ground)[1]<-"location"
sem.dat.ground<-subset(sem.vars.ground, select=-location)
sem.dat.scaled<-as.data.frame(apply(sem.dat.ground, 2, scale))
sem.dat.scaled<-cbind(sem.vars.ground$location, sem.dat.scaled)
names(sem.dat.scaled)[1]<- "location"

waka.psem<-psem(lm(form1, data=sem.dat.scaled), 
                lm(form2, data=sem.dat.scaled),
                lm(form3, data=sem.dat.scaled))


basisSet(waka.psem)
# NOTE: So far, there are no independence claims so basis set = 0

setwd(outdir)
txtname<-paste("stats_wakatobiSEM_", fish.col, ".txt", sep="")
sink(txtname)
print(summary(waka.psem, .progressBar = F))
sink()

# coefficients should already be standardized since data was already scaled
#coefs(waka.psem, standardize="scale")

### LEFT OFF HERE: how to formulate hierarchical model?
### Incorporate random effects by location
# With island-level random effects, here random intercepts are modeled for each island 

#wakarandom.psem<-psem(lme(form1, random = ~ 1 | location, data=sem.dat.scaled), 
#                lme(form2, random = ~ 1 | location, data=sem.dat.scaled),
#                lme(form3, random = ~ 1 | location, data=sem.dat.scaled))

#setwd(outdir)
#txtname<-paste("stats_wakatobiSEM_", fish.col, "_randomIslandIntercepts.txt", sep="")
#sink(txtname)
#print(summary(wakarandom.psem, .progressBar = F))
#sink()

### New Code: Still need to try: Random slopes and interecepts for all predictors
#wakarand_slope_intercept.psem<-psem(lme(form_a, random = ~ 1 + as.symbol(humanmetrics)[i] + All_HardCoral + Rugosity | location, data=sem.dat.scaled), 
#                      lme(form_b, random = ~ 1 + as.symbol(humanmetrics)[i] | location, data=sem.dat.scaled),
#                      lme(form_c, random = ~ 1 + All_HardCoral + as.symbol(humanmetrics)[i] | location, data=sem.dat.scaled))


#setwd(outdir)
#txtname<-paste("stats_wakatobiSEM_ecologyOnly-randomIslandSlopesAndIntercepts_", humanmetrics[i], ".txt", sep="")
#sink(txtname)
#print(summary(wakarand_slope_intercept.psem, .progressBar = F))
#sink()


