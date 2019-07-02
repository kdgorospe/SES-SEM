# PATH ANALYSIS of Wakatobi SES
rm(list=ls())

#update.packages(ask = FALSE, checkBuilt = TRUE)

library(piecewiseSEM)
library(nlme) # Version 3.1.122
library(ggplot2)
library(corrplot)
library(car) #VIF calculations


# Files live locally, not on GoogleDrive
indir<-"~/Analyses_notGit/fish-otakotak/indo-dat/Wakatobi"
outdir<-"~/Analyses_notGit/_RESULTS/fish-otakotak"


# FIRST: input / munge fish data
setwd(indir)
fishdat<-read.csv("_fishData/cleaned_wakatobi_fish_uvc.csv")

###### DONE: Removed site 17==Sombano from analysis

# Summarize total fish biomass at site level (site_id column)
fish.site<-aggregate(biomass_g ~ site_id, data=fishdat, FUN=sum)

# Is there a difference? First calculate average within transect, and then average across three transects per site
# no difference: average of averages equivalent to global average
#fish.site1<-aggregate(biomass_g ~ site_id + transect, data=fishdat, FUN=sum)
#fish.site<-aggregate(biomass_g ~ site_id, data=fish.site1, FUN=sum)

# Summarize fish biomass by functional group at site level
fishfun.tmp<-aggregate(biomass_g ~ site_id + trophic_group, data=fishdat, FUN=sum)
### OTHER POTENTIAL RESPONSE VARIABLES: DIVERSITY METRICS e.g., species richness and functional trait diversity

# For now, focus on total fish biomass for analysis: 

# NEXT: input / munge coral cover data
coraldat<-read.csv("_coralData/raw_coral_data_wakatobi_may_2018.csv")
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


# Process rugosity data
setwd(indir)
rugdat<-read.csv("_coralData/raw_rugosity_data_wakatobi_may_2018.csv")
rug.site<-aggregate(Rugosity ~ Site.Name, data=rugdat, FUN=mean)
setwd(outdir)
write.csv(rug.site, "data_wakatobi_benthicRugosity.csv", quote=FALSE, row.names=FALSE)


# READ-IN HUMAN POPULATION METRICS DATA:
setwd(indir)
#distWeighted.dat<-read.csv("_humanPopData/data_wakatobiHumans_distanceWeighted.csv") # LEAVE THESE OUT FOR NOW
#distToLandings.dat<-read.csv("_humanPopData/data_wakatobiHumans_distanceToLandingsSite.csv") # LEAVE THESE OUT FOR NOW
humanDensity.dat<-read.csv("_humanPopData/data_wakatobiHumans_areaWeightedDensityMetrics_5_km_buffer.csv") # weights each village's population density by its area to get "total population" within 5km buffer



# READ-IN MSEC oceanographic (and other) variables:
setwd(indir)
msec.dat<-read.csv("_MSECData/msec_out_5km.csv")

# READ-IN SST data 
setwd(indir)
sst.dat<-read.csv("_SSTsummaries/Wakatobi_2018_SSTExtract.csv")

# NEXT: Merge fish, benthic, rugosity, human population data using "site journal.xlsx" as site key
setwd(indir)
site.key<-read.csv("site journal-CLEANED-siteNames-removedsite17-decimalDegrees-meanVisibility.csv")
site.key<-subset(site.key, select=c("site_id", "Site.Name", "lat_dd", "long_dd", "exposed", "u_visibility", "type_reef", "location"))

# Merge all data: 
##### Do this in the following order: fish, MSEC, human pop data, rugosity, benthic cover

# Merge fish data
dat.tmp<-merge(site.key, fish.site, by="site_id")

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

# Insert 0 for missing human pop data
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
hist(alldat.site[,"biomass_g"],xlab="Total Biomass", main="Histogram of Site-Level Total Fish Biomass")
dev.off()

# Try log biomass
pdf(file="plot_histogram.LOGtotalbiomass.pdf")
hist(log10(alldat.site[,"biomass_g"]),xlab="log Total Biomass", main="Histogram of Site-Level log Fish Biomass")
dev.off()

# Include logbiomass as a potential response variable in data.frame
alldat.site$log_biomass_g<-log10(alldat.site[,"biomass_g"])

# QUESTION: how sensitive is SEM to (response) variable(s) having normal distribution?

# NEXT: create scatterplots

# Select all response + predictor variables + hierarchical variables
scatter.final<-alldat.site[,-c(1:4)]
loc.col<-grep("location", names(scatter.final))
bio.col<-grep("biomass_g", names(scatter.final))
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
  newfile=paste("plot_scatter_totalbiomass_v_", scatter.names[i], ".pdf", sep="")
  p<-ggplot(data=scatter.final, aes(x=get(scatter.names[i]), y=biomass_g)) + 
    geom_point(aes(x=get(scatter.names[i]), y=biomass_g, color=location, shape=type_reef), size=2) +
    labs(y="Total Biomass (g)", x=scatter.titles[i]) +  
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
  newfile=paste("plot_scatter_LOGtotalbiomass_v_", scatter.names[i], ".pdf", sep="")
  p<-ggplot(data=scatter.final, aes(x=get(scatter.names[i]), y=log_biomass_g)) + 
    geom_point(aes(x=get(scatter.names[i]), y=log_biomass_g, color=location, shape=type_reef), size=2) +
    labs(y="log Biomass (g)", x=scatter.titles[i]) +  
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


