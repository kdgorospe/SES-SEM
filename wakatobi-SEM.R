# PATH ANALYSIS of Wakatobi SES
rm(list=ls())
library(piecewiseSEM)
library(tidyverse)
library(nlme) # Version 3.1.122
library(ggplot2)
library(corrplot)
library(car) #VIF calculations
library(vegan)

# Read-in organized data
setwd("~/Analyses/_RESULTS/SES-SEM/__Organized Data")
#alldat.site<-read.csv("OriginalLandingsData/SiteLevelHumans-5km/data_wakatobi_allDataMerged.csv")
#alldat.site<-read.csv("LandingsDataPerKm2/SiteLevelHumans-5km/data_wakatobi_allDataMerged.csv")
#alldat.site<-read.csv("QC-LandingsData/SiteLevelHumans-5km/data_wakatobi_allDataMerged.csv")

# Alternatively, for aggregated analysis, read in fishing ground-level data
alldat.site<-read.csv("QC-LandingsData/SiteLevelHumans-5km/data_wakatobi_allDataMerged-fishingGroundLevel.csv")
alldat.site<-read.csv("OriginalLandingsData/SiteLevelHumans-5km/data_wakatobi_allDataMerged-fishingGroundLevel.csv")


# Set output directory
outdir<-"~/Analyses/_RESULTS/SES-SEM/"

################################################################################################################
## The following adds flexibility to choose different response variables (fish ecology metrics) for the SEM
################################################################################################################
responseDF<-as.data.frame(cbind(fish.response=c("log_biomass_g", "biomass_g", "avg_size", "no_of_species", "shannon", "invsimpson", "evenness"),
                                fish.title=c("log Total Biomass (g)", "Total Biomass (g)", "Average Size (cm)", "Richness", "Shannon Diversity (H')", "Inverse Simpson's Diversity (D2)", "Simpson's Evenness (E)" )))
## Identify fish response column name and title
## CHOICES:
## for biomass, set as biomass_g
## for avg size, set as avg_size
## for log biomass, set as log_biomass_g
## for species richness, set as no_of_species
## for shannon diversity, set as shannon
## for inverse simpson's, set as invsimpson
## for simpson's evenness, set as SimpsonEvenness
fish.col<-"biomass_g" # Set response here
fish.row<-responseDF$fish.response %in% fish.col
fish.title<-as.character(responseDF[fish.row, "fish.title"])

# NEXT: Create df for scatter plots (scatter.final) with all response + predictor variables + hierarchical variables
otherresponse<-responseDF$fish.response[!responseDF$fish.response %in% fish.col]
nonpredictors<-names(alldat.site) %in% c("site_name", as.character(otherresponse))
scatter.final<-alldat.site[,!nonpredictors]
loc.col<-grep("location", names(scatter.final))
response.col<-grep(fish.col, names(scatter.final))
scatter.names<-names(scatter.final)[-c(loc.col, response.col)]
########################################################################################
########################################################################################
########################################################################################
# Set graph names here: should match object scatter.names
scatter.titles<-c( # Site metadata
                  "Reef Type",  
                  
                  # Landings data
                  "Total Landings", "Total Personal",  "Total On-Island", "Total Papalele", "Total Pengumpul", "Total Market",
                  "Mean Landings", "Mean Personal", "Mean On-Island", "Mean Papalele", "Mean Pengumpul", "Mean Market", 
                  "Proportion Personal", "Proportion On-Island", "Proportion Papalele", "Proportion Pengumpul", "Proportion Market", 
                  
                  # Benthic cover data
                  "Rugosity",
                  "Dead Coral With Algae", "Macroalgae", "Rubble", "Rock", "Soft Coral", "Sand", 
                  "All Hard Coral", "All Abiotic Substrate",
                  
                  # Oceanography data (from MSEC)
                  "NPP (Mean)", "NPP (Min)", "NPP (Max)", "NPP (SD)", "NPP (Interannual SD)",
                  "Reef Area within 5km", 
                  "Wave Energy (Mean)", "Wave Energy (SD)", "Wave Energy (Interannual SD)", 
                  "Wind Fetch", 

                  # SST data
                  "SST SD", "SST 50th Percentile", "SST 98th Percentile", "SST 2nd Percentile", "SST Kurtosis", "SST Skewness",
                  
                  # Human Population data
                  "Population", "No of Fishers", "Rowboats", "Motorboats"
                  )

# Create scatterplots
setwd(outdir )
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
discreetvar<-match(c("wave_wind_fetch", "type_reef"), scatter.names)
corr.names<-scatter.names[-discreetvar]
corr.final<-scatter.final[corr.names]

##### Rename rows and columns for correlation plots
discreetvar.titles<-match(c("Wind Fetch", "Reef Type"), scatter.titles)
corr.titles<-scatter.titles[-discreetvar.titles]
colnames(corr.final)<-corr.titles


cor.dat<-cor(corr.final)
setwd(outdir)
write.csv(cor.dat, file="_Table_CorrelationsPearson.csv")
cor.test<-abs(cor.dat)>0.5
write.csv(cor.test, file="_Table_CorrelationsTestPearson.csv")

# Generate p values and confidence intervals for each correlation pair
pvals<-cor.mtest(corr.final, conf.level=0.95)

pdf(file="_Figure_CorrelationVisualPearson.pdf")
corrplot(cor.dat, method="color", tl.col="black", tl.cex=0.7, number.cex=0.4, p.mat=pvals$p, sig.level=0.05, insig="blank", cl.align.text="r", addgrid.col="grey")
dev.off()

pdf(file="_Figure_CorrelationVisualPearson-clustered.pdf")
corrplot(cor.dat, order="hclust", method="color", tl.col="black", tl.cex=0.7, number.cex=0.4, p.mat=pvals$p, sig.level=0.05, insig="blank", cl.align.text="r", addgrid.col="grey")
dev.off()


################################################################################
################################################################################
# STRUCTURAL EQUATION MODELS - SITE LEVEL ANALYSIS
################################################################################
################################################################################

### First, test for MULTICOLLINEARITY (calculate VIF)
### Equations for MULTICOLLINEARITY can be recycled as PSEM equations
analysis.col<-grep(fish.col, names(alldat.site))

### TWO DIFFERENT FRAMEWORK:
### (1) total landings variable nested between market and fish response variable
### (2) total landings split into two market vs nonmarket landings

### FRAMEWORK (1): 

### NOTES: predicting landings_sum_tot with fishing ground effects gives error: computationally singular (i.e., groupings perfectly predict response variable)
### Essentially, this is the problem https://stackoverflow.com/questions/25752259/error-in-nlme-repeated-measures
### Try: plot(alldat.site$landings_prop_market, alldat.site$All_HardCoral)
### i.e., - if a fishing ground variable is the response variable, can't use fishing ground as a random effect
form1a<-as.formula(paste(names(alldat.site)[analysis.col], " ~ All_HardCoral + reef_area_5km + landings_sum_tot", sep=""))
form1b<-as.formula("All_HardCoral ~  Population_2017")
form1c<-as.formula("reef_area_5km ~  Population_2017")
form1d<-as.formula("landings_sum_tot ~ landings_prop_market")

fit1a <- lm(form1a, data=alldat.site)
#fit1b <- lm(form1b, data=alldat.site)

# Note: only need to calculate vif for formulas with at least two predictors
vif(fit1a)
#vif(fit1b)

### SUBSET only model variables from alldat
all.forms_1<-ls(pattern="form1")
all.vars_1<-all.vars(get(all.forms_1[1]))
for (i in 2:length(all.forms_1))
{
  vars_i<-all.vars(get(all.forms_1[i]))
  all.vars_1<-append(all.vars_1, vars_i)    
}
model.vars_1<-unique(all.vars_1)

# Create site-level dataframe
sem.vars.site<-alldat.site[c("location", "type_reef", model.vars_1)]
sem.site.scaled<-as.data.frame(apply(sem.vars.site[model.vars_1], 2, scale))
sem.site.scaled<-cbind(sem.vars.site$location, sem.vars.site$type_reef, sem.site.scaled)
names(sem.site.scaled)[1:2]<- c("location", "reef_type")


### PSEM only fixed effects
waka.sitelevel.psem<-psem(lm(form1a, data=sem.site.scaled), 
                lm(form1b, data=sem.site.scaled),
                lm(form1c, data=sem.site.scaled),
                lm(form1d, data=sem.site.scaled))

setwd(outdir)
txtname<-paste("stats_wakatobiSEM_siteLevelData_", fish.col, ".txt", sep="")
sink(txtname)
print(summary(waka.sitelevel.psem, .progressBar = F))
sink()

### PSEM with fishing ground effects on fish response variable ONLY
waka.sitelevel.groundEffects.psem<-psem(lme(form1a, random = ~ 1 | location, data=sem.site.scaled), 
                                        lm(form1b,  data=sem.site.scaled),
                                        lm(form1c,  data=sem.site.scaled),
                                        lm(form1d, data=sem.site.scaled))

setwd(outdir)
txtname<-paste("stats_wakatobiSEM_siteLevelData_", fish.col, "_groundEffectsOnly.txt", sep="")
sink(txtname)
print(summary(waka.sitelevel.groundEffects.psem, .progressBar = F))
sink()


### PSEM with fishing ground effects on fish response variable, and reef type effects on all other variables
waka.sitelevel.groundAndReefEffects.psem<-psem(lme(form1a, random = ~ 1 | location, data=sem.site.scaled), 
                                   lme(form1b,  random = ~ 1 | reef_type, data=sem.site.scaled),
                                   lme(form1c,  random = ~ 1 | reef_type, data=sem.site.scaled),
                                   lme(form1d, random = ~ 1 | reef_type, data=sem.site.scaled))

setwd(outdir)
txtname<-paste("stats_wakatobiSEM_siteLevelData_", fish.col, "_groundAndReefEffects.txt", sep="")
sink(txtname)
print(summary(waka.sitelevel.groundAndReefEffects.psem, .progressBar = F))
sink()

### PSEM with NESTED effects on fish response (fishing ground nested within reef type), and no other random effects
### SYNTAX: if location is nested within reef type then random = ~ 1 | reef_type/location
waka.sitelevel.nestedEffects.psem<-psem(lme(form1a, random = ~ 1 | reef_type/location, data=sem.site.scaled), 
                                        lm(form1b, data=sem.site.scaled),
                                        lm(form1c, data=sem.site.scaled),
                                        lm(form1d, data=sem.site.scaled))

setwd(outdir)
txtname<-paste("stats_wakatobiSEM_siteLevelData_", fish.col, "_nestedEffectsOnly.txt", sep="")
sink(txtname)
print(summary(waka.sitelevel.nestedEffects.psem, .progressBar = F))
sink()



### PSEM with NESTED effects on fish response (fishing ground nested within reef type), and reef type effects on all other variables
### SYNTAX: if location is nested within reef type then random = ~ 1 | reef_type/location
waka.sitelevel.nestedAndReefEffects.psem<-psem(lme(form1a, random = ~ 1 | reef_type/location, data=sem.site.scaled), 
                                        lme(form1b, random = ~ 1 | reef_type, data=sem.site.scaled),
                                        lme(form1c, random = ~ 1 | reef_type, data=sem.site.scaled),
                                        lme(form1d, random = ~ 1 | reef_type, data=sem.site.scaled))

setwd(outdir)
txtname<-paste("stats_wakatobiSEM_siteLevelData_", fish.col, "_nestedAndReefEffects.txt", sep="")
sink(txtname)
print(summary(waka.sitelevel.nestedAndReefEffects.psem, .progressBar = F))
sink()




################################################################################
################################################################################
# STRUCTURAL EQUATION MODELS - Fishing ground level analysis
################################################################################
################################################################################


### First, test for MULTICOLLINEARITY (calculate VIF)
### Equations for MULTICOLLINEARITY can be recycled as PSEM equations
analysis.col<-grep(fish.col, names(alldat.site))

### TWO DIFFERENT FRAMEWORK:
### (1) total landings variable nested between market and fish response variable
### (2) total landings split into two market vs nonmarket landings

### FRAMEWORK (1): 

### NOTES: predicting landings_sum_tot with fishing ground effects gives error: computationally singular (i.e., groupings perfectly predict response variable)
### Essentially, this is the problem https://stackoverflow.com/questions/25752259/error-in-nlme-repeated-measures
### Try: plot(alldat.site$landings_prop_market, alldat.site$All_HardCoral)
### i.e., - if a fishing ground variable is the response variable, can't use fishing ground as a random effect
form1a<-as.formula(paste(names(alldat.site)[analysis.col], " ~ All_HardCoral + reef_area_5km + landings_sum_tot", sep=""))
form1b<-as.formula("All_HardCoral ~  Population_2017")
form1c<-as.formula("reef_area_5km ~  Population_2017")
form1d<-as.formula("landings_sum_tot ~ landings_prop_market")

fit1a <- lm(form1a, data=alldat.site)
#fit1b <- lm(form1b, data=alldat.site)

# Note: only need to calculate vif for formulas with at least two predictors
vif(fit1a)
#vif(fit1b)

### SUBSET only model variables from alldat
all.forms_1<-ls(pattern="form1")
all.vars_1<-all.vars(get(all.forms_1[1]))
for (i in 2:length(all.forms_1))
{
  vars_i<-all.vars(get(all.forms_1[i]))
  all.vars_1<-append(all.vars_1, vars_i)    
}
model.vars_1<-unique(all.vars_1)

# Create site-level dataframe
sem.vars.site<-alldat.site[c("location", "type_reef", model.vars_1)]
sem.site.scaled<-as.data.frame(apply(sem.vars.site[model.vars_1], 2, scale))
sem.site.scaled<-cbind(sem.vars.site$location, sem.vars.site$type_reef, sem.site.scaled)
names(sem.site.scaled)[1:2]<- c("location", "reef_type")


### PSEM only fixed effects
waka.sitelevel.psem<-psem(lm(form1a, data=sem.site.scaled), 
                          lm(form1b, data=sem.site.scaled),
                          lm(form1c, data=sem.site.scaled),
                          lm(form1d, data=sem.site.scaled))

setwd(outdir)
txtname<-paste("stats_wakatobiSEM_fishingGroundLevelData_", fish.col, ".txt", sep="")
sink(txtname)
print(summary(waka.sitelevel.psem, .progressBar = F))
sink()



### PSEM with reef type effects on fish response variable ONLY
waka.sitelevel.groundEffects.psem<-psem(lme(form1a, random = ~ 1 | reef_type, data=sem.site.scaled), 
                                        lm(form1b,  data=sem.site.scaled),
                                        lm(form1c,  data=sem.site.scaled),
                                        lm(form1d, data=sem.site.scaled))

setwd(outdir)
txtname<-paste("stats_wakatobiSEM_fishingGroundLevelData_", fish.col, "_reefEffectsOnly.txt", sep="")
sink(txtname)
print(summary(waka.sitelevel.groundEffects.psem, .progressBar = F))
sink()


### PSEM with reef type effects on all response variables
waka.sitelevel.groundEffects.psem<-psem(lme(form1a, random = ~ 1 | reef_type, data=sem.site.scaled), 
                                        lme(form1b, random = ~ 1 | reef_type, data=sem.site.scaled),
                                        lme(form1c, random = ~ 1 | reef_type, data=sem.site.scaled),
                                        lme(form1d, random = ~ 1 | reef_type, data=sem.site.scaled))

setwd(outdir)
txtname<-paste("stats_wakatobiSEM_fishingGroundLevelData_", fish.col, "_reefEffectsOnlyAllVars.txt", sep="")
sink(txtname)
print(summary(waka.sitelevel.groundEffects.psem, .progressBar = F))
sink()
