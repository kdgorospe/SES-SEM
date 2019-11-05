# PATH ANALYSIS of Wakatobi SES
rm(list=ls())
library(piecewiseSEM)
library(nlme) # Version 3.1.122
library(ggplot2)
library(corrplot)
library(car) #VIF calculations
library(vegan)

# Read-in organized data
setwd("~/Analyses/_RESULTS/SES-SEM/__Organized Data")
alldat.site<-read.csv("data_wakatobi_allDataMerged.csv")

# Set output directory
outdir<-"~/Analyses/_RESULTS/SES-SEM/"

################################################################################################################
## The following adds flexibility to choose different response variables (fish ecology metrics) for the SEM
################################################################################################################
responseDF<-as.data.frame(cbind(fish.response=c("log_biomass_g", "biomass_g", "no_of_species", "shannon", "invsimpson", "evenness"),
                                fish.title=c("log Total Biomass (g)", "Total Biomass (g)", "Richness", "Shannon Diversity (H')", "Inverse Simpson's Diversity (D2)", "Simpson's Evenness (E)" )))
## Identify fish response column name and title
## CHOICES:
## for biomass, set as biomass_g
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
                  "Latitude", "Longitude", "Exposure", "Visibility", "Reef Type",  
                  
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
                  "Population", "Fishers", "Rowboats", "Motorboats"
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
discreetvar<-match(c("wave_wind_fetch", "type_reef", "exposed"), scatter.names)
corr.names<-scatter.names[-discreetvar]
corr.final<-scatter.final[corr.names]

##### Rename rows and columns for correlation plots
discreetvar.titles<-match(c("Wind Fetch", "Reef Type", "Exposure"), scatter.titles)
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


### Test for MULTICOLLINEARITY (calculate VIF)
### Equations for MULTICOLLINEARITY can be recycled as PSEM equations
analysis.col<-grep(fish.col, names(alldat.site))

### FIRST, construct simple model that only uses site-level data (i.e., no fish landings data)
form1a<-as.formula(paste(names(alldat.site)[analysis.col], " ~ All_HardCoral + reef_area_5km", sep=""))
form1b<-as.formula("All_HardCoral ~ SST_98perc + Population_2017")

fit1a <- lm(form1a, data=alldat.site)
fit1b <- lm(form1b, data=alldat.site)

# Note: only need to calculate vif for formulas with at least two predictors
vif(fit1a)
vif(fit1b)

################################################################################
################################################################################
# STRUCTURAL EQUATION MODELS
################################################################################
################################################################################
# PSEM of site-level data only (i.e., no fish landings data)

### SUBSET only model variables from alldat
all.forms_1<-ls(pattern="form1")
all.vars_1<-all.vars(get(all.forms_1[1]))
for (i in 2:length(all.forms_1))
{
  vars_i<-all.vars(get(all.forms_1[i]))
  all.vars_1<-append(all.vars_1, vars_i)    
}
model.vars_1<-unique(all.vars_1)

sem.vars.site<-alldat.site[c("location", "type_reef", model.vars_1)]
sem.site.scaled<-as.data.frame(apply(sem.vars.site[model.vars_1], 2, scale))
sem.site.scaled<-cbind(sem.vars.site$location, sem.vars.site$type_reef, sem.site.scaled)
names(sem.site.scaled)[1:2]<- c("location", "reef_type")

waka.sitelevel.psem<-psem(lm(form1a, data=sem.site.scaled), 
                lm(form1b, data=sem.site.scaled))

setwd(outdir)
txtname<-paste("stats_wakatobiSEM_siteLevelData_", fish.col, ".txt", sep="")
sink(txtname)
print(summary(waka.sitelevel.psem, .progressBar = F))
sink()


## Site-level PSEM with reef_type hierarchy
waka.sitelevel.reeftype.psem<-psem(lme(form1a, random = ~ 1 | reef_type, data=sem.site.scaled), 
                      lme(form1b, random = ~ 1 | reef_type, data=sem.site.scaled))


setwd(outdir)
txtname<-paste("stats_wakatobiSEM_siteLevelData_", fish.col, "_reefTypeEffects.txt", sep="")
sink(txtname)
print(summary(waka.sitelevel.reeftype.psem, .progressBar = F))
sink()

################################################################################
################################################################################
# PSEM including fishing ground data (for now use total fish landings)
form2a<-as.formula(paste(names(alldat.site)[analysis.col], " ~ All_HardCoral + reef_area_5km + landings_sum_tot", sep=""))
form2b<-as.formula("All_HardCoral ~ SST_98perc + Population_2017")
form2c<-as.formula("landings_sum_tot ~ Population_2017")

# Note: only need to calculate vif for formulas with at least two predictors
vif(lm(form2a, data=alldat.site))
vif(lm(form2b, data=alldat.site))

### SUBSET only model variables from alldat
all.forms_2<-ls(pattern="form2")
all.vars_2<-all.vars(get(all.forms_2[1]))
for (i in 2:length(all.forms_2))
{
  vars_i<-all.vars(get(all.forms_2[i]))
  all.vars_2<-append(all.vars_2, vars_i)    
}
model.vars_2<-unique(all.vars_2)


sem.vars.catch<-alldat.site[c("location", "type_reef", model.vars_2)]
sem.catch.scaled<-as.data.frame(apply(sem.vars.catch[model.vars_2], 2, scale))
sem.catch.scaled<-cbind(sem.vars.catch$location, sem.vars.catch$type_reef, sem.catch.scaled)
names(sem.catch.scaled)[1:2]<- c("location", "reef_type")


waka.catch.psem<-psem(lm(form2a, data=sem.catch.scaled), 
                lm(form2b, data=sem.catch.scaled), 
                lm(form2c, data=sem.catch.scaled))


setwd(outdir)
txtname<-paste("stats_wakatobiSEM_withCatchData_", fish.col, ".txt", sep="")
sink(txtname)
print(summary(waka.catch.psem, .progressBar = F))
sink()

## Catch data PSEM with hierarchy
waka.catch.reeftype.psem<-psem(lme(form2a, random = ~ 1 | reef_type, data=sem.catch.scaled), 
                               lme(form2b, random = ~ 1 | reef_type, data=sem.catch.scaled),
                               lme(form2c, random = ~ 1 | reef_type, data=sem.catch.scaled))


setwd(outdir)
txtname<-paste("stats_wakatobiSEM_withCatchData_", fish.col, "_reefTypeEffects.txt", sep="")
sink(txtname)
print(summary(waka.catch.reeftype.psem, .progressBar = F))
sink()
### Notice: reef_type hierarchy here makes a difference

################################################################################
################################################################################
### PSEM with market data (i.e., substitute total fish landings with personal vs market-destined fish landings)
form3a<-as.formula(paste(names(alldat.site)[analysis.col], " ~ All_HardCoral + reef_area_5km + landings_mean_market + landings_mean_personal", sep=""))
form3b<-as.formula("All_HardCoral ~ Population_2017")
form3c<-as.formula("landings_mean_market ~ Population_2017")
form3d<-as.formula("landings_mean_personal ~ Population_2017")

vif(lm(form3a, data=alldat.site))
vif(lm(form3b, data=alldat.site))

### SUBSET only model variables from alldat
all.forms_3<-ls(pattern="form3")
all.vars_3<-all.vars(get(all.forms_3[1]))
for (i in 2:length(all.forms_3))
{
  vars_i<-all.vars(get(all.forms_3[i]))
  all.vars_3<-append(all.vars_3, vars_i)    
}
model.vars_3<-unique(all.vars_3)


sem.vars.market<-alldat.site[c("location", "type_reef", model.vars_3)]
sem.market.scaled<-as.data.frame(apply(sem.vars.market[model.vars_3], 2, scale))
sem.market.scaled<-cbind(sem.vars.market$location, sem.vars.market$type_reef, sem.market.scaled)
names(sem.market.scaled)[1:2]<- c("location", "reef_type")


waka.market.psem<-psem(lm(form3a, data=sem.market.scaled), 
                      lm(form3b, data=sem.market.scaled),
                      lm(form3c, data=sem.market.scaled),
                      lm(form3d, data=sem.market.scaled))


setwd(outdir)
txtname<-paste("stats_wakatobiSEM_withMarketData_", fish.col, ".txt", sep="")
sink(txtname)
print(summary(waka.market.psem, .progressBar = F))
sink()

## Market data PSEM with hierarchies
waka.market.reeftype.psem<-psem(lme(form3a, random = ~ 1 | reef_type, data=sem.market.scaled),
                                lme(form3b, random = ~ 1 | reef_type, data=sem.market.scaled),
                                lme(form3c, random = ~ 1 | reef_type, data=sem.market.scaled),
                                lme(form3d, random = ~ 1 | reef_type, data=sem.market.scaled))


setwd(outdir)
txtname<-paste("stats_wakatobiSEM_withMarketData_", fish.col, "_reefTypeEffects.txt", sep="")
sink(txtname)
print(summary(waka.market.reeftype.psem, .progressBar = F))
sink()

## In reality, there should be different groupings for different drivers: 
## market and catch variables should be grouped by fishing ground "location" 
## ecological variables should be grouped by reef type
## BUT, so far, can't get this to converge
## TO ME, the code below says that for form3a, there is a random effect on intercept based on reef_type groups, while the slopes for market and personal are based on location groups
waka.market.hierarch.psem<-psem(lme(form3a, random = ~ 1 + landings_mean_market + landings_mean_personal | reef_type/location/location, data=sem.market.scaled), 
                                lme(form3b, random = ~ 1 | reef_type, data=sem.market.scaled),
                                lme(form3c, random = ~ 1 | location, data=sem.market.scaled),
                                lme(form3d, random = ~ 1 | location, data=sem.market.scaled))


setwd(outdir)
txtname<-paste("stats_wakatobiSEM_withMarketData_", fish.col, "_reefAndLocationEffects.txt", sep="")
sink(txtname)
print(summary(waka.market.hierarch.psem, .progressBar = F))
sink()
