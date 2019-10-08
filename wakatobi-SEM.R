# PATH ANALYSIS of Wakatobi SES
rm(list=ls())
library(piecewiseSEM)
library(nlme) # Version 3.1.122
library(ggplot2)
library(corrplot)
library(car) #VIF calculations
library(vegan)

# Read-in organized data
setwd("~/Analyses/_RESULTS/SES-SEM")
alldat.site<-read.csv("data_wakatobi_allDataMerged.csv")

# Set output directory
outdir<-"~/Analyses/_RESULTS/SES-SEM/"

################################################################################################################
## The following adds flexibility to choose different response variables (fish ecology metrics) for the SEM
################################################################################################################
responseDF<-as.data.frame(cbind(fish.response=c("log_biomass_g", "biomass_g", "no_of_species", "shannon", "invsimpson", "SimpsonEvenness"),
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
nonpredictors<-names(alldat.site) %in% c("site_id", "Site.Name", as.character(otherresponse))
scatter.final<-alldat.site[,!nonpredictors]
loc.col<-grep("location", names(scatter.final))
response.col<-grep(fish.col, names(scatter.final))
scatter.names<-names(scatter.final)[-c(loc.col, response.col)]
########################################################################################
########################################################################################
########################################################################################
# Set graph names here: should match object scatter.names
scatter.titles<-c( # site journal columns
                  "Latitude", "Longitude", "Exposure", "Visibility", "Reef Type",  
                  
                  # Landings dat
                  "Landings per Trip", "Total Landings", 
                  "Total Personal", "Proportion Personal",
                  "Total Pengumpul", "Proportion Pengumpul", "Mean Pengumpul",
                  "Total Papalele", "Proportion Papalele",
                  "Total Market", "Proportion Market",
                  "Total On-Island", "Proportion On-Island", "Mean On-Island",
                  
                  # MSEC columns
                  "Mean NPP", "Min NPP", "Max NPP", "NPP SD", "Interannual NPP SD",
                  "Reef Area within 5km", 
                  "Mean Wave Energy", "SD of Wave Energy", "Interannual SD of Wave Energy", "Wind Fetch", 
                  
                  # getWakatobiPop columns
                  "Population within 2.5km", "Fishers within 2.5km", "Rowboats within 2.5km", "Motorboats within 2.5km",
                  
                  # benthic columns
                  "Rugosity",
                  "Dead Coral With Algae", "Macroalgae", "Rubble", "Rock", "Soft Coral", "Sand", 
                  "All Hard Coral", "All Abiotic Substrate",
                  
                  # SST columns
                  "SST SD", "SST 50th Percentile", "SST 98th Percentile", "SST 2nd Percentile", "SST Kurtosis", "SST Skewness"
                  )

# Create scatterplots
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

# Spearman's for all ranked, ordered vars (e.g., DACOR data)
cor.spear<-cor(corr.final, method="spearman")
write.csv(cor.spear, file="_Table_CorrelationsSpearman.csv")
cor.spear.test<-abs(cor.spear)>0.5
write.csv(cor.spear.test, file="_Table_CorrelationsTestSpearman.csv")

# Generate p values and confidence intervals for each correlation pair
pvals<-cor.mtest(corr.final, conf.level=0.95)

pdf(file="_Figure_CorrelationVisualPearson.pdf")
#corrplot.mixed(cor.dat, upper="circle", lower="number", tl.pos="lt", tl.col="black", tl.cex=0.7, lower.col="black", addCoefasPercent=TRUE, number.cex=0.7, p.mat=pvals$p, sig.level=0.05, insig="blank", diag="n")
corrplot(cor.dat, method="color", tl.col="black", tl.cex=0.7, number.cex=0.4, p.mat=pvals$p, sig.level=0.05, insig="blank", cl.align.text="r", addgrid.col="grey")
dev.off()

pdf(file="_Figure_CorrelationVisualSpearman.pdf")
#corrplot.mixed(cor.spear, upper="circle", lower="number", tl.pos="lt", tl.col="black", tl.cex=0.7, lower.col="black", addCoefasPercent=TRUE, number.cex=0.7, p.mat=pvals$p, sig.level=0.05, insig="blank", diag="y")
corrplot(cor.spear, method="circle", tl.col="black", tl.cex=0.7, number.cex=0.4, p.mat=pvals$p, sig.level=0.05, insig="blank")
dev.off()


### Test for MULTICOLLINEARITY (calculate VIF)
### Equations for MULTICOLLINEARITY can be recycled as PSEM equations
analysis.col<-grep(fish.col, names(alldat.site))

### FIRST, construct simple model that only uses site-level data (i.e., no fish landings data)
form1a<-as.formula(paste(names(alldat.site)[analysis.col], " ~ All_HardCoral + reef_area_5km", sep=""))
form1b<-as.formula("All_HardCoral ~ SST_98perc")

fit1a <- lm(form1a, data=alldat.site)
fit1b <- lm(form1b, data=alldat.site)

# Note: if some path equations only contain one predictor, VIF test below is invalid
vif.test1a<-vif(fit1a)
vif.test1b<-vif(fit1b)


################################################################################
################################################################################
################################################################################
# STRUCTURAL EQUATION MODEL: 
# Site-level PSEM (i.e., no fish landings data)

### SUBSET only model variables from alldat
vars1a<-all.vars(form1a)
vars1b<-all.vars(form1b)
model.vars_1<-unique(c(vars1a, vars1b))

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


# Now, include fishing ground (i.e., fish landings) data
form2a<-as.formula(paste(names(alldat.site)[analysis.col], " ~ All_HardCoral + reef_area_5km + landings_sum_tot", sep=""))
form2b<-as.formula("All_HardCoral ~ SST_98perc")

fit2a <- lm(form2a, data=alldat.site)
fit2b <- lm(form2b, data=alldat.site)

# Note: if some path equations only contain one predictor, VIF test below is invalid
vif.test2a<-vif(fit2a)
vif.test2b<-vif(fit2b)

### SUBSET only model variables from alldat
vars2a<-all.vars(form2a)
vars2b<-all.vars(form2b)
model.vars_2<-unique(c(vars2a, vars2b))


sem.vars.catch<-alldat.site[c("location", "type_reef", model.vars_2)]
sem.catch.scaled<-as.data.frame(apply(sem.vars.catch[model.vars_2], 2, scale))
sem.catch.scaled<-cbind(sem.vars.catch$location, sem.vars.catch$type_reef, sem.catch.scaled)
names(sem.catch.scaled)[1:2]<- c("location", "reef_type")


waka.catch.psem<-psem(lm(form2a, data=sem.catch.scaled), 
                lm(form2b, data=sem.catch.scaled))


setwd(outdir)
txtname<-paste("stats_wakatobiSEM_withCatchData_", fish.col, ".txt", sep="")
sink(txtname)
print(summary(waka.catch.psem, .progressBar = F))
sink()


#cbind(sem.catch.scaled, 

## Catch data PSEM with hierarchy
waka.catch.reeftype.psem<-psem(lme(form2a, random = ~ 1 | reef_type, data=sem.catch.scaled), 
                                   lme(form2b, random = ~ 1 | reef_type, data=sem.catch.scaled))


setwd(outdir)
txtname<-paste("stats_wakatobiSEM_withCatchData_", fish.col, "_reefTypeEffects.txt", sep="")
sink(txtname)
print(summary(waka.catch.reeftype.psem, .progressBar = F))
sink()

### Notice: reef_type hierarchy here makes a difference

### NEXT: Add market effects (i.e., divide total landings into personal vs market)
form3a<-as.formula(paste(names(alldat.site)[analysis.col], " ~ All_HardCoral + reef_area_5km + landings_sum_market + landings_sum_personal", sep=""))
form3b<-as.formula("All_HardCoral ~ SST_98perc")

fit3a <- lm(form3a, data=alldat.site)


# Note: if some path equations only contain one predictor, VIF test below is invalid
vif.test3a<-vif(fit3a)


### SUBSET only model variables from alldat
vars3a<-all.vars(form3a)
vars3b<-all.vars(form3b)
model.vars_3<-unique(c(vars3a, vars3b))


sem.vars.market<-alldat.site[c("location", "type_reef", model.vars_3)]
sem.market.scaled<-as.data.frame(apply(sem.vars.market[model.vars_3], 2, scale))
sem.market.scaled<-cbind(sem.vars.market$location, sem.vars.market$type_reef, sem.market.scaled)
names(sem.market.scaled)[1:2]<- c("location", "reef_type")


waka.market.psem<-psem(lm(form3a, data=sem.market.scaled), 
                      lm(form3b, data=sem.market.scaled))


setwd(outdir)
txtname<-paste("stats_wakatobiSEM_withMarketData_", fish.col, ".txt", sep="")
sink(txtname)
print(summary(waka.market.psem, .progressBar = F))
sink()

## Market data PSEM with hierarchies
waka.market.reeftype.psem<-psem(lme(form3a, random = ~ 1 | reef_type, data=sem.market.scaled), 
                       lme(form3b, random = ~ 1 | reef_type, data=sem.market.scaled))


setwd(outdir)
txtname<-paste("stats_wakatobiSEM_withMarketData_", fish.col, "_reefTypeEffects.txt", sep="")
sink(txtname)
print(summary(waka.market.reeftype.psem, .progressBar = F))
sink()


## LEFT OFF HERE - create Rmd file for Lefcheck?
## Shouldn't hierarchies be based on fishing ground "location" AND reef type?
waka.market.hierarch.psem<-psem(lme(form3a, random = ~ 1 + landings_sum_market + landings_sum_personal | reef_type/location/location, data=sem.market.scaled), 
                                lme(form3b, random = ~ 1 | reef_type, data=sem.market.scaled))


setwd(outdir)
txtname<-paste("stats_wakatobiSEM_withMarketData_", fish.col, "_locationEffects.txt", sep="")
sink(txtname)
print(summary(waka.market.hierarch.psem, .progressBar = F))
sink()


## Catch data PSEM with reef_type AND fishing ground "location" hierarchies?
#waka.catch.allgroups.psem<-psem(lme(form_a, random = ~ 1 + landings_sum_tot | reef_type/location, data=sem.catch.scaled), 
#                               lme(form_b, random = ~ 1 | reef_type, data=sem.catch.scaled))


#setwd(outdir)
#txtname<-paste("stats_wakatobiSEM_withCatchData_", fish.col, "_reefTypeAndLocationEffects.txt", sep="")
#sink(txtname)
#print(summary(waka.catch.allgroups.psem, .progressBar = F))
#sink()



# coefficients should already be standardized since data was already scaled
#coefs(waka.psem, standardize="scale")


#### RANDOM EFFECTS:
### Use site-level (unaggregated) data and incorporate random effects by location
# i.e., start with sem.dat.site


## For examples on hierarchical model specification see:
## http://www.rensenieuwenhuis.nl/r-sessions-21-multilevel-model-specification-nlme/

## Random intercepts for site-level + 
## predictor that is allowed to vary by reef type (e.g., effect of coral and landings on fish) + 
## group-level predictor (e.g., effect of total landings by location)
#wakarandom.psem<-psem(lme(form1, random = ~ 1 + All_HardCoral + landings_mean_onisland | reef_type, data=sem.site.scaled), 
#                      lme(form2, random = ~ 1 + Population_2017 | location, data=sem.site.scaled),
#                      lme(form3, random = ~ 1 + Population_2017 | reef_type, data=sem.site.scaled))


#setwd(outdir)
#txtname<-paste("stats_wakatobiSEM_", fish.col, "_randomIntercepts.txt", sep="")
#sink(txtname)
#print(summary(wakarandom.psem, .progressBar = F))
#sink()

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


